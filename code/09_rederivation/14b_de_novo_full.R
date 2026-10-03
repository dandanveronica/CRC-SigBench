#!/usr/bin/env Rscript
# =====================================================================
# 14b_de_novo_full.R —— 自产签名对照（修正版）
#
# 【为什么要重做 14 的 Part A】
# 初版有两个会被审稿人抓住的缺陷：
#   1. GSE39582（579 例，最大的队列）作为训练队列时静默失败了，
#      导致「自产签名」这一层系统性缺少最有希望的一批 —— 天平被压歪；
#   2. 40 条自产记录因为基因数超出已发表签名覆盖的随机基线范围被剔除，
#      剩下的 160 条存在选择偏差。
# 本版两个都修：
#   - 训练阶段显式释放内存 + 报错可见，保证 6 个队列全部参与
#   - 对【每一个】(验证队列, 基因数) 组合，单独重算随机零分布
#
# 公平性的核心约定（与 14 一致，务必保持）：
#   自产签名只借用 LASSO 的「选基因」能力；打分环节强制使用
#   z-score 等权重求和 + 方向定向，与已发表签名完全相同。
#
# 用法：Rscript 14b_de_novo_full.R [B_null] [cores]
# 输出：audit_out_clean/de_novo_full.csv      自产签名 + 配对随机基线
#       audit_out_clean/null_matched.csv      配对随机基线明细
# =====================================================================

suppressPackageStartupMessages({
  library(survival); library(data.table); library(glmnet); library(parallel)
})
args  <- commandArgs(trailingOnly = TRUE)
B     <- if (length(args) >= 1) as.integer(args[1]) else 4000L
CORES <- if (length(args) >= 2) as.integer(args[2]) else 4L
COH <- "cohorts"; OUT <- "audit_out_clean"
K_TARGETS <- c(3, 5, 8, 10, 15, 20, 30, 50)
PRE_N <- 500L

zscore <- function(m) {
  mu <- rowMeans(m, na.rm = TRUE); sd <- apply(m, 1, sd, na.rm = TRUE)
  sd[is.na(sd) | sd < 1e-8] <- 1; (m - mu) / sd
}
cz <- function(sv, xbuf) {
  cf <- try(concordancefit(sv, xbuf, reverse = TRUE, std.err = FALSE), silent = TRUE)
  if (inherits(cf, "try-error")) return(NA_real_)
  c0 <- cf$concordance
  if (is.na(c0)) return(NA_real_)
  if (c0 < 0.5) 1 - c0 else c0
}
oriented_C <- function(zz, genes, time, status) {
  g <- intersect(genes, rownames(zz))
  if (length(g) < 2) return(NA_real_)
  p <- as.numeric(colMeans(zz[g, , drop = FALSE]))
  if (length(unique(p)) < 2) return(NA_real_)
  cz(Surv(time, status), matrix(p, ncol = 1))
}

cox_score <- function(z, time, status) {
  ord <- order(time, decreasing = TRUE)
  X <- z[, ord, drop = FALSE]
  t_o <- time[ord]; s_o <- status[ord]
  r <- rle(t_o); ends <- cumsum(r$lengths)
  ne <- vapply(seq_along(ends), function(k) {
    ii <- if (k == 1L) 1:ends[1] else (ends[k-1] + 1):ends[k]
    sum(s_o[ii] == 1)
  }, numeric(1))
  keep <- ne > 0
  if (!any(keep)) return(setNames(rep(0, nrow(z)), rownames(z)))
  ends <- ends[keep]; ne <- ne[keep]
  CX <- t(apply(X, 1, cumsum))
  CX2 <- t(apply(X * X, 1, cumsum))
  XS <- X; XS[, s_o != 1] <- 0            # 注意：不能写 X * (s_o==1)
  SXblk <- t(apply(XS, 1, cumsum)); rm(XS)
  U <- numeric(nrow(z)); V <- numeric(nrow(z))
  for (i in seq_along(ends)) {
    eu <- ends[i]
    S1 <- CX[, eu] / eu
    Se <- SXblk[, eu] - if (i == 1L) 0 else SXblk[, ends[i-1]]
    S2 <- CX2[, eu] / eu
    U <- U + Se - ne[i] * S1
    V <- V + ne[i] * pmax(S2 - S1 * S1, 1e-8)
  }
  rm(CX, CX2, SXblk, X); invisible(gc(FALSE))
  zv <- ifelse(V > 1e-12, U / sqrt(V), 0); zv[is.na(zv)] <- 0
  setNames(zv, rownames(z))
}

cat("=== 载入队列 ===\n")
CL <- list()
for (cf in list.files(COH, pattern = "\\.rds$", full.names = TRUE)) {
  ch <- sub("\\.rds$", "", basename(cf))
  d <- readRDS(cf)
  z <- zscore(d$expr)
  cl <- d$clin; cl <- cl[match(colnames(z), cl$sample), , drop = FALSE]
  tm <- suppressWarnings(as.numeric(cl$time)); st <- suppressWarnings(as.integer(cl$status))
  kp <- !is.na(tm) & !is.na(st)
  CL[[ch]] <- list(z = z[, kp, drop = FALSE], time = tm[kp], status = st[kp])
  cat(sprintf("  %-24s %4d 样本 %3d 事件\n", ch, sum(kp), sum(st[kp])))
  rm(d, z, cl); invisible(gc(FALSE))
}
cohorts <- names(CL)

# ---------------------------------------------------------------- A 构建
cat("\n=== A. 构建自产签名（每队列都训练）===\n")
build_one <- function(T) {
  dd <- CL[[T]]; z <- dd$z; tm <- dd$time; st <- dd$status
  out <- list()
  sc <- try(cox_score(z, tm, st), silent = TRUE)
  if (inherits(sc, "try-error")) return(data.table(train = T, err = "score_failed"))
  g0 <- names(sort(abs(sc), decreasing = TRUE))[seq_len(min(PRE_N, length(sc)))]
  rm(sc); invisible(gc(FALSE))
  x <- t(z[g0, , drop = FALSE]); rm(g0); invisible(gc(FALSE))
  # GSE39582 里有 6 例随访时间 <= 0（其中 4 例是事件，推测为围手术期死亡）。
  # glmnet 的 Cox 族不接受非正生存时间，会直接报
  #   "Non-positive event times encountered" —— 这一条最初导致整个
  #   GSE39582（579 例，最大的队列）被静默排除，天平因此被压歪。
  # 处理：把非正时间压到最小正值的一半（该队列最小正值为 1 个月，
  #       取 0.5 可保持事件顺序且严格为正）；验证端仍用原始时间。
  tm_tr <- tm
  if (any(tm_tr <= 0)) {
    minpos <- min(tm_tr[tm_tr > 0])
    tm_tr[tm_tr <= 0] <- minpos / 2
    message("    [", T, "] 修正 ", sum(tm <= 0), " 例非正随访时间 -> ",
            minpos / 2, " 个月")
  }
  fit <- try(glmnet(x, Surv(tm_tr, st), family = "cox", alpha = 1, nlambda = 200,
                    standardize = TRUE, maxit = 5000), silent = TRUE)
  rm(x); invisible(gc(FALSE))
  if (inherits(fit, "try-error"))
    return(data.table(train = T, err = sub("\n.*", "", conditionMessage(attr(fit, "condition")))))
  beta <- as.matrix(fit$beta); nz <- colSums(beta != 0)
  for (K in K_TARGETS) {
    cand <- which(nz == K)
    if (!length(cand)) { d <- abs(nz - K); cand <- which(d == min(d)) }
    lam <- cand[which.max(fit$lambda[cand])]
    gs <- rownames(beta)[which(beta[, lam] != 0)]
    if (length(gs) < 2) next
    for (V in cohorts) {
      if (V == T) next
      out[[length(out) + 1]] <- data.table(
        train = T, test = V, K_target = K, n_gene = length(gs),
        genes = paste(gs, collapse = ";"),
        C = oriented_C(CL[[V]]$z, gs, CL[[V]]$time, CL[[V]]$status))
    }
  }
  if (!length(out)) return(data.table(train = T, err = "no_model"))
  rbindlist(out)
}
t0 <- Sys.time()
parts <- lapply(cohorts, build_one)   # 串行：mclapply 偶发死锁，第 113 个任务曾卡住不返回
res <- list()
for (i in seq_along(parts)) {
  p <- parts[[i]]
  if (!is.null(p) && "err" %in% names(p) && nrow(p) == 1 && !is.na(p$err[1]) &&
      p$err[1] != "") {
    cat("  [失败] ", cohorts[i], " -> ", p$err[1], "\n", sep = "")
  } else if (!is.null(p)) {
    res[[length(res) + 1]] <- p
    cat(sprintf("  %-24s 产出 %d 条\n", cohorts[i], nrow(p)))
  }
}
dnm <- rbindlist(res)
dnm <- dnm[!is.na(C)]
cat("  合计 ", nrow(dnm), " 条；参与训练的队列 ", uniqueN(dnm$train), " / ",
    length(cohorts), "\n", sep = "")
cat("  用时 ", round(as.numeric(Sys.time() - t0, units = "mins"), 1), " 分钟\n", sep = "")

# ------------------------------------------- B 为每个(队列,基因数)算配对随机基线
cat("\n=== B. 配对随机基线（B = ", B, "）===\n", sep = "")
need <- unique(dnm[, .(cohort = test, n_gene)])
cat("  需要 ", nrow(need), " 个（队列 × 基因数）组合\n", sep = "")
null_one <- function(i) {
  ch <- need$cohort[i]; ng <- need$n_gene[i]
  dd <- CL[[ch]]; z <- dd$z; tm <- dd$time; st <- dd$status
  pool <- rownames(z)
  if (length(pool) < ng * 2)
    return(data.table(cohort = ch, n_gene = ng, null_median = NA_real_,
                      null_mean = NA_real_, null_p95 = NA_real_))
  n <- length(tm); sv <- Surv(tm, st); xb <- matrix(0, n, 1)
  set.seed(990101L + i)
  v <- numeric(B)
  nk <- seq_len(nrow(z))
  # 提速要点：用【整数下标】而不是基因名取行，避开逐次的行名匹配。
  # 千万不要为了去行名而写 rownames(zm) <- NULL —— 那会整份拷贝矩阵
  # （最大队列约 100MB），循环上百次后 GC 抖到像是死锁。
  for (b in seq_len(B)) {
    xb[, 1] <- colMeans(z[sample(nk, ng), , drop = FALSE])
    v[b] <- cz(sv, xb)
  }
  cat(sprintf("    [%s|n=%d] %d/%d 完成\n", ch, ng, i, nrow(need)))
  data.table(cohort = ch, n_gene = ng,
             null_median = median(v, na.rm = TRUE),
             null_mean = mean(v, na.rm = TRUE),
             null_p95 = as.numeric(quantile(v, .95, na.rm = TRUE)))
}
t1 <- Sys.time()
nl <- rbindlist(lapply(seq_len(nrow(need)), null_one))
cat("  完成 ", nrow(nl), " 个组合，用时 ",
    round(as.numeric(Sys.time() - t1, units = "mins"), 1), " 分钟\n", sep = "")
fwrite(nl, file.path(OUT, "null_matched.csv"))

# ---------------------------------------------------------------- C 合并
dnm <- merge(dnm, nl, by.x = c("test", "n_gene"), by.y = c("cohort", "n_gene"),
             all.x = TRUE, sort = FALSE)
dnm[, dC := C - null_median]
dnm[, beat_random := as.integer(C > null_median)]
dnm[, exceed_p95 := as.integer(C > null_p95)]
fwrite(dnm, file.path(OUT, "de_novo_full.csv"))

cat("\n=== 汇总 ===\n")
cat("自产签名 ", nrow(dnm), " 条，随机基线缺失 ", sum(is.na(dnm$null_median)), " 条\n", sep = "")
cat("ΔC 中位 ", round(dnm[, median(dC, na.rm = TRUE)], 4),
    "  优于随机 ", dnm[, sum(beat_random == 1, na.rm = TRUE)], "/", nrow(dnm), "\n", sep = "")
cat("\n每个训练队列的表现：\n")
print(dnm[, .(记录 = .N, ΔC中位 = round(median(dC, na.rm = TRUE), 4),
              C中位 = round(median(C, na.rm = TRUE), 4),
              优于随机 = paste0(sum(beat_random == 1, na.rm = TRUE), "/", .N)),
          by = train][order(-ΔC中位)])
cat("\n训练队列规模与其产物质量（关键：验证「训练集越大越好」这一假设）\n")
info <- data.table(cohort = cohorts,
                   n = sapply(CL, function(x) length(x$time)),
                   events = sapply(CL, function(x) sum(x$status)))
print(merge(info, dnm[, .(ΔC中位 = round(median(dC, na.rm = TRUE), 4)), by = train],
            by.x = "cohort", by.y = "train", all.x = TRUE)[order(-n)])
