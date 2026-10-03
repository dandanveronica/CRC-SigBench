#!/usr/bin/env Rscript
# =====================================================================
# 07b_highres_null.R  (v2)
# 第二阶段：高分辨率零分布 —— 解决 B=1000 时经验 p 分辨率不足的问题
#
# 设计要点
#   1. 零分布只依赖 (队列, 实际参与打分的基因数)，与具体签名无关。
#      per_cohort.csv 里只有 138 个这样的组合 —— 所以【一次性把 138 个组合都
#      跑到高 B】，全部 320 个签名即同时获得高分辨率 p 值，
#      不需要"先筛候选再加密"的二阶段设计（后者会带来选择性推断偏倚）。
#   2. 经验 p 的下限 = 1/(B+1)。B=100,000 -> p_min ≈ 1.0e-5
#      全族检验数 1920，BH 后最小可达 = 1e-5 × 1920 = 0.019 < 0.05
#      => 这一次「有没有签名通过 FDR」才真正是有信息量的结论。
#   3. 两处提速（已验证，结果与 v1 完全一致）：
#      (a) 实测 C(-p) = 1 - C(p) 精确成立（多个基因数下和均为 1.000000），
#          所以 oriented C = max(C, 1-C)，不必调用两次 concordance；
#      (b) concordancefit(y, x, std.err=FALSE) 比 concordance(formula) 快 2.18×，
#          数值结果完全相同（差 = 0）。
#      合计提速 4.36×。
#   4. 按队列分批并行：同一批 worker 只碰同一个表达矩阵，显著缓解多进程
#      抢内存带宽导致的抖动（v1 实测并行开销高达 6 倍）。
#
# 用法：Rscript 07b_highres_null.R [B] [cores]
# 输出：audit_out_clean/highres_null.rds
#       audit_out_clean/highres_per_cohort.csv
#       audit_out_clean/highres_summary.csv
# =====================================================================

suppressPackageStartupMessages({
  library(survival); library(data.table); library(parallel)
})

args  <- commandArgs(trailingOnly = TRUE)
B     <- if (length(args) >= 1) as.integer(args[1]) else 100000L
CORES <- if (length(args) >= 2) as.integer(args[2]) else 7L
COH   <- "cohorts"
OUT   <- "audit_out_clean"
STAGE1<- file.path(OUT, "per_cohort.csv")
CACHE <- file.path(OUT, "highres_null.rds")

cat("=== 第二阶段：高分辨率零分布 v2 ===\n")
cat("B =", B, "  cores =", CORES, "\n")
cat("经验 p 分辨率 = 1/(B+1) =", sprintf("%.2e", 1 / (B + 1)), "\n\n")

zscore <- function(m) {
  mu <- rowMeans(m, na.rm = TRUE); sd <- apply(m, 1, sd, na.rm = TRUE)
  sd[is.na(sd) | sd < 1e-8] <- 1; (m - mu) / sd
}
# 定向 C-index：利用 C(-p) = 1 - C(p) 的恒等式，只算一次
cindex1 <- function(sv, xbuf) {
  cf <- try(concordancefit(sv, xbuf, reverse = TRUE, std.err = FALSE),
            silent = TRUE)
  if (inherits(cf, "try-error")) return(NA_real_)
  c0 <- cf$concordance
  if (is.na(c0)) return(NA_real_)
  if (c0 < 0.5) 1 - c0 else c0          # = max(C, 1-C)
}

# ---------------------------------------------------- 载入 stage1 结果
res0 <- fread(STAGE1); res0 <- res0[status_r == "ok"]
setnames(res0, "C", "C_stage1")
combos <- unique(res0[, .(cohort, n_measured)])[order(cohort, n_measured)]
cat("待算零分布组合数：", nrow(combos), "（队列 × 实测基因数）\n\n", sep = "")

# 断点续跑：已算过的组合直接复用
null_list <- if (file.exists(CACHE)) readRDS(CACHE) else list()
have <- function(ch, ng) {
  v <- null_list[[paste0(ch, "|", ng)]]
  !is.null(v) && length(v) >= B * 0.95
}
todo <- combos[!mapply(have, combos$cohort, combos$n_measured)]
cat("已完成 ", nrow(combos) - nrow(todo), " / ", nrow(combos),
    "，本次待算 ", nrow(todo), "\n\n", sep = "")

run_one <- function(ch, ng) {
  d   <- readRDS(file.path(COH, paste0(ch, ".rds")))
  z   <- zscore(d$expr)
  cl  <- d$clin; cl <- cl[match(colnames(z), cl$sample), , drop = FALSE]
  tm  <- suppressWarnings(as.numeric(cl$time))
  st  <- suppressWarnings(as.integer(cl$status))
  kp  <- !is.na(tm) & !is.na(st)
  z   <- z[, kp, drop = FALSE]; tm <- tm[kp]; st <- st[kp]
  rm(d, cl); invisible(gc(FALSE))

  pool <- rownames(z)
  n    <- length(tm)
  sv   <- Surv(tm, st)
  xbuf <- matrix(0, nrow = n, ncol = 1)
  set.seed(as.integer(20261001L + sum(utf8ToInt(ch)) * 1000L + ng))
  vals <- numeric(B)
  # 关键提速：抽样时用【整数下标】而非基因名，避开逐次的行名匹配，
  # 实测相差一个数量级。注意不要写 rownames(zm) <- NULL 去行名，
  # 那会整份拷贝矩阵（最大队列约 100MB），反复拷贝会把 GC 拖到像死锁。
  nk   <- seq_len(nrow(z))
  for (b in seq_len(B)) {
    xbuf[, 1] <- colMeans(z[sample(nk, ng), , drop = FALSE])
    vals[b] <- cindex1(sv, xbuf)
  }
  vals <- vals[!is.na(vals)]
  list(key = paste0(ch, "|", ng), vals = vals)
}

t00 <- Sys.time()
if (nrow(todo)) {
  by_coh <- split(todo, todo$cohort)
  for (ch in names(by_coh)) {
    sub <- by_coh[[ch]]
    cat(sprintf("[%s] %2d 个组合（基因数 %s）...\n", ch, nrow(sub),
                paste(sort(unique(sub$n_measured)), collapse = ",")))
    t0 <- Sys.time()
    out <- mclapply(seq_len(nrow(sub)), function(i) run_one(ch, sub$n_measured[i]),
                    mc.cores = min(CORES, nrow(sub)))
    for (o in out) {
      if (!is.null(o) && length(o$vals)) null_list[[o$key]] <- o$vals
    }
    saveRDS(null_list, CACHE)      # 增量落盘，崩溃可续跑
    cat(sprintf("    完成，用时 %.1f 分钟，累计 %.1f 分钟\n",
                as.numeric(Sys.time() - t0, units = "mins"),
                as.numeric(Sys.time() - t00, units = "mins")), flush = TRUE)
  }
}
cat("\n总耗时：", round(as.numeric(difftime(Sys.time(), t00, units = "mins")), 1), "分钟\n", sep = "")
cat("零分布组合总数：", length(null_list), " / ", nrow(combos), "\n\n", sep = "")

# ---------------------------------------------------- 重算经验 p
emp_p <- function(real, nv) {
  if (is.na(real) || is.null(nv) || length(nv) == 0) return(NA_real_)
  (sum(nv >= real) + 1) / (length(nv) + 1)
}
res0[, key := paste0(cohort, "|", n_measured)]
res0[, p_highres     := vapply(seq_len(.N), function(i) emp_p(C_stage1[i], null_list[[key[i]]]), numeric(1))]
res0[, null_median_hi:= vapply(seq_len(.N), function(i) { v <- null_list[[key[i]]]
                              if (is.null(v)) NA_real_ else median(v) }, numeric(1))]
res0[, null_mean_hi  := vapply(seq_len(.N), function(i) { v <- null_list[[key[i]]]
                              if (is.null(v)) NA_real_ else mean(v) }, numeric(1))]
res0[, null_sd_hi    := vapply(seq_len(.N), function(i) { v <- null_list[[key[i]]]
                              if (is.null(v)) NA_real_ else sd(v) }, numeric(1))]
res0[, z_vs_null     := (C_stage1 - null_mean_hi) / null_sd_hi]
# 全族校正：320 签名 × 6 队列 = 1920 条，作为一个 family
res0[, p_adj_highres := p.adjust(p_highres, method = "BH")]
fwrite(res0, file.path(OUT, "highres_per_cohort.csv"))

# ---------------------------------------------------- 签名级汇总
summ <- res0[, .(
  n_cohort    = .N,
  C_mean      = round(mean(C_stage1, na.rm = TRUE), 4),
  null_mean   = round(mean(null_mean_hi, na.rm = TRUE), 4),
  deltaC_mean = round(mean(C_stage1 - null_mean_hi, na.rm = TRUE), 4),
  z_mean      = round(mean(z_vs_null, na.rm = TRUE), 3),
  p_min       = round(min(p_highres, na.rm = TRUE), 6),
  p_adj_min   = round(min(p_adj_highres, na.rm = TRUE), 6),
  n_p05       = sum(p_highres < 0.05, na.rm = TRUE),
  n_fdr05     = sum(p_adj_highres < 0.05, na.rm = TRUE),
  n_coh_beats = sum(C_stage1 > null_median_hi, na.rm = TRUE)
), by = sig_id]
summ <- merge(summ, unique(res0[, .(sig_id, n_gene)]), by = "sig_id", all.x = TRUE)
setorder(summ, -n_fdr05, -deltaC_mean)
fwrite(summ, file.path(OUT, "highres_summary.csv"))

cat("=== 完成 ===\n")
cat("零分布：", length(null_list), " 个组合，每个有效样本数中位 = ",
    round(median(vapply(null_list, length, integer(1)))), "\n", sep = "")
cat("经验 p 最小可达：", sprintf("%.2e", 1 / (B + 1)),
    "   实测最小 p = ", sprintf("%.2e", min(res0$p_highres, na.rm = TRUE)), "\n", sep = "")
cat("全族检验数：", nrow(res0), "   BH 后最小可达 = ",
    sprintf("%.3f", 1 / (B + 1) * nrow(res0)), "\n", sep = "")
cat("原始 p<0.05：", sum(res0$p_highres < 0.05, na.rm = TRUE), " / ", nrow(res0),
    "（", round(100 * mean(res0$p_highres < 0.05, na.rm = TRUE), 1), "%）\n", sep = "")
cat("FDR<0.05 的签名×队列组合：", sum(res0$p_adj_highres < 0.05, na.rm = TRUE), "\n", sep = "")
cat("至少 1 个队列通过 FDR 的签名：", sum(summ$n_fdr05 >= 1), "\n", sep = "")
cat("至少 3 个队列通过 FDR 的签名：", sum(summ$n_fdr05 >= 3), "\n", sep = "")
