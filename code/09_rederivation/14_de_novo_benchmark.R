#!/usr/bin/env Rscript
# =====================================================================
# 14_de_novo_benchmark.R
# 核心新增分析：四层对照 —— 把结论从「排名」翻转成「天花板」
#
# 【为什么必须做这个】
#   Part 1 用随机基因集做零假设，但这个对照太弱：随机集是从两万个基因里
#   纯随机抽的，绝大多数与 CRC 预后毫无关系；而已发表签名是在别的数据集上
#   单因素 Cox 筛过的。真实签名赢过随机集是「应得的」，审稿人一句话就能
#   顶回来。真正公平的对手不是随机基因，而是：
#       【如果今天有人用同样的方法、同样的工作量，现场重做一个签名，
#         能做到什么水平？】
#
# 【四层证据】
#   L0 队列内 split-half      —— 同分布、无批次效应下的「技术天花板」
#   L1 自产签名 leave-one-out  —— 同等方法学努力的「跨站点可迁移上限」
#   L2 已发表签名外部验证      —— 真实世界表现（来自 07 per_cohort.csv）
#   L3 随机等大基因集          —— 零假设（来自 07 null_dist）
#
#   若 L1 ≈ L2  → 不是哪篇论文不行，是这个任务的跨队列信息上限就在这里
#   若 L1 > L2  → 已发表签名连标准流程都做不到位（更负面）
#   若 L2 > L1  → 已发表签名确有独到之处
#   无论哪种结果都是此前没人给过答案的新信息。
#
# 【保证对称的关键处理】
#   自产签名只借用 LASSO「选基因」的能力，打分环节强制使用与已发表签名
#   完全相同的做法：队列内 z-score 后等权重求和 + 方向定向。
#   否则天平又会歪向自产签名（它白拿了 LASSO 训练出来的系数）。
#
# 用法：Rscript 14_de_novo_benchmark.R [cores] [N_SPLIT]
# 输出：audit_out_clean/de_novo_自产签名.csv
#       audit_out_clean/de_novo_队列内上限.csv
#       audit_out_clean/de_novo_分层对比.csv
# =====================================================================

suppressPackageStartupMessages({
  library(survival); library(data.table); library(glmnet); library(parallel)
})

args    <- commandArgs(trailingOnly = TRUE)
CORES   <- if (length(args) >= 1) as.integer(args[1]) else 2L
N_SPLIT <- if (length(args) >= 2) as.integer(args[2]) else 50L
COH   <- "cohorts"; OUT <- "audit_out_clean"
K_TARGETS <- c(3, 5, 8, 10, 15, 20, 30, 50)
PRE_N     <- 1000L          # 预筛保留基因数
K_SPLIT   <- 10L            # 队列内 split-half 分析中固定的基因数

zscore <- function(m) {
  mu <- rowMeans(m, na.rm = TRUE); sd <- apply(m, 1, sd, na.rm = TRUE)
  sd[is.na(sd) | sd < 1e-8] <- 1; (m - mu) / sd
}
cz <- function(sv, xbuf) {                       # 定向 C（C(-p)=1-C(p)）
  cf <- try(concordancefit(sv, xbuf, reverse = TRUE, std.err = FALSE), silent = TRUE)
  if (inherits(cf, "try-error")) return(NA_real_)
  c0 <- cf$concordance
  if (is.na(c0)) return(NA_real_)
  if (c0 < 0.5) 1 - c0 else c0
}
oriented_C <- function(zz, genes, time, status) {
  genes <- intersect(genes, rownames(zz))
  if (length(genes) < 2) return(NA_real_)
  pred <- as.numeric(colMeans(zz[genes, , drop = FALSE]))
  if (length(unique(pred)) < 2) return(NA_real_)
  cz(Surv(time, status), matrix(pred, ncol = 1))
}

# ---------------------------------------------------------------
# 向量化单因素 Cox score 检验（对每个基因一次算出 score test z 值）
# 逐基因跑 coxph 会慢到不可行（2 万基因 × 多队列），这里用累积和
# 把「每个事件时间点的风险集均值」一次性算出来，复杂度 O(基因 × 样本)
# ---------------------------------------------------------------
cox_score <- function(z, time, status) {
  ord <- order(time, decreasing = TRUE)          # 降序：位置 <= j 即 t >= t_j
  X <- z[, ord, drop = FALSE]
  t_o <- time[ord]; s_o <- status[ord]
  r <- rle(t_o)
  ends <- cumsum(r$lengths)                      # 每个时间块的【末位置】
  ne <- vapply(seq_along(ends), function(k) {
    ii <- if (k == 1L) 1:ends[1] else (ends[k-1] + 1):ends[k]
    sum(s_o[ii] == 1)
  }, numeric(1))
  keep <- ne > 0
  if (!any(keep)) return(setNames(rep(0, nrow(z)), rownames(z)))
  ends <- ends[keep]; ne <- ne[keep]

  CX  <- t(apply(X, 1, cumsum))                  # 基因 × 样本 累积和
  CX2 <- t(apply(X * X, 1, cumsum))
  U <- numeric(nrow(z)); V <- numeric(nrow(z))
  # 事件块的 X 行和：同样用累积和差分得到
  # 注意：不能写成 X * (s_o == 1)。R 里「矩阵 × 向量」是按列优先存储顺序
  # 循环取值做逐元素乘，不是按列缩放 —— 那样会把行和列的对应关系彻底打乱
  # （实测 U 与暴力法相关性掉到 -0.11，而 V 不受影响，极具迷惑性）。
  # 正解：先把非事件的样本列整体置零。
  XS <- X; XS[, s_o != 1] <- 0
  SXblk <- t(apply(XS, 1, cumsum))
  for (i in seq_along(ends)) {
    eu <- ends[i]
    S1 <- CX[, eu] / eu
    Se <- SXblk[, eu] - if (i == 1L) 0 else SXblk[, ends[i-1]]
    S2 <- CX2[, eu] / eu
    U <- U + Se - ne[i] * S1
    V <- V + ne[i] * pmax(S2 - S1 * S1, 1e-8)
  }
  zval <- ifelse(V > 1e-12, U / sqrt(V), 0)
  zval[is.na(zval)] <- 0
  setNames(zval, rownames(z))
}

cat("=== 载入队列 ===\n")
CL <- list()
for (cf in list.files(COH, pattern = "\\.rds$", full.names = TRUE)) {
  ch <- sub("\\.rds$", "", basename(cf))
  d  <- readRDS(cf)
  z  <- zscore(d$expr)
  cl <- d$clin; cl <- cl[match(colnames(z), cl$sample), , drop = FALSE]
  tm <- suppressWarnings(as.numeric(cl$time))
  st <- suppressWarnings(as.integer(cl$status))
  kp <- !is.na(tm) & !is.na(st)
  CL[[ch]] <- list(z = z[, kp, drop = FALSE], time = tm[kp], status = st[kp])
  cat(sprintf("  %-24s %4d 样本 %3d 事件 %d 基因\n", ch, sum(kp), sum(st[kp]), nrow(z)))
  rm(d); invisible(gc(FALSE))
}
cohorts <- names(CL)

# ---------------------------------------------------------------
# A. 自产签名：在队列 T 构建 -> 在其余队列外部验证
#    只借用 LASSO 的「选基因」能力，打分仍用与已发表签名一致的等权重 z
# ---------------------------------------------------------------
cat("\n=== A. 自产签名 leave-one-cohort-out ===\n")
build_and_eval <- function(T, K_TARGETS, PRE_N) {
  dd <- CL[[T]]; z <- dd$z; tm <- dd$time; st <- dd$status
  sc <- cox_score(z, tm, st)
  genes0 <- names(sort(abs(sc), decreasing = TRUE))[seq_len(min(PRE_N, length(sc)))]
  x <- t(z[genes0, , drop = FALSE])              # 样本 × 基因
  fit <- try(glmnet(x, Surv(tm, st), family = "cox", alpha = 1, nlambda = 200,
                    standardize = TRUE, maxit = 5000), silent = TRUE)
  if (inherits(fit, "try-error")) return(NULL)
  beta <- as.matrix(fit$beta)
  nz   <- colSums(beta != 0)
  out  <- list()
  for (K in K_TARGETS) {
    cand <- which(nz == K)
    if (length(cand) == 0) {
      # 没有恰好 K 个的 lambda，退而取最接近的
      df <- abs(nz - K); cand <- which(df == min(df))
    }
    lam <- cand[which.max(fit$lambda[cand])]     # 同样 K 里取惩罚更重的
    gs  <- rownames(beta)[which(beta[, lam] != 0)]
    if (length(gs) < 2) next
    for (V in cohorts) {
      if (V == T) next
      out[[length(out) + 1]] <- data.table(
        source = "de_novo", train = T, test = V, K = K, n_gene = length(gs),
        genes = paste(gs, collapse = ";"),
        C = oriented_C(CL[[V]]$z, gs, CL[[V]]$time, CL[[V]]$status))
    }
  }
  rbindlist(out)
}
resA <- rbindlist(mclapply(cohorts, build_and_eval, K_TARGETS = K_TARGETS,
                           PRE_N = PRE_N, mc.cores = max(1L, CORES)))
resA <- resA[!is.na(C)]
cat("  自产签名记录：", nrow(resA), " 条\n", sep = "")
fwrite(resA, file.path(OUT, "de_novo_自产签名.csv"))

# ---------------------------------------------------------------
# B. 队列内 split-half：同分布、无批次效应的「技术天花板」
#    这部分回答「这个数据源本身的信噪比能撑到多少」
# ---------------------------------------------------------------
cat("\n=== B. 队列内随机 50/50 split（技术天花板）===\n")
split_once <- function(ch, rep_i) {
  dd <- CL[[ch]]; z <- dd$z; tm <- dd$time; st <- dd$status
  n <- length(tm)
  set.seed(20261001L + rep_i)
  tr <- sample(seq_len(n), floor(n / 2))
  te <- setdiff(seq_len(n), tr)
  if (sum(st[tr]) < 10 || sum(st[te]) < 5) return(NULL)
  sc <- cox_score(z[, tr, drop = FALSE], tm[tr], st[tr])
  genes0 <- names(sort(abs(sc), decreasing = TRUE))[seq_len(min(PRE_N, length(sc)))]
  x <- t(z[genes0, tr, drop = FALSE])
  fit <- try(glmnet(x, Surv(tm[tr], st[tr]), family = "cox", alpha = 1,
                    nlambda = 200, standardize = TRUE, maxit = 5000), silent = TRUE)
  if (inherits(fit, "try-error")) return(NULL)
  beta <- as.matrix(fit$beta); nz <- colSums(beta != 0)
  cand <- which(nz == K_SPLIT)
  if (!length(cand)) { df <- abs(nz - K_SPLIT); cand <- which(df == min(df)) }
  lam <- cand[which.max(fit$lambda[cand])]
  gs  <- rownames(beta)[which(beta[, lam] != 0)]
  if (length(gs) < 2) return(NULL)
  data.table(cohort = ch, rep = rep_i, K_actual = length(gs),
             C_val = oriented_C(z[, te, drop = FALSE], gs, tm[te], st[te]),
             n_train = length(tr), event_train = sum(st[tr]), n_test = length(te))
}
resB <- rbindlist(lapply(cohorts, function(ch)
  rbindlist(mclapply(seq_len(N_SPLIT), function(i) split_once(ch, i),
                     mc.cores = max(1L, CORES)))))
resB <- resB[!is.na(C_val)]
fwrite(resB, file.path(OUT, "de_novo_队列内上限.csv"))
aggB <- resB[, .(n = .N, C中位 = round(median(C_val, na.rm = TRUE), 3),
                 C四分位 = paste0(round(quantile(C_val, .25, na.rm = TRUE), 3), "~",
                                  round(quantile(C_val, .75, na.rm = TRUE), 3))), by = cohort]
cat("  各队列天花板（队列内 split-half）：\n")
print(aggB)
cat("  合并 ", nrow(resB), " 次重复，总体中位 C = ",
    round(resB[, median(C_val, na.rm = TRUE)], 3), "\n", sep = "")

# ---------------------------------------------------------------
# C. 分层对比：同一基因数 K 下，已发表 vs 自产 vs 随机
# ---------------------------------------------------------------
cat("\n=== C. 按基因数分层：已发表 vs 自产 vs 随机 ===\n")
per <- fread(file.path(OUT, "per_cohort.csv"))
per <- per[status_r == "ok"]
per[, dC := C - null_median]

# 把自产签名的 K 对齐到已发表签名的 n_gene 分箱
bins <- c(2, 3, 5, 8, 10, 15, 20, 30, 40, 100)
binof <- function(v) cut(v, bins, right = FALSE, labels = FALSE)
per[, Kbin := binof(n_measured)]
resA[, Kbin := binof(n_gene)]
labK <- sapply(seq_len(length(bins) - 1), function(i) paste0(bins[i], "~", bins[i+1] - 1))

cmp <- data.table()
for (kb in unique(c(per$Kbin, resA$Kbin))) {
  if (is.na(kb)) next
  a <- per[Kbin == kb, C]; b <- resA[Kbin == kb, C]; r <- per[Kbin == kb, null_median]
  if (length(a) < 5 || length(b) < 3) next
  pv <- if (length(a) >= 5 && length(b) >= 3)
    suppressWarnings(wilcox.test(a, b, alternative = "two.sided")$p.value) else NA_real_
  cmp <- rbind(cmp, data.table(
    基因数 = labK[kb],
    已发表_n = length(a), 已发表C中位 = round(median(a, na.rm = TRUE), 3),
    自产_n = length(b),   自产C中位   = round(median(b, na.rm = TRUE), 3),
    随机C中位 = round(median(r, na.rm = TRUE), 3),
    自产减已发表 = round(median(b, na.rm = TRUE) - median(a, na.rm = TRUE), 3),
    Wilcoxon_p = if (is.na(pv)) NA else signif(pv, 3)))
}
print(cmp)
fwrite(cmp, file.path(OUT, "de_novo_分层对比.csv"))

cat("\n=== D. 四层对照总览 ===\n")
allC <- resA[, median(C, na.rm = TRUE)]
pubC <- per[, median(C, na.rm = TRUE)]
rndC <- per[, median(null_median, na.rm = TRUE)]
inC  <- resB[, median(C_val, na.rm = TRUE)]
tab <- data.table(
  层级 = c("L0 队列内 split-half（同分布）", "L1 自产签名·跨队列",
           "L2 已发表签名·跨队列", "L3 随机等大基因集"),
  说明 = c("同一队列随机分半，无批次效应的技术天花板",
           "在队列 A 现场构建，丢到其余队列：同等方法学努力的可迁移上限",
           "原文在其自报训练集之外的独立队列上的表现",
           "纯零假设基线"),
  C中位数 = round(c(inC, allC, pubC, rndC), 4))
print(tab, nrows = 10)
fwrite(tab, file.path(OUT, "de_novo_四层对照.csv"))

cat("\n结论提示：\n")
d_pub_vs_de <- pubC - allC
if (abs(d_pub_vs_de) < 0.01) {
  cat("  已发表签名与自产签名几乎持平（差 ", round(d_pub_vs_de, 4),
      "）→ 支持【天花板】叙事：不是方法不行，是任务本身的跨队列信息上限如此。\n", sep = "")
} else if (d_pub_vs_de < -0.01) {
  cat("  自产签名反超已发表签名 ", round(-d_pub_vs_de, 4),
      " → 支持【已发表签名连标准流程都不如】的更负面结论。\n", sep = "")
} else {
  cat("  已发表签名高于自产签名 ", round(d_pub_vs_de, 4),
      " → 已发表签名确有独到之处，可直接作为正面结论。\n", sep = "")
}
cat("  队列内天花板 ", round(inC, 3), " vs 跨队列 ", round(allC, 3),
    "：同分布内的可迁移性远高于跨队列，落差 ",
    round(inC - allC, 3), " 即为批次/人群差异的真实代价。\n", sep = "")
