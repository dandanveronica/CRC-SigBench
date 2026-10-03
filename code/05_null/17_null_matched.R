#!/usr/bin/env Rscript
# =====================================================================
# 17_null_matched.R —— 只为「自产签名」重算配对的随机基线
#
# 为什么单独拆出来：
#   14b 里从内存列表 CL[[cohort]] 反复取大矩阵时，进程会在第 ~104 个组合
#   处明显卡住；而逐组合从磁盘 readRDS 再算，实测每个只要 1-2 秒。
#   拆开后 113 个组合并行约 1 分钟，且失败可单独重跑。
#
# 输入：audit_out_clean/de_novo_full.csv（只需 test 与 n_gene 两列）
# 输出：audit_out_clean/null_matched.csv
# =====================================================================
suppressPackageStartupMessages({
  library(survival); library(data.table); library(parallel)
})
args  <- commandArgs(trailingOnly = TRUE)
B     <- if (length(args) >= 1) as.integer(args[1]) else 5000L
CORES <- if (length(args) >= 2) as.integer(args[2]) else 4L
COH <- "cohorts"; OUT <- "audit_out_clean"

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

dnm <- fread(file.path(OUT, "de_novo_full.csv"))
need <- unique(dnm[, .(cohort = test, n_gene)])
cat("需要 ", nrow(need), " 个（验证队列 × 基因数）组合，B = ", B, "\n", sep = "")

null_one <- function(i) {
  ch <- need$cohort[i]; ng <- as.integer(need$n_gene[i])
  f <- file.path(COH, paste0(ch, ".rds"))
  d  <- readRDS(f)                       # 只读一次
  cl <- d$clin
  z  <- zscore(d$expr); rm(d)
  cl <- cl[match(colnames(z), cl$sample), , drop = FALSE]
  tm <- suppressWarnings(as.numeric(cl$time))
  st <- suppressWarnings(as.integer(cl$status))
  kp <- !is.na(tm) & !is.na(st)
  z <- z[, kp, drop = FALSE]; tm <- tm[kp]; st <- st[kp]
  nk <- seq_len(nrow(z)); n <- length(tm)
  sv <- Surv(tm, st); xb <- matrix(0, n, 1)
  set.seed(990101L + i)
  v <- numeric(B)
  for (b in seq_len(B)) {
    xb[, 1] <- colMeans(z[sample(nk, ng), , drop = FALSE])
    v[b] <- cz(sv, xb)
  }
  cat(sprintf("  [%s|n=%d] %d/%d 完成\n", ch, ng, i, nrow(need)))
  data.table(cohort = ch, n_gene = ng,
             null_median = median(v, na.rm = TRUE),
             null_mean = mean(v, na.rm = TRUE),
             null_p95 = as.numeric(quantile(v, .95, na.rm = TRUE)))
}

t0 <- Sys.time()
nl <- rbindlist(mclapply(seq_len(nrow(need)), null_one, mc.cores = CORES))
cat("完成 ", nrow(nl), " 个组合，用时 ",
    round(as.numeric(Sys.time() - t0, units = "mins"), 1), " 分钟\n", sep = "")
fwrite(nl, file.path(OUT, "null_matched.csv"))
cat("已写出 null_matched.csv\n")
