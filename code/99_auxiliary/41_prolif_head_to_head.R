#!/usr/bin/env Rscript
# ============================================================
# 41_prolif_head_to_head.R —— 增殖评分 P 单变量 vs 312 签名（新分析 ① 变体）
#
# 目的：量化 "312 个签名 collectively 打不过一个增殖评分"
#   每队列：C(P) oriented C-index（P = KEGG hsa04110 124 基因 z 均值，
#           与 39_proliferation_adjust.R 完全同口径）
#           + 对 size≈124 随机基因集 null 的增量（null 中位取自 S1 released，
#             n_measured∈[118,130] 行的中位，零新模拟）
#           + P 在该队列打赢多少 % 的签名（S1 C 分布）
# 附带：dump 每队列实测基因行名（union）→ universe 文件（供 42 ORA 用）
#
# 口径函数逐字复制 39_proliferation_adjust.R（其复制自 07_audit_signatures.R）
# 输出：audit_out_clean/pathway/
#   prolif_head_to_head.csv
#   measured_universe.csv   （gene, n_cohorts_measured）
# ============================================================
suppressPackageStartupMessages({
  library(survival); library(data.table)
})

BASE <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)))
if (!length(BASE)) BASE <- getwd()
ROOT <- normalizePath(file.path(BASE, "..", ".."))
COH  <- file.path(ROOT, "cohorts")
OUTD <- file.path(ROOT, "audit_out_clean", "pathway")
dir.create(OUTD, showWarnings = FALSE, recursive = TRUE)

# ---------------- 口径函数（39 原文复制） ----------------
zscore <- function(m) {
  mu <- rowMeans(m, na.rm = TRUE)
  sd <- apply(m, 1, sd, na.rm = TRUE)
  sd[is.na(sd) | sd < 1e-8] <- 1
  (m - mu) / sd
}
cindex <- function(time, status, pred) {
  ok <- !is.na(pred) & !is.na(time) & !is.na(status)
  if (sum(ok) < 20 || sum(status[ok]) < 5) return(NA_real_)
  if (length(unique(pred[ok])) < 2) return(NA_real_)
  cf <- try(survival::concordance(Surv(time[ok], status[ok]) ~ pred[ok],
                                  reverse = TRUE), silent = TRUE)
  if (inherits(cf, "try-error")) return(NA_real_)
  unname(cf$concordance)
}
oriented_cindex <- function(time, status, pred) {
  c0 <- cindex(time, status, pred)
  if (is.na(c0)) return(c(C = NA_real_, flipped = NA))
  c1 <- cindex(time, status, -pred)
  if (is.na(c1)) return(c(C = c0, flipped = 0))
  if (c1 > c0) c(C = c1, flipped = 1) else c(C = c0, flipped = 1 * 0 + 0)
}

# ---------------- 数据 ----------------
KEGG_FILE <- file.path(ROOT, "data", "proliferation", "kegg_hsa04110_symbols.txt")
kegg <- trimws(readLines(KEGG_FILE))
kegg <- kegg[nchar(kegg) > 0 & !startsWith(kegg, "#")]

s1 <- fread(file.path(ROOT, "results", "per_cohort_2496.csv"), encoding = "UTF-8")

measured <- list()   # 每队列实测基因（universe 用）
res_list <- list()

for (cf in sort(list.files(COH, pattern = "\\.rds$", full.names = TRUE))) {
  coh <- sub("\\.rds$", "", basename(cf))
  d <- readRDS(cf)
  z  <- zscore(d$expr)
  cl <- d$clin
  cl <- cl[match(colnames(z), cl$sample), , drop = FALSE]
  time <- suppressWarnings(as.numeric(cl$time))
  status <- suppressWarnings(as.integer(cl$status))
  keep <- !is.na(time) & !is.na(status)
  z <- z[, keep, drop = FALSE]; time <- time[keep]; status <- status[keep]

  measured[[coh]] <- rownames(z)

  kg <- intersect(kegg, rownames(z))
  P <- as.numeric(colMeans(z[kg, , drop = FALSE]))
  oc <- oriented_cindex(time, status, P)
  Cp <- oc[["C"]]

  # null：B=1000 个 n=length(kg) 随机基因集（S1 的 n_measured 是签名大小 3-57，
  # 无 124 基因组 null 行可借——39 的 null_plain_C 同款循环，单规模一次）
  pool <- rownames(z)
  set.seed(20261109)   # 固定种子：入库后重跑必须逐位可复现
  nv <- replicate(1000, {
    g <- sample(pool, length(kg))
    oriented_cindex(time, status, as.numeric(colMeans(z[g, , drop = FALSE])))[["C"]]
  })
  nv <- nv[!is.na(nv)]
  null_med <- median(nv)
  d_rand_P_alone <- Cp - null_med

  # P 打赢多少 % 签名（同队列 S1 C 分布，evaluation 级）
  sig_c <- s1[cohort == coh]$C
  pct_beaten <- 100 * mean(Cp > sig_c, na.rm = TRUE)

  res_list[[length(res_list) + 1]] <- data.table(
    cohort = coh, n_kegg_measured = length(kg), C_P = Cp,
    flipped = oc[["flipped"]], null_med_124 = null_med,
    d_rand_P_alone = d_rand_P_alone,
    median_sig_C = median(sig_c, na.rm = TRUE),
    pct_sigs_beaten = pct_beaten)
  message("  ", coh, ": C(P)=", sprintf("%.4f", Cp),
          " null=", sprintf("%.4f", null_med),
          " beats ", sprintf("%.1f%%", pct_beaten))
}

res <- rbindlist(res_list)
fwrite(res, file.path(OUTD, "prolif_head_to_head.csv"))

# ---------------- universe（8 队列并集＋计数） ----------------
allg <- unique(unlist(measured))
cnt <- sapply(allg, function(g) sum(sapply(measured, function(m) g %in% m)))
uni <- data.table(gene = names(cnt), n_cohorts_measured = as.integer(cnt))
setorder(uni, -n_cohorts_measured, gene)
fwrite(uni, file.path(OUTD, "measured_universe.csv"))
message("universe: ", nrow(uni), " 个唯一基因")

# ---------------- 汇总 ----------------
# 合并口径：逐评估比较 C_P(队列) vs C_sig —— 按 S1 行的 cohort 关联 C_P
s1m <- merge(s1, res[, .(cohort, C_P)], by = "cohort", all.x = TRUE)
stopifnot(nrow(s1m) == nrow(s1))
pooled_beat <- 100 * mean(s1m$C_P > s1m$C, na.rm = TRUE)
lines <- c(
  "════ 增殖评分 P 头对头（P vs 312 签名）════",
  sprintf("跨队列合并：C(P) median = %.4f | 签名 C median (evaluation 级) = %.4f",
          median(res$C_P), median(s1$C, na.rm = TRUE)),
  sprintf("逐评估合并：P 打赢 %.1f%% 的 2,496 次评估", pooled_beat),
  sprintf("null：B=1000 随机基因集（规模=各队列 KEGG 实测数 118-124）"),
  "",
  capture.output(print(res, digits = 4)))
writeLines(lines, file.path(OUTD, "prolif_head_to_head_summary.txt"))
cat(paste(lines, collapse = "\n"), "\n")
