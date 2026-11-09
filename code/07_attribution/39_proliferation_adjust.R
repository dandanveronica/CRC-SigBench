#!/usr/bin/env Rscript
# ============================================================
# 39_proliferation_adjust.R —— POST HOC 增殖混淆校正（不在注册内）
#
# 目的：检验签名的微小信号（ΔC_random +0.0107 / ΔC_clinical +0.0045）
#       在多大程度上由增殖/分期生物学解释（"little biology beyond the stage"）
#
# 设计：
#   P = 逐患者增殖评分 = 队列内 KEGG hsa04110（124 基因）z-score 等权平均
#       —— 与签名打分完全同口径（zscore -> colMeans）
#   主指标（新）：
#       dC_clin_P  = C(clin + P + S) - C(clin + P)   [超出分期+增殖的增量]
#   对照（同口径随机集，缓存按 cohort×n_used）：
#       dC_null_P  = C(clin + P + rand) - C(clin + P) 的中位数（B=100）
#       d_rand_P   = dC_clin_P - median(dC_null_P)
#   互证（必须与主管线一致，否则口径复刻失败）：
#       dC_clin    = delta_cindex() 逐行复制 07_audit_signatures.R
#       d_rand     = C - median(null C)             [B=1000，与原 B_PERM 一致]
#       逐行对 S1 的 delta_C / null_median 报最大绝对差
#
# 输出：audit_out_clean/proliferation/
#   per_row_prolif.csv       逐签名×逐队列
#   summary_prolif.txt       汇总（GEO / TCGA / 全部子集中位数）
#
# 🔴 本分析为 post hoc，不改变任何注册判据与主结果
# ============================================================
suppressPackageStartupMessages({
  library(survival); library(data.table)
})

BASE <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)))
if (!length(BASE)) BASE <- getwd()
ROOT <- normalizePath(file.path(BASE, "..", ".."))
COH  <- file.path(ROOT, "cohorts")
OUTD <- file.path(ROOT, "audit_out_clean", "proliferation")
dir.create(OUTD, showWarnings = FALSE, recursive = TRUE)

B_NULL_ADJ <- 100   # 调整后零分布重复数（post-hoc，她拍板 n=100）
B_NULL_PLN <- 1000  # 互证用，与主管线 B_PERM 一致
MIN_SAMPLES <- 40

# ---------------- KEGG 增殖基因集 ----------------
KEGG_FILE <- file.path(ROOT, "data", "proliferation", "kegg_hsa04110_symbols.txt")
kegg <- trimws(readLines(KEGG_FILE))
kegg <- kegg[nchar(kegg) > 0 & !startsWith(kegg, "#")]
message("KEGG hsa04110 基因集：", length(kegg), " 个符号")

# ---------------- 签名库 ----------------
sig <- fread(file.path(ROOT, "data", "signatures", "signature_library_312_genes.csv"),
             encoding = "UTF-8")
message("签名库：", nrow(sig))

# ---------------- 逐函数复制的口径（07_audit_signatures.R 原文） ----------------
`%||%` <- function(a, b) if (is.null(a)) b else a

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
  if (c1 > c0) c(C = c1, flipped = 1) else c(C = c0, flipped = 0)
}

risk_score <- function(z, genes, coef = NULL, direction = 1) {
  g <- intersect(genes, rownames(z))
  if (length(g) == 0) return(NULL)
  sub <- z[g, , drop = FALSE]
  if (is.null(coef)) {
    w <- rep(1, length(g))
  } else {
    w <- coef[match(g, genes)]
    w[is.na(w)] <- 1
  }
  as.numeric(direction * crossprod(w, sub)) / sum(abs(w))
}

# —— 与 07 逐字一致（互证 dC_clin 用）——
delta_cindex <- function(time, status, pred, clin_df) {
  d <- clin_df
  d$.t <- time; d$.s <- status; d$.p <- pred
  ok <- complete.cases(d[, c(".t", ".s", ".p")])
  d <- d[ok, , drop = FALSE]
  if (nrow(d) < 40 || sum(d$.s) < 5) return(NA_real_)
  vars <- intersect(c("stage", "age", "sex", "site"), colnames(d))
  vars <- vars[sapply(d[, vars, drop = FALSE], function(x) length(unique(x)) > 1)]
  if (length(vars) == 0) return(NA_real_)
  need <- c(".t", ".s", ".p", vars)
  keepcc <- stats::complete.cases(d[, need, drop = FALSE])
  if (sum(keepcc) < 40 || sum(d$.s[keepcc]) < 5) return(NA_real_)
  dd <- d[keepcc, , drop = FALSE]
  f_base <- as.formula(paste("Surv(.t, .s) ~", paste(vars, collapse = " + ")))
  f_full <- as.formula(paste("Surv(.t, .s) ~ .p +", paste(vars, collapse = " + ")))
  cb <- try(coxph(f_base, data = dd), silent = TRUE)
  cf <- try(coxph(f_full, data = dd), silent = TRUE)
  if (inherits(cb, "try-error") || inherits(cf, "try-error")) return(NA_real_)
  lp_b <- predict(cb, type = "lp"); lp_f <- predict(cf, type = "lp")
  cb_i <- cindex(dd$.t, dd$.s, lp_b)
  cf_i <- cindex(dd$.t, dd$.s, lp_f)
  if (is.na(cb_i) || is.na(cf_i)) return(NA_real_)
  cf_i - cb_i
}

# —— 新指标：base 与 full 都含 P，其余口径完全一致 ——
delta_cindex_P <- function(time, status, pred, clin_df, P) {
  d <- clin_df
  d$.t <- time; d$.s <- status; d$.p <- pred; d$.P <- P
  if (length(unique(d$.P[!is.na(d$.P)])) < 2) return(NA_real_)
  ok <- complete.cases(d[, c(".t", ".s", ".p", ".P")])
  d <- d[ok, , drop = FALSE]
  if (nrow(d) < 40 || sum(d$.s) < 5) return(NA_real_)
  vars <- intersect(c("stage", "age", "sex", "site"), colnames(d))
  vars <- vars[sapply(d[, vars, drop = FALSE], function(x) length(unique(x)) > 1)]
  base_rhs <- if (length(vars)) paste(c(".P", vars), collapse = " + ") else ".P"
  need <- c(".t", ".s", ".p", ".P", vars)
  keepcc <- stats::complete.cases(d[, need, drop = FALSE])
  if (sum(keepcc) < 40 || sum(d$.s[keepcc]) < 5) return(NA_real_)
  dd <- d[keepcc, , drop = FALSE]
  f_base <- as.formula(paste("Surv(.t, .s) ~", base_rhs))
  f_full <- as.formula(paste("Surv(.t, .s) ~ .p +", base_rhs))
  cb <- try(coxph(f_base, data = dd), silent = TRUE)
  cf <- try(coxph(f_full, data = dd), silent = TRUE)
  if (inherits(cb, "try-error") || inherits(cf, "try-error")) return(NA_real_)
  lp_b <- predict(cb, type = "lp"); lp_f <- predict(cf, type = "lp")
  cb_i <- cindex(dd$.t, dd$.s, lp_b)
  cf_i <- cindex(dd$.t, dd$.s, lp_f)
  if (is.na(cb_i) || is.na(cf_i)) return(NA_real_)
  cf_i - cb_i
}

# ---------------- 零分布缓存（按 cohort×n_used） ----------------
nc_plain <- new.env(hash = TRUE)   # C 零分布（B=1000，互证）
nc_adj   <- new.env(hash = TRUE)   # 调整后增量零分布（B=100）

null_plain_C <- function(coh, n_gene, z, time, status) {
  key <- paste0(coh, "|", n_gene)
  if (!is.null(nc_plain[[key]])) return(nc_plain[[key]])
  pool <- rownames(z)
  if (length(pool) < n_gene * 2) { v <- NULL; assign(key, v, nc_plain); return(v) }
  vals <- numeric(B_NULL_PLN)
  for (b in seq_len(B_NULL_PLN)) {
    g <- sample(pool, n_gene)
    p <- as.numeric(colMeans(z[g, , drop = FALSE]))
    vals[b] <- oriented_cindex(time, status, p)[["C"]]
  }
  vals <- vals[!is.na(vals)]
  assign(key, vals, nc_plain)
  vals
}

null_adj_incr <- function(coh, n_gene, z, time, status, cl2, P) {
  key <- paste0(coh, "|", n_gene)
  if (!is.null(nc_adj[[key]])) return(nc_adj[[key]])
  pool <- rownames(z)
  if (length(pool) < n_gene * 2) { v <- NULL; assign(key, v, nc_adj); return(v) }
  vals <- numeric(B_NULL_ADJ)
  for (b in seq_len(B_NULL_ADJ)) {
    g <- sample(pool, n_gene)
    p <- as.numeric(colMeans(z[g, , drop = FALSE]))
    oc <- oriented_cindex(time, status, p)
    if (is.na(oc[["C"]])) { vals[b] <- NA; next }
    p_o <- if (identical(as.numeric(oc[["flipped"]]), 1)) -p else p
    vals[b] <- delta_cindex_P(time, status, p_o, cl2, P)
  }
  vals <- vals[!is.na(vals)]
  assign(key, vals, nc_adj)
  vals
}

# ---------------- 主循环 ----------------
res_list <- list()
summary_lines <- c()
for (cf in sort(list.files(COH, pattern = "\\.rds$", full.names = TRUE))) {
  coh <- sub("\\.rds$", "", basename(cf))
  d <- readRDS(cf)
  z  <- zscore(d$expr)
  cl <- d$clin
  cl <- cl[match(colnames(z), cl$sample), , drop = FALSE]
  time <- suppressWarnings(as.numeric(cl$time))
  status <- suppressWarnings(as.integer(cl$status))
  keep <- !is.na(time) & !is.na(status)
  if (sum(keep) < MIN_SAMPLES) { message("[跳过] ", coh); next }
  z <- z[, keep, drop = FALSE]; time <- time[keep]; status <- status[keep]
  cl2 <- cl[keep, , drop = FALSE]

  # 增殖评分 P（与签名同口径：队列内 z-score 等权平均）
  kg <- intersect(kegg, rownames(z))
  message("  ", coh, "：样本 ", ncol(z), " | KEGG 测得 ", length(kg), "/", length(kegg))
  if (length(kg) < 10) { message("  !! KEGG 基因测得过少，跳过该队列"); next }
  P <- as.numeric(colMeans(z[kg, , drop = FALSE]))
  # 增殖 vs 分期（分期为数值化 I<II<III<IV）
  # GSE72970 无分期字段（数据事实），相关系数无配对时置 NA，不中断（07 原管线无此行）
  stg <- suppressWarnings(as.numeric(as.character(cl2$stage)))
  rho_stage <- NA_real_
  if (sum(complete.cases(stg)) >= 3 &&
      length(unique(stg[!is.na(stg)])) > 1 &&
      sum(complete.cases(P) & !is.na(stg)) >= 3) {
    rho_stage <- suppressWarnings(cor(P, stg, method = "spearman",
                                      use = "complete.obs"))
  }

  for (i in seq_len(nrow(sig))) {
    genes <- trimws(unlist(strsplit(sig$genes[i], "[;,]")))
    genes <- genes[genes != ""]
    p <- risk_score(z, genes, NULL, 1)
    if (is.null(p)) next
    miss <- sum(!genes %in% rownames(z))
    n_used <- length(genes) - miss
    oc  <- oriented_cindex(time, status, p)
    if (is.na(oc[["C"]])) next
    p_o <- if (identical(as.numeric(oc[["flipped"]]), 1)) -p else p

    dC_clin   <- delta_cindex(time, status, p_o, cl2)
    dC_clin_P <- delta_cindex_P(time, status, p_o, cl2, P)
    rho_SP    <- suppressWarnings(cor(p_o, P, method = "spearman",
                                      use = "complete.obs"))

    nv  <- null_plain_C(coh, n_used, z, time, status)
    nva <- null_adj_incr(coh, n_used, z, time, status, cl2, P)

    res_list[[length(res_list) + 1]] <- data.table(
      sig_id = sig$sig_id[i], cohort = coh, n_used = n_used,
      C = oc[["C"]],
      dC_clin = dC_clin, dC_clin_P = dC_clin_P,
      rho_SP = rho_SP,
      null_median = if (is.null(nv)) NA_real_ else median(nv),
      null_median_adj = if (is.null(nva)) NA_real_ else median(nva),
      d_rand = if (is.null(nv) || is.na(oc[["C"]])) NA_real_ else oc[["C"]] - median(nv),
      d_rand_P = if (is.null(nva) || is.na(dC_clin_P)) NA_real_ else dC_clin_P - median(nva),
      rho_P_stage = rho_stage)
  }
  summary_lines <- c(summary_lines,
    sprintf("%s: KEGG 测得 %d | rho(P, stage) = %s", coh, length(kg),
            sprintf("%.3f", rho_stage)))
  message("  -> ", coh, " 完成")
}

res <- rbindlist(res_list, fill = TRUE)
fwrite(res, file.path(OUTD, "per_row_prolif.csv"))

# ---------------- 汇总 ----------------
summ <- res[, .(
  n = .N,
  med_C = median(C, na.rm = TRUE),
  med_dC_clin = median(dC_clin, na.rm = TRUE),
  med_dC_clin_P = median(dC_clin_P, na.rm = TRUE),
  med_d_rand = median(d_rand, na.rm = TRUE),
  med_d_rand_P = median(d_rand_P, na.rm = TRUE),
  med_rho_SP = median(rho_SP, na.rm = TRUE),
  pct_rho_pos = 100 * mean(rho_SP > 0, na.rm = TRUE),
  pct_drandP_pos = 100 * mean(d_rand_P > 0, na.rm = TRUE)
), by = .(is_tcga = grepl("TCGA", cohort))]
setnames(summ, "is_tcga", "is_TCGA")

lines <- c(
  "════ POST HOC 增殖混淆校正汇总 ════",
  paste(summary_lines, collapse = "\n"),
  "",
  "口径：行 = 签名×队列；P = KEGG hsa04110 124 基因 z 均值",
  "      dC_clin   = C(clin+S) − C(clin)         [互证 S1.delta_C]",
  "      dC_clin_P = C(clin+P+S) − C(clin+P)     [超出分期+增殖]",
  "      d_rand    = C − median(null C), B=1000  [互证 S1.d_rand]",
  "      d_rand_P  = dC_clin_P − median(null 增量), B=100",
  "",
  capture.output(print(summ, digits = 4)),
  "",
  sprintf("全表中位数: med_C=%.4f med_dC_clin=%+.4f med_dC_clin_P=%+.4f | med_d_rand=%+.4f med_d_rand_P=%+.4f | med_rho_SP=%.3f (%.1f%%>0)",
          median(res$C, na.rm = TRUE), median(res$dC_clin, na.rm = TRUE),
          median(res$dC_clin_P, na.rm = TRUE), median(res$d_rand, na.rm = TRUE),
          median(res$d_rand_P, na.rm = TRUE), median(res$rho_SP, na.rm = TRUE),
          100 * mean(res$rho_SP > 0, na.rm = TRUE))
)
writeLines(lines, file.path(OUTD, "summary_prolif.txt"))
cat(paste(lines, collapse = "\n"), "\n")
message("输出：", OUTD)
