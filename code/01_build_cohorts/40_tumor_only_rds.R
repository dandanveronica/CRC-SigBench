#!/usr/bin/env Rscript
# ============================================================
# 40_tumor_only_rds.R —— 恢复 312-audit 底座：tumor-only rds
#
# 背景（2026-11-09 互证发现的归因）：
#   per_cohort_2496（312×8 audit）用的底座是 2026-10-16/17 会话做的
#   tumor-only 重建（Table 1 口径：GSE39582 566/556/187、GSE87211 203/190/28）。
#   本次 06b refresh（全样本 585/363）覆盖了这两个 rds，导致 39 增殖
#   校正在这两个队列行级偏移（其余 6 队列无 normal 组织，逐行一致）。
#
# 本脚本按当年口径从 series matrix characteristics 过滤：
#   GSE39582: dataset 字段 != "Non Tumoral"      （585 -> 566）
#   GSE87211: tissue 字段 == "rectal tumor"      （363 -> 203）
# 只动样本子集，expr/注释/临床列一律沿用现有 rds。
#
# 校验锚点（Table 1）：GSE39582 566（OS 可用 556、事件 187）；
#                     GSE87211 203（OS 可用 190、事件 28）
# ============================================================
suppressPackageStartupMessages(library(data.table))

BASE <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)))
if (!length(BASE)) BASE <- getwd()
ROOT <- normalizePath(file.path(BASE, "..", ".."))
RAW  <- file.path(ROOT, "cohorts_raw")
COH  <- file.path(ROOT, "cohorts")

# 从 series matrix 逐格提取指定键的 {GSM: value}
char_map <- function(gz, key_lower) {
  con <- gzfile(gz, "rt"); on.exit(close(con))
  gsm <- NULL; rows <- list()
  repeat {
    line <- readLines(con, n = 1, warn = FALSE)
    if (!length(line)) break
    if (startsWith(line, "!Sample_geo_accession")) {
      gsm <- trimws(gsub('"', '', strsplit(sub("^!Sample_geo_accession\t", "", line), "\t")[[1]]))
    } else if (startsWith(line, "!Sample_characteristics_ch1")) {
      rows[[length(rows) + 1]] <- trimws(gsub('"', '', strsplit(sub("^!Sample_characteristics_ch1\t", "", line), "\t")[[1]]))
    } else if (startsWith(line, "!series_matrix_table_begin")) {
      break
    }
  }
  out <- rep(NA_character_, length(gsm)); names(out) <- gsm
  for (row in rows) {
    for (j in seq_along(row)) {
      if (j > length(gsm)) break
      p <- regexpr(":", row[j], fixed = TRUE)
      if (p < 2) next
      k <- trimws(substr(row[j], 1, p - 1))
      if (tolower(k) == key_lower) out[gsm[j]] <- trimws(substr(row[j], p + 1, nchar(row[j])))
    }
  }
  out
}

rebuild <- function(gse, gz, key, keep_fn, expect_n, exp_ok, exp_ev) {
  lab <- char_map(file.path(RAW, gz), key)
  d <- readRDS(file.path(COH, paste0(gse, ".rds")))
  smp <- colnames(d$expr)
  lv <- lab[smp]
  keep <- !is.na(lv) & keep_fn(lv)
  cat(sprintf("[%s] 标签值: %s\n", gse, paste(names(table(lv, useNA = "ifany")),
      table(lv, useNA = "ifany"), collapse = " / ")))
  d$expr <- d$expr[, keep, drop = FALSE]
  d$clin <- d$clin[match(colnames(d$expr), d$clin$sample), , drop = FALSE]
  n <- ncol(d$expr)
  ok <- !is.na(d$clin$time) & !is.na(d$clin$status)
  ev <- sum(d$clin$status[ok] == 1, na.rm = TRUE)
  cat(sprintf("[%s] -> %d 样本 | OS 可用 %d | 事件 %d  (锚点: %d/%d/%d)\n",
              gse, n, sum(ok), ev, expect_n, exp_ok, exp_ev))
  stopifnot(n == expect_n, sum(ok) == exp_ok, ev == exp_ev)
  saveRDS(d, file.path(COH, paste0(gse, ".rds")))
  cat(sprintf("[%s] 已覆盖 %s\n", gse, file.path(COH, paste0(gse, ".rds"))))
}

rebuild2 <- function(gse, gz, key, keep_fn, expect_n, exp_ok, exp_ev,
                     drop_bad_time = FALSE) {
  lab <- char_map(file.path(RAW, gz), key)
  d <- readRDS(file.path(COH, paste0(gse, ".rds")))
  smp <- colnames(d$expr)
  lv <- lab[smp]
  keep <- !is.na(lv) & keep_fn(lv)
  cat(sprintf("[%s] 标签值: %s\n", gse, paste(names(table(lv, useNA = "ifany")),
      table(lv, useNA = "ifany"), collapse = " / ")))
  d$expr <- d$expr[, keep, drop = FALSE]
  d$clin <- d$clin[match(colnames(d$expr), d$clin$sample), , drop = FALSE]
  if (drop_bad_time) {
    # GSE39582 当年口径 = 分析集（OS 完整且 time>0）；GSE87211 保留全集
    bad <- is.na(d$clin$time) | is.na(d$clin$status) | d$clin$time <= 0
    cat(sprintf("[%s] 剔除无 OS/非正时间 %d 例\n", gse, sum(bad)))
    d$expr <- d$expr[, !bad, drop = FALSE]
    d$clin <- d$clin[!bad, , drop = FALSE]
  }
  n <- ncol(d$expr)
  ok <- !is.na(d$clin$time) & !is.na(d$clin$status)
  ev <- sum(d$clin$status[ok] == 1, na.rm = TRUE)
  cat(sprintf("[%s] -> %d 样本 | OS 可用 %d | 事件 %d  (锚点: %d/%d/%d)\n",
              gse, n, sum(ok), ev, expect_n, exp_ok, exp_ev))
  stopifnot(n == expect_n, sum(ok) == exp_ok, ev == exp_ev)
  saveRDS(d, file.path(COH, paste0(gse, ".rds")))
  cat(sprintf("[%s] 已覆盖 %s\n", gse, file.path(COH, paste0(gse, ".rds"))))
}

# Table1_定稿核验报告_20261017 口径：
#   GSE39582 rds = 556（分析集，OS 全完整、187 事件）
#   GSE87211 rds = 203 行（肿瘤-only 全集，13 例 time NA，190 OS / 28 事件）
rebuild2("GSE39582", "GSE39582_series_matrix.txt.gz", "dataset",
         function(v) v != "Non Tumoral", 556, 556, 187, drop_bad_time = TRUE)
rebuild2("GSE87211", "GSE87211_series_matrix.txt.gz", "tissue",
         function(v) v == "rectal tumor", 203, 190, 28)
cat("\nDONE — tumor-only 底座恢复\n")
