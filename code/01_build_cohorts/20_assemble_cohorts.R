#!/usr/bin/env Rscript
# ============================================================
# 20_assemble_cohorts.R
# 把 Python 侧构建好的队列（GSE87211 / TCGA）转成与 06b 完全一致的
# cohorts/*.rds 结构：list(expr = 基因×样本矩阵, clin = 临床 data.frame)
#
# 统一约定（与 06b 对齐，否则 07 会读错）：
#   expr : 行=基因符号，列=样本；已 log 尺度
#   clin : sample, time, status, stage, age, sex, site
# ============================================================
suppressPackageStartupMessages({library(data.table)})

BASE <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)))
if (!length(BASE)) BASE <- getwd()
RAW  <- file.path(BASE, "cohorts_raw")
COH  <- file.path(BASE, "cohorts")
dir.create(COH, showWarnings = FALSE)

normalise_clin <- function(cl) {
  # 强制补齐 07 需要的全部列，缺失填 NA
  need <- c("sample", "time", "status", "stage", "age", "sex", "site")
  for (k in need) if (!k %in% names(cl)) cl[[k]] <- NA
  cl$sample <- as.character(cl$sample)
  cl$time   <- suppressWarnings(as.numeric(cl$time))
  cl$status <- suppressWarnings(as.integer(cl$status))
  cl$stage  <- suppressWarnings(as.numeric(cl$stage))
  cl$age    <- suppressWarnings(as.numeric(cl$age))
  cl$sex    <- as.character(cl$sex)
  cl$site   <- as.character(cl$site)
  cl[, need, drop = FALSE]
}

build <- function(name) {
  fe <- file.path(RAW, paste0(name, "_expr.csv.gz"))
  fc <- file.path(RAW, paste0(name, "_clin.csv"))
  if (!file.exists(fe) || !file.exists(fc)) {
    cat(sprintf("[%s] 缺少输入，跳过\n", name)); return(invisible(NULL))
  }
  cat(sprintf("[%s] 读表达 ...\n", name))
  expr <- fread(fe, data.table = FALSE, showProgress = FALSE)
  gn   <- expr[[1]]
  expr <- as.matrix(expr[, -1, drop = FALSE])
  rownames(expr) <- gn
  storage.mode(expr) <- "double"

  cat(sprintf("[%s] 读临床 ...\n", name))
  cl <- fread(fc, data.table = FALSE, encoding = "UTF-8")
  cl <- normalise_clin(cl)
  cl <- cl[match(colnames(expr), cl$sample), , drop = FALSE]

  # 生存字段自检
  ok <- !is.na(cl$time) & !is.na(cl$status) & cl$time > 0
  cat(sprintf("[%s] 基因 %d × 样本 %d | 可用生存 %d 例，事件 %d\n",
              name, nrow(expr), ncol(expr), sum(ok), sum(cl$status[ok], na.rm = TRUE)))
  if (sum(ok) < 30 || sum(cl$status[ok], na.rm = TRUE) < 10) {
    cat(sprintf("[%s] !! 事件数过少，不写入\n", name)); return(invisible(NULL))
  }
  saveRDS(list(expr = expr, clin = cl), file.path(COH, paste0(name, ".rds")))
  cat(sprintf("[%s] 已写入 cohorts/%s.rds\n", name, name))
}

for (nm in c("GSE87211", "TCGA_CRC")) build(nm)
cat("\n完成。当前 cohorts 目录：\n")
print(list.files(COH, pattern = "\\.rds$"))
