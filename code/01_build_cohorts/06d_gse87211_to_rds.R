#!/usr/bin/env Rscript
# ============================================================
# 06d_gse87211_to_rds.R
# 把 06c（python）产出的 GSE87211_expr.csv.gz + GSE87211_clin.csv
# 转成 39_proliferation_adjust.R 消费的 cohorts/GSE87211.rds 格式
#   expr: 基因(行) × 样本(列) 矩阵；clin: data.frame(sample,time,status,stage,age,sex,...)
# ============================================================
suppressPackageStartupMessages(library(data.table))

BASE <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)))
if (!length(BASE)) BASE <- getwd()
ROOT <- normalizePath(file.path(BASE, "..", ".."))
RAW  <- file.path(ROOT, "cohorts_raw")
OUT  <- file.path(ROOT, "cohorts")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

e <- fread(file.path(RAW, "GSE87211_expr.csv.gz"))
stopifnot("V1" %in% names(e))
genes <- e$V1
mat   <- as.matrix(e[, -1, with = FALSE])
rownames(mat) <- genes
storage.mode(mat) <- "double"
message("expr: ", nrow(mat), " 基因 x ", ncol(mat), " 样本")

c <- fread(file.path(RAW, "GSE87211_clin.csv"))
# 列：sample,time,status,dfs_time,dfs_status,age,sex,stage
# 39 用 OS (time,status)；stage/age/sex 为协变量
c[, sample := as.character(sample)]
miss <- setdiff(colnames(mat), c$sample)
if (length(miss)) stop("clin 缺样本: ", head(miss, 5))
c2 <- c[match(colnames(mat), sample)]
out <- list(expr = mat,
            clin = data.frame(sample = c2$sample,
                              time = as.numeric(c2$time),
                              status = as.integer(c2$status),
                              stage = as.integer(c2$stage),
                              age = as.numeric(c2$age),
                              sex = as.character(c2$sex),
                              stringsAsFactors = FALSE))
n_ok <- sum(!is.na(out$clin$time) & !is.na(out$clin$status))
message("OS 可用: ", n_ok, " / ", nrow(out$clin), " 事件 ", sum(out$clin$status == 1, na.rm = TRUE))
saveRDS(out, file.path(OUT, "GSE87211.rds"))
message("-> ", file.path(OUT, "GSE87211.rds"))
