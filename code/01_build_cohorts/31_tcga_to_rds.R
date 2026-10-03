#!/usr/bin/env Rscript
# =====================================================================
# 31_tcga_to_rds.R —— 把 30_tcga_xena.py 的输出转为与 GEO 队列【同构】的 rds
#
# 结构契约（对照 cohorts/GSE87211.rds 实读结果）：
#   list(expr = matrix<double> 基因×样本, clin = data.frame(sample, time,
#        status, stage, age, sex, site))
#   expr rownames = HGNC symbol；time 单位 = 月；sex = F/M 大写；site 可 NA
#
# 三个必须对齐的点，错一个都会让 07 静默跑偏：
#   1. expr 必须是【基因×样本】（30 脚本为了好写落成样本×基因，这里要 t() 回来）
#   2. sex 必须折成 F/M（GEO 队列口径），male/female 会让七个队列的哑变量基准不一致
#   3. clin 行序必须与 expr 列序严格对齐（179 行 stopifnot 校验）
# =====================================================================
suppressPackageStartupMessages({ library(data.table) })

E <- "cohorts_raw/TCGA_expr.tsv.gz"
C <- "cohorts_raw/TCGA_clin.csv"
if (!file.exists(E)) stop("缺 ", E, "，请先跑 30_tcga_xena.py")
if (!file.exists(C)) stop("缺 ", C, "，请先跑 30_tcga_xena.py")

cat("① 读表达 ...\n")
# 首列是样本 barcode，直接交给 read.delim 当 row.names，比 fread 手动摘首列更稳
X <- read.delim(gzfile(E), row.names = 1, check.names = FALSE,
                stringsAsFactors = FALSE)
X <- as.matrix(X)
storage.mode(X) <- "double"
m <- t(X)                                  # -> 基因 x 样本
dimnames(m) <- list(colnames(X), rownames(X))
rm(X); invisible(gc())
cat("   表达", nrow(m), "基因 x", ncol(m), "样本\n")

cat("② 读临床 ...\n")
cl <- fread(C, data.table = FALSE)
cat("   临床", nrow(cl), "行\n")

# sex 折成 F/M（GEO 口径）
cl$sexM <- toupper(substr(as.character(cl$sex), 1, 1))
cl$sexM[!cl$sexM %in% c("F", "M")] <- NA_character_

out <- data.frame(
  sample = as.character(cl$sample),
  time   = suppressWarnings(as.numeric(cl$time)),
  status = suppressWarnings(as.integer(cl$status)),
  stage  = suppressWarnings(as.numeric(cl$stage)),
  age    = suppressWarnings(as.numeric(cl$age)),
  sex    = cl$sexM,
  site   = as.character(cl$site),
  stringsAsFactors = FALSE
)

cat("③ 对齐 ...\n")
idx <- match(colnames(m), out$sample)
stopifnot(!any(is.na(idx)))
out <- out[idx, , drop = FALSE]
rownames(out) <- NULL
stopifnot(nrow(out) == ncol(m))
stopifnot(identical(out$sample, colnames(m)))

ok <- !is.na(out$time) & !is.na(out$status)
cat("   有 OS 结局 ", sum(ok), " / ", ncol(m),
    "；事件 ", sum(out$status == 1, na.rm = TRUE), "\n", sep = "")
cat("   分期可用 ", sum(!is.na(out$stage)),
    "；性别可用 ", sum(out$sex %in% c("F", "M")),
    "；年龄可用 ", sum(!is.na(out$age)), "\n", sep = "")
cat("   随访中位 ", round(median(out$time, na.rm = TRUE), 1), " 月，最长 ",
    round(max(out$time, na.rm = TRUE), 0), " 月\n", sep = "")
cat("   >=60 月（5 年）且仍存活 ", sum(out$time >= 60 & out$status == 0, na.rm = TRUE),
    " 人 <- 做 5 年 AUC 的对照组，太少会不稳\n", sep = "")

dir.create("cohorts", showWarnings = FALSE)
saveRDS(list(expr = m, clin = out), "cohorts/TCGA_COADREAD.rds")
cat("\n④ 已写出 cohorts/TCGA_COADREAD.rds\n")

# 签名基因覆盖自检：tcga 在定稿前就该知道有多少签名基因测得到
libf <- "sigminer/output/signature_library_main.csv"
if (file.exists(libf)) {
  lib <- fread(libf, data.table = FALSE)
  gs <- unique(unlist(strsplit(lib$genes, ";")))
  cat("   签名基因 ", length(gs), " 个，其中在本队列测到 ",
      length(intersect(gs, rownames(m))), " 个 (",
      round(100 * length(intersect(gs, rownames(m))) / length(gs), 1), "%)\n", sep = "")
}
