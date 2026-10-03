#!/usr/bin/env Rscript
# =====================================================================
# 06b_build_cohorts.R
# 由 06b_download.sh 下载的 GEO series matrix 构建 07 所需的 cohorts/*.rds
#
# 与旧 06 的关键差别：
#   1) 不依赖 GEOquery（SOFT 文件大、易超时），直接解析 series matrix
#   2) 临床字段从 !Sample_characteristics_ch1 的 "key: value" 里提取，
#      并为每个队列配置字段优先级——旧 06 用固定 key 表，实测匹配不上
#      （GSE29621 的生存信息写在 "overall survival (os):" / "os event:" 里）
#   3) 探针->基因符号：优先用本地 GPL*.annot.gz，缺失时用 Bioconductor 注释包
#
# 用法：Rscript 06b_build_cohorts.R [原始目录] [输出目录]
# =====================================================================

suppressPackageStartupMessages({
  for (p in c("data.table")) if (!requireNamespace(p, quietly = TRUE))
    install.packages(p, repos = "https://cloud.r-project.org")
  library(data.table)
})

args  <- commandArgs(trailingOnly = TRUE)
RAW   <- ifelse(length(args) >= 1, args[1], "cohorts_raw")
OUT   <- ifelse(length(args) >= 2, args[2], "cohorts")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------- 临床字段映射表
# 每个队列自己一套 key 优先级（小写、去空格后匹配）。留空则用通用启发式。
CLIN_MAP <- list(
  "GSE29621" = list(
    time   = c("overall survival (os)", "os", "dfs"),
    status = c("os event", "event", "vital status"),
    stage  = c("ajcc staging", "stage", "tnm stage"),
    sex    = c("gender", "sex")
  ),
  # 实测字段：age / gender / ethnicity / ajcc_stage / grade /
  #   overall_event (death from any cause) / dss_event / dfs_event /
  #   "overall survival follow-up time" / dss_time / dfs_time
  # overall_event 取值 death / no death
  "GSE17536" = list(
    time   = c("overall survival follow-up time", "dss_time", "overall survival", "survival_time"),
    status = c("overall_event", "dss_event", "survival_status", "event", "vital status"),
    stage  = c("ajcc_stage", "stage", "tnm stage"),
    sex    = c("gender", "sex"), age = c("age"), site = c("location", "site")
  ),
  # 与 GSE17536 同结构，但 FIELD 是锯齿状（同一行不同样本键不同）
  "GSE17537" = list(
    time   = c("overall survival follow-up time", "overall survival", "survival_time"),
    status = c("overall_event", "survival_status", "event", "dfs_event"),
    stage  = c("ajcc_stage", "ajcc stage", "stage"),
    sex    = c("gender", "sex"), age = c("age"), site = c("location", "site")
  ),
  # GSE39582 实测字段：rfs.event / rfs.delay / os.event / "os.delay (months)" /
  #   tnm.stage / "age.at.diagnosis (year)" / Sex / tumor.location / dataset
  # 注意：time 必须写成 "os.delay"，否则 grepl("os") 会先命中 os.event（0/1）
  "GSE39582" = list(
    time   = c("os.delay", "os delay", "overall survival", "dfs.delay", "survival time"),
    status = c("os.event", "os event", "dfs.event", "event", "vital status"),
    stage  = c("tnm.stage", "tnm stage", "stage", "ajcc stage"),
    sex    = c("sex", "gender"),
    age    = c("age.at.diagnosis", "age at diagnosis", "age"),
    site   = c("tumor.location", "tumor location", "location", "sidedness")
  ),
  # 实测字段：ajcc_stage / dfs_event / dfs_time / dss_event("0 (no death)"/"1 (death…)")
  #           / "dss_time (disease specific survival time, months)"
  "GSE38832" = list(
    time   = c("dss_time", "dfs_time", "overall survival", "survival time"),
    status = c("dss_event", "dfs_event", "event", "vital status"),
    stage  = c("ajcc_stage", "stage", "tnm stage"),
    sex    = c("sex", "gender")
  ),
  # 实测：os / os censored / pfs / pfs censored / tumor location / Sex / age
  # ⚠ os censored 取值 1=删失 0=事件 —— 必须取反
  "GSE72970" = list(
    time   = c("os", "overall survival", "survival time", "pfs"),
    status = c("os censored", "os_censored", "os event", "vital status"),
    stage  = c("stage", "ajcc_stage", "tnm stage"),
    sex    = c("sex", "gender"), age = c("age"),
    site   = c("tumor location", "tumor.location", "location")
  ),
  "GSE87211" = list(
    time   = c("os", "overall survival", "survival time"),
    status = c("os event", "event", "vital status"),
    stage  = c("stage", "tnm stage", "ajcc stage"),
    sex    = c("sex", "gender")
  )
)

# 生存状态需要取反的队列（1=删失 0=事件）
INVERT_STATUS <- c("GSE72970")

GENERIC <- list(
  time   = c("os", "overall survival", "overall survival (os)", "survival time",
             "survival_time", "os_time", "dfs", "dfs_time", "time"),
  status = c("os event", "event", "vital status", "vital_status",
             "survival_status", "status", "death"),
  stage  = c("tnm stage", "stage", "ajcc stage", "ajcc_stage", "pathologic stage"),
  age    = c("age", "age at diagnosis", "age_at_diagnosis"),
  sex    = c("sex", "gender"),
  site   = c("site", "tumor site", "location", "sidedness", "primary site")
)

pick <- function(keys, value_by_key) {
  for (k in keys) {
    hit <- names(value_by_key)[tolower(trimws(names(value_by_key))) == k]
    if (length(hit)) return(value_by_key[[hit[1]]])
    # 退一步：包含匹配
    hit <- names(value_by_key)[grepl(k, tolower(names(value_by_key)), fixed = TRUE)]
    if (length(hit)) return(value_by_key[[hit[1]]])
  }
  NA_character_
}

# ------------------------------------------------------- series matrix 解析
read_series_matrix <- function(gz) {
  con <- gzfile(gz, "rt")
  lines <- readLines(con, warn = FALSE)
  close(con)
  meta <- grep("^!", lines, value = TRUE)
  beg  <- grep("^!series_matrix_table_begin", lines)
  end  <- grep("^!series_matrix_table_end", lines)
  if (!length(beg) || !length(end)) stop("series matrix 格式异常：", gz)
  tab  <- lines[(beg[1] + 1):(end[1] - 1)]
  pline <- grep("^!Series_platform_id", meta, value = TRUE)[1]
  platform <- if (is.na(pline)) "" else
    trimws(gsub('"', '', sub('^!Series_platform_id\t', '', pline)))
  list(meta = meta, tab = tab, platform = platform)
}

parse_expr <- function(sm) {
  m <- fread(text = paste(sm$tab, collapse = "\n"), sep = "\t", header = TRUE,
             quote = "\"", fill = TRUE, na.strings = c("", "NA", "null", "NULL"),
             colClasses = "character", showProgress = FALSE)
  idcol <- colnames(m)[1]
  probes <- as.character(m[[idcol]])
  num <- m[, -1, with = FALSE]
  num <- num[, lapply(.SD, as.numeric)]
  mat <- as.matrix(num)
  rownames(mat) <- probes
  mat
}

parse_chars <- function(sm) {
  # 关键修正：GEO 的 characteristics 是【锯齿状】的——同一行里不同样本的
  # key 可能不一样（GSE17537 第 5 行：样本1 是 overall_event，样本2 却是 grade）。
  # 按「行」取 key 会整体错位，必须逐格解析 "key: value"。
  ch <- grep("^!Sample_characteristics_ch1", sm$meta, value = TRUE)
  gs <- grep("^!Sample_geo_accession", sm$meta, value = TRUE)
  if (!length(gs)) return(NULL)
  gsm <- gsub('"', '', strsplit(sub("^!Sample_geo_accession\t", "", gs[1]), "\t")[[1]])
  fields <- list()
  for (ln in ch) {
    v <- gsub('"', '', strsplit(sub("^!Sample_characteristics_ch1\t", "", ln), "\t")[[1]])
    if (length(v) != length(gsm)) next
    for (j in seq_along(v)) {
      s <- trimws(v[j])
      if (!nzchar(s)) next
      p <- regexpr(":", s, fixed = TRUE)
      if (p < 2) next
      k   <- trimws(substr(s, 1, p - 1))
      val <- trimws(substr(s, p + 1, nchar(s)))
      if (!nzchar(k)) next
      if (is.null(fields[[k]])) fields[[k]] <- rep(NA_character_, length(gsm))
      fields[[k]][j] <- val
    }
  }
  list(gsm = gsm, fields = fields)
}

# ------------------------------------------------------- 探针 -> 基因符号
annot_path <- function(platform) file.path(RAW, paste0(platform, ".annot.gz"))

ensure_annot <- function(platform) {
  fp <- annot_path(platform)
  if (file.exists(fp) && file.size(fp) > 1e5) return(fp)
  url <- paste0("https://ftp.ncbi.nlm.nih.gov/geo/platforms/",
                substr(platform, 1, 3), "nnn/", platform, "/annot/", platform, ".annot.gz")
  message("  下载平台注释 ", platform, " ...")
  ok <- try(download.file(url, fp, mode = "wb", quiet = TRUE), silent = TRUE)
  if (inherits(ok, "try-error") || !file.exists(fp) || file.size(fp) < 1e5) {
    message("  [警告] 注释下载失败，改用 Bioconductor 注释包")
    return(NULL)
  }
  fp
}

parse_annot <- function(fp) {
  # GEO annot.gz 结构： "#" 开头是列说明 -> "!platform_table_begin" -> 真正的列名行 -> 数据
  con <- gzfile(fp, "rt")
  lines <- readLines(con, warn = FALSE)
  close(con)
  beg <- grep("^!platform_table_begin", lines)
  if (!length(beg)) return(NULL)
  hdr <- strsplit(lines[beg[1] + 1], "\t")[[1]]
  end <- grep("^!platform_table_end", lines)
  data_start <- beg[1] + 2
  data_end   <- if (length(end)) end[1] - 1 else length(lines)
  ic <- grep("^ID$", hdr, ignore.case = TRUE)
  sc <- grep("^Gene symbol$|^Symbol$", hdr, ignore.case = TRUE)
  if (!length(ic) || !length(sc)) return(NULL)
  ic <- min(ic); sc <- min(sc)
  d <- fread(text = paste(lines[data_start:data_end], collapse = "\n"),
             sep = "\t", header = FALSE, fill = TRUE, quote = "",
             select = c(ic, sc), colClasses = "character", showProgress = FALSE)
  setnames(d, c("ID", "SYMBOL"))
  # 一个探针对应多个基因（"MIR4640///DDR1"）时取第一个，并丢掉空值
  d[, SYMBOL := sub("///.*$", "", SYMBOL)]
  d[!is.na(SYMBOL) & SYMBOL != "" & SYMBOL != "---"]
}

# 解析一次注释要几十秒，缓存成 rds，7 个队列共用
annot_map <- function(platform) {
  cache <- file.path(RAW, paste0(platform, "_symbol.rds"))
  if (file.exists(cache)) return(readRDS(cache))
  fp <- ensure_annot(platform)
  mp <- if (!is.null(fp)) parse_annot(fp) else NULL
  if (is.null(mp) || !nrow(mp)) mp <- sym_from_bioc(platform, NULL)
  if (!is.null(mp) && nrow(mp)) { mp <- as.data.frame(mp); saveRDS(mp, cache) }
  mp
}

sym_from_bioc <- function(platform, probes) {
  pkg <- switch(platform,
                "GPL570"   = "hgu133plus2.db",
                "GPL96"    = "hgu133a.db",
                "GPL97"    = "hgu133b.db",
                "GPL8300"  = "hgu95av2.db",
                "GPL571"   = "hgu133a2.db",
                NULL)
  if (is.null(pkg)) return(NULL)
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message("  安装注释包 ", pkg, " ...")
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager", repos = "https://cloud.r-project.org")
    try(BiocManager::install(pkg, ask = FALSE, update = FALSE), silent = TRUE)
  }
  if (!requireNamespace(pkg, quietly = TRUE)) return(NULL)
  mp <- try(AnnotationDbi::select(get(pkg), keys = probes,
                                  columns = "SYMBOL", keytype = "PROBEID"), silent = TRUE)
  if (inherits(mp, "try-error")) return(NULL)
  setDT(mp)[!duplicated(PROBEID) & !is.na(SYMBOL)]
}

collapse_max <- function(mat, sym) {
  keep <- !is.na(sym) & sym != "" & sym != "---"
  mat <- mat[keep, , drop = FALSE]; sym <- sym[keep]
  if (!nrow(mat)) return(NULL)
  ord <- order(sym, -rowMeans(mat, na.rm = TRUE))
  mat <- mat[ord, , drop = FALSE]; sym <- sym[ord]
  dup <- duplicated(sym)
  m <- mat[!dup, , drop = FALSE]
  rownames(m) <- sym[!dup]
  m
}

# ------------------------------------------------------- 主流程
files <- list.files(RAW, pattern = "_series_matrix\\.txt\\.gz$", full.names = TRUE)
if (!length(files)) stop("原始目录里没有 series matrix：", RAW)

for (f in files) {
  gse <- sub("_series_matrix.*$", "", basename(f))
  outf <- file.path(OUT, paste0(gse, ".rds"))
  # 增量：已构建过的跳过（想重建就用 Rscript ... refresh，
  # 或删掉 cohorts/GSEXXXX.rds）
  if (file.exists(outf) && !("refresh" %in% args)) {
    d <- readRDS(outf)
    message("=== ", gse, " 已有 rds（", nrow(d$expr), " 基因 x ", ncol(d$expr),
            " 样本，事件 ", sum(d$clin$status == 1, na.rm = TRUE), "），跳过")
    next
  }
  message("=== ", gse, " ===")
  sm <- try(read_series_matrix(f), silent = TRUE)
  if (inherits(sm, "try-error")) { message("  [跳过] 解析失败"); next }

  expr <- parse_expr(sm)
  message("  探针数：", nrow(expr), " 样本数：", ncol(expr), " 平台：", sm$platform)

  mp <- try(annot_map(sm$platform), silent = TRUE)
  if (inherits(mp, "try-error")) mp <- NULL
  if (is.null(mp) || !nrow(mp)) { message("  [跳过] 拿不到基因符号，无法匹配签名"); next }

  sym <- mp$SYMBOL[match(rownames(expr), mp$ID)]
  m <- collapse_max(expr, sym)
  if (is.null(m)) { message("  [跳过] 注释后为空"); next }
  if (max(m, na.rm = TRUE) > 50) m <- log2(pmax(m, 0) + 1)
  message("  基因数：", nrow(m))

  # ---- 临床 ----
  pc <- parse_chars(sm)
  n <- ncol(m)
  gsm <- colnames(m)
  # 互斥词：避免 grepl 误命中语义相反的字段
  EXCL <- list(
    time   = "event|status|vital|organ|dataset|type|number|protein|exon|mutat",
    status = "delay|time|month|year|day|follow|age|stage",
    stage  = "event|status",
    age    = "event|status",
    sex    = "event|status|delay"
  )
  getv <- function(map, field) {
    if (is.null(pc)) return(rep(NA_character_, n))
    idx <- match(gsm, pc$gsm)
    v <- rep(NA_character_, n)
    keys <- names(pc$fields)
    lk <- tolower(trimws(keys))
    look <- if (!is.null(map[[field]]) && length(map[[field]])) map[[field]] else GENERIC[[field]]
    bad <- EXCL[[field]]
    for (k in look) {
      cand <- which(lk == k)
      if (!length(cand)) cand <- which(grepl(k, lk, fixed = TRUE))
      if (length(bad)) cand <- cand[!grepl(bad, lk[cand])]
      if (length(cand)) {
        v <- pc$fields[[keys[cand[1]]]][idx]
        break
      }
    }
    v
  }
  map <- CLIN_MAP[[gse]]
  if (is.null(map)) map <- GENERIC
  raw_t <- getv(map, "time");  raw_s <- getv(map, "status")
  raw_st <- getv(map, "stage"); raw_sex <- getv(map, "sex")
  raw_age <- getv(map, "age");  raw_site <- getv(map, "site")

  time <- suppressWarnings(as.numeric(raw_t))
  # 注意顺序：先判否定式再判肯定式。否则 "no death" 会被 'death' 命中成 1
  # （而且 "death" 不含 "dead"，只写 dead 会漏判，所以死亡类要写全）
  st <- tolower(trimws(raw_s))
  status <- rep(NA_integer_, length(st))
  # ① 数字编码最可靠，优先：
  #    "1 (death from cancer)" -> 1     "0 (no death)" -> 0
  #    注意不能先做关键词，否则 "0 (no death)" 会被 'death' 命中成 1（实测踩过）
  lead <- suppressWarnings(as.numeric(sub("^[^0-9]*([0-9]+).*$", "\\1", st)))
  code <- !is.na(lead) & lead %in% c(0, 1)
  status[code] <- as.integer(lead[code])
  # ② 剩余走文本：先否定后肯定
  remain <- is.na(status)
  if (any(remain)) {
    s2 <- st[remain]
    neg <- grepl("\\b(no|not|non|without)\\b.*(death|dead|recur|relaps|progress|metasta|event|disease)|^alive|^living|censored|^n$|^no$", s2)
    pos <- grepl("death|dead|deceased|died|recur|relaps|progress|metasta|^yes$|^y$", s2)
    v <- rep(NA_integer_, length(s2))
    v[neg] <- 0
    v[pos & !neg] <- 1
    status[which(remain)] <- v
  }
  if (gse %in% INVERT_STATUS) status <- 1L - status

  stage <- suppressWarnings(as.numeric(sub(".*?([0-9]+).*", "\\1", tolower(raw_st))))
  clin <- data.frame(sample = gsm, time = time, status = status,
                     stage = stage,
                     age = suppressWarnings(as.numeric(raw_age)),
                     sex = raw_sex, site = raw_site,
                     stringsAsFactors = FALSE)

  n_ok <- sum(!is.na(clin$time) & !is.na(clin$status))
  n_ev <- sum(clin$status == 1, na.rm = TRUE)
  message("  可用生存：", n_ok, " / ", n, "（事件 ", n_ev, "）")
  if (n_ok < 40) message("  [警告] 生存信息不足，07 会跳过该队列。请检查字段映射：",
                         paste(head(names(pc$fields), 20), collapse = " | "))

  saveRDS(list(expr = m, clin = clin), outf)
  message("  -> ", outf)
}

message("\n=== 就绪队列 ===")
for (f in list.files(OUT, pattern = "\\.rds$", full.names = TRUE)) {
  d <- readRDS(f)
  cat(sprintf("  %-24s %6d 基因 x %4d 样本，有生存结局 %d，事件 %d\n",
              basename(f), nrow(d$expr), ncol(d$expr),
              sum(!is.na(d$clin$time) & !is.na(d$clin$status)),
              sum(d$clin$status == 1, na.rm = TRUE)))
}
