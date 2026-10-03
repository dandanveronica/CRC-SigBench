#!/usr/bin/env Rscript
# =====================================================================
# 07c_clinical_benchmark.R
# 关键补充分析：签名在【已有临床变量】之上的增量价值（有基线的 ΔC-index）
#
# 为什么必须做这一步：
#   07 里算的 delta_C 只给了「加签名后 C 提升了多少」，
#   但读者/审稿人真正要问的是「临床变量本身已经能到多少 C、
#   加上签名之后是 0.66 -> 0.67 还是 0.52 -> 0.60」。
#   这两句话的含金量完全不同，必须给出基线模型的绝对 C 值。
#
# 同时加两个对照，让结论站得住：
#   (1) 随机等大基因集：这是本研究的零假设锚点
#   (2) 仅年龄/仅分期：看签名到底是被 stage 代理了多少
#
# 对每个队列输出：
#   C(stage+age+sex+site)  →  加入 Top 签名后的 C  →  绝对增量
# 并给出 permutation-free 的点估计（两侧波动用签名层面的分布展示）
# =====================================================================

suppressPackageStartupMessages({ library(survival); library(data.table) })

COH <- "cohorts"
PER <- "audit_out_clean/per_cohort.csv"
LIB <- "sigminer/output/signature_library_clean.csv"
OUT <- "audit_out_clean"
TOP_N <- 20

zscore <- function(m) {
  mu <- rowMeans(m, na.rm = TRUE); sd <- apply(m, 1, sd, na.rm = TRUE)
  sd[is.na(sd) | sd < 1e-8] <- 1; (m - mu) / sd
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
oriented <- function(time, status, pred) {
  c0 <- cindex(time, status, pred)
  if (is.na(c0)) return(list(C = NA_real_, p = pred))
  c1 <- cindex(time, status, -pred)
  if (!is.na(c1) && c1 > c0) list(C = c1, p = -pred) else list(C = c0, p = pred)
}

# ------------------------------------------------ 挑最能打的一批签名
res <- fread(PER); res <- res[status_r == "ok"]
cons <- res[, .(beat = sum(C > null_median, na.rm = TRUE),
                dC   = mean(C - null_median, na.rm = TRUE)), by = sig_id]
setorder(cons, -beat, -dC)
top <- head(cons$sig_id, TOP_N)
cat("纳入 Benchmark 的 Top 签名：", paste(top, collapse = ", "), "\n\n")

lib <- fread(LIB, sep = ",", header = TRUE, fill = TRUE,
             colClasses = "character", quote = "\"")
lib <- lib[sig_id %in% top]

ROWS <- list()
for (cf in list.files(COH, pattern = "\\.rds$", full.names = TRUE)) {
  coh <- sub("\\.rds$", "", basename(cf))
  d <- readRDS(cf)
  z  <- zscore(d$expr)
  cl <- d$clin; cl <- cl[match(colnames(z), cl$sample), , drop = FALSE]
  time <- suppressWarnings(as.numeric(cl$time))
  status <- suppressWarnings(as.integer(cl$status))
  keep <- !is.na(time) & !is.na(status)
  z <- z[, keep, drop = FALSE]; time <- time[keep]; status <- status[keep]
  cl2 <- cl[keep, , drop = FALSE]

  vars_all <- intersect(c("stage", "age", "sex", "site"), colnames(cl2))
  vars_all <- vars_all[sapply(cl2[, vars_all, drop = FALSE],
                              function(x) length(unique(x[!is.na(x)])) > 1)]
  dd <- cl2[, vars_all, drop = FALSE]
  dd$.t <- time; dd$.s <- status
  okK <- complete.cases(dd)
  dd <- dd[okK, , drop = FALSE]
  tt <- dd$.t; ss <- dd$.s
  zz <- z[, okK, drop = FALSE]
  cat(sprintf("[%s] n=%d 事件=%d  可用临床变量：%s\n",
              coh, length(tt), sum(ss), paste(vars_all, collapse = "+")))

  C_of <- function(v) {
    f <- as.formula(paste("Surv(.t, .s) ~", v))
    fit <- try(coxph(f, data = dd), silent = TRUE)
    if (inherits(fit, "try-error")) return(NA_real_)
    unname(cindex(tt, ss, predict(fit, type = "lp")))
  }
  base_full <- C_of(paste(vars_all, collapse = " + "))
  base_stage <- if ("stage" %in% vars_all) C_of("stage") else NA_real_
  base_age   <- if ("age"   %in% vars_all) C_of("age")   else NA_real_
  cat(sprintf("   基线 C(stage+age+sex+site)=%.3f  C(stage)=%.3f  C(age)=%.3f\n",
              base_full, base_stage, base_age))

  # ---- 逐个 Top 签名
  for (i in seq_len(nrow(lib))) {
    genes <- trimws(unlist(strsplit(lib$genes[i], "[;,]")))
    genes <- genes[genes != ""]
    g <- intersect(genes, rownames(zz))
    if (length(g) < 2) next
    rr <- oriented(tt, ss, as.numeric(colMeans(zz[g, , drop = FALSE])))
    dd2 <- dd; dd2$.p <- rr$p
    f_txt <- paste("Surv(.t, .s) ~ .p +", paste(vars_all, collapse = " + "))
    fit <- try(coxph(as.formula(f_txt), data = dd2), silent = TRUE)
    if (inherits(fit, "try-error")) next
    C_full <- cindex(tt, ss, predict(fit, type = "lp"))
    # 随机对照
    n_rand <- 200
    set.seed(1234L)
    dc_rand <- numeric(n_rand)
    for (b in seq_len(n_rand)) {
      gs <- sample(rownames(zz), length(g))
      rb <- oriented(tt, ss, as.numeric(colMeans(zz[gs, , drop = FALSE])))
      dd3 <- dd; dd3$.p <- rb$p
      fb <- try(coxph(as.formula(f_txt), data = dd3), silent = TRUE)
      if (inherits(fb, "try-error")) next
      dc_rand[b] <- cindex(tt, ss, predict(fb, type = "lp")) - base_full
    }
    dc_rand <- dc_rand[!is.na(dc_rand)]
    ROWS[[length(ROWS) + 1]] <- data.table(
      cohort = coh, sig_id = lib$sig_id[i], n_gene = length(genes),
      n_measured = length(g),
      C_clinical_only = round(base_full, 4),
      C_clinical_plus_sig = round(C_full, 4),
      abs_gain = round(C_full - base_full, 4),
      rand_gain_med = if (length(dc_rand)) round(median(dc_rand), 4) else NA_real_,
      rand_gain_p95 = if (length(dc_rand)) round(quantile(dc_rand, .95), 4) else NA_real_,
      beat_random = if (length(dc_rand)) as.integer(C_full - base_full > quantile(dc_rand, .95)) else NA
    )
  }
  rm(d); invisible(gc(FALSE))
}

out <- rbindlist(ROWS)
fwrite(out, file.path(OUT, "临床增量基准.csv"))

cat("\n=== 汇总：Top 签名在临床模型之上的绝对增益 ===\n")
agg <- out[, .(n = .N,
               C临床基线 = round(mean(C_clinical_only, na.rm = TRUE), 3),
               C加签名后 = round(mean(C_clinical_plus_sig, na.rm = TRUE), 3),
               增益中位 = round(median(abs_gain, na.rm = TRUE), 4),
               增益四分位 = paste0(round(quantile(abs_gain, .25, na.rm = TRUE), 4),
                                   " ~ ",
                                   round(quantile(abs_gain, .75, na.rm = TRUE), 4)),
               超过随机95分位 = paste0(sum(beat_random == 1, na.rm = TRUE), "/", .N)
), by = cohort]
print(agg)
cat("\n整体：临床基线 C 中位 ", round(out[, median(C_clinical_only, na.rm = TRUE)], 3),
    " -> 加签名后 ", round(out[, median(C_clinical_plus_sig, na.rm = TRUE)], 3),
    "，绝对增益中位 ", round(out[, median(abs_gain, na.rm = TRUE)], 4), "\n", sep = "")
fwrite(agg, file.path(OUT, "临床增量_队列汇总.csv"))
