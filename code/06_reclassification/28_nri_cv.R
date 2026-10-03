#!/usr/bin/env Rscript
# =====================================================================
# 28_nri_cv.R
# NRI / IDI 的【交叉验证版】（out-of-sample），用于校正 in-sample 乐观偏倚
#
# 为什么必须做：
#   27 的 in-sample NRI 里，连【随机基因集】都能拿到 +0.13 的 NRI——
#   因为新加一个协变量在样本内必然"改善"分类。这不是 bug，是 NRI 的固有偏倚。
#   审稿人会直接问："your NRI is in-sample, show me the cross-validated one."
#
# 做法：10 折 CV。每一折在 9/10 上估 Cox 系数与基线生存，
#       在剩下 1/10 上算 5 年风险概率。汇总全部折的 out-of-sample 概率后算 NRI/IDI。
#   注意：签名风险分 .p 是外部给定的（不依赖拟合），CV 只重新估 .p 的系数 + 基线。
#        这与"签名来自文献、系数已发表"的实际使用场景一致。
#
# 规模控制：只对 Top N 签名 + 随机基因集做（全量 ×10 折太慢）
# =====================================================================

suppressPackageStartupMessages({ library(survival); library(data.table) })

COH    <- "cohorts"
LIB    <- "sigminer/output/signature_library_clean.csv"
PER    <- "audit_out_clean/per_cohort.csv"
OUT    <- "audit_out_clean"
T0     <- 60
KFOLD  <- 10
TOP_N  <- 50        # 取 C 均值最高的 N 个签名
N_RAND <- 60        # 每队列随机基因集个数
SEED   <- 20261001L

zscore <- function(m) {
  mu <- rowMeans(m, na.rm = TRUE); sd <- apply(m, 1, sd, na.rm = TRUE)
  sd[is.na(sd) | sd < 1e-8] <- 1; (m - mu) / sd
}
cindex <- function(time, status, pred) {
  ok <- !is.na(pred) & !is.na(time) & !is.na(status)
  if (sum(ok) < 20 || sum(status[ok]) < 5) return(NA_real_)
  if (length(unique(pred[ok])) < 2) return(NA_real_)
  # 注意：concordancefit 在部分 survival 版本/输入形态下会报错，必须有 fallback，
  # 否则 cindex 全 NA → oriented 返回 C=NA → 所有签名被静默跳过（曾踩过：0 个签名完成）
  cf <- try(survival::concordancefit(Surv(time[ok], status[ok]) ~ matrix(pred[ok], ncol = 1),
                                     std.err = FALSE, reverse = TRUE), silent = TRUE)
  if (inherits(cf, "try-error")) {
    cf <- try(survival::concordance(Surv(time[ok], status[ok]) ~ pred[ok],
                                    reverse = TRUE), silent = TRUE)
    if (inherits(cf, "try-error")) return(NA_real_)
    return(unname(cf$concordance))
  }
  unname(cf$concordance)
}
oriented <- function(time, status, pred) {
  c0 <- cindex(time, status, pred)
  if (is.na(c0)) return(list(C = NA_real_, p = pred))
  c1 <- cindex(time, status, -pred)
  if (!is.na(c1) && c1 > c0) list(C = c1, p = -pred) else list(C = c0, p = pred)
}
risk5 <- function(fit, newdata) {
  sf <- try(survfit(fit, newdata = newdata), silent = TRUE)
  if (inherits(sf, "try-error")) return(rep(NA_real_, nrow(newdata)))
  s <- summary(sf, times = T0, extend = TRUE)$surv
  if (is.matrix(s)) s <- s[nrow(s), ]
  as.numeric(1 - s)
}
nri_idi <- function(p_old, p_new, case) {
  ctrl <- !case; n1 <- sum(case); n0 <- sum(ctrl)
  if (n1 < 5 || n0 < 5) return(NULL)
  d <- p_new - p_old
  up <- d > 1e-12; dn <- d < -1e-12
  pu1 <- sum(up[case])/n1; pd1 <- sum(dn[case])/n1
  pu0 <- sum(up[ctrl])/n0; pd0 <- sum(dn[ctrl])/n0
  nri <- (pu1 - pd1) - (pu0 - pd0)
  vn <- (pu1 + pd1 - (pu1 - pd1)^2)/n1 + (pu0 + pd0 - (pu0 - pd0)^2)/n0
  z_nri <- if (vn > 0) nri/sqrt(vn) else NA_real_
  idi <- (mean(p_new[case]) - mean(p_old[case])) - (mean(p_new[ctrl]) - mean(p_old[ctrl]))
  vi <- var(d[case])/n1 + var(d[ctrl])/n0
  z_idi <- if (vi > 0) idi/sqrt(vi) else NA_real_
  list(n_case = n1, n_ctrl = n0, NRI = nri, NRI_z = z_nri,
       NRI_p = if (is.na(z_nri)) NA_real_ else 2*(1-pnorm(abs(z_nri))),
       IDI = idi, IDI_z = z_idi,
       IDI_p = if (is.na(z_idi)) NA_real_ else 2*(1-pnorm(abs(z_idi))))
}

# ---- 选 Top 签名
res <- fread(PER); res <- res[status_r == "ok"]
cons <- res[, .(dC = mean(C - null_median, na.rm = TRUE)), by = sig_id]
setorder(cons, -dC)
top <- head(cons$sig_id, TOP_N)

lib <- fread(LIB, sep = ",", header = TRUE, fill = TRUE,
             colClasses = "character", quote = "\"")
lib <- lib[sig_id %in% top]

# CV 核心：给定 .p（风险分向量）与临床数据，返回 out-of-sample 的 p_old / p_new
cv_probs <- function(dd, .p, fold) {
  n <- nrow(dd)
  p_old <- rep(NA_real_, n); p_new <- rep(NA_real_, n)
  f0 <- as.formula(paste("Surv(.t, .s) ~", paste(vars_all, collapse = " + ")))
  f1 <- as.formula(paste("Surv(.t, .s) ~ .p +", paste(vars_all, collapse = " + ")))
  for (k in seq_len(KFOLD)) {
    te <- which(fold == k); tr <- which(fold != k)
    if (length(te) == 0 || length(tr) < 20) next
    dtr <- dd[tr, , drop = FALSE]; dte <- dd[te, , drop = FALSE]
    fit0 <- try(coxph(f0, data = dtr), silent = TRUE)
    dtr$.p <- .p[tr]; dte$.p <- .p[te]
    fit1 <- try(coxph(f1, data = dtr), silent = TRUE)
    if (!inherits(fit0, "try-error")) p_old[te] <- risk5(fit0, dte)
    if (!inherits(fit1, "try-error")) p_new[te] <- risk5(fit1, dte)
  }
  list(p_old = p_old, p_new = p_new)
}

ROWS <- list(); RND <- list()
for (cf in list.files(COH, pattern = "\\.rds$", full.names = TRUE)) {
  coh <- sub("\\.rds$", "", basename(cf))
  cat(sprintf("\n===== %s =====\n", coh))
  d <- readRDS(cf)
  z <- zscore(d$expr)
  cl <- d$clin; cl <- cl[match(colnames(z), cl$sample), , drop = FALSE]
  time <- suppressWarnings(as.numeric(cl$time)); status <- suppressWarnings(as.integer(cl$status))
  keep <- !is.na(time) & !is.na(status)
  z <- z[, keep, drop = FALSE]; time <- time[keep]; status <- status[keep]
  cl2 <- cl[keep, , drop = FALSE]

  vars_all <- intersect(c("stage", "age", "sex", "site"), colnames(cl2))
  vars_all <- vars_all[sapply(cl2[, vars_all, drop = FALSE],
                              function(x) length(unique(x[!is.na(x)])) > 1)]
  dd <- cl2[, vars_all, drop = FALSE]; dd$.t <- time; dd$.s <- status
  okK <- complete.cases(dd); dd <- dd[okK, , drop = FALSE]
  tt <- dd$.t; ss <- dd$.s; zz <- z[, okK, drop = FALSE]
  rm(z); invisible(gc(FALSE))

  case <- (tt <= T0 & ss == 1); ctrl <- (tt > T0); usable <- case | ctrl
  cat(sprintf("  n=%d case=%d ctrl=%d 可用=%d\n", length(tt), sum(case), sum(ctrl), sum(usable)))
  if (sum(case) < 10 || sum(ctrl) < 10) { rm(d); next }

  set.seed(SEED)
  fold <- sample(rep_len(seq_len(KFOLD), length.out = length(tt)))
  # 保证每折都有事件
  fold <- as.integer(fold)

  for (i in seq_len(nrow(lib))) {
    genes <- trimws(unlist(strsplit(lib$genes[i], "[;,]"))); genes <- genes[genes != ""]
    g <- intersect(genes, rownames(zz))
    if (length(g) < 2) next
    rr <- oriented(tt, ss, as.numeric(colMeans(zz[g, , drop = FALSE])))
    if (is.na(rr$C)) next
    cp <- cv_probs(dd, rr$p, fold)
    okp <- usable & !is.na(cp$p_old) & !is.na(cp$p_new)
    if (sum(case[okp]) < 5 || sum(ctrl[okp]) < 5) next
    r <- nri_idi(cp$p_old[okp], cp$p_new[okp], case[okp])
    if (is.null(r)) next
    ROWS[[length(ROWS)+1]] <- data.table(
      cohort = coh, sig_id = lib$sig_id[i], n_gene = length(genes),
      n_measured = length(g), n_case = r$n_case, n_ctrl = r$n_ctrl,
      NRI_cv = round(r$NRI, 4), NRI_cv_z = round(r$NRI_z, 4), NRI_cv_p = signif(r$NRI_p, 4),
      IDI_cv = round(r$IDI, 5), IDI_cv_z = round(r$IDI_z, 4), IDI_cv_p = signif(r$IDI_p, 4))
  }
  cat(sprintf("  %d 个签名完成 CV-NRI/IDI\n", length(ROWS)))

  set.seed(SEED + 1L)
  for (b in seq_len(N_RAND)) {
    ng <- sample(c(3L, 5L, 8L, 12L, 20L), 1L)
    gs <- sample(rownames(zz), ng)
    rb <- oriented(tt, ss, as.numeric(colMeans(zz[gs, , drop = FALSE])))
    if (is.na(rb$C)) next
    cp <- cv_probs(dd, rb$p, fold)
    okp <- usable & !is.na(cp$p_old) & !is.na(cp$p_new)
    if (sum(case[okp]) < 5 || sum(ctrl[okp]) < 5) next
    r <- nri_idi(cp$p_old[okp], cp$p_new[okp], case[okp])
    if (is.null(r)) next
    RND[[length(RND)+1]] <- data.table(cohort = coh, rep = b, n_gene = ng,
                                       NRI_cv = r$NRI, IDI_cv = r$IDI,
                                       NRI_cv_p = r$NRI_p, IDI_cv_p = r$IDI_p)
  }
  rm(d, zz, dd); invisible(gc(FALSE))
}

out <- rbindlist(ROWS); fwrite(out, file.path(OUT, "NRI_IDI_CV_签名.csv"))
rnd <- if (length(RND)) rbindlist(RND) else NULL
if (!is.null(rnd)) fwrite(rnd, file.path(OUT, "NRI_IDI_CV_随机基线.csv"))

cat("\n\n================ CV 版汇总 ================\n")
cat(sprintf("已发表签名（Top%d）CV-NRI 中位 %.4f  为正 %.1f%%  p<0.05 占 %.1f%%\n",
            TOP_N, out[, median(NRI_cv, na.rm = TRUE)],
            100*out[, mean(NRI_cv > 0, na.rm = TRUE)],
            100*out[, mean(NRI_cv_p < 0.05, na.rm = TRUE)]))
cat(sprintf("已发表签名（Top%d）CV-IDI 中位 %.5f  为正 %.1f%%\n",
            TOP_N, out[, median(IDI_cv, na.rm = TRUE)],
            100*out[, mean(IDI_cv > 0, na.rm = TRUE)]))
if (!is.null(rnd)) {
  cat(sprintf("随机基因集       CV-NRI 中位 %.4f  为正 %.1f%%   ← 应≈0\n",
              rnd[, median(NRI_cv, na.rm = TRUE)], 100*rnd[, mean(NRI_cv > 0, na.rm = TRUE)]))
  cat(sprintf("随机基因集       CV-IDI 中位 %.5f  为正 %.1f%%\n",
              rnd[, median(IDI_cv, na.rm = TRUE)], 100*rnd[, mean(IDI_cv > 0, na.rm = TRUE)]))
  cat(sprintf("\n净效应（签名 − 随机）: CV-NRI %.4f   CV-IDI %.5f\n",
              out[, median(NRI_cv, na.rm = TRUE)] - rnd[, median(NRI_cv, na.rm = TRUE)],
              out[, median(IDI_cv, na.rm = TRUE)] - rnd[, median(IDI_cv, na.rm = TRUE)]))
}
agg <- out[, .(n = .N, NRI中位 = round(median(NRI_cv, na.rm = TRUE), 4),
               IDI中位 = round(median(IDI_cv, na.rm = TRUE), 5)), by = cohort]
print(agg)
fwrite(agg, file.path(OUT, "NRI_IDI_CV_队列汇总.csv"))
cat("\n已存：NRI_IDI_CV_签名.csv / NRI_IDI_CV_随机基线.csv / NRI_IDI_CV_队列汇总.csv\n")
