#!/usr/bin/env Rscript
# =====================================================================
# 22_recompute_train_C.R
# 目的：复算协议合规签名在【训练集 GSE39582】上的 C-index（主稿 0.801），
#       并把结果落盘为 CSV，补齐该数字此前只存在于控制台输出的问题。
#
# 说明：
#   * 训练逻辑与 21_part3_compliant.R 的 ① 部分【逐行一致】，
#     包括同一个 set.seed(20261001)，以保证 LASSO/10 折 CV 可复现。
#   * 队列 rds 由 01_build_cohorts/06b_build_cohorts.R 构建。
#
# 用法：Rscript 22_recompute_train_C.R [cohorts目录] [输出csv]
# =====================================================================
suppressPackageStartupMessages({
  library(survival); library(data.table); library(glmnet)
})

args   <- commandArgs(trailingOnly = TRUE)
COH    <- ifelse(length(args) >= 1, args[1], "cohorts")
OUTCSV <- ifelse(length(args) >= 2, args[2], "results/training_performance.csv")
TRAIN  <- "GSE39582"
HORIZON <- 5

zscore <- function(m) {
  mu <- rowMeans(m, na.rm = TRUE); sd <- apply(m, 1, sd, na.rm = TRUE)
  sd[is.na(sd) | sd < 1e-8] <- 1
  (m - mu) / sd
}

oriented_c <- function(sv, pred) {
  xb <- matrix(as.numeric(pred), ncol = 1)
  cf <- try(concordancefit(sv, xb, reverse = TRUE, std.err = FALSE), silent = TRUE)
  if (inherits(cf, "try-error")) return(NA_real_)
  c0 <- cf$concordance
  if (is.na(c0)) return(NA_real_)
  if (c0 < 0.5) 1 - c0 else c0
}

auc5 <- function(time, status, pred, h = HORIZON * 12) {
  ok <- !is.na(pred) & !is.na(time) & !is.na(status)
  t <- time[ok]; s <- status[ok]; p <- as.numeric(pred[ok])
  ev <- s == 1 & t <= h
  no <- t > h
  if (sum(ev) < 5 || sum(no) < 5) return(NA_real_)
  y <- c(rep(1, sum(ev)), rep(0, sum(no)))
  x <- c(p[ev], p[no])
  r <- rank(x)
  n1 <- sum(y); n0 <- length(y) - n1
  (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

calib_slope <- function(time, status, pred) {
  dd <- data.frame(t = as.numeric(time), s = as.integer(status), p = as.numeric(pred))
  dd <- dd[is.finite(dd$t) & !is.na(dd$s) & is.finite(dd$p) & dd$t > 0, ]
  if (nrow(dd) < 30 || sum(dd$s) < 5) return(NA_real_)
  f <- try(coxph(Surv(t, s) ~ p, data = dd), silent = TRUE)
  if (inherits(f, "try-error")) return(NA_real_)
  unname(coef(f))
}

cox_score <- function(z, time, status) {
  o <- order(time); tm <- time[o]; s <- status[o]; X <- z[, o, drop = FALSE]
  X <- X[, tm > 0, drop = FALSE]; s <- s[tm > 0]
  n <- length(s); p <- nrow(X)
  w <- rev(cumsum(rev(s)))
  w[w <= 0] <- 1
  SX <- t(apply(X, 1, cumsum))
  SXr <- SX[, n] - SX + X
  SXb <- t(apply(X * rep(s, each = p), 1, function(v) rev(cumsum(rev(v)))))
  U <- rowSums((X - t(t(SXr) / w)) * rep(s, each = p))
  V <- rowSums((X^2 - t(t(SXr^2) / w^2)) * rep(s, each = p))
  num <- U; den <- sqrt(pmax(V, 1e-12))
  zz <- num / den; zz[!is.finite(zz)] <- 0
  setNames(zz, rownames(X))
}

# ---------------- 载入训练队列（与 21 完全一致的过滤） ----------------
f <- file.path(COH, paste0(TRAIN, ".rds"))
if (!file.exists(f)) stop("缺少队列文件：", f)
d  <- readRDS(f)
z  <- zscore(d$expr)
cl <- d$clin; cl <- cl[match(colnames(z), cl$sample), , drop = FALSE]
tm <- suppressWarnings(as.numeric(cl$time)); st <- suppressWarnings(as.integer(cl$status))
kp <- !is.na(tm) & !is.na(st) & tm > 0
z <- z[, kp, drop = FALSE]; tm <- tm[kp]; st <- st[kp]
cat(sprintf("训练集 %s：n=%d 事件=%d\n", TRAIN, length(tm), sum(st)))

zt <- z; tmt <- tm; stt <- st
sv <- Surv(tmt, stt)

# ① 单因素 Cox 预筛
sc <- cox_score(zt, tmt, stt)
presel <- names(sort(abs(sc), decreasing = TRUE))[seq_len(min(2000, length(sc)))]
cat("预筛后候选基因：", length(presel), "\n")

# ② LASSO-Cox（与 21 同一个随机种子）
x <- t(zt[presel, , drop = FALSE])
fit <- glmnet(x, sv, family = "cox", alpha = 1, nlambda = 200,
              standardize = TRUE, maxit = 10000)
set.seed(20261001)
cv  <- cv.glmnet(x, sv, family = "cox", alpha = 1, nfolds = 10,
                 standardize = TRUE, maxit = 10000)
co <- as.matrix(coef(fit, s = cv$lambda.min))
co <- co[co[, 1] != 0, , drop = FALSE]
genes <- rownames(co); genes <- genes[genes != "(Intercept)"]
beta <- co[genes, 1]
cat("LASSO 选中基因数：", length(genes), "（lambda.min =", signif(cv$lambda.min, 4), "）\n")

# ③ 多因素 Cox 定系数
dd <- data.frame(t = tmt, s = stt, t(zt[genes, , drop = FALSE]))
names(dd) <- c("t", "s", make.names(genes))
mfit <- try(coxph(as.formula(paste("Surv(t, s) ~", paste(make.names(genes), collapse = " + "))),
                  data = dd), silent = TRUE)
if (inherits(mfit, "try-error")) {
  cat("多因素 Cox 失败，回退 LASSO 系数\n"); beta_final <- beta
} else {
  beta_final <- coef(mfit); names(beta_final) <- genes
}

# ④ 训练集风险分与性能
risk_fun <- function(z, g, b) as.numeric(t(z[g, , drop = FALSE]) %*% b)
r_tr <- risk_fun(zt, genes, beta_final)
cutoff <- median(r_tr)
C_tr   <- oriented_c(sv, r_tr)
A5_tr  <- auc5(tmt, stt, r_tr)
CS_tr  <- calib_slope(tmt, stt, r_tr)

cat(sprintf("\n训练集 C-index  = %.4f\n", C_tr))
cat(sprintf("训练集 5年 AUC  = %.4f\n", A5_tr))
cat(sprintf("训练集校准斜率  = %.4f\n", CS_tr))
cat(sprintf("锁死截点(中位)  = %.5f\n", cutoff))

# ⑤ 落盘
out <- data.table(
  cohort        = TRAIN,
  n             = length(tmt),
  events        = sum(stt),
  n_gene        = length(genes),
  lambda_min    = cv$lambda.min,
  C_training    = C_tr,
  auc5_training = A5_tr,
  calib_slope_training = CS_tr,
  locked_cutoff = cutoff,
  seed          = 20261001L,
  script        = "22_recompute_train_C.R"
)
dir.create(dirname(OUTCSV), showWarnings = FALSE, recursive = TRUE)
fwrite(out, OUTCSV, encoding = "UTF-8")
cat("\n已写出：", OUTCSV, "\n")
