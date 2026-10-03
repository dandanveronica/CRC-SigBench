#!/usr/bin/env Rscript
# =====================================================================
# 21_part3_compliant.R
# Part 3：用归因分析得出的"最小合规流程"重做一个签名，再扔进同一套审判
#
# 从 Part 1/2 得到的关键规律
#   ① 决定成败的是【发现队列的事件数】，不是样本数
#      （事件 ≥70 的队列做出的签名能迁移；<30 的等于白做）
#   ② 方法学"合规标记"（外部验证/多因素校正/公开系数）与真实表现无关
#   ③ 原文声称值与独立重算值相差 +0.16（AUC5）/+0.19（C-index）
#
# 因此"最小合规流程"定义为：
#   A. 在【事件数最多】的队列上训练（本例 GSE39582，194 事件）
#   B. 单因素 Cox 预筛 -> LASSO-Cox 选基因 -> 多因素 Cox 定系数
#   C. 系数与【风险分截点】在训练集上一次确定后【锁死】，验证时不再调整
#   D. 报告 C-index + 5 年 AUC + 校准斜率 + DCA 净获益
#
# 关键对照（本脚本最有价值的部分）
#   合规版：验证队列用训练集锁死的截点
#   违规版：验证队列【就地重新优化】截点（模拟常见的"随行就市"做法）
#   -> 两者之差 = 不锁死截点能虚增多少，这是可量化的报告偏倚上界
# =====================================================================
suppressPackageStartupMessages({
  library(survival); library(data.table); library(glmnet)
})

BASE <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)))
if (!length(BASE)) BASE <- getwd()
setwd(BASE)
COH <- "cohorts"; OUT <- "audit_out_clean"
dir.create(OUT, showWarnings = FALSE)

TRAIN   <- "GSE39582"     # 事件数最多的队列
HORIZON <- 5              # 5 年
B_NULL  <- 3000           # 随机基因集重复次数

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
  # 时间依赖 AUC（简化：5 年内的事件 vs 非事件，剔除 5 年前删失）
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
  # 校准斜率：Cox 模型中风险分的回归系数，理想值 = 1
  dd <- data.frame(t = as.numeric(time), s = as.integer(status), p = as.numeric(pred))
  dd <- dd[is.finite(dd$t) & !is.na(dd$s) & is.finite(dd$p) & dd$t > 0, ]
  if (nrow(dd) < 30 || sum(dd$s) < 5) return(NA_real_)
  f <- try(coxph(Surv(t, s) ~ p, data = dd), silent = TRUE)
  if (inherits(f, "try-error")) return(NA_real_)
  unname(coef(f))
}

net_benefit <- function(time, status, pred, pt, h = HORIZON * 12) {
  # DCA 净获益（简化：以 h 时刻的可判定样本近似，剔除 h 前删失）
  ok <- !is.na(pred) & !is.na(time) & !is.na(status)
  t <- time[ok]; s <- status[ok]; p <- as.numeric(pred[ok])
  ev <- s == 1 & t <= h
  no <- t > h
  if (sum(ev) < 5 || sum(no) < 5) return(NA_real_)
  # 把风险分转成 0-1 的"预测概率"：用经验 logit，避免量纲问题
  y <- c(rep(1, sum(ev)), rep(0, sum(no)))
  x <- c(p[ev], p[no])
  fit <- suppressWarnings(try(glm(y ~ x, family = binomial()), silent = TRUE))
  if (inherits(fit, "try-error")) return(NA_real_)
  pr <- as.numeric(predict(fit, type = "response"))
  n <- length(y)
  tp <- sum(pr >= pt & y == 1); fp <- sum(pr >= pt & y == 0)
  tp / n - (fp / n) * (pt / (1 - pt))
}

# ---------------- 载入队列 ----------------
files <- list.files(COH, pattern = "\\.rds$", full.names = TRUE)
nm <- sub("\\.rds$", "", basename(files))
cat("可用队列：", paste(nm, collapse = " "), "\n\n")

CL <- list()
for (f in files) {
  d <- readRDS(f)
  z <- zscore(d$expr)
  cl <- d$clin; cl <- cl[match(colnames(z), cl$sample), , drop = FALSE]
  tm <- suppressWarnings(as.numeric(cl$time)); st <- suppressWarnings(as.integer(cl$status))
  kp <- !is.na(tm) & !is.na(st) & tm > 0        # 非正时间必须剔除
  z <- z[, kp, drop = FALSE]; tm <- tm[kp]; st <- st[kp]
  CL[[sub("\\.rds$", "", basename(f))]] <- list(z = z, time = tm, status = st)
  cat(sprintf("  %-11s n=%3d 事件 %3d\n", sub("\\.rds$", "", basename(f)), length(tm), sum(st)))
}

# ================= ① 训练（合规流程） =================
cat("\n=== ① 在", TRAIN, "上训练 ===\n")
tr <- CL[[TRAIN]]
zt <- tr$z; tmt <- tr$time; stt <- tr$status
sv <- Surv(tmt, stt)

# A. 单因素 Cox 预筛（向量化 score 检验）
cox_score <- function(z, time, status) {
  o <- order(time); tm <- time[o]; s <- status[o]; X <- z[, o, drop = FALSE]
  X <- X[, tm > 0, drop = FALSE]; s <- s[tm > 0]
  n <- length(s); p <- nrow(X)
  w <- rev(cumsum(rev(s)))                     # 风险集大小
  w[w <= 0] <- 1
  SX <- t(apply(X, 1, cumsum))
  SXr <- SX[, n] - SX + X                      # 风险集内的和
  SXb <- t(apply(X * rep(s, each = p), 1, function(v) rev(cumsum(rev(v)))))
  # 注意：必须用 rep(s, each = p) 做列缩放。
  # 写成 X * s 会被 R 按【列优先存储顺序】循环取值，结果全错。
  U <- rowSums((X - t(t(SXr) / w)) * rep(s, each = p))
  V <- rowSums((X^2 - t(t(SXr^2) / w^2)) * rep(s, each = p))
  num <- U; den <- sqrt(pmax(V, 1e-12))
  zz <- num / den; zz[!is.finite(zz)] <- 0
  setNames(zz, rownames(X))
}
sc <- cox_score(zt, tmt, stt)
presel <- names(sort(abs(sc), decreasing = TRUE))[seq_len(min(2000, length(sc)))]
cat("预筛后候选基因：", length(presel), "\n")

# B. LASSO-Cox
x <- t(zt[presel, , drop = FALSE])
fit <- glmnet(x, sv, family = "cox", alpha = 1, nlambda = 200, standardize = TRUE, maxit = 10000)
# 用 10 折 CV 选 lambda
set.seed(20261001)
cv <- cv.glmnet(x, sv, family = "cox", alpha = 1, nfolds = 10, standardize = TRUE, maxit = 10000)
co <- as.matrix(coef(fit, s = cv$lambda.min))
co <- co[co[, 1] != 0, , drop = FALSE]
genes <- rownames(co); genes <- genes[genes != "(Intercept)"]
beta <- co[genes, 1]
cat("LASSO 选中基因数：", length(genes), "（lambda.min =", signif(cv$lambda.min, 4), "）\n")
cat("基因：", paste(genes, collapse = ", "), "\n")

# C. 多因素 Cox 定系数（锁死）
dd <- data.frame(t = tmt, s = stt, t(zt[genes, , drop = FALSE]))
names(dd) <- c("t", "s", make.names(genes))
mfit <- try(coxph(as.formula(paste("Surv(t, s) ~", paste(make.names(genes), collapse = " + "))),
                  data = dd), silent = TRUE)
if (inherits(mfit, "try-error")) {
  cat("多因素 Cox 失败，回退到 LASSO 系数\n"); beta_final <- beta
} else {
  beta_final <- coef(mfit); names(beta_final) <- genes
}

# 训练集风险分 + 锁死截点
risk_fun <- function(z, g, b) as.numeric(t(z[g, , drop = FALSE]) %*% b)
r_tr <- risk_fun(zt, genes, beta_final)
cutoff <- median(r_tr)                        # ★ 锁死：训练集风险分中位数
cat("\n训练集风险分：中位", signif(cutoff, 5), "（此截点将锁死用于所有验证队列）\n")
cat("训练集 C-index：", signif(oriented_c(sv, r_tr), 4), "\n")
cat("训练集 5 年 AUC：", signif(auc5(tmt, stt, r_tr), 4), "\n")

# ================= ② 在其余队列验证 =================
cat("\n=== ② 独立验证（系数与截点全部锁死）===\n")
tests <- setdiff(names(CL), TRAIN)
res <- rbindlist(lapply(tests, function(ch) {
  d <- CL[[ch]]
  g_here <- intersect(genes, rownames(d$z))
  miss <- setdiff(genes, rownames(d$z))
  r <- risk_fun(d$z, g_here, beta_final[g_here])
  svh <- Surv(d$time, d$status)
  C <- oriented_c(svh, r)
  # 合规版：训练集锁死截点
  hi_lock <- r > cutoff
  hr_lock <- NA_real_; pv_lock <- NA_real_
  if (length(unique(hi_lock)) == 2 && sum(hi_lock) >= 5 && sum(!hi_lock) >= 5) {
    f <- try(coxph(svh ~ hi_lock), silent = TRUE)
    if (!inherits(f, "try-error")) {
      hr_lock <- unname(exp(coef(f)))
      pv_lock <- summary(f)$coefficients[1, 5]
    }
  }
  # 违规版：就地用验证队列中位数重新定截点
  cut_local <- median(r)
  hi_local <- r > cut_local
  hr_local <- NA_real_; pv_local <- NA_real_
  if (length(unique(hi_local)) == 2 && sum(hi_local) >= 5 && sum(!hi_local) >= 5) {
    f2 <- try(coxph(svh ~ hi_local), silent = TRUE)
    if (!inherits(f2, "try-error")) {
      hr_local <- unname(exp(coef(f2)))
      pv_local <- summary(f2)$coefficients[1, 5]
    }
  }
  # 违规版 2：就地搜索最优截点（最大化 log-rank）
  if (length(unique(r)) > 10) {
    qs <- quantile(r, probs = seq(.2, .8, by = .05))
    best <- NA_real_; bestp <- 1
    for (q in qs) {
      h <- r > q
      if (sum(h) < 5 || sum(!h) < 5) next
      f3 <- try(coxph(svh ~ h), silent = TRUE)
      if (inherits(f3, "try-error")) next
      pp <- summary(f3)$coefficients[1, 5]
      if (is.finite(pp) && pp < bestp) { bestp <- pp; best <- unname(exp(coef(f3))) }
    }
  } else { best <- NA_real_; bestp <- NA_real_ }
  data.table(cohort = ch, n = length(r), events = sum(d$status),
             n_gene = length(g_here), missing_gene = length(miss),
             C = C, auc5 = auc5(d$time, d$status, r),
             calib_slope = calib_slope(d$time, d$status, r),
             HR_locked = hr_lock, p_locked = pv_lock,
             HR_local_median = hr_local, p_local_median = pv_local,
             HR_bestcut = best, p_bestcut = bestp,
             nb_locked = net_benefit(d$time, d$status, r, 0.2))
}))
print(res, digits = 4)

# ================= ③ 随机基因集基线（同队列、同基因数） =================
cat("\n=== ③ 随机基因集基线（同队列、同基因数、等权重）===\n")
nulls <- rbindlist(lapply(tests, function(ch) {
  d <- CL[[ch]]; nk <- seq_len(nrow(d$z))
  ng <- length(intersect(genes, rownames(d$z)))
  svh <- Surv(d$time, d$status); xb <- matrix(0, length(d$time), 1)
  set.seed(4242L + sum(utf8ToInt(ch)))
  v <- numeric(B_NULL)
  for (b in seq_len(B_NULL)) {
    xb[, 1] <- colMeans(d$z[sample(nk, ng), , drop = FALSE])
    v[b] <- oriented_c(svh, xb)
  }
  data.table(cohort = ch, null_median = median(v, na.rm = TRUE),
             null_p95 = quantile(v, .95, na.rm = TRUE))
}))
res <- merge(res, nulls, by = "cohort")
res[, dC := C - null_median]

# ================= ④ 与 320 个已发表签名比较 =================
cat("\n=== ④ 与已发表签名的头对头 ===\n")
per <- tryCatch(fread(file.path(OUT, "per_cohort.csv"), data.table = FALSE), error = function(e) NULL)
comp <- NULL
if (!is.null(per)) {
  per <- per[per$status_r == "ok", ]
  pub <- per[per$cohort %in% tests, ]
  cat(sprintf("已发表签名在同样 %d 个验证队列上：中位 C = %.4f，优于随机 %.1f%%\n",
              length(tests), median(pub$C, na.rm = TRUE),
              100 * mean(pub$C > pub$null_median, na.rm = TRUE)))
  cat(sprintf("本合规签名：                    中位 C = %.4f，优于随机 %.1f%%\n",
              median(res$C, na.rm = TRUE), 100 * mean(res$C > res$null_median, na.rm = TRUE)))
  w <- try(wilcox.test(res$C, pub$C), silent = TRUE)
  if (!inherits(w, "try-error"))
    cat(sprintf("配对比较 Wilcoxon p = %.3g\n", w$p.value))
  comp <- data.table(对比 = c("合规签名", "已发表签名（320 个）"),
                     中位C = c(median(res$C, na.rm = TRUE), median(pub$C, na.rm = TRUE)),
                     中位ΔC = c(median(res$dC, na.rm = TRUE), median(pub$delta_C, na.rm = TRUE)),
                     优于随机比例 = c(mean(res$C > res$null_median, na.rm = TRUE),
                                 mean(pub$C > pub$null_median, na.rm = TRUE)))
  print(comp, digits = 4)
}

# ================= ⑤ 输出 =================
fwrite(res, file.path(OUT, "Part3_合规签名_验证.csv"), encoding = "UTF-8")
if (!is.null(comp)) fwrite(comp, file.path(OUT, "Part3_对比.csv"), encoding = "UTF-8")
sig_tbl <- data.table(gene = genes, coef = as.numeric(beta_final))
fwrite(sig_tbl, file.path(OUT, "Part3_签名基因与系数.csv"), encoding = "UTF-8")
cat("\n已写出：\n  Part3_合规签名_验证.csv\n  Part3_签名基因与系数.csv\n")
cat("\n★ 锁死截点 vs 就地重优化截点（违规做法能虚增多少）：\n")
sel <- res[, .(cohort, HR_locked, HR_local_median, HR_bestcut, p_locked, p_bestcut)]
print(sel, digits = 4)
