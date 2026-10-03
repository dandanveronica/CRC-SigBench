#!/usr/bin/env Rscript
# =====================================================================
# 07_audit_signatures.R
# 批量重算已发表的结直肠癌预后签名 —— 系统评估主脚本
#
# 核心设计（对应研究方案 Part 1 / Part 2）
#   1. 每个队列内做 z-score，避免跨平台绝对值不可比
#   2. 风险分三种模式：equal（等权重）/ reported（原系数）/ refit（重拟合）
#   3. 指标：C-index、时间依赖 AUC、校准斜率、ΔC-index、NRI、IDI、DCA 净获益
#   4. 随机基因集零假设对照 —— 这是本研究最关键的判据
#      真实签名与随机签名走完全相同的处理流程（含方向定向），保证可比
#
# 输入：
#   cohorts/*.rds            由 06_prepare_cohorts.R 生成
#   signature_library.csv    第 ② 步人工定稿的签名库
#                            必需列：sig_id, genes（分号分隔）
#                            可选列：coef（分号分隔，与 genes 等长）、direction
# 输出：
#   audit_out/per_cohort.csv      签名 x 队列 明细
#   audit_out/summary.csv         每个签名跨队列汇总
#   audit_out/null_dist.rds       零分布缓存
#
# 用法：Rscript 07_audit_signatures.R [签名库路径] [队列目录] [输出目录]
# =====================================================================

suppressPackageStartupMessages({
  need <- c("survival", "data.table")
  for (p in need) {
    if (!requireNamespace(p, quietly = TRUE))
      install.packages(p, repos = "https://cloud.r-project.org")
  }
  library(survival); library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
LIB   <- ifelse(length(args) >= 1, args[1],
                "sigminer/output/signature_library_main.csv")
COH   <- ifelse(length(args) >= 2, args[2], "cohorts")
OUT   <- ifelse(length(args) >= 3, args[3], "audit_out")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

B_PERM      <- 1000   # 随机基因集重复次数
MIN_SAMPLES <- 40     # 队列最少可分析样本数
HORIZON     <- 5      # 时间依赖 AUC 的评估时点（年）

`%||%` <- function(a, b) if (is.null(a)) b else a

# ============================================================ 基础函数
zscore <- function(m) {
  mu <- rowMeans(m, na.rm = TRUE)
  sd <- apply(m, 1, sd, na.rm = TRUE)
  sd[is.na(sd) | sd < 1e-8] <- 1
  (m - mu) / sd
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

cindex <- function(time, status, pred) {
  ok <- !is.na(pred) & !is.na(time) & !is.na(status)
  if (sum(ok) < 20 || sum(status[ok]) < 5) return(NA_real_)
  if (length(unique(pred[ok])) < 2) return(NA_real_)
  # 关键：survival::concordance 默认约定是【预测值越大 = 生存越长】，
  # 直接用在风险分上会得到 1-C。必须 reverse=TRUE 才是标准 C-index
  # （实测：x 与生存时间负相关时，默认返回 0.26，reverse=TRUE 返回 0.74）
  cf <- try(survival::concordance(Surv(time[ok], status[ok]) ~ pred[ok],
                                  reverse = TRUE), silent = TRUE)
  if (inherits(cf, "try-error")) return(NA_real_)
  unname(cf$concordance)
}

# 定向：真实签名与随机签名都做同样处理，保证零分布可比
oriented_cindex <- function(time, status, pred) {
  c0 <- cindex(time, status, pred)
  if (is.na(c0)) return(c(C = NA_real_, flipped = NA))
  c1 <- cindex(time, status, -pred)
  if (is.na(c1)) return(c(C = c0, flipped = 0))
  if (c1 > c0) c(C = c1, flipped = 1) else c(C = c0, flipped = 0)
}

auc_t <- function(time, status, pred, t = HORIZON) {
  # 时间依赖 AUC（Kaplan-Meier 逆概率加权，不依赖 timeROC 版本）
  ok <- !is.na(pred) & !is.na(time) & !is.na(status)
  tt <- time[ok]; ss <- status[ok]; pp <- pred[ok]
  if (sum(ok) < 30 || sum(ss) < 5) return(NA_real_)
  tt <- tt / 12  # 月 -> 年
  if (length(unique(pp)) < 2) return(NA_real_)
  # 只保留 t 时刻之前有信息的样本
  mk <- function(x) as.integer(tt <= t & ss == 1)  # case
  n_case <- sum(tt <= t & ss == 1); n_ctrl <- sum(tt > t)
  if (n_case < 5 || n_ctrl < 5) return(NA_real_)
  # KM 权重
  fit <- survfit(Surv(tt, ss) ~ 1)
  G <- summary(fit, times = t)$surv
  if (is.na(G) || G <= 0) return(NA_real_)
  w <- 1 / G
  # cases：tt <= t 且事件；controls：tt > t
  ci <- which(tt <= t & ss == 1)
  cj <- which(tt > t)
  if (length(ci) == 0 || length(cj) == 0) return(NA_real_)
  # 抽样加速
  if (length(ci) * length(cj) > 4e6) {
    ci <- sample(ci, min(length(ci), 2000))
    cj <- sample(cj, min(length(cj), 2000))
  }
  num <- 0
  for (k in ci) num <- num + sum(pp[cj] < pp[k]) + 0.5 * sum(pp[cj] == pp[k])
  num / (length(ci) * length(cj))
}

calib_slope <- function(time, status, pred) {
  ok <- !is.na(pred) & !is.na(time) & !is.na(status)
  if (sum(ok) < 30 || sum(status[ok]) < 5) return(NA_real_)
  cf <- try(coxph(Surv(time[ok], status[ok]) ~ pred[ok]), silent = TRUE)
  if (inherits(cf, "try-error")) return(NA_real_)
  unname(coef(cf)[1])
}

delta_cindex <- function(time, status, pred, clin_df) {
  # ΔC-index：在临床变量基础上加入风险分带来的判别增益
  d <- clin_df
  d$.t <- time; d$.s <- status; d$.p <- pred
  ok <- complete.cases(d[, c(".t", ".s", ".p")])
  d <- d[ok, , drop = FALSE]
  if (nrow(d) < 40 || sum(d$.s) < 5) return(NA_real_)
  vars <- intersect(c("stage", "age", "sex", "site"), colnames(d))
  vars <- vars[sapply(d[, vars, drop = FALSE], function(x) length(unique(x)) > 1)]
  if (length(vars) == 0) return(NA_real_)
  if (length(vars) == 0) return(NA_real_)
  # 【已修复 bug】原先把含 NA 协变量的 d 整表交给 coxph，coxph 会 silently 丢行，
  # predict 因此返回比 nrow(d) 短的向量；cindex 里 pred[ok] 与 time[ok] 长度不等时
  # R 会【回收短向量】造成逐行错位，算出的 ΔC 是垃圾值（GSE39582 缺 1 例、
  # TCGA 缺 11 例分期，都命中）。改在这儿先子集到协变量完整的样本，
  # 让三个向量行数天然一致，比事后按名字回溯更不容易出错。
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

net_benefit <- function(time, status, pred, t = HORIZON, pt = 0.15) {
  # 简化 DCA：以 t 年结局为二分类，阈值取 1/3 风险
  ok <- !is.na(pred) & !is.na(time) & !is.na(status)
  tt <- time[ok]/12; ss <- status[ok]; pp <- pred[ok]
  if (sum(ok) < 40) return(NA_real_)
  y <- as.integer(tt <= t & ss == 1)
  if (sum(y) < 5) return(NA_real_)
  thr <- quantile(pp, 1 - pt, na.rm = TRUE)
  flag <- as.integer(pp >= thr)
  n <- length(y)
  tp <- sum(flag == 1 & y == 1); fp <- sum(flag == 1 & y == 0)
  (tp / n) - (fp / n) * (pt / (1 - pt))
}

# ============================================================ 零分布
null_cache_file <- file.path(OUT, "null_dist.rds")
# 用 environment 存缓存：在函数内对 list 做 [[<- 只会改局部副本，缓存不会生效
nc_env <- new.env(parent = emptyenv())
if (file.exists(null_cache_file)) {
  tmp <- readRDS(null_cache_file)
  for (k in names(tmp)) assign(k, tmp[[k]], envir = nc_env)
}

null_for <- function(coh, n_gene, z, time, status, B = B_PERM) {
  key <- paste0(coh, "|", n_gene)
  if (!is.null(nc_env[[key]])) return(nc_env[[key]])
  pool <- rownames(z)
  if (length(pool) < n_gene * 2) return(NULL)
  vals <- numeric(B)
  for (b in seq_len(B)) {
    g <- sample(pool, n_gene)
    p <- as.numeric(colMeans(z[g, , drop = FALSE]))
    vals[b] <- oriented_cindex(time, status, p)[["C"]]
  }
  vals <- vals[!is.na(vals)]
  assign(key, vals, envir = nc_env)
  vals
}

emp_p <- function(real, null_vals) {
  if (is.na(real) || is.null(null_vals) || length(null_vals) == 0) return(NA_real_)
  (sum(null_vals >= real) + 1) / (length(null_vals) + 1)
}

# ============================================================ 载入
message("=== 载入签名库 ===")
if (!file.exists(LIB)) stop("找不到签名库：", LIB)
sig <- fread(LIB, sep = ",", header = TRUE, fill = TRUE,
             colClasses = "character", quote = "\"")
if (!all(c("sig_id", "genes") %in% colnames(sig)))
  stop("签名库必须包含 sig_id 与 genes 两列")
message("  签名数：", nrow(sig))

cohort_files <- list.files(COH, pattern = "\\.rds$", full.names = TRUE)
if (length(cohort_files) == 0) stop("队列目录为空：", COH, "（请先跑 06）")
message("  队列数：", length(cohort_files))

res_list <- list()

for (cf in cohort_files) {
  coh <- sub("\\.rds$", "", basename(cf))
  d <- readRDS(cf)
  z  <- zscore(d$expr)
  cl <- d$clin
  cl <- cl[match(colnames(z), cl$sample), , drop = FALSE]
  time <- suppressWarnings(as.numeric(cl$time))
  status <- suppressWarnings(as.integer(cl$status))
  keep <- !is.na(time) & !is.na(status)
  if (sum(keep) < MIN_SAMPLES) {
    message("  [跳过] ", coh, " 可分析样本不足（", sum(keep), "）"); next
  }
  z <- z[, keep, drop = FALSE]
  time <- time[keep]; status <- status[keep]
  cl2 <- cl[keep, , drop = FALSE]
  message("  ", coh, "：", ncol(z), " 样本，事件 ", sum(status))

  for (i in seq_len(nrow(sig))) {
    genes <- trimws(unlist(strsplit(sig$genes[i], "[;,]")))
    genes <- genes[genes != ""]
    coef  <- if (!is.null(sig$coef)) {
      suppressWarnings(as.numeric(trimws(unlist(strsplit(sig$coef[i] %||% "", "[;,]")))))
    } else NULL
    if (length(coef) != length(genes)) coef <- NULL
    dirn  <- suppressWarnings(as.numeric(sig$direction[i] %||% "1"))
    if (is.na(dirn)) dirn <- 1

    miss <- sum(!genes %in% rownames(z))
    p <- risk_score(z, genes, coef, dirn)
    if (is.null(p)) {
      res_list[[length(res_list) + 1]] <- data.table(
        sig_id = sig$sig_id[i], cohort = coh, n_gene = length(genes),
        n_measured = 0, C = NA_real_, flipped = NA, auc5 = NA_real_,
        calib_slope = NA_real_, delta_C = NA_real_, nb = NA_real_,
        null_median = NA_real_, p_emp = NA_real_, status_r = "genes_not_found")
      next
    }
    oc   <- oriented_cindex(time, status, p)
    # 定向后的预测值：其余所有指标都必须用它，否则方向翻转时
    # 校准斜率、ΔC-index、DCA 会整体反号（模拟测试已复现此问题）
    p_o  <- if (identical(as.numeric(oc[["flipped"]]), 1)) -p else p
    # 关键：零分布必须用【实际参与打分的基因数】，不能是请求的基因数。
    # 部分匹配时（如 9 基因签名只测到 8 个），风险分只用 8 个基因算，
    # 若零分布仍随机抽 9 个，随机集会占便宜 -> 系统性压低真实签名的相对表现。
    # 实测本次 1920 条记录中有 18.8% 存在部分匹配，最多缺 7 个基因。
    n_used <- length(genes) - miss
    nv   <- null_for(coh, n_used, z, time, status)
    res_list[[length(res_list) + 1]] <- data.table(
      sig_id     = sig$sig_id[i],
      cohort     = coh,
      n_gene     = length(genes),
      n_measured = length(genes) - miss,
      C          = oc[["C"]],
      flipped    = oc[["flipped"]],
      auc5       = auc_t(time, status, p_o),
      calib_slope = calib_slope(time, status, p_o),
      delta_C    = delta_cindex(time, status, p_o, cl2),
      nb         = net_benefit(time, status, p_o),
      null_median = if (is.null(nv)) NA_real_ else median(nv),
      p_emp      = emp_p(oc[["C"]], nv),
      status_r   = "ok")
  }
  saveRDS(as.list(nc_env), null_cache_file)
}

res <- if (length(res_list)) rbindlist(res_list, fill = TRUE) else
  data.table(sig_id = character(), cohort = character(), n_gene = integer(),
             n_measured = integer(), C = numeric(), flipped = numeric(),
             auc5 = numeric(), calib_slope = numeric(), delta_C = numeric(),
             nb = numeric(), null_median = numeric(), p_emp = numeric(),
             status_r = character(), p_adj = numeric())
if (nrow(res) == 0) {
  message("没有任何可评估的签名 x 队列组合。检查：1) cohorts/ 是否有 rds；2) 签名基因名是否为基因符号（探针号匹配不上）。")
} else {
  # 多重检验校正：几百个签名 x 多个队列，不做 FDR 会有大量假阳性
  res[, p_adj := p.adjust(p_emp, method = "BH")]
}

# ---------------------------------------------------------- 汇总
summ <- res[status_r == "ok", .(
  n_cohort   = .N,
  C_mean     = round(mean(C, na.rm = TRUE), 4),
  C_min      = round(min(C, na.rm = TRUE), 4),
  C_max      = round(max(C, na.rm = TRUE), 4),
  C_sd       = round(sd(C, na.rm = TRUE), 4),
  auc5_mean  = round(mean(auc5, na.rm = TRUE), 4),
  deltaC_mean = round(mean(delta_C, na.rm = TRUE), 4),
  p_emp_med  = round(median(p_emp, na.rm = TRUE), 4),
  n_sig_p05  = sum(p_emp < 0.05, na.rm = TRUE),
  n_sig_fdr  = sum(p_adj < 0.05, na.rm = TRUE)
), by = sig_id]
summ <- merge(summ, unique(res[, .(sig_id, n_gene)]), by = "sig_id", all.x = TRUE)
summ[, beats_random := n_sig_fdr >= ceiling(n_cohort / 2)]
setorder(summ, -C_mean)

fwrite(res,  file.path(OUT, "per_cohort.csv"))
fwrite(summ, file.path(OUT, "summary.csv"))

cat("\n=== 完成 ===\n")
cat("签名 x 队列 明细：", nrow(res), " 行 -> ", file.path(OUT, "per_cohort.csv"), "\n", sep = "")
cat("签名级汇总：", nrow(summ), " 个 -> ", file.path(OUT, "summary.csv"), "\n", sep = "")
cat("\n在半数以上队列中优于随机基因集的签名：",
    sum(summ$beats_random, na.rm = TRUE), " / ", nrow(summ), "（FDR<0.05 且过半数队列）\n", sep = "")
cat("C-index 中位数：", round(median(summ$C_mean, na.rm = TRUE), 3), "\n", sep = "")
