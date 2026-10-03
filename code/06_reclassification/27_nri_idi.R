#!/usr/bin/env Rscript
# =====================================================================
# 27_nri_idi.R
# 补齐 Part 1 承诺但未实现的两个指标：NRI 与 IDI
#
# 设计要点
# --------
# 1) 用【连续型 NRI（category-free NRI, NRI(>0)）】，不设风险切点。
#    理由：本研究的核心论点之一就是"就地搜最优切点会凭空造出显著性"
#    （见 22_cutpoint_manipulation.R）。若这里再用切点型 NRI，等于自相矛盾。
#    连续型 NRI 只看预测概率的【移动方向】，不需要任何切点。
#
# 2) 时间窗 = 5 年 = 60 个月（cohorts 的 time 单位是月，与 07 的 HORIZON=5 一致，
#    07 里 tt <- tt/12 转年）。
#
# 3) case / control 定义与 07 的 auc_t 完全一致：
#      case    = 5 年内发生事件（tt <= 60 & status == 1）
#      control = 随访超过 5 年仍存活（tt > 60）
#      5 年内删失者（tt <= 60 & status == 0）剔除——结局不确定，不能判类
#
# 4) 基线模型（old）每个队列只拟合一次：Surv ~ stage + age + sex + site
#    新模型（new）= 基线 + 签名风险分。风险分方向 oriented（与 07 一致）。
#
# 5) 每个队列另做 100 个随机等大基因集，给出 NRI/IDI 的随机基线，
#    与 C 的 null_median 平行，便于"是否优于瞎猜"的判断。
#
# 6) 内存：一次只加载一个队列（教训见 2026-10-01 日志——7 队列同时加载必爆 swap）
# =====================================================================

suppressPackageStartupMessages({ library(survival); library(data.table) })

COH  <- "cohorts"
LIB  <- "sigminer/output/signature_library_clean.csv"
OUT  <- "audit_out_clean"
T0   <- 60          # 5 年（月）
N_RAND <- 100       # 每队列随机基因集个数

zscore <- function(m) {
  mu <- rowMeans(m, na.rm = TRUE); sd <- apply(m, 1, sd, na.rm = TRUE)
  sd[is.na(sd) | sd < 1e-8] <- 1; (m - mu) / sd
}
cindex <- function(time, status, pred) {
  ok <- !is.na(pred) & !is.na(time) & !is.na(status)
  if (sum(ok) < 20 || sum(status[ok]) < 5) return(NA_real_)
  if (length(unique(pred[ok])) < 2) return(NA_real_)
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

# 5 年风险概率 P(事件 ≤ 5年) = 1 - S(60)
risk5 <- function(fit, newdata) {
  sf <- try(survfit(fit, newdata = newdata), silent = TRUE)
  if (inherits(sf, "try-error")) return(rep(NA_real_, nrow(newdata)))
  s <- summary(sf, times = T0, extend = TRUE)$surv
  if (is.matrix(s)) s <- s[nrow(s), ]          # 取最后一个时间点（= T0）
  as.numeric(1 - s)
}

# 连续型 NRI + IDI（Pencina et al. 2008 / 2011）
nri_idi <- function(p_old, p_new, case) {
  # case: 逻辑向量，TRUE = 5 年内发生事件
  ctrl <- !case
  n1 <- sum(case); n0 <- sum(ctrl)
  if (n1 < 5 || n0 < 5) return(NULL)
  d  <- p_new - p_old
  up <- d >  1e-12
  dn <- d < -1e-12
  pu1 <- sum(up[case]) / n1; pd1 <- sum(dn[case]) / n1
  pu0 <- sum(up[ctrl]) / n0; pd0 <- sum(dn[ctrl]) / n0
  nri <- (pu1 - pd1) - (pu0 - pd0)
  # Pencina 方差
  v1 <- (pu1 + pd1 - (pu1 - pd1)^2) / n1
  v0 <- (pu0 + pd0 - (pu0 - pd0)^2) / n0
  vn <- v1 + v0
  z_nri <- if (vn > 0) nri / sqrt(vn) else NA_real_
  p_nri <- if (!is.na(z_nri)) 2 * (1 - pnorm(abs(z_nri))) else NA_real_
  # IDI
  is_case <- mean(p_new[case]) - mean(p_old[case])
  is_ctrl <- mean(p_new[ctrl]) - mean(p_old[ctrl])
  idi <- is_case - is_ctrl
  vi <- var(d[case]) / n1 + var(d[ctrl]) / n0
  z_idi <- if (vi > 0) idi / sqrt(vi) else NA_real_
  p_idi <- if (!is.na(z_idi)) 2 * (1 - pnorm(abs(z_idi))) else NA_real_
  list(n_case = n1, n_ctrl = n0, NRI = nri, NRI_z = z_nri, NRI_p = p_nri,
       IDI = idi, IDI_z = z_idi, IDI_p = p_idi,
       p_up_case = pu1, p_down_case = pd1, p_up_ctrl = pu0, p_down_ctrl = pd0)
}

lib <- fread(LIB, sep = ",", header = TRUE, fill = TRUE,
             colClasses = "character", quote = "\"")

ROWS <- list(); RND <- list()
for (cf in list.files(COH, pattern = "\\.rds$", full.names = TRUE)) {
  coh <- sub("\\.rds$", "", basename(cf))
  cat(sprintf("\n===== %s =====\n", coh))
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
  rm(z); invisible(gc(FALSE))

  # case / control
  case  <- (tt <= T0 & ss == 1)
  ctrl  <- (tt >  T0)
  usable <- case | ctrl
  cat(sprintf("  n=%d  事件=%d  5年内事件(case)=%d  超5年存活(control)=%d  可用=%d  临床变量=%s\n",
              length(tt), sum(ss), sum(case), sum(ctrl), sum(usable),
              paste(vars_all, collapse = "+")))
  if (sum(case) < 5 || sum(ctrl) < 5) { rm(d); next }

  base_f <- as.formula(paste("Surv(.t, .s) ~", paste(vars_all, collapse = " + ")))
  fit0 <- try(coxph(base_f, data = dd), silent = TRUE)
  if (inherits(fit0, "try-error")) { rm(d); next }
  p_old_all <- risk5(fit0, dd)
  if (all(is.na(p_old_all))) { rm(d); next }
  C_base <- cindex(tt, ss, predict(fit0, type = "lp"))

  # ---- 逐签名
  ng_ok <- 0
  for (i in seq_len(nrow(lib))) {
    genes <- trimws(unlist(strsplit(lib$genes[i], "[;,]")))
    genes <- genes[genes != ""]
    g <- intersect(genes, rownames(zz))
    if (length(g) < 2) next
    rr <- oriented(tt, ss, as.numeric(colMeans(zz[g, , drop = FALSE])))
    if (is.na(rr$C)) next
    dd2 <- dd; dd2$.p <- rr$p
    f1 <- as.formula(paste("Surv(.t, .s) ~ .p +", paste(vars_all, collapse = " + ")))
    fit1 <- try(coxph(f1, data = dd2), silent = TRUE)
    if (inherits(fit1, "try-error")) next
    p_new <- risk5(fit1, dd2)
    if (all(is.na(p_new))) next
    okp <- usable & !is.na(p_old_all) & !is.na(p_new)
    r <- nri_idi(p_old_all[okp], p_new[okp], case[okp])
    if (is.null(r)) next
    ng_ok <- ng_ok + 1
    ROWS[[length(ROWS) + 1]] <- data.table(
      cohort = coh, sig_id = lib$sig_id[i], n_gene = length(genes),
      n_measured = length(g), n_case = r$n_case, n_ctrl = r$n_ctrl,
      C_clinical = round(C_base, 4),
      C_plus_sig = round(cindex(tt, ss, predict(fit1, type = "lp")), 4),
      NRI = round(r$NRI, 4), NRI_z = round(r$NRI_z, 4), NRI_p = signif(r$NRI_p, 4),
      IDI = round(r$IDI, 5), IDI_z = round(r$IDI_z, 4), IDI_p = signif(r$IDI_p, 4),
      p_up_case = round(r$p_up_case, 4), p_down_case = round(r$p_down_case, 4),
      p_up_ctrl = round(r$p_up_ctrl, 4), p_down_ctrl = round(r$p_down_ctrl, 4))
  }
  cat(sprintf("  %d 个签名完成 NRI/IDI\n", ng_ok))

  # ---- 随机基因集基线（同队列、随机等大）
  set.seed(20261001L)
  nr <- 0
  for (b in seq_len(N_RAND)) {
    ng <- sample(c(3L, 5L, 8L, 12L, 20L), 1L)
    gs <- sample(rownames(zz), ng)
    rb <- oriented(tt, ss, as.numeric(colMeans(zz[gs, , drop = FALSE])))
    if (is.na(rb$C)) next
    dd3 <- dd; dd3$.p <- rb$p
    fb <- try(coxph(as.formula(paste("Surv(.t, .s) ~ .p +",
                                     paste(vars_all, collapse = " + "))),
                    data = dd3), silent = TRUE)
    if (inherits(fb, "try-error")) next
    pn <- risk5(fb, dd3)
    okp <- usable & !is.na(p_old_all) & !is.na(pn)
    r <- nri_idi(p_old_all[okp], pn[okp], case[okp])
    if (is.null(r)) next
    nr <- nr + 1
    RND[[length(RND) + 1]] <- data.table(cohort = coh, rep = b, n_gene = ng,
                                         NRI = r$NRI, IDI = r$IDI,
                                         NRI_p = r$NRI_p, IDI_p = r$IDI_p)
  }
  cat(sprintf("  随机基因集基线：%d 个\n", nr))
  rm(d, zz, dd); invisible(gc(FALSE))
}

out <- rbindlist(ROWS)
fwrite(out, file.path(OUT, "per_cohort_NRI_IDI.csv"))
rnd <- if (length(RND)) rbindlist(RND) else NULL
if (!is.null(rnd)) fwrite(rnd, file.path(OUT, "NRI_IDI_随机基线.csv"))

cat("\n\n================ 汇总 ================\n")
agg <- out[, .(n_sig = .N,
               NRI中位 = round(median(NRI, na.rm = TRUE), 4),
               NRI四分位 = paste0(round(quantile(NRI, .25, na.rm = TRUE), 3), "~",
                                  round(quantile(NRI, .75, na.rm = TRUE), 3)),
               NRI为正比例 = paste0(round(100 * mean(NRI > 0, na.rm = TRUE), 1), "%"),
               NRI显著率 = paste0(round(100 * mean(NRI_p < 0.05, na.rm = TRUE), 1), "%"),
               IDI中位 = round(median(IDI, na.rm = TRUE), 5),
               IDI为正比例 = paste0(round(100 * mean(IDI > 0, na.rm = TRUE), 1), "%"),
               IDI显著率 = paste0(round(100 * mean(IDI_p < 0.05, na.rm = TRUE), 1), "%")), by = cohort]
print(agg)
cat("\n—— 全队列（记录级）——\n")
cat(sprintf("  NRI 中位 %.4f  为正 %.1f%%  p<0.05 占 %.1f%%\n",
            out[, median(NRI, na.rm = TRUE)], 100 * out[, mean(NRI > 0, na.rm = TRUE)],
            100 * out[, mean(NRI_p < 0.05, na.rm = TRUE)]))
cat(sprintf("  IDI 中位 %.5f  为正 %.1f%%  p<0.05 占 %.1f%%\n",
            out[, median(IDI, na.rm = TRUE)], 100 * out[, mean(IDI > 0, na.rm = TRUE)],
            100 * out[, mean(IDI_p < 0.05, na.rm = TRUE)]))
if (!is.null(rnd)) {
  cat("\n—— 随机基因集基线（应≈0）——\n")
  cat(sprintf("  NRI 中位 %.4f（%.1f%% 为正）  IDI 中位 %.5f（%.1f%% 为正）\n",
              rnd[, median(NRI, na.rm = TRUE)], 100 * rnd[, mean(NRI > 0, na.rm = TRUE)],
              rnd[, median(IDI, na.rm = TRUE)], 100 * rnd[, mean(IDI > 0, na.rm = TRUE)]))
  cat(sprintf("  随机集 NRI p<0.05 占 %.1f%%  ← 这是假阳性率参照\n",
              100 * rnd[, mean(NRI_p < 0.05, na.rm = TRUE)]))
}
fwrite(agg, file.path(OUT, "NRI_IDI_队列汇总.csv"))
cat("\n已存：per_cohort_NRI_IDI.csv / NRI_IDI_随机基线.csv / NRI_IDI_队列汇总.csv\n")
