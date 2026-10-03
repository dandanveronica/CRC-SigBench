#!/usr/bin/env Rscript
# =====================================================================
# 22_cutpoint_manipulation.R
# 量化"截点操纵"能凭空制造多少显著性
#
# 背景（Part 3 的发现）
#   同一个签名、同一份数据：
#     用训练集锁死的截点 -> GSE17536  HR=1.017, p=0.943（完全无效）
#     就地搜索最优截点   -> GSE17536  HR=1.637, p=0.049（"显著"）
#   如果这在 320 个签名里普遍成立，那么"显著"的门槛就是可以调出来的。
#
# 设计
#   对【每个签名 × 每个队列】：
#     (a) 预指定截点：风险分中位数（不利用结局信息）—— 合规
#     (b) 就地最优截点：在 20%~80% 分位间搜索最小 p 值 —— 违规
#   统计两种做法下的"显著率"（p<0.05），差值即截点操纵的收益。
#
#   关键校正：搜索 k 个截点会把假阳性率从 5% 抬到约 1-(1-0.05)^k_eff，
#   因此要与【随机基因集】在完全相同流程下的显著率对照，
#   这样才能区分"真的有信号"和"纯靠搜索制造"。
# =====================================================================
suppressPackageStartupMessages({
  library(survival); library(data.table)
})

BASE <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)))
if (!length(BASE)) BASE <- getwd(); setwd(BASE)
OUT <- "audit_out_clean"

zscore <- function(m) {
  mu <- rowMeans(m, na.rm = TRUE); sd <- apply(m, 1, sd, na.rm = TRUE)
  sd[is.na(sd) | sd < 1e-8] <- 1
  (m - mu) / sd
}

# 对一条风险分，比较预指定截点 vs 最优截点
cmp_cut <- function(time, status, r, probs = seq(.2, .8, by = .05)) {
  ok <- is.finite(r) & !is.na(time) & !is.na(status) & time > 0
  t <- time[ok]; s <- status[ok]; r <- r[ok]
  if (length(t) < 40 || sum(s) < 5) return(NULL)
  sv <- Surv(t, s)
  # (a) 预指定：中位数
  p_med <- NA_real_; hr_med <- NA_real_
  h <- r > median(r)
  if (length(unique(h)) == 2 && min(table(h)) >= 5) {
    f <- try(coxph(sv ~ h), silent = TRUE)
    if (!inherits(f, "try-error")) {
      p_med <- summary(f)$coefficients[1, 5]; hr_med <- unname(exp(coef(f)))
    }
  }
  # (b) 就地最优：搜索最小 p
  p_best <- NA_real_; hr_best <- NA_real_
  for (q in quantile(r, probs)) {
    h <- r > q
    if (length(unique(h)) != 2 || min(table(h)) < 5) next
    f <- try(coxph(sv ~ h), silent = TRUE)
    if (inherits(f, "try-error")) next
    pp <- summary(f)$coefficients[1, 5]
    if (!is.finite(pp)) next
    if (is.na(p_best) || pp < p_best) { p_best <- pp; hr_best <- unname(exp(coef(f))) }
  }
  data.table(p_median = p_med, HR_median = hr_med,
             p_best = p_best, HR_best = hr_best)
}

# ---------------- 载入 ----------------
lib <- fread("sigminer/output/signature_library_clean.csv", encoding = "UTF-8")
files <- list.files("cohorts", pattern = "\\.rds$", full.names = TRUE)
cohorts <- sub("\\.rds$", "", basename(files))
cat("签名：", nrow(lib), " 队列：", paste(cohorts, collapse = " "), "\n\n")

CL <- list()
for (f in files) {
  d <- readRDS(f); z <- zscore(d$expr)
  cl <- d$clin; cl <- cl[match(colnames(z), cl$sample), , drop = FALSE]
  tm <- suppressWarnings(as.numeric(cl$time)); st <- suppressWarnings(as.integer(cl$status))
  kp <- !is.na(tm) & !is.na(st) & tm > 0
  CL[[sub("\\.rds$", "", basename(f))]] <-
    list(z = z[, kp, drop = FALSE], time = tm[kp], status = st[kp])
}

# ---------------- 已发表签名 ----------------
cat("① 已发表签名：预指定截点 vs 就地最优截点\n")
rows <- rbindlist(lapply(cohorts, function(ch) {
  d <- CL[[ch]]; nk <- seq_len(nrow(d$z))
  rbindlist(lapply(seq_len(nrow(lib)), function(i) {
    g <- intersect(strsplit(lib$genes[i], ";")[[1]], rownames(d$z))
    if (length(g) < 2) return(NULL)
    r <- as.numeric(colMeans(d$z[g, , drop = FALSE]))
    o <- cmp_cut(d$time, d$status, r)
    if (is.null(o)) return(NULL)
    data.table(sig_id = lib$sig_id[i], cohort = ch, n_gene = length(g), o)
  }))
}))
cat("  记录数：", nrow(rows), "\n")

summ <- function(x, lab) {
  v1 <- x$p_median; v2 <- x$p_best
  v1 <- v1[!is.na(v1)]; v2 <- v2[!is.na(v2)]
  cat(sprintf("  %-14s n=%5d  预指定显著率 %.1f%%   最优截点显著率 %.1f%%   ×%.2f\n",
              lab, length(v2), 100 * mean(v1 < .05, na.rm = TRUE),
              100 * mean(v2 < .05, na.rm = TRUE),
              mean(v2 < .05, na.rm = TRUE) / max(mean(v1 < .05, na.rm = TRUE), 1e-9)))
  c(median(v1, na.rm = TRUE), median(v2, na.rm = TRUE))
}
cat("\n  —— 已发表签名 ——\n"); summ(rows, "全部队列")

# ---------------- 随机基因集对照（同样流程） ----------------
cat("\n② 随机基因集对照（完全相同的搜索流程）\n")
NREP <- 20     # 每个（队列x基因数）重复次数（内存吃紧时调小）
need <- unique(rows[, .(cohort, n_gene)])
set.seed(20261001)
rows_n <- rbindlist(lapply(seq_len(nrow(need)), function(j) {
  ch <- need$cohort[j]; ng <- need$n_gene[j]; d <- CL[[ch]]
  nk <- seq_len(nrow(d$z))
  rbindlist(lapply(1:NREP, function(k) {
    r <- as.numeric(colMeans(d$z[sample(nk, ng), , drop = FALSE]))
    o <- cmp_cut(d$time, d$status, r)
    if (is.null(o)) return(NULL)
    data.table(sig_id = paste0("NULL_", ch, "_", ng, "_", k), cohort = ch,
               n_gene = ng, o)
  }))
}))
cat("  记录数：", nrow(rows_n), "\n")
cat("  —— 随机基因集 ——\n"); summ(rows_n, "全部队列")

# ---------------- 汇总 ----------------
tab <- rbind(
  data.table(对象 = "已发表签名(320)", 预指定_显著率 = mean(rows$p_median < .05, na.rm = TRUE),
             最优截点_显著率 = mean(rows$p_best < .05, na.rm = TRUE),
             中位p_预指定 = median(rows$p_median, na.rm = TRUE),
             中位p_最优 = median(rows$p_best, na.rm = TRUE), n = nrow(rows)),
  data.table(对象 = "随机基因集", 预指定_显著率 = mean(rows_n$p_median < .05, na.rm = TRUE),
             最优截点_显著率 = mean(rows_n$p_best < .05, na.rm = TRUE),
             中位p_预指定 = median(rows_n$p_median, na.rm = TRUE),
             中位p_最优 = median(rows_n$p_best, na.rm = TRUE), n = nrow(rows_n))
)
cat("\n=== 汇总 ===\n")
print(tab[, .(对象, n, 预指定_显著率 = round(100 * 预指定_显著率, 1),
             最优截点_显著率 = round(100 * 最优截点_显著率, 1),
             中位p_预指定 = round(中位p_预指定, 3), 中位p_最优 = round(中位p_最优, 3))])

gain_pub <- mean(rows$p_best < .05, na.rm = TRUE) - mean(rows$p_median < .05, na.rm = TRUE)
gain_nul <- mean(rows_n$p_best < .05, na.rm = TRUE) - mean(rows_n$p_median < .05, na.rm = TRUE)
cat(sprintf("\n★ 截点操纵的收益：已发表 +%.1f 个百分点，随机基因集 +%.1f 个百分点\n",
            100 * gain_pub, 100 * gain_nul))
cat(sprintf("  即：最优截点带来的额外'显著'中，约 %.0f%% 在随机基因集上同样能得到（纯搜索产物）\n",
            100 * gain_nul / max(gain_pub, 1e-9)))

fwrite(rows, file.path(OUT, "截点操纵_已发表签名.csv"), encoding = "UTF-8")
fwrite(rows_n, file.path(OUT, "截点操纵_随机基因集.csv"), encoding = "UTF-8")
fwrite(tab, file.path(OUT, "截点操纵_汇总.csv"), encoding = "UTF-8")
cat("\n已写出：截点操纵_已发表签名.csv / 截点操纵_随机基因集.csv / 截点操纵_汇总.csv\n")
