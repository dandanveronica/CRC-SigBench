#!/usr/bin/env Rscript
# =====================================================================
# 34_clinical_baseline_final.R —— 验证：队列的临床基线越强，签名的 ΔC_临床 越低
#
# 待解释现象（8 队列，delta_C 已修复错位 bug 后）：
#   TCGA(RNA-seq)   ΔC_随机 +0.0201（8 队最高）   ΔC_临床 +0.0018（偏低）
#   7 个芯片队列    ΔC_随机 +0.0105              ΔC_临床 +0.0113
# 同一队列上，签名相对随机基因集更突出，相对临床变量却更没增量。
#
# 假说：这些差异主要由【临床基线模型的强度】决定，而不是签名本身。
# 判据：逐队列算纯临床模型 C-index，再与 ΔC_临床 做队列级相关（n=8）。
#
# 口径必须与 07 的 delta_cindex 完全一致：
#   vars = intersect(stage, age, sex, site)，再用 length(unique(x)) > 1 过滤，
#   然后 complete.cases 子集 —— 差一步就不可比。
# =====================================================================
suppressPackageStartupMessages({ library(survival) })

cindex <- function(time, status, pred) {
  ok <- is.finite(pred) & is.finite(time) & !is.na(status)
  time <- time[ok]; status <- status[ok]; pred <- pred[ok]
  if (length(time) < 20 || sum(status) < 5) return(NA_real_)
  if (length(unique(pred)) < 2) return(NA_real_)
  cf <- try(survival::concordance(Surv(time, status) ~ pred, reverse = TRUE),
            silent = TRUE)
  if (inherits(cf, "try-error")) return(NA_real_)
  unname(cf$concordance)
}

files <- sort(list.files("cohorts", pattern = "\\.rds$", full.names = TRUE))
res <- data.frame()
for (f in files) {
  d <- readRDS(f)
  cl <- d$clin
  cl <- cl[match(colnames(d$expr), cl$sample), , drop = FALSE]
  ok <- !is.na(cl$time) & !is.na(cl$status) & is.finite(cl$time) & cl$time > 0
  cl2 <- cl[ok, , drop = FALSE]
  cl2$.t <- cl2$time; cl2$.s <- cl2$status

  # ==== 与 07::delta_cindex 完全同口径 ====
  vars <- intersect(c("stage", "age", "sex", "site"), colnames(cl2))
  vars <- vars[sapply(cl2[, vars, drop = FALSE], function(x) length(unique(x)) > 1)]
  if (length(vars) == 0) next
  keep <- stats::complete.cases(cl2[, c(".t", ".s", vars), drop = FALSE])
  dd <- cl2[keep, , drop = FALSE]
  if (nrow(dd) < 40 || sum(dd$.s) < 5) next
  fm <- as.formula(paste("Surv(.t, .s) ~", paste(vars, collapse = " + ")))
  cb <- try(coxph(fm, data = dd), silent = TRUE)
  if (inherits(cb, "try-error")) next
  C_clin <- cindex(dd$.t, dd$.s, predict(cb, type = "lp"))

  res <- rbind(res, data.frame(
    cohort = sub("\\.rds$", "", basename(f)),
    n = nrow(dd), events = sum(dd$.s), C_clin = C_clin,
    vars = paste(vars, collapse = "+"),
    med_fu = median(dd$.t), stringsAsFactors = FALSE))
}

cat("=== 纯临床模型（口径与 07 的 delta_C 完全一致）===\n")
cat(sprintf("%-16s %6s %6s %9s %10s  %s\n", "队列", "n", "事件", "C_临床", "随访中位", "使用变量"))
cat(strrep("-", 92), "\n", sep = "")
for (i in seq_len(nrow(res)))
  cat(sprintf("%-16s %6d %6d %9.4f %10.1f  %s\n", res$cohort[i], res$n[i],
              res$events[i], res$C_clin[i], res$med_fu[i], res$vars[i]))

# 合并 per_cohort 的 ΔC_临床 / ΔC_随机 中位数
pc <- read.csv("audit_out_8/per_cohort.csv", fileEncoding = "UTF-8")
pc$d_rand <- pc$C - pc$null_median
agg <- stats::aggregate(cbind(d_clin = pc$delta_C, d_rand = pc$d_rand, C_sig = pc$C),
                        by = list(cohort = pc$cohort),
                        FUN = function(z) median(z, na.rm = TRUE))
M <- merge(res, agg, by = "cohort", all.x = TRUE)

cat("\n=== 队列级关联：临床基线越强 -> 签名还能榨出多少增量？（n=%d）===\n",
    nrow(M))
cat(sprintf("%-16s %9s %9s %9s\n", "队列", "C_临床", "ΔC_临床", "ΔC_随机"))
for (i in seq_len(nrow(M)))
  cat(sprintf("%-16s %9.4f %9.4f %9.4f\n", M$cohort[i], M$C_clin[i],
              M$d_clin[i], M$d_rand[i]))

cd <- cor.test(M$C_clin, M$d_clin, method = "spearman")
cr <- cor.test(M$C_clin, M$d_rand, method = "spearman")
cat(sprintf("\n  C_临床 vs ΔC_临床:  Spearman rho = %+.3f   p = %.4g\n",
            cd$estimate, cd$p.value))
cat(sprintf("  C_临床 vs ΔC_随机:  Spearman rho = %+.3f   p = %.4g\n",
            cr$estimate, cr$p.value))

rna <- M[M$cohort == "TCGA_COADREAD", ]
cat(sprintf("\n  TCGA 临床基线 %.4f，7 个 GEO 中位 %.4f，差 %+.4f\n",
            rna$C_clin, median(M$C_clin[M$cohort != "TCGA_COADREAD"]),
            rna$C_clin - median(M$C_clin[M$cohort != "TCGA_COADREAD"])))
cat(sprintf("  TCGA ΔC_临床 %.4f，7 个 GEO 中位 %.4f，差 %+.4f\n",
            rna$d_clin, median(M$d_clin[M$cohort != "TCGA_COADREAD"]),
            rna$d_clin - median(M$d_clin[M$cohort != "TCGA_COADREAD"])))

write.csv(M, "audit_out_clean/队列临床基线_合并.csv", row.names = FALSE,
          fileEncoding = "UTF-8")
cat("\n已存: audit_out_clean/队列临床基线_合并.csv\n")
