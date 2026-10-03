# 00_check_env.R — verify the environment before running the pipeline.
# Usage:  Rscript code/00_check_env.R
# Run this from the repository root.

cat("R version:", R.version.string, "\n")
cat("Platform :", R.version$platform, "\n\n")

r_pkgs <- c("survival", "glmnet", "data.table", "parallel", "timeROC")
for (p in r_pkgs) {
  ok <- requireNamespace(p, quietly = TRUE)
  v  <- if (ok) tryCatch(as.character(packageVersion(p)), error = function(e) "?") else "-"
  cat(sprintf("  [%-3s] %-12s %s\n", if (ok) "OK" else "MISS", p, v))
}

cat("\nPython side (used by the .py scripts):\n")
py <- Sys.which("python3")
if (nzchar(py)) {
  cat("  interpreter:", py, "\n")
  ver <- tryCatch(system2(py, "--version", stdout = TRUE, stderr = TRUE),
                  error = function(e) "?")
  cat("  version    :", paste(ver, collapse = " "), "\n")
  cat("\n  Check the Python packages with:\n")
  cat("    python3 -c \"import pandas, numpy, scipy, matplotlib, openpyxl, docx, requests; print('python deps OK')\"\n")
} else {
  cat("  python3 not found on PATH\n")
}

cat("\nWorking directory must be the repository root; run_all.sh handles this.\n")
cat("Public inputs (GEO, UCSC Xena, cBioPortal) are downloaded on first run.\n")
