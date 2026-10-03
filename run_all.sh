#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# CRC-SigBench — top-level driver
#
#   bash run_all.sh                  # stage B + C (fully reproducible from public inputs)
#   bash run_all.sh --with-extraction # also stage A (needs publisher subscriptions; see README)
#
# Stage A extracts claimed performance values from publisher full texts. Those
# full texts are not redistributable, so stage A cannot be re-run by a third
# party; its OUTPUT is deposited in data/signatures/ instead.
# ---------------------------------------------------------------------------
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$ROOT"

WITH_EXTRACTION=0
for a in "$@"; do [ "$a" = "--with-extraction" ] && WITH_EXTRACTION=1; done

# --- path compatibility -----------------------------------------------------
# The R scripts resolve their inputs through the relative paths  cohorts/  and
# audit_out_clean/  (some from getwd(), some from dirname(script)). Running them
# from sub-directories would therefore point at the wrong place, so we expose
# both directories from inside every code/*/ folder as symlinks.
mkdir -p cohorts audit_out_clean
for d in code/*/; do
  ln -sfn "$ROOT/cohorts"         "$d/cohorts"
  ln -sfn "$ROOT/audit_out_clean" "$d/audit_out_clean"
done

run_r()  { echo; echo "==> R      $1"; Rscript "$1"; }
run_py() { echo; echo "==> python $1"; python3 "$1"; }

echo "CRC-SigBench driver | root=$ROOT"

# --- stage A: literature extraction (opt-in) ---------------------------------
if [ "$WITH_EXTRACTION" = "1" ]; then
  echo; echo "### STAGE A — literature extraction (requires publisher full-text access)"
  run_py code/02_signatures/09_metadata.py
  run_py code/02_signatures/10_fetch_fulltext.py
  run_py code/02_signatures/18_claim_extract.py
  run_py code/02_signatures/24_export_claim_sample.py
else
  echo; echo "### STAGE A skipped — using deposited data/signatures/claimed_values_*.csv"
  echo "###   (re-run with --with-extraction if you hold the publisher licences)"
fi

# --- stage B: cohorts, scoring, evaluation ----------------------------------
echo; echo "### STAGE B — cohort construction and evaluation (reproducible)"
run_r  code/01_build_cohorts/06b_build_cohorts.R
run_py code/01_build_cohorts/06c_build_gse87211.py
run_py code/01_build_cohorts/30_tcga_xena.py
run_r  code/01_build_cohorts/31_tcga_to_rds.R
run_r  code/01_build_cohorts/20_assemble_cohorts.R

run_r  code/03_score/07_audit_signatures.R

run_r  code/05_null/07b_highres_null.R
run_r  code/05_null/17_null_matched.R

run_r  code/04_evaluate/07c_clinical_benchmark.R
run_r  code/04_evaluate/34_clinical_baseline_final.R
run_py code/04_evaluate/32_tcga_stratified.py

run_r  code/06_reclassification/27_nri_idi.R
run_r  code/06_reclassification/28_nri_cv.R
run_py code/06_reclassification/29_nri_summary.py

run_py code/07_attribution/19_part2_attribution.py
run_py code/07_attribution/26_usability_and_trend.py
run_py code/07_attribution/25d_home_final.py
run_py code/07_attribution/15_compare_tiers.py

run_r  code/08_cutpoint/22_cutpoint_manipulation.R

run_r  code/09_rederivation/14_de_novo_benchmark.R
run_r  code/09_rederivation/14b_de_novo_full.R
run_r  code/09_rederivation/21_part3_compliant.R

# --- stage C: figures and the numbers quoted in the paper -------------------
echo; echo "### STAGE C — figures and manuscript numbers"
run_py code/10_figures/37_make_figures.py
run_py code/10_figures/38_make_prisma.py
run_py code/99_auxiliary/35_numbers_for_paper.py

echo; echo "Done. Fresh outputs are under audit_out_clean/;"
echo "the deposited snapshot of those outputs is results/."
