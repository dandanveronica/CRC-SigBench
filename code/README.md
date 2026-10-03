# Code: execution order and conventions

## How to run

Always drive the pipeline from `run_all.sh` in the repository root:

```bash
bash run_all.sh                    # stage B + C
bash run_all.sh --with-extraction  # also stage A
```

Do not run individual scripts by hand unless you know what you are doing — several of them
assume that earlier stages have already written their outputs.

## Directory conventions

The scripts were written against a flat working directory that exposes two relative paths:

```
cohorts/           eight standardised per-cohort objects, built by code/01_build_cohorts/
audit_out_clean/   every result table
```

Some scripts resolve these from `getwd()`, others from `dirname(script)`. `run_all.sh`
therefore creates both directories at the repository root **and** symlinks them into every
`code/*/` folder, so either resolution style lands in the same place. Both directories are
git-ignored; the committed `results/` folder is the deposited snapshot of `audit_out_clean/`.

## Execution order

| Stage | Script | Role |
|---|---|---|
| A | `02_signatures/09_metadata.py` | signature metadata (needs publisher access) |
| A | `02_signatures/10_fetch_fulltext.py` | fetch full texts (needs publisher access) |
| A | `02_signatures/18_claim_extract.py` | extract claimed performance values |
| A | `02_signatures/24_export_claim_sample.py` | export the 30-item manual-verification sample |
| B | `01_build_cohorts/06b_build_cohorts.R` | GEO download, normalisation, clinical harmonisation |
| B | `01_build_cohorts/06c_build_gse87211.py` | GSE87211 construction |
| B | `01_build_cohorts/30_tcga_xena.py` | TCGA expression (Xena) + clinical (cBioPortal) |
| B | `01_build_cohorts/31_tcga_to_rds.R` | TCGA → standardised object |
| B | `01_build_cohorts/20_assemble_cohorts.R` | assemble and assert cross-cohort isomorphism |
| B | `03_score/07_audit_signatures.R` | signature reconstruction, scoring, C-index, AUC, calibration |
| B | `05_null/07b_highres_null.R` | random-gene-set null, B = 40,000 |
| B | `05_null/17_null_matched.R` | size-matched null |
| B | `04_evaluate/07c_clinical_benchmark.R` | clinical-variable Cox benchmark |
| B | `04_evaluate/34_clinical_baseline_final.R` | clinical baseline C-index per cohort |
| B | `04_evaluate/32_tcga_stratified.py` | TCGA-stratified analyses |
| B | `06_reclassification/27_nri_idi.R` | category-free NRI/IDI |
| B | `06_reclassification/28_nri_cv.R` | NRI/IDI under 10-fold CV |
| B | `06_reclassification/29_nri_summary.py` | NRI/IDI summaries |
| B | `07_attribution/19_part2_attribution.py` | methodological attribution |
| B | `07_attribution/26_usability_and_trend.py` | usability and time trends |
| B | `07_attribution/25d_home_final.py` | home-cohort advantage |
| B | `07_attribution/15_compare_tiers.py` | four-tier comparison |
| B | `08_cutpoint/22_cutpoint_manipulation.R` | cut-point search experiment |
| B | `09_rederivation/14_de_novo_benchmark.R` | de novo comparator |
| B | `09_rederivation/14b_de_novo_full.R` | de novo full grid |
| B | `09_rederivation/21_part3_compliant.R` | minimal-compliant signature, locked validation |
| C | `10_figures/37_make_figures.py` | manuscript figures |
| C | `10_figures/38_make_prisma.py` | PRISMA flow diagram |
| C | `99_auxiliary/35_numbers_for_paper.py` | the exact numbers quoted in the manuscript |

`99_auxiliary/` also holds supporting and superseded scripts (for example `11_part2_attribution.py`,
superseded by `19_part2_attribution.py`) retained so that every step we actually took is visible.
The manuscript-producing script (`36_write_imrad_docx.py`) is not included: it operates on the
manuscript text rather than on results.

## Known defect

See `04_evaluate/DEFECT_LOG.md`. The pre-fix snapshot of the affected script is kept as
`04_evaluate/pre_fix_backup_07_audit_signatures.R` for verification.
