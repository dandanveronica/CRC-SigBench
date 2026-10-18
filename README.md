# CRC-SigBench — a cross-cohort benchmark of 320 colorectal cancer prognostic transcriptomic signatures

https://github.com/dandanveronica/CRC-SigBench · OSF registration: DOI 10.17605/OSF.IO/6PF43 · Code: MIT · Data and results: CC0 1.0

> **Version note (17 October 2026).** Releases **v1.0.0–v1.0.2** are the *pre-audit snapshot*: they
> hold the 320-signature library and the 2,560-row result matrix exactly as analysed at the time of
> OSF registration (3 October 2026). A post-hoc eligibility audit (17 October 2026) re-read all 338
> candidate publications in full text and removed 8 that did not meet the eligibility criteria,
> leaving **312 signatures** and **2,496** external validations. Release **v1.1.0** deposits the
> audited artefacts — `data/signatures/signature_library_312_genes.csv` and
> `results/per_cohort_2496.csv` — next to the 320-signature snapshot; both catalogues share an
> identical schema, so the audited one can be substituted straight into stage B. No decision rule,
> threshold or cohort selection was altered by the audit, and the audited matrices are also supplied
> as Additional files 1–3 of the companion article.
>
> **Pipeline note.** The audited matrices come from a pipeline that additionally restricts every
> cohort to primary tumour samples (19 mucosa samples excluded from GSE39582, 160 from GSE87211).
> That restriction is documented in Section 2.2 of the article's supplementary material but is
> *not yet* in the code at v1.0.2, so re-running stage B on the 312-signature catalogue will not
> reproduce `per_cohort_2496.csv` exactly. The pipeline upgrade and an end-to-end re-run
> verification are scheduled for **v1.2.0**.

## What this is

A reproducible benchmark that re-evaluates **320** published protein-coding mRNA prognostic
signatures for colorectal cancer across **8** fully independent cohorts
(7 GEO microarray cohorts + TCGA COADREAD RNA-seq), giving **2,560** independent external
validations.

Each signature is reconstructed strictly from its published gene list and coefficients, scored
within-cohort after gene-wise z-standardisation, and benchmarked against:

1. a **size-matched random-gene-set null distribution** (B = 1,000; a high-resolution B = 40,000 run is also provided),
2. a **Cox model of routine clinical variables** (stage, age, sex, tumour site),
3. a **de novo comparator set** built by us under a pre-specified minimal-compliant protocol.

## Repository layout

```
CRC-SigBench/
├── README.md
├── LICENSE                       # MIT (code)
├── .gitignore
├── run_all.sh                    # top-level driver (see "Reproducibility scope" below)
├── code/
│   ├── 00_check_env.R            # prints R version, verifies required packages
│   ├── 01_build_cohorts/         # GEO download, TCGA via Xena/cBioPortal, normalisation, clinical harmonisation
│   ├── 02_signatures/            # parse published gene lists + coefficients -> machine-readable catalogue
│   ├── 03_score/                 # signature reconstruction, z-scoring, risk score
│   ├── 04_evaluate/              # C-index, 5-year AUC, calibration slope, clinical benchmark
│   │   └── DEFECT_LOG.md         # the row-misalignment defect we found and fixed
│   ├── 05_null/                  # random-gene-set null distribution (B = 1,000 / 40,000; size-matched)
│   ├── 06_reclassification/      # category-free NRI/IDI, in-sample and 10-fold CV
│   ├── 07_attribution/           # methodological attribution, time trends, home-cohort advantage
│   ├── 08_cutpoint/              # cut-point search experiment
│   ├── 09_rederivation/          # Part 3: minimal-compliant de novo signature, locked validation
│   ├── 10_figures/               # manuscript figures
│   ├── 11_search/                # PRISMA 2020 search strategy for the signature corpus
│   └── 99_auxiliary/             # supporting and superseded scripts retained for transparency
├── data/
│   └── signatures/               # 320 (v1.0.x) + 312 (v1.1.0) signatures: gene lists, metadata, claimed values
└── results/
    ├── README.md                 # data dictionary for every result file
    ├── per_cohort_2560.csv       # ONE ROW PER SIGNATURE x COHORT — pre-audit core deliverable (v1.0.x)
    └── per_cohort_2496.csv       # post-audit core deliverable, 312 signatures (v1.1.0)
```

### The signature catalogue (what gets scored)

`data/signatures/signature_library_320_genes.csv` is the machine-readable catalogue that stage B
consumes: one row per audited signature with `sig_id`, year, journal, title, PMCID, source URL,
gene count and the **actual gene list** (semicolon-separated). All 320 signatures are present and
all 320 carry a non-empty gene list (7.8 genes on average).

Release v1.1.0 adds `data/signatures/signature_library_312_genes.csv` — the same schema, restricted
to the 312 signatures that survived the post-hoc eligibility audit of 17 October 2026 (the 8 removed
signatures are listed, with reasons, in Table A8 of the article's supplementary material).

> **Scoring note (honest disclosure).** The catalogue contains **gene lists only** — no risk-score
> coefficients. `code/03_score/` therefore scores every signature as the *unweighted* mean of
> member-gene z-scores (`coef = NULL` in `07_audit_signatures.R`). Published coefficient sets were
> not systematically extracted from the source papers, so for signatures that were originally
> coefficient-weighted our reconstruction is a standardised approximation rather than a verbatim
> reimplementation. The direction of the resulting bias is not established and is disclosed as a
> limitation in the manuscript.

`data/cohorts/` is intentionally **not** populated. The eight standardised per-cohort objects are
rebuilt locally by `code/01_build_cohorts/`; they are not redistributed because the underlying GEO
and TCGA data remain subject to their own terms of use.

## Reproducing the results

```bash
git clone [repository URL]
cd CRC-SigBench
Rscript code/00_check_env.R     # prints R version + verifies required packages
bash run_all.sh                 # end-to-end; writes everything under results/
```

Requirements: R ≥ 4.4 (survival, glmnet, data.table, parallel, timeROC);
Python ≥ 3.9 (pandas, numpy, scipy, matplotlib, openpyxl, python-docx, requests).
Input cohorts are downloaded automatically from GEO, UCSC Xena and cBioPortal; no private data
are required.

### Reproducibility scope (read before judging us)

`run_all.sh` is organised in three stages, and only stage B is fully self-contained:

- **Stage A — literature extraction (not reproducible from this repository).**
  `code/02_signatures/10_fetch_fulltext.py` and `18_claim_extract.py` operated on publisher
  full-text XML obtained under our institutional subscriptions. Those files are **not
  redistributable**, so this stage cannot be re-run by a third party. What we deposit instead is
  its *output*: `data/signatures/claimed_values_extracted.csv` and
  `claimed_values_manual30.csv`. Every published performance value we compare against is in
  those two files, with the source sentence and the manual-verification verdict.
- **Stage B — cohort construction and evaluation (reproducible).**
  From the deposited signature catalogue and the public cohorts, every number in the manuscript
  is regenerated by script, with no manual step. This is where the 2,560-row result matrix,
  the null distributions, NRI/IDI, calibration and the de novo comparator are produced.
- **Stage C — figures and numbers-for-paper (reproducible from `results/`).**
  `code/10_figures/` and `code/99_auxiliary/35_numbers_for_paper.py` read the deposited result
  tables and emit the manuscript figures and the exact figures quoted in the text.

So: every *computed* result in the paper is reproducible from public inputs by script; the one
step that cannot be re-run by others is the extraction of published claims from paywalled full
texts, and for that step we deposit the extracted table rather than the inputs.

## The core result

After Benjamini–Hochberg correction across all 2,560 tests, **no individual signature was
significantly better than a size-matched random gene set** (minimum adjusted p = 0.077).
Reported 5-year AUC exceeded our recomputation by a median of **+0.176** across 103 matched
pairs (99% in the same direction); **7.8%** of pairs reproduced within ±0.05.

### Post-audit (v1.1.0, 312 signatures / 2,496 evaluations)

Removing the 8 ineligible publications leaves the conclusion unchanged and slightly stronger.
After Benjamini–Hochberg correction across all 2,496 tests, **no signature was significantly
better than the size-matched random null** (minimum adjusted p = 0.071). The median signature
reached C = 0.555 and a median ΔC_random of **+0.0107** (+0.0045 over routine clinical variables);
299 of 2,496 evaluations were nominally positive and none survived correction. Median calibration
slope across signatures was 0.363, and **3 of 312** signatures fell inside the usable region
(C ≥ 0.60 with slope ≥ 0.70) — 2 of them under the joint criterion locked in the OSF registration.
The claimed-versus-recomputed AUC comparison covers **100** matched pairs (median difference
+0.161; 99.0% in the same direction; 10.0% reproduced within ±0.05).

## Known and corrected defect

An early version of the scoring code contained a silent row-misalignment defect: `coxph` drops
rows with missing covariates, and `predict()` returned a short vector that R recycled without
warning, misaligning predictions with survival times. It affected ΔC-index in the two cohorts
with missing stage (GSE39582: 1 sample; TCGA: 11 samples) and nothing else. It was found and
fixed before any result was reported; the fix, its numerical impact and the pre-fix script
snapshot are documented in `code/04_evaluate/DEFECT_LOG.md` and reported in full in the
supplementary material of the companion article (Section 2.11).

## Citation

[Authors, title, journal, year — TO BE COMPLETED]
