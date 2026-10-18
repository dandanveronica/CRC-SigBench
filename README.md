# CRC-SigBench — a cross-cohort benchmark of colorectal cancer prognostic transcriptomic signatures

**Current release: v1.1.0 (17 October 2026) — 312 signatures, 2,496 external validations, after the
post-hoc eligibility audit. The 320-signature pre-audit snapshot analysed at the time of OSF
registration (3 October 2026) is preserved unchanged at v1.0.0–v1.0.2 — see the Version note below.**

https://github.com/dandanveronica/CRC-SigBench · OSF registration: DOI 10.17605/OSF.IO/6PF43 · Code: MIT · Data and results: CC0 1.0

> **Version note (17 October 2026).** Releases **v1.0.0–v1.0.2** are the *pre-audit snapshot*: they
> hold the 320-signature library and the 2,560-row result matrix exactly as analysed at the time of
> OSF registration (3 October 2026). A post-hoc eligibility audit (17 October 2026) re-examined all
> 320 publications in that library — every title and abstract, with full-text verification for every
> flagged or borderline item — against the four pre-specified eligibility criteria, and removed 8
> that did not meet them (listed, with reasons, in Table A8 of the article's supplementary
> material), leaving **312 signatures** and **2,496** external validations. Together with the 18
> candidates removed at standardisation (13 with ≤2 genes, 1 treatment-response outcome,
> 4 non-reproducible gene lists), the 338 primary-analysis candidates were reduced to the 312
> audited signatures. Release **v1.1.0** deposits the audited artefacts —
> `data/signatures/signature_library_312_genes.csv` and `results/per_cohort_2496.csv` — next to the
> 320-signature snapshot; both catalogues share an identical schema, so the audited one can be
> substituted straight into stage B. No decision rule, threshold or cohort selection was altered by
> the audit, and the audited catalogue and matrices are also supplied as Additional files 1–3 of the
> companion article (S1 result matrix, S5 signature catalogue, S8 gene list).
>
> **Pipeline note.** The audited matrices come from a pipeline that additionally restricts every
> cohort to primary tumour samples (19 mucosa samples excluded from GSE39582, 160 from GSE87211).
> That restriction is documented in Section 2.2 of the article's supplementary material but is
> *not yet* in the code at v1.0.2, so re-running stage B on the 312-signature catalogue will not
> reproduce `per_cohort_2496.csv` exactly. The pipeline upgrade and an end-to-end re-run
> verification are scheduled for **v1.2.0**.

## What this is

A reproducible benchmark that re-evaluates published protein-coding mRNA prognostic signatures for
colorectal cancer across **8** fully independent cohorts (7 GEO microarray cohorts + TCGA COADREAD
RNA-seq). The **current audited release (v1.1.0)** covers **312 signatures** and gives **2,496**
signature × cohort evaluations; the **pre-audit v1.0.x snapshot** covers 320 signatures and 2,560
evaluations. Both are kept side by side and share an identical schema.

Each signature is reconstructed from its published gene list under one uniform scoring rule — the
unweighted mean of member-gene z-scores within cohort, orientation chosen by the higher of the two
concordance directions — and benchmarked against:

1. a **size-matched random-gene-set null distribution** (B = 1,000; a high-resolution B = 40,000 run is also provided),
2. a **Cox model of routine clinical variables** (stage, age, sex, tumour site),
3. a **de novo comparator set** built by us under a pre-specified minimal-compliant protocol.

Published risk-score coefficients were deliberately not carried into the evaluation, so that
published signatures, random gene sets and the clinical comparator are all scored on one footing;
the approximation this imposes on signatures that were originally coefficient-weighted is disclosed
as a limitation in the companion article.

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
│   └── signatures/               # 312 (v1.1.0) + 320 (v1.0.x) signatures: gene lists, metadata, claimed values
└── results/
    ├── README.md                 # data dictionary for every result file
    ├── per_cohort_2560.csv       # ONE ROW PER SIGNATURE x COHORT — pre-audit core deliverable (v1.0.x)
    └── per_cohort_2496.csv       # post-audit core deliverable, 312 signatures (v1.1.0)
```

### The signature catalogue (what gets scored)

`data/signatures/signature_library_312_genes.csv` (v1.1.0, current) and
`data/signatures/signature_library_320_genes.csv` (v1.0.x, pre-audit snapshot) are the
machine-readable catalogues that stage B consumes: one row per audited signature with `sig_id`,
year, journal, title, PMCID, source URL, gene count and the **actual gene list**
(semicolon-separated). Every signature is present and every signature carries a non-empty gene list
(7.8 genes on average in the 320-signature snapshot, 7.7 in the audited 312).

The v1.1.0 catalogue has the same schema, restricted to the 312 signatures that survived the
post-hoc eligibility audit of 17 October 2026 (the 8 removed signatures are listed, with reasons, in
Table A8 of the article's supplementary material).

> **Scoring note (honest disclosure).** The catalogue contains **gene lists only** — no risk-score
> coefficients. `code/03_score/` therefore scores every signature as the *unweighted* mean of
> member-gene z-scores (`coef = NULL` in `07_audit_signatures.R`). Published regression coefficients
> were deliberately not carried into the evaluation, because transporting each model's own weights
> would test every signature under a different rule and make the 312 mutually incomparable; for
> signatures that were originally coefficient-weighted our reconstruction is therefore a standardised
> approximation rather than a verbatim reimplementation. The direction of that error is unknown and
> is stated as a limitation (Section 4) of the companion article.

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
- **Stage B — cohort construction and evaluation (reproducible for the v1.0.x analysis).**
  From the deposited signature catalogue and the public cohorts, every number of the v1.0.x
  (320-signature / 2,560-evaluation) analysis is regenerated by script, with no manual step. This
  is where the 2,560-row result matrix, the null distributions, NRI/IDI, calibration and the de
  novo comparator are produced. The post-audit v1.1.0 matrices additionally apply the
  primary-tumour restriction, which is documented in the companion article but not yet in the code
  — see the Pipeline note above; an exact re-run is scheduled for v1.2.0.
- **Stage C — figures and numbers-for-paper (reproducible from `results/`).**
  `code/10_figures/` and `code/99_auxiliary/35_numbers_for_paper.py` read the deposited result
  tables and emit the manuscript figures and the exact figures quoted in the text.

So: every *computed* result of the v1.0.x analysis is reproducible from public inputs by script;
the one step that cannot be re-run by others is the extraction of published claims from paywalled
full texts, and for that step we deposit the extracted table rather than the inputs. For the
post-audit v1.1.0 matrices, the honest position is the Pipeline note: artefacts deposited,
end-to-end code verification scheduled for v1.2.0.

## The core result

### Pre-audit (v1.0.x: 320 signatures / 2,560 evaluations)

After Benjamini–Hochberg correction across all 2,560 tests, **no individual signature was
significantly better than a size-matched random gene set** (minimum adjusted p = 0.077).
Reported 5-year AUC exceeded our recomputation by a median of **+0.176** across 103 matched
pairs (99% in the same direction); **7.8%** of pairs reproduced within ±0.05.

### Post-audit (v1.1.0, 312 signatures / 2,496 evaluations)

Removing the 8 ineligible publications changed no qualitative conclusion.
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

Yang D, Chen P, Wang J, Zhao L, Xiang T, Wu Z. Colorectal cancer prognostic signatures carry a
real but unusable signal: a meta-research re-evaluation of 312 published signatures across 2,496
external validations. *Journal of Clinical Epidemiology* — manuscript in submission, October 2026.
