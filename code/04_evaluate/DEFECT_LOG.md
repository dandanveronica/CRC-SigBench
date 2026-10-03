# DEFECT LOG — silent row misalignment in ΔC-index

**Status:** found and fixed before any result was reported. Reported in full in the manuscript
(Section 2.11). The pre-fix script snapshot is kept alongside this file as
`pre_fix_backup_07_audit_signatures.R`.

## What happened

In `07_audit_signatures.R`, the ΔC-index routine fitted a Cox model and predicted linear
predictors like this:

```r
cb   <- coxph(f_base, data = d)
lp_b <- predict(cb, type = "lp")
cindex(d$.t, d$.s, lp_b)
```

`coxph()` silently drops rows with missing covariates, so `predict()` returned a vector *shorter*
than `nrow(d)`. Inside `cindex()`, when `pred[ok]` and `time[ok]` differ in length, R **recycles
the shorter vector without any warning**, producing a row-by-row misalignment: predictions were
compared against the wrong patients' survival times.

## How it was found

The eighth cohort (TCGA, RNA-seq) was added last. On TCGA, ΔC_clinical was +0.00017 with only
53.8% of signatures positive — an apparent finding ("signatures add nothing over routine clinical
variables on RNA-seq"). It was an artefact: TCGA is one of the two cohorts with missing stage.

## The fix

Subset to complete cases **before** fitting, so that time, status and linear predictor are
aligned by construction rather than repaired afterwards:

```r
d_cc <- d[complete.cases(d[, covars]), ]
cb   <- coxph(f_base, data = d_cc)
lp_b <- predict(cb, type = "lp")
cindex(d_cc$.t, d_cc$.s, lp_b)
```

## Impact, measured exactly

Only the cohorts with missing covariates changed:

| Cohort | ΔC_clinical before | after | difference | missing |
|---|---|---|---|---|
| GSE17536 | +0.00357 | +0.00357 | 0.00000 | – |
| GSE17537 | +0.00568 | +0.00568 | 0.00000 | – |
| GSE29621 | +0.00992 | +0.00992 | 0.00000 | – |
| GSE38832 | +0.00868 | +0.00868 | 0.00000 | – |
| GSE39582 | +0.00086 | +0.00133 | +0.00037 | 1 sample without stage |
| GSE72970 | +0.03164 | +0.03164 | 0.00000 | – |
| GSE87211 | +0.00206 | +0.00206 | 0.00000 | – |
| TCGA COADREAD | +0.00017 | +0.00179 | +0.00174 | 11 samples without stage |

Six cohorts changed by exactly 0.00000, which is what confirms the localisation of the defect.
C-index, 5-year AUC and calibration slope were unaffected to six decimal places; the null
distribution moved only within random-seed noise. NRI/IDI was unaffected because `27_nri_idi.R`
already subset on complete cases before fitting.

All results reported in the manuscript were computed after the fix.
