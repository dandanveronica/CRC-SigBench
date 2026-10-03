# Supplementary Table S9. Search strategy and reproducible query

- **Database**: Europe PMC — sole literature source for signature identification
- **Interface**: Europe PMC REST API
- **Date of search**: Executed in September 2026; record set locked 29 September 2026
- **Date of re-execution for this table**: 3 October 2026
- **Date limits**: First publication date 1 January 2020 to search date
- **Language limit**: English-language full text required; applied at full-text assessment, not in the query
- **Records returned at re-execution**: 6086
- **Coverage check**: 332/338 (98.2%) of the pre-standardisation library retrieved directly
- **Screening register**: 1,241 records in an auditable register, released with the benchmark

## Query string

```
(TITLE_ABS:"colorectal" OR TITLE_ABS:"colorectal cancer" OR TITLE_ABS:"colon cancer" OR TITLE_ABS:"rectal cancer" OR TITLE_ABS:"colon adenocarcinoma" OR TITLE_ABS:"rectal adenocarcinoma" OR TITLE_ABS:"colorectal carcinoma" OR TITLE_ABS:"CRC" OR TITLE_ABS:"colon" OR TITLE_ABS:"rectal" OR TITLE_ABS:"colonic") AND (TITLE_ABS:"signature" OR TITLE_ABS:"signatures" OR TITLE_ABS:"gene signature" OR TITLE_ABS:"prognostic signature" OR TITLE_ABS:"prognostic model" OR TITLE_ABS:"risk score" OR TITLE_ABS:"risk model" OR TITLE_ABS:"nomogram" OR TITLE_ABS:"multigene" OR TITLE_ABS:"multi-gene" OR TITLE_ABS:"gene set" OR TITLE_ABS:"gene panel") AND (TITLE_ABS:"survival" OR TITLE_ABS:"prognosis" OR TITLE_ABS:"prognostic" OR TITLE_ABS:"overall survival" OR TITLE_ABS:"recurrence" OR TITLE_ABS:"outcome" OR TITLE_ABS:"predict" OR TITLE_ABS:"biomarker") AND (FIRST_PDATE:[2020-01-01 TO 2026-12-31])
```
