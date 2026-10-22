# legacy/v1.0.x — 320-signature snapshot (superseded)

本目录保存 **v1.0.0–v1.0.2 的 320 库快照**，仅用于复现已发布 tag 的历史结果。
**不要**用它复现 v1.1.0 之后（含投稿）的任何数字。

| 文件 | 说明 |
|---|---|
| `signatures_320.csv` | 剔除前的 320 条签名目录（v1.1.0 起为 312 条，见 `results/submission_tables/S5_312签名清单.csv`） |
| `signature_library_320_genes.csv` | 320 库的基因构成（v1.1.0 起为 `data/signatures/signature_library_312_genes.csv`） |
| `per_cohort_2560.csv` | 320 × 8 = 2,560 次外部验证逐条结果（v1.1.0 起为 312 × 8 = 2,496，见 `results/per_cohort_2496.csv` 与 `results/submission_tables/S1_2496次评估逐条结果.csv`） |

## 为什么移到这里

v1.1.0（2026-10-18）复审后剔除 8 条签名，主口径改为 **312 签名 / 2,496 次评估**。
320 库文件留在 `data/` 与 `results/` 下会被误读为当前口径，故隔离到本目录。
文件本身未作任何改动，逐字节等同 v1.0.2。

## 仍依赖 320 口径的代码

`code/07_attribution/` 下的归因脚本（`15_compare_tiers.py`、`19_part2_attribution.py`）
读 `签名320清单.csv`，属 v1.0.x 归因链，**刻意保持 320 口径不动**；
Stage B 的 320 口径可复现性声明范围见仓库 README。
