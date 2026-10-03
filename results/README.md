# results/ — data dictionary

Every file below is a deposited snapshot produced by `run_all.sh`.
Columns are documented per file inside each header row; the "original file" column gives the
name used by the pipeline, kept so that any intermediate can be traced back on request.
Licence: CC0 1.0 Universal.

| Group | File | Original file | Description |
|---|---|---|---|
| Core matrix | `per_cohort_2560.csv` | `per_cohort_8队列.csv` | CORE: 320 签名 x 8 队列，一行一次外部验证 |
| Core matrix | `signature_level_summary.csv` | `签名级_8队列汇总.csv` | 签名级跨队列汇总 |
| Core matrix | `signature_level_usability.csv` | `签名级可用性判决.csv` | 可用性判决（C / ΔC / 校准三条件） |
| Core matrix | `gene_count_vs_deltaC.csv` | `签名级_基因数vsΔC.csv` | 基因数 vs ΔC |
| Cohorts and null | `cohort_summary.csv` | `队列汇总_8队列.csv` | 8 队列汇总 |
| Cohorts and null | `cohort_characteristics.csv` | `队列特征表.csv` | 队列基线特征 |
| Cohorts and null | `cohort_summary_with_tests.csv` | `队列级汇总_含检验.csv` | 队列级汇总含检验 |
| Cohorts and null | `null_highres_per_cohort.csv` | `highres_per_cohort.csv` | 高分辨率零分布 B=40000 逐队列 |
| Cohorts and null | `null_highres_summary.csv` | `highres_summary.csv` | 高分辨率零分布汇总 |
| Cohorts and null | `null_matched.csv` | `null_matched.csv` | 基因数匹配零分布 |
| Reclassification | `nri_idi.csv` | `NRI_IDI_队列汇总.csv` | NRI/IDI 队列汇总（样本内） |
| Reclassification | `nri_idi_cv.csv` | `NRI_IDI_CV_队列汇总.csv` | NRI/IDI 10 折交叉验证 |
| Reclassification | `nri_idi_cv_signature.csv` | `NRI_IDI_CV_签名.csv` | NRI/IDI CV 逐签名 |
| Reclassification | `nri_idi_random_baseline.csv` | `NRI_IDI_随机基线.csv` | NRI/IDI 随机基因集基线 |
| Reclassification | `nri_idi_cv_random_baseline.csv` | `NRI_IDI_CV_随机基线.csv` | NRI/IDI CV 随机基线 |
| Reclassification | `per_cohort_nri_idi.csv` | `per_cohort_NRI_IDI.csv` | 逐签名 x 队列 NRI/IDI |
| Attribution and bias | `four_tier_comparison.csv` | `四层对照.csv` | 四层对照主表 |
| Attribution and bias | `four_tier_stratified.csv` | `四层对照_分层.csv` | 四层对照分层 |
| Attribution and bias | `methodological_attribution.csv` | `Part2_归因检验.csv` | 方法学归因检验 |
| Attribution and bias | `publication_bias_paired.csv` | `原文声称值vs重算值.csv` | 声称值 vs 重算值（配对） |
| Attribution and bias | `publication_bias_pairs.csv` | `复现成功率_配对明细.csv` | 复现配对明细（103 对 AUC5 + C 对） |
| Attribution and bias | `probast_assessment.csv` | `PROBAST_自动评估.csv` | PROBAST 逐签名自动评估 |
| Attribution and bias | `methodology_extraction.csv` | `全文方法学提取.csv` | 全文方法学要素提取 |
| Attribution and bias | `manual_verification_30.csv` | `人工核对结果_30条.csv` | 30 条人工抽查核对结果 |
| Attribution and bias | `manual_verification_30.xlsx` | `人工核对30条_已核对.xlsx` | 同一份 30 条核对底稿（Excel 版，含逐项人工判定与备注） |
| Cut-point and trends | `cutpoint_summary.csv` | `截点操纵_汇总.csv` | 截点操纵实验汇总 |
| Cut-point and trends | `cutpoint_published.csv` | `截点操纵_已发表签名.csv` | 截点操纵：已发表签名 |
| Cut-point and trends | `cutpoint_random.csv` | `截点操纵_随机基因集.csv` | 截点操纵：随机基因集 |
| Cut-point and trends | `time_trends.csv` | `时间趋势.csv` | 年份趋势 |
| De novo and home advantage | `rederivation_validation.csv` | `Part3_合规签名_验证.csv` | Part3 自产合规签名验证 |
| De novo and home advantage | `rederivation_coefficients.csv` | `Part3_签名基因与系数.csv` | Part3 签名基因与系数 |
| De novo and home advantage | `part3_comparison.csv` | `Part3_对比.csv` | Part3 落差对比 |
| De novo and home advantage | `de_novo_signatures.csv` | `de_novo_自产签名.csv` | 自产签名（de novo） |
| De novo and home advantage | `de_novo_with_random_baseline.csv` | `自产签名_含随机基线.csv` | 自产签名含随机基线 |
| De novo and home advantage | `home_advantage.csv` | `回巢优势.csv` | 回巢优势主结果 |
| De novo and home advantage | `home_advantage_lenient.csv` | `回巢优势_终版_宽松.csv` | 回巢优势（宽松口径） |
| De novo and home advantage | `home_advantage_highconf.csv` | `回巢优势_终版_高置信.csv` | 回巢优势（高置信口径） |
| De novo and home advantage | `home_advantage_permutation.csv` | `回巢_置换检验.csv` | 回巢优势置换检验 |
| De novo and home advantage | `home_vs_unfamiliar.csv` | `回巢_自家vs陌生.csv` | 自家 vs 陌生队列 |
| Cross-platform | `top20_robust.csv` | `最稳健签名_Top20_8队列.csv` | 8 队列最稳健签名 Top20 |
| Cross-platform | `cross_platform_pairing.csv` | `跨平台配对_8队列.csv` | RNA-seq vs 芯片跨平台配对 |
| Cross-platform | `tcga_stratified.csv` | `TCGA分层.csv` | TCGA 分层结果 |
| Cross-platform | `training_source_tcga.csv` | `训练来源判定_含TCGA.csv` | 训练来源判定（含 TCGA） |
| Cross-platform | `clinical_baseline_cindex.csv` | `队列临床基线_合并.csv` | 各队列纯临床模型 C-index |

**Not deposited on purpose:** the eight standardised per-cohort expression objects
(`cohorts/*.rds`) and the publisher full-text XML used for claim extraction. The former are
rebuilt from public sources by `code/01_build_cohorts/`; the latter are not redistributable.
Working intermediates not listed here are available from the authors on request.
