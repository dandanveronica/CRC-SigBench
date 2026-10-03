# -*- coding: utf-8 -*-
"""① 临床可用性判决：有没有任何签名过线  ② 时间趋势：2018→2026 这个领域有没有进步"""
import numpy as np, pandas as pd
from scipy import stats

pc = pd.read_csv('audit_out_8/per_cohort.csv')     # 8 队列（含 TCGA），delta_C 已修复
KEEP = set(pd.read_csv('audit_out_clean/签名320清单.csv').sig_id)
pc = pc[pc.sig_id.isin(KEEP)]                      # 统一到已定稿的 320 条口径
meta = pd.read_csv('audit_out_clean/签名元数据_归因变量.csv')
claim = pd.read_csv('audit_out_clean/原文声称值vs重算值.csv')

pc['d_rand'] = pc.C - pc.null_median
S = pc.groupby('sig_id').agg(
    n_cohort=('C','size'), C_mean=('C','mean'), C_med=('C','median'), C_min=('C','min'),
    d_rand_mean=('d_rand','mean'), d_rand_med=('d_rand','median'),
    d_clin_mean=('delta_C','mean'),
    calib=('calib_slope','mean'),
    n_pos_rand=('d_rand', lambda x:(x>0).sum())).reset_index()
ft = pd.read_csv('audit_out_clean/全文方法学提取.csv')[['sig_id','是否外部验证','是否给出风险公式','是否多因素校正']]
S = S.merge(meta[['sig_id','year','journal','n_gene','方法_单细胞/空间','方法_LASSO']], on='sig_id', how='left')
S = S.merge(ft, on='sig_id', how='left')
S = S.merge(claim[['sig_id','偏倚_AUC5','偏倚_C']], on='sig_id', how='left')

print('=== ① 临床可用性判决（签名级，跨 8 队列聚合，n=%d）===' % len(S))
print(f'  C 均值分布:  中位 {S.C_mean.median():.3f}  最高 {S.C_mean.max():.3f}  ≥0.60 的 {100*(S.C_mean>=0.60).mean():.1f}%  ≥0.65 的 {100*(S.C_mean>=0.65).mean():.1f}%')
print(f'  ΔC(随机之上) 均值: 中位 {S.d_rand_mean.median():+.4f}  最高 {S.d_rand_mean.max():+.4f}')
print(f'  校准斜率均值: 中位 {S.calib.median():.3f}  ≥0.7 的 {100*(S.calib>=0.7).mean():.1f}%  （理想=1）')
print()
gates = [
    ('门槛A 严格  C均值≥0.65 且 ΔC_rand≥0.03', (S.C_mean>=0.65) & (S.d_rand_mean>=0.03)),
    ('门槛B 中等  C均值≥0.60 且 ΔC_rand≥0.02', (S.C_mean>=0.60) & (S.d_rand_mean>=0.02)),
    ('门槛C 宽松  C均值≥0.58 且 ΔC_rand>0',    (S.C_mean>=0.58) & (S.d_rand_mean>0)),
    ('门槛D 一致  ΔC_rand>0 且 ≥5/8 队列为正',  (S.d_rand_mean>0) & (S.n_pos_rand>=5)),
    ('门槛E 校准  ΔC_rand>0 且 校准斜率≥0.7',   (S.d_rand_mean>0) & (S.calib>=0.7)),
]
for lab, m in gates:
    m = m.fillna(False)
    print(f'  {lab}: 通过 {int(m.sum())} / {len(S)}  ({100*m.mean():.1f}%)')
best = S.nlargest(8, 'C_mean')[['sig_id','year','n_gene','C_mean','d_rand_mean','d_clin_mean','calib','n_pos_rand']]
print('\n  C 均值最高的 8 个签名:')
print(best.round(4).to_string(index=False))
print('\n  注：ΔC_rand 最高者 与 C 最高者 是否同一批 →', 
      '是' if set(S.nlargest(8,'C_mean').sig_id) & set(S.nlargest(8,'d_rand_mean').sig_id) else '否（两个榜单几乎不重合）')

print('\n=== ② 时间趋势 2018→2026 ===')
S['yr'] = pd.to_numeric(S.year, errors='coerce')
S0 = S.dropna(subset=['yr'])
S0 = S0[S0.yr >= 2018]
g = S0.groupby(S0.yr.astype(int))
tab = g.agg(n=('sig_id','size'), C_mean=('C_mean','median'), d_rand=('d_rand_mean','median'),
            d_clin=('d_clin_mean','median'), calib=('calib','median'),
            外部验证率=('是否外部验证','mean'), 公开公式率=('是否给出风险公式','mean'),
            单细胞率=('方法_单细胞/空间','mean'),
            LASSO率=('方法_LASSO','mean'), 偏倚AUC=('偏倚_AUC5','median'),
            基因数=('n_gene','median')).round(4)
print(tab.to_string())
for col in ['d_rand_mean','d_clin_mean','C_mean','偏倚_AUC5','是否外部验证','n_gene','calib']:
    sub = S.dropna(subset=['yr', col])
    if len(sub) > 20:
        r, p = stats.spearmanr(sub.yr, sub[col])
        print(f'  年份 vs {col}: Spearman ρ={r:+.3f}  p={p:.3g}')
S.to_csv('audit_out_clean/签名级可用性判决.csv', index=False, encoding='utf-8-sig')
tab.to_csv('audit_out_clean/时间趋势.csv', encoding='utf-8-sig')
print('\n已存: 签名级可用性判决.csv / 时间趋势.csv')
