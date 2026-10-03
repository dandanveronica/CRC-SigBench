# -*- coding: utf-8 -*-
"""NRI/IDI 结果深挖：签名 vs 随机基因集基线，配对比较 + 显著率对比 + 敏感性分析"""
import numpy as np, pandas as pd
from scipy import stats

o = pd.read_csv('audit_out_clean/per_cohort_NRI_IDI.csv')
r = pd.read_csv('audit_out_clean/NRI_IDI_随机基线.csv')
pc = pd.read_csv('audit_out_clean/per_cohort.csv')
pc['d_rand'] = pc.C - pc.null_median
dr = pc.groupby('sig_id').d_rand.mean()

print('='*70)
print('一、按队列：签名 NRI vs 随机基因集 NRI（同队列配对比较）')
print('='*70)
print(f"{'队列':<12}{'n_sig':>6}{'签名NRI中位':>11}{'随机NRI中位':>11}{'净NRI':>9}{'MW-p':>10}{'签名显著率':>10}{'随机显著率':>10}")
rows=[]
for c, s in o.groupby('cohort'):
    rr = r[r.cohort == c]
    if len(rr) == 0: continue
    u = stats.mannwhitneyu(s.NRI.dropna(), rr.NRI.dropna(), alternative='greater')
    sr = 100*s.NRI_p.lt(0.05).mean(); rr_ = 100*rr.NRI_p.lt(0.05).mean()
    print(f'{c:<12}{len(s):>6}{s.NRI.median():>11.4f}{rr.NRI.median():>11.4f}{s.NRI.median()-rr.NRI.median():>9.4f}{u.pvalue:>10.3g}{sr:>9.1f}%{rr_:>9.1f}%')
    rows.append((c, s.NRI.median()-rr.NRI.median(), u.pvalue))
net = np.median([x[1] for x in rows])
print(f'\n净 NRI（签名−随机）跨队列中位: {net:+.4f}')

print('\n'+'='*70)
print('二、显著率对比：用 NRI 宣称"显著改善"到底有多少是白捡的')
print('='*70)
a = int((o.NRI_p < 0.05).sum()); na = len(o.NRI_p.dropna())
b = int((r.NRI_p < 0.05).sum()); nb = len(r.NRI_p.dropna())
print(f'  已发表签名: {a}/{na} = {100*a/na:.1f}% 的 NRI 达 p<0.05')
print(f'  随机基因集: {b}/{nb} = {100*b/nb:.1f}% 的 NRI 达 p<0.05   ← 名义应为 5%')
tab = [[a, na-a],[b, nb-b]]
chi2, p, dof, _ = stats.chi2_contingency(tab)
print(f'  卡方检验 p = {p:.3g}')
print(f'  → 随机基因集的假阳性率是名义值的 {100*b/nb/5:.1f} 倍（反保守，NRI 检验在样本内失效）')
print(f'  → 签名比随机多的部分: {100*a/na - 100*b/nb:+.1f} 个百分点')

print('\n'+'='*70)
print('三、NRI 与我们主指标的关系')
print('='*70)
o2 = o.merge(pc.groupby('sig_id').d_rand.mean().rename('d_rand_sig'), on='sig_id', how='left')
for cc in ['NRI','IDI']:
    s = o2.dropna(subset=[cc,'d_rand_sig'])
    rho, pv = stats.spearmanr(s[cc], s.d_rand_sig)
    print(f'  {cc} vs ΔC(随机之上，签名级): Spearman ρ={rho:+.3f} p={pv:.3g}')
s = o2.dropna(subset=['NRI','IDI'])
rho, pv = stats.spearmanr(s.NRI, s.IDI)
print(f'  NRI vs IDI: Spearman ρ={rho:+.3f} p={pv:.3g}')

print('\n'+'='*70)
print('四、敏感性分析：剔除 control 数太少的小队列')
print('='*70)
nc = o.groupby('cohort').n_ctrl.median()
print('  各队列 control 数（超5年存活）:'); print('   ', nc.round(0).to_dict())
big = [c for c in nc.index if nc[c] >= 30]
ob = o[o.cohort.isin(big)]
print(f'\n  保留 control≥30 的队列: {sorted(big)}  (n={len(ob)} 条)')
print(f'  NRI 中位 {ob.NRI.median():+.4f}  为正 {100*(ob.NRI>0).mean():.1f}%  显著率 {100*(ob.NRI_p<0.05).mean():.1f}%')
rb = r[r.cohort.isin(big)]
print(f'  同队列随机基线 NRI 中位 {rb.NRI.median():+.4f}  显著率 {100*(rb.NRI_p<0.05).mean():.1f}%')
print(f'  净 NRI {ob.NRI.median()-rb.NRI.median():+.4f}')

print('\n'+'='*70)
print('五、IDI 汇总')
print('='*70)
print(f'  签名 IDI 中位 {o.IDI.median():+.5f}  为正 {100*(o.IDI>0).mean():.1f}%  显著率 {100*(o.IDI_p<0.05).mean():.1f}%')
print(f'  随机 IDI 中位 {r.IDI.median():+.5f}  为正 {100*(r.IDI>0).mean():.1f}%  显著率 {100*(r.IDI_p<0.05).mean():.1f}%')
print(f'  净 IDI {o.IDI.median()-r.IDI.median():+.5f}')
print('\n  参考口径：IDI 0.01 通常被视为"小但有意义"，0.005 以下基本无临床价值')
