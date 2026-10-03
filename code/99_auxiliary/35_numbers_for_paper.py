# -*- coding: utf-8 -*-
"""一次性提取 IMRaD 全文所需的全部数字，避免凭记忆写错。"""
import pandas as pd, numpy as np, json
from scipy import stats

B='audit_out_clean'
KEEP = set(pd.read_csv(f'{B}/签名320清单.csv').sig_id)
pc = pd.read_csv('audit_out_8/per_cohort.csv')
pc = pc[pc.sig_id.isin(KEEP)].copy()
pc['d_rand'] = pc.C - pc.null_median
OUT={}
def p(*a): print(*a)

p('='*72); p('【1】主结果 · 记录级（320 签名 × 8 队列 = %d 条）'%len(pc)); p('='*72)
for lab,s in [('全部 8 队列',pc),('仅 7 个 GEO 芯片',pc[pc.cohort!='TCGA_COADREAD']),('仅 TCGA RNA-seq',pc[pc.cohort=='TCGA_COADREAD'])]:
    w=stats.wilcoxon(s.d_rand.dropna()); w2=stats.wilcoxon(s.delta_C.dropna())
    p(f'{lab}: ΔC_rand 中位 {s.d_rand.median():+.4f} ({100*(s.d_rand>0).mean():.1f}%>0) p={w.pvalue:.3g} | '
      f'ΔC_clin 中位 {s.delta_C.median():+.4f} ({100*(s.delta_C>0).mean():.1f}%>0) p={w2.pvalue:.3g} | C中位 {s.C.median():.4f}')
a=pc[pc.cohort!='TCGA_COADREAD'].d_rand.dropna(); b=pc[pc.cohort=='TCGA_COADREAD'].d_rand.dropna()
p(f'RNA-seq vs 芯片 ΔC_rand: MW p={stats.mannwhitneyu(b,a,alternative="greater").pvalue:.3g}  倍率 {b.median()/a.median():.2f}x')
a2=pc[pc.cohort!='TCGA_COADREAD'].delta_C.dropna(); b2=pc[pc.cohort=='TCGA_COADREAD'].delta_C.dropna()
p(f'RNA-seq vs 芯片 ΔC_clin: MW p={stats.mannwhitneyu(b2,a2,alternative="less").pvalue:.3g}  倍率 {b2.median()/a2.median():.2f}x')

# 逐签名折叠
S = pc.groupby('sig_id').agg(C=('C','mean'),drand=('d_rand','mean'),dclin=('delta_C','mean'),
                             calib=('calib_slope','mean'),auc5=('auc5','mean'),
                             npos=('d_rand',lambda x:(x>0).sum()),
                             pmin=('p_adj','min')).reset_index()
p(f'\n签名级聚合 (n={len(S)}): C均值中位 {S.C.median():.4f} [IQR {S.C.quantile(.25):.3f}-{S.C.quantile(.75):.3f}]')
p(f'  ΔC_rand 签名级均值中位 {S.drand.median():+.4f}, {(S.drand>0).mean()*100:.1f}% 为正, p={stats.wilcoxon(S.drand).pvalue:.3g}')
p(f'  ΔC_clin 签名级均值中位 {S.dclin.median():+.4f}, {(S.dclin>0).mean()*100:.1f}% 为正, p={stats.wilcoxon(S.dclin).pvalue:.3g}')
p(f'  校准斜率均值中位 {S.calib.median():.4f} [IQR {S.calib.quantile(.25):.3f}-{S.calib.quantile(.75):.3f}], >=0.7 占 {100*(S.calib>=0.7).mean():.1f}%')
p(f'  5年AUC 均值中位 {S.auc5.median():.4f}' if S.auc5.notna().any() else '  无 auc5')
p(f'  至少 1 队列 p_adj<0.05: {(S.pmin<0.05).sum()} 个签名; >=5/8 队列 ΔC_rand>0: {(S.npos>=5).sum()} 个')
p(f'  全部 8 队列 ΔC_rand>0: {(S.npos==8).sum()} 个')

p('\n'+'='*72); p('【2】四层对照 / 自产签名'); p('='*72)
try:
    t=pd.read_csv(f'{B}/四层对照.csv'); p(t.to_string(index=False))
    t2=pd.read_csv(f'{B}/四层对照_分层.csv'); p('\n按事件数分层:'); p(t2.to_string(index=False))
except Exception as e: p('ERR',e)

p('\n'+'='*72); p('【3】Part 2 归因 + 发表偏倚'); p('='*72)
try:
    n=pd.read_csv(f'{B}/Part2_归因检验.csv'); p(n.to_string(index=False))
except Exception as e: p('ERR',e)

p('\n'+'='*72); p('【4】Part 3 合规签名'); p('='*72)
p(pd.read_csv(f'{B}/Part3_合规签名_验证.csv').to_string(index=False))
p('\n'+pd.read_csv(f'{B}/Part3_对比.csv').to_string(index=False))

p('\n'+'='*72); p('【5】NRI / IDI'); p('='*72)
p(pd.read_csv(f'{B}/NRI_IDI_队列汇总.csv').to_string(index=False))
p('\nCV:'); p(pd.read_csv(f'{B}/NRI_IDI_CV_队列汇总.csv').to_string(index=False))

p('\n'+'='*72); p('【6】可用性判决（8队列版）'); p('='*72)
U=pd.read_csv(f'{B}/签名级可用性判决.csv')
p(f'C均值>=0.65: {(U.C_mean>=0.65).sum()}  | >=0.60: {(U.C_mean>=0.60).sum()} ({100*(U.C_mean>=0.60).mean():.1f}%)  | 校准>=0.7: {(U.calib>=0.7).sum()}')
p('\nTop6 by C_mean:'); p(U.nlargest(6,'C_mean')[['sig_id','year','n_gene','C_mean','d_rand_mean','d_clin_mean','calib','n_pos_rand']].round(4).to_string(index=False))
p('\n最高校准:'); p(U.nlargest(4,'calib')[['sig_id','n_gene','C_mean','calib','d_rand_mean']].round(4).to_string(index=False))

p('\n'+'='*72); p('【7】跨平台一致性'); p('='*72)
X=pd.read_csv(f'{B}/跨平台配对_8队列.csv')
p(X.head(3).to_string(index=False)) if len(X.columns)>2 else p(X.to_string(index=False))

p('\n'+'='*72); p('【8】回巢优势（8队列版）'); p('='*72)
for f in ['回巢优势_终版_宽松.csv','回巢优势_终版_高置信.csv']:
    H=pd.read_csv(f'{B}/{f}')
    good=H.dropna(subset=['d_home','d_away'])
    p(f'{f}: n={len(H)}  Δhome {good.d_home.median():+.4f} ({100*(good.d_home>0).mean():.1f}%)  Δaway {good.d_away.median():+.4f} ({100*(good.d_away>0).mean():.1f}%)  '
      f'配对差 {good.delta.median():+.4f} p={stats.wilcoxon(good.delta).pvalue:.3g}')
