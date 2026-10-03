# -*- coding: utf-8 -*-
"""回巢优势 · 终版
① 置换检验改用「队列难度已校正」的量 d_rand = C − null_median（随机基线已吸收队列难度）
② 高置信判定子集做敏感性分析（自动判定误分类会把效应稀释向 0，高置信子集应更强）
"""
import os, re, glob
import numpy as np, pandas as pd
from scipy import stats

OUR = ['GSE39582','GSE17536','GSE17537','GSE29621','GSE38832','GSE72970','GSE87211']
OURSET = set(OUR)
GEO = re.compile(r'GSE\d{4,6}')
TW = re.compile(r'\b(train|training|discovery|derivation|deriving|construction|constructed|developed)\b', re.I)
TW_TIGHT = re.compile(r'\b(train|training|discovery|derivation)\b', re.I)
VW = re.compile(r'\b(valid|validation|testing|external|verified)\b', re.I)

def strip(fn): return re.sub(r'<[^>]+>', ' ', open(fn, encoding='utf-8', errors='ignore').read())

def near(win, rel, pat):
    p = [x.start() for x in pat.finditer(win)]
    return min([abs(x - rel) for x in p], default=10**9)

def judge(txt):
    loose, tight = set(), set()
    for m in GEO.finditer(txt):
        g, pos = m.group(), m.start()
        lo, hi = max(0, pos-120), min(len(txt), pos+120)
        win = txt[lo:hi]; rel = pos - lo
        dt, dv = near(win, rel, TW), near(win, rel, VW)
        if dt < dv: loose.add(g)
        dtt = near(win, rel, TW_TIGHT)
        if dtt < 40 and dtt < dv: tight.add(g)      # 40 字符内紧邻 train 词
    return loose, tight

meta = pd.read_csv('audit_out_clean/全文方法学提取.csv')[['sig_id','pmcid']]
rows = []
for fn in sorted(glob.glob('fulltext_xml/*')):
    pmc = re.sub(r'\.(xml|txt|html)$', '', os.path.basename(fn))
    t = strip(fn)
    lo, ti = judge(t)
    rows.append(dict(pmc=pmc, home_loose=';'.join(sorted(g for g in lo if g in OURSET)),
                     home_tight=';'.join(sorted(g for g in ti if g in OURSET))))
J = pd.DataFrame(rows)
def norm(s): return re.sub(r'[^A-Z0-9]','',str(s)).upper()
J['k'] = J.pmc.map(norm); meta['k'] = meta.pmcid.map(norm)
M = meta.merge(J.drop(columns=['pmc']), on='k', how='left')
M[['home_loose','home_tight']] = M[['home_loose','home_tight']].fillna('')
print('判定情况：宽松集', (M.home_loose!='').sum(), ' 高置信集', (M.home_tight!='').sum(), ' 总计', len(M))

pc = pd.read_csv('audit_out_clean/per_cohort.csv')
pc['d_rand'] = pc.C - pc.null_median
rng = np.random.default_rng(11)

def run(col, label):
    hm = {}
    for _, r in M.iterrows():
        hs = [g for g in str(r[col]).split(';') if g]
        if hs: hm[r.sig_id] = hs
    out = []
    for sid, hs in hm.items():
        s = pc[pc.sig_id == sid]
        h = s[s.cohort.isin(hs)]; a = s[~s.cohort.isin(hs)]
        if len(h) == 0 or len(a) < 2: continue
        dh, da = h.d_rand.mean(), a.d_rand.median()
        av = a.d_rand.values
        sim = np.array([av[rng.integers(0, len(av), len(h))].mean() for _ in range(500)])
        out.append(dict(sig_id=sid, d_home=dh, d_away=da, delta=dh-da,
                        z=(dh - sim.mean())/sim.std(),
                        p_perm=(np.sum(sim >= dh)+1)/501))
    O = pd.DataFrame(out)
    print(f'\n===== {label}（n={len(O)}）=====')
    print(f'  自家队列 d_rand 中位 {O.d_home.median():+.4f}  {100*(O.d_home>0).mean():.1f}%>0  p={stats.wilcoxon(O.d_home.dropna()).pvalue:.3g}')
    print(f'  陌生队列 d_rand 中位 {O.d_away.median():+.4f}  {100*(O.d_away>0).mean():.1f}%>0  p={stats.wilcoxon(O.d_away.dropna()).pvalue:.3g}')
    print(f'  回巢优势（配对差）中位 {O.delta.median():+.4f}  {100*(O.delta>0).mean():.1f}%>0  p={stats.wilcoxon(O.delta.dropna()).pvalue:.3g}')
    print(f'  置换 z 中位 {O.z.median():+.3f}  z>0 占 {100*(O.z>0).mean():.1f}%  p={stats.wilcoxon(O.z.dropna()).pvalue:.3g}')
    print(f'  p_perm<0.05: {(O.p_perm<0.05).sum()}/{len(O)}')
    print(f'  衰减倍数: 自家/陌生 = {O.d_home.median()/O.d_away.median():.2f}×' if O.d_away.median() != 0 else '')
    return O

O1 = run('home_loose', '宽松判定（全量）')
O2 = run('home_tight', '高置信判定（train 词 40 字符内紧邻）')
O1.to_csv('audit_out_clean/回巢优势_终版_宽松.csv', index=False, encoding='utf-8-sig')
O2.to_csv('audit_out_clean/回巢优势_终版_高置信.csv', index=False, encoding='utf-8-sig')
print('\n已存: 回巢优势_终版_宽松.csv / 回巢优势_终版_高置信.csv')
