# -*- coding: utf-8 -*-
"""
32_tcga_stratified.py —— 第 8 队列 TCGA 到位后的分层分析

要解决的核心问题：TCGA 上 ΔC_随机 = +0.0205（8 个队列里最高），
但 ΔC_临床 = +0.0002（几乎归零）。这两个数字可以被两种机制解释：
  A. 回巢效应：多数签名本就是用 TCGA 训练的，在"自家"队列上天然占优
  B. 临床协变量更强：TCGA 的分期/年龄/性别完整度高，临床模型本身就好
两者政策含义完全不同，必须拆开 —— 办法是按【是否用 TCGA 训练】分层。

同时产出：
  ① 8 队列合并的整体指标更新
  ② RNA-seq vs 芯片的平台效应
"""
import os, re, glob
import numpy as np, pandas as pd
from scipy import stats

BASE = "audit_out_clean"
OUR = ['GSE39582', 'GSE17536', 'GSE17537', 'GSE29621',
       'GSE38832', 'GSE72970', 'GSE87211']

allq = pd.read_csv("audit_out_8/per_cohort.csv")   # 8 队列、delta_C 已修复版本
pc7 = allq[allq.cohort != "TCGA_COADREAD"].copy()
tcg = allq[allq.cohort == "TCGA_COADREAD"].copy()

for d in (pc7, tcg):
    d['d_rand'] = d.C - d.null_median
    d['d_clin'] = d.delta_C

# ---- 签名口径与已定稿的 320 条保持一致 ----
# 统一到已定稿的 320 条口径（签名库 338 条含后被剔除的 18 条）
keep = set(pd.read_csv(f"{BASE}/签名320清单.csv").sig_id)
pc7 = pc7[pc7.sig_id.isin(keep)].copy()
tcg = tcg[tcg.sig_id.isin(keep)].copy()
print(f"定稿签名 {len(keep)} 条，TCGA {len(tcg)} 条，每队列 {len(tcg)} 条\n")

pc8 = pd.concat([pc7, tcg], ignore_index=True)
pc8.to_csv(f"{BASE}/per_cohort_8队列.csv", index=False, encoding="utf-8-sig")


def blk(t): print("\n" + "=" * 72 + f"\n{t}\n" + "=" * 72)


# ============================ ① 8 队列整体 ============================
blk("① 8 队列整体（%d 签名 x 8 队列 = %d 条）" % (pc8.sig_id.nunique(), len(pc8)))
print(f"{'队列':<16}{'n':>5}{'C中位':>8}{'ΔC随机中位':>11}{'%>0':>7}"
      f"{'ΔC临床中位':>11}{'%>0':>7}{'校准':>7}{'事件数':>7}")
rows = []
for c, s in pc8.groupby('cohort'):
    rows.append(dict(cohort=c, n=len(s), C=s.C.median(),
                     d_rand=s.d_rand.median(), pos_rand=100 * (s.d_rand > 0).mean(),
                     d_clin=s.d_clin.median(), pos_clin=100 * (s.d_clin > 0).mean(),
                     calib=s.calib_slope.median()))
T = pd.DataFrame(rows).sort_values('d_rand')
for _, r in T.iterrows():
    print(f"{r.cohort:<16}{int(r.n):>5}{r.C:>8.4f}{r.d_rand:>+11.4f}{r.pos_rand:>6.1f}%"
          f"{r.d_clin:>+11.4f}{r.pos_clin:>6.1f}%{r.calib:>7.3f}{'':>7}")
print(f"\n全部 8 队列: ΔC随机 中位 {pc8.d_rand.median():+.4f}，"
      f"{100*(pc8.d_rand>0).mean():.1f}% 为正")
print(f"            ΔC临床 中位 {pc8.d_clin.median():+.4f}，"
      f"{100*(pc8.d_clin>0).mean():.1f}% 为正")
print(f"（7 队列版对照）ΔC随机 {pc7.d_rand.median():+.4f}，ΔC临床 {pc7.d_clin.median():+.4f}")

only_t = pc8[pc8.cohort == 'TCGA_COADREAD']
print(f"\nTCGA 在 8 队列中的 ΔC随机 排名: "
      f"{list(T.sort_values('d_rand').cohort).index('TCGA_COADREAD')+1} / 8（从低到高）")
print(f"Wilcoxon ΔC随机 TCGA vs 其余: "
      f"p={stats.mannwhitneyu(only_t.d_rand.dropna(), pc8[pc8.cohort!='TCGA_COADREAD'].d_rand.dropna(), alternative='greater').pvalue:.3g}")


# ============================ ② 训练来源判定 ============================
blk("② 判定每个签名是否【用 TCGA 训练】")
GEO = re.compile(r'GSE\d{4,6}')
TW = re.compile(r'\b(train|training|discovery|derivation|deriving|construction|'
                r'constructed|developed)\b', re.I)
VW = re.compile(r'\b(valid|validation|testing|test set|external|verified)\b', re.I)

meta = pd.read_csv(f"{BASE}/全文方法学提取.csv")[['sig_id', 'pmcid']]


def norm(s): return re.sub(r'[^A-Z0-9]', '', str(s)).upper()


def near(pat, txt, pos, win=150):
    lo, hi = max(0, pos - win), min(len(txt), pos + win)
    seg = txt[lo:hi]; rel = pos - lo
    p = [m.start() for m in pat.finditer(seg)]
    return min([abs(x - rel) for x in p], default=10 ** 9)


rec = []
for fn in sorted(glob.glob('fulltext_xml/*')):
    pmc = re.sub(r'\.(xml|txt|html)$', '', os.path.basename(fn))
    t = re.sub(r'<[^>]+>', ' ', open(fn, encoding='utf-8', errors='ignore').read())
    tc_train = None
    # TCGA 出现处（含全称）：看最近的是 train 词还是 valid 词
    for m in re.finditer(r'TCGA|The Cancer Genome Atlas', t, re.I):
        dt, dv = near(TW, t, m.start()), near(VW, t, m.start())
        if dt < dv:
            tc_train = 1
            break
        if dv < dt:
            tc_train = 0
    # GEO 训练队列（沿用 25d 的字符距离法）
    homes = set()
    for m in GEO.finditer(t):
        g = m.group()
        if g not in OUR:
            continue
        dt, dv = near(TW, t, m.start()), near(VW, t, m.start())
        if dt < dv:
            homes.add(g)
    rec.append(dict(pmcid_raw=pmc, tcga_train=tc_train,
                    home=';'.join(sorted(homes))))
J = pd.DataFrame(rec)
J['k'] = J.pmcid_raw.map(norm)
meta['k'] = meta.pmcid.map(norm)
M = meta.merge(J.drop(columns=['pmcid_raw']), on='k', how='left')
M['home'] = M.home.fillna('')
M.to_csv(f"{BASE}/训练来源判定_含TCGA.csv", index=False, encoding='utf-8-sig')

tc = M.set_index('sig_id').tcga_train
print(f"能判定的签名 {int(tc.notna().sum())} / {len(M)}")
print(f"  判定为【用 TCGA 训练】: {int((tc==1).sum())}")
print(f"  判定为【未用 TCGA 训练】: {int((tc==0).sum())}")
print(f"  无法判定（全文没提 TCGA 或其角色不明）: {int(tc.isna().sum())}")

only_t = only_t.copy()
only_t['tcga_train'] = only_t.sig_id.map(tc)
only_t['n_home'] = only_t.sig_id.map(
    M.set_index('sig_id').home.fillna('').apply(lambda s: len([x for x in s.split(';') if x])))


# ============================ ③ 回巢效应 ============================
blk("③ 关键：TCGA 上的表现，按【是否用 TCGA 训练】分层")
sub = only_t.dropna(subset=['tcga_train', 'd_rand'])
print(f"可分层签名 {len(sub)} 条\n")
print(f"{'':<22}{'n':>5}{'C中位':>8}{'ΔC随机中位':>11}{'%>0':>7}"
      f"{'ΔC临床中位':>11}{'%>0':>7}")
for k, lab in [(1, '用 TCGA 训练（"回家考"）'), (0, '未用 TCGA 训练（陌生队列）')]:
    s = sub[sub.tcga_train == k]
    if len(s) < 5:
        continue
    print(f"{lab:<22}{len(s):>5}{s.C.median():>8.4f}{s.d_rand.median():>+11.4f}"
          f"{100*(s.d_rand>0).mean():>6.1f}%{s.d_clin.median():>+11.4f}"
          f"{100*(s.d_clin>0).mean():>6.1f}%")
a = sub[sub.tcga_train == 1].d_rand.dropna()
b = sub[sub.tcga_train == 0].d_rand.dropna()
if len(a) > 5 and len(b) > 5:
    u = stats.mannwhitneyu(a, b, alternative='greater')
    print(f"\n  ΔC随机 差值中位 {a.median()-b.median():+.4f}   Mann-Whitney p={u.pvalue:.3g}")
    ca = sub[sub.tcga_train == 1].d_clin.dropna()
    cb = sub[sub.tcga_train == 0].d_clin.dropna()
    u2 = stats.mannwhitneyu(ca, cb, alternative='greater')
    print(f"  ΔC临床 差值中位 {ca.median()-cb.median():+.4f}   Mann-Whitney p={u2.pvalue:.3g}")

# 三重对照：TCGA 训练组在【GEO 队列】上是否也占优（若不占优，才是真的回巢效应）
geo = pc8[pc8.cohort != 'TCGA_COADREAD'].copy()
geo['tcga_train'] = geo.sig_id.map(tc)
g1 = geo[geo.tcga_train == 1].d_rand.dropna()
g0 = geo[geo.tcga_train == 0].d_rand.dropna()
print(f"\n  【对照】同样按 TCGA 训练分层，在 7 个 GEO 队列上:")
print(f"    用 TCGA 训练: ΔC随机中位 {g1.median():+.4f}（n={len(g1)}）")
print(f"    未用 TCGA 训练: ΔC随机中位 {g0.median():+.4f}（n={len(g0)}）")
if len(g1) > 5 and len(g0) > 5:
    print(f"    差 {g1.median()-g0.median():+.4f}  p={stats.mannwhitneyu(g1, g0, alternative='greater').pvalue:.3g}")
    print("    → 若此处≈0 而 TCGA 上明显>0，说明优势不是签名本身更好，而是回家考效应")


# ============================ ④ 平台效应 ============================
blk("④ 平台效应：RNA-seq（TCGA）vs 芯片（7 个 GEO）")
rna = pc8[pc8.cohort == 'TCGA_COADREAD']
geo_all = pc8[pc8.cohort != 'TCGA_COADREAD']
S = pc8.groupby('sig_id').agg(
    n=('C', 'size'), C_rna=('C', 'mean'),
    d_rand_all=('d_rand', 'mean'), d_clin_all=('d_clin', 'mean')).reset_index()
r = pc8[pc8.cohort == 'TCGA_COADREAD'].set_index('sig_id')
g = pc8[pc8.cohort != 'TCGA_COADREAD'].groupby('sig_id').agg(
    C_geo=('C', 'mean'), d_rand_geo=('d_rand', 'mean'), d_clin_geo=('d_clin', 'mean'))
D = g.join(r[['C', 'd_rand', 'd_clin']].rename(
    columns={'C': 'C_rna', 'd_rand': 'd_rand_rna', 'd_clin': 'd_clin_rna'})).dropna()
print(f"跨两平台均可评估的签名 {len(D)} 条\n")
for a, b_, lab in [('C_rna', 'C_geo', 'C-index'),
                   ('d_rand_rna', 'd_rand_geo', 'ΔC随机'),
                   ('d_clin_rna', 'd_clin_geo', 'ΔC临床')]:
    rho, pv = stats.spearmanr(D[a], D[b_])
    w = stats.wilcoxon(D[a], D[b_])
    print(f"  {lab:<8} RNA-seq 中位 {D[a].median():+.4f}  芯片中位 {D[b_].median():+.4f}  "
          f"配对 Wilcoxon p={w.pvalue:.3g}")
    print(f"  {'':<8} 跨平台 Spearman ρ={rho:+.3f}  p={pv:.3g}"
          f"   <-ρ 低说明签名在哪平台好用基本不可预测")

D.to_csv(f"{BASE}/跨平台配对_8队列.csv", encoding="utf-8-sig")
sub.to_csv(f"{BASE}/TCGA分层.csv", index=False, encoding="utf-8-sig")
T.to_csv(f"{BASE}/队列汇总_8队列.csv", index=False, encoding="utf-8-sig")
print("\n已存: per_cohort_8队列.csv / 队列汇总_8队列.csv / TCGA分层.csv / "
      "跨平台配对_8队列.csv / 训练来源判定_含TCGA.csv")
