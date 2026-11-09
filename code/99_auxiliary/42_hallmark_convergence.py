# -*- coding: utf-8 -*-
"""42_hallmark_convergence.py —— 312 签名 Hallmark 通路收敛（新分析 ①，2026-11-09）

设计：
  基因集  data/proliferation/hallmark_2020_symbols.txt（MSigDB Hallmark 2020，50 条，
          Enrichr API 2026-11-09 抓取，留痕见 README_hallmark.txt）
  universe audit_out_clean/pathway/measured_universe.csv（41 脚本 dump 的 8 队列
          并集实测基因；ORA 在 universe 内做，K/n 均先与 universe 求交）
  检验    每签名 × 每 Hallmark 通路超几何 ORA（右尾），签名内 BH FDR<0.05 记"富集"
  收敛指标 每通路被多少签名富集（n_signatures_enriched / 312）
  null 对照 size-matched 随机基因集（每签名 1 条、B=100 轮，从 universe 抽同样大小），
          得通路的 null 富集签名数分布 → 经验 p（收敛是否超出随机）＋ fold

输出：audit_out_clean/pathway/
  hallmark_ora_per_signature.csv   逐签名×逐通路（仅 FDR<0.05 行）
  hallmark_convergence.csv         每通路汇总（富集签名数、null 中位、fold、经验 p）
  convergence_summary.txt          汇总文本（top 通路、增殖类通路专项）
"""
import os
import csv
import math
import random
import numpy as np

random.seed(20261109)
np.random.seed(20261109)

ROOT = '/Users/yangdan/Desktop/阳丹/文章课题/文章/CRC-SigBench'
OUTD = f'{ROOT}/audit_out_clean/pathway'
os.makedirs(OUTD, exist_ok=True)

# ---------------- 输入 ----------------
rows = list(csv.DictReader(open(f'{ROOT}/data/signatures/signature_library_312_genes.csv',
                                encoding='utf-8-sig')))
assert len(rows) == 312
sig_sets = {}
for r in rows:
    genes = set(g.strip() for g in r['genes'].replace(',', ';').split(';') if g.strip())
    sig_sets[r['sig_id']] = genes

universe_rows = list(csv.DictReader(open(f'{OUTD}/measured_universe.csv', encoding='utf-8')))
universe = set(r['gene'] for r in universe_rows if int(r['n_cohorts_measured']) >= 1)
N = len(universe)

hallmark = {}
for line in open(f'{ROOT}/data/proliferation/hallmark_2020_symbols.txt', encoding='utf-8'):
    parts = [p.strip() for p in line.rstrip('\n').split('\t')]
    name, genes = parts[0], [g for g in parts[2:] if g]
    hallmark[name] = set(genes) & universe          # K = 通路∩universe
assert len(hallmark) == 50, len(hallmark)

# Enrichr 版 Hallmark 2020 无 "G2M Checkpoint"（实测 50 条命名），增殖专项用实际存在的 4 条
PROLIF_PATHS = ['E2F Targets', 'Myc Targets V1', 'Myc Targets V2', 'Mitotic Spindle']


def hypergeom_sf(k, N, K, n):
    """P(X >= k)，X~Hypergeom(N, K, n)。pmf 逐项 lgamma 累加（n 小，数值安全）"""
    from math import lgamma, exp
    if k > min(K, n):
        return 1.0
    def log_pmf(i):
        return (lgamma(K+1)-lgamma(i+1)-lgamma(K-i+1)
                + lgamma(N-K+1)-lgamma(n-i+1)-lgamma(N-K-n+i+1)
                - (lgamma(N+1)-lgamma(n+1)-lgamma(N-n+1)))
    k_eff = max(k, 1)  # k=0 时 sf=1
    return min(1.0, sum(exp(log_pmf(i)) for i in range(k_eff, min(K, n)+1)))


def bh_fdr(pvals):
    """BH 校正，返回同序 fdr 数组"""
    n = len(pvals)
    order = np.argsort(pvals)
    ranked = np.asarray(pvals)[order]
    q = ranked * n / (np.arange(n) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.clip(q, 0, 1)
    return out


def ora_one(genes):
    """单签名对 50 通路 ORA，返回 {pathway: (k, K, n, p, fdr)}，fdr 为签名内 BH"""
    in_uni = genes & universe
    n = len(in_uni)
    res = {}
    ps = []
    keys = []
    for name, H in hallmark.items():
        K = len(H)
        k = len(in_uni & H)
        p = hypergeom_sf(k, N, K, n) if k >= 1 else 1.0
        res[name] = (k, K, n, p)
        ps.append(p)
        keys.append(name)
    fdrs = bh_fdr(ps)
    return {keys[i]: (*res[keys[i]], fdrs[i]) for i in range(len(keys))}, in_uni


# ---------------- 主计算：312 签名 ----------------
sig_order = [r['sig_id'] for r in rows]
per_sig_hits = {}          # sig -> set(enriched pathways)
ora_rows = []
for sid in sig_order:
    ora, _ = ora_one(sig_sets[sid])
    hits = set()
    for name, (k, K, n, p, fdr) in ora.items():
        if fdr < 0.05:
            hits.add(name)
            ora_rows.append([sid, name, k, K, n, f'{p:.3e}', f'{fdr:.3e}'])
    per_sig_hits[sid] = hits

# ---------------- null：size-matched 随机基因集，B=100 轮 ----------------
B = 100
sizes = {sid: len(sig_sets[sid] & universe) for sid in sig_order}
null_counts = {name: np.zeros(B, dtype=int) for name in hallmark}
for b in range(B):
    rnd_hits = {name: 0 for name in hallmark}
    for sid in sig_order:
        rnd = set(np.random.choice(sorted(universe), sizes[sid], replace=False))
        ora, _ = ora_one(rnd)
        for name, (k, K, n, p, fdr) in ora.items():
            if fdr < 0.05:
                rnd_hits[name] += 1
    for name in hallmark:
        null_counts[name][b] = rnd_hits[name]

# ---------------- 汇总 ----------------
sig_n_enriched = [len(per_sig_hits[s]) for s in sig_order]
conv_rows = []
for name, H in hallmark.items():
    obs = sum(1 for s in sig_order if name in per_sig_hits[s])
    nd = null_counts[name]
    emp_p = (np.sum(nd >= obs) + 1) / (B + 1)
    fold = obs / max(0.5, np.median(nd))
    conv_rows.append([name, obs, f'{100*obs/312:.1f}',
                      f'{np.median(nd):.1f}', f'{fold:.1f}', f'{emp_p:.4f}'])
conv_rows.sort(key=lambda r: -int(r[1]))

with open(f'{OUTD}/hallmark_ora_per_signature.csv', 'w', newline='', encoding='utf-8') as f:
    w = csv.writer(f)
    w.writerow(['sig_id', 'pathway', 'k_overlap', 'K_pathway_in_universe',
                'n_sig_in_universe', 'p_hyper', 'fdr_bh_within_signature'])
    w.writerows(ora_rows)

with open(f'{OUTD}/hallmark_convergence.csv', 'w', newline='', encoding='utf-8') as f:
    w = csv.writer(f)
    w.writerow(['pathway', 'n_signatures_enriched', 'pct_of_312',
                'null_median_B100', 'fold_over_null', 'empirical_p'])
    w.writerows(conv_rows)

sig_hits_series = sorted(sig_n_enriched)
none_hits = sum(1 for x in sig_n_enriched if x == 0)
prolif_stats = {}
for name in PROLIF_PATHS:
    obs = sum(1 for s in sig_order if name in per_sig_hits[s])
    prolif_stats[name] = obs

lines = [
    '════ 312 签名 Hallmark 通路收敛 ════',
    f'universe = {N} 基因（8 队列并集实测）；Hallmark {len(hallmark)} 通路；'
    f'signature 内 BH FDR<0.05；null = size-matched 随机集 B=100',
    f'每签名富集通路数: median {sig_hits_series[len(sig_hits_series)//2]}, '
    f'range {sig_hits_series[0]}-{sig_hits_series[-1]}；无任何富集的签名 {none_hits}/312',
    '',
    'Top 15 收敛通路（按富集签名数）:',
]
for r in conv_rows[:15]:
    lines.append(f"  {r[0]:35s} {int(r[1]):>3d}/312 ({r[2]:>5s}%)  null {r[3]:>5s}  "
                 f"fold {r[4]:>5s}  p {r[5]}")
lines.append('')
lines.append('增殖类通路专项:')
for name in PROLIF_PATHS:
    obs = prolif_stats[name]
    nd = null_counts[name]
    emp_p = (np.sum(nd >= obs) + 1) / (B + 1)
    lines.append(f'  {name:20s} {obs}/312 ({100*obs/312:.1f}%)  '
                 f'null median {np.median(nd):.0f}  emp p {emp_p:.4f}')
open(f'{OUTD}/convergence_summary.txt', 'w', encoding='utf-8').write('\n'.join(lines))
print('\n'.join(lines))
