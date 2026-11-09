# -*- coding: utf-8 -*-
"""S19-A 基因复现频谱 —— 312 签名基因频率（新分析 ③，2026-11-09）

输入：data/signatures/signature_library_312_genes.csv（genes 列，';' 分隔，HGNC symbols）
输出：audit_out_clean/pathway/
  gene_recurrence.csv        每基因出现于多少签名（全量）
  recurrence_summary.txt     频谱汇总（top 基因、1/2/3+ 签名基因占比等）

对标：TNBC 62 签名 meta 分析（bioRxiv 2026-07）的 gene-level divergence 呈现
口径：pure_mRNA 与全部签名均按库原样（与主管线打分一致，genes 列逐字拆分）
"""
import os
import csv
from collections import Counter

ROOT = ('/Users/yangdan/Desktop/阳丹/文章课题/文章/CRC-SigBench')
OUTD = os.path.join(ROOT, 'audit_out_clean', 'pathway')
os.makedirs(OUTD, exist_ok=True)

rows = list(csv.DictReader(open(f'{ROOT}/data/signatures/signature_library_312_genes.csv',
                                encoding='utf-8-sig')))
assert len(rows) == 312, len(rows)

gene2sigs = Counter()
n_per_sig = []
for r in rows:
    genes = [g.strip() for g in r['genes'].replace(',', ';').split(';') if g.strip()]
    n_per_sig.append(len(genes))
    gene2sigs.update(set(genes))          # set：同一签名内重复基因只计一次

uniq_genes = len(gene2sigs)
total_mentions = sum(gene2sigs.values())
singles = sum(1 for v in gene2sigs.values() if v == 1)
multi2 = sum(1 for v in gene2sigs.values() if v >= 2)
multi3 = sum(1 for v in gene2sigs.values() if v >= 3)

top = gene2sigs.most_common(20)
max_overlap = top[0][1] if top else 0

# 覆盖率：k 个最常见基因累计覆盖多少签名（≥1 个 top-k 基因）
sig_gene_sets = []
for r in rows:
    genes = set(g.strip() for g in r['genes'].replace(',', ';').split(';') if g.strip())
    sig_gene_sets.append(genes)
cov = []
for k in (1, 5, 10, 20, 50):
    topk = [g for g, _ in top[:k]]
    n_cov = sum(1 for s in sig_gene_sets if s & set(topk))
    cov.append((k, n_cov, 100 * n_cov / 312))

with open(f'{OUTD}/gene_recurrence.csv', 'w', newline='', encoding='utf-8') as f:
    w = csv.writer(f)
    w.writerow(['gene', 'n_signatures', 'pct_of_312'])
    for g, v in sorted(gene2sigs.items(), key=lambda x: (-x[1], x[0])):
        w.writerow([g, v, f'{100*v/312:.2f}'])

lines = [
    '════ 312 签名基因复现频谱 ════',
    f'签名数 312 | 基因提及总次数 {total_mentions} | 唯一基因 {uniq_genes}',
    f'每签名基因数: median {sorted(n_per_sig)[156]}, range {min(n_per_sig)}-{max(n_per_sig)}',
    f'只出现 1 次的基因: {singles} ({100*singles/uniq_genes:.1f}%)',
    f'出现 ≥2 次的基因: {multi2} ({100*multi2/uniq_genes:.1f}%)',
    f'出现 ≥3 次的基因: {multi3} ({100*multi3/uniq_genes:.1f}%)',
    f'最常复现基因: {top[0][0]} ({max_overlap}/312)',
    '',
    'Top 20 基因（出现签名数）:',
]
for g, v in top:
    lines.append(f'  {g:12s} {v:3d}/312 ({100*v/312:.1f}%)')
lines.append('')
lines.append('Top-k 基因覆盖率（≥1 个 top-k 基因的签名占比）:')
for k, n, pct in cov:
    lines.append(f'  top-{k:2d}: {n}/312 ({pct:.1f}%)')
open(f'{OUTD}/recurrence_summary.txt', 'w', encoding='utf-8').write('\n'.join(lines))
print('\n'.join(lines))
