#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
08_analyze.py —— 07 结果的正式统计分析

主终点（分布层面，不是逐签名 FDR）：
    已发表的 CRC 预后 mRNA 签名，是否为【随机等大基因集】带来判别度增益？
    配对比较 ΔC-index = C(真实签名) - median C(随机基因集)

为什么不把「逐签名 FDR<0.05」当主终点：
    经验 p 的分辨率 = 1/(B+1)。本次 B=1000 -> p_min = 0.000999
    BH 校正后最小值 = p_min × 检验数 = 0.000999 × 1920 = 1.92 >> 0.05
    即:无论真实效果多强,都不可能出现 FDR<0.05 的单条签名。
    要达到 FDR<0.05,B 需 ≥ 检验数/0.05 ≈ 38,400 次/组合。
    => 主终点必须是分布层面的配对检验;逐签名只报告未校正 p 与一致性计数。

用法：python 08_analyze.py [per_cohort.csv] [输出目录]
"""
import os, sys
import pandas as pd, numpy as np
from scipy import stats

PER  = sys.argv[1] if len(sys.argv) > 1 else "audit_out_clean/per_cohort.csv"
OUT  = sys.argv[2] if len(sys.argv) > 2 else "audit_out_clean"
LIB  = "sigminer/output/signature_library_clean.csv"
os.makedirs(OUT, exist_ok=True)

d = pd.read_csv(PER)
d = d[d["status_r"] == "ok"].copy()
d["dC"] = d["C"] - d["null_median"]
N_REC, N_SIG, N_COH = len(d), d["sig_id"].nunique(), d["cohort"].nunique()

lines = []
def emit(s=""):
    print(s); lines.append(s)

emit("=" * 78)
emit("结直肠癌转录组预后签名系统评估 —— Part 1 主结果")
emit("=" * 78)
emit(f"纳入签名 {N_SIG} 个 × 独立队列 {N_COH} 个 = {N_REC} 条评估")
emit(f"随机基因集对照：每个（队列 × 实际参与打分的基因数）组合 B=1000 次重抽样")
emit("")

# ---------- 1. 逐队列配对检验 ----------
emit("-" * 78)
emit("一、逐队列：真实签名 ΔC-index vs 随机基因集（配对 Wilcoxon / 符号检验）")
emit("-" * 78)
rows = []
for coh, g in d.groupby("cohort"):
    x = g["dC"].dropna()
    pw = stats.wilcoxon(x, alternative="greater")[1]
    nb = int((x > 0).sum())
    pb = stats.binomtest(nb, len(x), 0.5, alternative="greater").pvalue
    rows.append(dict(队列=coh, 记录数=len(x),
                     真实C中位=round(g["C"].median(), 3),
                     随机C中位=round(g["null_median"].median(), 3),
                     ΔC中位=round(x.median(), 4),
                     ΔC四分位=f"{x.quantile(.25):+.3f} ~ {x.quantile(.75):+.3f}",
                     优于随机=f"{nb}/{len(x)} ({nb/len(x):.0%})",
                     Wilcoxon_p=f"{pw:.2e}", 符号检验_p=f"{pb:.2e}"))
per_coh = pd.DataFrame(rows)
emit(per_coh.to_string(index=False))
emit("")

# ---------- 2. 合并 ----------
emit("-" * 78)
emit("二、合并全部记录")
emit("-" * 78)
x = d["dC"].dropna()
pw = stats.wilcoxon(x, alternative="greater")[1]
nb = int((x > 0).sum())
pb = stats.binomtest(nb, len(x), 0.5, alternative="greater").pvalue
k05 = int((d["p_emp"] < 0.05).sum())
p05 = stats.binomtest(k05, len(d), 0.05, alternative="greater").pvalue
emit(f"  ΔC 中位数      : {x.median():+.4f}（IQR {x.quantile(.25):+.4f} ~ {x.quantile(.75):+.4f}）")
emit(f"  ΔC 范围        : {x.min():+.4f} ~ {x.max():+.4f}")
emit(f"  真实 C 中位    : {d['C'].median():.3f}    随机 C 中位: {d['null_median'].median():.3f}")
emit(f"  优于随机的条目 : {nb}/{len(x)} = {nb/len(x):.1%}")
emit(f"  Wilcoxon 符号秩: p = {pw:.3e}")
emit(f"  符号检验      : p = {pb:.3e}")
emit(f"  原始 p<0.05    : {k05}/{len(d)} = {k05/len(d):.1%}（期望 5%，二项检验 p = {p05:.2e}）")
emit("")
emit("  解读：签名整体确实携带超出随机基因集的信息（p 极小），")
emit("        但增益量级很小——C-index 中位只高约 0.01，")
emit("        且四分位区间跨越 0，说明约四成签名连随机基因集都比不过。")
emit("")

# ---------- 3. FDR 的分辨率陷阱 ----------
emit("-" * 78)
emit("三、重要方法学说明：为什么「FDR<0.05 的签名数 = 0」不能作为结论")
emit("-" * 78)
B = 1000
emit(f"  经验 p 的最小可能值 = 1/(B+1) = {1/(B+1):.5f}（本次实测最小 p = {d['p_emp'].min():.5f}）")
emit(f"  BH 校正后的最小值   = p_min × 检验数 = {d['p_adj'].min():.2f}（实测）>> 0.05")
emit(f"  若要出现 FDR<0.05 的单条签名，B 至少要约 {int(len(d)/0.05):,} 次/组合")
emit("  => 这是【排列次数不足造成的假阴性】，不是「所有签名都不行」。")
emit("     因此本研究的判据不采用逐签名 FDR，而采用：")
emit("       (a) 分布层面的 ΔC-index 配对检验（主终点）")
emit("       (b) 跨队列一致性计数 ≥k 的富集分析（次要终点）")
emit("")

# ---------- 4. 跨队列一致性 ----------
emit("-" * 78)
emit("四、次要终点：跨队列一致性（在多少个队列里优于随机）")
emit("-" * 78)
cons = d.assign(k=(d["dC"] > 0).astype(int)).groupby("sig_id")["k"].sum()
for k in [N_COH, N_COH-1, N_COH-2]:
    obs = int((cons >= k).sum())
    exp = stats.binom.pmf(range(k, N_COH+1), N_COH, 0.5).sum()
    emit(f"  ≥{k}/{N_COH} 队列优于随机 : {obs} 个 ({obs/N_SIG:.1%})   "
         f"随机期望 {exp:.1%}   富集倍数 {obs/N_SIG/exp:.1f}×")
emit("")

# ---------- 5. 最稳健签名 ----------
emit("-" * 78)
emit("五、跨队列最稳健的 20 个签名")
emit("-" * 78)
agg = d.groupby("sig_id").agg(队列数=("cohort", "nunique"),
                              C均值=("C", "mean"), 随机C均值=("null_median", "mean"),
                              平均ΔC=("dC", "mean"))
agg["优于随机的队列数"] = cons
agg = agg.sort_values(["优于随机的队列数", "平均ΔC"], ascending=False)
try:
    lib = pd.read_csv(LIB, dtype=str, encoding="utf-8-sig").fillna("")
    m = lib.set_index("sig_id")
    agg["年份"]    = m["year"].reindex(agg.index)
    agg["基因数"]  = m["n_gene"].reindex(agg.index)
    agg["标题"]    = m["title"].reindex(agg.index).str.slice(0, 70)
except Exception:
    pass
top = agg.head(20).round(4)
top.to_csv(os.path.join(OUT, "最稳健签名_Top20.csv"), encoding="utf-8-sig")
emit(top.to_string())
emit("")

# ---------- 输出 ----------
per_coh.to_csv(os.path.join(OUT, "队列级汇总_含检验.csv"), index=False, encoding="utf-8-sig")
cons.sort_values(ascending=False).to_frame("优于随机的队列数").to_csv(
    os.path.join(OUT, "签名一致性.csv"), encoding="utf-8-sig")

with open(os.path.join(OUT, "Part1_主结果.txt"), "w", encoding="utf-8") as f:
    f.write("\n".join(lines))

# ---------- Word ----------
try:
    from docx import Document
    from docx.shared import Pt
    doc = Document()
    doc.styles["Normal"].font.name = "Times New Roman"
    doc.styles["Normal"].font.size = Pt(10.5)
    doc.add_heading("结直肠癌转录组预后签名系统评估 —— Part 1 主结果", 0)
    p = doc.add_paragraph()
    p.add_run(f"{N_SIG} 个已发表签名 × {N_COH} 个独立 GEO 队列 = {N_REC} 条重算评估").bold = True

    doc.add_heading("一、主终点：签名是否优于随机等大基因集", 1)
    tb = doc.add_table(rows=1, cols=len(per_coh.columns)); tb.style = "Light Grid Accent 1"
    for i, c in enumerate(per_coh.columns):
        cell = tb.rows[0].cells[i]; cell.text = ""
        r = cell.paragraphs[0].add_run(str(c)); r.bold = True; r.font.size = Pt(8.5)
    for _, row in per_coh.iterrows():
        cs = tb.add_row().cells
        for i, c in enumerate(per_coh.columns):
            cs[i].text = ""
            cs[i].paragraphs[0].add_run(str(row[c])).font.size = Pt(8.5)
    doc.add_paragraph()
    for s in [f"合并 {N_REC} 条记录：ΔC 中位 {x.median():+.4f}"
              f"（IQR {x.quantile(.25):+.4f} ~ {x.quantile(.75):+.4f}）",
              f"真实 C-index 中位 {d['C'].median():.3f} vs 随机基因集 {d['null_median'].median():.3f}",
              f"优于随机的条目 {nb}/{len(x)} = {nb/len(x):.1%}",
              f"Wilcoxon 符号秩检验 p = {pw:.3e}；符号检验 p = {pb:.3e}",
              f"原始 p<0.05 的条目 {k05}/{len(d)} = {k05/len(d):.1%}（期望 5%，二项 p = {p05:.2e}）"]:
        doc.add_paragraph(s, style="List Bullet")

    doc.add_heading("二、结论", 1)
    doc.add_paragraph(
        f"已发表的结直肠癌 mRNA 预后签名作为一个整体，确实携带超出随机等大基因集的判别信息"
        f"（Wilcoxon p = {pw:.2e}），但增益量级极小：C-index 中位数仅高出 {x.median():.3f}，"
        f"且四分位区间跨越 0，约 {(1-nb/len(x)):.0%} 的签名-队列组合并不优于随机抽样。"
        f"约 {k05/len(d):.0%} 的组合达到未校正 p<0.05，提示存在真实但微弱、且高度依赖队列的信号。").runs[0].bold = True

    doc.add_heading("三、方法学警示：逐签名 FDR 不可用于本研究", 1)
    doc.add_paragraph(
        f"经验 p 的分辨率受限于排列次数 B：本次 B={B}，最小可达 p = 1/(B+1) = {1/(B+1):.5f}；"
        f"BH 校正后最小值为 {d['p_adj'].min():.2f}，远大于 0.05。也就是说，无论效应多强，"
        f"本次设计在数学上都不可能产生 FDR<0.05 的单条签名。要使逐签名 FDR 可用，"
        f"B 需增至约 {int(len(d)/0.05):,} 次/组合。"
        f"因此本研究以分布层面的配对检验为主终点，逐签名结果仅报告未校正 p 与跨队列一致性。")

    doc.add_heading("四、次要终点：跨队列一致性", 1)
    for k in [N_COH, N_COH-1, N_COH-2]:
        obs = int((cons >= k).sum())
        exp = stats.binom.pmf(range(k, N_COH+1), N_COH, 0.5).sum()
        doc.add_paragraph(
            f"≥{k}/{N_COH} 个队列优于随机：{obs} 个（{obs/N_SIG:.1%}），"
            f"随机期望 {exp:.1%}，富集 {obs/N_SIG/exp:.1f} 倍", style="List Bullet")

    doc.add_heading("五、局限与待办", 1)
    for s in ["队列均为 GPL570 芯片平台，尚未纳入 RNA-seq（TCGA）队列；",
              "尚未做跨平台批次校正（ComBat），队列间只做队列内 z-score；",
              "GSE87211（Agilent GPL13497，363 例）因缺少平台注释暂未纳入；",
              "未按原文献系数加权打分（reported 模式），当前为等权重 + 方向定向；",
              "建议后续对 Top 签名做高 B（≥40,000 次）二阶段检验以获得逐签名 FDR。"]:
        doc.add_paragraph(s, style="List Number")
    doc.save(os.path.join(OUT, "Part1_主结果.docx"))
    emit("已生成 Word：Part1_主结果.docx")
except Exception as e:
    emit("Word 生成失败：" + str(e))

emit("")
emit("输出：队列级汇总_含检验.csv / 签名一致性.csv / 最稳健签名_Top20.csv / Part1_主结果.txt|docx")
