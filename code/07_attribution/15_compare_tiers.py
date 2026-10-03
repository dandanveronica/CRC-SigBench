#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
15_compare_tiers.py —— 四层对照的正式比较与检验

核心问题：已发表签名到底只是「挪用了运气」，还是真的比同等努力做出来的好？

四层（都在同一批外部队列上、用同一个定向 C-index 口径评估）：
    L0 队列内 split-half    同分布下的技术天花板（无批次效应）
    L1 自产签名             在队列 A 现场构建 → 丢到其余队列
    L2 已发表签名           原文训练集之外的独立队列
    L3 随机等大基因集       零假设基线

关键公平性处理：
    不能用「C 的绝对值」直接比，因为随机基线本身随【队列】和【基因数】变化。
    必须先把每一层都换算成 ΔC = C - 该(队列, 基因数)组合的随机基线中位数，
    再在同一刻度上比较。这一步错了，整个对照就是假的。

输入：audit_out_clean/per_cohort.csv       （L2 + L3 基线）
      audit_out_clean/de_novo_自产签名.csv （L1）
      audit_out_clean/de_novo_队列内上限.csv（L0，可选：14 跑完才有）
输出：audit_out_clean/四层对照.csv /.txt /.docx
"""
import os, sys
import numpy as np, pandas as pd
from scipy import stats

OUT = "audit_out_clean"
PER = os.path.join(OUT, "per_cohort.csv")
DNF = os.path.join(OUT, "de_novo_full.csv")      # 修正版：自带配对随机基线
DNO = os.path.join(OUT, "de_novo_自产签名.csv")  # 初版：只有 5 个训练队列，已弃用
SP  = os.path.join(OUT, "de_novo_队列内上限.csv")

lines = []
def emit(s=""):
    print(s); lines.append(s)

per = pd.read_csv(PER)
# 统一到已定稿的 320 条签名口径（签名库有 338 条，含后被剔除的 18 条）。
# 不统一会出现「同一篇文章里两处记录数不同」：此处 338×8=2704，正文其余部分 320×8=2560。
try:
    _keep = set(pd.read_csv(os.path.join(OUT, "签名320清单.csv"))["sig_id"])
    if per["sig_id"].nunique() > len(_keep):
        per = per[per["sig_id"].isin(_keep)].copy()
except FileNotFoundError:
    pass
per = per[per["status_r"] == "ok"].copy()
per["dC"] = per["C"] - per["null_median"]

# ---- 建立 (队列, 基因数) -> 随机基线中位数 的查表 ----
nullmed = per.groupby(["cohort", "n_measured"])["null_median"].median()

emit("=" * 78)
emit("四层对照：已发表签名 vs 同等努力的自产签名 vs 随机基因集")
emit("=" * 78)
emit(f"随机基线查表：{len(nullmed)} 个（队列 × 基因数）组合")

# ---- 优先用修正版：它为每个(验证队列, 基因数)单独重算了随机基线 ----
if os.path.exists(DNF):
    dn = pd.read_csv(DNF)
    # 统一用 17_null_matched.R 单独重算的配对基线（14b 内联算的那份曾全部为 NA）
    nm = os.path.join(OUT, "null_matched.csv")
    if os.path.exists(nm):
        nmt = pd.read_csv(nm)
        dn = dn.drop(columns=[c for c in ["null_median", "null_mean", "null_p95",
                                          "dC", "beat_random", "exceed_p95"]
                              if c in dn.columns])
        dn = dn.merge(nmt, left_on=["test", "n_gene"],
                      right_on=["cohort", "n_gene"], how="left")
    dn["dC"] = dn["C"] - dn["null_median"]
    dn["beat_random"] = (dn["C"] > dn["null_median"]).astype(int)
    dn_ok = dn.dropna(subset=["dC"]).copy()
    emit("来源：de_novo_full.csv + null_matched.csv"
         "（每个验证队列×基因数均单独重算随机基线，6 个队列全部参与训练）")
else:
    dn = pd.read_csv(DNO)
    dn["null"] = [nullmed.get((c, int(n)), np.nan)
                  for c, n in zip(dn["test"], dn["n_gene"])]
    dn["dC"] = dn["C"] - dn["null"]
    miss = dn["dC"].isna().sum()
    if miss:
        emit(f"（{miss} 条自产签名记录的基因数超出已发表签名覆盖的随机基线范围，"
             f"已在分层比较中剔除）")
    dn_ok = dn.dropna(subset=["dC"]).copy()
    dn_ok = dn_ok.rename(columns={"test": "test_c"})
    emit("警告：使用的是初版 de_novo_自产签名.csv（缺 GSE39582 训练，存在偏倚）")
dn_ok = dn_ok.rename(columns={"train": "train_c"})
dn_ok["train"] = dn_ok["train_c"]
emit(f"自产签名可比记录：{len(dn_ok)}")
emit("")

emit("-" * 78)
emit("一、对齐随机基线后的 ΔC-index（这才是可比的量）")
emit("-" * 78)
rows = []
tiers = [("L2 已发表签名", per["dC"].dropna()), ("L1 自产签名", dn_ok["dC"].dropna())]
if os.path.exists(SP):
    sp = pd.read_csv(SP)
    sp["null"] = [nullmed.get((c, int(n)), np.nan)
                  for c, n in zip(sp["cohort"], sp["K_actual"])]
    sp["dC"] = sp["C_val"] - sp["null"]
    sp_ok = sp.dropna(subset=["dC"])
    tiers.insert(0, ("L0 队列内 split-half", sp_ok["dC"].dropna()))
for name, x in tiers:
    w = stats.wilcoxon(x, alternative="greater")[1]
    rows.append(dict(层级=name, n=len(x), ΔC中位=round(x.median(), 4),
                     四分位=f"{x.quantile(.25):+.4f} ~ {x.quantile(.75):+.4f}",
                     优于随机比例=f"{(x>0).mean():.1%}", 对0检验p=f"{w:.3g}"))
tab = pd.DataFrame(rows)
emit(tab.to_string(index=False))
emit("")

emit("-" * 78)
emit("二、头对头：已发表 vs 自产（同一批外部队列、同一 Randomness 基线）")
emit("-" * 78)
a, b = dn_ok["dC"].dropna(), per["dC"].dropna()
u, p = stats.mannwhitneyu(a, b, alternative="two-sided")
cliffs = 2 * u / (len(a) * len(b)) - 1
emit(f"  自产 ΔC 中位 {a.median():+.4f}   已发表 ΔC 中位 {b.median():+.4f}")
emit(f"  差值（自产 - 已发表）= {a.median()-b.median():+.4f}")
emit(f"  Mann-Whitney p = {p:.3g}   Cliff's δ = {cliffs:+.3f}")
emit("")

emit("-" * 78)
emit("三、按基因数分层（控制工作量差異后的比较）")
emit("-" * 78)
bins = [2, 4, 6, 9, 12, 18, 25, 40, 100]
per["kb"] = pd.cut(per["n_measured"], bins, right=False)
dn_ok["kb"] = pd.cut(dn_ok["n_gene"], bins, right=False)
strat = []
for kb in per["kb"].cat.categories:
    pa = per.loc[per["kb"] == kb, "dC"].dropna()
    db = dn_ok.loc[dn_ok["kb"] == kb, "dC"].dropna()
    if len(pa) < 10 or len(db) < 3:
        continue
    pv = stats.mannwhitneyu(pa, db, alternative="two-sided")[1]
    strat.append(dict(基因数=str(kb), 已发表n=len(pa), 已发表ΔC=round(pa.median(), 4),
                      自产n=len(db), 自产ΔC=round(db.median(), 4),
                      差=round(db.median() - pa.median(), 4), p=f"{pv:.3g}"))
st = pd.DataFrame(strat)
emit(st.to_string(index=False))
emit("")

emit("-" * 78)
emit("四、逐自产签名：最好的那几个做到了多少")
emit("-" * 78)
kc = "K_target" if "K_target" in dn_ok.columns else "K"
top = dn_ok.sort_values("dC", ascending=False).head(12)[
    ["train", "test", kc, "n_gene", "C", "null_median", "dC"]].round(4)
emit(top.to_string(index=False))
emit("")

# 训练/验证配对的稳定性：同一train→不同test 的一致性
emit("-" * 78)
emit("五、每个训练队列自产签名在其他队列上的平均表现")
emit("-" * 78)
bt = dn_ok.groupby("train").agg(记录数=("C", "size"), ΔC中位=("dC", "median"),
                                C中位=("C", "median")).round(4)
emit(bt.to_string())
emit("")
bte = dn_ok.groupby("test").agg(记录数=("C", "size"), ΔC中位=("dC", "median")).round(4)
emit("作为验证队列时被打击的情况：")
emit(bte.to_string())

tab_out = tab.copy()
tab_out.to_csv(os.path.join(OUT, "四层对照.csv"), index=False, encoding="utf-8-sig")
st.to_csv(os.path.join(OUT, "四层对照_分层.csv"), index=False, encoding="utf-8-sig")
dn.to_csv(os.path.join(OUT, "自产签名_含随机基线.csv"), index=False, encoding="utf-8-sig")
with open(os.path.join(OUT, "四层对照.txt"), "w", encoding="utf-8") as f:
    f.write("\n".join(lines))

try:
    from docx import Document
    from docx.shared import Pt
    doc = Document()
    doc.styles["Normal"].font.name = "Times New Roman"
    doc.styles["Normal"].font.size = Pt(10.5)
    doc.add_heading("四层对照：同等努力的自产签名能做得多好", 0)
    doc.add_paragraph(
        f"对照 {per['sig_id'].nunique()} 个已发表签名 × {per['cohort'].nunique()} 个队列，"
        f"与 {dn_ok['train'].nunique()} 个训练队列现场生成的 {len(dn_ok)} 个自产签名。"
        f"所有层级均换算为 ΔC = C − 该(队列, 基因数)组合的随机基线中位数。").runs[0].bold = True
    doc.add_heading("一、四层对照结果", 1)
    t = doc.add_table(rows=1, cols=len(tab.columns)); t.style = "Light Grid Accent 1"
    for i, c in enumerate(tab.columns):
        cell = t.rows[0].cells[i]; cell.text = ""
        r = cell.paragraphs[0].add_run(str(c)); r.bold = True; r.font.size = Pt(8.5)
    for _, row in tab.iterrows():
        cs = t.add_row().cells
        for i, c in enumerate(tab.columns):
            cs[i].text = ""; cs[i].paragraphs[0].add_run(str(row[c])).font.size = Pt(8.5)
    doc.add_heading("二、头对头检验", 1)
    doc.add_paragraph(f"自产 ΔC 中位 {a.median():+.4f} vs 已发表 {b.median():+.4f}；"
                      f"差 {a.median()-b.median():+.4f}，Mann-Whitney p = {p:.3g}，"
                      f"Cliff's δ = {cliffs:+.3f}")
    doc.add_heading("三、按基因数分层", 1)
    t2 = doc.add_table(rows=1, cols=len(st.columns)); t2.style = "Light Grid Accent 1"
    for i, c in enumerate(st.columns):
        cell = t2.rows[0].cells[i]; cell.text = ""
        r = cell.paragraphs[0].add_run(str(c)); r.bold = True; r.font.size = Pt(8.5)
    for _, row in st.iterrows():
        cs = t2.add_row().cells
        for i, c in enumerate(st.columns):
            cs[i].text = ""; cs[i].paragraphs[0].add_run(str(row[c])).font.size = Pt(8.5)
    doc.add_paragraph()
    doc.add_paragraph("所有层级都做了方向定向并采用统一的随机基线查表，"
                      "因此 L1 与 L2 的 ΔC 可直接相减。").runs[0].italic = True
    doc.save(os.path.join(OUT, "四层对照.docx"))
    emit("已生成 Word：四层对照.docx")
except Exception as e:
    emit("Word 生成失败：" + str(e))

emit("\n输出：四层对照.csv / 四层对照_分层.csv / 自产签名_含随机基线.csv / 四层对照.txt|docx")
