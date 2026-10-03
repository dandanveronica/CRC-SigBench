#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
11_part2_attribution.py —— Part 2：什么样的签名更能打？

结局（签名层面，n=320）：
    Y1  跨 6 队列的平均 ΔC-index（连续，主结局）
    Y2  优于随机的队列数 0~6（计数，次要结局）
    Y3  最大标准化位置 max z(C vs 零分布)（次要结局）

解释变量（+ 为 09/10 抽取的元数据）：
    文献层面：发表年份、期刊
    设计层面：基因数(log)、基因复用度、独有基因占比、GPL570 检出率
    方法层面：LASSO / 单因素Cox / 多因素Cox / 逐步 / 机器学习 / WGCNA /
              列线图 / 单细胞  +  训练集样本量(log) / 是否外部验证 /
              是否多因素校正 / 是否公开系数
    主题层面：生物学主主题（免疫、代谢、铁死亡、CAF…）

统计：
    单因素  连续 -> Spearman；二分类 -> Mann-Whitney；多分类 -> Kruskal-Wallis
    多重性  BH 校正（全部检验作为一个 family，报告 p 与 q）
    多因素  OLS + HC1 稳健标准误；另做 L1(LASSO) 无关——只用稳健 OLS

输入：audit_out_clean/per_cohort.csv
      audit_out_clean/签名元数据_归因变量.csv
      audit_out_clean/全文方法学提取.csv（可选）
输出：audit_out_clean/Part2_归因分析.txt / .docx
      audit_out_clean/Part2_单因素.csv / Part2_多因素.csv
"""
import os, re, sys
import numpy as np, pandas as pd
from scipy import stats

PER   = "audit_out_clean/per_cohort.csv"
META  = "audit_out_clean/签名元数据_归因变量.csv"
FT    = "audit_out_clean/全文方法学提取.csv"
OUT   = "audit_out_clean"
os.makedirs(OUT, exist_ok=True)

lines = []
def emit(s=""):
    print(s); lines.append(s)

# ------------------------------------------------------------------ 结局
d = pd.read_csv(PER)
d = d[d["status_r"] == "ok"].copy()
d["dC"] = d["C"] - d["null_median"]
Y = d.groupby("sig_id").agg(
    mean_dC   = ("dC", "mean"),
    med_dC    = ("dC", "median"),
    beats     = ("dC", lambda s: int((s > 0).sum())),
    C_mean    = ("C", "mean"),
).reset_index()

# ------------------------------------------------------------- 解释变量
m = pd.read_csv(META)
m["year"] = pd.to_numeric(m["year"], errors="coerce")
ft = None
if os.path.exists(FT):
    ft = pd.read_csv(FT)
    for c in ["是否外部验证", "是否多因素校正", "是否给出风险公式"]:
        if c in ft.columns:
            ft[c] = pd.to_numeric(ft[c], errors="coerce")
    ft["训练集样本量"] = pd.to_numeric(ft.get("训练集样本量"), errors="coerce")
    m = m.merge(ft[["sig_id", "训练集样本量", "是否外部验证",
                    "是否多因素校正", "是否给出风险公式"]].drop_duplicates("sig_id"),
                on="sig_id", how="left")

df = m.merge(Y, on="sig_id", how="inner")
df["log_n_gene"]  = np.log(df["n_gene"])
df["log_n_train"] = np.log(pd.to_numeric(df.get("训练集样本量"), errors="coerce"))
emit("=" * 78)
emit("Part 2：什么样的签名更能打？")
emit("=" * 78)
emit(f"可分析的签名：{len(df)} 个（全部完成了 6 个队列的重算）")
if ft is not None and "训练集样本量" in df.columns:
    nn = df["log_n_train"].notna().sum()
    emit(f"其中 {nn} 个可从全文中提取到训练集样本量")

CONT = [("year", "发表年份"), ("log_n_gene", "基因数(log)"),
        ("gene_popularity_mean", "基因复用度"),
        ("pct_unique_gene", "独有基因占比"),
        ("pct_hot_gene", "热点基因占比"),
        ("检出率均值", "GPL570检出率")]
if "log_n_train" in df.columns and df["log_n_train"].notna().sum() >= 50:
    CONT.append(("log_n_train", "训练集样本量(log)"))

BIN = [c for c in df.columns if c.startswith("方法_")] + \
      [c for c in ["是否外部验证", "是否多因素校正", "是否给出风险公式"] if c in df.columns]
THEME_COL = None
if "主主题" in df.columns:
    THEME_COL = "主主题"

MAIN_Y = "mean_dC"
emit(f"\n主结局 Y = {MAIN_Y}（跨 6 队列平均 ΔC-index）")
emit("")

# ------------------------------------------------------------- 单因素检验
defuni = []
for col, lab in [(c, c) for c, _ in CONT]:
    sub = df[[col, MAIN_Y]].dropna()
    if len(sub) < 30:
        continue
    rho, p = stats.spearmanr(sub[col], sub[MAIN_Y])
    defuni.append(dict(变量=CONT[[c for c, _ in CONT].index(col)][1], 类型="连续",
                        n=len(sub), 效应量=f"ρ={rho:+.3f}", p=p))

for c in BIN:
    sub = df[[c, MAIN_Y]].dropna()
    sub = sub[sub[c].isin([0, 1])]
    if sub[c].nunique() < 2 or len(sub) < 10:
        continue
    a = sub.loc[sub[c] == 1, MAIN_Y]; b = sub.loc[sub[c] == 0, MAIN_Y]
    if min(len(a), len(b)) < 5:
        continue
    u, p = stats.mannwhitneyu(a, b, alternative="two-sided")
    # Cliff's delta
    cd = 2 * u / (len(a) * len(b)) - 1
    defuni.append(dict(变量=c.replace("方法_", ""), 类型="二分类",
                        n=f"{len(a)} vs {len(b)}", 效应量=f"δ={cd:+.3f} "
                        f"(中位 {a.median():+.4f} vs {b.median():+.4f})", p=p))

if THEME_COL:
    g = df[[THEME_COL, MAIN_Y]].dropna()
    keep = g[THEME_COL].value_counts()
    keep = keep[keep >= 8].index                # 太小的类别不做检验
    g = g[g[THEME_COL].isin(keep)]
    groups = [x[MAIN_Y].values for _, x in g.groupby(THEME_COL)]
    if len(groups) >= 3:
        h, p = stats.kruskal(*groups)
        meds = g.groupby(THEME_COL)[MAIN_Y].median().sort_values(ascending=False)
        defuni.append(dict(变量="生物学主主题", 类型=f"多分类({len(groups)}类)",
                            n=len(g), 效应量="; ".join(f"{k}:{v:+.3f}" for k, v in meds.head(5).items()),
                            p=p))

uni = pd.DataFrame(defuni)
if len(uni):
    uni = uni.dropna(subset=["p"])
    uni = uni[uni["p"].notna()]
    uni["q"] = stats.false_discovery_control(uni["p"].values, method="bh")
    uni = uni.sort_values("p")
    uni.to_csv(os.path.join(OUT, "Part2_单因素.csv"), index=False, encoding="utf-8-sig")

emit("-" * 78)
emit("一、单因素关联（BH 校正后，q<0.05 为显著）")
emit("-" * 78)
show = uni[["变量", "类型", "n", "效应量", "p", "q"]].copy()
show["p"] = show["p"].map(lambda v: f"{v:.3g}")
show["q"] = show["q"].map(lambda v: f"{v:.3g}")
emit(show.to_string(index=False))
sigz = uni[uni["q"] < 0.05]
emit("")
if len(sigz):
    emit(f"  显著变量 {len(sigz)} 个：" + "、".join(sigz["变量"].tolist()))
else:
    emit("  没有任何变量在 BH 校正后仍然显著 —— 这本身就是一个值得报告的发现：")
    emit("  已发表签名的跨队列稳健性，无法用常见的文献/设计/方法特征预测。")
emit("")

# ------------------------------------------------------------- 多因素 OLS
emit("-" * 78)
emit("二、多因素模型（OLS + HC1 稳健标准误）")
emit("-" * 78)
CAND = [c for c, _ in CONT] + BIN
rows = []
for c in CAND:
    s = df[[c, MAIN_Y]].dropna()
    if len(s) < 100:
        continue
    try:
        rho, p = stats.spearmanr(s[c], s[MAIN_Y])
    except Exception:
        continue
    if p < 0.15 and not np.isnan(p):
        rows.append((c, abs(rho)))
rows.sort(key=lambda x: -x[1])
SEL = [c for c, _ in rows[:6]]
emit(f"纳入变量（单因素 p<0.15 中效应最强的 6 个）：{SEL}")

Xdf = df[SEL + [MAIN_Y]].dropna()
emit(f"完整病例（无缺失）：{len(Xdf)} / {len(df)}")
multi = []
if len(Xdf) >= 60 and len(SEL) >= 2:
    try:
        import statsmodels.api as sm
        X = sm.add_constant(Xdf[SEL].astype(float))
        fit = sm.OLS(Xdf[MAIN_Y].astype(float), X).fit(cov_type="HC1")
        for v in fit.params.index:
            if v == "const":
                continue
            multi.append(dict(变量=v, 系数=fit.params[v], 稳健SE=fit.bse[v],
                              t=fit.tvalues[v], p=fit.pvalues[v]))
        multi = pd.DataFrame(multi)
        if len(multi):
            multi["q"] = stats.false_discovery_control(multi["p"].values, method="bh")
            multi.to_csv(os.path.join(OUT, "Part2_多因素.csv"), index=False, encoding="utf-8-sig")
            mm = multi.copy()
            for c in ["系数", "稳健SE", "t", "p", "q"]:
                mm[c] = mm[c].map(lambda v: f"{v:.4g}")
            emit(mm.to_string(index=False))
        emit(f"\n  模型整体 R² = {fit.rsquared:.3f}（调整后 {fit.rsquared_adj:.3f}），"
             f"F 检验 p = {fit.f_pvalue:.3g}")
    except Exception as e:
        emit("  多因素模型失败：" + str(e))
else:
    emit("  样本不足以拟合多因素模型，跳过。")
emit("")

# ------------------------------------------------- 三、最稳健 vs 最差的描述
emit("-" * 78)
emit("三、最稳健（6/6）与最不稳（≤1/6）两组签名的特征对比")
emit("-" * 78)
hi = df[df["beats"] == 6]
lo = df[df["beats"] <= 1]
emit(f"  6/6 组 {len(hi)} 个    ≤1/6 组 {len(lo)} 个")
cmp_rows = []
for col, lab in CONT:
    if col not in df.columns:
        continue
    a = hi[col].dropna(); b = lo[col].dropna()
    if len(a) < 5 or len(b) < 5:
        continue
    _, p = stats.mannwhitneyu(a, b, alternative="two-sided")
    cmp_rows.append(dict(特征=lab, 稳健组中位=f"{a.median():.3g}",
                         不稳组中位=f"{b.median():.3g}", p=f"{p:.3g}"))
for c in BIN:
    if c not in df.columns:
        continue
    a = hi[c].dropna(); b = lo[c].dropna()
    if a.nunique() < 2 and b.nunique() < 2:
        continue
    try:
        _, p = stats.mannwhitneyu(a, b, alternative="two-sided")
    except Exception:
        continue
    cmp_rows.append(dict(特征=c.replace("方法_", ""),
                         稳健组中位=f"{a.mean():.1%}", 不稳组中位=f"{b.mean():.1%}",
                         p=f"{p:.3g}"))
cr = pd.DataFrame(cmp_rows)
if len(cr):
    cr["q"] = [f"{v:.3g}" for v in stats.false_discovery_control(
        pd.to_numeric(cr["p"], errors="coerce").values, method="bh")]
    emit(cr.to_string(index=False))

with open(os.path.join(OUT, "Part2_归因分析.txt"), "w", encoding="utf-8") as f:
    f.write("\n".join(lines))

# ---------------------------------------------------------------- Word
try:
    from docx import Document
    from docx.shared import Pt
    doc = Document()
    doc.styles["Normal"].font.name = "Times New Roman"
    doc.styles["Normal"].font.size = Pt(10.5)
    doc.add_heading("Part 2：什么样的签名更能打？——虚拟人群归因分析", 0)
    p = doc.add_paragraph()
    p.add_run(f"{len(df)} 个签名 × 6 个独立队列；结局为跨队列平均 ΔC-index").bold = True

    doc.add_heading("一、单因素关联", 1)
    tb = doc.add_table(rows=1, cols=6); tb.style = "Light Grid Accent 1"
    for i, c in enumerate(["变量", "类型", "n", "效应量", "p", "q(BH)"]):
        cell = tb.rows[0].cells[i]; cell.text = ""
        r = cell.paragraphs[0].add_run(c); r.bold = True; r.font.size = Pt(8.5)
    for _, row in uni.iterrows():
        cs = tb.add_row().cells
        vals = [row["变量"], row["类型"], str(row["n"]), row["效应量"],
                f"{row['p']:.3g}", f"{row['q']:.3g}"]
        for i, v in enumerate(vals):
            cs[i].text = ""
            cs[i].paragraphs[0].add_run(str(v)).font.size = Pt(8.5)
    doc.add_paragraph()
    if len(sigz):
        doc.add_paragraph(f"BH 校正后仍显著的变量：" + "、".join(sigz["变量"].tolist()))
    else:
        doc.add_paragraph(
            "没有任何文献/设计/方法层面的特征在 BH 校正后与跨队列稳健性显著相关。"
            "这意味着：无法在阅读论文时预判一个签名能否复现，"
            "唯一可靠的判断依据仍然是外部数据上的实际重算。").runs[0].bold = True

    doc.add_heading("二、多因素模型", 1)
    for s in lines[lines.index("二、多因素模型（OLS + HC1 稳健标准误）")
                   if "二、多因素模型（OLS + HC1 稳健标准误）" in lines else 0:][:14]:
        doc.add_paragraph(s)
    doc.save(os.path.join(OUT, "Part2_归因分析.docx"))
    emit("\n已生成 Word：Part2_归因分析.docx")
except Exception as e:
    emit("Word 生成失败：" + str(e))

emit("\n输出：Part2_单因素.csv / Part2_多因素.csv / Part2_归因分析.txt|docx")
