# -*- coding: utf-8 -*-
"""
19_part2_attribution.py
=======================
Part 2：方法学归因 —— 哪些"做法"与独立队列上的真实表现相关？

结局变量（签名级，跨队列汇总）——【两个都要报，含义不同】
------------------------------------------------------
  dC_rand   签名 C − 随机基因集基线      = "比瞎猜强多少"
  dC_clin   (临床+签名) − 仅临床         = "有 stage 之后还值多少"  ← 审稿人必问
  C_mean    重算 C-index 均值
  beat_rate 优于随机基线的队列比例（0~1）
  n_cohort  有效队列数

  注意：per_cohort.csv 里的 `delta_C` 列是【临床之上增量】，不是"减随机基线"。
  两者 Spearman ρ 只有 0.44，不能互相替代，必须分开报告。

解释变量（来自 320 篇全文 + 签名库）
------------------------------------
  主题约束（主主题是否明确）、单细胞衍生与否、算法、是否声明外部验证、
  是否公开风险公式、是否多因素校正、基因数、发表年份、训练集样本量、
  热门基因占比、是否报告 C-index。

统计
----
  二分组：Mann–Whitney U
  多分组：Kruskal–Wallis
  连续变量：Spearman
  BH 校正（跨全部检验）
"""
import os, re
import numpy as np
import pandas as pd
from scipy import stats

BASE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(BASE, "audit_out_clean")
os.makedirs(OUT, exist_ok=True)

lines = []


def emit(s=""):
    print(s)
    lines.append(str(s))


# ---------------- 载入 ----------------
per = pd.read_csv(os.path.join(OUT, "per_cohort.csv"))
per = per[per.status_r == "ok"].copy()
per["beat"] = (per["C"] > per["null_median"]).astype(int)

sig = per.groupby("sig_id").agg(
    dC_clin=("delta_C", "mean"),                       # 临床之上的增量
    C_mean=("C", "mean"),
    C_median=("C", "median"),
    beat_rate=("beat", "mean"),
    n_cohort=("C", "count"),
).reset_index()
# 随机基线之上的增量：逐记录相减后再按签名汇总
per["d_rand"] = per["C"] - per["null_median"]
sig = sig.merge(per.groupby("sig_id").d_rand.mean().rename("dC_rand").reset_index(),
                on="sig_id", how="left")

meta = pd.read_csv(os.path.join(OUT, "签名元数据_归因变量.csv"))
ft = pd.read_csv(os.path.join(OUT, "全文方法学提取.csv"))[
    ["sig_id", "训练集样本量", "模型构建方法", "是否外部验证",
     "是否多因素校正", "是否给出风险公式"]]
claim = pd.read_csv(os.path.join(OUT, "原文声称值vs重算值.csv"))

d = sig.merge(meta, on="sig_id", how="left").merge(ft, on="sig_id", how="left")
if "偏倚_AUC5" in claim.columns:
    d = d.merge(claim[["sig_id", "偏倚_AUC5", "偏倚_C", "report_c", "report_auc"]],
                on="sig_id", how="left")
d["n_cohort_ok"] = d["n_cohort"] >= 4          # 至少 4 个有效队列才纳入
emit("=" * 78)
emit("Part 2  方法学归因分析")
emit("=" * 78)
emit("可分析签名：%d 个（有效队列 ≥4 的：%d 个）" % (len(d), d.n_cohort_ok.sum()))
d = d[d.n_cohort_ok].copy()

# ---------------- 构造解释变量 ----------------
THEME_NONE = {"", "无", "其他", "未分类", "nan"}


def has_theme(x):
    s = str(x).strip()
    return 0 if (s in THEME_NONE or s.lower() == "nan") else 1


d["X_主题约束"] = d["主主题"].apply(has_theme)
d["X_单细胞衍生"] = d["方法_单细胞/空间"].fillna(0).astype(int)
d["X_LASSO"] = d["方法_LASSO"].fillna(0).astype(int)
d["X_机器学习"] = d["方法_机器学习"].fillna(0).astype(int)
d["X_WGCNA"] = d["方法_WGCNA"].fillna(0).astype(int)
d["X_声明外部验证"] = d["是否外部验证"].fillna(0).astype(int)
d["X_公开风险公式"] = d["是否给出风险公式"].fillna(0).astype(int)
d["X_多因素校正"] = d["是否多因素校正"].fillna(0).astype(int)
d["X_列线图"] = d["方法_列线图"].fillna(0).astype(int)
# 注意：二分类切分要先看分布，中位数切分在取值高度集中时会退化成 100%/0%，
# 检验毫无意义（pct_hot_gene 大量并列、year 集中在 2022+）。这两个改用连续变量。
d["X_大基因集"] = (d["n_gene"] >= 10).astype(int)
d["X_报告Cindex"] = d.get("report_c", pd.Series(np.nan, index=d.index)).fillna(0).astype(int)

emit("\n暴露 prevalence：")
for c in [c for c in d.columns if c.startswith("X_")]:
    emit("  %-18s n=%3d (%.1f%%)" % (c[2:], int(d[c].sum()), 100 * d[c].mean()))

# ---------------- 检验 ----------------
emit("\n" + "-" * 78)
emit("单因素归因（签名级，两个结局分别做）")
emit("-" * 78)
for lab, ycol in [("ΔC_临床之上（临床+签名 − 仅临床）", "dC_clin"),
                  ("ΔC_随机之上（签名 − 随机基线）", "dC_rand")]:
    v = d[ycol].dropna()
    emit("  【%s】n=%d  中位 %+.4f  四分位 %+.4f ~ %+.4f  >0 占 %.1f%%  p=%.2g"
         % (lab, len(v), v.median(), v.quantile(.25), v.quantile(.75),
            100 * (v > 0).mean(), stats.wilcoxon(v).pvalue))
emit("")
bins = [c for c in d.columns if c.startswith("X_")]

rows = []
for YMAIN, ylab in [("dC_clin", "临床之上"), ("dC_rand", "随机之上")]:
    for c in bins:
        g1 = d.loc[d[c] == 1, YMAIN].dropna()
        g0 = d.loc[d[c] == 0, YMAIN].dropna()
        if len(g1) < 8 or len(g0) < 8:
            continue
        u, p = stats.mannwhitneyu(g1, g0, alternative="two-sided")
        rows.append(dict(结局=ylab, 变量=c[2:], n1=len(g1), n0=len(g0),
                         中位1=g1.median(), 中位0=g0.median(),
                         差值=g1.median() - g0.median(), p=p, 类型="二分"))

# 多分组（只打印一次明细）
YMAIN = "dC_clin"
for c, lab in [("模型构建方法", "算法组合")]:
    sub = d[[c, YMAIN]].dropna()
    top = sub[c].value_counts()
    keep = top[top >= 12].index
    sub = sub[sub[c].isin(keep)]
    if sub[c].nunique() >= 3:
        groups = [g[YMAIN].values for _, g in sub.groupby(c)]
        h, p = stats.kruskal(*groups)
        rows.append(dict(结局="临床之上", 变量=lab, n1=len(sub), n0=int(sub[c].nunique()),
                         中位1=np.nan, 中位0=np.nan, 差值=np.nan,
                         p=p, 类型="多组"))
        emit("\n  【%s】各组中位 ΔC（临床之上）：" % lab)
        tb = sub.groupby(c)[YMAIN].agg(["count", "median"]).sort_values("median", ascending=False)
        for k, r in tb.iterrows():
            emit("     %-34s n=%3d  ΔC中位 %+.4f" % (str(k)[:34], int(r["count"]), r["median"]))

# 连续变量
for c, lab in [("n_gene", "基因数"), ("year", "发表年份"),
               ("训练集样本量", "训练集样本量"),
               ("pct_hot_gene", "热门基因占比"),
               ("gene_popularity_mean", "基因平均流行度")]:
    if c not in d.columns:
        continue
    sub = d[[c, YMAIN]].dropna()
    if len(sub) < 30:
        continue
    sub2 = d[[c, "dC_rand"]].dropna()
    rho, p = stats.spearmanr(sub[c], sub[YMAIN])
    rows.append(dict(结局="临床之上", 变量=lab, n1=len(sub), n0=np.nan, 中位1=np.nan,
                     中位0=np.nan, 差值=rho, p=p, 类型="连续(Spearman ρ)"))
    if len(sub2) >= 30:
        rho2, p2 = stats.spearmanr(sub2[c], sub2["dC_rand"])
        rows.append(dict(结局="随机之上", 变量=lab, n1=len(sub2), n0=np.nan, 中位1=np.nan,
                         中位0=np.nan, 差值=rho2, p=p2, 类型="连续(Spearman ρ)"))

res = pd.DataFrame(rows)
if len(res):
    # BH 校正
    pp = res["p"].values
    order = np.argsort(pp)
    ranked = pp[order]
    m = len(pp)
    q = ranked * m / (np.arange(m) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    outq = np.empty(m)
    outq[order] = np.clip(q, 0, 1)
    res["q_BH"] = outq
    res = res.sort_values("p")
    emit("\n  检验汇总（按 p 排序）：")
    for _, r in res.iterrows():
        star = "***" if r["p"] < .001 else ("**" if r["p"] < .01 else ("*" if r["p"] < .05 else ""))
        if r["类型"] == "连续(Spearman ρ)":
            emit("    [%-6s] %-14s ρ=%+.3f  n=%3d  p=%.3g q=%.3g %s"
                 % (r["结局"], r["变量"], r["差值"], r["n1"], r["p"], r["q_BH"], star))
        elif r["类型"] == "多组":
            emit("    [%-6s] %-14s Kruskal %d组 n=%3d  p=%.3g q=%.3g %s"
                 % (r["结局"], r["变量"], r["n0"], r["n1"], r["p"], r["q_BH"], star))
        else:
            emit("    [%-6s] %-14s 有 n=%3d(中位%+.4f) vs 无 n=%3d(中位%+.4f) 差%+.4f  p=%.3g q=%.3g %s"
                 % (r["结局"], r["变量"], r["n1"], r["中位1"], r["n0"], r["中位0"],
                    r["差值"], r["p"], r["q_BH"], star))

# ---------------- 发表偏倚的归因 ----------------
if "偏倚_AUC5" in d.columns:
    emit("\n" + "-" * 78)
    emit("发表偏倚的归因（结局：原文声称 AUC5 − 本研究重算 AUC5）")
    emit("-" * 78)
    v = d["偏倚_AUC5"].dropna()
    emit("可配对：%d 个签名   中位 %+.4f   >0 占 %.1f%%"
         % (len(v), v.median(), 100 * (v > 0).mean()))
    r2 = []
    for c in bins:
        g1 = d.loc[d[c] == 1, "偏倚_AUC5"].dropna()
        g0 = d.loc[d[c] == 0, "偏倚_AUC5"].dropna()
        if len(g1) < 6 or len(g0) < 6:
            continue
        u, p = stats.mannwhitneyu(g1, g0, alternative="two-sided")
        r2.append((c[2:], len(g1), g1.median(), len(g0), g0.median(),
                   g1.median() - g0.median(), p))
    r2 = sorted(r2, key=lambda t: t[6])
    emit("\n  %-16s %-28s %-28s %s" % ("变量", "有该特征", "无该特征", "p"))
    for name, n1, m1, n0, m0, diff, p in r2:
        star = "***" if p < .001 else ("**" if p < .01 else ("*" if p < .05 else ""))
        emit("  %-16s n=%3d 偏倚%+.4f      n=%3d 偏倚%+.4f      差%+.4f p=%.3g %s"
             % (name, n1, m1, n0, m0, diff, p, star))

# ---------------- 输出 ----------------
res.to_csv(os.path.join(OUT, "Part2_归因检验.csv"), index=False, encoding="utf-8-sig")
p = os.path.join(OUT, "Part2_归因分析.txt")
open(p, "w", encoding="utf-8").write("\n".join(lines))
emit("\n已写出：%s" % p)
