#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
09_metadata.py —— 构建 Part 2 归因分析用的签名元数据

从标题 / 摘要证据句 / 基因列表抽取可解释变量：
  A. 方法学：LASSO、单因素 Cox、多因素 Cox、逐步回归、机器学习、WGCNA、
             列线图、单细胞/空间组学
  B. 生物学主题：NETs、缺氧、脂代谢、铁死亡、铜死亡、自噬、焦亡、m6A、
                 CAF、免疫、糖酵解、EMT、线粒体、内质网应激、血管生成、
                 外泌体、衰老、干细胞性、失巢凋亡、氨基酸、泛素化等
  C. 基因层面：基因数、基因复用度（该签名基因在全部签名中被使用的频率）
  D. 平台可测量性：跨队列平均检出率

输入：sigminer/output/signature_library_clean.csv
      audit_out_clean/per_cohort.csv（用于检出率）
输出：audit_out_clean/签名元数据_归因变量.csv
"""
import os, re, sys
import pandas as pd, numpy as np

LIB = "sigminer/output/signature_library_clean.csv"
PER = sys.argv[1] if len(sys.argv) > 1 else "audit_out_clean/per_cohort.csv"
OUT = "audit_out_clean"
os.makedirs(OUT, exist_ok=True)

lib = pd.read_csv(LIB, dtype=str, encoding="utf-8-sig").fillna("")
lib["year"] = pd.to_numeric(lib["year"], errors="coerce")
lib["n_gene"] = pd.to_numeric(lib["n_gene"], errors="coerce")

lib["txt"] = (lib["title"] + " || " + lib["证据句"]).str.lower()

# ------------------------------------------------------------------ A 方法学
METHOD = {
    "LASSO":            r"\blasso\b|least absolute shrinkage",
    "单因素Cox":         r"univariate cox|univariable cox",
    "多因素Cox":         r"multivariate cox|multivariable cox|multivariate analysis",
    "逐步回归":          r"stepwise|step-wise|stepwise regression",
    "机器学习":          r"random forest|support vector|\bsvm\b|xgboost|gradient boost|"
                        r"machine learning|deep learning|neural network|boruta",
    "WGCNA":            r"\bwgcna\b|weighted gene co-?expression network",
    "列线图":           r"\bnomogram\b",
    "单细胞/空间":       r"single-?cell|scRNA|spatial transcriptom|single cell RNA",
    "独立外部验证":      r"external validation|independent (?:validation )?cohort|"
                        r"validated in an independent|external (?:test )?(?:set|cohort)|"
                        r"validation cohort|external cohorts",
    "多中心":           r"multicenter|multi-?cent(?:er|re)|multi-?institutional",
    "免疫浸润分析":      r"immune infiltration|immune cell infiltration|CIBERSORT|"
                        r"immune microenvironment|tumor immune",
    "药物敏感性/治疗反应": r"drug sensitivity|immunotherapy response|therapeutic response|"
                          r"treatment response|chemosensitivity|checkpoint",
}

# -------------------------------------------------------------- B 生物学主题
THEME = {
    "NETs中性粒细胞胞外诱捕网": r"neutrophil extracellular trap|\bNETs?\b",
    "缺氧":                  r"hypoxi",
    "脂代谢":                r"lipid metabolism|lipid-?related|fatty acid|cholesterol metabolism",
    "铁死亡":                r"ferroptosis",
    "铜死亡":                r"cuproptosis",
    "自噬":                  r"autophag",
    "焦亡":                  r"pyroptosis",
    "失巢凋亡":              r"anoikis",
    "m6A甲基化":             r"\bm6a\b|n6-methyladenosine|methyladenosine",
    "泛素化/蛋白酶体":        r"ubiquitin|sumoylat|degron",
    "CAF成纤维细胞":         r"cancer-?associated fibroblast|\bCAF\b|fibroblast",
    "免疫/炎症":             r"immun|inflammat|macrophage|t-?cell|b-?cell|lymphocyte|"
                            r"myeloid|neutrophil|cytokine|chemokine",
    "糖酵解/代谢重编程":      r"glycolys|metabolic reprogram|metabolism-?related|metabolic",
    "线粒体":                r"mitochondri",
    "内质网应激":            r"endoplasmic reticulum stress|er stress|unfolded protein",
    "血管生成":              r"angiogenesis|angiogenic",
    "外泌体":                r"exosome|exosomal",
    "细胞衰老":              r"senescen",
    "干细胞性":             r"stemness|stem cell|cancer stem|\bCSC\b",
    "EMT/转移":             r"epithelial-?mesenchymal|\bEMT\b|metastas",
    "氨基酸代谢":            r"glutamine|amino acid|tryptophan|arginine|lysine|glycine|serine",
    "细胞周期/增殖":         r"cell cycle|proliferat|mitotic",
    "DNA损伤修复":          r"dna damage|dna repair|homologous recombination",
    "表观/染色质":           r"chromatin|histone|epigenetic|methylation",
    "内质网/高尔基体":        r"golgi|endoplasmic",
    "迁移体/新细胞器":        r"migrasome|organelle",
    "昼夜/节律":             r"circadian",
}

def flag(col_map, row):
    return {k: int(bool(re.search(v, row))) for k, v in col_map.items()}


# ---------------------------------------------------------- C 基因复用度
gene_lists = [sorted({g.strip() for g in re.split(r"[;,]", s) if g.strip()})
              for s in lib["genes"]]
from collections import Counter
cnt = Counter()
for gl in gene_lists:
    cnt.update(gl)
n_sig = len(lib)

rows = []
for i, row in lib.iterrows():
    gl = gene_lists[i]
    if not gl:
        continue
    freq = np.array([cnt.get(g, 1) for g in gl], dtype=float)
    # 归一化为「出现在多少个其它签名里」的比例（0~1）
    reuse = (freq - 1) / max(n_sig - 1, 1)
    rows.append(dict(
        sig_id=row["sig_id"],
        gene_popularity_mean=round(reuse.mean(), 4),        # 平均复用度
        gene_popularity_max=round(reuse.max(), 4),          # 最热基因的复用度
        pct_hot_gene=round(float((freq >= 10).mean()), 4),  # 高频基因(≥10个签名用过)占比
        pct_unique_gene=round(float((freq <= 2).mean()), 4),# 近乎独有基因占比
    ))
meta_g = pd.DataFrame(rows)
lib = lib.merge(meta_g, on="sig_id", how="left")

# ------------------------------------------------------ D 打方法/主题标签
mm = pd.DataFrame([flag(METHOD, t) for t in lib["txt"]])
mm.columns = ["方法_" + c for c in mm.columns]
tt = pd.DataFrame([flag(THEME, t) for t in lib["txt"]])
tt.columns = ["主题_" + c for c in tt.columns]
met = pd.concat([lib[["sig_id", "year", "journal", "title", "n_gene"]].reset_index(drop=True),
                 mm, tt,
                 lib[["gene_popularity_mean", "gene_popularity_max",
                      "pct_hot_gene", "pct_unique_gene"]].reset_index(drop=True)], axis=1)

# 主题互斥化：一个签名可能命中多个主题，另设「主主题」= 第一个命中的（按表序）
theme_cols = [c for c in met.columns if c.startswith("主题_")]
def main_theme(r):
    for c in theme_cols:
        if r[c] == 1:
            return c.replace("主题_", "")
    return "未分类"
met["主主题"] = met.apply(main_theme, axis=1)
met["主题命中数"] = met[theme_cols].sum(axis=1)

# ------------------------------------------------ E 跨平台可测量性（来自重算）
try:
    per = pd.read_csv(PER)
    per = per[per["status_r"] == "ok"]
    meas = per.groupby("sig_id").agg(
        检出率均值=("n_measured", lambda s: np.mean(s / per.loc[s.index, "n_gene"])),
        平台数=("cohort", "nunique"),
    ).reset_index()
    met = met.merge(meas, on="sig_id", how="left")
except Exception as e:
    print("跳过检出率合并：", e)

# ------------------------------------------------------------------ 输出
met.to_csv(os.path.join(OUT, "签名元数据_归因变量.csv"), index=False, encoding="utf-8-sig")
print(f"签名元数据：{met.shape[0]} 行 × {met.shape[1]} 列")
print("\n【方法学要素覆盖率】")
for c in [x for x in met.columns if x.startswith("方法_")]:
    print(f"  {c.replace('方法_',''):<16s} {met[c].sum():>4d} ({met[c].mean():.1%})")
print("\n【生物学主题 Top15】")
print(met[theme_cols].sum().sort_values(ascending=False).head(15).to_string())
print("\n【主主题分布 Top12】")
print(met["主主题"].value_counts().head(12).to_string())
print("\n【基因复用度】")
print(met[["gene_popularity_mean", "gene_popularity_max",
           "pct_hot_gene", "pct_unique_gene"]].describe().round(3).to_string())
print(f"\n检出率：<100% 的签名 {(met['检出率均值'] < 1).sum()} 个"
      f"（即存在部分 Sigpli 基因在 GPL570 上测不到）"
      if "检出率均值" in met.columns else "")
print("\n已写出：", os.path.join(OUT, "签名元数据_归因变量.csv"))
