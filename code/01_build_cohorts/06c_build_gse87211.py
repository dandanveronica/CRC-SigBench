# -*- coding: utf-8 -*-
"""
06c_build_gse87211.py
=====================
构建 GSE87211 队列（Agilent GPL13497，直肠癌新辅助放化疗，363 例）。

这个队列之前搁置的原因：GEO 上没有 GPL13497 的 annot.gz（下载到的是 0 字节）。
改用 NCBI 的 platform SOFT（view=full）拿 ID -> GENE_SYMBOL 映射解决。

临床字段（逐格解析，因为 characteristics 是锯齿状的）
----------------------------------------------------
  OS : survival time (month)      + death due to tumor (1=事件)
  DFS: disease free time (month)  + cancer recurrance after surgery (1=事件)
  注意：disease free time 存在负值（-4.08），必须过滤，否则 Cox 会报错。

输出：cohorts_raw/GSE87211_expr.csv.gz + GSE87211_clin.csv（供 R 转 rds）
"""
import gzip, os, re
import numpy as np
import pandas as pd

BASE = os.path.dirname(os.path.abspath(__file__))
RAW = os.path.join(BASE, "cohorts_raw")
MAT = os.path.join(RAW, "GSE87211_series_matrix.txt.gz")
GPL = os.path.join(RAW, "GPL13497_platform.txt")


# ---------- 1. 平台注释：ID -> GENE_SYMBOL ----------
def load_gpl():
    if not os.path.exists(GPL) or os.path.getsize(GPL) < 1000:
        raise SystemExit("缺少平台注释文件：%s（需先从 NCBI 下载 GPL13497 view=full）" % GPL)
    cols = None
    rows = []
    started = False
    with open(GPL, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.rstrip("\n").rstrip("\r")
            if line.startswith("!platform_table_begin"):
                started = True
                continue
            if line.startswith("!platform_table_end"):
                break
            if not started:
                continue
            parts = line.split("\t")
            if cols is None:
                cols = parts
                continue
            rows.append(parts[:len(cols)])
    t = pd.DataFrame(rows, columns=cols)
    print("平台表：%d 探针，列：%s" % (len(t), [c for c in cols][:8]))
    if "GENE_SYMBOL" not in t.columns:
        raise SystemExit("平台表缺少 GENE_SYMBOL 列")
    t = t[["ID", "GENE_SYMBOL"]].copy()
    t["GENE_SYMBOL"] = t["GENE_SYMBOL"].astype(str).str.strip()
    t = t[(t.GENE_SYMBOL != "") & (t.GENE_SYMBOL.str.lower() != "nan")]
    # Agilent 常出现一个探针对多个 symbol，取第一个
    t["GENE_SYMBOL"] = t["GENE_SYMBOL"].str.split(r"///").str[0].str.strip()
    t = t[t.GENE_SYMBOL != ""]
    print("有基因符号的探针：%d，唯一基因：%d" % (len(t), t.GENE_SYMBOL.nunique()))
    return dict(zip(t.ID, t.GENE_SYMBOL))


# ---------- 2. 读 series matrix ----------
def read_matrix():
    meta = {}
    chars = []
    with gzip.open(MAT, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("!Sample_characteristics_ch1"):
                chars.append(line.rstrip("\n").split("\t")[1:])
            elif line.startswith("!Sample_geo_accession"):
                meta["sample"] = [x.strip('"') for x in line.rstrip("\n").split("\t")[1:]]
            elif line.startswith("!Sample_platform_id"):
                meta["platform"] = [x.strip('"') for x in line.rstrip("\n").split("\t")[1:]]
            elif line.startswith("!series_matrix_table_begin"):
                break
    return meta, chars


def read_expr():
    """series matrix 的表达矩阵部分：跳过所有以 ! 开头的行。"""
    rows = []
    header = None
    with gzip.open(MAT, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("!"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            if header is None:
                # 表头必须是列数最多的那一行（前面可能残留短行）
                header = [p.strip('"') for p in parts]
                ncol = len(header)
                continue
            if len(parts) != ncol:
                parts = (parts + [""] * ncol)[:ncol]
            rows.append(parts)
    df = pd.DataFrame(rows, columns=header).set_index(header[0])
    df.index = df.index.astype(str).str.strip('"')
    df.columns = [c.strip('"') for c in df.columns]
    return df


def main():
    p2g = load_gpl()
    meta, chars = read_matrix()
    samples = meta["sample"]
    print("样本数：%d，平台：%s" % (len(samples), meta.get("platform", ["?"])[0]))

    expr = read_expr()
    print("原始矩阵：%d 探针 × %d 样本" % expr.shape)
    expr = expr.reindex(columns=samples)
    expr = expr.apply(pd.to_numeric, errors="coerce")

    # 探针 -> 基因符号
    sym = pd.Series([p2g.get(i, "") for i in expr.index], index=expr.index)
    expr = expr[sym != ""]
    expr.index = sym[sym != ""].values
    # 多探针对应同一基因：取均值
    expr = expr.groupby(level=0).mean()
    print("映射后：%d 个基因 × %d 样本" % expr.shape)

    # ---- 逐格解析临床 ----
    kv = {}
    for row in chars:
        for j, cell in enumerate(row):
            if j >= len(samples):
                continue
            cell = str(cell).strip().strip('"')
            if ":" not in cell:
                continue
            k, v = cell.split(":", 1)
            kv.setdefault(k.strip().lower(), {})[samples[j]] = v.strip()
    print("\n可用临床键：%s" % list(kv.keys()))

    def get(key):
        # 必须返回 list（位置对齐），不能返回以 sample 为索引的 Series：
        # clin 的索引是 0..n-1，按索引赋值会全部对不上变成 NaN。
        d = kv.get(key, {})
        return [d.get(s, np.nan) for s in samples]

    clin = pd.DataFrame({"sample": samples})
    clin["time"] = pd.to_numeric(get("survival time (month)"), errors="coerce")
    clin["status"] = pd.to_numeric(get("death due to tumor"), errors="coerce")
    clin["dfs_time"] = pd.to_numeric(get("disease free time (month)"), errors="coerce")
    clin["dfs_status"] = pd.to_numeric(get("cancer recurrance after surgery"), errors="coerce")
    clin["age"] = pd.to_numeric(get("age"), errors="coerce")
    g = pd.Series(get("gender")).astype(str).str.upper().str[0]
    clin["sex"] = g.map({"M": "M", "F": "F"}).fillna("NA")
    # 分期：用 rct 前的浸润深度 + 淋巴结 + 转移合成 TNM
    t = pd.to_numeric(get("depth of invasion before rct"), errors="coerce")
    n = pd.to_numeric(get("lymph node metastasis before rct"), errors="coerce")
    m = pd.to_numeric(get("metastasis before rct"), errors="coerce")
    stage = np.where(pd.Series(m).fillna(0) >= 1, 4,
                     np.where(pd.Series(n).fillna(0) >= 1, 3,
                              np.where(pd.Series(t).fillna(0) >= 2, 2, 1)))
    clin["stage"] = pd.Series(stage)

    # ---- 清理：非正生存时间必须剔除（glmnet/coxph 不接受）----
    for tt, ss in [("time", "status"), ("dfs_time", "dfs_status")]:
        bad = (pd.to_numeric(clin[tt], errors="coerce") <= 0) & clin[ss].notna()
        if bad.sum():
            print("  %s：剔除 %d 例非正时间" % (tt, int(bad.sum())))
            clin.loc[bad, [tt, ss]] = np.nan
        clin.loc[clin[ss].isna(), tt] = np.nan

    ok = clin.time.notna() & clin.status.notna()
    print("\nOS 可用：%d 例，事件 %d" % (ok.sum(), int(clin.loc[ok, "status"].sum())))
    ok2 = clin.dfs_time.notna() & clin.dfs_status.notna()
    print("DFS可用：%d 例，事件 %d" % (ok2.sum(), int(clin.loc[ok2, "dfs_status"].sum())))

    expr = expr.loc[:, clin["sample"].values]
    # log2 转换判断：Agilent 数据若已 log 化则跳过
    mx = np.nanmax(expr.values)
    if mx > 50:
        print("数值最大 %.1f -> 判定为线性信号，做 log2(x+1)" % mx)
        expr = np.log2(expr.clip(lower=0) + 1)
    else:
        print("数值最大 %.1f -> 已在 log 尺度，不再转换" % mx)

    out_e = os.path.join(RAW, "GSE87211_expr.csv.gz")
    out_c = os.path.join(RAW, "GSE87211_clin.csv")
    expr.to_csv(out_e, compression="gzip")
    clin.to_csv(out_c, index=False, encoding="utf-8-sig")
    print("\n已写出：\n  %s\n  %s" % (out_e, out_c))


if __name__ == "__main__":
    main()
