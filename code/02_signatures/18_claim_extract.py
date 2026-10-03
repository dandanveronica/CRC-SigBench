# -*- coding: utf-8 -*-
"""
18_claim_extract.py
====================
Part 2 的核心：从 320 篇全文 XML 中抽取【原文声称的预测性能】，用于量化发表偏倚。

与 v1 的关键改动
----------------
v1 用「句子」做窗口，召回只有 71/320 —— 因为大量文章把数值写在表格或下一句。
v2 改为**块级窗口**（以指标词为中心 ±320 字符），并同时抽取 AUC 及其年份，
这样可与本研究重算的 auc5（5 年时间依赖 AUC）配对，样本量从 16 提升到百余。

输出的主分析指标
----------------
  claim_c        原文 C-index（取所有候选的中位）
  claim_auc5     原文 5 年 AUC（与本研究 auc5 直接配对）—— 发表偏倚主终点
  claim_auc3     原文 3 年 AUC
  claim_auc1     原文 1 年 AUC
  claim_auc_any  原文任意 AUC

发表偏倚定义
------------
  Δ = 原文声称值 − 本研究在独立外部队列上的重算值
"""
import re, os, glob
import pandas as pd
import numpy as np

BASE = os.path.dirname(os.path.abspath(__file__))
XMLD = os.path.join(BASE, "fulltext_xml")
OUT = os.path.join(BASE, "audit_out_clean")
os.makedirs(OUT, exist_ok=True)

WIN = 320          # 指标词前后各取多少字符作为上下文窗口
LO, HI = 0.45, 0.995


def protect(t):
    """把小数整体占位，避免被句子切分/正则误伤。"""
    store = {}

    def rp(m):
        k = "NUM%dNUM" % len(store)
        store[k] = m.group(0)
        return k
    t = re.sub(r"\d+\.\d+", rp, t)
    return t, store


def restore(s, store):
    for k, v in store.items():
        s = s.replace(k, v)
    return s


CI_WORD = r"(?:C-?index|C\s+index|concordance\s+index|Harrell'?s?\s+C|C-?statistic)"
AUC_WORD = r"(?:AUC|area\s+under\s+(?:the\s+)?(?:ROC|curve))"
TRAIN_W = r"(?:training|discovery|derivation|internal|TCGA|entire|whole|primary\s+cohort)"
VALID_W = r"(?:validation|validating|testing|external|independent|verification|GEO|GSE\d+|ICGC)"
NUMP = re.compile(r"NUM\d+NUM")
# 年份标记：既覆盖 "5-year"，也覆盖缩写连排 "1-, 3-, and 5-year" 中的 "3-"
YEAR_FULL = re.compile(r"(\d)\s*[-–]\s*(?:year|yr)", re.I)
YEAR_SEQ = re.compile(
    r"(\d)\s*[-–]\s*,?\s*(?:(\d)\s*[-–]\s*,?\s*)?(?:and\s*)?(\d)\s*[-–]\s*(?:year|yr)", re.I)


def _unpack_seq(m):
    return [int(g) for g in m.groups() if g]


SPAN = 250   # 年份标记与数值之间允许的最大间隔


def bind_years(win, store):
    """
    把窗口内的 AUC 数值按年份绑定，返回 (dict{年份: 值}, 绑定方式)。

    最常见也最容易错的句式：
      "1-, 3-, and 5-year AUCs of 0.698, 0.711, and 0.645, respectively"
      -> 必须【按位置顺序一一对应】；取中位数会错把 3 年值当成 5 年值。
      "0.6215, 0.6313, and 0.6715 for 1-, 3-, and 5-year survival"（数字在前）
      -> 同样按序对应。
      "5-year AUC = 0.743"  -> 单年份显式绑定。

    绑定成立的条件：年份标记 ±SPAN 字符内的数值个数 **严格等于** 年份个数。
    数量对不上就放弃（宁可漏，不可错）。
    """
    vals = []
    for m in NUMP.finditer(win):
        v = store.get(m.group(0))
        if v is None:
            continue
        try:
            x = float(v)
        except (TypeError, ValueError):
            continue
        if LO <= x <= HI:
            vals.append((m.start(), m.end(), x))
    if not vals:
        return {}, "none"

    # --- 年份标记 ---
    # (a) 缩写连排："1-, 3-, and 5-year"
    marks = []            # (start, end, [年份...])
    spans = []
    for m in YEAR_SEQ.finditer(win):
        ys = _unpack_seq(m)
        if len(ys) >= 2:
            marks.append((m.start(), m.end(), ys))
            spans.append((m.start(), m.end()))
    # (b) 全称连排："1-year, 3-year, and 5-year"
    #     每个都是独立匹配，必须先把相邻的合并成一个序列，否则会各自绑到同一个数。
    singles = [(m.start(), m.end(), [int(m.group(1))])
               for m in YEAR_FULL.finditer(win)
               if not any(s <= m.start() <= e for s, e in spans)]
    if singles:
        groups = [[singles[0]]]
        for cur in singles[1:]:
            if cur[0] - groups[-1][-1][1] < 80:
                groups[-1].append(cur)
            else:
                groups.append([cur])
        for g in groups:
            ys = [y for _, _, yy in g for y in yy]
            marks.append((g[0][0], g[-1][1], ys))
    if not marks:
        return {}, "noyear"
    marks.sort(key=lambda t: t[0])

    out, hows = {}, []
    for ps, pe, ys in marks:
        after = [x for a, b, x in vals if pe <= a <= pe + SPAN]
        if len(after) == len(ys):
            for y, x in zip(ys, after):
                out.setdefault(y, x)
            hows.append("seq_post")
            continue
        before = [x for a, b, x in vals if ps - SPAN <= b <= ps]
        if len(before) == len(ys):
            for y, x in zip(ys, before):
                out.setdefault(y, x)
            hows.append("seq_pre")
            continue
        if len(ys) == 1:
            cand = [(abs(a - pe), x) for a, b, x in vals if abs(a - pe) <= 150]
            if cand:
                out.setdefault(ys[0], min(cand)[1])
                hows.append("one_nearest")
    return out, ("+".join(sorted(set(hows))) if hows else "unbound")


def ctx_of(win_low):
    mt = re.search(TRAIN_W, win_low)
    mv = re.search(VALID_W, win_low)
    if mt and mv:
        return "train" if mt.start() < mv.start() else "valid"
    if mt:
        return "train"
    if mv:
        return "valid"
    return "unk"


def _scan(seg, store, conf):
    """在一个片段（句子或窗口）里抽性能指标。conf 标记置信度。"""
    recs = []
    for kind, word in (("C-index", CI_WORD), ("AUC", AUC_WORD)):
        for m in re.finditer(word, seg, re.I):
            if kind == "AUC":
                ymap, how = bind_years(seg, store)
                for y, x in ymap.items():
                    pos = [nm.start() for nm in NUMP.finditer(seg)
                           if store.get(nm.group(0))
                           and abs(float(store[nm.group(0)]) - x) < 1e-9]
                    dist = min([abs(p - m.start()) for p in pos], default=999)
                    recs.append(dict(metric="AUC", value=float(x),
                                     ctx=ctx_of(seg.lower()), year=str(y),
                                     bind=how, dist=dist, conf=conf,
                                     sent=restore(seg, store)[:240]))
            else:
                for nm in NUMP.finditer(seg):
                    v = store.get(nm.group(0))
                    if v is None:
                        continue
                    try:
                        x = float(v)
                    except (TypeError, ValueError):
                        continue
                    if not (LO <= x <= HI):
                        continue
                    recs.append(dict(metric="C-index", value=x,
                                     ctx=ctx_of(seg.lower()), year="",
                                     bind="", dist=abs(nm.start() - m.start()),
                                     conf=conf,
                                     sent=restore(seg, store)[:240]))
    return recs


def parse_one(path, want_detail=False):
    raw = open(path, encoding="utf-8", errors="ignore").read()
    txt = re.sub(r"<[^>]+>", " ", raw)
    txt = re.sub(r"&[a-z]+;", " ", txt)
    txt = re.sub(r"\s+", " ", txt)
    txt, store = protect(txt)

    # ---- 第一轮：句子级（置信度高）----
    # 数字已占位，句子切分不会把 "0.783" 切坏。
    recs = []
    for s in re.split(r"(?<=[.;:])\s+", txt):
        s = s.strip()
        if len(s) < 8:
            continue
        if re.search(r"\bp\s*[<>=]\s*NUM\d+NUM", s, re.I):
            continue          # p 值句整句丢弃
        if not (re.search(CI_WORD, s, re.I) or re.search(AUC_WORD, s, re.I)):
            continue
        recs.extend(_scan(s, store, "sent"))

    # ---- 第二轮：仅当句子级一无所获时，才用跨句窗口兜底（置信度低）----
    if not recs:
        for m in re.finditer(CI_WORD + "|" + AUC_WORD, txt, re.I):
            a = max(0, m.start() - WIN)
            b = min(len(txt), m.end() + WIN)
            win = txt[a:b]
            if re.search(r"\bp\s*[<>=]\s*NUM\d+NUM", win, re.I):
                continue
            recs.extend(_scan(win, store, "win"))
    return recs


def best(recs, metric, ctx=None, year=None):
    """在同类中挑最可信的一个：优先距离指标词最近的。"""
    sel = [r for r in recs if r["metric"] == metric]
    if ctx:
        sel = [r for r in sel if r["ctx"] == ctx]
    if year is not None:
        sel = [r for r in sel if year in str(r.get("year", ""))]
    if not sel:
        return np.nan
    sel.sort(key=lambda r: r["dist"])
    return float(sel[0]["value"])


def _vals(recs, metric, ctx=None, year=None):
    """优先取【同句】命中的候选；同句没有时才退回跨句窗口。"""
    sel = [r for r in recs if r["metric"] == metric]
    if ctx:
        sel = [r for r in sel if r["ctx"] == ctx]
    if year is not None:
        sel = [r for r in sel if str(r.get("year", "")) == str(year)]
    hi = [r["value"] for r in sel if r.get("conf") == "sent"]
    if hi:
        return hi
    return [r["value"] for r in sel]


def median_of(recs, metric, ctx=None, year=None):
    v = _vals(recs, metric, ctx, year)
    return float(np.median(v)) if len(v) else np.nan


def max_of(recs, metric, ctx=None, year=None):
    v = _vals(recs, metric, ctx, year)
    return float(max(v)) if len(v) else np.nan


def summarize(recs):
    if not recs:
        return dict(claim_c=np.nan, claim_c_train=np.nan, claim_c_valid=np.nan,
                    claim_auc5=np.nan, claim_auc3=np.nan, claim_auc1=np.nan,
                    claim_auc_any=np.nan, n_ci=0, n_auc=0,
                    report_c=0, report_auc=0)
    ci = [r for r in recs if r["metric"] == "C-index"]
    au = [r for r in recs if r["metric"] == "AUC"]
    return dict(
        claim_c=median_of(recs, "C-index"),
        claim_c_train=median_of(recs, "C-index", ctx="train"),
        claim_c_valid=median_of(recs, "C-index", ctx="valid"),
        claim_auc5=median_of(recs, "AUC", year="5"),
        claim_auc3=median_of(recs, "AUC", year="3"),
        claim_auc1=median_of(recs, "AUC", year="1"),
        claim_auc_any=median_of(recs, "AUC"),
        claim_auc5_max=max_of(recs, "AUC", year="5"),
        claim_auc5_valid=median_of(recs, "AUC", year="5", ctx="valid"),
        n_ci=len(ci), n_auc=len(au),
        report_c=int(len(ci) > 0), report_auc=int(len(au) > 0),
        conf=("sent" if any(r.get("conf") == "sent" for r in recs) else "win"),
    )


def main():
    meta = pd.read_csv(os.path.join(OUT, "全文方法学提取.csv"))
    m = meta[["sig_id", "pmcid"]].drop_duplicates("sig_id")

    rows, det = [], []
    files = sorted(glob.glob(os.path.join(XMLD, "*.xml")))
    print("全文 XML：%d 篇" % len(files))
    for f in files:
        pmc = os.path.basename(f).replace(".xml", "")
        try:
            recs = parse_one(f)
        except Exception:
            recs = []
        s = summarize(recs)
        s["pmcid"] = pmc
        rows.append(s)
        for r in recs[:12]:
            r = dict(r); r["pmcid"] = pmc
            det.append(r)

    df = pd.DataFrame(rows)
    print("报告 C-index 的篇数：%d (%.1f%%)" % (df.report_c.sum(), 100 * df.report_c.mean()))
    print("报告 AUC    的篇数：%d (%.1f%%)" % (df.report_auc.sum(), 100 * df.report_auc.mean()))

    df = m.merge(df, on="pmcid", how="left").drop(columns=["pmcid"])

    per = pd.read_csv(os.path.join(OUT, "per_cohort.csv"))
    per = per[per.status_r == "ok"]
    recalc = per.groupby("sig_id").agg(
        重算C中位=("C", "median"), 重算C均值=("C", "mean"),
        重算AUC5中位=("auc5", "median"), 重算AUC5均值=("auc5", "mean"),
        重算AUC5最大=("auc5", "max"), 重算AUC5最小=("auc5", "min"),
        重算校准斜率=("calib_slope", "median"),
        有效队列数=("C", "count"),
    ).reset_index()

    out = df.merge(recalc, on="sig_id", how="left")
    out["偏倚_C"] = out["claim_c"] - out["重算C中位"]
    out["偏倚_AUC5"] = out["claim_auc5"] - out["重算AUC5中位"]
    out["偏倚_AUC任意"] = out["claim_auc_any"] - out["重算AUC5中位"]

    p = os.path.join(OUT, "原文声称值vs重算值.csv")
    out.to_csv(p, index=False, encoding="utf-8-sig")
    dd = pd.DataFrame(det)
    if len(dd):
        dd = m.merge(dd, on="pmcid", how="right")
        dd.to_csv(os.path.join(OUT, "原文声称值_抽取明细.csv"), index=False, encoding="utf-8-sig")

    print("\n=== 发表偏倚 ===")
    for lab, col in [("C-index（原文 vs 重算）", "偏倚_C"),
                     ("5年AUC（原文 vs 重算）", "偏倚_AUC5"),
                     ("任意AUC（原文 vs 重算5年）", "偏倚_AUC任意")]:
        v = out[col].dropna()
        if not len(v):
            print("%-24s 无可配对" % lab); continue
        from scipy import stats
        try:
            w = stats.wilcoxon(v)
            pv = w.pvalue
        except Exception:
            pv = np.nan
        print("%-24s n=%3d  中位 %+.4f  四分位 %+.4f~%+.4f  >0 占 %.1f%%  Wilcoxon p=%.3g"
              % (lab, len(v), v.median(), v.quantile(.25), v.quantile(.75),
                 100 * (v > 0).mean(), pv))
    print("\n已写出：%s" % p)
    return out


if __name__ == "__main__":
    main()
