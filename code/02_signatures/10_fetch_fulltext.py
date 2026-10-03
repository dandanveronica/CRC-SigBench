#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
10_fetch_fulltext.py —— 从 Europe PMC 抓取全文并抽取方法学元数据

为什么需要：原先的归因变量来自摘要片段，方法学标签召回率明显偏低
（LASSO 仅 16.6%、单因素 Cox 仅 6.9%），不足以支撑 Part 2 的方法学归因。
全文才能回答这些问题：
    - 训练集样本量到底是多少
    - 有没有独立的外部验证队列（还是只做 TCGA 内部随机拆分）
    - 用什么方法筛基因建模（LASSO / 逐步 Cox / 随机森林 / RSF）
    - 有没有做多因素 Cox 校正（即"独立预后因子"的证据强度）
    - 有没有公开完整风险计算公式

输入：sigminer/output/signature_library_clean.csv（需 pmcid 列）
输出：fulltext_xml/{pmcid}.xml            原始全文（缓存，避免重复请求）
      audit_out_clean/全文方法学提取.csv    每条签名的结构化方法学字段
"""
import os, re, sys, time, json, concurrent.futures as cf
import xml.etree.ElementTree as ET
import pandas as pd
import requests

LIB = "sigminer/output/signature_library_clean.csv"
XML_DIR = "fulltext_xml"
OUT = "audit_out_clean"
BASE = "https://www.ebi.ac.uk/europepmc/webservices/rest/{}/fullTextXML"
os.makedirs(XML_DIR, exist_ok=True)
os.makedirs(OUT, exist_ok=True)

S = requests.Session()
S.headers.update({"User-Agent": "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) "
                                "AppleWebKit/537.36 (KHTML, like Gecko) Chrome/120 Safari/537.36"})

lib = pd.read_csv(LIB, dtype=str, encoding="utf-8-sig").fillna("")


def pmc_ids(row):
    """尽量把 pmcid 规整成 PMCxxxxxxx"""
    v = str(row.get("pmcid", "") or "").strip()
    if not v:
        m = re.search(r"PMC\d+", str(row.get("url", "")) or "", re.I)
        v = m.group(0) if m else ""
    if not v:
        return ""
    v = v.upper()
    if not v.startswith("PMC"):
        v = "PMC" + v.lstrip("0") if v.isdigit() else v
    return v


def fetch(pid, tries=2):
    fp = os.path.join(XML_DIR, pid + ".xml")
    if os.path.exists(fp) and os.path.getsize(fp) > 5000:
        return pid, "cached", fp
    for t in range(tries):
        try:
            r = S.get(BASE.format(pid), timeout=40)
            if r.status_code == 200 and len(r.content) > 5000 and "<article" in r.text[:2000]:
                with open(fp, "wb") as f:
                    f.write(r.content)
                return pid, "ok", fp
            if r.status_code in (403, 404):
                return pid, f"http{r.status_code}", ""
        except Exception as e:
            pass
        time.sleep(1.5 + 2 * t)
    return pid, "fail", ""


def xml_text(fp):
    """抽 <abstract> + <body>，尽量保留 <title> 层级信息"""
    try:
        tree = ET.parse(fp)
    except Exception:
        try:
            txt = open(fp, encoding="utf-8", errors="ignore").read()
            return re.sub(r"<[^>]+>", " ", txt)
        except Exception:
            return ""
    root = tree.getroot()
    parts = []
    for tag in ("abstract", "body"):
        for node in root.iter(tag):
            for s in node.itertext():
                s = s.strip()
                if s:
                    parts.append(s)
    return re.sub(r"\s+", " ", " ".join(parts))


# ---------------------------------------------------- 从全文抽取结构化字段
N_PAT = [
    r"(?:training|discovery|derivation)\s+(?:cohort|set|group)[^.]{0,120}?\(?\s*n\s*[=:]\s*(\d{2,4})",
    r"(?:TCGA|training)\s+(?:cohort|dataset|database)[^.]{0,80}?(\d{2,4})\s+(?:patients|samples|cases)",
    r"a total of\s+(\d{2,4})\s+(?:CRC|colorectal[^.]{0,30})?\s*(?:patients|samples|cases)",
    r"(\d{2,4})\s+(?:patients|samples|cases)\s+(?:were|was)\s+(?:enrolled|included|recruited|collected)",
    r"(?:comprised|consisted of|included)\s+(\d{2,4})\s+(?:patients|samples|cases)",
]
GEO_PAT = re.compile(r"GSE\d{4,6}")


def extract(rec):
    pid, status, fp = rec
    out = dict(pmcid=pid, 取全文状态=status)
    keys = ["训练集样本量", "模型构建方法", "是否外部验证", "是否多因素校正",
            "是否给出风险公式", "是否有非负对照队列", "GEO队列提及数", "全文字符数"]
    for k in keys:
        out[k] = "" if k != "GEO队列提及数" else 0
    if not fp:
        return out
    txt = xml_text(fp)
    if not txt:
        return out
    low = txt.lower()
    out["全文字符数"] = len(txt)

    # 1) 训练集样本量
    n_found = []
    for p in N_PAT:
        for m in re.finditer(p, low, re.I):
            v = int(m.group(1))
            if 15 <= v <= 5000:
                n_found.append(v)
    if n_found:
        # 取出现最早的（通常是训练集/整体队列）
        out["训练集样本量"] = min(n_found) if len(set(n_found)) == 1 else max(set(n_found), key=n_found.count)

    # 2) 建模方法
    meths = []
    if re.search(r"\blasso\b|least absolute shrinkage", low):            meths.append("LASSO")
    if re.search(r"univariate cox|univariable cox", low):                 meths.append("单因素Cox")
    if re.search(r"multivariate cox|multivariable cox", low):             meths.append("多因素Cox")
    if re.search(r"stepwise|step-wise", low):                             meths.append("逐步回归")
    if re.search(r"random (?:survival )?forest|rsf\b", low):              meths.append("随机森林")
    if re.search(r"support vector|\bsvm\b|xgboost|gradient boost|boruta|"
                 r"neural network|deep learning", low):                   meths.append("机器学习-其他")
    if re.search(r"\bwgcna\b", low):                                      meths.append("WGCNA")
    if re.search(r"pca\b|principal component", low):                      meths.append("PCA")
    out["模型构建方法"] = "+".join(meths) if meths else "未说明"

    # 3) 外部验证：提到 TCGA 之外的独立数据集（GEO、独立医院队列）
    geos = {g.upper() for g in GEO_PAT.findall(txt)}
    out["GEO队列提及数"] = len(geos)
    indep = bool(re.search(
        r"external validation|independent (?:validation |external )?(?:cohort|set|population)|"
        r"validated (?:in|using) (?:an )?(?:independent|external)|"
        r"independent cohort from|prospectively collected|"
        r"our (?:own|institution)(?:al)? cohort|patients from our hospital|"
        r"clinical samples? (?:from|were obtained)", low))
    out["是否外部验证"] = int(indep or len(geos) >= 2)

    # 4) 多因素 Cox 校正
    out["是否多因素校正"] = int(bool(re.search(
        r"multivariate cox|multivariable cox|multivariate analysis", low)))

    # 5) 是否给出可用系数/风险公式
    out["是否给出风险公式"] = int(bool(re.search(
        r"risk score\s*=|riskscore\s*=|\u03a3\s*\(?\s*(?:coef|expression|exp)|"
        r"coefficients?\s+(?:from|of)\s+(?:the\s+)?(?:lasso|multivariate|cox)|"
        r"β\s*(?:coefficient|value)|beta coefficient", low)))
    return out


if __name__ == "__main__":
    sig_ids, pid_list = [], []
    for _, row in lib.iterrows():
        pid = pmc_ids(row)
        sig_ids.append(row["sig_id"])
        pid_list.append(pid)
    have = sum(1 for p in pid_list if p)
    print(f"待抓取 {len(pid_list)} 条，其中有 PMCID 的 {have} 条")

    todo = [(sid, p) for sid, p in zip(sig_ids, pid_list) if p]
    print(f"开始并发抓取（Europe PMC 全文 XML）...")

    t0 = time.time()
    results = {}
    with cf.ThreadPoolExecutor(max_workers=5) as ex:
        futs = {ex.submit(fetch, p): sid for sid, p in todo}
        done = 0
        for fu in cf.as_completed(futs):
            sid = futs[fu]
            results[sid] = fu.result()
            done += 1
            if done % 25 == 0 or done == len(todo):
                el = time.time() - t0
                print(f"  {done}/{len(todo)}  用时 {el:.0f}s  预计剩余 "
                      f"{el/done*(len(todo)-done):.0f}s", flush=True)

    # 抽取结构化字段
    rows = []
    for i, (sid, p) in enumerate(zip(sig_ids, pid_list)):
        if not p:
            continue
        rec = results.get(sid, (p, "fail", ""))
        r = extract(rec)
        r["sig_id"] = sid
        rows.append(r)
    df = pd.DataFrame(rows)
    if len(df):
        for c in ["全文字符数", "是否外部验证", "是否多因素校正",
                  "是否给出风险公式", "GEO队列提及数"]:
            if c in df.columns:
                df[c] = pd.to_numeric(df[c], errors="coerce").fillna(0).astype(int)
        df = df[["sig_id"] + [c for c in df.columns if c != "sig_id"]]
        df.to_csv(os.path.join(OUT, "全文方法学提取.csv"), index=False, encoding="utf-8-sig")

        stat = df["取全文状态"].value_counts()
        print("\n【全文获取】")
        print(stat.to_string())
        ok = df[df["全文字符数"] > 5000]
        print(f"\n可用全文 {len(ok)} / {len(df)} = {len(ok)/len(df):.1%}")
        if len(ok):
            print("\n【方法构成】")
            print(ok["模型构建方法"].value_counts().head(12).to_string())
            print(f"\n有外部独立验证队列：{ok['是否外部验证'].sum()} ({ok['是否外部验证'].mean():.1%})")
            print(f"做了多因素 Cox 校正：{ok['是否多因素校正'].sum()} ({ok['是否多因素校正'].mean():.1%})")
            print(f"公开了可用系数/公式：{ok['是否给出风险公式'].sum()} ({ok['是否给出风险公式'].mean():.1%})")
            nn = pd.to_numeric(ok["训练集样本量"], errors="coerce").dropna()
            print(f"\n训练集样本量：可提取 {len(nn)} 篇，中位 {nn.median():.0f}，"
                  f"范围 {nn.min():.0f}~{nn.max():.0f}，四分位 {nn.quantile(.25):.0f}/{nn.quantile(.75):.0f}")
    print("\n已写出：", os.path.join(OUT, "全文方法学提取.csv"))
