#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
30_tcga_xena.py —— 构建第 8 个队列 TCGA_COADREAD（唯一 RNA-seq 平台）

【为什么放弃 cBioPortal API 逐批拉表达】
实测：10 例 155 秒 / 响应 90MB；25 例 427 秒 / 响应 216MB。
594 例要 3 小时以上，且 JSON 每条记录都重复携带完整 gene 对象 + profileId +
studyId，冗余到近乎荒谬（有用载荷不到十分之一）。

【现在的方案：表达走 Xena S3 静态文件，临床走 cBioPortal API】
  表达  https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/
        TCGA.COADREAD.sampleMap/HiSeqV2.gz       23MB，整个矩阵一个文件
  临床  cBioPortal REST（3 秒返回全 594 例 PETIENT 级临床，含 OS_MONTHS 本来就是【月】）

【两个数据源合并的关键】样品号口径不一致：
  Xena   列名为 sample barcode  TCGA-3L-AA1B-01A（末尾带 A）
  cBio   sampleId 为           TCGA-3L-AA1B-01
  -> 统一规约到 12 位 TCGA-3L-AA1B，再各自映射回去

产出：cohorts_raw/TCGA_expr.tsv.gz（基因×样本）、cohorts_raw/TCGA_clin.csv
      再由 31_tcga_to_rds.R 转成 cohorts/TCGA_COADREAD.rds
"""
import os, sys, time, gzip, threading, re
import requests
import numpy as np
import pandas as pd

OUT = "cohorts_raw"
PARTS = os.path.join(OUT, ".parts_tcga")
os.makedirs(OUT, exist_ok=True)

EXPR_URL = ("https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/"
            "TCGA.COADREAD.sampleMap/HiSeqV2.gz")
TARGET = os.path.join(OUT, "TCGA_HiSeqV2.gz")

S = requests.Session()
S.headers.update({"User-Agent": "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) "
                                "AppleWebKit/537.36 Chrome/120 Safari/537.36"})

# ============================ 一、分段下载 ============================
_lock = threading.Lock()
_stat = {"n": 0}


def head_len(url, tries=3):
    for i in range(tries):
        try:
            r = S.head(url, timeout=60, allow_redirects=True)
            if r.status_code == 200:
                return int(r.headers.get("Content-Length", 0)), r.headers.get("Accept-Ranges")
        except Exception:
            time.sleep(2 + i * 2)
    return 0, None


def fetch_part(url, a, b, path):
    """下载 [a,b] 字节写入 path，支持续传。"""
    while True:
        have = os.path.getsize(path) if os.path.exists(path) else 0
        if have >= b - a + 1:
            return True
        try:
            r = S.get(url, headers={"Range": f"bytes={a+have}-{b}"},
                      timeout=120, stream=True)
            if r.status_code not in (200, 206):
                time.sleep(3); continue
            with open(path, "ab") as f:
                for ck in r.iter_content(1 << 18):
                    f.write(ck)
                    with _lock:
                        _stat["n"] += len(ck)
        except Exception:
            time.sleep(3)


def download(nw=8):
    if os.path.exists(TARGET) and os.path.getsize(TARGET) > 1_000_000 and gzip_ok(TARGET):
        print(f"  [跳过] 已存在 {TARGET} {os.path.getsize(TARGET)/1048576:.1f} MB", flush=True)
        return True
    L, ar = head_len(EXPR_URL)
    if not L:
        print("  取不到文件长度", flush=True); return False
    print(f"  HiSeqV2.gz {L/1048576:.1f} MB, Accept-Ranges={ar}", flush=True)
    os.makedirs(PARTS, exist_ok=True)

    if ar != "bytes":
        # 不支持 Range -> 单线程流式
        print("  服务端不支持断点，改单线程流式下载", flush=True)
        t0 = time.time()
        r = S.get(EXPR_URL, stream=True, timeout=180)
        with open(TARGET, "wb") as f:
            for ck in r.iter_content(1 << 20):
                f.write(ck)
                with _lock:
                    _stat["n"] += len(ck)
        print(f"  完成 {os.path.getsize(TARGET)/1048576:.1f} MB / {time.time()-t0:.0f}s", flush=True)
        return gzip_ok(TARGET)

    seg = np.linspace(0, L, num=nw + 1, dtype=int)
    jobs = [(seg[i], seg[i+1] - 1, os.path.join(PARTS, f"seg_{i}.part"))
            for i in range(nw)]
    t0 = time.time()
    last = 0
    ths = []
    for a, b, p in jobs:
        t = threading.Thread(target=fetch_part, args=(EXPR_URL, a, b, p), daemon=True)
        t.start(); ths.append(t)
    while any(t.is_alive() for t in ths):
        time.sleep(3)
        n = _stat["n"] / 1048576
        dt = time.time() - t0
        if n > last:
            print(f"    {n:.1f}/{L/1048576:.1f} MB  平均 {n/dt:.2f} MB/s", flush=True)
            last = n
    for t in ths:
        t.join()

    total = sum(os.path.getsize(p) for _, _, p in jobs if os.path.exists(p))
    if total < L:
        print(f"  不完整 {total/1048576:.1f}/{L/1048576:.1f} MB", flush=True)
        return False
    with open(TARGET, "wb") as fo:
        for a, b, p in jobs:
            with open(p, "rb") as fi:
                while True:
                    blk = fi.read(1 << 20)
                    if not blk: break
                    fo.write(blk)
    ok = gzip_ok(TARGET)
    print(f"  拼装 {'OK' if ok else 'FAIL'} {os.path.getsize(TARGET)/1048576:.1f} MB "
          f"/ {time.time()-t0:.0f}s（平均 {L/1048576/(time.time()-t0):.2f} MB/s）", flush=True)
    if ok:
        for _, _, p in jobs:
            try: os.remove(p)
            except OSError: pass
    return ok


def gzip_ok(path):
    try:
        with gzip.open(path, "rb") as f:
            while f.read(1 << 20):
                pass
        return True
    except Exception:
        return False


# ============================ 二、拉临床 ============================
CB = "https://www.cbioportal.org/api"
STUDY = "coadread_tcga_pan_can_atlas_2018"
CLIN_ATTRS = ["OS_MONTHS", "OS_STATUS", "AGE", "SEX",
              "AJCC_PATHOLOGIC_TUMOR_STAGE", "DFS_MONTHS", "DFS_STATUS"]


def get_json(url, params=None, tries=3):
    for t in range(tries):
        try:
            r = S.get(url, params=params, timeout=120)
            if r.status_code == 200:
                return r.json()
        except Exception:
            pass
        time.sleep(2 + 3 * t)
    return None


def fetch_clinical():
    """临床走 GET /studies/{id}/clinical-data。

    【踩坑】原打算 POST /clinical-data/fetch 带 patientIds 批处理，但这个 API
    版本不认该字段名（试 identifiers 也报 "error in the JSON format"），且会
    静默返回 0 条而不报错——上一次就是因此把空 DataFrame 直接送进了 pivot 才崩的。
    GET 端点不仅可用，还直接给出明文 patientId（省掉 base64 解码 uniquePatientKey
    这一步），一次请求拿全部 ATTRIBUTE，比批处理还快。
    """
    print("\n② cBioPortal 拉临床 ...", flush=True)
    rows, off, PAGE = [], 0, 1000000
    while True:
        d = get_json(f"{CB}/studies/{STUDY}/clinical-data",
                     {"clinicalDataType": "PATIENT", "projection": "SUMMARY",
                      "pageSize": PAGE, "pageNumber": off})
        if not d:
            break
        rows.extend(d)
        if len(d) < PAGE:
            break
        off += 1
    print(f"   GET 拿到 {len(rows)} 条临床记录", flush=True)
    cdf = pd.DataFrame(rows)
    if cdf.empty:
        raise RuntimeError("临床记录为空，检查 cBioPortal 接口")
    piv = cdf.pivot_table(index="patientId", columns="clinicalAttributeId",
                          values="value", aggfunc="first")
    piv.columns.name = None
    print(f"   临床覆盖 {piv.shape[0]} 人；可用属性 {len(piv.columns)} 个", flush=True)
    for a in CLIN_ATTRS:
        print(f"     {a:<32} 有值 {int(piv[a].notna().sum()) if a in piv.columns else 0}", flush=True)

    # 注：分期的解析函数见下方 clin 组装之后（要等 clin 建好才能读 AJCC 原值）
    # Xena 的列名是 12 位病人 barcode（TCGA-3L-AA1B），cBio 600 多个 patientId 也是同一口径
    clin = pd.DataFrame({"pid12": piv.index})
    clin.columns = ["pid12"]
    clin = clin.join(piv.reset_index(drop=True))
    for c in CLIN_ATTRS:
        if c not in clin.columns:
            clin[c] = np.nan
    clin["time"] = pd.to_numeric(clin["OS_MONTHS"], errors="coerce")
    clin["status"] = pd.to_numeric(
        clin["OS_STATUS"].astype(str).str.extract(r"^([01])", expand=False), errors="coerce")
    clin["age"] = pd.to_numeric(clin["AGE"], errors="coerce")
    clin["sex"] = clin["SEX"].astype(str).str.strip().str.lower()
    # 【踩坑】这里原本写死 {"Stage IIA": 2, ...}，而 cBioPortal GET 返回的是全大写
    # "STAGE IIA"，一个大小写之差让 580 个分期全部静默变成 NA。stage 是 delta_C
    # （签名相对临床变量的增量）的核心协变量，缺了它 delta_C 会被系统性高估，
    # 且与那 7 个 GEO 队列不可比。改为大小写无关的正则 + 罗马数字查表。
    ROMAN = {"0": 0, "I": 1, "IA": 1, "IB": 1, "II": 2, "IIA": 2, "IIB": 2,
             "IIC": 2, "III": 3, "IIIA": 3, "IIIB": 3, "IIIC": 3,
             "IV": 4, "IVA": 4, "IVB": 4, "IVC": 4}
    _RE_STAGE = re.compile(r"^\s*(?:stage\s*)?([0IV]{1,3}[ABC]?)\s*$", re.I)

    def parse_stage(v):
        if v is None or (isinstance(v, float) and np.isnan(v)):
            return np.nan
        m = _RE_STAGE.match(str(v).strip())
        return ROMAN.get(m.group(1).upper(), np.nan) if m else np.nan

    raw_stage = clin["AJCC_PATHOLOGIC_TUMOR_STAGE"]
    print("   分期原始取值:", dict(raw_stage.value_counts().head(6)), flush=True)
    clin["stage"] = raw_stage.map(parse_stage)
    unmapped = raw_stage.notna() & clin["stage"].isna()
    if unmapped.any():
        print(f"   ⚠ 未被识别的分期 {int(unmapped.sum())} 例: "
              f"{dict(raw_stage[unmapped].value_counts().head(5))}", flush=True)
    clin["dfs_time"] = pd.to_numeric(clin["DFS_MONTHS"], errors="coerce")
    clin["dfs_status"] = pd.to_numeric(
        clin["DFS_STATUS"].astype(str).str.extract(r"^([01])", expand=False), errors="coerce")
    clin["site"] = clin["pid12"].str.extract(r"^TCGA-([A-Z0-9]{2})-", expand=False)
    print(f"   有 OS 随访 {int(clin.time.notna().sum())}/{len(clin)}，"
          f"事件 {int((clin.status==1).sum())}", flush=True)
    print(f"   分期 {int(clin.stage.notna().sum())}，"
          f"性别 {int(clin.sex.isin(['male','female']).sum())}，"
          f"年龄 {int(clin.age.notna().sum())}，"
          f"DFS 事件 {int((clin.dfs_status==1).sum())}", flush=True)
    return clin


# ============================ 三、读表达并对齐 ============================
def bar12(x):
    """TCGA-3L-AA1B-01A / TCGA-3L-AA1B-01 -> TCGA-3L-AA1B（12 位）"""
    m = re.match(r"^(TCGA-[A-Z0-9]{2}-[A-Z0-9]{4})", str(x))
    return m.group(1) if m else None


def load_expr():
    print("\n③ 读 HiSeqV2 表达矩阵 ...", flush=True)
    with gzip.open(TARGET, "rt") as f:
        mat = pd.read_csv(f, sep="\t", index_col=0)
    print(f"   原始 {mat.shape[0]} 基因 x {mat.shape[1]} 样本", flush=True)

    # 只看原发灶：barcode 位点 15-16 == 01
    def is_tumor(s):
        p = str(s).split("-")
        return len(p) >= 4 and p[3][:2] == "01"
    keep = [c for c in mat.columns if is_tumor(c)]
    print(f"   原发灶样本（01）: {len(keep)} / {mat.shape[1]}", flush=True)
    if len(keep) < mat.shape[1]:
        print(f"   剔除的例（多为正常组织 -11）: "
              f"{[c for c in mat.columns if c not in keep][:5]}", flush=True)
    mat = mat[keep]

    # 重复 Sample \barcode ->12 位去重（同一个病人既有 01A 又有 01B 时取均值容易失真，取第一个）
    key = pd.Series(keep).map(bar12).values
    mat.columns = key
    dup = pd.Index(mat.columns).duplicated().sum()
    if dup:
        print(f"   去重前的重复病人 {dup} 个 -> 按列去重（保留首列）", flush=True)
        mat = mat.loc[:, ~mat.columns.duplicated()]

    vmax = float(np.nanmax(mat.values))
    print(f"   数值范围 {np.nanmin(mat.values):.2f} ~ {vmax:.2f}", flush=True)
    if vmax > 60:
        print("   -> 判定线性尺度，做 log2(x+1)", flush=True)
        mat = np.log2(mat + 1)
    else:
        print("   -> 判定已是对数尺度（RSEM log2），不变换", flush=True)
    mat = mat.astype("float32")
    return mat


if __name__ == "__main__":
    t0 = time.time()
    print("① 下载表达矩阵 ...", flush=True)
    if not download(int(sys.argv[1]) if len(sys.argv) > 1 else 8):
        sys.exit(1)

    mat = load_expr()
    clin = fetch_clinical()

    print("\n④ 样本对齐（表达 Xena <-> 临床 cBio）...", flush=True)
    # mat.columns 在 load_expr 里已规约成 12 位，clin 的主键也是 12 位，口径一致
    m12 = set(mat.columns)
    c12 = set(clin["pid12"].dropna())
    common = sorted(m12 & c12)
    print(f"   表达 {len(m12)} 人，临床 {len(c12)} 人，交集 {len(common)} 人", flush=True)
    print(f"   有表达但缺临床 {len(m12 - c12)} 人；有临床但缺表达 {len(c12 - m12)} 人",
          flush=True)
    mat = mat[common]

    cl = clin.set_index("pid12").loc[common].reset_index(drop=True)
    cl.insert(0, "sample", common)
    keep_cols = ["sample", "time", "status", "stage", "age", "sex", "site",
                 "dfs_time", "dfs_status"]
    cl = cl[keep_cols]

    ok = cl.time.notna() & cl.status.notna()
    mat = mat[[s for s in mat.columns]]
    X = mat.loc[:, ok.values]
    cl = cl[ok.values].reset_index(drop=True)
    print(f"   有 OS 结局的最终样本 {cl.shape[0]}，事件 {int((cl.status==1).sum())}",
          flush=True)
    print(f"   分期 {int(cl.stage.notna().sum())}，"
          f"性别 {int(cl.sex.isin(['male','female']).sum())}，"
          f"年龄 {int(cl.age.notna().sum())}，"
          f"DFS 事件 {int((cl.dfs_status==1).sum())}", flush=True)
    print(f"   随访中位数 {cl.time.median():.1f} 月，最长 {cl.time.max():.0f} 月",
          flush=True)

    print("\n⑤ 落盘 ...", flush=True)
    # 过滤低表达基因：与 13_tcga_to_rds.R 的口径一致（>20% 样本有表达）
    expressed = (X > 0).mean(axis=1)
    X = X[expressed >= 0.2]
    print(f"   过滤低表达后 {X.shape[0]} 基因 x {X.shape[1]} 样本", flush=True)
    X.T.to_csv(os.path.join(OUT, "TCGA_expr.tsv.gz"), sep="\t",
               compression="gzip", float_format="%.4f")
    cl.to_csv(os.path.join(OUT, "TCGA_clin.csv"), index=False, encoding="utf-8-sig")
    print(f"   {OUT}/TCGA_expr.tsv.gz + TCGA_clin.csv")
    print(f"   总耗时 {(time.time()-t0)/60:.1f} 分钟", flush=True)
