#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
07_down.py —— GEO series matrix 多线程断点续传下载器

为什么不用 curl + xargs：
  NCBI 在被并发打狠时会返回 XML 错误页（HTTP 200），
  shell 版 worker 把错误页当成数据追加写进分段，越下越脏。
  这里每收到一块就校验前 5 字节是否为 '<?xml'，是就丢弃并退避重试。

用法：
  python 07_down.py              # 下载全部缺失队列
  python 07_down.py --workers 8  # 指定线程数
"""
import os, sys, time, threading, urllib.request, urllib.error

RAW = os.path.join(os.path.dirname(os.path.abspath(__file__)), "cohorts_raw")
PARTS = os.path.join(RAW, ".parts")
ALL = ["GSE39582", "GSE17536", "GSE72970", "GSE38832", "GSE87211",
       "GSE17537", "GSE29621"]

def url_of(g):
    return ("https://ftp.ncbi.nlm.nih.gov/geo/series/%snnn/%s/matrix/"
            "%s_series_matrix.txt.gz" % (g[:5], g, g))

def head_len(url, tries=3):
    for i in range(tries):
        try:
            req = urllib.request.Request(url, method="HEAD")
            with urllib.request.urlopen(req, timeout=30) as r:
                return int(r.headers.get("Content-Length", 0))
        except Exception:
            time.sleep(2 + i * 2)
    return 0

def gzip_ok(path):
    """能解压到最后一块且无异常即认为完整。"""
    import gzip
    try:
        with gzip.open(path, "rb") as f:
            while f.read(1 << 20):
                pass
        return True
    except Exception:
        return False

class Fetcher(threading.Thread):
    def __init__(self, g, start, end, stop, stats):
        super().__init__(daemon=True)
        self.g, self.s, self.e = g, start, end
        self.stop = stop
        self.stats = stats
        self.path = os.path.join(PARTS, "%s.%d_%d.part" % (g, start, end))
        self.url = url_of(g)

    def run(self):
        lock = self.stats["lock"]
        while self.s <= self.e and not self.stop.is_set():
            have = os.path.getsize(self.path) if os.path.exists(self.path) else 0
            # 已写盘的数据如果起始部分是错误页，清空
            if have > 0 and self.is_xml_head(self.path):
                try: os.remove(self.path)
                except OSError: pass
                have = 0
            want_s = self.s + have
            if want_s > self.e:
                break
            data = self.get(want_s, self.e)
            if data is None:            # 失败 -> 退避后重试
                if self.stop.is_set(): break
                time.sleep(3)
                continue
            if data[:5] == b"<?xml":    # NCBI 限流错误页
                time.sleep(6)
                continue
            try:
                with open(self.path, "ab") as f:
                    f.write(data)
                with lock:
                    self.stats["bytes"] += len(data)
            except OSError:
                time.sleep(2)

    def is_xml_head(self, p):
        try:
            with open(p, "rb") as f:
                return f.read(5) == b"<?xml"
        except OSError:
            return False

    def get(self, a, b):
        hdr = {"Range": "bytes=%d-%d" % (a, b)}
        try:
            req = urllib.request.Request(self.url, headers=hdr)
            with urllib.request.urlopen(req, timeout=45) as r:
                return r.read()
        except Exception:
            return None

def coverage(g):
    """返回 [(start, size), ...] 已覆盖区间"""
    out = []
    pref = g + "."
    if not os.path.isdir(PARTS):
        return out
    for fn in os.listdir(PARTS):
        if not fn.startswith(pref) or not fn.endswith(".part"):
            continue
        mid = fn[len(pref):-5]           # start_end
        st = mid.split("_")[0]
        if not st.isdigit():
            continue
        sz = os.path.getsize(os.path.join(PARTS, fn))
        if sz > 0:
            out.append((int(st), sz))
    out.sort()
    return out

def gaps_of(cov, L):
    pos, gaps = 0, []
    for st, sz in cov:
        if st > pos:
            gaps.append((pos, st - 1))
        e = st + sz
        if e > pos:
            pos = e
    if pos < L:
        gaps.append((pos, L - 1))
    return gaps

def assemble(g, L):
    cov = coverage(g)
    pos, ok = 0, True
    for st, sz in cov:
        if st > pos:
            ok = False
        e = st + sz
        if e > pos:
            pos = e
    if pos < L:
        ok = False
    if not ok:
        return False, "已覆盖 %d/%d MB" % (pos >> 20, L >> 20)
    out = os.path.join(RAW, "%s_series_matrix.txt.gz" % g)
    tmp = out + ".asm"
    with open(tmp, "wb") as fo:
        pos = 0
        for st, sz in cov:
            files = [f for f in os.listdir(PARTS)
                     if f.startswith("%s.%d_" % (g, st)) and f.endswith(".part")]
            if not files:
                continue
            p = os.path.join(PARTS, files[0])
            with open(p, "rb") as fi:
                if st < pos:
                    fi.seek(pos - st)
                while True:
                    blk = fi.read(1 << 20)
                    if not blk:
                        break
                    fo.write(blk)
            pos = st + sz
    if gzip_ok(tmp):
        os.replace(tmp, out)
        # 清理该队列分段
        for f in os.listdir(PARTS):
            if f.startswith(g + ".") and f.endswith(".part"):
                try: os.remove(os.path.join(PARTS, f))
                except OSError: pass
        return True, "完成 %.1f MB" % (os.path.getsize(out) / 1048576)
    try: os.remove(tmp)
    except OSError: pass
    return False, "gzip 校验失败"

def main():
    nw = 8
    if "--workers" in sys.argv:
        nw = int(sys.argv[sys.argv.index("--workers") + 1])
    os.makedirs(PARTS, exist_ok=True)
    todo = []
    print("=== 检查现有文件 ===")
    for g in ALL:
        out = os.path.join(RAW, "%s_series_matrix.txt.gz" % g)
        if os.path.exists(out) and os.path.getsize(out) > 1000000 and gzip_ok(out):
            print("  [跳过] %s 已完整 %.1f MB" % (g, os.path.getsize(out) / 1048576))
            continue
        L = head_len(url_of(g))
        if not L:
            print("  [警告] %s 取不到大小" % g)
            continue
        cov = coverage(g)
        gaps = gaps_of(cov, L)
        need = sum(b - a + 1 for a, b in gaps)
        print("  [待下] %s %d/%d MB" % (g, need >> 20, L >> 20))
        if need:
            todo.append((g, L, gaps))
    if not todo:
        print("\n全部完成。")
        return

    stop = threading.Event()
    stats = {"bytes": 0, "lock": threading.Lock()}
    for g, L, gaps in todo:
        jobs = []
        for a, b in gaps:
            sz = b - a + 1
            k = max(1, (sz + (4 << 20) - 1) // (4 << 20))   # 每段 ≤4 MB
            cs = (sz + k - 1) // k
            for i in range(k):
                s = a + i * cs
                e = min(s + cs - 1, b)
                if s > b:
                    break
                jobs.append(Fetcher(g, s, e, stop, stats))
        print("\n>>> %s 切 %d 段，线程 %d" % (g, len(jobs), nw))
        t0 = time.time()
        for i in range(0, len(jobs), nw):
            batch = jobs[i:i + nw]
            for t in batch: t.start()
            for t in batch: t.join()
        d = time.time() - t0
        print("    %s 本轮 %.1f MB / %.0f 秒" % (g, stats["bytes"] / 1048576, d))

    print("\n=== 拼装与校验 ===")
    for g, L, _ in todo:
        ok, msg = assemble(g, L)
        print("  %-11s %s %s" % (g, "✓" if ok else "✗", msg))

if __name__ == "__main__":
    main()
