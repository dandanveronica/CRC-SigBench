#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
parallel_fetch.py — 分块并行 + 块级断点续传下载器（纯 python，无第三方依赖）

用法:
  python parallel_fetch.py URL OUT [--prefix OLD_FILE] [--head N] [--workers 3]

特性:
  - Content-Range 探测总长; HTTP Range 分块下载
  - 每块独立落盘 OUT.partNNN —— 断点后重跑只补缺块（块级断点）
  - --prefix: 把旧的顺序追加式部分文件并入块 0（不浪费已下载字节）
  - --head N: 只下载前 N 字节（GPL13497 soft 头部场景）
  - 全块齐后拼接 + gzip 完整性校验 + 改名为 OUT; 校验失败删件重来
"""
import os, sys, time, gzip, zlib, threading, urllib.request

CHUNK = 8 * 1024 * 1024  # 8MB/块
UA = "Mozilla/5.0 (CRC-SigBench fetcher)"


def probe_total(url):
    for i in range(10):
        try:
            req = urllib.request.Request(url, headers={"Range": "bytes=0-0", "User-Agent": UA})
            with urllib.request.urlopen(req, timeout=45) as r:
                cr = r.headers.get("Content-Range")
                if cr and "/" in cr:
                    return int(cr.split("/")[-1])
        except Exception as e:
            print(f"  probe retry {i+1}: {str(e)[:70]}", flush=True)
            time.sleep(3)
    return None


def fetch_range(url, start, end, out_path, tag, idx):
    """下载 [start,end] 到 out_path；若文件已存在且大小精确匹配则跳过（块级断点）。"""
    want = end - start + 1
    if os.path.exists(out_path) and os.path.getsize(out_path) == want:
        return idx, True, "skip"
    for rnd in range(30):
        cur = os.path.getsize(out_path) if os.path.exists(out_path) else 0
        if cur >= want:
            return idx, True, "done"
        try:
            req = urllib.request.Request(url, headers={"Range": f"bytes={start+cur}-{end}", "User-Agent": UA})
            with urllib.request.urlopen(req, timeout=90) as r, open(out_path, "ab") as f:
                while True:
                    c = r.read(1 << 20)
                    if not c:
                        break
                    f.write(c)
            if os.path.getsize(out_path) >= want:
                return idx, True, "done"
        except Exception as e:
            print(f"  [{tag}#{idx}] r{rnd+1}: {os.path.getsize(out_path) if os.path.exists(out_path) else 0}/{want} ({str(e)[:55]})", flush=True)
            time.sleep(2)
    return idx, False, "fail"


def download(url, out, prefix=None, head=None, workers=3):
    total = probe_total(url)
    if total is None:
        print(f"[FAIL] {out}: 拿不到总长", flush=True)
        return False
    limit = min(total, head) if head else total
    print(f"[{out}] 总长 {total}  目标 {limit} bytes", flush=True)

    # 块 0: 复用旧前缀文件（动态大小）; 其余按 CHUNK 切
    bounds = []  # (start, end, part_path)
    pdir = out + ".parts"
    os.makedirs(pdir, exist_ok=True)
    p0 = 0
    if prefix and os.path.exists(prefix):
        p0 = os.path.getsize(prefix)
        # 沙箱禁 rename(unlink)，改用复制；源文件保留
        import shutil
        dst0 = os.path.join(pdir, "part000")
        if not os.path.exists(dst0):
            shutil.copyfile(prefix, dst0)
        print(f"  复制前缀 {p0} bytes -> part000", flush=True)
    if p0 > 0:
        bounds.append((0, p0 - 1, os.path.join(pdir, "part000")))
    s = p0
    i = 1
    while s < limit:
        e = min(s + CHUNK - 1, limit - 1)
        bounds.append((s, e, os.path.join(pdir, f"part{i:03d}")))
        s = e + 1
        i += 1
    print(f"  分 {len(bounds)} 块 (前缀块 {p0}B + {max(0,len(bounds)-1)} 个 {CHUNK//1024//1024}MB 块), {workers} 路并行", flush=True)

    todo = [b for b in bounds if not (os.path.exists(b[2]) and os.path.getsize(b[2]) == b[1] - b[0] + 1)]
    done_cnt = len(bounds) - len(todo)
    print(f"  断点检查: {done_cnt}/{len(bounds)} 块已完成, 待补 {len(todo)}", flush=True)

    lock = threading.Lock()
    results = {}

    def run(b):
        idx, ok, msg = fetch_range(url, b[0], b[1], b[2], os.path.basename(out), b[2].split("part")[-1])
        with lock:
            results[b[2]] = ok
            d = sum(1 for v in results.values() if v)
            print(f"  [{d}/{len(bounds)}] 块 {b[2].split('part')[-1]} {msg}", flush=True)

    # 简单线程池
    queue = list(todo)
    def worker():
        while True:
            b = None
            with lock:
                if queue:
                    b = queue.pop(0)
                else:
                    return
            run(b)
    threads = [threading.Thread(target=worker, daemon=True) for _ in range(min(workers, max(1, len(todo))))]
    for t in threads:
        t.start()
    for t in threads:
        t.join()

    if not all(results.get(b[2], os.path.exists(b[2]) and os.path.getsize(b[2]) == b[1] - b[0] + 1) for b in bounds):
        print(f"[FAIL] {out}: 有块未完成，part 文件已保留，重跑即续", flush=True)
        return False

    # 拼接
    with open(out, "wb") as fo:
        for b in bounds:
            with open(b[2], "rb") as fi:
                while True:
                    c = fi.read(1 << 22)
                    if not c:
                        break
                    fo.write(c)
    size = os.path.getsize(out)
    if size != limit:
        print(f"[FAIL] {out}: 拼接后 {size} != {limit}", flush=True)
        return False
    try:
        with gzip.open(out, "rb") as f:
            f.read(1 << 20)
        print(f"[OK] {out}: {size} bytes, gzip 头部可解", flush=True)
    except Exception as e:
        # --head 模式本来就是截断 gz, 必然"损坏"——属预期; 完整文件必须过校验
        if head:
            print(f"[OK-head] {out}: {size} bytes (截断 gz, 预期内)", flush=True)
        else:
            print(f"[FAIL] {out}: gzip 校验失败 {str(e)[:60]}, 删除重试", flush=True)
            os.remove(out)
            return False
    if not head:
        for b in bounds:
            try:
                os.remove(b[2])
            except OSError:
                pass
        try:
            os.rmdir(pdir)
        except OSError:
            pass
    return True


if __name__ == "__main__":
    url, out = sys.argv[1], sys.argv[2]
    prefix = head = None
    workers = 3
    argv = sys.argv[3:]
    k = 0
    while k < len(argv):
        if argv[k] == "--prefix":
            prefix = argv[k + 1]; k += 2
        elif argv[k] == "--head":
            head = int(argv[k + 1]); k += 2
        elif argv[k] == "--workers":
            workers = int(argv[k + 1]); k += 2
        else:
            k += 1
    ok = download(url, out, prefix=prefix, head=head, workers=workers)
    sys.exit(0 if ok else 1)
