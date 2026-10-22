# -*- coding: utf-8 -*-
"""Figure 5（Seven years without improvement）— 312/2,496 口径生成脚本 · 2026-10-18

底本：CRC-SigBench/code/10_figures/37_make_figures.py 的 fig4()（10-14 版，措辞与布局）
尺寸：figsize (7.2, 5.36) + dpi 1400 —— 成图 176.58 x 130.08 mm（JCE 限 <=190 mm）

数据源（全部实算，不硬编码）：
  A  S1_2496 次评估 + S5_312 签名清单  → record/signature 级时间趋势 Spearman
  B  S5_312                              → 每年 usable-formula / external-validation 率
  C  S1_2496                             → 7 GEO 均值 vs TCGA RNA-seq 跨平台 Spearman
  D  S1_2496                             → 3/312 在可用区（SIG095/SIG186/SIG018）

修复/补齐：
  1. 产出矢量 PDF（此前 Figure 5 只有 PNG）
  2. 按内容裁边（_trim）+ 同步裁 PDF（_trim_pdf），PNG 与 PDF 逐毫米一致
  3. 措辞恢复 10-14 旧脚本标题，数字按 312 口径更新
  4. D 面板标注 3 个可用区签名（caption 已点名 SIG018 + Table 4 两个）
  5. 新增渲染器级几何自检块

# H-42/H-43 备注（2026-10-18 第五轮核验）：
#   docstring 首行由 "Eight years" 更正为 "Seven years"（2020-2026 为 7 个日历年）。
#   _trim_pdf 在 bbox_inches='tight' 之后实际不产生额外裁切（content_box 加 pad 后
#   恰好覆盖全图），故 PDF 保留 matplotlib 的 MediaBox；实测 MediaBox == CropBox，
#   尺寸 176.58 x 130.08 mm，与 PNG 一致，无内容被裁。
#   ⚠ 勿用 pymupdf 的 page.rect / page.cropbox 判断 Box 偏移——其返回值被内部归一化
#     （会把 y0 报成 0），据此误判过"顶部被裁 2.1 mm"。判断须用 pypdf 读原始对象，
#     或用 CoreGraphics(CGPDFPageGetBoxRect) 交叉验证。
输出：results/figures/Figure_5.png + Figure_5.pdf
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import rcParams
from scipy import stats

rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['Arial', 'Helvetica']
rcParams['font.size'] = 8.5
rcParams['axes.linewidth'] = 0.8
rcParams['xtick.major.width'] = 0.8
rcParams['ytick.major.width'] = 0.8
rcParams['pdf.fonttype'] = 42
rcParams['ps.fonttype'] = 42

GREY, BLUE, RED, ORANGE, TEAL = '#5B6770', '#2E6DA4', '#C0392B', '#D68910', '#178C7F'

# --- repository-relative paths (v1.2.0) -------------------------------------
# 原脚本硬编码本机投稿包绝对路径；收编进仓库后改为由脚本自身位置推导，
# clone 后无需改任何路径即可复现投稿图。
#   SUP : 投稿用补充材料 CSV（从投稿包 supplementary 目录原样搬入）
#   OUT : 图产物（投稿包 figure 目录的同名定稿图）
_HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(_HERE, os.pardir, os.pardir))
SUP = os.path.join(REPO, 'results', 'submission_tables')
OUT = os.path.join(REPO, 'results', 'figures')
os.makedirs(OUT, exist_ok=True)
# ---------------------------------------------------------------------------
OUT_PNG = os.path.join(OUT, 'Figure_5.png')
OUT_PDF = os.path.join(OUT, 'Figure_5.pdf')
DPI = 1400


def _content_box(path, pad=24):
    from PIL import Image
    im = Image.open(path)
    bg = Image.new('RGB', im.size, (255, 255, 255))
    bg.paste(im, mask=im.split()[3] if im.mode == 'RGBA' else None)
    a = np.asarray(bg, int)
    ys, xs = np.where(a.sum(2) < 750)
    return (max(0, xs.min() - pad), max(0, ys.min() - pad),
            min(a.shape[1], xs.max() + 1 + pad), min(a.shape[0], ys.max() + 1 + pad))


def _trim(path, dpi, pad=24):
    from PIL import Image
    crop = _content_box(path, pad)
    im = Image.open(path)
    bg = Image.new('RGB', im.size, (255, 255, 255))
    bg.paste(im, mask=im.split()[3] if im.mode == 'RGBA' else None)
    bg.crop(crop).save(path, dpi=(dpi, dpi))
    return crop


def _trim_pdf(pdf_path, crop, size):
    from pypdf import PdfReader, PdfWriter
    x0, y0, x1, y1 = crop
    W0, H0 = size
    r = PdfReader(pdf_path)
    page = r.pages[0]
    ox, oy = float(page.mediabox.left), float(page.mediabox.bottom)
    pw, ph = float(page.mediabox.width), float(page.mediabox.height)
    left = ox + (x0 / W0) * pw
    right = ox + (x1 / W0) * pw
    bottom = oy + ph - (y1 / H0) * ph
    top = oy + ph - (y0 / H0) * ph
    for box in ('mediabox', 'cropbox'):
        getattr(page, box).lower_left = (left, bottom)
        getattr(page, box).upper_right = (right, top)
    w = PdfWriter()
    w.add_page(page)
    w.write(pdf_path)


def build():
    # ---------------- 数据 ----------------
    s1 = pd.read_csv(f'{SUP}/S1_2496次评估逐条结果.csv', encoding='utf-8-sig')
    s5 = pd.read_csv(f'{SUP}/S5_312签名清单.csv', encoding='utf-8-sig')
    yrmap = dict(zip(s5.sig_id, s5.publication_year))
    s1['year'] = s1.sig_id.map(yrmap)

    fig, axs = plt.subplots(2, 2, figsize=(7.2, 5.36))

    # ---------------- A ----------------
    ax = axs[0, 0]
    rec = s1.dropna(subset=['year', 'd_rand'])
    r_rec, p_rec = stats.spearmanr(rec.year, rec.d_rand)
    sig = s1.groupby('sig_id').agg(year=('year', 'first'), d_rand=('d_rand', 'mean')).reset_index()
    r_sig, p_sig = stats.spearmanr(sig.year, sig.d_rand)
    years = sorted(rec['year'].dropna().unique())
    box_dat = [rec[rec.year == y].d_rand.values for y in years]
    bp = ax.boxplot(box_dat, positions=years, widths=.55, patch_artist=True,
                    medianprops=dict(color=RED, lw=1.4), showfliers=False)
    for b in bp['boxes']:
        b.set(facecolor='#C8DCEC', edgecolor=GREY, lw=.8)
    ax.axhline(0, color='k', lw=.9)
    ax.text(.02, .94,
            f'record level: ρ = {r_rec:+.3f}, p = {p_rec:.3f}\n'
            f'(signature level: ρ = {r_sig:+.3f}, p = {p_sig:.3f})',
            transform=ax.transAxes, va='top', fontsize=7.4,
            bbox=dict(fc='white', ec=GREY, lw=.6, pad=3))
    ax.set_xlabel('Publication year')
    ax.set_ylabel('ΔC_random')
    # H-26 修复：抬高上界，使 2023–2026 上须帽（0.131–0.138）落在信息框（白底）下方
    ax.set_ylim(-0.075, 0.21)
    ax.set_title('A  No trend in true incremental value', loc='left',
                 fontsize=9, fontweight='bold')

    # ---------------- B ----------------
    ax = axs[0, 1]
    byv, bye = {}, {}
    for _, r in s5.iterrows():
        y = int(r.publication_year)
        byv.setdefault(y, []).append(int(r.risk_formula_published))
        bye.setdefault(y, []).append(int(r.external_validation_declared))
    vr = sorted((y, sum(v) / len(v) * 100, len(v)) for y, v in byv.items() if len(v) >= 5)
    er = sorted((y, sum(v) / len(v) * 100, len(v)) for y, v in bye.items() if len(v) >= 5)
    ax.plot([y for y, _, _ in er], [v for _, v, _ in er],
            'o-', color=BLUE, lw=1.6, ms=4.5, label='claims external validation')
    ax.plot([y for y, _, _ in vr], [v for _, v, _ in vr],
            's--', color=RED, lw=1.6, ms=4.5, label='publishes usable coefficients')
    r_trend, p_trend = stats.spearmanr(s5.publication_year, s5.external_validation_declared)
    # H-33 修复：ρ 的样本是 312 个签名（per-signature 二分变量），与面板上 6 个年度聚合点
    # （若据其计算 ρ 恰为 +1.000）不是同一集合——必须在图面标出 n = 312，否则读者会误读。
    # H-35 修复：面板仅含 ≥5 篇的年份（2021–2026），故无 2020 点；须在图面披露阈值。
    # H-37 修复：3 行加长后信息框右端与右下图例互叠（B 面板过窄放不下两个白底框）——
    # 图例两行并入信息框（分色），不再使用独立 legend。
    ax.text(.02, .06,
            '● claims external validation\n'
            '■ publishes usable coefficients\n'
            f'trend (per signature, n = {len(s5)};\n'
            'years with n ≥ 5):\n'
            f'ρ = {r_trend:+.3f}, p = {p_trend:.3f}\n'
            f'verifiability {vr[0][1]:.1f}% → {vr[-1][1]:.1f}%',
            transform=ax.transAxes, fontsize=7.4, va='bottom', linespacing=1.35,
            bbox=dict(fc='white', ec=GREY, lw=.6, pad=3))
    ax.set_ylim(0, 105)
    ax.set_xlabel('Publication year')
    ax.set_ylabel('% of signatures')
    # H-27 原修复（图例移右下）已被 H-37 取代：图例并入信息框，独立 legend 删除
    ax.set_title('B  Reported standards rise, verifiability falls', loc='left',
                 fontsize=9, fontweight='bold')

    # ---------------- C ----------------
    ax = axs[1, 0]
    geo = s1[s1.cohort != 'TCGA_COADREAD'].groupby('sig_id')['d_rand'].mean()
    rna = s1[s1.cohort == 'TCGA_COADREAD'].groupby('sig_id')['d_rand'].mean()
    common = sorted(set(geo.index) & set(rna.index))
    xg = [geo[k] for k in common]
    yrna = [rna[k] for k in common]
    r_c, p_c = stats.spearmanr(xg, yrna)
    ax.scatter(xg, yrna, s=13, c=BLUE, alpha=.45, edgecolors='none')
    m, b = np.polyfit(xg, yrna, 1)
    ax.plot([min(xg), max(xg)], [m * min(xg) + b, m * max(xg) + b],
            color=RED, lw=1.0)
    ax.axvline(0, color='k', lw=.7)
    ax.axhline(0, color='k', lw=.7)
    # H-32 修复：红线是 OLS 最小二乘拟合（Pearson r = +0.310），与所报 Spearman ρ = +0.272
    # 并非同一统计量，图注未解释该线——读者无法得知红线含义。图面补注，行数不变避免新遮挡。
    # H-38 修复：不用缩写 "OLS"（主稿 Abbreviations 无此条目，图面不应出现未定义缩写），
    # 改为完整写法 "least-squares fit"。
    ax.text(.02, .95,
            f'ρ = {r_c:+.3f}\np = {p_c:.1e}\n'
            f'n = {len(common)}; red line, least-squares fit',
            transform=ax.transAxes, va='top', fontsize=7.6,
            bbox=dict(fc='white', ec=GREY, lw=.6, pad=3))
    # H-30 修复：抬高上界，信息框覆盖的数据区（y>0.17）无散点，不再压 SIG013/069/081/251/304/SIG244
    # H-38 修复：信息框第 3 行改为完整写法 "least-squares fit" 后变宽，上界再抬至 0.28
    ax.set_ylim(-0.054, 0.28)
    ax.set_xlabel('Mean ΔC_random in 7 GEO microarray cohorts')
    ax.set_ylabel('ΔC_random in TCGA RNA-seq')
    ax.set_title('C  Cross-platform reproducibility is weak', loc='left',
                 fontsize=9, fontweight='bold')

    # ---------------- D ----------------
    ax = axs[1, 1]
    means = s1.groupby('sig_id').agg(C_mean=('C', 'mean'),
                                     calib_mean=('calib_slope', 'mean')).reset_index()
    hit = means[(means.C_mean >= 0.60) & (means.calib_mean >= 0.70)].copy()
    ax.scatter(means.C_mean, means.calib_mean, s=13, c=GREY, alpha=.45,
               edgecolors='none', zorder=3)
    ax.scatter(hit.C_mean, hit.calib_mean, s=20, c=TEAL, alpha=.85,
               edgecolors='none', zorder=5)
    ax.axvline(0.60, color=ORANGE, ls='--', lw=1.0)
    ax.axhline(0.70, color=ORANGE, ls='--', lw=1.0)
    # H-31 修复：信息框拆 3 行短行，右缘缩至 C≈0.56，避开 0.60 竖虚线与全部高 slope 散点
    ax.text(.025, .93,
            'usable region:\nC ≥ 0.60, slope ≥ 0.70\n'
            f'{len(hit)} of {len(means)} signatures',
            transform=ax.transAxes, va='top', fontsize=7.4,
            bbox=dict(fc='white', ec=RED, lw=.6, pad=3))
    # 标注可用区签名；手动偏移避免重叠
    offsets = {
        'SIG095': (-6, -7),   # H-36 修复：原 (5, 5) 使文字贴面板右缘（含 pad 即越界）
        'SIG186': (7, -2),    # H-29 修复：原 (5,-12) 使文字被 0.70 虚线穿过；H-36 微调避 SIG095
        'SIG018': (6, -16),   # H-36 修复：原 (-40, 2) 压 3 个散点，改置点右下后仅叠 2 个淡点
    }
    for _, row in hit.iterrows():
        dx, dy = offsets.get(row.sig_id, (5, 5))
        ax.annotate(row.sig_id, (row.C_mean, row.calib_mean),
                    fontsize=7.0, color=RED, zorder=6,   # H-28：6.6 → 7.0 pt（最低字号下限）
                    xytext=(dx, dy), textcoords='offset points',
                    ha='left' if dx > 0 else 'right')
    ax.set_xlabel('Mean C-index across 8 cohorts')
    ax.set_ylabel('Mean calibration slope (ideal = 1.0)')
    ax.set_title('D  The clinically usable quadrant is nearly empty', loc='left',
                 fontsize=9, fontweight='bold')

    fig.tight_layout(rect=[0, 0, 1, 1])

    # ---- 几何自检：数字/注记/图元不得越出所属面板或画布 ----
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    canvas = FigureCanvasAgg(fig)
    canvas.draw()
    rend = canvas.get_renderer()
    fb = fig.bbox
    print('=== 几何自检（越界图元）===')
    n_bad = 0

    def _flag(name, what, s, detail=''):
        nonlocal n_bad
        n_bad += 1
        print(f'  ⚠ {name} {what} {s}  {detail}')

    for name, ax in (('A', axs[0, 0]), ('B', axs[0, 1]),
                     ('C', axs[1, 0]), ('D', axs[1, 1])):
        ab = ax.get_window_extent(renderer=rend)
        for t in ax.texts:
            tb = t.get_window_extent(renderer=rend)
            if tb.width == 0 and tb.height == 0:
                continue
            out_ax = (tb.x0 < ab.x0 - 1 or tb.x1 > ab.x1 + 1 or
                      tb.y0 < ab.y0 - 1 or tb.y1 > ab.y1 + 1)
            out_fig = (tb.x0 < 0 or tb.x1 > fb.width or tb.y0 < 0 or tb.y1 > fb.height)
            if out_ax or out_fig:
                _flag(name, '文字', f'"{t.get_text()[:34]}"',
                      ('越出面板' if out_ax else '') +
                      ('/' if out_ax and out_fig else '') +
                      ('越出画布' if out_fig else '') +
                      f'  (text x {tb.x0:.0f}-{tb.x1:.0f} y {tb.y0:.0f}-{tb.y1:.0f}, '
                      f'axes x {ab.x0:.0f}-{ab.x1:.0f} y {ab.y0:.0f}-{ab.y1:.0f})')
        xl, yl = ax.get_xlim(), ax.get_ylim()
        from matplotlib.collections import LineCollection
        for coll in ax.collections:
            if isinstance(coll, LineCollection):
                continue
            offs = coll.get_offsets()
            if offs is None or len(offs) == 0:
                continue
            for x, y in offs[:, :2]:
                if x < xl[0] or x > xl[1] or y < yl[0] or y > yl[1]:
                    _flag(name, '数据点', f'({x:.3f}, {y:.3f})',
                          f'越出坐标范围 x[{xl[0]:.2f},{xl[1]:.2f}] y[{yl[0]:.2f},{yl[1]:.2f}]')
        for lab, t, axis in (('title', ax.title, 'x'),
                             ('xlabel', ax.xaxis.label, 'x'),
                             ('ylabel', ax.yaxis.label, 'y')):
            tb = t.get_window_extent(renderer=rend)
            if tb.width == 0 and tb.height == 0:
                continue
            if (tb.x0 < -1 or tb.x1 > fb.width + 1 or tb.y0 < -1 or tb.y1 > fb.height + 1):
                _flag(name, lab, f'"{t.get_text()[:30]}"', '越出画布')
            elif axis == 'x' and (tb.x0 < ab.x0 - 1 or tb.x1 > ab.x1 + 1):
                _flag(name, lab, f'"{t.get_text()[:30]}"', '横向越出面板')
            elif axis == 'y' and (tb.y0 < ab.y0 - 1 or tb.y1 > ab.y1 + 1):
                _flag(name, lab, f'"{t.get_text()[:30]}"', '纵向越出面板')
        lg = ax.get_legend()
        if lg is not None:
            lb = lg.get_window_extent(renderer=rend)
            if (lb.x0 < -1 or lb.x1 > fb.width + 1 or lb.y0 < -1 or lb.y1 > fb.height + 1):
                _flag(name, '图例', '', f'越出画布 (x {lb.x0:.0f}-{lb.x1:.0f}, '
                                        f'y {lb.y0:.0f}-{lb.y1:.0f})')
    print(f'  → 越界图元合计 {n_bad} 个' + ('（须修）' if n_bad else '（无）'))

    fig.savefig(OUT_PNG, dpi=DPI, bbox_inches='tight', facecolor='white')
    fig.savefig(OUT_PDF, dpi=DPI, bbox_inches='tight', facecolor='white')
    from PIL import Image
    size = Image.open(OUT_PNG).size
    crop = _trim(OUT_PNG, DPI)
    _trim_pdf(OUT_PDF, crop, size)
    plt.close(fig)

    print('=== Figure 5 复算自检 ===')
    print(f' A record ρ={r_rec:+.3f} p={p_rec:.3f}; signature ρ={r_sig:+.3f} p={p_sig:.3f}')
    print(f' B external-validation trend ρ={r_trend:+.3f} p={p_trend:.3f}; '
          f'verifiability {vr[0][1]:.1f}% → {vr[-1][1]:.1f}%')
    print(f' C cross-platform ρ={r_c:+.3f} p={p_c:.1e} n={len(common)}')
    print(f' D usable region {len(hit)}/{len(means)}: {list(hit.sig_id)}')


if __name__ == '__main__':
    build()
    from PIL import Image
    with Image.open(OUT_PNG) as im:
        w, h = im.size
        dpi = im.info.get('dpi', (DPI, DPI))[0]
    print()
    print('=== 交付规格 ===')
    print(f' PNG {w} x {h} px @ {dpi:.0f} dpi = {w / dpi * 25.4:.1f} x {h / dpi * 25.4:.1f} mm')
