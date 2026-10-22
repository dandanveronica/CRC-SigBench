# -*- coding: utf-8 -*-
"""Figure 3（增量判别力 vs 两个参照）— 312/2,496 口径生成脚本 · 2026-10-18

底本：08_核验底稿_勿上传/重跑产物_肿瘤only_20261017/figure3_redo.py（10-17 重绘，视觉布局）
对照：CRC-SigBench/code/10_figures/37_make_figures.py 的 fig2()（10-14 权威措辞）

本轮修复（H-19 及 Figure 3 遗留项）：
  1. 数据源换成交付源 S1_2496 次评估逐条结果.csv（312 x 8 = 2,496）
  2. A 面板标题数字改为由数据动态生成（旧脚本硬编码 "2,560"，重跑会退回旧口径）
  3. B 面板标题恢复 10-14 原文 "Increment above routine clinical variables"
     （10-17 redo 省略了 routine；正文 idx26 有 "routine clinical variables" 出处）
  4. D 面板队列名恢复 "TCGA COADREAD"（10-14 用 replace('_',' ')；
     与 Table 1 一致，旧 redo 保留下划线与表格不一致）
  5. C 面板**保留 within-signature mean**（77.6% / 98.7%）——该口径为 10-14 既定、
     且 recheck D47 锁定；与 S1 §3.17 的 median 惯例差异改由 §3.17 例外清单披露
  6. 按内容裁边（_trim）+ 同步裁矢量 PDF（_trim_pdf），消除 H-18d 类 PNG/PDF 不一致

数据源：results/submission_tables/S1_2496次评估逐条结果.csv
输出：results/figures/Figure_3.png + Figure_3.pdf
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import rcParams
from scipy import stats

# Elsevier artwork instructions: 只允许 Arial/Helvetica、Courier、Symbol、Times
rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['Arial', 'Helvetica']
rcParams['font.size'] = 8.5
# 面板标题字号 8.5。B 标题恢复 10-14 原文 "routine" 后为四面板最长，
# 单行时右边界顶到画布（OCR 实测 x1=0.999，末字被裁），故拆为两行；其余三面板单行。
rcParams['axes.titlesize'] = 8.5
rcParams['axes.labelsize'] = 8
rcParams['xtick.labelsize'] = 7.5
rcParams['ytick.labelsize'] = 7.5
rcParams['axes.linewidth'] = 0.8
rcParams['xtick.major.width'] = 0.8
rcParams['ytick.major.width'] = 0.8
rcParams['pdf.fonttype'] = 42
rcParams['ps.fonttype'] = 42

BLUE1 = '#2E6DA4'
BLUE2 = '#5B9BD5'
TEAL1 = '#178C7F'
TEAL2 = '#76D7C4'
RED = '#C0392B'
GREY = '#5B6770'

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
SRC = os.path.join(SUP, 'S1_2496次评估逐条结果.csv')
OUT_PNG = os.path.join(OUT, 'Figure_3.png')
OUT_PDF = os.path.join(OUT, 'Figure_3.pdf')
DPI = 1000


def _content_box(path, pad=24):
    """返回按实际内容（非白像素）计算的裁剪框。"""
    from PIL import Image
    im = Image.open(path)
    bg = Image.new('RGB', im.size, (255, 255, 255))
    bg.paste(im, mask=im.split()[3] if im.mode == 'RGBA' else None)
    a = np.asarray(bg, int)
    ys, xs = np.where(a.sum(2) < 750)
    crop = (max(0, xs.min() - pad), max(0, ys.min() - pad),
            min(a.shape[1], xs.max() + 1 + pad), min(a.shape[0], ys.max() + 1 + pad))
    return crop


def _trim(path, dpi, pad=24):
    """按内容裁 PNG（savefig 的白边不依赖 bbox_inches='tight'，避免与 constrained_layout 冲突）。"""
    from PIL import Image
    crop = _content_box(path, pad)
    im = Image.open(path)
    bg = Image.new('RGB', im.size, (255, 255, 255))
    bg.paste(im, mask=im.split()[3] if im.mode == 'RGBA' else None)
    bg.crop(crop).save(path, dpi=(dpi, dpi))
    return crop


def _trim_pdf(pdf_path, crop, size):
    """把 PDF 裁到与 PNG 相同的边界（H-18d：PIL 左上 y 向下 → PDF 左下 y 向上）。"""
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
    pc = pd.read_csv(SRC, encoding='utf-8-sig')
    # d_rand 已在源文件中；d_clin 与 delta_C 完全相等（10-17 专项核查：相关系数 1.000000）
    pc['d_clin'] = pc['delta_C']
    n_eval = len(pc)
    assert n_eval == 2496, f'期望 2,496 行，实际 {n_eval}'
    assert pc.sig_id.nunique() == 312 and pc.cohort.nunique() == 8

    S = pc.groupby('sig_id').agg(drand=('d_rand', 'mean'),
                                 dclin=('d_clin', 'mean')).reset_index()

    fig = plt.figure(figsize=(7.48, 6.0), constrained_layout=True)
    # w_pad 由 0.35 收到 0.25：给最长的 B 面板标题让出右侧空间（标题从 axes 左边界
    # 起，可溢出到 figure padding；trim 时另有 pad=24 px 补回边距）
    fig.set_constrained_layout_pads(w_pad=0.25, h_pad=0.35, hspace=0.28, wspace=0.28)
    gs = fig.add_gridspec(2, 2)

    # ---------------- A ----------------
    ax1 = fig.add_subplot(gs[0, 0])
    vals = pc.d_rand.dropna()
    ax1.hist(vals, bins=60, color=BLUE2, edgecolor='white', linewidth=0.3)
    med = vals.median()
    ax1.axvline(0, color='black', linewidth=0.9)
    ax1.axvline(med, color=RED, linestyle='--', linewidth=1.2)
    pos = (vals > 0).mean() * 100
    p = stats.wilcoxon(vals).pvalue
    txt = f'median {med:+.4f}\n{pos:.1f}% positive\np = {p:.1e}'.replace('e-0', 'e-')
    ax1.text(0.03, 0.97, txt, transform=ax1.transAxes, ha='left', va='top',
             fontsize=8.0, bbox=dict(boxstyle='round,pad=0.3', facecolor='white',
                                     edgecolor='grey', linewidth=0.5))
    ax1.set_xlabel('ΔC_random  (signature C − matched random-gene-set median)')
    ax1.set_ylabel('Number of evaluations')
    # 数字由数据动态生成，避免硬编码退回旧口径
    ax1.set_title(f'A  {n_eval:,} evaluations vs random gene sets', loc='left', fontweight='bold')

    # ---------------- B ----------------
    ax2 = fig.add_subplot(gs[0, 1])
    vals2 = pc.d_clin.dropna()
    ax2.hist(vals2, bins=60, color=TEAL1, edgecolor='white', linewidth=0.3)
    med2 = vals2.median()
    ax2.axvline(0, color='black', linewidth=0.9)
    ax2.axvline(med2, color=RED, linestyle='--', linewidth=1.2)
    pos2 = (vals2 > 0).mean() * 100
    p2 = stats.wilcoxon(vals2).pvalue
    txt2 = f'median {med2:+.4f}\n{pos2:.1f}% positive\np = {p2:.1e}'.replace('e-0', 'e-')
    ax2.text(0.03, 0.97, txt2, transform=ax2.transAxes, ha='left', va='top',
             fontsize=8.0, bbox=dict(boxstyle='round,pad=0.3', facecolor='white',
                                     edgecolor='grey', linewidth=0.5))
    ax2.set_xlabel('ΔC_clinical  (clinical + signature, minus clinical alone)')
    ax2.set_ylabel('Number of evaluations')
    # 恢复 10-14 原文 "routine"（正文 idx26 "against routine clinical variables"）
    ax2.set_title('B  Increment above routine\nclinical variables', loc='left', fontweight='bold')

    # ---------------- C ----------------
    ax3 = fig.add_subplot(gs[1, 0])
    cats = ['random:\nrecord', 'random:\nsignature', 'clinical:\nrecord', 'clinical:\nsignature']
    heights = [100 * (pc.d_rand > 0).mean(), 100 * (S.drand > 0).mean(),
               100 * (pc.d_clin > 0).mean(), 100 * (S.dclin > 0).mean()]
    colors = [BLUE1, BLUE2, TEAL1, TEAL2]
    bars = ax3.bar(range(4), heights, color=colors, edgecolor='white', width=0.7)
    ax3.axhline(50, color=RED, linestyle=':', linewidth=0.9)
    ax3.text(3.35, 51, 'chance', fontsize=7.5, color=RED, va='bottom')
    for b, h in zip(bars, heights):
        ax3.text(b.get_x() + b.get_width() / 2, h + 2, f'{h:.1f}%',
                 ha='center', va='bottom', fontsize=8.0)
    ax3.set_xticks(range(4))
    ax3.set_xticklabels(cats)
    ax3.set_ylim(0, 110)
    ax3.set_ylabel('% of evaluations / signatures positive')
    ax3.set_title('C  Evaluation noise resolves into direction once averaged', loc='left',
                  fontweight='bold')

    # ---------------- D ----------------
    ax4 = fig.add_subplot(gs[1, 1])
    g = pc.groupby('cohort').agg(m=('d_rand', 'median'),
                                 q1=('d_rand', lambda x: x.quantile(0.25)),
                                 q3=('d_rand', lambda x: x.quantile(0.75))).reset_index()
    g = g.sort_values('m', ascending=True).reset_index(drop=True)
    y = np.arange(len(g))
    assert (g.m >= 0).all(), '存在负中位数队列，D 标题/图注需改写'
    ax4.axvline(0, color='black', linewidth=0.8)
    for i, row in g.iterrows():
        col = TEAL1 if row.cohort == 'TCGA_COADREAD' else ('#E67E22' if row.m < 0 else GREY)
        ax4.errorbar(row.m, i, xerr=[[row.m - row.q1], [row.q3 - row.m]], fmt='o',
                     markersize=6, color=col, ecolor=col, elinewidth=2, capsize=0)
    ax4.set_yticks(y)
    # 恢复 10-14：去下划线，与 Table 1 "TCGA COADREAD" 一致
    ax4.set_yticklabels([c.replace('_', ' ') for c in g.cohort])
    ax4.set_xlabel('Median ΔC_random with interquartile range')
    ax4.set_title('D  Per-cohort median ΔC_random', loc='left', fontweight='bold')

    # ---------------- 保存 ----------------
    tmp_pdf = OUT_PDF
    fig.savefig(OUT_PNG, dpi=DPI, facecolor='white')
    fig.savefig(tmp_pdf, dpi=DPI, facecolor='white')
    size = plt.gcf().canvas.get_width_height()[::-1] if False else None
    from PIL import Image
    size = Image.open(OUT_PNG).size
    crop = _trim(OUT_PNG, DPI)
    _trim_pdf(tmp_pdf, crop, size)
    plt.close(fig)

    # ---------------- 复算自检 ----------------
    print('=== Figure 3 复算自检 ===')
    print(f' n evaluations      = {n_eval:,}')
    print(f' A median ΔC_random = {med:+.4f}  ({pos:.1f}% positive, p = {p:.2e})')
    print(f' B median ΔC_clinic = {med2:+.4f}  ({pos2:.1f}% positive, p = {p2:.2e})')
    print(f' C record/signature = {heights[0]:.1f}% / {heights[1]:.1f}% '
          f'(random), {heights[2]:.1f}% / {heights[3]:.1f}% (clinical)')
    print(f' D per-cohort median min = {g.m.min():+.4f} ({g.cohort[g.m.idxmin()]})，'
          f'负中位数队列数 = {(g.m < 0).sum()}')
    # 尾部截断量（图注未披露，需确认是否可接受）
    lo_a, hi_a = vals.min(), vals.max()
    print(f' A 数据范围 [{lo_a:.4f}, {hi_a:.4f}]，轴外（< -0.05 或 > 0.20）= '
          f'{int(((vals < -0.05) | (vals > 0.20)).sum())}')
    lo_b, hi_b = vals2.min(), vals2.max()
    print(f' B 数据范围 [{lo_b:.4f}, {hi_b:.4f}]，轴外（< -0.05 或 > 0.15）= '
          f'{int(((vals2 < -0.05) | (vals2 > 0.15)).sum())}')
    return fig


if __name__ == '__main__':
    build()
    from PIL import Image
    with Image.open(OUT_PNG) as im:
        w, h = im.size
        dpi = im.info.get('dpi', (DPI, DPI))[0]
    print()
    print('=== 交付规格 ===')
    print(f' PNG {w} x {h} px @ {dpi:.0f} dpi = {w / dpi * 25.4:.1f} x {h / dpi * 25.4:.1f} mm')
