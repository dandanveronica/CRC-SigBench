# -*- coding: utf-8 -*-
"""Figure 4（三种制造显著性的机制）— 312/2,496 口径生成脚本 · 2026-10-18

底本：CRC-SigBench/code/10_figures/37_make_figures.py 的 fig3()（10-14 版，布局与措辞）
尺寸：figsize (7.2, 5.36) + dpi 1400 —— 沿用 10-15 第四轮终版（成品 176 mm、最小印刷字号 7.54 pt）

面板数据源（全部实算，不硬编码；仅 B 面板因逐行源文件缺失而用 JCE_S1 Table A2 定稿值）：
  A  S10_截点搜索_已发表2184与随机3280.csv    → 12.4%→37.8% / 7.8%→30.7%，n=2,184 / 3,280
  B  A2 定稿值 + S11 实算随机基线            → in-sample 18.5% vs 15.0%（+3.5 pp, p=0.037）
                                              10-fold CV 25.2% vs 16.4%（+8.8 pp）
  C  S7_复现率配对明细.csv                   → n=100, median gap +0.161, 99% above, p=4.3e-18
  D  S12b（训练 apparent）+ S13（7 队列验证） → 0.902 → 0.597, drop 0.305；literature gap 同源 C

本轮修复/补齐：
  1. 产出矢量 PDF（此前 Figure 4 只有 PNG，无矢量版）
  2. 按内容裁边（_trim）+ 同步裁 PDF（_trim_pdf），PNG 与 PDF 逐毫米一致
  3. D 面板标题沿用 10-17 版 "Compliant signature: a larger drop"
     （10-14 原文为 "the same drop"，但 312 口径下 drop 0.305 > literature gap 0.161，
       "larger" 才与数据相符；10-17 的改法正确，予以保留）
  4. 全部数字由源数据渲染，杜绝 10-15 §3.2 那类「数据源更新了、图上数字没跟着变」

输出：results/figures/Figure_4.png + Figure_4.pdf
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib.transforms import blended_transform_factory
from scipy import stats

# Elsevier artwork instructions: 只允许 Arial/Helvetica、Courier、Symbol、Times
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
OUT_PNG = os.path.join(OUT, 'Figure_4.png')
OUT_PDF = os.path.join(OUT, 'Figure_4.pdf')
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
    """把 PDF 裁到与 PNG 相同边界（PIL 左上 y 向下 → PDF 左下 y 向上）。"""
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
    s10 = pd.read_csv(f'{SUP}/S10_截点搜索_已发表2184与随机3280.csv', encoding='utf-8-sig')
    pub = s10[s10.arm == 'published_signature']
    rnd = s10[s10.arm == 'size_matched_random']
    assert len(pub) == 2184 and len(rnd) == 3280, (len(pub), len(rnd))

    s11 = pd.read_csv(f'{SUP}/S11_NRI与IDI_随机基因集基线700.csv', encoding='utf-8-sig')
    rnd_in = 100 * (s11.NRI_p < .05).mean()          # 15.0%（105/700）
    # B 面板其余三项来自 JCE_S1 Table A2 定稿（逐行源文件未随 312 重跑保留）
    A2 = {'sig_in': 18.5, 'rnd_in': 15.0, 'p_in': 0.037,
          'sig_cv': 25.2, 'rnd_cv': 16.4}
    assert abs(rnd_in - A2['rnd_in']) < 0.05, f'S11 实算 {rnd_in:.2f}% 与 A2 {A2["rnd_in"]}% 不符'

    s7 = pd.read_csv(f'{SUP}/S7_复现率配对明细.csv', encoding='utf-8-sig')
    cl = s7.dropna(subset=['claim_auc5', '重算AUC5中位'])
    gap = (cl.claim_auc5 - cl['重算AUC5中位']).median()
    above = 100 * ((cl.claim_auc5 - cl['重算AUC5中位']) > 0).mean()
    p_c = stats.wilcoxon(cl.claim_auc5 - cl['重算AUC5中位']).pvalue

    s12b = pd.read_csv(f'{SUP}/S12b_合规签名_训练集apparent表现.csv', encoding='utf-8-sig')
    s13 = pd.read_csv(f'{SUP}/S13_合规签名_队列间验证7队列.csv', encoding='utf-8-sig')
    c_train = s12b.C_train.iloc[0]
    v_med = s13.C.median()
    drop = c_train - v_med

    fig, axs = plt.subplots(2, 2, figsize=(7.2, 5.36))

    # ---------------- A ----------------
    ax = axs[0, 0]
    pre = [100 * (pub.p_median < .05).mean(), 100 * (rnd.p_median < .05).mean()]
    best = [100 * (pub.p_best < .05).mean(), 100 * (rnd.p_best < .05).mean()]
    x = np.arange(2)
    w = .36
    ax.bar(x - w / 2, pre, w, color=GREY, edgecolor='k', lw=.6, label='pre-specified cut-point')
    ax.bar(x + w / 2, best, w, color=RED, edgecolor='k', lw=.6, label='best cut-point by search')
    ax.axhline(5, color=BLUE, ls='--', lw=1.1)
    # H-22b：旧版 x=1.48 是数据坐标，但 xlim 上界只有 1.446（bar 自动留边），
    # ha='right' 会把右边界顶到面板外 4 px 被裁 → 改用 axes 横向比例定位，恒在框内
    trans = blended_transform_factory(ax.transAxes, ax.transData)
    ax.text(0.995, 6, 'nominal 5%', color=BLUE, ha='right', va='bottom', fontsize=7,
            transform=trans)
    for i, v in enumerate(pre):
        ax.text(i - w / 2, v + 1.4, f'{v:.1f}%', ha='center', fontsize=7.4, zorder=6)
    for i, v in enumerate(best):
        ax.text(i + w / 2, v + 1.4, f'{v:.1f}%', ha='center', fontsize=7.4, zorder=6)
    ax.set_xticks(x)
    ax.set_xticklabels([f'Published\nsignatures (n={len(pub):,})',
                        f'Random\ngene sets (n={len(rnd):,})'], fontsize=7.4)
    ax.set_ylabel('% reaching p < 0.05')
    ax.set_ylim(0, 50)
    ax.legend(fontsize=7, frameon=True, loc='upper center', ncol=2, framealpha=.95)
    ax.set_title('A  Shared grey-to-red gain is the search itself', loc='left',
                 fontsize=9, fontweight='bold')

    # ---------------- B ----------------
    ax = axs[0, 1]
    sets = [('In-sample', A2['sig_in'], A2['rnd_in']), ('10-fold CV', A2['sig_cv'], A2['rnd_cv'])]
    x = np.arange(2)
    ax.bar(x - w / 2, [s[1] for s in sets], w, color=TEAL, edgecolor='k', lw=.6,
           label='published signatures')
    ax.bar(x + w / 2, [s[2] for s in sets], w, color='#C9CBCD', edgecolor='k', lw=.6,
           label='random gene sets')
    ax.axhline(5, color=BLUE, ls='--', lw=1.1)
    # H-22b：同 A 面板，数据坐标 1.48 越出 xlim 上界 → 改 axes 横向比例定位
    trans = blended_transform_factory(ax.transAxes, ax.transData)
    ax.text(0.995, 6, 'nominal 5%', color=BLUE, ha='right', va='bottom', fontsize=7,
            transform=trans)
    for i, (_, a_, b_) in enumerate(sets):
        ax.text(i - w / 2, a_ + .8, f'{a_:.1f}%', ha='center', fontsize=7.4)
        ax.text(i + w / 2, b_ + .8, f'{b_:.1f}%', ha='center', fontsize=7.4)
    # 注记上移并 va='bottom'（10-15 §3.3：原固定 y 会与柱值标签叠印）
    g_in = A2['sig_in'] - A2['rnd_in']
    g_cv = A2['sig_cv'] - A2['rnd_cv']
    ax.annotate(f'gap {g_in:.1f} pp\n(p = {A2["p_in"]:.3f})', xy=(-w / 2, 22.4),
                fontsize=7, color=RED, ha='left', va='bottom')
    ax.annotate(f'gap {g_cv:.1f} pp', xy=(1 - w / 2, 27.4),
                fontsize=7, color=RED, ha='left', va='bottom')
    ax.set_xticks(x)
    ax.set_xticklabels([s[0] for s in sets], fontsize=7.4)
    ax.set_ylabel('% of NRI tests reaching p < 0.05')
    ax.set_ylim(0, 36)
    ax.legend(fontsize=7, loc='upper left')
    ax.set_title(f'B  In-sample NRI is anti-conservative ({A2["rnd_in"]/5:.1f}x nominal)',
                 loc='left', fontsize=9, fontweight='bold')

    # ---------------- C ----------------
    ax = axs[1, 0]
    ax.scatter(cl['重算AUC5中位'], cl.claim_auc5, s=16, c=RED, alpha=.55, edgecolors='none')
    lo, hi = 0.3, 1.0
    ax.plot([lo, hi], [lo, hi], color='k', lw=.9)
    ax.plot([lo, hi - gap], [lo + gap, hi], color=RED, lw=1.3, ls='--')
    ax.text(.33, .95, f'median gap +{gap:.3f}\n{above:.0f}% above the diagonal\n'
                      f'n = {len(cl)},  p = {p_c:.1e}',
            fontsize=7.6, va='top', bbox=dict(fc='white', ec=GREY, lw=.6, pad=3))
    ax.set_xlabel('Recomputed 5-year AUC (this study)')
    ax.set_ylabel('Published claimed 5-year AUC')
    ax.set_xlim(lo, hi)
    ax.set_ylim(lo, hi)
    ax.set_title('C  Claimed versus recomputed', loc='left', fontsize=9, fontweight='bold')

    # ---------------- D ----------------
    ax = axs[1, 1]
    v = s13.sort_values('C')
    nm = [({'TCGA_COADREAD': 'TCGA'}.get(c, c)) for c in v.cohort]
    yp = np.arange(len(v))
    ax.scatter([c_train], [len(v)], s=70, marker='D', c=RED, zorder=5)
    for k, c in enumerate(v.C):
        ax.scatter(c, k, s=34, c=BLUE, zorder=4)
    ax.plot([c_train, v.C.median()], [len(v), len(v) / 2], color=RED, lw=1.7, ls='--')
    ax.hlines(len(v) / 2, v.C.min(), c_train, ls=':', lw=.9, color=GREY)
    # 文字随 training C 走（旧版固定 x=0.62 是为 training=0.801 设计的，
    # 现 training=0.902，固定值会让注记与点分离；
    # H-22c： baseline y 须下调，避免两行文本顶出面板顶部
    # H-23： text 必须避开 red diamond（zorder 提升至 marker 之上，并整体左移）
    ax.text(c_train - .05, len(v) - .25, f'training\n{c_train:.3f}',
            fontsize=7.4, ha='right', va='center', color=RED, zorder=6)
    # validation median 注记右移并提高 zorder，避免被蓝色 dot 遮挡
    ax.text(v.C.median() + .05, len(v) / 2 + .25, f'validation median {v.C.median():.3f}',
            fontsize=7.4, ha='left', va='center', color=BLUE, zorder=6)
    # 三行注记整体左移，避免右端越出面板（旧版 x=0.815 在新 xlim 下右端会冲出 76 px）
    ax.text(0.70, 1.0, f'drop {drop:.3f}\nliterature gap {gap:.3f}\ncompliant drop is larger',
            fontsize=7.4, color=RED, ha='left', va='bottom', zorder=6)
    ax.set_yticks(list(yp) + [len(v)])
    ax.set_yticklabels(list(nm) + ['training\nGSE39582'], fontsize=7.2)
    # xlim 上界必须容得下 training C：10-14 时 training=0.801 故上界 0.86 够用；
    # 312 口径 training=0.902，沿用 0.86 会把 training 点画到面板外被裁掉（H-22）
    ax.set_xlabel('C-index')
    ax.set_xlim(0.50, 0.95)
    ax.grid(axis='x', ls=':', lw=.4, color='#BBBBBB')
    ax.set_title('D  Compliant signature: a larger drop', loc='left', fontsize=9, fontweight='bold')

    fig.tight_layout(rect=[0, 0, 1, 1])

    # ---- 几何自检：数字/注记/图元不得越出所属面板（axes）或画布（figure）----
    # 判据分层（避免误报）：
    #   · ax.texts（数据锚定的数字与注记）→ 必须完整落在 axes 内
    #   · collections（散点/柱）        → 必须落在 xlim/ylim 内，否则被 axes 裁掉
    #   · title / xlabel / ylabel       → 本就在 axes 之外，只做「越出画布」+ 横向/纵向对齐检查
    #   · legend                        → 必须完整落在画布内
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

    for name, ax in (('A', axs[0, 0]), ('B', axs[0, 1]), ('C', axs[1, 0]), ('D', axs[1, 1])):
        ab = ax.get_window_extent(renderer=rend)
        # 1) 数据锚定文字：必须完全在 axes 内
        for t in ax.texts:
            tb = t.get_window_extent(renderer=rend)
            if tb.width == 0 and tb.height == 0:
                continue
            out_ax = (tb.x0 < ab.x0 - 1 or tb.x1 > ab.x1 + 1 or
                      tb.y0 < ab.y0 - 1 or tb.y1 > ab.y1 + 1)
            out_fig = (tb.x0 < 0 or tb.x1 > fb.width or tb.y0 < 0 or tb.y1 > fb.height)
            if out_ax or out_fig:
                _flag(name, '文字', f'"{t.get_text()[:34]}"',
                      ('越出面板' if out_ax else '') + ('/' if out_ax and out_fig else '') +
                      ('越出画布' if out_fig else '') +
                      f'  (text x {tb.x0:.0f}-{tb.x1:.0f} y {tb.y0:.0f}-{tb.y1:.0f}, '
                      f'axes x {ab.x0:.0f}-{ab.x1:.0f} y {ab.y0:.0f}-{ab.y1:.0f})')
        # 2) 数据点：必须落在坐标范围内（否则被 axes 裁掉）
        xl, yl = ax.get_xlim(), ax.get_ylim()
        from matplotlib.collections import LineCollection
        for coll in ax.collections:
            # LineCollection（如 hlines）的 offsets 不是散点数据，跳过
            if isinstance(coll, LineCollection):
                continue
            offs = coll.get_offsets()
            if offs is None or len(offs) == 0:
                continue
            for x, y in offs[:, :2]:
                if x < xl[0] or x > xl[1] or y < yl[0] or y > yl[1]:
                    _flag(name, '数据点', f'({x:.3f}, {y:.3f})',
                          f'越出坐标范围 x[{xl[0]:.2f},{xl[1]:.2f}] y[{yl[0]:.2f},{yl[1]:.2f}]')
        # 3) 标题 / 轴标签：不得越出画布；标题与 xlabel 横向不得越过面板，ylabel 纵向不得越过面板
        for lab, t, axis in ((('title'), ax.title, 'x'),
                             (('xlabel'), ax.xaxis.label, 'x'),
                             (('ylabel'), ax.yaxis.label, 'y')):
            tb = t.get_window_extent(renderer=rend)
            if tb.width == 0 and tb.height == 0:
                continue
            if (tb.x0 < -1 or tb.x1 > fb.width + 1 or tb.y0 < -1 or tb.y1 > fb.height + 1):
                _flag(name, lab, f'"{t.get_text()[:30]}"', '越出画布')
            elif axis == 'x' and (tb.x0 < ab.x0 - 1 or tb.x1 > ab.x1 + 1):
                _flag(name, lab, f'"{t.get_text()[:30]}"', '横向越出面板')
            elif axis == 'y' and (tb.y0 < ab.y0 - 1 or tb.y1 > ab.y1 + 1):
                _flag(name, lab, f'"{t.get_text()[:30]}"', '纵向越出面板')
        # 4) 图例：必须完整落在画布内（允许比 axes 宽，因为顶置图例常故意跨面板宽度）
        lg = ax.get_legend()
        if lg is not None:
            lb = lg.get_window_extent(renderer=rend)
            if (lb.x0 < -1 or lb.x1 > fb.width + 1 or lb.y0 < -1 or lb.y1 > fb.height + 1):
                _flag(name, '图例', '', f'越出画布 (x {lb.x0:.0f}-{lb.x1:.0f}, '
                                        f'y {lb.y0:.0f}-{lb.y1:.0f}, canvas {fb.width:.0f}x{fb.height:.0f})')
    print(f'  → 越界图元合计 {n_bad} 个' + ('（须修）' if n_bad else '（无）'))

    fig.savefig(OUT_PNG, dpi=DPI, bbox_inches='tight', facecolor='white')
    fig.savefig(OUT_PDF, dpi=DPI, bbox_inches='tight', facecolor='white')
    from PIL import Image
    size = Image.open(OUT_PNG).size
    crop = _trim(OUT_PNG, DPI)
    _trim_pdf(OUT_PDF, crop, size)
    plt.close(fig)

    print('=== Figure 4 复算自检 ===')
    print(f' A published {pre[0]:.1f}% -> {best[0]:.1f}%  (n={len(pub):,})')
    print(f' A random    {pre[1]:.1f}% -> {best[1]:.1f}%  (n={len(rnd):,})')
    print(f' B in-sample {A2["sig_in"]}% vs {A2["rnd_in"]}%  gap {g_in:.1f} pp (p={A2["p_in"]})')
    print(f' B 10-fold   {A2["sig_cv"]}% vs {A2["rnd_cv"]}%  gap {g_cv:.1f} pp')
    print(f' C n={len(cl)}  median gap {gap:+.4f}  {above:.0f}% above  p={p_c:.2e}')
    print(f' D training {c_train:.4f} -> validation median {v_med:.4f}  drop {drop:.4f} '
          f'(literature gap {gap:.3f})')


if __name__ == '__main__':
    build()
    from PIL import Image
    with Image.open(OUT_PNG) as im:
        w, h = im.size
        dpi = im.info.get('dpi', (DPI, DPI))[0]
    print()
    print('=== 交付规格 ===')
    print(f' PNG {w} x {h} px @ {dpi:.0f} dpi = {w / dpi * 25.4:.1f} x {h / dpi * 25.4:.1f} mm')
