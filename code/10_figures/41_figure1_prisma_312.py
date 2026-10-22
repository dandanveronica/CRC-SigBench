# -*- coding: utf-8 -*-
"""Figure 1 — PRISMA 2020 flow diagram, adapted to signature-level identification.

Renumbered 2026-10-14: this diagram is now Figure 1 (it is cited first, in 3.1);
the study-design figure became Figure 2.

All counts are taken verbatim from the auditable screening register
(sigminer/output/PRISMA_counts.txt) and from the frozen library files:
    signature_library_main.csv        n = 338  (primary-analysis library)
    signature_library_clean.csv       n = 320  (registration-time snapshot, v1.0.x)
    signature_library_312_genes.csv   n = 312  (audited benchmark, v1.1.0)
No count is estimated.

2026-10-18 (H-18 round) — content updated 320 -> 312 after the post-hoc eligibility
audit (8 signatures removed; JCE_S1, Table A8). The 18 prespecified standardisation
removals (13 <=2 genes, 1 treatment-response outcome, 4 non-reproducible gene list)
and the 8 post-hoc removals are shown together as 26 in the final-standardisation
box, exactly as in the main text ("26 were removed at standardisation (13 <=2 genes,
1 treatment-response outcome, 4 non-reproducible gene list, 8 ineligible on post-hoc
audit; JCE_S1, Table A8)").
Rendering restored to 1500 dpi + bbox_inches='tight' + _trim: the 2026-10-17
regeneration had dropped bbox_inches='tight', which left 20-41 mm of white margin
around a 131 mm diagram (fonts would have printed at ~6.9 pt in a 165.9 mm column).
"""
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import rcParams

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

os.makedirs(OUT, exist_ok=True)

# Elsevier artwork instructions: 只允许 Arial/Helvetica、Courier、Symbol、Times；
# 用非标准字体会被替换为 Elsevier 标准字体，可能导致缺符号或叠字。故弃用默认 DejaVu Sans。
rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['Arial', 'Helvetica']
rcParams['font.size'] = 8.0
rcParams['pdf.fonttype'] = 42
rcParams['ps.fonttype'] = 42      # 矢量导出时嵌入 TrueType，满足 “Embed the font”

GREY = '#5B6770'
BLUE_FILL, BLUE_EDGE = '#DCE6F1', '#2E6DA4'
GREEN_FILL, GREEN_EDGE = '#D5F0E3', '#178C7F'
GREY_FILL = '#F2F2F0'

LINE = 2.1   # vertical space per body line, in data units


def fig_prisma():
    fig, ax = plt.subplots(figsize=(7.4, 10.6))
    ax.axis('off')
    ax.set_xlim(0, 100)
    ax.set_ylim(10, 114)

    LX, LW = 7, 46          # left (flow) column
    RX, RW = 57, 39         # right (exclusion) column

    def box(x, y, w, h, title, body, fc, ec, tl=1, fs_t=8.3, fs_b=7.3):
        """tl = number of lines the title occupies (affects body offset)."""
        ax.add_patch(plt.Rectangle((x, y), w, h, fc=fc, ec=ec, lw=1.0,
                                   joinstyle='round', zorder=2))
        if body:
            ax.text(x + w / 2, y + h - 1.4, title, ha='center', va='top',
                    fontsize=fs_t, fontweight='bold', zorder=3)
            body_y = y + h - 1.4 - 2.5 * tl - 0.5
            ax.text(x + w / 2, body_y, body, ha='center', va='top',
                    fontsize=fs_b, zorder=3, linespacing=1.42)
        else:
            ax.text(x + w / 2, y + h / 2, title, ha='center', va='center',
                    fontsize=fs_t, fontweight='bold', zorder=3)

    def down(y1, y2):
        ax.annotate('', xy=(LX + LW / 2, y2), xytext=(LX + LW / 2, y1),
                    arrowprops=dict(arrowstyle='-|>', lw=1.15, color=GREY,
                                    shrinkA=2, shrinkB=2))

    def across(y_from, y_to):
        ax.annotate('', xy=(RX, y_to), xytext=(LX + LW, y_from),
                    arrowprops=dict(arrowstyle='-|>', lw=0.95, color=GREY,
                                    shrinkA=2, shrinkB=2))

    # ---------------- Identification ----------------
    box(LX, 100, LW, 9.5, 'Records identified from database search',
        'Europe PMC, English language, 2020-2026\nn = 1,241',
        BLUE_FILL, BLUE_EDGE)
    box(RX, 96, RW, 12, 'Records removed before screening',
        'Full text not retrievable\n(no open-access version)\nn = 641',
        GREY_FILL, GREY)
    across(104, 104)
    down(100, 95)

    # ---------------- Screening ----------------
    box(LX, 87, LW, 8, 'Records screened at title level', 'n = 600',
        BLUE_FILL, BLUE_EDGE)
    box(RX, 71.5, RW, 22.5, 'Records excluded at title level',
        'n = 70  (all 70 verified by hand,\nnone reinstated)\n'
        '\u2022 imaging / digital pathology                    11\n'
        '\u2022 nutrition-inflammation score / nomogram  12\n'
        '\u2022 screening or incidence risk, not prognosis  17\n'
        '\u2022 non-coding RNA (lncRNA / miRNA)             14\n'
        '\u2022 methylation-only model                          5\n'
        '\u2022 blood / stool / CTC / microbiome            11',
        GREY_FILL, GREY)
    across(91, 91)
    down(87, 82)

    # ---------------- Eligibility ----------------
    box(LX, 74, LW, 8, 'Full texts assessed for eligibility', 'n = 530',
        BLUE_FILL, BLUE_EDGE)
    box(RX, 57, RW, 12.5, 'Gene list not extractable\nfrom the published full text',
        'n = 150', GREY_FILL, GREY, tl=2)
    across(76, 66)
    down(74, 70)

    box(LX, 58.5, LW, 11.5, 'Gene lists extracted',
        'n = 380   (226 automatic at high\nconfidence, 154 requiring\nmanual confirmation)',
        BLUE_FILL, BLUE_EDGE)
    box(RX, 36.5, RW, 18.5, 'Removed during verification\nand symbol standardisation',
        'n = 42\n'
        '\u2022 non-coding transcripts\n'
        '\u2022 splice-event identifiers\n'
        '\u2022 not mappable to official HGNC\n'
        '\u2022 judged non-mRNA on review',
        GREY_FILL, GREY, tl=2)
    across(62, 50)
    down(58.5, 53.5)

    # ---------------- Included ----------------
    box(LX, 45.5, LW, 8, 'Primary-analysis library', 'n = 338',
        BLUE_FILL, BLUE_EDGE)
    # n = 26 = 18 prespecified standardisation removals (13 / 1 / 4)
    #       + 8 ineligible on the post-hoc eligibility audit (JCE_S1, Table A8)
    box(RX, 18.5, RW, 16, 'Removed during final standardisation',
        'n = 26\n'
        '\u2022 \u22642 genes (fails multigene criterion)   13\n'
        '\u2022 treatment-response outcome               1\n'
        '\u2022 gene list not reproducible on review   4\n'
        '\u2022 ineligible on post-hoc re-read           8',
        GREY_FILL, GREY)
    across(49.5, 30)
    down(45.5, 40.5)

    box(LX, 31, LW, 9.5, 'Signatures included in the benchmark',
        'n = 312 from 312 publications\n1,463 unique HGNC symbols',
        GREEN_FILL, GREEN_EDGE)
    down(31, 25.5)

    box(LX, 18, LW, 7.5,
        'External evaluations\n312 signatures \u00d7 8 cohorts = 2,496',
        None, GREEN_FILL, GREEN_EDGE, fs_t=8.0)

    # ---------------- section labels ----------------
    for name, ymid in [('Identification', 104), ('Screening', 90),
                       ('Eligibility', 65), ('Included', 33)]:
        ax.text(2.2, ymid, name, ha='center', va='center',
                rotation=90, fontsize=8.6, fontweight='bold', color=GREY)

    # 不烧录图题：Elsevier 规定 “A caption should consist of a brief title
    # (not displayed on the figure itself)”，图注已在正文 P56 单独提供。
    # 1500 dpi：本图为纯线条/文字流程，属 line art；全书阔(190 mm)排版时
    # 最低需 7480 px (=1000 dpi x 7.48 in)。7.4 in @1500 dpi 经 tight+trim 后约 8200 px，
    # 留约 10% 余量；实测像素宽见 08_核验底稿。
    _p = f'{OUT}/Figure_1.png'
    _pdf = f'{OUT}/Figure_1.pdf'
    fig.savefig(_p, dpi=1500, bbox_inches='tight', facecolor='white')
    # 同时输出矢量备份（Elsevier 对 line art 优先推荐矢量；仅在需要时提供）
    fig.savefig(_pdf, bbox_inches='tight', facecolor='white')
    # 2026-10-18 补：PNG 与 PDF 用同一内容框裁边。tight 对两者输出同一 bbox，
    # 但会把不可见坐标轴背景算入，PDF 因此比按内容裁过的 PNG 大一圈
    # （实测 150.7x212.4 vs 139.1x183.6 mm），投稿若提交矢量版会被多缩一次、字号跌破 7 pt。
    _crop, _tight_size = _content_box(_p)
    _trim(_p, 1500)
    _trim_pdf(_pdf, _crop, _tight_size)
    plt.close(fig)
    print('Figure 1 (PRISMA) ok ->', _p)


def _content_box(path, pad=24):
    """返回（裁剪框, tight 输出尺寸）；裁剪框按实际内容 + pad 计算。"""
    from PIL import Image
    import numpy as np
    im = Image.open(path)
    bg = Image.new('RGB', im.size, (255, 255, 255))
    bg.paste(im, mask=im.split()[3] if im.mode == 'RGBA' else None)
    a = np.asarray(bg, int)
    ys, xs = np.where(a.sum(2) < 750)
    crop = (max(0, xs.min()-pad), max(0, ys.min()-pad),
            min(a.shape[1], xs.max()+1+pad), min(a.shape[0], ys.max()+1+pad))
    return crop, im.size


def _trim_pdf(pdf_path, crop, size):
    """把 PDF 裁到与 PNG 相同边界（PIL 原点左上 y 向下 → PDF 原点左下 y 向上）。"""
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


def _trim(path, dpi, pad=24):
    """bbox_inches='tight' 会把不可见的坐标轴背景算入，导致上下留 ~10 mm 死白；
    按实际内容重新裁边，使图幅恰好包住图形。"""
    from PIL import Image
    import numpy as np
    im = Image.open(path)
    bg = Image.new('RGB', im.size, (255, 255, 255))
    bg.paste(im, mask=im.split()[3] if im.mode == 'RGBA' else None)
    a = np.asarray(bg, int)
    ys, xs = np.where(a.sum(2) < 750)
    box_ = (max(0, xs.min()-pad), max(0, ys.min()-pad),
            min(a.shape[1], xs.max()+1+pad), min(a.shape[0], ys.max()+1+pad))
    bg.crop(box_).save(path, dpi=(dpi, dpi))


if __name__ == '__main__':
    fig_prisma()
