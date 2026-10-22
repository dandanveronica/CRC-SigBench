# -*- coding: utf-8 -*-
"""Figure 2（研究设计与四个分析臂）312 口径重绘 · 2026-10-18

修复 10-17 重绘引入的三处回归：
  H-18a 白边回归：丢了 bbox_inches='tight' + _trim，内容仅 139.6x88.1 mm/画布 190x132.1 mm，
          按 165.9 mm 栏宽排印字号跌至 6.4 pt（<7 pt 下限）
  H-18b 措辞回归："coefficients as published" 与主稿 P31 "Published coefficients were
          deliberately not transported" 直接矛盾；恢复旧版（10-14）措辞
  H-18c 措辞回归：Arm B "Four-tier comparison" 全文无出处（A1 题注正式名称为
          "Head-to-head comparison"）；恢复旧版措辞

数字全部经 Table 1 / S1 / A1 复核：
  312 signatures / 312 publications / 1,463 unique HGNC symbols（312 库实测）
  7 GEO microarrays n=1,289（556+177+55+124+122+65+190）+ TCGA 358 = 1,647 patients
  events 473（393+80）；312 x 8 = 2,496 evaluations；B = 1,000；16 variables（S2）
"""
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import rcParams

# Elsevier artwork instructions: 只允许 Arial/Helvetica、Courier、Symbol、Times
rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['Arial', 'Helvetica']
rcParams['font.size'] = 8.5
rcParams['pdf.fonttype'] = 42
rcParams['ps.fonttype'] = 42      # 矢量导出嵌入 TrueType

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

os.makedirs(OUT, exist_ok=True)


def _content_box(path, pad=24):
    """返回（裁剪框, tight 输出尺寸）；裁剪框按实际内容 + pad 计算。"""
    from PIL import Image
    import numpy as np
    im = Image.open(path)
    bg = Image.new('RGB', im.size, (255, 255, 255))
    bg.paste(im, mask=im.split()[3] if im.mode == 'RGBA' else None)
    a = np.asarray(bg, int)
    ys, xs = np.where(a.sum(2) < 750)
    crop = (max(0, xs.min() - pad), max(0, ys.min() - pad),
            min(a.shape[1], xs.max() + 1 + pad), min(a.shape[0], ys.max() + 1 + pad))
    return crop, im.size


def _trim(path, dpi, pad=24):
    """bbox_inches='tight' 会把不可见的坐标轴背景算入，按实际内容重新裁边。"""
    from PIL import Image
    crop, _ = _content_box(path, pad)
    im = Image.open(path)
    bg = Image.new('RGB', im.size, (255, 255, 255))
    bg.paste(im, mask=im.split()[3] if im.mode == 'RGBA' else None)
    bg.crop(crop).save(path, dpi=(dpi, dpi))
    return crop


def _trim_pdf(pdf_path, crop, size):
    """把 PDF 裁到与 PNG 相同的边界。

    bbox_inches='tight' 对 PNG 与 PDF 输出同一个 bbox，但会把不可见的坐标轴背景算入，
    使 PDF 比按内容裁剪后的 PNG 大一圈（本图实测 154.7x114.6 vs 142.8x95.7 mm）。
    若投稿提交矢量版，排版时会被多缩一次，字号跌破 7 pt 下限，故须同步裁边。
    坐标换算：PIL 原点左上（y 向下），PDF 原点左下（y 向上）。
    """
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
    fig, ax = plt.subplots(figsize=(7.6, 5.6))
    ax.axis('off')
    ax.set_xlim(0, 100)
    ax.set_ylim(26, 95)

    def box(x, y, w, h, title, body, fc, fs_t=8.4, fs_b=7.8):
        ax.add_patch(plt.Rectangle((x, y), w, h, fc=fc, ec=GREY, lw=1.0,
                                   joinstyle='round', zorder=2))
        ax.text(x + w / 2, y + h - (1.6 if body else h / 2), title,
                ha='center', va='top' if body else 'center',
                fontsize=fs_t, fontweight='bold', zorder=3)
        if body:
            ax.text(x + w / 2, y + 1.2, body, ha='center', va='bottom',
                    fontsize=fs_b, zorder=3, linespacing=1.35)

    def arrow(x1, y1, x2, y2):
        ax.annotate('', xy=(x2, y2), xytext=(x1, y1),
                    arrowprops=dict(arrowstyle='-|>', lw=1.1, color=GREY,
                                    shrinkA=2, shrinkB=2))

    CX, CW = 15, 70
    # 顶部三框（检索 / 筛选条件 / 候选）已删除：其全部细节见 Figure 1 的 PRISMA 流程图。
    # 本图只负责 PRISMA 不覆盖的部分——研究库、队列、统一评分与四个分析臂。
    box(CX, 80, CW, 9.5, '312 signatures   |   312 publications',
        '1,463 unique HGNC symbols\n'
        'identified from 1,241 records; screening flow in Figure 1', '#D5F0E3')
    arrow(50, 80, 50, 76.2)
    box(CX, 66, CW, 9.5, '8 external cohorts',
        '7 GEO microarrays (n = 1,289) + TCGA COADREAD RNA-seq (n = 358)\n'
        'total 1,647 patients, 473 events', '#DCE6F1')
    arrow(50, 66, 50, 63.7)
    # h=11.5：三行正文实占 5.9 data 单位，9.5 高的框会与标题叠印（renderer 实测），
    # 加高后留 0.66 单位间隙（10-14 版实测结论，本次沿用）
    box(CX, 52, CW, 11.5, 'Uniform reconstruction and scoring',
        'within-cohort z-standardisation  |  oriented C-index\n'
        'uniform unweighted z-mean over measured genes; published\n'
        'coefficients deliberately not carried into the evaluation', '#EAF0F7')

    labels = [
        ('Arm A', 'Primary evaluation\n312 \u00d7 8 = 2,496\nexternal validations\n'
                  'C, AUC, calibration,\nDCA', '#DCE6F1'),
        ('Arm B', 'Head-to-head comparison\nvs de novo comparators:\nsize-matched random\n'
                  'null (B = 1,000) +\nclinical model', '#F7E5D8'),
        ('Arm C', 'Attribution & bias\n16 prespecified\nmethodological variables;\n'
                  'claimed vs recomputed', '#D5F0E3'),
        ('Arm D', 'Re-derivation\nminimal-compliant\nprotocol; own\nsignature, same\n'
                  'judgement', '#FDF2CC'),
    ]
    ws = 22.0
    for k, (tag, body, col) in enumerate(labels):
        x = 2.5 + k * (ws + 2.2)
        ax.annotate('', xy=(x + ws / 2, 45), xytext=(50, 52),
                    arrowprops=dict(arrowstyle='-|>', lw=0.9, color=GREY,
                                    connectionstyle='arc3,rad=0',
                                    shrinkA=2, shrinkB=1))
        box(x, 30, ws, 15, tag, body, col, fs_t=8.6, fs_b=7.0)

    # 不烧录图题：Elsevier 规定 "A caption should consist of a brief title
    # (not displayed on the figure itself)"，图注已在正文单独提供。
    # 1400 dpi：本图四臂并排，须全书阔(190 mm)排版，线条图全宽最低 7480 px
    # (=1000 dpi x 7.48 in)；裁边后仍 >7480 px。
    png = os.path.join(OUT, 'Figure_2.png')
    pdf = os.path.join(OUT, 'Figure_2.pdf')
    fig.savefig(png, dpi=1400, bbox_inches='tight', facecolor='white')
    fig.savefig(pdf, bbox_inches='tight', facecolor='white')
    # 先按 tight 输出量出内容框，PNG 与 PDF 用同一比例裁边，保证两者物理尺寸一致
    crop, tight_size = _content_box(png)
    _trim(png, 1400)
    _trim_pdf(pdf, crop, tight_size)
    return fig, ax   # 不在 build() 内 close，便于校验脚本复用同一 figure 做几何检测


if __name__ == '__main__':
    fig, _ = build()
    plt.close(fig)
    print('Figure_2.png / Figure_2.pdf written to', OUT)
