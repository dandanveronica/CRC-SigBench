# -*- coding: utf-8 -*-
"""生成 IMRaD 全文的 4 张 Figure（1400 dpi PNG；组合图全幅 190 mm 有效 >1000 dpi）。
数据全部来自 audit_out_clean/ 下已定稿的 8 队列结果。
2026-10-14：全部改 Arial、删烧录图题、按重排后编号输出（fig1→Figure_2, fig2→Figure_3, fig3→Figure_4, fig4→Figure_5）。

DEPRECATED (v1.2.0, 2026-10-22)：本脚本**不再用于生成投稿图**，保留仅供溯源。原因：
  1. 依赖 audit_out_clean/ 与 audit_out_8/ 下未入库的中间数据（.gitignore 已排除），
     clone 后无法运行；
  2. 签名过滤名单为 320 库（签名320清单.csv），与 v1.1.0 起的 312 口径不符；
  3. fig3()/fig4() 未做投稿级改造（figsize 仍 9.4x7.0 > 190 mm、无 _trim、字号未达标）。
投稿定稿图（Figure_1–Figure_5）请改用同目录下 41–45 号脚本：
  41_figure1_prisma_312.py / 42_figure2_design_312.py / 43_figure3_increments_312.py /
  44_figure4_mechanisms_312.py / 45_figure5_trends_312.py
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import rcParams
from scipy import stats

B = 'audit_out_clean'
OUT = os.path.join(B, 'figures')
os.makedirs(OUT, exist_ok=True)
# Elsevier artwork instructions: 只允许 Arial/Helvetica、Courier、Symbol、Times，
# 用非标准字体会被替换为 Elsevier 标准字体，可能导致缺符号或叠字。故弃用默认 DejaVu Sans。
rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['Arial', 'Helvetica']
rcParams['font.size'] = 8.5
rcParams['axes.linewidth'] = 0.8
rcParams['xtick.major.width'] = 0.8
rcParams['ytick.major.width'] = 0.8
rcParams['pdf.fonttype'] = 42
rcParams['ps.fonttype'] = 42      # 矢量导出时嵌入 TrueType，满足 “Embed the font”

GREY, BLUE, RED, ORANGE, TEAL = '#5B6770', '#2E6DA4', '#C0392B', '#D68910', '#178C7F'

KEEP = set(pd.read_csv(f'{B}/签名320清单.csv').sig_id)   # DEPRECATED: 320 库名单，v1.1.0 起为 312
pc = pd.read_csv('audit_out_8/per_cohort.csv')
pc = pc[pc.sig_id.isin(KEEP)].copy()
pc['d_rand'] = pc.C - pc.null_median
GEO = pc[pc.cohort != 'TCGA_COADREAD']
RNA = pc[pc.cohort == 'TCGA_COADREAD']

# ============================================================ Figure 1: flow
def fig1():
    fig, ax = plt.subplots(figsize=(7.6, 5.6))
    ax.axis('off')
    ax.set_xlim(0, 100); ax.set_ylim(26, 95)

    def box(x, y, w, h, title, body, fc, fs_t=8.4, fs_b=7.8):
        ax.add_patch(plt.Rectangle((x, y), w, h, fc=fc, ec=GREY, lw=1.0,
                                   joinstyle='round', zorder=2))
        ax.text(x + w/2, y + h - (1.6 if body else h/2), title,
                ha='center', va='top' if body else 'center',
                fontsize=fs_t, fontweight='bold', zorder=3)
        if body:
            ax.text(x + w/2, y + 1.2, body, ha='center', va='bottom',
                    fontsize=fs_b, zorder=3, linespacing=1.35)

    def arrow(x1, y1, x2, y2):
        ax.annotate('', xy=(x2, y2), xytext=(x1, y1),
                    arrowprops=dict(arrowstyle='-|>', lw=1.1, color=GREY,
                                    shrinkA=2, shrinkB=2))

    CX, CW = 15, 70
    # 顶部三框（检索 / 筛选条件 / 候选）已删除：其全部细节见 Figure 1 的 PRISMA 流程图。
    # 本图只负责 PRISMA 不覆盖的部分——研究库、队列、统一评分与四个分析臂。
    # 2026-10-18 更新为 312 口径（H-18 轮）：
    #   - 数字：312/312/1,463；7 GEO n=1,289（Table 1 黏膜修正后）+ TCGA 358 = 1,647/473；312x8=2,496
    #   - 恢复 10-14 措辞：published coefficients deliberately not carried into the evaluation
    #     （10-17 重绘曾错改为 "coefficients as published"，与主稿 P31 矛盾）
    #   - 恢复 10-14 措辞：Arm B "Head-to-head comparison"（A1 题注正式名称；
    #     10-17 重绘曾错改为无出处的 "Four-tier comparison"）
    #   - fs_b 6.6 -> 7.0：自然尺寸下最小字号 >=7 pt（Elsevier 下限）
    box(CX, 80, CW, 9.5, '312 signatures   |   312 publications',
        '1,463 unique HGNC symbols\n'
        'identified from 1,241 records; screening flow in Figure 1', '#D5F0E3')
    arrow(50, 80, 50, 76.2)
    box(CX, 66, CW, 9.5, '8 external cohorts',
        '7 GEO microarrays (n = 1,289) + TCGA COADREAD RNA-seq (n = 358)\n'
        'total 1,647 patients, 473 events', '#DCE6F1')
    arrow(50, 66, 50, 63.7)
    # h=11.5：三行正文(7.0pt)实占 5.9 data 单位，9.5 高的框会与标题叠印 7.7pt（renderer 实测），加高后留 0.66 单位间隙
    box(CX, 52, CW, 11.5, 'Uniform reconstruction and scoring',
        'within-cohort z-standardisation  |  oriented C-index\n'
        'uniform unweighted z-mean over measured genes; published\n'
        'coefficients deliberately not carried into the evaluation', '#EAF0F7')

    labels = [('Arm A', 'Primary evaluation\n312 \u00d7 8 = 2,496\nexternal validations\n'
                        'C, AUC, calibration,\nDCA', '#DCE6F1'),
              ('Arm B', 'Head-to-head comparison\nvs de novo comparators:\nsize-matched random\n'
                        'null (B = 1,000) +\nclinical model', '#F7E5D8'),
              ('Arm C', 'Attribution & bias\n16 prespecified\nmethodological variables;\n'
                        'claimed vs recomputed', '#D5F0E3'),
              ('Arm D', 'Re-derivation\nminimal-compliant\nprotocol; own\nsignature, same\n'
                        'judgement', '#FDF2CC')]
    ws = 22.0
    for k, (tag, body, col) in enumerate(labels):
        x = 2.5 + k * (ws + 2.2)
        ax.annotate('', xy=(x + ws/2, 45), xytext=(50, 52),
                    arrowprops=dict(arrowstyle='-|>', lw=0.9, color=GREY,
                                    connectionstyle='arc3,rad=0',
                                    shrinkA=2, shrinkB=1))
        box(x, 30, ws, 15, tag, body, col, fs_t=8.6, fs_b=7.0)

    # 不烧录图题：Elsevier 规定 “A caption should consist of a brief title
    # (not displayed on the figure itself)”，图注已在正文 P56 单独提供。
    # 1400 dpi：本图四臂并排，须全书阔(190 mm)排版，线条图全宽最低 7480 px
    # (=1000 dpi x 7.48 in)；裁边后仍 >7480 px。
    # 输出名即投稿文件名（图序重排 2026-10-14 后：研究设计图 = Figure 2）
    _p = f'{OUT}/Figure_2.png'
    fig.savefig(_p, dpi=1400, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    _trim(_p, 1400)
    print('Figure 1 ok')


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

# ============================================================ Figure 3: increments
def fig2():
    # figsize 由 (9.4, 7.0) 收到 (7.2, 5.36)：原尺寸 tight 后成品物理宽 232.8 mm，超过期刊
    # 全幅 190 mm，排版时被缩小 0.82x，最小字号掉到 5.71 pt（Elsevier 要求 >=7 pt）。
    # 收到 178 mm 后排版放大 1.07x，最小印刷字号 7.46 pt，达标。fig4/fig5 同病，需同法处理。
    fig, axs = plt.subplots(2, 2, figsize=(7.2, 5.36))
    # 不烧录图题（Elsevier：图注不得显示在图上，题注在正文单独提供）；重排后本图 = Figure 3

    # A
    ax = axs[0, 0]
    ax.hist(pc.d_rand.dropna(), bins=70, color=BLUE, alpha=.85, edgecolor='white', lw=.2)
    ax.axvline(0, color='k', lw=.9)
    ax.axvline(pc.d_rand.median(), color=RED, lw=1.4, ls='--')
    ax.text(.02, .93, f'median {pc.d_rand.median():+.4f}\n{100*(pc.d_rand>0).mean():.1f}% positive\np = 1.5e-48',
            transform=ax.transAxes, va='top', fontsize=7.6,
            bbox=dict(fc='white', ec=GREY, lw=.6, pad=3))
    # 注意：轴标签不得用 matplotlib mathtext($...$)。mathtext.fontset 默认 dejavusans，
    # $\Delta$ 会由 DejaVu Sans 渲染，违反 Elsevier「只允许 Arial/Courier/Symbol/Times」。
    # Arial 含 U+0394 字形（已实测），故直接用 Unicode Δ。
    ax.set_xlabel('ΔC_random  (signature C - matched random-gene-set median)')
    ax.set_ylabel('Number of evaluations')
    ax.set_title('A  2,496 evaluations vs random gene sets', loc='left', fontsize=9, fontweight='bold')

    # B
    ax = axs[0, 1]
    ax.hist(pc.delta_C.dropna(), bins=70, color=TEAL, alpha=.85, edgecolor='white', lw=.2)
    ax.axvline(0, color='k', lw=.9)
    ax.axvline(pc.delta_C.median(), color=RED, lw=1.4, ls='--')
    ax.text(.02, .93, f'median {pc.delta_C.median():+.4f}\n{100*(pc.delta_C>0).mean():.1f}% positive\np = 4.2e-244',
            transform=ax.transAxes, va='top', fontsize=7.6,
            bbox=dict(fc='white', ec=GREY, lw=.6, pad=3))
    ax.set_xlabel('ΔC_clinical  (clinical + signature, minus clinical alone)')
    ax.set_ylabel('Number of evaluations')
    ax.set_title('B  Increment above routine clinical variables', loc='left', fontsize=9, fontweight='bold')

    # C
    ax = axs[1, 0]
    groups = [r'random: record', r'random: signature', r'clinical: record', r'clinical: signature']
    vals = [100 * (pc.d_rand > 0).mean()]
    S = pc.groupby('sig_id').agg(dr=('d_rand', 'mean'), dc=('delta_C', 'mean'))
    vals += [100 * (S.dr > 0).mean(), 100 * (pc.delta_C > 0).mean(), 100 * (S.dc > 0).mean()]
    cols = [BLUE, '#7FB3DB', TEAL, '#7FD1BF']
    bars = ax.bar(range(4), vals, color=cols, edgecolor=GREY, lw=.7, width=.62)
    for b, v in zip(bars, vals):
        ax.text(b.get_x() + b.get_width() / 2, v + 1.2, f'{v:.1f}%', ha='center', fontsize=7.8)
    ax.axhline(50, color=RED, ls=':', lw=1.0)
    ax.text(3.55, 52, 'chance', color=RED, ha='right', fontsize=7)
    # 18 度足够：相邻标签锚点间距 70.7 px(dpi100)，平行标签垂直间距 = 70.7*sin18 = 21.9 px
    # > 字高 11.8 px(8.5pt)，故不重叠。注意：旋转标签的 AABB 必然相交，用 AABB 判叠印会假阳性。
    ax.set_xticks(range(4)); ax.set_xticklabels(groups, rotation=18, ha='right')
    ax.set_ylim(0, 108); ax.set_ylabel('% of evaluations / signatures positive')
    ax.set_title('C  Evaluation noise resolves into direction once averaged', loc='left', fontsize=9, fontweight='bold')

    # D
    ax = axs[1, 1]
    rows = []
    for c, s in pc.groupby('cohort'):
        rows.append((c, s.d_rand.median(), s.d_rand.quantile(.25), s.d_rand.quantile(.75), s.C.median()))
    rows.sort(key=lambda r: r[1])
    # 先去掉队列名里的下划线（TCGA_COADREAD -> TCGA COADREAD），再换行，避免行尾残留 "TCGA_"
    names = [r[0].replace('_', ' ').replace(' COADREAD', '\nCOADREAD') for r in rows]
    ypos = np.arange(len(rows))
    for k, (nm, md, lo, hi, cm) in enumerate(rows):
        col = ORANGE if md < 0 else (BLUE if nm.startswith('TCGA') else GREY)
        ax.plot([lo, hi], [k, k], color=col, lw=1.6, solid_capstyle='round')
        ax.plot(md, k, 'o', ms=4.6, color=col)
    ax.axvline(0, color='k', lw=.9)
    ax.set_yticks(ypos); ax.set_yticklabels(names, fontsize=7.4)
    ax.set_xlabel('Median ΔC_random with interquartile range')
    ax.set_title('D  Per-cohort: GSE87211 is the only negative cohort', loc='left', fontsize=9, fontweight='bold')
    ax.grid(axis='y', ls=':', lw=.4, color='#BBBBBB')

    fig.tight_layout(rect=[0, 0, 1, 1])
    _p = f'{OUT}/Figure_3.png'
    fig.savefig(_p, dpi=1400, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    _trim(_p, 1400)
    print('Figure 3 ok')

# ============================================================ Figure 4: three mechanisms
def fig3():
    fig, axs = plt.subplots(2, 2, figsize=(9.4, 7.0))
    # 不烧录图题；重排后本图 = Figure 4

    # A cutpoint
    ax = axs[0, 0]
    cuts = pd.read_csv(f'{B}/截点操纵_汇总.csv')
    pub = cuts.iloc[0]
    rnd = cuts.iloc[1]
    x = np.arange(2); w = .36
    ax.bar(x - w / 2, [100 * pub['预指定_显著率'], 100 * rnd['预指定_显著率']], w,
           color=GREY, edgecolor='k', lw=.6, label='pre-specified cut-point')
    ax.bar(x + w / 2, [100 * pub['最优截点_显著率'], 100 * rnd['最优截点_显著率']], w,
           color=RED, edgecolor='k', lw=.6, label='best cut-point by search')
    ax.axhline(5, color=BLUE, ls='--', lw=1.1)
    ax.text(1.48, 6, 'nominal 5%', color=BLUE, ha='right', fontsize=7)
    for i, v in enumerate([100 * pub['预指定_显著率'], 100 * rnd['预指定_显著率']]):
        ax.text(i - w / 2, v + 1.4, f'{v:.1f}%', ha='center', fontsize=7.4, zorder=6)
    for i, v in enumerate([100 * pub['最优截点_显著率'], 100 * rnd['最优截点_显著率']]):
        ax.text(i + w / 2, v + 1.4, f'{v:.1f}%', ha='center', fontsize=7.4, zorder=6)
    ax.set_xticks(x); ax.set_xticklabels(['Published\nsignatures (n=2,240)', 'Random\ngene sets (n=3,280)'])
    ax.set_ylabel('% reaching p < 0.05'); ax.set_ylim(0, 50)
    ax.legend(fontsize=7, frameon=True, loc='upper center', ncol=2, framealpha=.95)
    ax.set_title('A  Shared grey-to-red gain is the search itself', loc='left', fontsize=9, fontweight='bold')

    # B NRI significance
    ax = axs[0, 1]
    cv = pd.read_csv(f'{B}/NRI_IDI_CV_签名.csv'); cvr = pd.read_csv(f'{B}/NRI_IDI_CV_随机基线.csv')
    issig = pd.read_csv(f'{B}/per_cohort_NRI_IDI.csv'); isrnd = pd.read_csv(f'{B}/NRI_IDI_随机基线.csv')
    sets = [('In-sample', 100 * (issig.NRI_p < .05).mean(), 100 * (isrnd.NRI_p < .05).mean()),
            ('10-fold CV', 100 * (cv.NRI_cv_p < .05).mean(), 100 * (cvr.NRI_cv_p < .05).mean())]
    x = np.arange(2); w = .36
    ax.bar(x - w / 2, [s[1] for s in sets], w, color=TEAL, edgecolor='k', lw=.6, label='published signatures')
    ax.bar(x + w / 2, [s[2] for s in sets], w, color='#C9CBCD', edgecolor='k', lw=.6, label='random gene sets')
    ax.axhline(5, color=BLUE, ls='--', lw=1.1)
    ax.text(1.48, 6, 'nominal 5%', color=BLUE, ha='right', fontsize=7)
    for i, (_, a, b) in enumerate(sets):
        ax.text(i - w / 2, a + .8, f'{a:.1f}%', ha='center', fontsize=7.4)
        ax.text(i + w / 2, b + .8, f'{b:.1f}%', ha='center', fontsize=7.4)
    ax.annotate('gap 1.6 pp\n(p = 0.356)', xy=(-w / 2, 21.5), fontsize=7, color=RED, ha='left')
    ax.annotate('gap 9.7 pp', xy=(1 - w / 2, 26.6), fontsize=7, color=RED, ha='left')
    ax.set_xticks(x); ax.set_xticklabels([s[0] for s in sets])
    ax.set_ylabel('% of NRI tests reaching p < 0.05'); ax.set_ylim(0, 36)
    ax.legend(fontsize=7, loc='upper left')
    ax.set_title('B  In-sample NRI is anti-conservative (3.2x nominal)', loc='left', fontsize=9, fontweight='bold')

    # C claimed vs recomputed
    ax = axs[1, 0]
    # C 面板数据源用投稿包定稿口径 S7（人工核对后 claim 值，median gap = +0.176）。
    # 初版 原文声称值vs重算值.csv 后来被自动化重跑覆盖成 +0.174，与正文/复算锁定的 0.176 矛盾，
    # 故直接读 S7_复现率配对明细.csv（n=103，同列名 claim_auc5 / 重算AUC5中位）。
    # v1.2.0：原为本机投稿包绝对路径（不可复现），改为仓库内相对路径
    cl = pd.read_csv(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..',
                                  'results', 'submission_tables', 'S7_复现率配对明细.csv'))
    d = cl.dropna(subset=['claim_auc5', '重算AUC5中位'])
    ax.scatter(d['重算AUC5中位'], d.claim_auc5, s=16, c=RED, alpha=.55, edgecolors='none')
    lo, hi = 0.3, 1.0
    ax.plot([lo, hi], [lo, hi], color='k', lw=.9)
    gap = (d.claim_auc5 - d['重算AUC5中位']).median()
    ax.plot([lo, hi - gap], [lo + gap, hi], color=RED, lw=1.3, ls='--')
    ax.text(.33, .95, f'median gap +{gap:.3f}\n{100*((d.claim_auc5-d["重算AUC5中位"])>0).mean():.0f}% above the diagonal\nn = {len(d)},  p = 1.4e-18',
            fontsize=7.6, va='top', bbox=dict(fc='white', ec=GREY, lw=.6, pad=3))
    ax.set_xlabel('Recomputed 5-year AUC (this study)')
    ax.set_ylabel('Published claimed 5-year AUC')
    ax.set_xlim(lo, hi); ax.set_ylim(lo, hi)
    ax.set_title('C  Claimed versus recomputed', loc='left', fontsize=9, fontweight='bold')

    # D training-to-validation
    ax = axs[1, 1]
    P = pd.read_csv(f'{B}/Part3_合规签名_验证.csv')
    v = P[P.cohort != 'GSE39582'].sort_values('C')
    nm = [({'TCGA_COADREAD': 'TCGA'}.get(c, c)) for c in v.cohort]
    yp = np.arange(len(v))
    ax.scatter([0.8008] * 1, [len(v)], s=70, marker='D', c=RED, zorder=5, label='training C = 0.801')
    for k, c in enumerate(v.C):
        ax.scatter(c, k, s=34, c=BLUE, zorder=4)
    ax.plot([0.8008, v.C.median()], [len(v), len(v) / 2], color=RED, lw=1.7, ls='--')
    ax.hlines(len(v) / 2, v.C.min(), 0.8008, ls=':', lw=.9, color=GREY)
    ax.text(0.62, len(v) - .35, 'training\n0.801', fontsize=7.4, ha='center', color=RED)
    ax.text(v.C.median() + .012, len(v) / 2 + .28, f'validation median {v.C.median():.3f}',
            fontsize=7.4, ha='left', color=BLUE)
    ax.text(0.815, 1.0, 'drop 0.223\nliterature gap 0.174\n- same number', fontsize=7.4,
            color=RED, ha='left', va='bottom')
    ax.set_yticks(list(yp) + [len(v)])
    ax.set_yticklabels(list(nm) + ['training\nGSE39582'], fontsize=7.2)
    ax.set_xlabel('C-index'); ax.set_xlim(0.50, 0.86)
    ax.grid(axis='x', ls=':', lw=.4, color='#BBBBBB')
    ax.set_title('D  Our own compliant signature: the same drop', loc='left', fontsize=9, fontweight='bold')

    fig.tight_layout(rect=[0, 0, 1, 1])
    _p = f'{OUT}/Figure_4.png'
    fig.savefig(_p, dpi=1400, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    _trim(_p, 1400)
    print('Figure 4 ok')

# ============================================================ Figure 5: seven years
def fig4():
    fig, axs = plt.subplots(2, 2, figsize=(9.4, 7.0))
    # 不烧录图题；重排后本图 = Figure 5

    lib = pd.read_csv('sigminer/output/signature_library_main.csv')
    lib = lib[lib.sig_id.isin(KEEP)].copy()
    U = pd.read_csv(f'{B}/签名级可用性判决.csv')

    # A
    ax = axs[0, 0]
    yr = lib.groupby('year').size()
    dr = pc.assign(y=pc.sig_id.map(dict(zip(lib.sig_id, lib.year)))) \
           .merge(lib[['sig_id', 'year']], on='sig_id', how='left', suffixes=('', '_y2'))
    years = sorted(dr['year'].dropna().unique())
    box_dat = [dr[dr['year'] == y].d_rand.dropna().values for y in years]
    bp = ax.boxplot(box_dat, tick_labels=[int(y) for y in years], widths=.55, patch_artist=True,
                    medianprops=dict(color=RED, lw=1.4), showfliers=False)
    for b in bp['boxes']:
        b.set(facecolor='#C8DCEC', edgecolor=GREY, lw=.8)
    ax.axhline(0, color='k', lw=.9)
    r, p = stats.spearmanr(dr.dropna(subset=['d_rand'])['year'], dr.dropna(subset=['d_rand']).d_rand)
    ax.text(.02, .94, f'record level: rho = {r:+.3f}, p = {p:.2f}\n(signature level: rho = +0.014, p = 0.806)',
            transform=ax.transAxes, va='top', fontsize=7.4,
            bbox=dict(fc='white', ec=GREY, lw=.6, pad=3))
    ax.set_xlabel('Publication year'); ax.set_ylabel('ΔC_random')
    ax.set_title('A  No trend in true incremental value', loc='left', fontsize=9, fontweight='bold')

    # B
    ax = axs[0, 1]
    tr = pd.read_csv(f'{B}/时间趋势.csv')
    ax.plot(tr.yr, 100 * tr['外部验证率'], 'o-', color=BLUE, lw=1.6, ms=4.5, label='claims external validation')
    ax.plot(tr.yr, 100 * tr['公开公式率'], 's--', color=RED, lw=1.6, ms=4.5, label='publishes usable coefficients')
    ax.plot(tr.yr, 100 * tr['单细胞率'], '^:', color=TEAL, lw=1.4, ms=4.2, label='single-cell derived')
    _, pv = stats.spearmanr(tr.yr, tr['外部验证率'])
    ax.text(.02, .06, f'claimed validation rho = +0.136, p = 0.015\nverifiability falling 83.3% -> 35.3%',
            transform=ax.transAxes, fontsize=7.4, va='bottom',
            bbox=dict(fc='white', ec=GREY, lw=.6, pad=3))
    ax.set_ylim(0, 105); ax.set_xlabel('Publication year'); ax.set_ylabel('% of publications')
    ax.legend(fontsize=7, loc='center left')
    ax.set_title('B  Reported standards rise, verifiability falls', loc='left', fontsize=9, fontweight='bold')

    # C
    ax = axs[1, 0]
    X = pd.read_csv(f'{B}/跨平台配对_8队列.csv').dropna(subset=['d_rand_geo', 'd_rand_rna'])
    ax.scatter(X.d_rand_geo, X.d_rand_rna, s=13, c=BLUE, alpha=.45, edgecolors='none')
    r, p = stats.spearmanr(X.d_rand_geo, X.d_rand_rna)
    ax.axvline(0, color='k', lw=.7); ax.axhline(0, color='k', lw=.7)
    ax.text(.02, .95, f'Spearman rho = {r:+.3f}\np = {p:.1e}\nn = 312 signatures',
            transform=ax.transAxes, va='top', fontsize=7.6,
            bbox=dict(fc='white', ec=GREY, lw=.6, pad=3))
    ax.set_xlabel('ΔC_random averaged over 7 microarray cohorts')
    ax.set_ylabel('ΔC_random on TCGA RNA-seq')
    ax.set_title('C  Cross-platform reproducibility is weak', loc='left', fontsize=9, fontweight='bold')

    # D
    ax = axs[1, 1]
    ax.scatter(U.C_mean, U.calib, s=13, c=GREY, alpha=.45, edgecolors='none')
    ax.axhline(0.7, color=ORANGE, ls='--', lw=1.0)
    ax.axvline(0.60, color=ORANGE, ls='--', lw=1.0)
    ax.text(.015, .93, f'usable region (C>=0.60, slope>=0.70):\n{int(((U.C_mean>=0.60)&(U.calib>=0.70)).sum())} of 312 signatures',
            transform=ax.transAxes, va='top', fontsize=7.4,
            bbox=dict(fc='white', ec=RED, lw=.6, pad=3))
    top = U.nlargest(5, 'C_mean').sort_values('calib')
    # 标签错位排布：点彼此靠近时（SIG064 与 SIG284 的 C 与 slope 都几乎相同）
    # 统一 (3,3) 偏移会让两段文字叠印成不可读的一团，故按 calib 升序逐个判断，
    # 与前一个标签足够近则把本标签上抬一级。仅影响标注位置，不改任何数据点。
    prev = None
    lvl = 0
    for _, t in top.iterrows():
        if prev is not None and abs(t.C_mean - prev[0]) < 0.012 and abs(t.calib - prev[1]) < 0.06:
            lvl += 1
        else:
            lvl = 0
        ax.annotate(t.sig_id, (t.C_mean, t.calib), fontsize=6.6, xytext=(3, 3 + 11 * lvl),
                    textcoords='offset points', color=RED)
        prev = (t.C_mean, t.calib)
    ax.set_xlabel('Mean C-index across 8 cohorts')
    ax.set_ylabel('Mean calibration slope (ideal = 1.0)')
    ax.set_title('D  The clinically usable quadrant is nearly empty', loc='left', fontsize=9, fontweight='bold')

    fig.tight_layout(rect=[0, 0, 1, 1])
    _p = f'{OUT}/Figure_5.png'
    fig.savefig(_p, dpi=1400, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    _trim(_p, 1400)
    print('Figure 5 ok')


if __name__ == '__main__':
    fig1(); fig2(); fig3(); fig4()
    for f in sorted(os.listdir(OUT)):
        print(' ', f, round(os.path.getsize(f'{OUT}/{f}') / 1024), 'KB')
