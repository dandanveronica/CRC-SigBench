# -*- coding: utf-8 -*-
"""Figure 5 — PRISMA 2020 flow diagram, adapted to signature-level identification.

All counts are taken verbatim from the auditable screening register
(sigminer/output/PRISMA_counts.txt) and from the two frozen library files:
    signature_library_main.csv   n = 338
    signature_library_clean.csv  n = 320
No count is estimated.
"""
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import rcParams

OUT = 'audit_out_clean/figures'
os.makedirs(OUT, exist_ok=True)

rcParams['font.family'] = 'DejaVu Sans'
rcParams['font.size'] = 8.0
rcParams['pdf.fonttype'] = 42

GREY = '#5B6770'
BLUE_FILL, BLUE_EDGE = '#DCE6F1', '#2E6DA4'
GREEN_FILL, GREEN_EDGE = '#D5F0E3', '#178C7F'
GREY_FILL = '#F2F2F0'

LINE = 2.1   # vertical space per body line, in data units


def fig5():
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
    box(RX, 20, RW, 14, 'Removed during final standardisation',
        'n = 18\n'
        '\u2022 \u22642 genes (fails multigene criterion)   13\n'
        '\u2022 treatment-response outcome               1\n'
        '\u2022 gene list not reproducible on review   4',
        GREY_FILL, GREY)
    across(49.5, 30)
    down(45.5, 40.5)

    box(LX, 31, LW, 9.5, 'Signatures included in the benchmark',
        'n = 320 from 320 publications\n1,514 unique HGNC symbols',
        GREEN_FILL, GREEN_EDGE)
    down(31, 25.5)

    box(LX, 18, LW, 7.5,
        'External evaluations\n320 signatures \u00d7 8 cohorts = 2,560',
        None, GREEN_FILL, GREEN_EDGE, fs_t=8.0)

    # ---------------- section labels ----------------
    for name, ymid in [('Identification', 104), ('Screening', 90),
                       ('Eligibility', 65), ('Included', 33)]:
        ax.text(2.2, ymid, name, ha='center', va='center',
                rotation=90, fontsize=8.6, fontweight='bold', color=GREY)

    ax.text(50, 112.6, 'Figure 5. PRISMA 2020 flow diagram (signature-level identification)',
            ha='center', fontsize=9.6, fontweight='bold')

    fig.savefig(f'{OUT}/Figure5_prisma.png', dpi=300, bbox_inches='tight',
                facecolor='white')
    plt.close(fig)
    print('Figure 5 (PRISMA) ok ->', f'{OUT}/Figure5_prisma.png')


if __name__ == '__main__':
    fig5()
