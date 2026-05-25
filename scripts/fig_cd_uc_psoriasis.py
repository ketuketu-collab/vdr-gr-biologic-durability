#!/usr/bin/env python3
"""
CD vs UC vs Psoriasis focused figure
疾患病態とVDR-domain/GR-domain仮説の対応を可視化
"""

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.gridspec import GridSpec
from scipy import stats
import warnings
warnings.filterwarnings('ignore')

RESULTS = '/Volumes/M4_SSD/projects/tlr_chipseq/results/'
COL_VDR  = '#2166AC'
COL_GR   = '#D6604D'
COL_APP  = '#1a9641'
COL_FAIL = '#d7191c'

# ── データ読み込み ──────────────────────────────────────────────────
df = pd.read_csv(RESULTS + 'gene_disease_phase_expanded.csv')

gene_dis = (df.groupby(['gene','disease','category'])
            .agg(
                best_status = ('mol_status', lambda x:
                               'approved' if 'approved' in x.values else 'failed'),
                VDR_score   = ('VDR_score','first'),
                GR_score    = ('GR_score','first'),
                VDR_dom     = ('VDR_dominant','first'),
            ).reset_index())

TARGET_DIS = ['Crohn Disease', 'Colitis, Ulcerative', 'Psoriasis']

# 病態アノテーション
DISEASE_META = {
    'Crohn Disease': {
        'short': 'Crohn\nDisease',
        'organ': 'Small Intestine',
        'vdr_tpm': 38.1,
        'immunotype': 'Th1 / Th17',
        'steroid_resp': 'Partial /\nResistant',
        'domain': 'VDR-domain',
        'color': '#2166AC',
        'note': 'Transmural; VitD deficiency\nassociated with relapse risk',
    },
    'Colitis, Ulcerative': {
        'short': 'Ulcerative\nColitis',
        'organ': 'Colon',
        'vdr_tpm': 29.7,
        'immunotype': 'Th2 / Treg',
        'steroid_resp': 'Responsive',
        'domain': 'GR-domain',
        'color': '#D6604D',
        'note': 'Mucosal; GR-driven\nsteroid efficacy high',
    },
    'Psoriasis': {
        'short': 'Psoriasis',
        'organ': 'Skin',
        'vdr_tpm': 45.2,
        'immunotype': 'Th17 / Th1',
        'steroid_resp': 'Partial\n(topical)',
        'domain': 'VDR-domain',
        'color': '#4393C3',
        'note': 'Epidermal; VDR highest\nin skin (TPM=45.2)',
    },
}

# Fisher's test 計算
def fisher_disease(grp):
    n_va = int((grp['VDR_dom'] & (grp['best_status']=='approved')).sum())
    n_vf = int((grp['VDR_dom'] & (grp['best_status']=='failed')).sum())
    n_ga = int(((~grp['VDR_dom']) & (grp['best_status']=='approved')).sum())
    n_gf = int(((~grp['VDR_dom']) & (grp['best_status']=='failed')).sum())

    # p-value: always use actual Fisher (scipy handles zero cells)
    _, p_v = stats.fisher_exact([[n_va, n_vf],[n_ga, n_gf]], alternative='greater')

    # OR and CI: use 0.5 zero-cell correction throughout for numerical stability
    a, b, c, d = n_va+0.5, n_vf+0.5, n_ga+0.5, n_gf+0.5
    or_v = (a*d)/(b*c)
    log_or = np.log(or_v)
    se = np.sqrt(1/a + 1/b + 1/c + 1/d)
    ci_lo = np.exp(log_or - 1.96*se)
    ci_hi = np.exp(log_or + 1.96*se)

    return dict(n_va=n_va, n_vf=n_vf, n_ga=n_ga, n_gf=n_gf,
                or_v=or_v, ci_lo=ci_lo, ci_hi=ci_hi, p=p_v)

stats_dict = {}
for dis in TARGET_DIS:
    grp = gene_dis[gene_dis['disease'] == dis]
    stats_dict[dis] = fisher_disease(grp)
    s = stats_dict[dis]
    print(f"{dis}: n_va={s['n_va']} n_vf={s['n_vf']} n_ga={s['n_ga']} n_gf={s['n_gf']} "
          f"OR={s['or_v']:.2f} ({s['ci_lo']:.2f}–{s['ci_hi']:.2f}) p={s['p']:.3f}")

# ══════════════════════════════════════════════════════════════════════
# FIGURE
# Layout:
#   Row 0 (tall): 3 x 承認率バー (CD, UC, Psoriasis) + annotation
#   Row 1 (mid):  3 x 主要承認遺伝子 VDR/GR scatter
#   Row 2 (short): OR forest comparison strip
# ══════════════════════════════════════════════════════════════════════
fig = plt.figure(figsize=(15, 13))
gs  = GridSpec(3, 3, figure=fig,
               height_ratios=[3.5, 3, 2],
               hspace=0.55, wspace=0.38)

# ──────────────────────────────────────────────────────────────────────
# ROW 0: 承認率バー + disease annotation
# ──────────────────────────────────────────────────────────────────────
for col_idx, dis in enumerate(TARGET_DIS):
    ax = fig.add_subplot(gs[0, col_idx])
    meta = DISEASE_META[dis]
    s    = stats_dict[dis]
    grp  = gene_dis[gene_dis['disease'] == dis]

    vdr_app_rate = s['n_va']/(s['n_va']+s['n_vf']) if (s['n_va']+s['n_vf'])>0 else 0
    gr_app_rate  = s['n_ga']/(s['n_ga']+s['n_gf']) if (s['n_ga']+s['n_gf'])>0 else 0

    bars = ax.bar([0, 1], [gr_app_rate*100, vdr_app_rate*100],
                  color=[COL_GR, COL_VDR], alpha=0.85, width=0.6,
                  edgecolor='white', linewidth=0.5)
    for bar, rate, n in zip(bars,
                             [gr_app_rate, vdr_app_rate],
                             [s['n_ga']+s['n_gf'], s['n_va']+s['n_vf']]):
        ax.text(bar.get_x()+bar.get_width()/2, rate*100 + 1.5,
                f"{rate*100:.0f}%\n(n={n})",
                ha='center', va='bottom', fontsize=9.5, fontweight='bold')

    # OR annotation box
    p_str = f"p={s['p']:.3f}" if s['p'] < 1 else "p=n.s."
    if s['p'] < 0.1: p_str += " †"
    if s['p'] < 0.05: p_str += "*"
    or_text = (f"OR = {s['or_v']:.1f}\n"
               f"95%CI [{s['ci_lo']:.1f}–{s['ci_hi']:.1f}]\n"
               f"{p_str}")
    ax.text(0.97, 0.97, or_text, transform=ax.transAxes,
            ha='right', va='top', fontsize=8.5,
            bbox=dict(boxstyle='round,pad=0.3', facecolor='#f0f0f0', alpha=0.8))

    # domain badge
    dom_col = '#2166AC' if meta['domain']=='VDR-domain' else '#D6604D'
    ax.text(0.03, 0.97, meta['domain'], transform=ax.transAxes,
            ha='left', va='top', fontsize=9, fontweight='bold', color=dom_col,
            bbox=dict(boxstyle='round,pad=0.3', facecolor=dom_col, alpha=0.15))

    ax.set_xticks([0, 1])
    ax.set_xticklabels(['GR-dom', 'VDR-dom'], fontsize=10)
    ax.set_ylabel('Approval rate (%)', fontsize=10)
    ax.set_ylim(0, 115)
    ax.set_title(f"{meta['short']}\n[{meta['organ']}, VDR TPM={meta['vdr_tpm']}]",
                 fontsize=11, fontweight='bold')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    # annotation footer
    annot = (f"Immunotype: {meta['immunotype']}\n"
             f"Steroid: {meta['steroid_resp']}\n"
             f"{meta['note']}")
    ax.text(0.5, -0.38, annot, transform=ax.transAxes,
            ha='center', va='top', fontsize=7.8, color='#444444',
            style='italic', multialignment='center')

# ──────────────────────────────────────────────────────────────────────
# ROW 1: approved drug target VDR vs GR scatter per disease
# ──────────────────────────────────────────────────────────────────────
for col_idx, dis in enumerate(TARGET_DIS):
    ax = fig.add_subplot(gs[1, col_idx])
    grp = gene_dis[gene_dis['disease'] == dis]
    app  = grp[grp['best_status']=='approved']
    fail = grp[grp['best_status']=='failed']

    # scatter: 失敗=灰, 承認=色
    ax.scatter(fail['GR_score'], fail['VDR_score'],
               color='#bbbbbb', s=40, alpha=0.7, zorder=2, label='Failed')
    ax.scatter(app['GR_score'],  app['VDR_score'],
               color=COL_APP, s=70, alpha=0.9, zorder=3, label='Approved',
               edgecolors='white', linewidths=0.5)

    # 対角線 VDR=GR
    lim_max = max(grp[['VDR_score','GR_score']].max().max(), 20) * 1.05
    ax.plot([0, lim_max], [0, lim_max], color='#aaaaaa', lw=1, ls='--', alpha=0.6)
    ax.fill_between([0, lim_max], [0, lim_max], [lim_max, lim_max],
                    color=COL_VDR, alpha=0.04)
    ax.text(lim_max*0.05, lim_max*0.9, 'VDR > GR', fontsize=7.5,
            color=COL_VDR, alpha=0.7)

    # 主要遺伝子ラベル
    for _, row in app.iterrows():
        ax.annotate(row['gene'],
                    xy=(row['GR_score'], row['VDR_score']),
                    xytext=(4, 4), textcoords='offset points',
                    fontsize=7.5, color=COL_APP, fontweight='bold')

    ax.set_xlabel('GR score', fontsize=9)
    ax.set_ylabel('VDR score', fontsize=9)
    ax.set_xlim(-2, lim_max)
    ax.set_ylim(-2, lim_max)
    ax.set_title(f"Approved targets: VDR vs GR\n({dis.split(',')[0]})",
                 fontsize=9.5, fontweight='bold')
    if col_idx == 0:
        ax.legend(fontsize=8, loc='lower right')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

# ──────────────────────────────────────────────────────────────────────
# ROW 2: OR comparison strip (forest-style horizontal)
# ──────────────────────────────────────────────────────────────────────
ax_forest = fig.add_subplot(gs[2, :])

y_positions = [2, 1, 0]
y_labels    = [DISEASE_META[d]['short'].replace('\n',' ') for d in TARGET_DIS]
domain_cols = [DISEASE_META[d]['color'] for d in TARGET_DIS]

for i, (dis, yp) in enumerate(zip(TARGET_DIS, y_positions)):
    s    = stats_dict[dis]
    meta = DISEASE_META[dis]
    dom  = meta['domain']
    col  = meta['color']
    ci_hi_plot = min(s['ci_hi'], 25)

    ax_forest.plot([s['ci_lo'], ci_hi_plot], [yp, yp],
                   color=col, lw=3, solid_capstyle='round')
    if s['ci_hi'] > 25:
        ax_forest.annotate('', xy=(26, yp), xytext=(25, yp),
                           arrowprops=dict(arrowstyle='->', color=col, lw=1.5))

    marker = 'D' if s['p'] < 0.1 else 'o'
    ax_forest.plot(s['or_v'], yp, marker, color=col, ms=11, zorder=5)

    # label
    p_str = f"p={s['p']:.3f}"
    if s['p'] < 0.1: p_str += " †"
    dom_label = f"[{dom}]"
    dom_c = '#2166AC' if dom=='VDR-domain' else '#D6604D'
    n_total = (gene_dis[gene_dis['disease']==dis]).shape[0]
    ax_forest.text(27, yp,
                   f"OR={s['or_v']:.1f}  {p_str}",
                   va='center', fontsize=9.5, color='black')
    ax_forest.text(-1.5, yp,
                   f"n={n_total}  {dom_label}",
                   va='center', ha='right', fontsize=9,
                   color=dom_c, fontweight='bold')

ax_forest.axvline(1, color='black', lw=1.2, ls='--', alpha=0.6)
ax_forest.set_yticks(y_positions)
ax_forest.set_yticklabels(y_labels, fontsize=11)
ax_forest.set_xlabel('Odds Ratio (VDR-dominant target: Approved vs Failed)', fontsize=10)
ax_forest.set_xlim(-2, 35)
ax_forest.set_ylim(-0.6, 2.6)
ax_forest.set_title(
    'Disease Etiology Determines VDR-domain Advantage\n'
    '◆ = p<0.10 trend, VDR-domain diseases (CD, Psoriasis) show OR>1; GR-domain disease (UC) shows OR≈1',
    fontsize=10, fontweight='bold')
ax_forest.spines['top'].set_visible(False)
ax_forest.spines['right'].set_visible(False)
ax_forest.spines['left'].set_visible(False)

# 全体タイトル
fig.suptitle(
    "VDR-domain vs GR-domain Disease: Crohn's Disease, Ulcerative Colitis, and Psoriasis\n"
    "VDR-dominant drug targets succeed selectively in VDR-domain diseases where steroids underperform",
    fontsize=13, fontweight='bold', y=1.01)

out = RESULTS + 'fig_cd_uc_psoriasis'
fig.savefig(out+'.pdf', dpi=300, bbox_inches='tight', facecolor='white')
fig.savefig(out+'.png', dpi=200, bbox_inches='tight', facecolor='white')
print(f"\nSaved: {out}.pdf / .png")
print("DONE")
