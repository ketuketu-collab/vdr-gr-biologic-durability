#!/usr/bin/env python3
"""Publication-quality figure: ReMap2022 reproducibility analysis"""

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

df = pd.read_csv('/Volumes/M4_SSD/projects/tlr_chipseq/results/remap_breadth_per_gene.csv')

# Group definitions
df['dominant'] = np.where(df['VDR_score'] > df['GR_score'], 'VDR-dominant', 'GR-dominant')
vdr_dom = df[df['dominant'] == 'VDR-dominant']
gr_dom  = df[df['dominant'] == 'GR-dominant']
approved = df[df['status'] == 'approved']
failed   = df[df['status'] == 'failed']

# Stats
_, p_vdr = stats.mannwhitneyu(vdr_dom['VDR_repro'], gr_dom['VDR_repro'], alternative='greater')
_, p_gr  = stats.mannwhitneyu(gr_dom['GR_repro'],   vdr_dom['GR_repro'], alternative='greater')
r_vdr_s, p_vdr_s = stats.spearmanr(df['VDR_score'], df['VDR_repro'])
r_gr_s,  p_gr_s  = stats.spearmanr(df['GR_score'],  df['GR_repro'])

print(f"VDR: repro dom vs non-dom p={p_vdr:.3e}")
print(f"GR:  repro dom vs non-dom p={p_gr:.3e}")
print(f"VDR score-repro r={r_vdr_s:.3f} p={p_vdr_s:.1e}")
print(f"GR  score-repro r={r_gr_s:.3f}  p={p_gr_s:.1e}")

# Colors
col_vdr      = '#1565C0'
col_gr       = '#E53935'
col_approved = '#2E7D32'
col_failed   = '#C62828'
col_ongoing  = '#E65100'
col_other    = '#BDBDBD'
cmap_status  = {'approved': col_approved, 'failed': col_failed,
                'ongoing': col_ongoing, 'other': col_other}

fig = plt.figure(figsize=(18, 6))
fig.patch.set_facecolor('white')
gs = GridSpec(1, 3, figure=fig, wspace=0.38)

# ===== Panel A: VDR repro and GR repro by dominance group =====
ax_a = fig.add_subplot(gs[0])

groups = ['VDR-dominant\n(n=%d)' % len(vdr_dom),
          'GR-dominant\n(n=%d)'  % len(gr_dom)]
vdr_data = [vdr_dom['VDR_repro'].values, gr_dom['VDR_repro'].values]
gr_data  = [vdr_dom['GR_repro'].values,  gr_dom['GR_repro'].values]

x = np.array([1, 2])
w = 0.3

# Box plots side by side
bp1 = ax_a.boxplot(vdr_data, positions=x - w/2, widths=w*0.85,
                   patch_artist=True, medianprops=dict(color='white', lw=2),
                   whiskerprops=dict(lw=1.2), capprops=dict(lw=1.2),
                   flierprops=dict(marker='', alpha=0))
bp2 = ax_a.boxplot(gr_data,  positions=x + w/2, widths=w*0.85,
                   patch_artist=True, medianprops=dict(color='white', lw=2),
                   whiskerprops=dict(lw=1.2), capprops=dict(lw=1.2),
                   flierprops=dict(marker='', alpha=0))
for patch in bp1['boxes']: patch.set_facecolor(col_vdr); patch.set_alpha(0.8)
for patch in bp2['boxes']: patch.set_facecolor(col_gr);  patch.set_alpha(0.8)

# Jitter
for i, (vdat, gdat) in enumerate(zip(vdr_data, gr_data)):
    xv = np.random.normal(x[i] - w/2, 0.05, len(vdat))
    xg = np.random.normal(x[i] + w/2, 0.05, len(gdat))
    ax_a.scatter(xv, vdat, c=col_vdr, s=18, alpha=0.5, zorder=3)
    ax_a.scatter(xg, gdat, c=col_gr,  s=18, alpha=0.5, zorder=3)

# Significance brackets
y_max = max(df['VDR_repro'].max(), df['GR_repro'].max()) + 0.02
def sig_star(p):
    if p < 0.001: return '***'
    elif p < 0.01: return '**'
    elif p < 0.05: return '*'
    return 'n.s.'

# VDR repro bracket (above left box pair → left box pair comparison)
y1 = 0.65
ax_a.plot([x[0]-w/2, x[1]-w/2], [y1, y1], 'k-', lw=1)
ax_a.text((x[0]+x[1])/2 - w/2, y1+0.01,
          f'p={p_vdr:.3f} {sig_star(p_vdr)}', ha='center', fontsize=8.5,
          color=col_vdr, fontweight='bold')

# GR repro bracket
y2 = 0.75
ax_a.plot([x[0]+w/2, x[1]+w/2], [y2, y2], 'k-', lw=1)
ax_a.text((x[0]+x[1])/2 + w/2, y2+0.01,
          f'p={p_gr:.4f} {sig_star(p_gr)}', ha='center', fontsize=8.5,
          color=col_gr, fontweight='bold')

ax_a.set_xticks(x)
ax_a.set_xticklabels(groups, fontsize=9)
ax_a.set_ylabel('Binding reproducibility\n(fraction of independent experiments)', fontsize=10)
ax_a.set_ylim(0, 0.85)
ax_a.set_title('A   Score reflects reproducibility across experiments', fontsize=11, fontweight='bold', loc='left')

legend_el = [
    mpatches.Patch(facecolor=col_vdr, alpha=0.8, label='VDR binding reproducibility'),
    mpatches.Patch(facecolor=col_gr,  alpha=0.8, label='GR binding reproducibility'),
]
ax_a.legend(handles=legend_el, fontsize=8.5, loc='upper right')
ax_a.spines['top'].set_visible(False)
ax_a.spines['right'].set_visible(False)

# ===== Panel B: VDR score vs VDR repro (per gene, approval colored) =====
ax_b = fig.add_subplot(gs[1])

size_map   = {'approved': 70, 'failed': 70, 'ongoing': 40, 'other': 18}
zorder_map = {'approved': 5, 'failed': 5, 'ongoing': 3, 'other': 2}

for status, grp in df.groupby('status'):
    ax_b.scatter(np.log1p(grp['VDR_score']), grp['VDR_repro'],
                 c=cmap_status.get(status, col_other),
                 s=size_map.get(status, 18),
                 alpha=0.75, edgecolors='white', lw=0.3,
                 zorder=zorder_map.get(status, 2))

# Regression line
xs = np.log1p(df['VDR_score'].values)
slope, intercept, _, _, _ = stats.linregress(xs, df['VDR_repro'].values)
x_line = np.linspace(xs.min(), xs.max(), 100)
ax_b.plot(x_line, slope*x_line + intercept, color=col_vdr, lw=1.5, ls='--', alpha=0.7)

ax_b.text(0.05, 0.97,
          f'Spearman ρ = {r_vdr_s:.3f}\np = {p_vdr_s:.1e}\nn = {len(df)}',
          transform=ax_b.transAxes, fontsize=9.5, va='top', color=col_vdr, fontweight='bold',
          bbox=dict(boxstyle='round,pad=0.35', facecolor='#E3F2FD', alpha=0.9, edgecolor=col_vdr))

# Label key genes
key_genes = {'IL4R','CCR3','TNF','IL17A','TSLP','IL5','VEGFA','LAG3','CD274',
             'DUSP1','TNFAIP3','TLR2','MAPK14','MMP9'}
for _, row in df.iterrows():
    if row['gene'] in key_genes:
        ax_b.annotate(row['gene'], (np.log1p(row['VDR_score']), row['VDR_repro']),
                      fontsize=7, xytext=(4,3), textcoords='offset points', color='#212121')

ax_b.set_xlabel('VDR ChIP-seq score (log₁₊ₓ)', fontsize=10)
ax_b.set_ylabel('VDR binding reproducibility', fontsize=10)
ax_b.set_title('B   VDR score validates against reproducibility', fontsize=11, fontweight='bold', loc='left')

legend_el2 = [
    mpatches.Patch(facecolor=col_approved, label='Approved biologic'),
    mpatches.Patch(facecolor=col_failed,   label='Phase II/III failure'),
    mpatches.Patch(facecolor=col_ongoing,  label='In development'),
    mpatches.Patch(facecolor=col_other,    label='Other targets'),
]
ax_b.legend(handles=legend_el2, fontsize=8, loc='lower right')
ax_b.spines['top'].set_visible(False)
ax_b.spines['right'].set_visible(False)

# ===== Panel C: Per-gene heatmap of VDR vs GR repro for key genes =====
ax_c = fig.add_subplot(gs[2])

# Select approved + failed genes with non-zero scores
key_df = df[df['status'].isin(['approved','failed'])].copy()
key_df = key_df[key_df[['VDR_score','GR_score']].max(axis=1) > 0].copy()
key_df['repro_diff'] = key_df['VDR_repro'] - key_df['GR_repro']
key_df = key_df.sort_values(['status', 'repro_diff'], ascending=[True, False])

y_pos = np.arange(len(key_df))
colors = [col_approved if s == 'approved' else col_failed for s in key_df['status']]
bar_vals = key_df['repro_diff'].values
bar_cols = [col_vdr if v > 0 else col_gr for v in bar_vals]

ax_c.barh(y_pos, bar_vals, color=bar_cols, alpha=0.8, height=0.7, edgecolor='white', lw=0.3)
ax_c.axvline(0, color='#424242', lw=0.8)

# Gene labels colored by approval
for i, (_, row) in enumerate(key_df.iterrows()):
    col = col_approved if row['status'] == 'approved' else col_failed
    ax_c.text(-0.005, i, row['gene'], ha='right', va='center', fontsize=7.5,
              color=col, fontweight='bold' if row['status'] == 'failed' else 'normal')

ax_c.set_yticks([])
ax_c.set_xlabel('VDR repro − GR repro\n(positive = VDR more reproducible)', fontsize=10)
ax_c.set_title('C   Reproducibility balance (approved vs failed)', fontsize=11, fontweight='bold', loc='left')
ax_c.set_xlim(-0.45, 0.55)

# Annotations
ax_c.text(0.25, len(key_df)-1.5, '← VDR-favored', fontsize=8, color=col_vdr, style='italic', ha='center')
ax_c.text(-0.22, len(key_df)-1.5, 'GR-favored →', fontsize=8, color=col_gr, style='italic', ha='center')

legend_el3 = [
    mpatches.Patch(facecolor=col_approved, label='Approved biologic'),
    mpatches.Patch(facecolor=col_failed,   label='Phase II/III failure'),
]
ax_c.legend(handles=legend_el3, fontsize=8.5, loc='lower right')
ax_c.spines['top'].set_visible(False)
ax_c.spines['right'].set_visible(False)

plt.tight_layout()
out = '/Volumes/M4_SSD/projects/tlr_chipseq/results/fig_remap_breadth_pubquality'
plt.savefig(out + '.pdf', dpi=300, bbox_inches='tight', facecolor='white')
plt.savefig(out + '.png', dpi=300, bbox_inches='tight', facecolor='white')
print(f"Saved: {out}.pdf / .png")
