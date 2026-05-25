"""
THP-1 vs Organoid VDR comparison figure
Shows that VDR binding at inflammatory disease target loci is immune cell-specific
"""

import pandas as pd
import numpy as np
import matplotlib
matplotlib.rcParams['font.family'] = 'DejaVu Sans'
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.patches as mpatches
from scipy import stats

# ── Load data ─────────────────────────────────────────────────────────────────

thp1 = pd.read_csv('/Volumes/M4_SSD/projects/tlr_chipseq/results/thp1_ls180_scores.csv')
org  = pd.read_csv('/Volumes/M4_SSD/projects/tlr_chipseq/results/organoid_vdr_scores.csv')

# Merge on gene
df = pd.merge(thp1[['gene','category','thp1_vdr','thp1_vdr_n','thp1_pattern']],
              org[['gene','organoid_vitD','organoid_veh','vitD_enrichment']],
              on='gene', how='inner')

print(f"Common genes: {len(df)}")

# ── Define groups ─────────────────────────────────────────────────────────────

# THP-1 pattern: VDR vs GR (already classified)
# For coloring: use category
cat_color = {
    'IBD/UC/CD':    '#2196F3',
    'Asthma/Atopy': '#FF9800',
    'RA':           '#9C27B0',
    'SLE/Renal':    '#E91E63',
    'Cancer/irAE':  '#4CAF50',
    'TLR/innate':   '#795548',
}
# Map Japanese → English categories
cat_map = {
    '喘息/アトピー': 'Asthma/Atopy',
    'SLE/腎疾患':   'SLE/Renal',
    'がん免疫/irAE':'Cancer/irAE',
}

df['category'] = df['category'].replace(cat_map)
df['color'] = df['category'].map(cat_color).fillna('#9E9E9E')

# Classify by THP-1 VDR score
HIGH_VDR = 50   # score threshold
df['thp1_class'] = np.where(df['thp1_vdr'] >= HIGH_VDR, 'VDR-high (≥50)', 'VDR-low (<50)')

# ── Figure layout ──────────────────────────────────────────────────────────────

fig = plt.figure(figsize=(16, 12))
gs = gridspec.GridSpec(2, 2, figure=fig, hspace=0.42, wspace=0.38)

ax_scatter  = fig.add_subplot(gs[0, 0])
ax_bar      = fig.add_subplot(gs[0, 1])
ax_org_rank = fig.add_subplot(gs[1, 0])
ax_heatmap  = fig.add_subplot(gs[1, 1])

# ── Panel A: scatter THP-1 VDR vs Organoid VDR ───────────────────────────────

ax = ax_scatter
ax.scatter(df['thp1_vdr'], df['organoid_vitD'],
           c=df['color'], s=60, alpha=0.8, edgecolors='none', zorder=3)

# Annotation for notable genes
label_genes = {'IL23A','TNF','IL6R','NOD2','OSM','CD274','LAG3','CCR9','TLR10','NLRP3'}
for _, row in df.iterrows():
    if row['gene'] in label_genes:
        ax.annotate(row['gene'],
                    xy=(row['thp1_vdr'], row['organoid_vitD']),
                    xytext=(5, 3), textcoords='offset points',
                    fontsize=7, color='#333333')

ax.axvline(HIGH_VDR, color='gray', ls='--', lw=0.8, alpha=0.6)
ax.set_xlabel('THP-1 VDR MACS2 score (immune monocyte)', fontsize=10)
ax.set_ylabel('Organoid VDR signal\n(colonic epithelium, +VitD)', fontsize=10)
ax.set_title('A  THP-1 vs Colonic Organoid\nVDR binding at target gene loci', fontsize=11, fontweight='bold')

# correlation
r, p = stats.pearsonr(df['thp1_vdr'], df['organoid_vitD'])
ax.text(0.97, 0.97, f'r = {r:.2f}\np = {p:.3f}',
        transform=ax.transAxes, ha='right', va='top', fontsize=9,
        bbox=dict(fc='white', ec='gray', alpha=0.7, boxstyle='round,pad=0.3'))

# legend
handles = [mpatches.Patch(facecolor=c, label=cat) for cat, c in cat_color.items()]
ax.legend(handles=handles, fontsize=7, loc='upper right',
          bbox_to_anchor=(1.0, 0.82), framealpha=0.7)

# ── Panel B: Mean organoid signal by THP-1 class ─────────────────────────────

ax = ax_bar
grp = df.groupby('thp1_class')['organoid_vitD'].agg(['mean','sem','count']).reset_index()
grp = grp.sort_values('mean', ascending=False)

colors_bar = ['#1565C0' if '≥50' in g else '#90CAF9' for g in grp['thp1_class']]
bars = ax.bar(grp['thp1_class'], grp['mean'],
              yerr=grp['sem'], capsize=5,
              color=colors_bar, edgecolor='none', alpha=0.85)

for bar, (_, row) in zip(bars, grp.iterrows()):
    ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + row['sem'] + 0.01,
            f'n={int(row["count"])}', ha='center', va='bottom', fontsize=9)

ax.set_ylabel('Mean organoid VDR signal (+VitD)', fontsize=10)
ax.set_title('B  Organoid VDR signal\nby THP-1 VDR binding class', fontsize=11, fontweight='bold')
ax.set_ylim(0, ax.get_ylim()[1] * 1.2)
ax.axhline(1.0, color='gray', ls=':', lw=0.8)
ax.text(0.97, 0.97, 'Dashed = background level',
        transform=ax.transAxes, ha='right', va='top', fontsize=8, color='gray')

# t-test
if len(grp) == 2:
    g1 = df[df['thp1_class'].str.contains('≥50')]['organoid_vitD']
    g2 = df[df['thp1_class'].str.contains('<50')]['organoid_vitD']
    t, p_t = stats.ttest_ind(g1, g2)
    ax.text(0.5, 0.05, f't-test p = {p_t:.3f}',
            transform=ax.transAxes, ha='center', fontsize=9, color='gray')

# ── Panel C: Organoid VitD enrichment ranked ─────────────────────────────────

ax = ax_org_rank
df_sorted = df.sort_values('vitD_enrichment', ascending=False)

bar_colors = ['#E53935' if v > 0 else '#1E88E5' for v in df_sorted['vitD_enrichment']]
ax.barh(range(len(df_sorted)), df_sorted['vitD_enrichment'],
        color=bar_colors, edgecolor='none', alpha=0.8)
ax.set_yticks(range(len(df_sorted)))
ax.set_yticklabels(df_sorted['gene'], fontsize=7)
ax.axvline(0, color='black', lw=0.8)
ax.set_xlabel('VitD enrichment\n(organoid vitD − vehicle)', fontsize=10)
ax.set_title('C  VitD-induced VDR enrichment\nin colonic organoids', fontsize=11, fontweight='bold')

# mark THP-1 VDR-high genes
for i, (_, row) in enumerate(df_sorted.iterrows()):
    if row['thp1_vdr'] >= HIGH_VDR:
        ax.text(ax.get_xlim()[1] * 0.95, i, '*', ha='right', va='center',
                fontsize=10, color='#1565C0', fontweight='bold')

ax.text(0.98, 0.02, '* = THP-1 VDR-high (score ≥50)',
        transform=ax.transAxes, ha='right', va='bottom',
        fontsize=7.5, color='#1565C0')

# ── Panel D: Conceptual summary heatmap ───────────────────────────────────────

ax = ax_heatmap

# Select top genes (THP-1 VDR high) and show THP-1 vs organoid side by side
top_genes = df.nlargest(20, 'thp1_vdr')['gene'].tolist()
df_top = df[df['gene'].isin(top_genes)].sort_values('thp1_vdr', ascending=False)

# Normalize to 0-1 for heatmap display
thp1_norm = df_top['thp1_vdr'] / df_top['thp1_vdr'].max()
org_norm  = (df_top['organoid_vitD'] - df_top['organoid_vitD'].min()) / \
            (df_top['organoid_vitD'].max() - df_top['organoid_vitD'].min())

data_matrix = np.column_stack([thp1_norm.values, org_norm.values])

cmap = LinearSegmentedColormap.from_list('blue_red', ['#EEF5FF', '#1565C0'])
im = ax.imshow(data_matrix, aspect='auto', cmap=cmap, vmin=0, vmax=1)

ax.set_yticks(range(len(df_top)))
ax.set_yticklabels(df_top['gene'], fontsize=8)
ax.set_xticks([0, 1])
ax.set_xticklabels(['THP-1\n(monocyte)', 'Organoid\n(colonocyte)'], fontsize=10)
ax.set_title('D  Top THP-1 VDR targets:\nimmune vs epithelial context', fontsize=11, fontweight='bold')

cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
cbar.set_label('Normalized VDR binding', fontsize=8)

# ── Overall title ─────────────────────────────────────────────────────────────

fig.suptitle('VDR binding at inflammatory disease target loci is immune cell-specific\n'
             'THP-1 monocytes vs human colonic organoids (GSE206176)',
             fontsize=13, fontweight='bold', y=1.01)

# ── Save ─────────────────────────────────────────────────────────────────────

out_png = '/Volumes/M4_SSD/projects/tlr_chipseq/results/fig_thp1_vs_organoid.png'
out_pdf = '/Volumes/M4_SSD/projects/tlr_chipseq/results/fig_thp1_vs_organoid.pdf'
plt.savefig(out_png, dpi=200, bbox_inches='tight')
plt.savefig(out_pdf, bbox_inches='tight')
print(f"Saved: {out_png}")
print(f"Saved: {out_pdf}")

# ── Print summary stats ───────────────────────────────────────────────────────

print("\n=== SUMMARY ===")
print(f"Common genes analyzed: {len(df)}")
print(f"\nTHP-1 VDR-high (score ≥{HIGH_VDR}): n={df['thp1_vdr_n_x'].sum() if 'thp1_vdr_n_x' in df else (df['thp1_vdr']>=HIGH_VDR).sum()}")
high = df[df['thp1_vdr'] >= HIGH_VDR]
low  = df[df['thp1_vdr'] <  HIGH_VDR]
print(f"  Organoid signal (mean±SEM): {high['organoid_vitD'].mean():.3f} ± {high['organoid_vitD'].sem():.3f}")
print(f"THP-1 VDR-low  (score <{HIGH_VDR}): n={(df['thp1_vdr']<HIGH_VDR).sum()}")
print(f"  Organoid signal (mean±SEM): {low['organoid_vitD'].mean():.3f} ± {low['organoid_vitD'].sem():.3f}")

print(f"\nPearson r (THP-1 vs Organoid): {r:.3f}, p={p:.4f}")

print(f"\nVitD enrichment in organoids:")
print(f"  Genes with positive enrichment: {(df['vitD_enrichment']>0).sum()}/{len(df)}")
print(f"  Median enrichment: {df['vitD_enrichment'].median():.3f}")

print("\nTop 5 organoid VitD-enriched genes:")
print(df.nlargest(5, 'vitD_enrichment')[['gene','thp1_vdr','organoid_vitD','vitD_enrichment']].to_string(index=False))
