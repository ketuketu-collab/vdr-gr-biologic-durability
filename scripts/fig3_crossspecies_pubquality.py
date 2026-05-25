#!/usr/bin/env python3
"""Publication-quality Figure 3: Cross-species VDR/GR conservation validation"""

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.gridspec import GridSpec
from scipy import stats
import openpyxl
import warnings
warnings.filterwarnings('ignore')

# --- Load tableS4 approval status ---
wb = openpyxl.load_workbook('/Volumes/M4_SSD/projects/tlr_chipseq/results/tableS4_phase1_vdr_gr_scores.xlsx')
ws = wb['全体 VDR降順']

approval = {}
for row in ws.iter_rows(min_row=2, values_only=True):
    if len(row) < 8:
        continue
    gene, result = row[1], row[7]
    if gene and result:
        approval[str(gene).strip()] = str(result).strip()

# --- Load cross-species data ---
cs = pd.read_csv('/Volumes/M4_SSD/projects/tlr_chipseq/results/crossspecies_validation.csv')

# log1p for display
cs['log_vdr_h'] = np.log1p(cs['VDR_human'])
cs['log_vdr_m'] = np.log1p(cs['VDR_mouse'])
eps = 1.0
cs['ratio_human'] = np.log2((cs['VDR_human'] + eps) / (cs['GR_human'] + eps))
cs['ratio_mouse'] = np.log2((cs['VDR_mouse'] + eps) / (cs['GR_mouse'] + eps))

# Map approval status
def get_status(gene):
    s = approval.get(str(gene), '')
    if '承認済' in s:
        return 'approved'
    elif '失敗' in s or '中止' in s:
        return 'failed'
    elif 'Phase' in s or '対照' in s or '研究段階' in s or '未達' in s:
        return 'ongoing'
    return 'other'

cs['status'] = cs['human_gene'].apply(get_status)
print(cs['status'].value_counts())

# Spearman correlations
r_vdr, p_vdr = stats.spearmanr(cs['VDR_human'], cs['VDR_mouse'])
r_gr,  p_gr  = stats.spearmanr(cs['GR_human'],  cs['GR_mouse'])
r_ratio, p_ratio = stats.spearmanr(cs['ratio_human'], cs['ratio_mouse'])
print(f"VDR r={r_vdr:.3f}, p={p_vdr:.2e}")
print(f"GR  r={r_gr:.3f},  p={p_gr:.2e}")

# --- Figure setup ---
fig = plt.figure(figsize=(18, 7))
fig.patch.set_facecolor('white')
gs = GridSpec(1, 3, figure=fig, wspace=0.38)

# Color palette
col_approved = '#2E7D32'
col_failed   = '#C62828'
col_ongoing  = '#E65100'
col_other    = '#BDBDBD'

color_map  = {'approved': col_approved, 'failed': col_failed, 'ongoing': col_ongoing, 'other': col_other}
size_map   = {'approved': 70, 'failed': 70, 'ongoing': 40, 'other': 25}
zorder_map = {'approved': 5, 'failed': 5, 'ongoing': 3, 'other': 2}

# ===== Panel A: VDR score cross-species scatter =====
ax_a = fig.add_subplot(gs[0])

for status, grp in cs.groupby('status'):
    ax_a.scatter(grp['log_vdr_h'], grp['log_vdr_m'],
                 c=color_map[status], s=size_map[status],
                 alpha=0.75, edgecolors='white', linewidths=0.4,
                 zorder=zorder_map[status])

# Regression line
slope, intercept, _, _, _ = stats.linregress(cs['log_vdr_h'], cs['log_vdr_m'])
x_line = np.linspace(cs['log_vdr_h'].min(), cs['log_vdr_h'].max(), 100)
ax_a.plot(x_line, slope * x_line + intercept, color='#1565C0', lw=1.5, ls='--', alpha=0.7, zorder=3)

# Stat box
ax_a.text(0.05, 0.97,
          f"Spearman ρ = {r_vdr:.3f}\np = {p_vdr:.1e}\nn = {len(cs)}",
          transform=ax_a.transAxes, fontsize=9.5, va='top',
          color='#1565C0', fontweight='bold',
          bbox=dict(boxstyle='round,pad=0.35', facecolor='#E3F2FD', alpha=0.9, edgecolor='#1565C0', lw=1.2))

# Label key genes
label_genes = {'IL4R','IL13','CCR3','IL23A','TNF','VEGFA','TLR10','IL17A','TSLP',
               'LAG3','IL1B','IL5','IL33','ITGB7','CD274'}
texts = []
try:
    from adjustText import adjust_text
    for _, row in cs.iterrows():
        if row['human_gene'] in label_genes:
            t = ax_a.annotate(row['human_gene'],
                              xy=(row['log_vdr_h'], row['log_vdr_m']),
                              fontsize=7.5, color='#212121',
                              xytext=(5, 4), textcoords='offset points')
            texts.append(t)
    adjust_text(texts, ax=ax_a, arrowprops=dict(arrowstyle='-', color='#9E9E9E', lw=0.5))
except ImportError:
    for _, row in cs.iterrows():
        if row['human_gene'] in label_genes:
            ax_a.annotate(row['human_gene'],
                          xy=(row['log_vdr_h'], row['log_vdr_m']),
                          fontsize=7.5, color='#212121',
                          xytext=(5, 4), textcoords='offset points')

ax_a.set_xlabel('VDR ChIP-seq score — Human (hg38, log₁₊ₓ)', fontsize=10)
ax_a.set_ylabel('VDR ChIP-seq score — Mouse (mm10, log₁₊ₓ)', fontsize=10)
ax_a.set_title('A   VDR binding conserved across species', fontsize=11, fontweight='bold', loc='left')

legend_elements = [
    mpatches.Patch(facecolor=col_approved, label='Approved biologic target'),
    mpatches.Patch(facecolor=col_failed,   label='Phase II/III failure'),
    mpatches.Patch(facecolor=col_ongoing,  label='In clinical development'),
    mpatches.Patch(facecolor=col_other,    label='Other targets'),
]
ax_a.legend(handles=legend_elements, fontsize=8, loc='upper left', framealpha=0.9)
ax_a.spines['top'].set_visible(False)
ax_a.spines['right'].set_visible(False)

# ===== Panel B: Asymmetric conservation bar chart =====
ax_b = fig.add_subplot(gs[1])

bars_r = [r_vdr, r_gr]
bar_labels = ['VDR\n(hg38 vs mm10)', 'GR/NR3C1\n(hg38 vs mm10)']
bar_colors = ['#1565C0', '#E53935']
p_vals = [p_vdr, p_gr]

bars = ax_b.bar(bar_labels, bars_r, color=bar_colors, width=0.45, alpha=0.85,
                edgecolor='white', linewidth=1.2, zorder=3)
ax_b.axhline(0, color='#424242', lw=0.8, zorder=2)

for i, (r_val, p_val) in enumerate(zip(bars_r, p_vals)):
    if p_val < 0.001:
        sig = '***'
    elif p_val < 0.01:
        sig = '**'
    elif p_val < 0.05:
        sig = '*'
    else:
        sig = 'n.s.'
    offset = 0.015 if r_val >= 0 else -0.04
    ax_b.text(i, r_val + offset, sig, ha='center', va='bottom' if r_val >= 0 else 'top',
              fontsize=12, fontweight='bold', color='#212121')
    ax_b.text(i, r_val / 2, f'ρ = {r_val:.3f}', ha='center', va='center',
              fontsize=9.5, color='white', fontweight='bold')

# ENCODE A549 inset
ax_b_inset = ax_b.inset_axes([0.58, 0.52, 0.38, 0.40])
r_encode = 0.202; p_encode = 3.73e-04
ax_b_inset.bar(['GR\nA549'], [r_encode], color='#EF6C00', alpha=0.85, width=0.5, edgecolor='white')
ax_b_inset.axhline(0, color='#424242', lw=0.6)
ax_b_inset.set_ylim(-0.05, 0.45)
ax_b_inset.text(0, r_encode + 0.02, '***', ha='center', fontsize=9, fontweight='bold')
ax_b_inset.text(0, r_encode / 2, f'{r_encode:.2f}', ha='center', va='center',
                fontsize=7.5, color='white', fontweight='bold')
ax_b_inset.set_title('ENCODE\nA549', fontsize=7.5)
ax_b_inset.tick_params(labelsize=6)
ax_b_inset.set_ylabel('ρ', fontsize=7)
ax_b_inset.spines['top'].set_visible(False)
ax_b_inset.spines['right'].set_visible(False)

ax_b.set_ylabel("Spearman's ρ (human vs. mouse)", fontsize=10)
ax_b.set_ylim(-0.05, 0.62)
ax_b.set_title('B   Asymmetric evolutionary conservation', fontsize=11, fontweight='bold', loc='left')
ax_b.text(0.5, 0.98,
          f'n = {len(cs)} orthologous genes\n(294/308 mapped)',
          transform=ax_b.transAxes, ha='center', va='top', fontsize=8.5,
          bbox=dict(boxstyle='round,pad=0.3', facecolor='#F5F5F5', edgecolor='#BDBDBD'))
ax_b.spines['top'].set_visible(False)
ax_b.spines['right'].set_visible(False)

# ===== Panel C: Evolutionary timeline concept =====
ax_c = fig.add_subplot(gs[2])
ax_c.set_xlim(0, 100)
ax_c.set_ylim(0, 100)
ax_c.axis('off')
ax_c.set_title('C   Evolutionary model', fontsize=11, fontweight='bold', loc='left')

# Timeline arrow
ax_c.annotate('', xy=(95, 52), xytext=(5, 52),
              arrowprops=dict(arrowstyle='->', color='#424242', lw=2))
ax_c.text(50, 45, 'Evolutionary time (~80 million years)', ha='center', fontsize=9, color='#616161')

# Species markers
for x, label in [(10, 'Common\nAncestor'), (50, 'Mus\nmusculus'), (88, 'Homo\nsapiens')]:
    ax_c.plot(x, 52, 'o', ms=11, color='#424242', zorder=5)
    ax_c.text(x, 61, label, ha='center', fontsize=8, color='#212121')

# VDR conservation track
vdr_y = 38
ax_c.annotate('', xy=(88, vdr_y), xytext=(10, vdr_y),
              arrowprops=dict(arrowstyle='->', color='#1565C0', lw=3.5))
ax_c.text(49, vdr_y + 5, 'VDR binding sites', ha='center', fontsize=9.5,
          color='#1565C0', fontweight='bold')
ax_c.text(49, vdr_y - 2,
          f'CONSERVED   ρ = {r_vdr:.3f},  p = {p_vdr:.0e}',
          ha='center', fontsize=8.5, color='#1565C0',
          bbox=dict(boxstyle='round,pad=0.3', facecolor='#E3F2FD', edgecolor='#1565C0', lw=1.3))

# GR divergence track
gr_y = 19
ax_c.plot([10, 88], [gr_y, gr_y], lw=2.5, color='#E53935', ls='--', alpha=0.5)
np.random.seed(42)
xs = np.linspace(10, 88, 50)
ys = gr_y + np.random.randn(50) * 2.0
ax_c.plot(xs, ys, lw=2, color='#E53935', alpha=0.65)
ax_c.text(49, gr_y + 6, 'GR binding sites', ha='center', fontsize=9.5,
          color='#C62828', fontweight='bold')
ax_c.text(49, gr_y - 2,
          f'DIVERGED   ρ = {r_gr:.3f},  n.s.',
          ha='center', fontsize=8.5, color='#C62828',
          bbox=dict(boxstyle='round,pad=0.3', facecolor='#FFEBEE', edgecolor='#C62828', lw=1.3))

# Interpretation
ax_c.text(50, 7,
          '"VDR = homeostatic chronic regulation\n  GR = context-dependent acute response"',
          ha='center', va='bottom', fontsize=8.5, style='italic', color='#212121',
          bbox=dict(boxstyle='round,pad=0.45', facecolor='#FFFDE7', edgecolor='#F9A825', lw=1.3))

# --- Save ---
out_base = '/Volumes/M4_SSD/projects/tlr_chipseq/results/fig3_crossspecies_validation_pubquality'
plt.savefig(out_base + '.pdf', dpi=300, bbox_inches='tight', facecolor='white')
plt.savefig(out_base + '.png', dpi=300, bbox_inches='tight', facecolor='white')
print(f"\nSaved: {out_base}.pdf / .png")
