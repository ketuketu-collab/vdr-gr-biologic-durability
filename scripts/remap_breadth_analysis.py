#!/usr/bin/env python3
"""
ReMap2022 per-gene breadth/reproducibility analysis
For each target gene TSS±5kb:
  - How many unique cell types show VDR binding?
  - How many unique cell types show GR binding?
  - Reproducibility score = n_experiments_with_binding / total_experiments_for_TF
"""

import gzip, re
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy import stats
import openpyxl
import warnings
warnings.filterwarnings('ignore')

DATA = '/Volumes/M4_SSD/projects/tlr10_chipseq_vdr_gr/data/'
VDR_BED  = DATA + 'remap2022_VDR_all_macs2_hg38.bed.gz'
GR_BED   = DATA + 'remap2022_NR3C1_all_macs2_hg38.bed.gz'
TARGETS  = '/Volumes/M4_SSD/projects/tlr_chipseq/results/steroid_targets_vdr_gr_scores.csv'
TSS_F    = '/Volumes/M4_SSD/projects/tlr_chipseq/data/gr_monocyte_primary/tss_hg19.tsv'
HALF     = 5000

# ---------- 1. TSS coordinates (hg38 via Ensembl REST API) ----------
print("Fetching TSS coordinates from Ensembl REST API...")
import urllib.request, urllib.parse, json, time

def fetch_tss_ensembl_batch(genes):
    url  = 'https://rest.ensembl.org/lookup/symbol/homo_sapiens'
    data = json.dumps({'symbols': genes}).encode()
    req  = urllib.request.Request(url, data=data,
           headers={'Content-Type': 'application/json', 'Accept': 'application/json'})
    try:
        with urllib.request.urlopen(req, timeout=60) as r:
            d = json.loads(r.read())
        tss = {}
        for gene, info in d.items():
            if not info: continue
            chrom = info.get('seq_region_name', '')
            if re.match(r'^\d+$|^X$|^Y$', str(chrom)):
                chrom = 'chr' + str(chrom)
                strand = info.get('strand', 1)
                tss_pos = info['start'] if strand == 1 else info['end']
                tss[gene] = (chrom, tss_pos)
        return tss
    except Exception as e:
        print(f"  Ensembl REST error: {e}")
        return {}

targets_df = pd.read_csv(TARGETS)
genes = targets_df['gene'].tolist()

tss_map = {}
for i in range(0, len(genes), 50):
    chunk = genes[i:i+50]
    t = fetch_tss_ensembl_batch(chunk)
    tss_map.update(t)
    print(f"  {min(i+50, len(genes))}/{len(genes)} genes processed ({len(tss_map)} with TSS)")
    time.sleep(0.5)

print(f"TSS found for {len(tss_map)}/{len(genes)} genes")

# ---------- 2. Build interval index from BED ----------
def build_index(bed_gz):
    idx = {}  # chrom -> list of (start, end, celltype, gse_id)
    with gzip.open(bed_gz, 'rt') as f:
        for line in f:
            p = line.rstrip().split('\t')
            if len(p) < 4: continue
            chrom, start, end, name = p[0], int(p[1]), int(p[2]), p[3]
            parts = name.split('.')
            gse_id   = parts[0] if len(parts) > 0 else 'unknown'
            celltype = parts[2] if len(parts) > 2 else 'unknown'
            # Normalize cell type (strip condition suffixes for base cell type)
            base_cell = re.sub(r'_[A-Z0-9_]+$', '', celltype, flags=re.IGNORECASE)
            if chrom not in idx:
                idx[chrom] = []
            idx[chrom].append((start, end, celltype, base_cell, gse_id))
    # Sort by start
    for chrom in idx:
        idx[chrom].sort(key=lambda x: x[0])
    return idx

print("Building VDR index...")
vdr_idx = build_index(VDR_BED)
print("Building GR index...")
gr_idx  = build_index(GR_BED)

# Total experiments per TF
vdr_total_exp = len(set(
    item[4] for items in vdr_idx.values() for item in items
))
gr_total_exp = len(set(
    item[4] for items in gr_idx.values() for item in items
))
vdr_total_cell = len(set(
    item[2] for items in vdr_idx.values() for item in items
))
gr_total_cell = len(set(
    item[2] for items in gr_idx.values() for item in items
))
print(f"VDR: {vdr_total_exp} experiments, {vdr_total_cell} cell types")
print(f"GR:  {gr_total_exp} experiments, {gr_total_cell} cell types")

# ---------- 3. Per-gene overlap ----------
def query_tss(idx, chrom, tss_pos):
    if chrom not in idx: return [], [], []
    start = max(0, tss_pos - HALF)
    end   = tss_pos + HALF
    ovlp  = [x for x in idx[chrom] if x[0] < end and x[1] > start]
    celltypes = set(x[2] for x in ovlp)
    base_cells = set(x[3] for x in ovlp)
    gse_ids    = set(x[4] for x in ovlp)
    return celltypes, base_cells, gse_ids

results = []
for _, row in targets_df.iterrows():
    gene = row['gene']
    if gene not in tss_map:
        continue
    chrom, tss = tss_map[gene]

    vdr_ct, vdr_bc, vdr_gse = query_tss(vdr_idx, chrom, tss)
    gr_ct,  gr_bc,  gr_gse  = query_tss(gr_idx,  chrom, tss)

    results.append({
        'gene':            gene,
        'VDR_score':       row['VDR'],
        'GR_score':        row['GR'],
        'VDR_celltypes':   len(vdr_ct),
        'GR_celltypes':    len(gr_ct),
        'VDR_base_cells':  len(vdr_bc),
        'GR_base_cells':   len(gr_bc),
        'VDR_experiments': len(vdr_gse),
        'GR_experiments':  len(gr_gse),
        # Reproducibility: fraction of available experiments showing binding
        'VDR_repro':       len(vdr_gse) / vdr_total_exp,
        'GR_repro':        len(gr_gse)  / gr_total_exp,
    })

df = pd.DataFrame(results)
df['VDR_dominant'] = df['VDR_score'] > df['GR_score']
df['repro_ratio']  = np.log2((df['VDR_repro'] + 0.01) / (df['GR_repro'] + 0.01))

# ---------- 4. Load approval status ----------
wb = openpyxl.load_workbook('/Volumes/M4_SSD/projects/tlr_chipseq/results/tableS4_phase1_vdr_gr_scores.xlsx')
ws = wb['全体 VDR降順']
approval = {}
for row in ws.iter_rows(min_row=2, values_only=True):
    if len(row) >= 8 and row[1] and row[7]:
        approval[str(row[1]).strip()] = str(row[7]).strip()

def get_status(gene):
    s = approval.get(str(gene), '')
    if '承認済' in s: return 'approved'
    elif '失敗' in s or '中止' in s: return 'failed'
    elif 'Phase' in s: return 'ongoing'
    return 'other'

df['status'] = df['gene'].apply(get_status)

# Save
out_csv = '/Volumes/M4_SSD/projects/tlr_chipseq/results/remap_breadth_per_gene.csv'
df.to_csv(out_csv, index=False)
print(f"\nSaved: {out_csv}")
print(df[['gene','VDR_score','GR_score','VDR_celltypes','GR_celltypes',
          'VDR_repro','GR_repro','status']].head(20).to_string())

# ---------- 5. Group statistics ----------
approved = df[df['status'] == 'approved']
failed   = df[df['status'] == 'failed']

print(f"\n=== VDR reproducibility ===")
print(f"Approved  mean={approved['VDR_repro'].mean():.3f}  median={approved['VDR_repro'].median():.3f}")
print(f"Failed    mean={failed['VDR_repro'].mean():.3f}  median={failed['VDR_repro'].median():.3f}")
stat, p = stats.mannwhitneyu(approved['VDR_repro'], failed['VDR_repro'], alternative='greater')
print(f"Mann-Whitney U (approved > failed): p={p:.3e}")

print(f"\n=== GR reproducibility ===")
print(f"Approved  mean={approved['GR_repro'].mean():.3f}  median={approved['GR_repro'].median():.3f}")
print(f"Failed    mean={failed['GR_repro'].mean():.3f}  median={failed['GR_repro'].median():.3f}")
stat2, p2 = stats.mannwhitneyu(failed['GR_repro'], approved['GR_repro'], alternative='greater')
print(f"Mann-Whitney U (failed > approved): p={p2:.3e}")

# Spearman: VDR_score vs VDR_repro
r_vdr_s, p_vdr_s = stats.spearmanr(df['VDR_score'], df['VDR_repro'])
r_gr_s,  p_gr_s  = stats.spearmanr(df['GR_score'],  df['GR_repro'])
print(f"\nVDR score vs VDR reproducibility: r={r_vdr_s:.3f}, p={p_vdr_s:.2e}")
print(f"GR  score vs GR  reproducibility: r={r_gr_s:.3f},  p={p_gr_s:.2e}")

# ---------- 6. Figure ----------
fig, axes = plt.subplots(1, 3, figsize=(16, 6))
fig.patch.set_facecolor('white')

col_approved = '#2E7D32'
col_failed   = '#C62828'
col_ongoing  = '#E65100'
col_other    = '#BDBDBD'
cmap_s = {'approved': col_approved, 'failed': col_failed, 'ongoing': col_ongoing, 'other': col_other}

# Panel A: VDR reproducibility by status (boxplot + swarm)
ax = axes[0]
order = ['approved', 'failed', 'ongoing']
data_by_status = {s: df[df['status']==s]['VDR_repro'].values for s in order}
labels_n = [f"{s}\n(n={len(data_by_status[s])})" for s in order]
bp = ax.boxplot([data_by_status[s] for s in order], labels=labels_n,
                patch_artist=True, medianprops=dict(color='white', lw=2),
                whiskerprops=dict(lw=1.2), capprops=dict(lw=1.2),
                flierprops=dict(marker='', alpha=0))
for patch, s in zip(bp['boxes'], order):
    patch.set_facecolor(cmap_s[s])
    patch.set_alpha(0.7)
# Jitter
for i, s in enumerate(order):
    y = data_by_status[s]
    x = np.random.normal(i+1, 0.06, len(y))
    ax.scatter(x, y, c=cmap_s[s], s=30, alpha=0.7, zorder=3, edgecolors='white', lw=0.3)

# p-value bracket
y_max = df['VDR_repro'].max() * 1.05
ax.plot([1, 2], [y_max, y_max], 'k-', lw=1)
ax.text(1.5, y_max * 1.02, f'p={p:.3f}', ha='center', fontsize=9)
ax.set_ylabel('VDR binding reproducibility\n(fraction of 10 experiments)', fontsize=10)
ax.set_title('A  VDR reproducibility by drug outcome', fontsize=11, fontweight='bold', loc='left')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# Panel B: Scatter VDR repro vs GR repro per gene
ax = axes[1]
for status, grp in df.groupby('status'):
    sm = {'approved': 70, 'failed': 70, 'ongoing': 40, 'other': 20}
    zo = {'approved': 5, 'failed': 5, 'ongoing': 3, 'other': 2}
    ax.scatter(grp['VDR_repro'], grp['GR_repro'],
               c=cmap_s.get(status, col_other), s=sm.get(status, 20),
               alpha=0.75, edgecolors='white', lw=0.3, zorder=zo.get(status, 2))

# Diagonal
lim = max(df[['VDR_repro','GR_repro']].max())
ax.plot([0, lim], [0, lim], 'k--', lw=0.8, alpha=0.5)
ax.set_xlabel('VDR binding reproducibility', fontsize=10)
ax.set_ylabel('GR binding reproducibility', fontsize=10)
ax.set_title('B  Per-gene VDR vs GR reproducibility', fontsize=11, fontweight='bold', loc='left')

# Label key genes
key = {'IL4R','IL13','TNF','TSLP','IL17A','IL5','DUSP1','TNFAIP3','TLR2','CCR3','VEGFA'}
for _, row in df.iterrows():
    if row['gene'] in key:
        ax.annotate(row['gene'], (row['VDR_repro'], row['GR_repro']),
                    fontsize=7, xytext=(4,3), textcoords='offset points')
legend_elements = [
    mpatches.Patch(facecolor=col_approved, label='Approved biologic'),
    mpatches.Patch(facecolor=col_failed,   label='Phase II/III failure'),
    mpatches.Patch(facecolor=col_ongoing,  label='In development'),
]
ax.legend(handles=legend_elements, fontsize=8, loc='upper right')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# Panel C: Cell type breadth — VDR vs GR per gene (approved genes only)
ax = axes[2]
app_df = df[df['status']=='approved'].copy()
app_df = app_df.sort_values('VDR_celltypes', ascending=True)
y_pos = np.arange(len(app_df))
ax.barh(y_pos - 0.2, app_df['VDR_celltypes'], height=0.38,
        color=col_approved, alpha=0.8, label='VDR')
ax.barh(y_pos + 0.2, app_df['GR_celltypes'],  height=0.38,
        color='#E53935', alpha=0.8, label='GR')
ax.set_yticks(y_pos)
ax.set_yticklabels(app_df['gene'], fontsize=8)
ax.set_xlabel('Number of cell types with binding\n(TSS ± 5 kb)', fontsize=10)
ax.set_title('C  Cell type breadth (approved targets)', fontsize=11, fontweight='bold', loc='left')
ax.legend(fontsize=9)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.tight_layout()
out = '/Volumes/M4_SSD/projects/tlr_chipseq/results/fig_remap_breadth_analysis'
plt.savefig(out + '.pdf', dpi=300, bbox_inches='tight', facecolor='white')
plt.savefig(out + '.png', dpi=300, bbox_inches='tight', facecolor='white')
print(f"\nFigure saved: {out}.pdf / .png")
