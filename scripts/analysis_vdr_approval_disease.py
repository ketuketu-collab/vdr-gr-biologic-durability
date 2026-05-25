#!/usr/bin/env python3
"""
① Fisher's exact test: VDR-dominant enrichment in approved drugs
② Disease mapping: VDR-dominant 81 genes via OpenTargets GraphQL API
"""

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.gridspec import GridSpec
import matplotlib.colors as mcolors
from scipy import stats
import urllib.request, json, time
import warnings
warnings.filterwarnings('ignore')

RESULTS = '/Volumes/M4_SSD/projects/tlr_chipseq/results/'
CSV     = RESULTS + 'remap_breadth_per_gene.csv'

df = pd.read_csv(CSV)
print(f"Loaded: {len(df)} genes")

# ── colours ────────────────────────────────────────────────────────────────
COL_VDR  = '#2166AC'
COL_GR   = '#D6604D'
COL_BOTH = '#762a83'
COL_APP  = '#1a9641'
COL_FAIL = '#d7191c'

# ══════════════════════════════════════════════════════════════════════════
# PART 1: Fisher's exact test
# ══════════════════════════════════════════════════════════════════════════
print("\n=== PART 1: Fisher's exact test ===")

# --- Test A: approved vs failed (strict) ---
sub = df[df['status'].isin(['approved','failed'])].copy()
ct_strict = pd.crosstab(sub['VDR_dominant'], sub['status'])
print("\nStrict (approved vs failed):")
print(ct_strict)

# VDR-dom: approved=a, failed=b; GR-dom: approved=c, failed=d
n_va = int(ct_strict.loc[True,  'approved'])   # VDR-dom approved
n_vf = int(ct_strict.loc[True,  'failed'])     # VDR-dom failed
n_ga = int(ct_strict.loc[False, 'approved'])   # GR-dom approved
n_gf = int(ct_strict.loc[False, 'failed'])     # GR-dom failed
table_strict = [[n_va, n_vf],[n_ga, n_gf]]
or_strict, p_strict = stats.fisher_exact(table_strict, alternative='greater')
print(f"OR={or_strict:.2f}  p={p_strict:.4f}")

# 95% CI via log-OR ± 1.96*SE
log_or = np.log(or_strict)
se = np.sqrt(1/n_va + 1/n_vf + 1/n_ga + 1/n_gf)
ci_lo = np.exp(log_or - 1.96*se)
ci_hi = np.exp(log_or + 1.96*se)
print(f"95%CI: [{ci_lo:.2f}, {ci_hi:.2f}]")

# --- Test B: approved vs failed+other (broad) ---
df['binary'] = df['status'].apply(lambda x: 'approved' if x=='approved' else 'non-approved')
ct_broad = pd.crosstab(df['VDR_dominant'], df['binary'])
print("\nBroad (approved vs non-approved):")
print(ct_broad)
a2 = int(ct_broad.loc[True,  'approved'])
b2 = int(ct_broad.loc[True,  'non-approved'])
c2 = int(ct_broad.loc[False, 'approved'])
d2 = int(ct_broad.loc[False, 'non-approved'])
or_broad, p_broad = stats.fisher_exact([[a2,b2],[c2,d2]], alternative='greater')
log_or2 = np.log(or_broad)
se2 = np.sqrt(1/a2 + 1/b2 + 1/c2 + 1/d2)
ci_lo2 = np.exp(log_or2 - 1.96*se2)
ci_hi2 = np.exp(log_or2 + 1.96*se2)
print(f"OR={or_broad:.2f}  p={p_broad:.4f}  95%CI: [{ci_lo2:.2f}, {ci_hi2:.2f}]")

# ══════════════════════════════════════════════════════════════════════════
# PART 2: OpenTargets disease mapping for VDR-dominant genes
# ══════════════════════════════════════════════════════════════════════════
print("\n=== PART 2: OpenTargets disease mapping ===")

vdr_dom_genes = df[df['VDR_dominant']==True]['gene'].tolist()
print(f"VDR-dominant genes: {len(vdr_dom_genes)}")

OT_URL = "https://api.platform.opentargets.org/api/v4/graphql"

def ot_query(gene_symbol):
    query = """
    query($sym: String!) {
      search(queryString: $sym, entityNames: ["target"]) {
        hits {
          id
          object {
            ... on Target {
              approvedSymbol
              associatedDiseases(page: {index: 0, size: 5}) {
                rows {
                  disease {
                    id
                    name
                    therapeuticAreas { name }
                  }
                  score
                }
              }
            }
          }
        }
      }
    }
    """
    payload = json.dumps({"query": query, "variables": {"sym": gene_symbol}}).encode()
    req = urllib.request.Request(OT_URL, data=payload,
          headers={"Content-Type": "application/json", "Accept": "application/json"})
    try:
        with urllib.request.urlopen(req, timeout=15) as r:
            return json.loads(r.read())
    except Exception as e:
        print(f"  ERROR {gene_symbol}: {e}")
        return None

disease_rows = []
print("Querying OpenTargets...")
for i, gene in enumerate(vdr_dom_genes):
    res = ot_query(gene)
    if res and res.get('data'):
        hits = res['data']['search']['hits']
        # find exact symbol match
        for hit in hits:
            obj = hit.get('object', {})
            if obj.get('approvedSymbol','').upper() == gene.upper():
                for row in obj.get('associatedDiseases', {}).get('rows', []):
                    dis = row['disease']
                    areas = [ta['name'] for ta in dis.get('therapeuticAreas', [])]
                    disease_rows.append({
                        'gene': gene,
                        'disease_id': dis['id'],
                        'disease': dis['name'],
                        'therapeutic_areas': '; '.join(areas),
                        'score': row['score'],
                        'status': df[df['gene']==gene]['status'].values[0]
                    })
                break
    if (i+1) % 10 == 0:
        print(f"  {i+1}/{len(vdr_dom_genes)} done")
    time.sleep(0.15)

dis_df = pd.DataFrame(disease_rows)
print(f"\nTotal disease associations: {len(dis_df)}")

if len(dis_df) > 0:
    dis_df.to_csv(RESULTS + 'vdr_dominant_disease_map.csv', index=False)
    print("Saved: vdr_dominant_disease_map.csv")

    # top therapeutic areas
    all_areas = []
    for ta_str in dis_df['therapeutic_areas'].dropna():
        all_areas.extend([x.strip() for x in ta_str.split(';') if x.strip()])
    area_counts = pd.Series(all_areas).value_counts()
    print("\nTop therapeutic areas:")
    print(area_counts.head(15))

# ══════════════════════════════════════════════════════════════════════════
# FIGURE: 3-panel publication-quality figure
# ══════════════════════════════════════════════════════════════════════════
print("\n=== Generating figures ===")

fig = plt.figure(figsize=(18, 14))
gs = GridSpec(2, 3, figure=fig, hspace=0.45, wspace=0.4)
ax_a = fig.add_subplot(gs[0, 0])   # 2x2 heatmap (strict)
ax_b = fig.add_subplot(gs[0, 1])   # forest plot (both tests)
ax_c = fig.add_subplot(gs[0, 2])   # scatter: VDR vs GR score, colored by status
ax_d = fig.add_subplot(gs[1, :])   # disease area bar (full width)

# ── Panel A: 2x2 contingency heatmap (strict) ─────────────────────────
mat = np.array([[n_va, n_vf],[n_ga, n_gf]], dtype=float)
mat_pct = mat / mat.sum(axis=1, keepdims=True) * 100
im = ax_a.imshow(mat_pct, cmap='Blues', vmin=0, vmax=100)
for ri in range(2):
    for ci in range(2):
        ax_a.text(ci, ri, f"{int(mat[ri,ci])}\n({mat_pct[ri,ci]:.0f}%)",
                  ha='center', va='center', fontsize=11, fontweight='bold',
                  color='white' if mat_pct[ri,ci]>55 else 'black')
ax_a.set_xticks([0,1]); ax_a.set_xticklabels(['Approved','Failed'], fontsize=10)
ax_a.set_yticks([0,1]); ax_a.set_yticklabels(['VDR-dominant','GR-dominant'], fontsize=10)
ax_a.set_title(f'A   2×2 contingency (strict)\nOR={or_strict:.1f}, p={p_strict:.3f}',
               fontsize=11, fontweight='bold', loc='left')
plt.colorbar(im, ax=ax_a, label='Row %', shrink=0.8)

# ── Panel B: Forest plot ────────────────────────────────────────────────
tests = ['Strict\n(approved vs failed)', 'Broad\n(approved vs non-approved)']
ors   = [or_strict, or_broad]
ci_los= [ci_lo, ci_lo2]
ci_his= [ci_hi, ci_hi2]
ps    = [p_strict, p_broad]
y = [1, 0]
for i,(yi,o,lo,hi,p) in enumerate(zip(y,ors,ci_los,ci_his,ps)):
    ax_b.plot([lo,hi],[yi,yi], color=COL_VDR, lw=2)
    ax_b.plot(o, yi, 'o', color=COL_VDR, ms=10, zorder=5)
    ax_b.text(hi+0.05, yi, f"OR={o:.1f}\np={p:.3f}", va='center', fontsize=9)
ax_b.axvline(1, color='gray', ls='--', lw=1)
ax_b.set_yticks(y); ax_b.set_yticklabels(tests, fontsize=10)
ax_b.set_xlabel('Odds ratio (VDR-dominant: approved vs not)', fontsize=10)
ax_b.set_title('B   Forest plot', fontsize=11, fontweight='bold', loc='left')
ax_b.set_xlim(0, max(ci_his)*1.3)

# ── Panel C: VDR vs GR scatter, colored by status ──────────────────────
status_col = {'approved':COL_APP, 'failed':COL_FAIL, 'ongoing':'#f59b00', 'other':'#aaaaaa'}
status_z   = {'approved':4, 'failed':3, 'ongoing':3, 'other':1}
for st, grp in df.groupby('status'):
    ax_c.scatter(np.log1p(grp['VDR_score']), np.log1p(grp['GR_score']),
                 c=status_col[st], alpha=0.7, s=40 if st=='other' else 70,
                 zorder=status_z[st], label=st, edgecolors='none')

# label approved VDR-dominant genes
for _, row in df[(df['status']=='approved') & (df['VDR_dominant']==True)].iterrows():
    ax_c.annotate(row['gene'],
                  (np.log1p(row['VDR_score']), np.log1p(row['GR_score'])),
                  fontsize=7, xytext=(3,3), textcoords='offset points', color=COL_APP)

diag = max(np.log1p(df['VDR_score'].max()), np.log1p(df['GR_score'].max()))
ax_c.plot([0,diag],[0,diag], 'k--', lw=0.8, alpha=0.4)
ax_c.set_xlabel('VDR score (log₁₊ₓ)', fontsize=10)
ax_c.set_ylabel('GR score (log₁₊ₓ)', fontsize=10)
ax_c.set_title('C   VDR vs GR score by drug status', fontsize=11, fontweight='bold', loc='left')
ax_c.legend(fontsize=8, markerscale=1.2)

# ── Panel D: Disease therapeutic areas ─────────────────────────────────
if len(dis_df) > 0:
    top_areas = area_counts.head(12)
    colors_d = plt.cm.tab20(np.linspace(0, 1, len(top_areas)))
    bars = ax_d.barh(range(len(top_areas)), top_areas.values[::-1],
                     color=colors_d[::-1], edgecolor='white', height=0.7)
    ax_d.set_yticks(range(len(top_areas)))
    ax_d.set_yticklabels(top_areas.index[::-1], fontsize=10)
    ax_d.set_xlabel('Number of VDR-dominant gene associations', fontsize=11)
    ax_d.set_title('D   Therapeutic areas enriched in VDR-dominant genes (OpenTargets)',
                   fontsize=11, fontweight='bold', loc='left')
    for bar, val in zip(bars, top_areas.values[::-1]):
        ax_d.text(val+0.3, bar.get_y()+bar.get_height()/2,
                  str(val), va='center', fontsize=9)
    ax_d.spines['top'].set_visible(False)
    ax_d.spines['right'].set_visible(False)
else:
    ax_d.text(0.5, 0.5, 'OpenTargets query failed\n(no network or API error)',
              ha='center', va='center', transform=ax_d.transAxes, fontsize=12, color='gray')
    ax_d.set_title('D   Therapeutic areas (OpenTargets)', fontsize=11, fontweight='bold', loc='left')

for ax in [ax_a,ax_b,ax_c]:
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

fig.suptitle('VDR-dominant genes are enriched for approved drugs\nand span major chronic inflammatory disease areas',
             fontsize=13, fontweight='bold', y=0.98)

out = RESULTS + 'fig_vdr_approval_disease'
fig.savefig(out + '.pdf', dpi=300, bbox_inches='tight', facecolor='white')
fig.savefig(out + '.png', dpi=200, bbox_inches='tight', facecolor='white')
print(f"\nSaved: {out}.pdf / .png")

# ── Summary stats to console ────────────────────────────────────────────
print("\n=== SUMMARY ===")
print(f"VDR-dominant genes: {df['VDR_dominant'].sum()} / {len(df)}")
print(f"Approved VDR-dominant: {n_va} / {df['VDR_dominant'].sum()} ({n_va/df['VDR_dominant'].sum()*100:.1f}%)")
print(f"Approved GR-dominant:  {n_ga} / {(~df['VDR_dominant']).sum()} ({n_ga/(~df['VDR_dominant']).sum()*100:.1f}%)")
print(f"\nStrict Fisher's: OR={or_strict:.2f} (95%CI {ci_lo:.2f}–{ci_hi:.2f}), p={p_strict:.4f}")
print(f"Broad  Fisher's: OR={or_broad:.2f} (95%CI {ci_lo2:.2f}–{ci_hi2:.2f}), p={p_broad:.4f}")
if len(dis_df) > 0:
    print(f"\nTop disease area: {area_counts.index[0]} (n={area_counts.iloc[0]})")
    print(f"Genes with disease data: {dis_df['gene'].nunique()} / {len(vdr_dom_genes)}")

print("\nDONE")
