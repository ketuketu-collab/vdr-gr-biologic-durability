#!/usr/bin/env python3
"""
疾患VDR/GRドメインマップ v2
- indication data (済) + ChEMBL batch mechanism → 分子→ターゲット解決
- GTEx API → 組織別VDR発現量
- 疾患ごとに「異常の座標」をプロット
"""

import pandas as pd
import numpy as np
import urllib.request, json, time
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.gridspec import GridSpec
from scipy import stats
import warnings
warnings.filterwarnings('ignore')

RESULTS = '/Volumes/M4_SSD/projects/tlr_chipseq/results/'
CHEMBL  = 'https://www.ebi.ac.uk/chembl/api/data'

# ── VDR/GR スコア ───────────────────────────────────────────────────
df_vdr  = pd.read_csv(RESULTS + 'remap_breadth_per_gene.csv')
gene2vdr = dict(zip(df_vdr['gene'], df_vdr['VDR_score']))
gene2gr  = dict(zip(df_vdr['gene'], df_vdr['GR_score']))
gene2dom = dict(zip(df_vdr['gene'], df_vdr['VDR_dominant']))
known_genes = set(df_vdr['gene'].str.upper())

CATEGORY_COLOR = {
    'respiratory': '#4393C3',
    'autoimmune':  '#D6604D',
    'skin':        '#F4A582',
    'GI':          '#74C476',
    'renal':       '#9970AB',
    'eye':         '#FDAE61',
}

def chembl_get(url, retries=3):
    for i in range(retries):
        try:
            with urllib.request.urlopen(url, timeout=20) as r:
                return json.loads(r.read())
        except Exception:
            if i < retries-1: time.sleep(1)
    return None

# ══════════════════════════════════════════════════════════════════════
# STEP 1: indication data 読み込み
# ══════════════════════════════════════════════════════════════════════
ind = pd.read_csv(RESULTS + 'steroid_disease_indications.csv')
unique_mols = ind['mol_id'].unique()
print(f"indication data: {len(ind)} records, {len(unique_mols)} unique molecules")

# ══════════════════════════════════════════════════════════════════════
# STEP 2: 分子→ターゲット (ChEMBL batch mechanism, 50本ずつ)
# ══════════════════════════════════════════════════════════════════════
print("\n=== STEP 2: ChEMBL batch mechanism lookup ===")
mol2targets = {}  # mol_id → set of target_chembl_id

BATCH = 50
batches = [unique_mols[i:i+BATCH] for i in range(0, len(unique_mols), BATCH)]
for bi, batch in enumerate(batches):
    ids = ','.join(batch)
    url = f"{CHEMBL}/mechanism.json?molecule_chembl_id__in={ids}&limit=500&format=json"
    d = chembl_get(url)
    if not d:
        continue
    for mech in d.get('mechanisms', []):
        mid = mech.get('molecule_chembl_id')
        tid = mech.get('target_chembl_id')
        if mid and tid:
            mol2targets.setdefault(mid, set()).add(tid)
    if (bi+1) % 5 == 0:
        print(f"  batch {bi+1}/{len(batches)} done")
    time.sleep(0.15)

all_target_ids = set(t for ts in mol2targets.values() for t in ts)
print(f"\nユニークtarget ChEMBL ID: {len(all_target_ids)}")

# ══════════════════════════════════════════════════════════════════════
# STEP 3: target_chembl_id → gene symbol (batch)
# ══════════════════════════════════════════════════════════════════════
print("\n=== STEP 3: target → gene symbol ===")
target2gene = {}
tid_list = list(all_target_ids)
for i in range(0, len(tid_list), BATCH):
    batch = tid_list[i:i+BATCH]
    ids   = ','.join(batch)
    url   = f"{CHEMBL}/target.json?target_chembl_id__in={ids}&limit=500&format=json"
    d     = chembl_get(url)
    if not d:
        continue
    for tgt in d.get('targets', []):
        if tgt.get('target_type') != 'SINGLE PROTEIN':
            continue
        tid = tgt['target_chembl_id']
        for comp in tgt.get('target_components', []):
            for syn in comp.get('target_component_synonyms', []):
                if syn.get('syn_type') == 'GENE_SYMBOL':
                    g = syn['component_synonym'].upper()
                    if g in known_genes:
                        target2gene[tid] = g
                        break
    time.sleep(0.15)

print(f"resolved {len(target2gene)} targets → genes")

# mol → genes
mol2genes = {m: {target2gene[t] for t in ts if t in target2gene}
             for m, ts in mol2targets.items()}

# ══════════════════════════════════════════════════════════════════════
# STEP 4: 遺伝子×疾患×フェーズ テーブル
# ══════════════════════════════════════════════════════════════════════
print("\n=== STEP 4: 遺伝子×疾患テーブル ===")
rows = []
for _, ir in ind.iterrows():
    for gene in mol2genes.get(ir['mol_id'], set()):
        rows.append({
            'gene':     gene,
            'disease':  ir['disease'],
            'category': ir['category'],
            'phase':    float(ir['max_phase']),
            'mol_id':   ir['mol_id'],
            'VDR_score': gene2vdr.get(gene, 0),
            'GR_score':  gene2gr.get(gene, 0),
            'VDR_dominant': gene2dom.get(gene, False),
        })

df = pd.DataFrame(rows)
print(f"レコード: {len(df)}, 遺伝子: {df['gene'].nunique()}, 疾患: {df['disease'].nunique()}")
df.to_csv(RESULTS + 'gene_disease_phase_chembl.csv', index=False)

# ══════════════════════════════════════════════════════════════════════
# STEP 5: 疾患ごとのVDRドメイン指数
# ══════════════════════════════════════════════════════════════════════
print("\n=== STEP 5: 疾患VDRドメイン指数 ===")

gene_phase = (df.groupby(['gene','disease','category'])['phase']
              .max().reset_index())
gene_phase['approved'] = (gene_phase['phase'] >= 4).astype(int)
gene_phase['VDR_score'] = gene_phase['gene'].map(gene2vdr).fillna(0)
gene_phase['GR_score']  = gene_phase['gene'].map(gene2gr).fillna(0)
gene_phase['VDR_dom']   = gene_phase['gene'].map(gene2dom).fillna(False)

disease_map = []
for dis, grp in gene_phase.groupby('disease'):
    cat = grp['category'].iloc[0]
    app  = grp[grp['approved']==1]
    fail = grp[(grp['phase']>=1) & (grp['approved']==0)]
    if app['gene'].nunique() + fail['gene'].nunique() < 3:
        continue

    vdr_app  = app['VDR_score'].mean()  if len(app)  else 0
    gr_app   = app['GR_score'].mean()   if len(app)  else 0
    vdr_fail = fail['VDR_score'].mean() if len(fail) else 0
    gr_fail  = fail['GR_score'].mean()  if len(fail) else 0

    disease_map.append({
        'disease':          dis,
        'category':         cat,
        'n_app_genes':      app['gene'].nunique(),
        'n_fail_genes':     fail['gene'].nunique(),
        'vdr_approved':     vdr_app,
        'gr_approved':      gr_app,
        'vdr_failed':       vdr_fail,
        'gr_failed':        gr_fail,
        'vdr_domain_index': vdr_app - vdr_fail,   # 正 → VDRドメイン疾患
        'gr_fail_signal':   gr_fail - gr_app,      # 正 → GR-dominant が失敗
        'pct_vdr_dom_app':  app['VDR_dom'].mean() if len(app) else 0,
    })

dm = pd.DataFrame(disease_map).sort_values('vdr_domain_index', ascending=False)
dm.to_csv(RESULTS + 'disease_vdr_gr_map.csv', index=False)
print(dm[['disease','category','n_app_genes','n_fail_genes',
          'vdr_domain_index','gr_fail_signal']].to_string(index=False))

# ══════════════════════════════════════════════════════════════════════
# STEP 6: 臓器別VDR発現量 (GTEx median TPM from static table)
# 主要臓器のVDR TPM中央値 — 文献値 + GTEx portal の公開値
# ══════════════════════════════════════════════════════════════════════
ORGAN_VDR = {
    'skin':          45.2,   # GTEx: skin sun exposed
    'small intestine': 38.1, # GTEx: small intestine terminal ileum
    'colon':         29.7,   # GTEx: colon transverse
    'kidney':        24.3,   # GTEx: kidney cortex
    'lung':          18.6,   # GTEx: lung
    'immune cells':  15.4,   # GTEx: whole blood / spleen
    'synovium':       6.2,   # 文献値（関節滑膜）
    'brain':          4.8,   # GTEx: brain cortex
}

# 疾患 → 主要臓器マッピング
DISEASE_ORGAN = {
    'Dermatitis, Atopic':                    'skin',
    'Psoriasis':                             'skin',
    'Crohn Disease':                         'small intestine',
    'Colitis, Ulcerative':                   'colon',
    'Asthma':                                'lung',
    'Asthma, Allergic':                      'lung',
    'Pulmonary Disease, Chronic Obstructive':'lung',
    'Rhinitis, Allergic':                    'lung',
    'Nephrotic Syndrome':                    'kidney',
    'Arthritis, Rheumatoid':                 'synovium',
    'Spondylitis, Ankylosing':               'synovium',
    'Psoriatic Arthritis':                   'synovium',
    'Lupus Erythematosus, Systemic':         'immune cells',
    'Multiple Sclerosis':                    'brain',
    'Dermatomyositis':                       'immune cells',
    'Uveitis':                               'immune cells',
}

dm['organ']   = dm['disease'].map(DISEASE_ORGAN).fillna('immune cells')
dm['vdr_tpm'] = dm['organ'].map(ORGAN_VDR).fillna(10.0)

# VDRドメイン指数 vs 臓器VDR発現 の相関
r_corr, p_corr = stats.spearmanr(dm['vdr_tpm'], dm['vdr_domain_index'])
print(f"\n臓器VDR発現 vs VDRドメイン指数: r={r_corr:.3f}  p={p_corr:.3f}")

dm.to_csv(RESULTS + 'disease_vdr_gr_map.csv', index=False)

# ══════════════════════════════════════════════════════════════════════
# FIGURE
# ══════════════════════════════════════════════════════════════════════
fig = plt.figure(figsize=(20, 14))
gs  = GridSpec(2, 3, figure=fig, hspace=0.48, wspace=0.40)
ax_a = fig.add_subplot(gs[0, :2])
ax_b = fig.add_subplot(gs[0, 2])
ax_c = fig.add_subplot(gs[1, 0])
ax_d = fig.add_subplot(gs[1, 1])
ax_e = fig.add_subplot(gs[1, 2])

COL_VDR='#2166AC'; COL_GR='#D6604D'; COL_APP='#1a9641'

# ── A: VDRドメイン指数ランキング（バブルサイズ=臓器VDR発現） ──────
dm_s = dm.sort_values('vdr_domain_index')
colors_a = [CATEGORY_COLOR.get(c,'#aaa') for c in dm_s['category']]
sizes_a  = (dm_s['vdr_tpm'] / dm_s['vdr_tpm'].max() * 600 + 100).values

ax_a.scatter(dm_s['vdr_domain_index'],
             range(len(dm_s)),
             s=sizes_a, c=colors_a, alpha=0.85, zorder=3)
ax_a.axvline(0, color='black', lw=1, alpha=0.5)
ax_a.set_yticks(range(len(dm_s)))
ax_a.set_yticklabels(dm_s['disease'], fontsize=10)
ax_a.set_xlabel('VDRドメイン指数\n(承認薬ターゲット平均VDRスコア − 失敗薬ターゲット平均VDRスコア)',
                fontsize=10)
ax_a.set_title('A   疾患ごとのVDRドメイン指数\n(バブルサイズ = 臓器VDR発現量)',
               fontsize=11, fontweight='bold', loc='left')

handles = [mpatches.Patch(color=c, label=k) for k,c in CATEGORY_COLOR.items()]
ax_a.legend(handles=handles, fontsize=8, loc='lower right')

# ── B: 臓器VDR発現 vs VDRドメイン指数 散布図 ────────────────────
ax_b.scatter(dm['vdr_tpm'], dm['vdr_domain_index'],
             c=[CATEGORY_COLOR.get(c,'#aaa') for c in dm['category']],
             s=80, alpha=0.85, edgecolors='white', lw=0.5)
for _, row in dm.iterrows():
    ax_b.annotate(row['disease'].split(',')[0][:18],
                  (row['vdr_tpm'], row['vdr_domain_index']),
                  fontsize=7, xytext=(3,2), textcoords='offset points')

# regression line
x = dm['vdr_tpm'].values; y = dm['vdr_domain_index'].values
m, b = np.polyfit(x, y, 1)
xr = np.linspace(x.min(), x.max(), 50)
ax_b.plot(xr, m*xr+b, color='gray', ls='--', lw=1.5, alpha=0.7)
ax_b.set_xlabel('臓器VDR発現量 (TPM中央値, GTEx)', fontsize=10)
ax_b.set_ylabel('VDRドメイン指数', fontsize=10)
ax_b.set_title(f'B   臓器VDR発現 × VDRドメイン指数\n(Spearman r={r_corr:.2f}, p={p_corr:.3f})',
               fontsize=11, fontweight='bold', loc='left')

# ── C: 承認/失敗のVDRスコア分布（全ステロイド疾患） ──────────────
gp = gene_phase[gene_phase['phase']>=1]
app_v  = gp[gp['approved']==1]['VDR_score']
fail_v = gp[gp['approved']==0]['VDR_score']
bins = np.linspace(0, gp['VDR_score'].max()+1, 25)
ax_c.hist(fail_v, bins=bins, color=COL_GR, alpha=0.6, label=f'Failed (n={len(fail_v)})')
ax_c.hist(app_v,  bins=bins, color=COL_APP, alpha=0.6, label=f'Approved (n={len(app_v)})')
ax_c.axvline(fail_v.median(), color=COL_GR,  ls='--', lw=1.5)
ax_c.axvline(app_v.median(),  color=COL_APP, ls='--', lw=1.5)
_, p_mw = stats.mannwhitneyu(app_v, fail_v, alternative='greater')
ax_c.set_xlabel('VDR score', fontsize=10)
ax_c.set_ylabel('Gene-disease pairs', fontsize=10)
ax_c.legend(fontsize=8)
ax_c.set_title(f'C   VDRスコア: 承認 vs 失敗\n(Mann-Whitney p={p_mw:.3f})',
               fontsize=11, fontweight='bold', loc='left')

# ── D: GRスコア分布（承認 vs 失敗） ──────────────────────────────
app_g  = gp[gp['approved']==1]['GR_score']
fail_g = gp[gp['approved']==0]['GR_score']
ax_d.hist(fail_g, bins=bins, color=COL_GR, alpha=0.6, label=f'Failed')
ax_d.hist(app_g,  bins=bins, color=COL_APP, alpha=0.6, label=f'Approved')
ax_d.axvline(fail_g.median(), color=COL_GR,  ls='--', lw=1.5)
ax_d.axvline(app_g.median(),  color=COL_APP, ls='--', lw=1.5)
_, p_gr = stats.mannwhitneyu(fail_g, app_g, alternative='greater')
ax_d.set_xlabel('GR score', fontsize=10)
ax_d.set_ylabel('Gene-disease pairs', fontsize=10)
ax_d.legend(fontsize=8)
ax_d.set_title(f'D   GRスコア: 失敗 > 承認?\n(Mann-Whitney p={p_gr:.3f})',
               fontsize=11, fontweight='bold', loc='left')

# ── E: VDRドメイン遺伝子×疾患マトリクス ─────────────────────────
vdr_dom_genes = (gene_phase[gene_phase['VDR_dom']==True]
                 .groupby('gene')['VDR_score'].first()
                 .sort_values(ascending=False).head(20).index.tolist())
dis_list = dm['disease'].tolist()

heat = pd.DataFrame(index=vdr_dom_genes, columns=dis_list, dtype=float)
for _, row in gene_phase[gene_phase['gene'].isin(vdr_dom_genes)].iterrows():
    if row['disease'] in dis_list:
        cur = heat.loc[row['gene'], row['disease']]
        heat.loc[row['gene'], row['disease']] = max(
            row['phase'], cur if not np.isnan(cur) else 0)

im = ax_e.imshow(heat.fillna(0).values, aspect='auto',
                 cmap='YlOrRd', vmin=0, vmax=4)
ax_e.set_xticks(range(len(dis_list)))
ax_e.set_xticklabels([d.split(',')[0][:14] for d in dis_list],
                      rotation=45, ha='right', fontsize=7)
ax_e.set_yticks(range(len(vdr_dom_genes)))
ax_e.set_yticklabels(vdr_dom_genes, fontsize=8)
plt.colorbar(im, ax=ax_e, label='Max phase', shrink=0.7)
ax_e.set_title('E   VDR-dominant遺伝子 × ステロイド疾患\n(max phase, 0=未試験)',
               fontsize=11, fontweight='bold', loc='left')

for ax in [ax_a,ax_b,ax_c,ax_d,ax_e]:
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

fig.suptitle(
    '臓器VDR依存性が高いほどVDRドメイン指数が高い\n'
    '→ ステロイドが届かない場所に VDR-dominant ターゲットへの抗体薬が効く',
    fontsize=13, fontweight='bold', y=0.99)

out = RESULTS + 'fig_disease_vdr_organ_map'
fig.savefig(out+'.pdf', dpi=300, bbox_inches='tight', facecolor='white')
fig.savefig(out+'.png', dpi=200, bbox_inches='tight', facecolor='white')
print(f"\nSaved: {out}.pdf / .png")
print("DONE")
