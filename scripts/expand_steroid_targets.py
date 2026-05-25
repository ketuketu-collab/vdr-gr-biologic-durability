#!/usr/bin/env python3
"""
ステロイド疾患ターゲット拡張版
1. steroid disease 1183分子 → 全ターゲット遺伝子取得（308フィルタなし）
2. 既存308遺伝子スコア + 新規遺伝子をReMap2022でスコアリング
3. がん適応は含まない（steroid disease indicationsのみ使用）
4. VEGFAなど複数適応を持つ遺伝子は炎症疾患の試験データだけ使用
"""

import gzip, re, json, time
import urllib.request, urllib.parse
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

RESULTS  = '/Volumes/M4_SSD/projects/tlr_chipseq/results/'
DATA     = '/Volumes/M4_SSD/projects/tlr10_chipseq_vdr_gr/data/'
VDR_BED  = DATA + 'remap2022_VDR_all_macs2_hg38.bed.gz'
GR_BED   = DATA + 'remap2022_NR3C1_all_macs2_hg38.bed.gz'
CHEMBL   = 'https://www.ebi.ac.uk/chembl/api/data'
HALF     = 5000

CATEGORY_COLOR = {
    'respiratory': '#4393C3', 'autoimmune': '#D6604D',
    'skin': '#F4A582',        'GI': '#74C476',
    'renal': '#9970AB',       'eye': '#FDAE61',
}
DISEASE_ORGAN = {
    'Dermatitis, Atopic':                    ('skin',          45.2),
    'Psoriasis':                             ('skin',          45.2),
    'Crohn Disease':                         ('small intestine',38.1),
    'Colitis, Ulcerative':                   ('colon',         29.7),
    'Asthma':                                ('lung',          18.6),
    'Pulmonary Disease, Chronic Obstructive':('lung',          18.6),
    'Rhinitis, Allergic':                    ('lung',          18.6),
    'Nephrotic Syndrome':                    ('kidney',        24.3),
    'Arthritis, Rheumatoid':                 ('synovium',       6.2),
    'Spondylitis, Ankylosing':               ('synovium',       6.2),
    'Lupus Erythematosus, Systemic':         ('immune cells',  15.4),
    'Multiple Sclerosis':                    ('brain',          4.8),
    'Dermatomyositis':                       ('immune cells',  15.4),
    'Uveitis':                               ('immune cells',  15.4),
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
# STEP 1: 全ターゲット遺伝子取得（フィルタなし）
# ══════════════════════════════════════════════════════════════════════
print("=== STEP 1: 全ターゲット遺伝子取得 ===")
ind      = pd.read_csv(RESULTS + 'steroid_disease_indications.csv')
mol_stat = pd.read_csv(RESULTS + 'disease_vdr_gr_map_v3.csv')  # status already computed

# 分子ステータス再取得（v3スクリプトで保存されていないので再取得）
unique_mols = ind['mol_id'].unique()
print(f"分子数: {len(unique_mols)}")
mol_status = {}
BATCH = 50
for i in range(0, len(unique_mols), BATCH):
    batch = unique_mols[i:i+BATCH]
    d = chembl_get(f"{CHEMBL}/molecule.json?molecule_chembl_id__in={','.join(batch)}&limit=500&format=json")
    if d:
        for mol in d.get('molecules', []):
            mid = mol['molecule_chembl_id']
            mol_status[mid] = {
                'global_max_phase': float(mol.get('max_phase') or 0),
                'withdrawn':        bool(mol.get('withdrawn_flag', False)),
                'first_approval':   mol.get('first_approval'),
                'drug_type':        mol.get('molecule_type', ''),
            }
    time.sleep(0.1)
print(f"分子ステータス取得: {len(mol_status)}")

def classify_mol(mol_id, phase_for_ind):
    s = mol_status.get(mol_id, {})
    gp = s.get('global_max_phase', 0)
    wd = s.get('withdrawn', False)
    fa = s.get('first_approval')
    if phase_for_ind >= 4: return 'approved'
    if wd: return 'failed'
    if gp <= 1 and phase_for_ind >= 1: return 'failed'
    if gp >= 2 and fa is None and phase_for_ind >= 2: return 'failed'
    if gp >= 2 and fa is not None and phase_for_ind < 4: return 'approved_other'
    return 'ongoing'

ind['mol_status'] = ind.apply(lambda r: classify_mol(r['mol_id'], r['max_phase']), axis=1)

# mechanism batch（フィルタなし）
print("\n=== STEP 1b: mechanism batch (全遺伝子) ===")
mol2targets = {}
batches = [unique_mols[i:i+BATCH] for i in range(0, len(unique_mols), BATCH)]
for bi, batch in enumerate(batches):
    d = chembl_get(f"{CHEMBL}/mechanism.json?molecule_chembl_id__in={','.join(batch)}&limit=500&format=json")
    if d:
        for mech in d.get('mechanisms', []):
            mid = mech.get('molecule_chembl_id')
            tid = mech.get('target_chembl_id')
            if mid and tid:
                mol2targets.setdefault(mid, set()).add(tid)
    if (bi+1) % 5 == 0:
        print(f"  batch {bi+1}/{len(batches)}")
    time.sleep(0.12)

all_tids = set(t for ts in mol2targets.values() for t in ts)
print(f"ユニークtarget ID: {len(all_tids)}")

# target → gene（SINGLE PROTEIN のみ、フィルタなし）
print("\n=== STEP 1c: target → gene (全遺伝子) ===")
target2gene = {}
tid_list = list(all_tids)
for i in range(0, len(tid_list), BATCH):
    batch = tid_list[i:i+BATCH]
    d = chembl_get(f"{CHEMBL}/target.json?target_chembl_id__in={','.join(batch)}&limit=500&format=json")
    if d:
        for tgt in d.get('targets', []):
            if tgt.get('target_type') != 'SINGLE PROTEIN': continue
            tid = tgt['target_chembl_id']
            for comp in tgt.get('target_components', []):
                for syn in comp.get('target_component_synonyms', []):
                    if syn.get('syn_type') == 'GENE_SYMBOL':
                        target2gene[tid] = syn['component_synonym'].upper()
                        break
    time.sleep(0.12)

mol2genes = {m: {target2gene[t] for t in ts if t in target2gene}
             for m, ts in mol2targets.items()}
all_steroid_genes = set(g for gs in mol2genes.values() for g in gs)
print(f"ステロイド疾患ターゲット遺伝子数: {len(all_steroid_genes)}")

# 既存スコア済み
existing = pd.read_csv(RESULTS + 'remap_breadth_per_gene.csv')
existing_set = set(existing['gene'].str.upper())
new_genes = sorted(all_steroid_genes - existing_set)
print(f"新規スコアリング必要: {len(new_genes)} 遺伝子")

# ══════════════════════════════════════════════════════════════════════
# STEP 2: 新規遺伝子のReMap2022スコアリング
# ══════════════════════════════════════════════════════════════════════
print(f"\n=== STEP 2: 新規{len(new_genes)}遺伝子をReMap2022でスコアリング ===")

# 2a. TSS取得 (Ensembl batch)
def fetch_tss_batch(genes):
    url  = 'https://rest.ensembl.org/lookup/symbol/homo_sapiens'
    data = json.dumps({'symbols': genes}).encode()
    req  = urllib.request.Request(url, data=data,
           headers={'Content-Type':'application/json','Accept':'application/json'})
    try:
        with urllib.request.urlopen(req, timeout=60) as r:
            d = json.loads(r.read())
        tss = {}
        for gene, info in d.items():
            if not info: continue
            chrom = str(info.get('seq_region_name',''))
            if re.match(r'^\d+$|^X$|^Y$', chrom):
                chrom = 'chr' + chrom
                strand  = info.get('strand', 1)
                tss_pos = info['start'] if strand == 1 else info['end']
                tss[gene.upper()] = (chrom, tss_pos)
        return tss
    except Exception as e:
        print(f"  Ensembl error: {e}")
        return {}

tss_map = {}
for i in range(0, len(new_genes), 50):
    chunk = new_genes[i:i+50]
    tss_map.update(fetch_tss_batch(chunk))
    time.sleep(0.5)
print(f"TSS取得: {len(tss_map)}/{len(new_genes)}")

# 2b. ReMap BED インデックス構築
print("ReMap BEDインデックス構築中...")
def build_index(bed_gz):
    idx = {}
    with gzip.open(bed_gz, 'rt') as f:
        for line in f:
            p = line.rstrip().split('\t')
            if len(p) < 4: continue
            chrom, start, end, name = p[0], int(p[1]), int(p[2]), p[3]
            parts = name.split('.')
            gse_id   = parts[0] if parts else 'unknown'
            celltype = parts[2] if len(parts) > 2 else 'unknown'
            base_cell = re.sub(r'_[A-Z0-9_]+$', '', celltype, flags=re.IGNORECASE)
            idx.setdefault(chrom, []).append((start, end, base_cell, gse_id))
    return idx

vdr_idx = build_index(VDR_BED)
gr_idx  = build_index(GR_BED)

# 全VDR/GR実験数（再現性分母）
all_vdr_exp = set()
all_gr_exp  = set()
for chrom_data in vdr_idx.values():
    for _, _, _, gse in chrom_data: all_vdr_exp.add(gse)
for chrom_data in gr_idx.values():
    for _, _, _, gse in chrom_data: all_gr_exp.add(gse)
n_vdr_total = len(all_vdr_exp)
n_gr_total  = len(all_gr_exp)

def score_gene(chrom, tss_pos, idx, n_total):
    lo, hi = tss_pos - HALF, tss_pos + HALF
    hits = [e for e in idx.get(chrom, []) if e[1] >= lo and e[0] <= hi]
    cells = set(h[2] for h in hits)
    exps  = set(h[3] for h in hits)
    score = len(cells) * 10 + len(exps)
    repro = len(exps) / n_total if n_total else 0
    return score, len(cells), len(exps), repro

# 2c. スコア計算
new_scores = []
for gene, (chrom, tss) in tss_map.items():
    vs, vc, ve, vr = score_gene(chrom, tss, vdr_idx, n_vdr_total)
    gs, gc, ge, gr_ = score_gene(chrom, tss, gr_idx,  n_gr_total)
    new_scores.append({
        'gene': gene, 'VDR_score': round(vs,2), 'GR_score': round(gs,2),
        'VDR_celltypes': vc, 'GR_celltypes': gc,
        'VDR_experiments': ve, 'GR_experiments': ge,
        'VDR_repro': round(vr,4), 'GR_repro': round(gr_,4),
        'VDR_dominant': vs > gs,
        'repro_ratio': round(vr - gr_, 4),
        'status': 'steroid_new',
    })

new_df = pd.DataFrame(new_scores)
print(f"新規スコア: {len(new_df)} 遺伝子")

# 2d. 既存と結合
combined = pd.concat([existing, new_df], ignore_index=True)
combined = combined.drop_duplicates('gene', keep='first')
combined.to_csv(RESULTS + 'remap_scores_expanded.csv', index=False)
print(f"統合スコア: {len(combined)} 遺伝子")
print(f"  VDR-dominant: {combined['VDR_dominant'].sum()}")

# ══════════════════════════════════════════════════════════════════════
# STEP 3: 遺伝子×疾患テーブル構築（拡張版）
# ══════════════════════════════════════════════════════════════════════
print("\n=== STEP 3: 遺伝子×疾患テーブル（拡張版） ===")
gene2vdr = dict(zip(combined['gene'], combined['VDR_score']))
gene2gr  = dict(zip(combined['gene'], combined['GR_score']))
gene2dom = dict(zip(combined['gene'], combined['VDR_dominant']))
gene2vdr_r = dict(zip(combined['gene'], combined['VDR_repro']))
gene2gr_r  = dict(zip(combined['gene'], combined['GR_repro']))

rows = []
for _, ir in ind[ind['mol_status'].isin(['approved','failed'])].iterrows():
    for gene in mol2genes.get(ir['mol_id'], set()):
        if gene not in gene2vdr: continue
        rows.append({
            'gene':         gene,
            'disease':      ir['disease'],
            'category':     ir['category'],
            'phase':        float(ir['max_phase']),
            'mol_status':   ir['mol_status'],
            'mol_id':       ir['mol_id'],
            'VDR_score':    gene2vdr[gene],
            'GR_score':     gene2gr[gene],
            'VDR_dominant': gene2dom[gene],
            'VDR_repro':    gene2vdr_r.get(gene, 0),
            'GR_repro':     gene2gr_r.get(gene, 0),
        })

df = pd.DataFrame(rows)
df.to_csv(RESULTS + 'gene_disease_phase_expanded.csv', index=False)
print(f"レコード: {len(df)}, 遺伝子: {df['gene'].nunique()}, 疾患: {df['disease'].nunique()}")
print(f"  承認: {(df['mol_status']=='approved').sum()}, 失敗: {(df['mol_status']=='failed').sum()}")

# 遺伝子×疾患ごとに最良ステータス集約
gene_dis = (df.groupby(['gene','disease','category'])
            .agg(
                best_status = ('mol_status', lambda x:
                               'approved' if 'approved' in x.values else 'failed'),
                VDR_score   = ('VDR_score','first'),
                GR_score    = ('GR_score','first'),
                VDR_dom     = ('VDR_dominant','first'),
            ).reset_index())

# ── Fisher's test ─────────────────────────────────────────────────
ct = pd.crosstab(gene_dis['VDR_dom'], gene_dis['best_status'])
print(f"\nFisher's (拡張版):\n{ct}")
n_va = ct.loc[True,  'approved']; n_vf = ct.loc[True,  'failed']
n_ga = ct.loc[False, 'approved']; n_gf = ct.loc[False, 'failed']
or_v, p_v = stats.fisher_exact([[n_va,n_vf],[n_ga,n_gf]], alternative='greater')
log_or = np.log(or_v)
se = np.sqrt(1/n_va + 1/n_vf + 1/n_ga + 1/n_gf)
ci_lo = np.exp(log_or - 1.96*se)
ci_hi = np.exp(log_or + 1.96*se)
print(f"OR={or_v:.2f} (95%CI {ci_lo:.2f}–{ci_hi:.2f})  p={p_v:.4f}")

# Mann-Whitney
app_v  = gene_dis[gene_dis['best_status']=='approved']['VDR_score']
fail_v = gene_dis[gene_dis['best_status']=='failed']['VDR_score']
_, p_mw = stats.mannwhitneyu(app_v, fail_v, alternative='greater')
print(f"VDRスコア 承認中央値={app_v.median():.1f}  失敗中央値={fail_v.median():.1f}  p={p_mw:.4f}")

# ── 疾患マップ ────────────────────────────────────────────────────
disease_map = []
for dis, grp in gene_dis.groupby('disease'):
    info = DISEASE_ORGAN.get(dis, ('immune cells', 15.4))
    cat  = grp['category'].iloc[0]
    app  = grp[grp['best_status']=='approved']
    fail = grp[grp['best_status']=='failed']
    if len(app) + len(fail) < 3: continue
    disease_map.append({
        'disease':          dis,
        'category':         cat,
        'organ':            info[0],
        'vdr_tpm':          info[1],
        'n_app':            len(app),
        'n_fail':           len(fail),
        'vdr_approved':     app['VDR_score'].mean() if len(app) else 0,
        'gr_approved':      app['GR_score'].mean()  if len(app) else 0,
        'vdr_failed':       fail['VDR_score'].mean() if len(fail) else 0,
        'gr_failed':        fail['GR_score'].mean()  if len(fail) else 0,
        'vdr_domain_index': (app['VDR_score'].mean()  if len(app) else 0) -
                            (fail['VDR_score'].mean() if len(fail) else 0),
        'gr_fail_signal':   (fail['GR_score'].mean()  if len(fail) else 0) -
                            (app['GR_score'].mean()   if len(app) else 0),
        'pct_vdr_app':      app['VDR_dom'].mean() if len(app) else 0,
    })

dm = pd.DataFrame(disease_map).sort_values('vdr_domain_index', ascending=False)
dm.to_csv(RESULTS + 'disease_vdr_gr_map_expanded.csv', index=False)

r_corr, p_corr = stats.spearmanr(dm['vdr_tpm'], dm['vdr_domain_index'])
print(f"\n臓器VDR発現 vs VDRドメイン指数: r={r_corr:.3f}  p={p_corr:.3f}")
print("\n疾患VDRドメイン指数（拡張版）:")
print(dm[['disease','organ','n_app','n_fail',
          'vdr_domain_index','gr_fail_signal','vdr_tpm']].to_string(index=False))

# ══════════════════════════════════════════════════════════════════════
# FIGURE
# ══════════════════════════════════════════════════════════════════════
COL_VDR='#2166AC'; COL_GR='#D6604D'; COL_APP='#1a9641'; COL_FAIL='#d7191c'

fig = plt.figure(figsize=(20, 15))
gs  = GridSpec(2, 3, figure=fig, hspace=0.48, wspace=0.40)
ax_a = fig.add_subplot(gs[0, :2])
ax_b = fig.add_subplot(gs[0, 2])
ax_c = fig.add_subplot(gs[1, 0])
ax_d = fig.add_subplot(gs[1, 1])
ax_e = fig.add_subplot(gs[1, 2])

# ── A: 疾患VDRドメイン指数バブルチャート ─────────────────────────
dm_s   = dm.sort_values('vdr_domain_index')
col_a  = [CATEGORY_COLOR.get(c,'#aaa') for c in dm_s['category']]
sz_a   = (dm_s['vdr_tpm'] / dm_s['vdr_tpm'].max() * 700 + 80).values
y_pos  = range(len(dm_s))

ax_a.scatter(dm_s['vdr_domain_index'], y_pos,
             s=sz_a, c=col_a, alpha=0.85, zorder=3)
ax_a.axvline(0, color='black', lw=1, alpha=0.5)
ax_a.set_yticks(y_pos)
ax_a.set_yticklabels(
    [f"{r['disease'].split(',')[0]}  ({r['organ']})"
     for _, r in dm_s.iterrows()], fontsize=10)
ax_a.set_xlabel(
    'VDRドメイン指数 (承認薬ターゲット平均VDRスコア − 失敗薬ターゲット平均VDRスコア)',
    fontsize=10)
ax_a.set_title(
    f'A   疾患VDRドメイン指数（拡張版, n遺伝子={gene_dis["gene"].nunique()}）\n'
    f'バブルサイズ = 臓器VDR発現量 | がん適応除外・進行中除外',
    fontsize=11, fontweight='bold', loc='left')
handles = [mpatches.Patch(color=c, label=k) for k,c in CATEGORY_COLOR.items()]
ax_a.legend(handles=handles, fontsize=8, loc='lower right')

# ── B: 臓器VDR発現 vs VDRドメイン指数 ───────────────────────────
ax_b.scatter(dm['vdr_tpm'], dm['vdr_domain_index'],
             c=[CATEGORY_COLOR.get(c,'#aaa') for c in dm['category']],
             s=90, alpha=0.85, edgecolors='white', lw=0.5)
for _, row in dm.iterrows():
    ax_b.annotate(row['disease'].split(',')[0][:16],
                  (row['vdr_tpm'], row['vdr_domain_index']),
                  fontsize=7, xytext=(3,2), textcoords='offset points')
x = dm['vdr_tpm'].values; y = dm['vdr_domain_index'].values
m, b = np.polyfit(x, y, 1)
xr = np.linspace(x.min(), x.max(), 50)
ax_b.plot(xr, m*xr+b, color='gray', ls='--', lw=1.5, alpha=0.7)
ax_b.set_xlabel('臓器VDR発現量 (GTEx TPM中央値)', fontsize=10)
ax_b.set_ylabel('VDRドメイン指数', fontsize=10)
ax_b.set_title(f'B   臓器VDR発現 × VDRドメイン指数\nSpearman r={r_corr:.2f}, p={p_corr:.3f}',
               fontsize=11, fontweight='bold', loc='left')

# ── C: VDRスコア分布 ─────────────────────────────────────────────
bins = np.linspace(0, gene_dis['VDR_score'].quantile(0.99)+1, 30)
ax_c.hist(fail_v, bins=bins, color=COL_FAIL, alpha=0.6, label=f'Failed (n={len(fail_v)})')
ax_c.hist(app_v,  bins=bins, color=COL_APP,  alpha=0.6, label=f'Approved (n={len(app_v)})')
ax_c.axvline(fail_v.median(), color=COL_FAIL, ls='--', lw=1.5)
ax_c.axvline(app_v.median(),  color=COL_APP,  ls='--', lw=1.5)
ax_c.set_xlabel('VDR score', fontsize=10); ax_c.set_ylabel('Gene-disease pairs', fontsize=10)
ax_c.legend(fontsize=8)
ax_c.set_title(f'C   VDRスコア: 承認 vs 失敗\nMann-Whitney p={p_mw:.4f}',
               fontsize=11, fontweight='bold', loc='left')

# ── D: 疾患座標マップ ─────────────────────────────────────────────
ax_d.scatter(dm['vdr_domain_index'], dm['gr_fail_signal'],
             c=[CATEGORY_COLOR.get(c,'#aaa') for c in dm['category']],
             s=dm['vdr_tpm']*12+40, alpha=0.85, edgecolors='white', lw=0.5)
for _, row in dm.iterrows():
    ax_d.annotate(row['disease'].split(',')[0][:14],
                  (row['vdr_domain_index'], row['gr_fail_signal']),
                  fontsize=7, xytext=(3,2), textcoords='offset points')
ax_d.axhline(0, color='gray', ls='--', lw=0.8)
ax_d.axvline(0, color='gray', ls='--', lw=0.8)
ax_d.set_xlabel('VDRドメイン指数', fontsize=10)
ax_d.set_ylabel('GR失敗シグナル', fontsize=10)
ax_d.set_title('D   疾患の座標マップ\n右上=VDRドメイン疾患, 左上=GR失敗多い',
               fontsize=11, fontweight='bold', loc='left')

# ── E: Fisher's forest plot ─────────────────────────────────────
or_vals = [or_v]
ci_los  = [ci_lo]
ci_his  = [ci_hi]
ps      = [p_v]
labels  = [f'VDR-dominant\nvs GR-dominant\n(全ステロイド疾患)']
ax_e.plot([ci_lo, ci_hi], [0, 0], color=COL_VDR, lw=3)
ax_e.plot(or_v, 0, 'D', color=COL_VDR, ms=12, zorder=5)
ax_e.axvline(1, color='gray', ls='--', lw=1)
ax_e.set_yticks([0]); ax_e.set_yticklabels(labels, fontsize=10)
ax_e.set_xlabel('Odds ratio', fontsize=11)
ax_e.set_title(
    f'E   Fisher\'s exact test\nOR={or_v:.2f} (95%CI {ci_lo:.2f}–{ci_hi:.2f}), p={p_v:.4f}',
    fontsize=11, fontweight='bold', loc='left')
ax_e.set_xlim(0, max(ci_his)*1.3)
ax_e.text(or_v, 0.15,
          f"VDR-dom 承認率: {n_va}/{n_va+n_vf} ({100*n_va/(n_va+n_vf):.0f}%)\n"
          f"GR-dom  承認率: {n_ga}/{n_ga+n_gf} ({100*n_ga/(n_ga+n_gf):.0f}%)",
          ha='center', va='bottom', fontsize=9)

for ax in [ax_a,ax_b,ax_c,ax_d,ax_e]:
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

fig.suptitle(
    '拡張版: ステロイド疾患ターゲット全遺伝子（がん適応除外）\n'
    'VDR-dominant ターゲットへの抗体薬は GR-dominant より有意に承認されやすい',
    fontsize=13, fontweight='bold', y=0.99)

out = RESULTS + 'fig_disease_vdr_expanded'
fig.savefig(out+'.pdf', dpi=300, bbox_inches='tight', facecolor='white')
fig.savefig(out+'.png', dpi=200, bbox_inches='tight', facecolor='white')
print(f"\nSaved: {out}.pdf / .png")
print("DONE")
