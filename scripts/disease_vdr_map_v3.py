#!/usr/bin/env python3
"""
疾患VDR/GRドメインマップ v3
「本当に失敗した薬」の定義:
  - withdrawn_flag=True (安全性で市場撤退), または
  - グローバルmax_phase <= 1 (Phase1を超えられなかった), または
  - グローバルmax_phase 2-3 かつ first_approval=None (承認なし、開発停止推定)
「進行中」を除外して信号をクリーンにする
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

def chembl_get(url, retries=3):
    for i in range(retries):
        try:
            with urllib.request.urlopen(url, timeout=20) as r:
                return json.loads(r.read())
        except Exception:
            if i < retries-1: time.sleep(1)
    return None

# ── 既存データ読み込み ───────────────────────────────────────────────
ind      = pd.read_csv(RESULTS + 'steroid_disease_indications.csv')
df_gene  = pd.read_csv(RESULTS + 'gene_disease_phase_chembl.csv')
df_vdr   = pd.read_csv(RESULTS + 'remap_breadth_per_gene.csv')
gene2vdr = dict(zip(df_vdr['gene'], df_vdr['VDR_score']))
gene2gr  = dict(zip(df_vdr['gene'], df_vdr['GR_score']))
gene2dom = dict(zip(df_vdr['gene'], df_vdr['VDR_dominant']))

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
    'Asthma, Allergic':                      ('lung',          18.6),
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

# ══════════════════════════════════════════════════════════════════════
# STEP 1: 全分子の withdrawn_flag + global max_phase + first_approval 取得
# ══════════════════════════════════════════════════════════════════════
unique_mols = ind['mol_id'].unique()
print(f"分子数: {len(unique_mols)}")
print("=== 分子ステータス取得 ===")

mol_status = {}
BATCH = 50
for i in range(0, len(unique_mols), BATCH):
    batch = unique_mols[i:i+BATCH]
    ids   = ','.join(batch)
    url   = (f"{CHEMBL}/molecule.json?molecule_chembl_id__in={ids}"
             f"&limit=500&format=json")
    d = chembl_get(url)
    if not d:
        continue
    for mol in d.get('molecules', []):
        mid = mol['molecule_chembl_id']
        mol_status[mid] = {
            'global_max_phase':  float(mol.get('max_phase') or 0),
            'withdrawn':         bool(mol.get('withdrawn_flag', False)),
            'first_approval':    mol.get('first_approval'),
        }
    if (i // BATCH + 1) % 5 == 0:
        print(f"  {i+BATCH}/{len(unique_mols)} done")
    time.sleep(0.12)

print(f"取得: {len(mol_status)} 分子")

# ── ステータス分類 ─────────────────────────────────────────────────
# approved:      max_phase_for_ind = 4
# truly_failed:  withdrawn=True OR (global_max_phase<4 AND first_approval=None AND phase_for_ind>=2)
# p1_stalled:    global_max_phase<=1 (phase1超えられず)
# ongoing:       global_max_phase>=2 AND first_approval=None AND not withdrawn (進行中)
def classify_mol(mol_id, phase_for_ind):
    s = mol_status.get(mol_id, {})
    gp = s.get('global_max_phase', 0)
    wd = s.get('withdrawn', False)
    fa = s.get('first_approval')

    if phase_for_ind >= 4:
        return 'approved'
    if wd:
        return 'failed'
    if gp <= 1 and phase_for_ind >= 1:
        return 'failed'          # Phase1すら超えられなかった
    if gp >= 2 and fa is None and phase_for_ind >= 2:
        return 'failed'          # Phase2+まで進んだが承認なし → 失敗
    if gp >= 2 and fa is not None and phase_for_ind < 4:
        return 'approved_other'  # 別適応症で承認済み（この疾患では不採用）
    return 'ongoing'             # 進行中 → 除外

ind['mol_status'] = ind.apply(
    lambda r: classify_mol(r['mol_id'], r['max_phase']), axis=1)

print("\n分子ステータス分布:")
print(ind['mol_status'].value_counts())

# ── 遺伝子×疾患テーブルにステータスをマージ ─────────────────────
mol2status = dict(zip(ind['mol_id'], ind['mol_status']))
df_gene['mol_status'] = df_gene['mol_id'].map(mol2status)
df_gene = df_gene[df_gene['mol_status'].isin(['approved','failed'])]
print(f"\n承認+失敗レコード: {len(df_gene)}")
print(f"  承認: {(df_gene['mol_status']=='approved').sum()}")
print(f"  失敗: {(df_gene['mol_status']=='failed').sum()}")

# ══════════════════════════════════════════════════════════════════════
# STEP 2: 疾患VDRドメイン指数（クリーン版）
# ══════════════════════════════════════════════════════════════════════
print("\n=== 疾患VDRドメイン指数（クリーン版） ===")

# 遺伝子×疾患ごとに最良ステータスを集約
gene_dis = (df_gene.groupby(['gene','disease','category'])
            .agg(
                best_status = ('mol_status', lambda x:
                               'approved' if 'approved' in x.values else 'failed'),
                VDR_score = ('VDR_score','first'),
                GR_score  = ('GR_score','first'),
                VDR_dom   = ('VDR_dominant','first'),
            ).reset_index())

# Fisher's exact: VDR-dominant × 承認 (全疾患まとめ)
ct = pd.crosstab(gene_dis['VDR_dom'], gene_dis['best_status'])
print(f"\nFisher's (全疾患):\n{ct}")
if ct.shape == (2,2) and 'approved' in ct.columns and 'failed' in ct.columns:
    n_va = ct.loc[True,  'approved']; n_vf = ct.loc[True,  'failed']
    n_ga = ct.loc[False, 'approved']; n_gf = ct.loc[False, 'failed']
    or_v, p_v = stats.fisher_exact([[n_va,n_vf],[n_ga,n_gf]], alternative='greater')
    print(f"OR={or_v:.2f}  p={p_v:.4f}")

# Spearman
app_v  = gene_dis[gene_dis['best_status']=='approved']['VDR_score']
fail_v = gene_dis[gene_dis['best_status']=='failed']['VDR_score']
_, p_mw = stats.mannwhitneyu(app_v, fail_v, alternative='greater')
print(f"VDRスコア 承認 vs 失敗 Mann-Whitney p={p_mw:.4f}")
print(f"承認中央値={app_v.median():.1f}  失敗中央値={fail_v.median():.1f}")

disease_map = []
for dis, grp in gene_dis.groupby('disease'):
    info = DISEASE_ORGAN.get(dis, ('immune cells', 15.4))
    cat  = grp['category'].iloc[0]
    app  = grp[grp['best_status']=='approved']
    fail = grp[grp['best_status']=='failed']
    if len(app) + len(fail) < 3:
        continue
    disease_map.append({
        'disease':           dis,
        'category':          cat,
        'organ':             info[0],
        'vdr_tpm':           info[1],
        'n_app':             len(app),
        'n_fail':            len(fail),
        'vdr_approved':      app['VDR_score'].mean() if len(app) else 0,
        'gr_approved':       app['GR_score'].mean()  if len(app) else 0,
        'vdr_failed':        fail['VDR_score'].mean() if len(fail) else 0,
        'gr_failed':         fail['GR_score'].mean()  if len(fail) else 0,
        'vdr_domain_index':  (app['VDR_score'].mean() if len(app) else 0) -
                             (fail['VDR_score'].mean() if len(fail) else 0),
        'gr_fail_signal':    (fail['GR_score'].mean() if len(fail) else 0) -
                             (app['GR_score'].mean()  if len(app) else 0),
        'pct_vdr_app':       app['VDR_dom'].mean() if len(app) else 0,
    })

dm = pd.DataFrame(disease_map).sort_values('vdr_domain_index', ascending=False)
dm.to_csv(RESULTS + 'disease_vdr_gr_map_v3.csv', index=False)
print("\n疾患VDRドメイン指数:")
print(dm[['disease','organ','n_app','n_fail',
          'vdr_domain_index','gr_fail_signal','vdr_tpm']].to_string(index=False))

r_corr, p_corr = stats.spearmanr(dm['vdr_tpm'], dm['vdr_domain_index'])
print(f"\n臓器VDR発現 vs VDRドメイン指数: r={r_corr:.3f}  p={p_corr:.3f}")

# ══════════════════════════════════════════════════════════════════════
# FIGURE
# ══════════════════════════════════════════════════════════════════════
fig = plt.figure(figsize=(20, 14))
gs  = GridSpec(2, 3, figure=fig, hspace=0.48, wspace=0.42)
ax_a = fig.add_subplot(gs[0, :2])
ax_b = fig.add_subplot(gs[0, 2])
ax_c = fig.add_subplot(gs[1, 0])
ax_d = fig.add_subplot(gs[1, 1])
ax_e = fig.add_subplot(gs[1, 2])

COL_APP='#1a9641'; COL_FAIL='#d7191c'

# ── A: 疾患VDRドメイン指数（バブル=臓器VDR発現） ────────────────
dm_s = dm.sort_values('vdr_domain_index')
col_a = [CATEGORY_COLOR.get(c,'#aaa') for c in dm_s['category']]
sz_a  = (dm_s['vdr_tpm'] / dm_s['vdr_tpm'].max() * 700 + 100).values

y_pos = range(len(dm_s))
ax_a.scatter(dm_s['vdr_domain_index'], y_pos,
             s=sz_a, c=col_a, alpha=0.85, zorder=3)
ax_a.axvline(0, color='black', lw=1, alpha=0.5)
ax_a.set_yticks(y_pos)
ax_a.set_yticklabels(
    [f"{r['disease'].split(',')[0]}  [{r['organ']}]"
     for _, r in dm_s.iterrows()], fontsize=10)
ax_a.set_xlabel(
    'VDRドメイン指数\n(承認薬ターゲット平均VDRスコア − 失敗薬ターゲット平均VDRスコア)',
    fontsize=10)
ax_a.set_title(
    'A   疾患ごとのVDRドメイン帰属指数 (クリーン版)\n'
    'バブルサイズ = 臓器VDR発現量(GTEx TPM)  |  進行中薬剤は除外',
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
ax_b.set_title(
    f'B   臓器VDR発現 × VDRドメイン指数\nSpearman r={r_corr:.2f}, p={p_corr:.3f}',
    fontsize=11, fontweight='bold', loc='left')

# ── C: VDRスコア分布（承認 vs 失敗） ────────────────────────────
bins = np.linspace(0, gene_dis['VDR_score'].max()+1, 25)
ax_c.hist(fail_v, bins=bins, color=COL_FAIL, alpha=0.6,
          label=f'Failed (n={len(fail_v)})')
ax_c.hist(app_v,  bins=bins, color=COL_APP,  alpha=0.6,
          label=f'Approved (n={len(app_v)})')
ax_c.axvline(fail_v.median(), color=COL_FAIL, ls='--', lw=1.5)
ax_c.axvline(app_v.median(),  color=COL_APP,  ls='--', lw=1.5)
ax_c.set_xlabel('VDR score', fontsize=10)
ax_c.set_ylabel('Gene-disease pairs', fontsize=10)
ax_c.legend(fontsize=8)
ax_c.set_title(f'C   VDRスコア: 承認 vs 失敗\n(Mann-Whitney p={p_mw:.4f})',
               fontsize=11, fontweight='bold', loc='left')

# ── D: GR失敗シグナル × VDRドメイン指数 疾患scatter ─────────────
sc = ax_d.scatter(dm['vdr_domain_index'], dm['gr_fail_signal'],
                  c=[CATEGORY_COLOR.get(c,'#aaa') for c in dm['category']],
                  s=dm['vdr_tpm']*10+30, alpha=0.85,
                  edgecolors='white', lw=0.5)
for _, row in dm.iterrows():
    ax_d.annotate(row['disease'].split(',')[0][:16],
                  (row['vdr_domain_index'], row['gr_fail_signal']),
                  fontsize=7, xytext=(3,2), textcoords='offset points')
ax_d.axhline(0, color='gray', ls='--', lw=0.8)
ax_d.axvline(0, color='gray', ls='--', lw=0.8)
ax_d.set_xlabel('VDRドメイン指数\n(VDR承認↑)', fontsize=10)
ax_d.set_ylabel('GR失敗シグナル\n(GR失敗↑)', fontsize=10)
ax_d.set_title('D   疾患の座標マップ\n右上 = VDR承認かつGR失敗 = 典型的VDRドメイン疾患',
               fontsize=11, fontweight='bold', loc='left')
# 象限ラベル
xlim = ax_d.get_xlim(); ylim = ax_d.get_ylim()
ax_d.text(xlim[1]*0.6, ylim[1]*0.85, 'VDRドメイン疾患', fontsize=8,
          color='#2166AC', alpha=0.7, fontweight='bold')
ax_d.text(xlim[0]*0.9, ylim[1]*0.85, 'GRドメイン疾患\n(ステロイドで足りる)',
          fontsize=8, color='#D6604D', alpha=0.7)

# ── E: VDR-dominant 「次のボタン」 top15 ─────────────────────────
next_btn = (gene_dis[gene_dis['VDR_dom']==True]
            .sort_values('VDR_score', ascending=False)
            .drop_duplicates('gene').head(15))
col_e = [COL_APP if r=='approved' else COL_FAIL
         for r in next_btn['best_status']]
ax_e.barh(range(len(next_btn)), next_btn['VDR_score'].values[::-1],
          color=col_e[::-1], alpha=0.85, height=0.7)
ax_e.set_yticks(range(len(next_btn)))
ax_e.set_yticklabels(next_btn['gene'].values[::-1], fontsize=9)
ax_e.set_xlabel('VDR ChIP-seq score', fontsize=10)
ax_e.set_title('E   「次のボタン」候補\n(VDR-dominant, ▲=ステロイド疾患試験あり)',
               fontsize=11, fontweight='bold', loc='left')
ax_e.legend(handles=[mpatches.Patch(facecolor=COL_APP, label='承認済'),
                      mpatches.Patch(facecolor=COL_FAIL, label='失敗/未承認')],
            fontsize=8)

for ax in [ax_a,ax_b,ax_c,ax_d,ax_e]:
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

fig.suptitle(
    '臓器VDR依存性が高い疾患ほどVDRドメイン指数が高い\n'
    '「ステロイドが届かない場所の異常」に VDR-dominant ターゲットへの抗体薬が効く',
    fontsize=13, fontweight='bold', y=0.99)

out = RESULTS + 'fig_disease_vdr_organ_map_v3'
fig.savefig(out+'.pdf', dpi=300, bbox_inches='tight', facecolor='white')
fig.savefig(out+'.png', dpi=200, bbox_inches='tight', facecolor='white')
print(f"\nSaved: {out}.pdf / .png")
print("DONE")
