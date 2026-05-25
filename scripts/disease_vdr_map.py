#!/usr/bin/env python3
"""
疾患ごとのVDR/GRドメイン帰属マップ
- OpenTargets: 遺伝子 → 薬 → 疾患(phase) を逆引き
- ステロイド標準治療疾患に絞る
- 疾患ごとに「承認薬ターゲットのVDRスコア平均」を計算 → VDRドメイン指数
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
OT_URL  = "https://api.platform.opentargets.org/api/v4/graphql"

# ── VDR/GR スコア ───────────────────────────────────────────────────
df_vdr = pd.read_csv(RESULTS + 'remap_breadth_per_gene.csv')
gene2vdr  = dict(zip(df_vdr['gene'], df_vdr['VDR_score']))
gene2gr   = dict(zip(df_vdr['gene'], df_vdr['GR_score']))
gene2dom  = dict(zip(df_vdr['gene'], df_vdr['VDR_dominant']))
ALL_GENES = df_vdr['gene'].tolist()

# ── ステロイド疾患 EFO ID ────────────────────────────────────────────
# EFO/MONDO ID: OpenTargets で疾患を特定するために必要
STEROID_EFO = {
    'MONDO_0004979': ('Asthma',                     'respiratory'),
    'MONDO_0005002': ('COPD',                        'respiratory'),
    'EFO_0000676':   ('Rheumatoid Arthritis',        'autoimmune'),
    'MONDO_0007915': ('Lupus (SLE)',                 'autoimmune'),
    'EFO_0003885':   ('Psoriasis',                   'skin'),
    'EFO_0000676':   ('Rheumatoid Arthritis',        'autoimmune'),
    'MONDO_0005096': ('Psoriatic Arthritis',         'autoimmune'),
    'EFO_0000384':   ('Crohn Disease',               'GI'),
    'EFO_0000729':   ('Ulcerative Colitis',          'GI'),
    'EFO_0000274':   ('Atopic Dermatitis',           'skin'),
    'EFO_0003767':   ('Inflammatory Bowel Disease',  'GI'),
    'MONDO_0005301': ('Multiple Sclerosis',          'autoimmune'),
    'EFO_0000685':   ('Ankylosing Spondylitis',      'autoimmune'),
    'MONDO_0004975': ('Allergic Asthma',             'respiratory'),
}

CATEGORY_COLOR = {
    'respiratory': '#4393C3',
    'autoimmune':  '#D6604D',
    'skin':        '#F4A582',
    'GI':          '#74C476',
}

COL_VDR = '#2166AC'
COL_GR  = '#D6604D'
COL_APP = '#1a9641'
COL_FAIL= '#d7191c'

def ot_drugs_for_gene(gene_symbol):
    """遺伝子に対する薬剤+疾患情報をOpenTargetsから取得"""
    query = """
    query($sym: String!) {
      search(queryString: $sym, entityNames: ["target"]) {
        hits {
          object {
            ... on Target {
              approvedSymbol
              knownDrugs(size: 100) {
                rows {
                  phase
                  status
                  drug { name drugType }
                  disease { id name therapeuticAreas { id name } }
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
        with urllib.request.urlopen(req, timeout=20) as r:
            return json.loads(r.read())
    except:
        return None

# ── ステロイド疾患に含まれるか判定 ─────────────────────────────────
STEROID_KEYWORDS = {
    'asthma', 'chronic obstructive', 'copd', 'rhinitis',
    'rheumatoid', 'lupus', 'psoriasis', 'psoriatic',
    "crohn", 'colitis', 'inflammatory bowel',
    'atopic', 'dermatitis', 'eczema',
    'multiple sclerosis', 'ankylosing', 'spondylitis',
    'dermatomyositis', 'polymyositis', 'uveitis',
    'nephrotic', 'vasculitis', 'sarcoidosis',
}

def is_steroid_disease(disease_name):
    n = disease_name.lower()
    return any(kw in n for kw in STEROID_KEYWORDS)

def get_category(disease_name):
    n = disease_name.lower()
    if any(k in n for k in ['asthma','copd','obstructive','rhinitis','pulmonary','bronch']):
        return 'respiratory'
    if any(k in n for k in ['rheumatoid','lupus','sclerosis','ankylosing','spondylitis',
                             'dermatomyositis','uveitis','vasculitis','sarcoidosis']):
        return 'autoimmune'
    if any(k in n for k in ['psoriasis','psoriatic','atopic','dermatitis','eczema']):
        return 'skin'
    if any(k in n for k in ["crohn","colitis","bowel","ileitis"]):
        return 'GI'
    return 'other'

# ══════════════════════════════════════════════════════════════════════
# QUERY: 全308遺伝子 → 薬 → 疾患
# ══════════════════════════════════════════════════════════════════════
print("=== OpenTargets: 遺伝子→薬→疾患 クエリ ===")
rows = []
for i, gene in enumerate(ALL_GENES):
    res = ot_drugs_for_gene(gene)
    if res and res.get('data'):
        for hit in res['data']['search']['hits']:
            obj = hit.get('object', {})
            if obj.get('approvedSymbol','').upper() != gene.upper():
                continue
            for drug_row in obj.get('knownDrugs', {}).get('rows', []):
                dis = drug_row.get('disease')
                if not dis:
                    continue
                dis_name = dis.get('name', '')
                if not is_steroid_disease(dis_name):
                    continue
                # drugType filter: 抗体・生物製剤を優先タグ
                drug_info = drug_row.get('drug', {})
                drug_type = drug_info.get('drugType', '')
                is_biologic = drug_type in ('Antibody', 'Protein', 'Enzyme',
                                             'Oligonucleotide', 'Cell therapy',
                                             'Gene therapy', 'Oligosaccharide')
                rows.append({
                    'gene':       gene,
                    'disease':    dis_name,
                    'category':   get_category(dis_name),
                    'phase':      drug_row.get('phase', 0),
                    'status':     drug_row.get('status', ''),
                    'drug_name':  drug_info.get('name', ''),
                    'drug_type':  drug_type,
                    'is_biologic': is_biologic,
                    'VDR_score':  gene2vdr.get(gene, 0),
                    'GR_score':   gene2gr.get(gene, 0),
                    'VDR_dominant': gene2dom.get(gene, False),
                    'dis_id':     dis.get('id', ''),
                })
            break
    if (i+1) % 50 == 0:
        print(f"  {i+1}/{len(ALL_GENES)} done  ({len(rows)} records)")
    time.sleep(0.12)

df = pd.DataFrame(rows)
print(f"\n総レコード: {len(df)}")
print(f"疾患数: {df['disease'].nunique()}")
print(f"遺伝子数: {df['gene'].nunique()}")
df.to_csv(RESULTS + 'gene_disease_phase_ot.csv', index=False)

# ══════════════════════════════════════════════════════════════════════
# DISEASE MAP: 疾患ごとVDR/GRドメイン指数
# ══════════════════════════════════════════════════════════════════════
print("\n=== 疾患マップ計算 ===")
df['approved'] = (df['phase'] >= 4).astype(int)

# 疾患ごと:
#   vdr_domain_score = 承認薬ターゲットの平均VDRスコア
#   gr_failure_score = 失敗薬ターゲットの平均GRスコア（失敗 = phase<4 かつ phase>=1）
disease_map = []
for dis, grp in df.groupby('disease'):
    cat = grp['category'].iloc[0]
    approved  = grp[grp['phase'] >= 4]
    failed    = grp[(grp['phase'] >= 1) & (grp['phase'] < 4)]
    n_genes_app  = approved['gene'].nunique()
    n_genes_fail = failed['gene'].nunique()
    if n_genes_app + n_genes_fail < 3:
        continue

    vdr_app  = approved['VDR_score'].mean()  if len(approved)  else 0
    gr_app   = approved['GR_score'].mean()   if len(approved)  else 0
    vdr_fail = failed['VDR_score'].mean()    if len(failed)    else 0
    gr_fail  = failed['GR_score'].mean()     if len(failed)    else 0

    # VDRドメイン指数 = 承認薬のVDR平均 - 失敗薬のVDR平均
    vdr_domain_idx = vdr_app - vdr_fail
    # GR-failure signal = 失敗薬のGR平均 - 承認薬のGR平均
    gr_fail_signal = gr_fail - gr_app

    # 抗体薬限定でも計算
    bio_app  = approved[approved['is_biologic']]
    bio_fail = failed[failed['is_biologic']]
    vdr_bio_app  = bio_app['VDR_score'].mean()  if len(bio_app)  else np.nan
    gr_bio_fail  = bio_fail['GR_score'].mean()  if len(bio_fail) else np.nan

    disease_map.append({
        'disease':         dis,
        'category':        cat,
        'n_approved_genes': n_genes_app,
        'n_failed_genes':   n_genes_fail,
        'vdr_approved':    vdr_app,
        'gr_approved':     gr_app,
        'vdr_failed':      vdr_fail,
        'gr_failed':       gr_fail,
        'vdr_domain_index': vdr_domain_idx,
        'gr_failure_signal': gr_fail_signal,
        'vdr_bio_approved': vdr_bio_app,
        'gr_bio_failed':    gr_bio_fail,
    })

dm = pd.DataFrame(disease_map).sort_values('vdr_domain_index', ascending=False)
print(dm[['disease','category','n_approved_genes','n_failed_genes',
          'vdr_domain_index','gr_failure_signal']].to_string())
dm.to_csv(RESULTS + 'disease_vdr_gr_map.csv', index=False)

# ══════════════════════════════════════════════════════════════════════
# FIGURE: 疾患VDR/GRマップ
# ══════════════════════════════════════════════════════════════════════
fig = plt.figure(figsize=(18, 14))
gs  = GridSpec(2, 3, figure=fig, hspace=0.50, wspace=0.42)
ax_a = fig.add_subplot(gs[0, :2])   # 疾患マップ（横棒）
ax_b = fig.add_subplot(gs[0, 2])    # 疾患×VDR/GR scatter
ax_c = fig.add_subplot(gs[1, :])    # 遺伝子×疾患ヒートマップ

# ── Panel A: 疾患VDRドメイン指数ランキング ─────────────────────────
dm_sorted = dm.sort_values('vdr_domain_index')
colors_bar = [CATEGORY_COLOR.get(c, '#aaa') for c in dm_sorted['category']]
bars = ax_a.barh(range(len(dm_sorted)), dm_sorted['vdr_domain_index'],
                  color=colors_bar, alpha=0.85, height=0.7)
ax_a.set_yticks(range(len(dm_sorted)))
ax_a.set_yticklabels(dm_sorted['disease'], fontsize=10)
ax_a.axvline(0, color='black', lw=1)
ax_a.set_xlabel('VDRドメイン指数\n(承認薬ターゲットの平均VDRスコア − 失敗薬ターゲットの平均VDRスコア)',
                fontsize=10)
ax_a.set_title('A   疾患ごとのVDRドメイン帰属指数\n(正 = VDRドメイン疾患 → VDR-dominant ターゲットが正解)',
               fontsize=11, fontweight='bold', loc='left')

# カテゴリ凡例
handles = [mpatches.Patch(color=c, label=k)
           for k, c in CATEGORY_COLOR.items() if k != 'other']
ax_a.legend(handles=handles, fontsize=9, loc='lower right')

# ── Panel B: 疾患scatter (VDR承認 vs GR失敗) ─────────────────────────
sc = ax_b.scatter(dm['vdr_approved'], dm['gr_failed'],
                  c=[CATEGORY_COLOR.get(c,'#aaa') for c in dm['category']],
                  s=dm['n_approved_genes']*15+30, alpha=0.8, edgecolors='white', lw=0.5)
for _, row in dm.iterrows():
    ax_b.annotate(row['disease'].split(',')[0].split('(')[0].strip()[:20],
                  (row['vdr_approved'], row['gr_failed']),
                  fontsize=7, xytext=(3, 2), textcoords='offset points')
ax_b.set_xlabel('承認薬ターゲットの平均VDRスコア', fontsize=10)
ax_b.set_ylabel('失敗薬ターゲットの平均GRスコア', fontsize=10)
ax_b.set_title('B   疾患座標\n(右上 = VDR承認×GR失敗 = 典型的VDRドメイン疾患)',
               fontsize=11, fontweight='bold', loc='left')
ax_b.spines['top'].set_visible(False)
ax_b.spines['right'].set_visible(False)

# ── Panel C: 遺伝子×疾患ヒートマップ（VDR-dominant 遺伝子） ─────────
vdr_dom_genes = df[df['VDR_dominant']==True]['gene'].unique()
diseases_list  = dm['disease'].tolist()

# 各遺伝子の各疾患でのmax phase
heat_data = []
for g in vdr_dom_genes:
    row_data = {'gene': g, 'VDR': gene2vdr.get(g,0)}
    gdf = df[df['gene']==g]
    for dis in diseases_list:
        ddf = gdf[gdf['disease']==dis]
        row_data[dis] = ddf['phase'].max() if len(ddf) else np.nan
    heat_data.append(row_data)

heat_df = (pd.DataFrame(heat_data)
           .set_index('gene')
           .sort_values('VDR', ascending=False)
           .drop(columns='VDR'))

# 試験データがある遺伝子のみ
heat_df = heat_df.dropna(how='all')
if len(heat_df) > 25:
    heat_df = heat_df.head(25)

mat = heat_df.values.astype(float)
im  = ax_c.imshow(mat, aspect='auto', cmap='YlOrRd', vmin=0, vmax=4,
                  interpolation='nearest')

ax_c.set_xticks(range(len(diseases_list)))
ax_c.set_xticklabels([d.split(',')[0][:22] for d in diseases_list],
                      rotation=40, ha='right', fontsize=8)
ylabels = [f"▲ {g}" if gene2dom.get(g, False) else g for g in heat_df.index]
ax_c.set_yticks(range(len(heat_df)))
ax_c.set_yticklabels(ylabels, fontsize=8)
plt.colorbar(im, ax=ax_c, label='Max clinical phase', shrink=0.6, aspect=15)
ax_c.set_title('C   VDR-dominant 遺伝子 × ステロイド疾患マップ (max phase, ▲=VDR-dominant)',
               fontsize=11, fontweight='bold', loc='left')

for ax in [ax_a, ax_b]:
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

fig.suptitle('疾患ごとのVDR/GRドメイン帰属マップ\n'
             '「ステロイドが届かない場所」に VDR-dominant ターゲットへの抗体薬が効く',
             fontsize=13, fontweight='bold', y=0.99)

out = RESULTS + 'fig_disease_vdr_gr_map'
fig.savefig(out + '.pdf', dpi=300, bbox_inches='tight', facecolor='white')
fig.savefig(out + '.png', dpi=200, bbox_inches='tight', facecolor='white')
print(f"\nSaved: {out}.pdf / .png")

# ── コンソールサマリー ─────────────────────────────────────────────
print("\n=== 疾患VDRドメイン指数ランキング ===")
print(dm[['disease','category','vdr_domain_index','gr_failure_signal',
          'n_approved_genes','n_failed_genes']].to_string(index=False))
print("\nDONE")
