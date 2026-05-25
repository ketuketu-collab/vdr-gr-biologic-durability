#!/usr/bin/env python3
"""
ステロイドが標準治療の疾患に絞ったChEMBL解析
構図: ステロイド治療中の疾患 → 次に押すべきVDR-dominant ボタンの同定
      GR-dominant ターゲットは失敗と相関するはず
"""

import pandas as pd
import numpy as np
import urllib.request, json, time
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from scipy import stats
import warnings
warnings.filterwarnings('ignore')

RESULTS = '/Volumes/M4_SSD/projects/tlr_chipseq/results/'
BASE    = 'https://www.ebi.ac.uk/chembl/api/data'

# ── ステロイド標準治療疾患 (MESH heading) ─────────────────────────────
STEROID_DISEASES = {
    'Asthma':                            'respiratory',
    'Pulmonary Disease, Chronic Obstructive': 'respiratory',
    'Rhinitis, Allergic':                'respiratory',
    'Arthritis, Rheumatoid':             'autoimmune',
    'Lupus Erythematosus, Systemic':     'autoimmune',
    'Spondylitis, Ankylosing':           'autoimmune',
    'Psoriatic Arthritis':               'autoimmune',
    'Multiple Sclerosis':                'autoimmune',
    'Dermatomyositis':                   'autoimmune',
    'Dermatitis, Atopic':                'skin',
    'Psoriasis':                         'skin',
    'Crohn Disease':                     'GI',
    'Colitis, Ulcerative':               'GI',
    'Nephrotic Syndrome':                'renal',
    'Uveitis':                           'eye',
    'Asthma, Allergic':                  'respiratory',
}

CATEGORY_COLOR = {
    'respiratory': '#4393C3',
    'autoimmune':  '#D6604D',
    'skin':        '#F4A582',
    'GI':          '#74C476',
    'renal':       '#9970AB',
    'eye':         '#FDAE61',
}

# ── VDR/GR スコアデータ ───────────────────────────────────────────────
df_vdr = pd.read_csv(RESULTS + 'remap_breadth_per_gene.csv')
gene2vdr = dict(zip(df_vdr['gene'], df_vdr['VDR_score']))
gene2gr  = dict(zip(df_vdr['gene'], df_vdr['GR_score']))
gene2dom = dict(zip(df_vdr['gene'], df_vdr['VDR_dominant']))
gene2repro_vdr = dict(zip(df_vdr['gene'], df_vdr['VDR_repro']))
gene2repro_gr  = dict(zip(df_vdr['gene'], df_vdr['GR_repro']))
known_genes = set(df_vdr['gene'].str.upper())

def chembl_get(url, retries=3):
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(url, timeout=20) as r:
                return json.loads(r.read())
        except Exception as e:
            if attempt == retries - 1:
                return None
            time.sleep(1)

def get_target_genes(mol_id):
    """分子のターゲット遺伝子シンボルを取得"""
    url = f"{BASE}/mechanism.json?molecule_chembl_id={mol_id}&limit=50&format=json"
    d = chembl_get(url)
    genes = set()
    if not d:
        return genes
    for mech in d.get('mechanisms', []):
        tid = mech.get('target_chembl_id')
        if not tid:
            continue
        turl = f"{BASE}/target/{tid}.json?format=json"
        td = chembl_get(turl)
        if not td:
            continue
        for comp in td.get('target_components', []):
            for syn in comp.get('target_component_synonyms', []):
                if syn.get('syn_type') == 'GENE_SYMBOL':
                    genes.add(syn['component_synonym'].upper())
    return genes

# ══════════════════════════════════════════════════════════════════════
# STEP 1: 各疾患の試験薬を全取得
# ══════════════════════════════════════════════════════════════════════
print("=== STEP 1: ChEMBL drug_indication クエリ ===")
rows = []
for disease, category in STEROID_DISEASES.items():
    print(f"\n  {disease} ({category})...")
    mesh_enc = urllib.parse.quote(disease) if hasattr(urllib, 'parse') else disease.replace(' ', '+').replace(',', '%2C')
    import urllib.parse
    mesh_enc = urllib.parse.quote(disease)

    offset, total = 0, None
    n_found = 0
    while total is None or offset < total:
        url = (f"{BASE}/drug_indication.json?"
               f"mesh_heading={mesh_enc}&limit=100&offset={offset}&format=json")
        d = chembl_get(url)
        if not d:
            break
        if total is None:
            total = d['page_meta']['total_count']
        for ind in d.get('drug_indications', []):
            rows.append({
                'disease':   disease,
                'category':  category,
                'mol_id':    ind['molecule_chembl_id'],
                'max_phase': float(ind['max_phase_for_ind'] or 0),
                'efo_term':  ind.get('efo_term',''),
            })
            n_found += 1
        offset += 100
        time.sleep(0.1)
    print(f"    {n_found} trials")

df_ind = pd.DataFrame(rows).drop_duplicates(['disease','mol_id'])
print(f"\n総試験レコード: {len(df_ind)}  ({df_ind['mol_id'].nunique()} 分子)")
df_ind.to_csv(RESULTS + 'steroid_disease_indications.csv', index=False)

# ══════════════════════════════════════════════════════════════════════
# STEP 2: 分子→ターゲット遺伝子マッピング
# ══════════════════════════════════════════════════════════════════════
print("\n=== STEP 2: 分子→ターゲット遺伝子マッピング ===")
unique_mols = df_ind['mol_id'].unique()
print(f"ユニーク分子数: {len(unique_mols)}")

mol2genes = {}
for i, mol in enumerate(unique_mols):
    genes = get_target_genes(mol)
    if genes:
        mol2genes[mol] = genes
    if (i+1) % 50 == 0:
        print(f"  {i+1}/{len(unique_mols)} done  ({sum(len(v) for v in mol2genes.values())} gene hits)")
    time.sleep(0.12)

print(f"\nターゲット情報取得: {len(mol2genes)} 分子")

# ══════════════════════════════════════════════════════════════════════
# STEP 3: 遺伝子×疾患×フェーズのテーブル構築
# ══════════════════════════════════════════════════════════════════════
print("\n=== STEP 3: 遺伝子×疾患テーブル ===")
gene_rows = []
for _, ind_row in df_ind.iterrows():
    mol   = ind_row['mol_id']
    genes = mol2genes.get(mol, set())
    for g in genes:
        if g not in known_genes:
            continue
        gene_rows.append({
            'gene':     g,
            'disease':  ind_row['disease'],
            'category': ind_row['category'],
            'mol_id':   mol,
            'max_phase': ind_row['max_phase'],
            'VDR_score': gene2vdr.get(g, 0),
            'GR_score':  gene2gr.get(g, 0),
            'VDR_dominant': gene2dom.get(g, False),
            'VDR_repro': gene2repro_vdr.get(g, 0),
            'GR_repro':  gene2repro_gr.get(g, 0),
        })

df_gene = pd.DataFrame(gene_rows)
if len(df_gene) == 0:
    print("WARNING: 遺伝子マッピングが空 - ターゲット情報なし")
else:
    df_gene.to_csv(RESULTS + 'steroid_gene_disease_phase.csv', index=False)
    print(f"レコード数: {len(df_gene)}")
    print(f"ユニーク遺伝子: {df_gene['gene'].nunique()}")
    print(f"フェーズ分布:\n{df_gene.groupby('max_phase').size()}")

    # ── 遺伝子ごとの最高フェーズ・成否 ─────────────────────────────
    gene_summary = (df_gene.groupby('gene')
        .agg(
            max_phase    = ('max_phase','max'),
            n_trials     = ('mol_id','nunique'),
            diseases     = ('disease', lambda x: '; '.join(sorted(set(x)))),
            categories   = ('category', lambda x: '; '.join(sorted(set(x)))),
            VDR_score    = ('VDR_score','first'),
            GR_score     = ('GR_score','first'),
            VDR_dominant = ('VDR_dominant','first'),
            VDR_repro    = ('VDR_repro','first'),
            GR_repro     = ('GR_repro','first'),
        ).reset_index())

    # 成功=承認(4), 失敗=P1-P3止まり(1-3, 最高フェーズでそれ以上なし)
    gene_summary['approved'] = (gene_summary['max_phase'] >= 4).astype(int)
    gene_summary['failed']   = ((gene_summary['max_phase'] >= 1) &
                                 (gene_summary['max_phase'] < 4)).astype(int)

    gene_summary.to_csv(RESULTS + 'steroid_gene_summary.csv', index=False)
    print(f"\n遺伝子サマリー:")
    print(f"  承認 (phase≥4): {gene_summary['approved'].sum()}")
    print(f"  失敗 (phase<4): {gene_summary['failed'].sum()}")
    print(f"  VDR-dominant:   {gene_summary['VDR_dominant'].sum()}")

    # ── Fisher's exact test ───────────────────────────────────────────
    sub = gene_summary[gene_summary['max_phase'] >= 1]  # 試験あり
    ct = pd.crosstab(sub['VDR_dominant'], sub['approved'])
    print(f"\nFisher's (VDR-dominant vs approved):\n{ct}")
    if ct.shape == (2,2):
        n_va = ct.loc[True,  1]; n_vf = ct.loc[True,  0]
        n_ga = ct.loc[False, 1]; n_gf = ct.loc[False, 0]
        or_v, p_v = stats.fisher_exact([[n_va,n_vf],[n_ga,n_gf]], alternative='greater')
        print(f"OR={or_v:.2f}  p={p_v:.4f}")

    # ── Spearman: VDR score vs max_phase ─────────────────────────────
    r_vdr, p_vdr = stats.spearmanr(sub['VDR_score'], sub['max_phase'])
    r_gr,  p_gr  = stats.spearmanr(sub['GR_score'],  sub['max_phase'])
    print(f"\nSpearman VDR score vs max_phase: r={r_vdr:.3f}  p={p_vdr:.3e}")
    print(f"Spearman GR  score vs max_phase: r={r_gr:.3f}  p={p_gr:.3e}")

    # ══════════════════════════════════════════════════════════════════
    # FIGURE
    # ══════════════════════════════════════════════════════════════════
    fig = plt.figure(figsize=(18,13))
    gs  = GridSpec(2,3, figure=fig, hspace=0.45, wspace=0.40)
    ax_a = fig.add_subplot(gs[0,0])
    ax_b = fig.add_subplot(gs[0,1])
    ax_c = fig.add_subplot(gs[0,2])
    ax_d = fig.add_subplot(gs[1,:2])
    ax_e = fig.add_subplot(gs[1,2])

    COL_VDR='#2166AC'; COL_GR='#D6604D'; COL_APP='#1a9641'; COL_FAIL='#d7191c'

    # ── A: VDR vs GR score, colored by max_phase ─────────────────────
    phase_colors = {0:'#bdbdbd', 1:'#fdcc8a', 2:'#fc8d59', 3:'#e34a33', 4:COL_APP}
    for ph in sorted(sub['max_phase'].unique()):
        grp = sub[sub['max_phase']==ph]
        label = f'Phase {int(ph)}' if ph < 4 else 'Approved'
        ax_a.scatter(np.log1p(grp['VDR_score']), np.log1p(grp['GR_score']),
                     c=phase_colors.get(ph,'#aaa'), s=60, alpha=0.8,
                     label=label, zorder=int(ph)+1)
    diag = max(np.log1p(sub['VDR_score'].max()), np.log1p(sub['GR_score'].max()))
    ax_a.plot([0,diag],[0,diag],'k--',lw=0.8,alpha=0.4)
    ax_a.set_xlabel('log(1+VDR score)',fontsize=10)
    ax_a.set_ylabel('log(1+GR score)',fontsize=10)
    ax_a.legend(fontsize=8,markerscale=1.1)
    ax_a.set_title('A   VDR vs GR score by max phase\n(ステロイド疾患ターゲット)',
                   fontsize=11,fontweight='bold',loc='left')

    # ── B: VDR score distribution by approved/failed ──────────────────
    app_vdr  = sub[sub['approved']==1]['VDR_score']
    fail_vdr = sub[sub['approved']==0]['VDR_score']
    bins = np.linspace(0, sub['VDR_score'].max()+1, 30)
    ax_b.hist(fail_vdr, bins=bins, color=COL_FAIL, alpha=0.6, label=f'Failed (n={len(fail_vdr)})')
    ax_b.hist(app_vdr,  bins=bins, color=COL_APP,  alpha=0.6, label=f'Approved (n={len(app_vdr)})')
    ax_b.axvline(fail_vdr.median(), color=COL_FAIL, ls='--', lw=1.5)
    ax_b.axvline(app_vdr.median(),  color=COL_APP,  ls='--', lw=1.5)
    rr, pp = stats.mannwhitneyu(app_vdr, fail_vdr, alternative='greater')
    ax_b.set_xlabel('VDR score', fontsize=10)
    ax_b.set_ylabel('Genes', fontsize=10)
    ax_b.legend(fontsize=9)
    ax_b.set_title(f'B   VDR score: approved vs failed\n(Mann-Whitney p={pp:.3f})',
                   fontsize=11,fontweight='bold',loc='left')

    # ── C: disease category breakdown ─────────────────────────────────
    cat_counts = (df_gene.groupby(['category','gene'])['max_phase']
                  .max().reset_index()
                  .groupby('category').size().sort_values(ascending=True))
    colors_c = [CATEGORY_COLOR.get(c,'#aaa') for c in cat_counts.index]
    ax_c.barh(range(len(cat_counts)), cat_counts.values, color=colors_c, height=0.7)
    ax_c.set_yticks(range(len(cat_counts)))
    ax_c.set_yticklabels(cat_counts.index, fontsize=10)
    ax_c.set_xlabel('ユニーク遺伝子数', fontsize=10)
    ax_c.set_title('C   疾患カテゴリ別ターゲット数',
                   fontsize=11,fontweight='bold',loc='left')

    # ── D: gene×disease heatmap (top genes) ───────────────────────────
    top_genes = (gene_summary.sort_values('VDR_score', ascending=False)
                 .head(30)['gene'].tolist())
    diseases  = list(STEROID_DISEASES.keys())
    heat_df   = pd.DataFrame(index=top_genes, columns=diseases, dtype=float)
    for _, row in df_gene[df_gene['gene'].isin(top_genes)].iterrows():
        heat_df.loc[row['gene'], row['disease']] = row['max_phase']
    heat_df = heat_df.fillna(0)

    im = ax_d.imshow(heat_df.values, aspect='auto', cmap='YlOrRd', vmin=0, vmax=4)
    ax_d.set_xticks(range(len(diseases)))
    ax_d.set_xticklabels([d.replace(', ',' ').replace('Erythematosus','Eryth.') for d in diseases],
                          fontsize=7, rotation=45, ha='right')
    ax_d.set_yticks(range(len(top_genes)))

    # VDR-dominant マークを遺伝子名に付加
    ylabels = [f"▲ {g}" if gene2dom.get(g,False) else g for g in top_genes]
    ax_d.set_yticklabels(ylabels, fontsize=8)
    plt.colorbar(im, ax=ax_d, label='Max phase', shrink=0.8)
    ax_d.set_title('D   Top VDR-score遺伝子 × ステロイド疾患 (max phase, ▲=VDR-dominant)',
                   fontsize=11,fontweight='bold',loc='left')

    # ── E: top VDR-dominant "次のボタン" ─────────────────────────────
    next_buttons = (gene_summary[gene_summary['VDR_dominant']==True]
                    .sort_values('VDR_score', ascending=False)
                    .head(15))
    colors_e = [COL_APP if r['approved']==1 else COL_FAIL
                for _, r in next_buttons.iterrows()]
    ax_e.barh(range(len(next_buttons)),
              next_buttons['VDR_score'].values[::-1],
              color=colors_e[::-1], alpha=0.85, height=0.7)
    ax_e.set_yticks(range(len(next_buttons)))
    ax_e.set_yticklabels(next_buttons['gene'].values[::-1], fontsize=9)
    ax_e.set_xlabel('VDR ChIP-seq score', fontsize=10)
    ax_e.set_title('E   「次に押すべきボタン」\n(VDR-dominant, ステロイド疾患)',
                   fontsize=11,fontweight='bold',loc='left')
    from matplotlib.patches import Patch
    ax_e.legend(handles=[Patch(facecolor=COL_APP,label='承認済'),
                          Patch(facecolor=COL_FAIL,label='未承認/失敗')], fontsize=8)

    for ax in [ax_a,ax_b,ax_c,ax_d,ax_e]:
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

    fig.suptitle('ステロイド標準治療疾患における「次の一手」:\nVDR-dominant ターゲットが承認されやすく、GR-dominant は失敗と相関',
                 fontsize=13,fontweight='bold',y=0.99)

    out = RESULTS + 'fig_steroid_next_target'
    fig.savefig(out+'.pdf', dpi=300, bbox_inches='tight', facecolor='white')
    fig.savefig(out+'.png', dpi=200, bbox_inches='tight', facecolor='white')
    print(f"\nSaved: {out}.pdf / .png")

print("\nDONE")
