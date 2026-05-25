#!/usr/bin/env python3
"""
疾患別サブグループ解析
- 各疾患ごとに Fisher's exact test (VDR-dominant vs 承認)
- Forest plot + ヒートマップ + 個別散布図
"""

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
from scipy import stats
import warnings
warnings.filterwarnings('ignore')

RESULTS = '/Volumes/M4_SSD/projects/tlr_chipseq/results/'
COL_VDR = '#2166AC'; COL_GR = '#D6604D'
COL_APP = '#1a9641'; COL_FAIL = '#d7191c'
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

# ── データ読み込み ──────────────────────────────────────────────────
df = pd.read_csv(RESULTS + 'gene_disease_phase_expanded.csv')
gene_dis = (df.groupby(['gene','disease','category'])
            .agg(
                best_status = ('mol_status', lambda x:
                               'approved' if 'approved' in x.values else 'failed'),
                VDR_score   = ('VDR_score','first'),
                GR_score    = ('GR_score','first'),
                VDR_dom     = ('VDR_dominant','first'),
            ).reset_index())

# ══════════════════════════════════════════════════════════════════════
# 疾患別 Fisher's exact test
# ══════════════════════════════════════════════════════════════════════
results = []
for dis, grp in gene_dis.groupby('disease'):
    app  = grp[grp['best_status']=='approved']
    fail = grp[grp['best_status']=='failed']
    cat  = grp['category'].iloc[0]
    organ, vdr_tpm = DISEASE_ORGAN.get(dis, ('immune cells', 15.4))

    n_va = int((grp['VDR_dom'] &  (grp['best_status']=='approved')).sum())
    n_vf = int((grp['VDR_dom'] &  (grp['best_status']=='failed')).sum())
    n_ga = int(((~grp['VDR_dom']) & (grp['best_status']=='approved')).sum())
    n_gf = int(((~grp['VDR_dom']) & (grp['best_status']=='failed')).sum())

    # Fisher (n が小さい疾患はOR計算のみ、CI幅広い)
    if n_va + n_vf < 2 or n_ga + n_gf < 2 or n_vf == 0 or n_gf == 0:
        # ゼロセル補正 (+0.5)
        a, b, c, d = n_va+0.5, n_vf+0.5, n_ga+0.5, n_gf+0.5
        or_v = (a*d)/(b*c)
        log_or = np.log(or_v)
        se = np.sqrt(1/a+1/b+1/c+1/d)
        ci_lo = np.exp(log_or - 1.96*se)
        ci_hi = np.exp(log_or + 1.96*se)
        p_v = 1.0
        reliable = False
    else:
        or_v, p_v = stats.fisher_exact([[n_va,n_vf],[n_ga,n_gf]], alternative='greater')
        if or_v == 0 or np.isinf(or_v):
            or_v = (n_va+0.5)*(n_gf+0.5)/((n_vf+0.5)*(n_ga+0.5))
        log_or = np.log(max(or_v, 1e-6))
        se = np.sqrt(1/max(n_va,0.5)+1/max(n_vf,0.5)+
                     1/max(n_ga,0.5)+1/max(n_gf,0.5))
        ci_lo = np.exp(log_or - 1.96*se)
        ci_hi = np.exp(log_or + 1.96*se)
        reliable = (n_va+n_vf >= 3) and (n_ga+n_gf >= 3)

    # VDRスコア Mann-Whitney
    if len(app) >= 2 and len(fail) >= 2:
        _, p_mw = stats.mannwhitneyu(app['VDR_score'], fail['VDR_score'],
                                     alternative='greater')
    else:
        p_mw = np.nan

    results.append({
        'disease': dis, 'category': cat, 'organ': organ, 'vdr_tpm': vdr_tpm,
        'n_total': len(grp), 'n_app': len(app), 'n_fail': len(fail),
        'n_va': n_va, 'n_vf': n_vf, 'n_ga': n_ga, 'n_gf': n_gf,
        'or': or_v, 'ci_lo': ci_lo, 'ci_hi': ci_hi, 'p': p_v, 'p_mw': p_mw,
        'vdr_app_rate': n_va/(n_va+n_vf) if (n_va+n_vf)>0 else np.nan,
        'gr_app_rate':  n_ga/(n_ga+n_gf) if (n_ga+n_gf)>0 else np.nan,
        'vdr_mean_app': app[app['VDR_dom']==True]['VDR_score'].mean() if len(app) else np.nan,
        'vdr_mean_fail':fail[fail['VDR_dom']==True]['VDR_score'].mean() if len(fail) else np.nan,
        'reliable': reliable,
    })

    print(f"{dis[:35]:35s} n={len(grp):3d} "
          f"VDR-dom({n_va}app/{n_vf}fail) GR-dom({n_ga}app/{n_gf}fail) "
          f"OR={or_v:.2f} p={p_v:.3f}{'*' if p_v<0.05 else ''}"
          f"{'  [small n]' if not reliable else ''}")

res = pd.DataFrame(results)
res.to_csv(RESULTS + 'subgroup_by_disease.csv', index=False)

# ══════════════════════════════════════════════════════════════════════
# FIGURE 1: Forest plot (メイン)
# ══════════════════════════════════════════════════════════════════════
fig1, ax = plt.subplots(figsize=(14, 9))

# カテゴリ順に並び替え
cat_order = ['respiratory','skin','GI','renal','autoimmune','eye']
res['cat_ord'] = res['category'].map({c:i for i,c in enumerate(cat_order)}).fillna(99)
res_sorted = res.sort_values(['cat_ord','or'], ascending=[True, False]).reset_index(drop=True)

y_pos = list(range(len(res_sorted)))
for i, (_, row) in enumerate(res_sorted.iterrows()):
    col   = CATEGORY_COLOR.get(row['category'], '#aaa')
    alpha = 1.0 if row['reliable'] else 0.45
    lw    = 2.5 if row['reliable'] else 1.5
    ls    = '-'  if row['reliable'] else '--'

    # CI line
    ci_hi_plot = min(row['ci_hi'], 20)  # 表示上限
    ax.plot([row['ci_lo'], ci_hi_plot], [i, i],
            color=col, lw=lw, ls=ls, alpha=alpha, solid_capstyle='round')
    # OR点
    marker = 'D' if row['p'] < 0.05 else 'o'
    ax.plot(row['or'], i, marker, color=col, ms=10 if row['reliable'] else 7,
            alpha=alpha, zorder=5)
    # p値ラベル
    if row['reliable']:
        p_label = f"p={row['p']:.3f}" if row['p'] >= 0.001 else "p<0.001"
        if row['p'] < 0.05: p_label += ' *'
        ax.text(ci_hi_plot + 0.3, i, p_label, va='center', fontsize=8,
                color='black' if row['p']<0.05 else '#666')
    # サンプルサイズ
    ax.text(-0.5, i,
            f"n={row['n_total']} ({row['n_app']}app/{row['n_fail']}fail)",
            va='center', ha='right', fontsize=8, color='#444')

ax.axvline(1, color='black', lw=1.2, ls='--', alpha=0.7)
ax.set_yticks(y_pos)
ax.set_yticklabels(
    [f"{r['disease'].split(',')[0]}  [{r['organ']}]"
     for _, r in res_sorted.iterrows()], fontsize=10)
ax.set_xlabel('Odds ratio (VDR-dominant: 承認 vs 失敗)', fontsize=11)
ax.set_xlim(-0.8, 22)
ax.set_title(
    '疾患別サブグループ解析: VDR-dominant ターゲットの承認優位性\n'
    '◆=p<0.05, 実線=信頼できるn数, 破線=小サンプル(0.5補正)',
    fontsize=12, fontweight='bold')

# カテゴリ凡例
handles = [mpatches.Patch(color=c, label=k) for k,c in CATEGORY_COLOR.items()
           if k in res['category'].values]
ax.legend(handles=handles, fontsize=9, loc='lower right')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# カテゴリ区切り線
prev_cat = None
for i, (_, row) in enumerate(res_sorted.iterrows()):
    if row['category'] != prev_cat and i > 0:
        ax.axhline(i - 0.5, color='#cccccc', lw=1, ls='-')
    prev_cat = row['category']

fig1.tight_layout()
out1 = RESULTS + 'fig_subgroup_forest'
fig1.savefig(out1+'.pdf', dpi=300, bbox_inches='tight', facecolor='white')
fig1.savefig(out1+'.png', dpi=200, bbox_inches='tight', facecolor='white')
print(f"\nSaved: {out1}.pdf/.png")

# ══════════════════════════════════════════════════════════════════════
# FIGURE 2: 疾患別詳細 (2x7 サブプロット — 承認率比較 + VDRスコア分布)
# ══════════════════════════════════════════════════════════════════════
reliable_dis = res_sorted[res_sorted['reliable']]['disease'].tolist()
n_dis = len(reliable_dis)
ncols = 4
nrows = int(np.ceil(n_dis / ncols))

fig2, axes = plt.subplots(nrows, ncols, figsize=(ncols*4, nrows*3.5))
axes = axes.flatten()

for idx, dis in enumerate(reliable_dis):
    ax = axes[idx]
    grp = gene_dis[gene_dis['disease'] == dis]
    cat = grp['category'].iloc[0]
    col = CATEGORY_COLOR.get(cat, '#aaa')
    organ, _ = DISEASE_ORGAN.get(dis, ('?', 0))

    # VDR-dominant vs GR-dominant 承認率
    vdr_app = grp[grp['VDR_dom']==True]['best_status'].eq('approved').mean()
    gr_app  = grp[grp['VDR_dom']==False]['best_status'].eq('approved').mean()
    n_vdr   = grp['VDR_dom'].sum()
    n_gr    = (~grp['VDR_dom']).sum()

    x     = [0, 1]
    rates = [gr_app, vdr_app]
    ns    = [n_gr, n_vdr]
    cols  = [COL_GR, COL_VDR]
    bars  = ax.bar(x, rates, color=cols, alpha=0.8, width=0.55, edgecolor='white', lw=0.5)
    for bar, rate, n in zip(bars, rates, ns):
        ax.text(bar.get_x()+bar.get_width()/2, rate+0.02,
                f"{rate*100:.0f}%\n(n={n})",
                ha='center', va='bottom', fontsize=8.5, fontweight='bold')

    row_r = res[res['disease']==dis].iloc[0]
    p_str = f"OR={row_r['or']:.2f}, p={row_r['p']:.3f}"
    if row_r['p'] < 0.05: p_str += ' *'

    ax.set_xticks([0,1])
    ax.set_xticklabels(['GR-dominant', 'VDR-dominant'], fontsize=9)
    ax.set_ylabel('承認率', fontsize=9)
    ax.set_ylim(0, min(1.0, max(vdr_app, gr_app)*1.5 + 0.15))
    short = dis.split(',')[0].replace('Pulmonary Disease','COPD')
    ax.set_title(f"{short}\n[{organ}]  {p_str}", fontsize=9, fontweight='bold')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

# 残りのaxを非表示
for idx in range(n_dis, len(axes)):
    axes[idx].set_visible(False)

fig2.suptitle('疾患別: VDR-dominant vs GR-dominant ターゲットの承認率',
              fontsize=13, fontweight='bold', y=1.01)
fig2.tight_layout()
out2 = RESULTS + 'fig_subgroup_approval_rate'
fig2.savefig(out2+'.pdf', dpi=300, bbox_inches='tight', facecolor='white')
fig2.savefig(out2+'.png', dpi=200, bbox_inches='tight', facecolor='white')
print(f"Saved: {out2}.pdf/.png")

# ══════════════════════════════════════════════════════════════════════
# FIGURE 3: 疾患×遺伝子 詳細ヒートマップ
# ══════════════════════════════════════════════════════════════════════
# VDR-dominant 遺伝子 × 疾患マトリクス（承認=4, 失敗=赤, 未試験=白）
vdr_genes = (gene_dis[gene_dis['VDR_dom']==True]
             .groupby('gene')['VDR_score'].first()
             .sort_values(ascending=False).head(30).index.tolist())
dis_list = res_sorted['disease'].tolist()

heat_val  = pd.DataFrame(np.nan, index=vdr_genes, columns=dis_list)
heat_stat = pd.DataFrame('',     index=vdr_genes, columns=dis_list)

for _, row in gene_dis[gene_dis['gene'].isin(vdr_genes)].iterrows():
    if row['disease'] not in dis_list: continue
    heat_val.loc[row['gene'], row['disease']] = (
        4 if row['best_status']=='approved' else 1)
    heat_stat.loc[row['gene'], row['disease']] = row['best_status']

# カスタムカラーマップ: nan=白, 1=赤(失敗), 4=緑(承認)
import matplotlib.colors as mcolors
cmap = mcolors.ListedColormap(['#f5f5f5', COL_FAIL, '#aaaaaa', '#aaaaaa', COL_APP])
bounds = [0, 0.5, 1.5, 2.5, 3.5, 4.5]
norm = mcolors.BoundaryNorm(bounds, cmap.N)

fig3, ax3 = plt.subplots(figsize=(20, 10))
mat = heat_val.fillna(0).values
im  = ax3.imshow(mat, aspect='auto', cmap=cmap, norm=norm, interpolation='nearest')

ax3.set_xticks(range(len(dis_list)))
ax3.set_xticklabels(
    [d.split(',')[0].replace('Pulmonary Disease','COPD')[:18] for d in dis_list],
    rotation=40, ha='right', fontsize=9)
ax3.set_yticks(range(len(vdr_genes)))
ax3.set_yticklabels(
    [f"{g}  (VDR={int(gene_dis[gene_dis['gene']==g]['VDR_score'].iloc[0])})"
     for g in vdr_genes], fontsize=9)

# 凡例
handles3 = [mpatches.Patch(facecolor=COL_APP, label='承認'),
            mpatches.Patch(facecolor=COL_FAIL, label='失敗'),
            mpatches.Patch(facecolor='#f5f5f5', edgecolor='gray', label='未試験')]
ax3.legend(handles=handles3, fontsize=10, loc='lower right',
           bbox_to_anchor=(1.01, 0))

ax3.set_title(
    'VDR-dominant 遺伝子 × ステロイド疾患マップ (上位30遺伝子, VDRスコア順)\n'
    '緑=承認, 赤=失敗(discontinued/withdrawn), 白=未試験',
    fontsize=12, fontweight='bold')
ax3.spines['top'].set_visible(False); ax3.spines['right'].set_visible(False)

fig3.tight_layout()
out3 = RESULTS + 'fig_subgroup_heatmap'
fig3.savefig(out3+'.pdf', dpi=300, bbox_inches='tight', facecolor='white')
fig3.savefig(out3+'.png', dpi=200, bbox_inches='tight', facecolor='white')
print(f"Saved: {out3}.pdf/.png")

# ── サマリー出力 ───────────────────────────────────────────────────
print("\n=== 疾患別サブグループ サマリー ===")
print(res_sorted[['disease','organ','n_total','n_app','n_fail',
                  'or','ci_lo','ci_hi','p','reliable']]
      .to_string(index=False))

sig = res_sorted[res_sorted['p']<0.05]
print(f"\n有意（p<0.05）: {len(sig)} / {len(res_sorted[res_sorted['reliable']])} 疾患")
for _, r in sig.iterrows():
    print(f"  {r['disease']}: OR={r['or']:.2f}, p={r['p']:.3f}")
print("\nDONE")
