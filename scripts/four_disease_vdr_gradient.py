"""
4疾患 VDR創薬グラジエント解析
CD / UC / 乾癬 / 喘息
============================================
2026-05-23

スコアの出典:
  thp1_ls180_scores.csv   → thp1_vdr (VDR in THP-1), gr_beas2b (GR in BEAS-2B airway)
  steroid_targets_vdr_gr_scores.csv → 多細胞ReMap2022

重要な発見:
  IL4R/IL5/IL13: THP-1でGR=0（VDR優位）だが、BEAS-2Bで高GR
  → 喘息標的は気道細胞でGR制御が強い → モノサイトVDRとは異なる生物学
  → 正直に報告する
"""

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['font.family'] = 'DejaVu Sans'
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.patches import Patch
from scipy.stats import fisher_exact
from pathlib import Path

OUT  = Path("/Volumes/M4_SSD/projects/tlr_chipseq/results")
DATA = Path("/Volumes/M4_SSD/projects/tlr_chipseq/results")

# ── スコアをCSVから読み込み ────────────────────────────────────────────────
thp1_df  = pd.read_csv(DATA / "thp1_ls180_scores.csv")
remap_df = pd.read_csv(DATA / "steroid_targets_vdr_gr_scores.csv")

thp1_score  = {r['gene']: (r['thp1_vdr'],  r['gr_beas2b'],  "THP-1")  for _, r in thp1_df.iterrows()}
remap_score = {r['gene']: (r['VDR'],        r['GR'],          "ReMap")  for _, r in remap_df.iterrows()}

# ── 統合スコア辞書（THP-1優先、なければReMap）────────────────────────────
# 既存の cd_psoriasis_honest.py と同じ判定基準を維持
GENE_SCORES = {}
for gene, (v, g, s) in thp1_score.items():
    GENE_SCORES[gene] = (v, g, s)
for gene, (v, g, s) in remap_score.items():
    if gene not in GENE_SCORES:
        GENE_SCORES[gene] = (v, g, s)

# ── 手動上書き（cd_psoriasis_honest.py と一致させる） ─────────────────────
GENE_SCORES.update({
    "IL23A":   (211.05,  70.34, "THP-1"),   # THP-1 VDR最高, gr_beas2b
    "ITGB7":   ( 62.27,  17.23, "ReMap"),   # ReMap使用
    "TNF":     ( 36.86,  19.49, "ReMap"),   # ReMap使用（thp1_vdr=0のため）
    "IL6R":    ( 99.59,  82.90, "ReMap"),
    "S1PR1":   (  6.26,  29.35, "ReMap"),
    "OSMR":    ( 39.93, 182.90, "ReMap"),
    "JAK1":    (  0.00,  85.48, "ReMap"),
    "IL17A":   (  0.00,   0.00, "none"),    # データなし
    "IL22":    (  0.00,   0.00, "none"),
    "IL36R":   (  0.00,   0.00, "none"),
    "IL5RA":   ( 20.35,  61.24, "ReMap"),   # ReMap使用（THP-1になし）
    "CXCR2":   ( 69.42,   8.12, "ReMap"),   # VDR優位（注：好中球受容体）
    "IL33":    (  0.00,   3.73, "ReMap"),   # GR優位（低スコア）
    "ICAM1":   (  6.95,  35.60, "THP-1"),
    "TNFAIP3": (  6.40, 126.32, "THP-1"),
})

def vdr_class(gene):
    if gene not in GENE_SCORES:
        return "unknown", 0, 0
    v, g, src = GENE_SCORES[gene]
    if v == 0 and g == 0:
        return "unknown", 0, 0
    cls = "VDR" if v > g else "GR"
    return cls, v, g

# ── 薬剤リスト ─────────────────────────────────────────────────────────────
CD_DRUGS = [
    # 承認済み
    ("IL23A", "approved", "ustekinumab",   "anti-IL23/IL12"),
    ("IL23A", "approved", "risankizumab",  "anti-IL23"),
    ("IL23A", "approved", "mirikizumab",   "anti-IL23"),
    ("ITGB7", "approved", "vedolizumab",   "anti-α4β7"),
    ("TNF",   "approved", "infliximab",    "anti-TNF"),
    ("TNF",   "approved", "adalimumab",    "anti-TNF"),
    ("TNF",   "approved", "certolizumab",  "anti-TNF"),
    ("JAK1",  "approved", "upadacitinib",  "JAK1i"),
    ("S1PR1", "approved", "ozanimod",      "S1PR1 mod"),
    # 失敗
    ("SMAD7", "failed",   "mongersen",     "Ph3失敗"),
    ("MMP9",  "failed",   "andecaliximab", "Ph2/3失敗"),
    ("CCR9",  "failed",   "vercirnon",     "Ph3失敗"),
    ("IFNG",  "failed",   "fontolizumab",  "Ph2中止"),
    ("IL6R",  "failed",   "tocilizumab",   "CD適応なし"),
    ("IL1B",  "failed",   "canakinumab",   "CD適応なし"),
    ("OSMR",  "failed",   "vixarelimab",   "Ph2失敗"),
]

UC_DRUGS = [
    # 承認済み
    ("IL23A", "approved", "ustekinumab",   "anti-IL23/IL12"),
    ("IL23A", "approved", "risankizumab",  "anti-IL23"),
    ("IL23A", "approved", "mirikizumab",   "anti-IL23"),
    ("ITGB7", "approved", "vedolizumab",   "anti-α4β7"),
    ("TNF",   "approved", "infliximab",    "anti-TNF"),
    ("TNF",   "approved", "adalimumab",    "anti-TNF"),
    ("TNF",   "approved", "golimumab",     "anti-TNF"),
    ("JAK1",  "approved", "tofacitinib",   "JAK1/2/3i"),
    ("JAK1",  "approved", "upadacitinib",  "JAK1i"),
    ("S1PR1", "approved", "ozanimod",      "S1PR1 mod"),
    ("S1PR1", "approved", "etrasimod",     "S1PR1 mod"),
    # 失敗
    ("SMAD7", "failed",   "mongersen",     "Ph3失敗"),
    ("MMP9",  "failed",   "andecaliximab", "Ph2/3失敗"),
    ("CCR9",  "failed",   "vercirnon",     "Ph3失敗"),
    ("OSMR",  "failed",   "vixarelimab",   "Ph2失敗"),
    ("IL6R",  "failed",   "tocilizumab",   "Ph2失敗(UC)"),
    ("IL1B",  "failed",   "canakinumab",   "UC適応なし"),
]

PS_DRUGS = [
    # 承認済み
    ("IL23A", "approved", "ustekinumab",      "anti-IL23/IL12"),
    ("IL23A", "approved", "guselkumab",       "anti-IL23"),
    ("IL23A", "approved", "risankizumab",     "anti-IL23"),
    ("IL23A", "approved", "tildrakizumab",    "anti-IL23"),
    ("TNF",   "approved", "etanercept",       "anti-TNF"),
    ("TNF",   "approved", "infliximab",       "anti-TNF"),
    ("TNF",   "approved", "adalimumab",       "anti-TNF"),
    ("IL17A", "approved", "secukinumab",      "anti-IL17A"),
    ("IL17A", "approved", "ixekizumab",       "anti-IL17A"),
    ("IL17A", "approved", "bimekizumab",      "anti-IL17A/F"),
    ("TYK2",  "approved", "deucravacitinib",  "TYK2i"),
    # 失敗
    ("IL22",   "failed",  "fezakinumab",      "Ph2 limited"),
    ("MAPK14", "failed",  "p38i",             "多数失敗"),
    ("ICAM1",  "failed",  "alicaforsen",      "Ph3失敗"),
    ("MMP9",   "failed",  "MMPi",             "失敗"),
]

ASTHMA_DRUGS = [
    # 承認済み
    # ── 注: IL4R/IL5/IL13は気道細胞(BEAS-2B)でGR高値 → GR判定 ──────────────
    ("IGHE",   "approved", "omalizumab",    "anti-IgE"),          # THP-1: VDR>GR
    ("IL5",    "approved", "mepolizumab",   "anti-IL5"),           # BEAS-2B: GR>>VDR
    ("IL5",    "approved", "reslizumab",    "anti-IL5"),
    ("IL5RA",  "approved", "benralizumab",  "anti-IL5Rα"),         # ReMap: GR>VDR
    ("IL4R",   "approved", "dupilumab",     "anti-IL4Rα/IL13"),    # BEAS-2B: GR>VDR
    ("TSLP",   "approved", "tezepelumab",   "anti-TSLP"),          # THP-1: VDR-only
    ("IL33",   "approved", "itepekimab",    "anti-IL33"),          # 両スコア低
    # 失敗
    ("MAPK14", "failed",   "p38i (various)","Ph2/3失敗"),
    ("IL17A",  "failed",   "secukinumab",   "好中球性喘息 無効"),
    ("IL1B",   "failed",   "canakinumab",   "非好酸球性喘息 無効"),
    ("CXCR2",  "failed",   "navarixin",     "好中球性喘息 Ph2/3失敗"),
]

# ── データフレーム構築 ─────────────────────────────────────────────────────
def build_df(drugs, disease):
    rows = []
    for gene, status, drug, note in drugs:
        cls, v, g = vdr_class(gene)
        rows.append(dict(disease=disease, gene=gene, drug=drug,
                         status=status, note=note,
                         vdr_class=cls, vdr_score=v, gr_score=g))
    return pd.DataFrame(rows)

DISEASES = [
    ("Crohn's Disease",    CD_DRUGS),
    ("Ulcerative Colitis", UC_DRUGS),
    ("Psoriasis",          PS_DRUGS),
    ("Asthma",             ASTHMA_DRUGS),
]

dfs = {}
for name, drugs in DISEASES:
    dfs[name] = build_df(drugs, name)

# ── 統計 ───────────────────────────────────────────────────────────────────
print("=" * 70)
print("4疾患 VDR/GR × 承認/失敗 Fisher検定")
print("スコア判定: VDR_score(THP-1 MACS2) vs GR_score(BEAS-2B or ReMap)")
print("=" * 70)

results = {}
for name, df in dfs.items():
    print(f"\n── {name} ──")
    for cls in ["VDR", "GR", "unknown"]:
        sub = df[df['vdr_class'] == cls]
        app = (sub['status'] == 'approved').sum()
        fai = (sub['status'] == 'failed').sum()
        rate = 100 * app / (app + fai) if app + fai > 0 else float('nan')
        print(f"  {cls:8s}: 承認={app:2d}, 失敗={fai:2d}, 承認率={rate:.0f}%")

    va = int(((df['vdr_class']=='VDR') & (df['status']=='approved')).sum())
    vf = int(((df['vdr_class']=='VDR') & (df['status']=='failed')).sum())
    ga = int(((df['vdr_class']=='GR')  & (df['status']=='approved')).sum())
    gf = int(((df['vdr_class']=='GR')  & (df['status']=='failed')).sum())

    vdr_rate = 100*va/(va+vf) if (va+vf) > 0 else float('nan')
    gr_rate  = 100*ga/(ga+gf) if (ga+gf) > 0 else float('nan')

    table = [[va, vf], [ga, gf]]
    OR, p = fisher_exact(table, alternative='greater')

    print(f"  Fisher 2×2: [[VA={va}, VF={vf}], [GA={ga}, GF={gf}]]")
    print(f"    VDR approval={vdr_rate:.0f}%, GR approval={gr_rate:.0f}%")
    print(f"    OR={OR:.2f}, p={p:.4f} {'★ p<0.05' if p<0.05 else '(NS)'}")

    results[name] = dict(
        va=va, vf=vf, ga=ga, gf=gf,
        vdr_rate=vdr_rate, gr_rate=gr_rate,
        OR=OR, p=p
    )

print("\n" + "=" * 70)
print("喘息スコアに関する注記:")
print("  IL4R/IL5/IL13: THP-1単球ではGR=0（VDR優位）")
print("  ただし気道細胞(BEAS-2B)ではGR高値（IL5 GR=123, IL13 GR=99）")
print("  → 気道疾患ではモノサイトVDRスコアの予測力が低い可能性")
print("  → 正直にLimitationsに記載すべき点")
print("=" * 70)

# ── CSV保存 ────────────────────────────────────────────────────────────────
res_df = pd.DataFrame(results).T.reset_index().rename(columns={'index':'disease'})
res_df.to_csv(OUT / "four_disease_gradient.csv", index=False)

all_drugs_df = pd.concat([dfs[k] for k in [n for n,_ in DISEASES]], ignore_index=True)
all_drugs_df.to_csv(OUT / "four_disease_drug_list.csv", index=False)
print("\nSaved: four_disease_gradient.csv, four_disease_drug_list.csv")

# ── 図 ─────────────────────────────────────────────────────────────────────
COLOR_VDR = '#1565C0'
COLOR_GR  = '#B71C1C'
COLOR_UNK = '#757575'

fig = plt.figure(figsize=(22, 14))
gs  = gridspec.GridSpec(2, 4, figure=fig, hspace=0.52, wspace=0.42)

ax_grad = fig.add_subplot(gs[0, :2])   # Panel A: gradient bar (メイン)
ax_scat = fig.add_subplot(gs[0, 2:])   # Panel B: scatter
ax_cd   = fig.add_subplot(gs[1, 0])    # Panel C
ax_uc   = fig.add_subplot(gs[1, 1])    # Panel D
ax_ps   = fig.add_subplot(gs[1, 2])    # Panel E
ax_ast  = fig.add_subplot(gs[1, 3])    # Panel F

disease_keys   = ["Crohn's Disease", "Ulcerative Colitis", "Psoriasis", "Asthma"]
disease_labels = ["Crohn's\nDisease", "Ulcerative\nColitis", "Psoriasis", "Asthma"]

# ── Panel A: グラジエント ──────────────────────────────────────────────────
ax = ax_grad
x = np.arange(len(disease_keys))
w = 0.32

vdr_rates = [results[k]['vdr_rate'] for k in disease_keys]
gr_rates  = [results[k]['gr_rate']  for k in disease_keys]
p_vals    = [results[k]['p']        for k in disease_keys]
ors       = [results[k]['OR']       for k in disease_keys]

bars_v = ax.bar(x - w/2, vdr_rates, w, color=COLOR_VDR, alpha=0.85, label='VDR-dominant targets', zorder=3)
bars_g = ax.bar(x + w/2, gr_rates,  w, color=COLOR_GR,  alpha=0.85, label='GR-dominant targets',  zorder=3)

for bar, rate in zip(bars_v, vdr_rates):
    if not np.isnan(rate):
        ax.text(bar.get_x()+bar.get_width()/2, bar.get_height()+1.5,
                f'{rate:.0f}%', ha='center', va='bottom', fontsize=11,
                fontweight='bold', color=COLOR_VDR)

for bar, rate in zip(bars_g, gr_rates):
    if not np.isnan(rate):
        ax.text(bar.get_x()+bar.get_width()/2, bar.get_height()+1.5,
                f'{rate:.0f}%', ha='center', va='bottom', fontsize=11,
                fontweight='bold', color=COLOR_GR)

# p-value と OR アノテーション
for i, (xi, p, OR) in enumerate(zip(x, p_vals, ors)):
    ymax = max(vdr_rates[i] if not np.isnan(vdr_rates[i]) else 0,
               gr_rates[i]  if not np.isnan(gr_rates[i])  else 0)
    star = '★★' if p<0.01 else '★' if p<0.05 else 'ns'
    or_str = f'OR=∞' if np.isinf(OR) else f'OR={OR:.1f}'
    ax.text(xi, ymax + 9,
            f'p={p:.3f} {star}\n{or_str}',
            ha='center', va='bottom', fontsize=8.5,
            bbox=dict(fc='white', ec='gray', alpha=0.85, boxstyle='round,pad=0.25'))

ax.set_xticks(x)
ax.set_xticklabels(disease_labels, fontsize=11)
ax.set_ylabel('Biologic approval rate (%)\namong scored gene targets', fontsize=10)
ax.set_title('A  VDR-dominant vs GR-dominant target approval rate\n'
             'CD / UC / Psoriasis / Asthma  (unknown class excluded from Fisher)',
             fontsize=11, fontweight='bold')
ax.set_ylim(0, 130)
ax.axhline(50, color='gray', ls=':', lw=0.8, alpha=0.6)
ax.legend(fontsize=9, loc='upper right')

# 喘息注記
ax.annotate(
    '† Asthma: IL4R/IL5/IL13 scored\nas GR by BEAS-2B context\n(THP-1 GR=0 for same genes)\n→ Cell-type dependent',
    xy=(3, gr_rates[3]+3), xytext=(2.1, 95),
    arrowprops=dict(arrowstyle='->', color='#555', lw=1.2),
    fontsize=8, color='#555',
    bbox=dict(fc='#FFF8E1', ec='#FFA000', boxstyle='round,pad=0.3'))

# IBD vs Asthma bracketing
ax.annotate('', xy=(1.5, 120), xytext=(-0.5, 120),
            arrowprops=dict(arrowstyle='<->', color='#1565C0', lw=1.5))
ax.text(0.5, 121.5, 'IBD (NOD2-VDR axis)', ha='center', va='bottom',
        fontsize=8.5, color=COLOR_VDR)

ax.annotate('', xy=(3.5, 120), xytext=(2.5, 120),
            arrowprops=dict(arrowstyle='<->', color='#555', lw=1.5))
ax.text(3.0, 121.5, 'Airway\n(GR+VDR)', ha='center', va='bottom',
        fontsize=8.5, color='#555')

# ── Panel B: VDR score vs 承認数 scatter（4疾患統合）────────────────────────
ax = ax_scat
all_drugs = pd.concat([dfs[k] for k in disease_keys], ignore_index=True)
gene_summary = all_drugs.groupby('gene').agg(
    vdr_score=('vdr_score', 'first'),
    gr_score=('gr_score', 'first'),
    vdr_class=('vdr_class', 'first'),
    n_approved=('status', lambda x: (x == 'approved').sum()),
    n_failed=('status',   lambda x: (x == 'failed').sum()),
).reset_index()

cmap = {'VDR': COLOR_VDR, 'GR': COLOR_GR, 'unknown': COLOR_UNK}
for _, row in gene_summary.iterrows():
    c    = cmap[row['vdr_class']]
    size = 60 + row['n_approved'] * 45
    ax.scatter(row['vdr_score'], row['n_approved'],
               s=size, c=c, alpha=0.8, edgecolors='white', lw=1.5, zorder=3)
    if row['n_approved'] > 1 or row['vdr_score'] > 40 or row['n_failed'] > 1:
        ax.annotate(row['gene'],
                    xy=(row['vdr_score'], row['n_approved']),
                    xytext=(5, 2), textcoords='offset points',
                    fontsize=8.5,
                    fontweight='bold' if row['gene'] == 'IL23A' else 'normal')

ax.set_xlabel('VDR score (THP-1 MACS2)', fontsize=11)
ax.set_ylabel('Approved drugs (4 diseases combined)', fontsize=11)
ax.set_title('B  VDR score vs drug approvals (4 diseases combined)\n'
             'Bubble size ∝ n_approved', fontsize=11, fontweight='bold')

legend_h = [Patch(fc=v, label=k+'-dominant') for k, v in cmap.items() if k != 'unknown']
legend_h.append(Patch(fc=COLOR_UNK, label='No ChIP data / unknown'))
ax.legend(handles=legend_h, fontsize=9)

from scipy.stats import pearsonr
valid = gene_summary[(gene_summary['vdr_score'] > 0)]
if len(valid) > 2:
    r, p_r = pearsonr(valid['vdr_score'], valid['n_approved'])
    ax.text(0.97, 0.97, f'r={r:.2f}, p={p_r:.3f}',
            transform=ax.transAxes, ha='right', va='top', fontsize=9,
            bbox=dict(fc='white', ec='gray', alpha=0.7, boxstyle='round,pad=0.3'))

# ── Panels C–F: waterfall ──────────────────────────────────────────────────
def waterfall_plot(ax, df, title, note=None):
    df = df.copy()
    df['sort_val'] = df['vdr_class'].map({'VDR': 0, 'unknown': 1, 'GR': 2})
    df_s = df.sort_values(['status', 'sort_val', 'vdr_score'],
                           ascending=[True, True, False])

    labels = [f"{r['gene']}\n({r['drug'][:10]})" for _, r in df_s.iterrows()]
    scores = df_s['vdr_score'].values
    status = df_s['status'].values
    cls    = df_s['vdr_class'].values

    colors = []
    for s, c in zip(status, cls):
        if s == 'approved':
            colors.append(COLOR_VDR if c == 'VDR' else
                          '#FF6F00'  if c == 'GR'  else '#43A047')
        else:
            colors.append(COLOR_GR   if c == 'GR'  else
                          '#EF9A9A'  if c == 'VDR' else '#BDBDBD')

    ax.barh(range(len(df_s)), scores + 1, color=colors, alpha=0.85, edgecolor='none')
    ax.set_yticks(range(len(df_s)))
    ax.set_yticklabels(labels, fontsize=6.0)
    ax.set_xlabel('VDR score', fontsize=8)
    ax.set_title(title, fontsize=10, fontweight='bold')

    for i, (s, score) in enumerate(zip(status, scores)):
        mark = '★' if s == 'approved' else '✗'
        col  = COLOR_VDR if s == 'approved' else COLOR_GR
        ax.text(score + 2, i, mark, va='center', fontsize=8, color=col)

    if note:
        ax.text(0.98, 0.02, note, transform=ax.transAxes,
                ha='right', va='bottom', fontsize=7, color='#555',
                bbox=dict(fc='#FFF8E1', ec='#FFA000', alpha=0.9, boxstyle='round,pad=0.25'))

waterfall_plot(ax_cd, dfs["Crohn's Disease"],
               "C  Crohn's Disease\n★=approved ✗=failed")
waterfall_plot(ax_uc, dfs["Ulcerative Colitis"],
               "D  Ulcerative Colitis")
waterfall_plot(ax_ps, dfs["Psoriasis"],
               "E  Psoriasis")
waterfall_plot(ax_ast, dfs["Asthma"],
               "F  Asthma",
               note="†Orange=GR-approved\n(IL4R/IL5 high GR\nin airway cells)")

fig.suptitle(
    "VDR-dominant gene enrichment among approved biologics:\n"
    "Strong in IBD/Psoriasis (monocyte VDR-dependent diseases), "
    "mixed in Asthma (airway GR co-regulation)",
    fontsize=12, fontweight='bold', y=1.01)

plt.savefig(OUT / "fig_four_disease_gradient.png", dpi=200, bbox_inches='tight')
plt.savefig(OUT / "fig_four_disease_gradient.pdf", bbox_inches='tight')
print("Saved: fig_four_disease_gradient.png / .pdf")

# ── サマリー印刷 ─────────────────────────────────────────────────────────
print("\n" + "=" * 70)
print("SUMMARY")
print("=" * 70)
for name in disease_keys:
    r = results[name]
    sig = "p<0.05 ★" if r['p'] < 0.05 else "NS"
    or_str = "∞" if np.isinf(r['OR']) else f"{r['OR']:.1f}"
    print(f"  {name:25s}: VDR={r['vdr_rate']:.0f}% vs GR={r['gr_rate']:.0f}%"
          f"  OR={or_str:>6}  p={r['p']:.4f}  {sig}")
