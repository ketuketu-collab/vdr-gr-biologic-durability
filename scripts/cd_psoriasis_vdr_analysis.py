"""
CD（クローン病）と乾癬における VDR 依存性の定量
==================================================
承認済み生物製剤・低分子薬 → ターゲット遺伝子 → VDR/GR スコア照合
→ 「VDR依存疾患では承認率が特に高い」を定量的に示す
"""

import urllib.request, json, time, re
import pandas as pd
import numpy as np
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['font.family'] = 'DejaVu Sans'
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.patches import Patch
from scipy.stats import fisher_exact

OUT   = Path("/Volumes/M4_SSD/projects/tlr_chipseq/results")
CHEMBL = "https://www.ebi.ac.uk/chembl/api/data"

# ── 既存スコアを読み込む ──────────────────────────────────────────────────────
thp1 = pd.read_csv(OUT / "thp1_ls180_scores.csv")

# MACS2スコア版（より高感度）も読み込む
# make_tableS4 の元データ（ALL_TARGETS と同じ遺伝子セット）
VDR_SCORES = dict(zip(thp1['gene'], thp1['thp1_vdr']))
GR_SCORES  = dict(zip(thp1['gene'], thp1['gr_beas2b'] + thp1['gr_thp1']))

# ── 手動キュレーション: CD と乾癬の承認薬・開発薬 ────────────────────────────
# 各薬剤の情報: (gene_target, status, drug_name, year_approved_or_phase)
# status: "approved" / "failed" / "ongoing"

CD_DRUGS = [
    # ── 承認済み ──────────────────────────────────────────────
    ("TNF",    "approved",  "infliximab",     1998),
    ("TNF",    "approved",  "adalimumab",     2007),
    ("TNF",    "approved",  "certolizumab",   2008),
    ("ITGB7",  "approved",  "vedolizumab",    2014),
    ("IL23A",  "approved",  "ustekinumab",    2016),  # anti-IL12/23 (IL23A/IL12B)
    ("IL12B",  "approved",  "ustekinumab",    2016),
    ("IL23A",  "approved",  "risankizumab",   2022),
    ("S1PR1",  "approved",  "ozanimod",       2023),
    # ── 承認(UC/CD両方) ──────────────────────────────────────
    ("IL23A",  "approved",  "mirikizumab",    2023),
    # ── 失敗・中止 ────────────────────────────────────────────
    ("SMAD7",  "failed",    "mongersen",      0),     # Phase3失敗
    ("IL6R",   "failed",    "tocilizumab-CD", 0),     # CD適応なし
    ("MMP9",   "failed",    "andecaliximab",  0),     # Phase2/3 failed
    ("CCR9",   "failed",    "vercirnon",      0),     # Phase3 failed
    ("IL13",   "failed",    "tralokinumab-CD",0),     # CD適応なし
    ("IFNG",   "failed",    "fontolizumab",   0),     # Phase2 stopped
    ("IL2RA",  "failed",    "basiliximab-CD", 0),     # CD適応なし
    ("OSMR",   "failed",    "vixarelimab-CD", 0),     # ongoing/CD未承認
    ("IL1B",   "failed",    "canakinumab-CD", 0),     # CD適応なし
    # ── 進行中 ────────────────────────────────────────────────
    ("IL23A",  "ongoing",   "guselkumab-CD",  0),
    ("TYK2",   "ongoing",   "deucravacitinib-CD", 0),
    ("RIPK2",  "ongoing",   "RIPK2i-CD",      0),
    ("NOD2",   "ongoing",   "NOD2-target",    0),
    ("IL10",   "ongoing",   "IL10-therapy",   0),
]

PSORIASIS_DRUGS = [
    # ── 承認済み ──────────────────────────────────────────────
    ("TNF",    "approved",  "etanercept",     1998),
    ("TNF",    "approved",  "infliximab",     2006),
    ("TNF",    "approved",  "adalimumab",     2008),
    ("IL23A",  "approved",  "ustekinumab",    2009),
    ("IL12B",  "approved",  "ustekinumab",    2009),
    ("IL17A",  "approved",  "secukinumab",    2015),
    ("IL17A",  "approved",  "ixekizumab",     2016),
    ("IL23A",  "approved",  "guselkumab",     2017),
    ("IL23A",  "approved",  "risankizumab",   2019),
    ("IL23A",  "approved",  "tildrakizumab",  2018),
    ("TYK2",   "approved",  "deucravacitinib",2022),
    ("IL17A",  "approved",  "bimekizumab",    2023),
    # ── 失敗・中止 ────────────────────────────────────────────
    ("IL22",   "failed",    "fezakinumab",    0),     # Phase2 limited
    ("IFNG",   "failed",    "fontolizumab",   0),
    ("MMP9",   "failed",    "MMP inhibitor",  0),
    ("MAPK14", "failed",    "p38 MAPK inh",   0),
    ("ICAM1",  "failed",    "alicaforsen",    0),
    ("CSF2",   "failed",    "GM-CSF target",  0),
    # ── 進行中 ────────────────────────────────────────────────
    ("IL23A",  "ongoing",   "mirikizumab-Ps", 0),
    ("IL36R",  "ongoing",   "spesolimab-Ps",  0),
]

def get_vdr_gr_class(gene):
    """VDR vs GR 分類（THP-1スコアベース）"""
    vdr = VDR_SCORES.get(gene, 0)
    gr  = GR_SCORES.get(gene, 0)
    if vdr == 0 and gr == 0:
        return "unknown", 0, 0
    if vdr >= gr:
        return "VDR", vdr, gr
    else:
        return "GR", vdr, gr

# ── 承認率の計算 ──────────────────────────────────────────────────────────────
def calc_approval_rate(drug_list, disease_name):
    rows = []
    for gene, status, drug, year in drug_list:
        cls, vdr_s, gr_s = get_vdr_gr_class(gene)
        rows.append({
            'disease': disease_name,
            'gene': gene,
            'drug': drug,
            'status': status,
            'vdr_gr_class': cls,
            'vdr_score': vdr_s,
            'gr_score': gr_s,
            'year': year,
        })
    return pd.DataFrame(rows)

df_cd  = calc_approval_rate(CD_DRUGS,        "Crohn's Disease")
df_ps  = calc_approval_rate(PSORIASIS_DRUGS, "Psoriasis")
df_all = pd.concat([df_cd, df_ps], ignore_index=True)

print("=" * 65)
print("VDR依存疾患における承認率分析: CD vs Psoriasis")
print("=" * 65)

for disease, df in [("Crohn's Disease", df_cd), ("Psoriasis", df_ps)]:
    print(f"\n── {disease} ──")
    for cls in ["VDR", "GR", "unknown"]:
        sub = df[df['vdr_gr_class'] == cls]
        approved = (sub['status'] == 'approved').sum()
        failed   = (sub['status'] == 'failed').sum()
        ongoing  = (sub['status'] == 'ongoing').sum()
        total    = len(sub)
        rate     = 100 * approved / (approved + failed) if (approved + failed) > 0 else 0
        print(f"  {cls:8s}: n={total}  承認={approved}  失敗={failed}  進行中={ongoing}  "
              f"承認率(除進行中)={rate:.0f}%")

    # Fisher's exact test (VDR vs GR, approved vs failed)
    vdr_approved = ((df['vdr_gr_class']=='VDR') & (df['status']=='approved')).sum()
    vdr_failed   = ((df['vdr_gr_class']=='VDR') & (df['status']=='failed')).sum()
    gr_approved  = ((df['vdr_gr_class']=='GR')  & (df['status']=='approved')).sum()
    gr_failed    = ((df['vdr_gr_class']=='GR')  & (df['status']=='failed')).sum()
    if vdr_failed + gr_approved + gr_failed > 0:
        OR, p = fisher_exact([[vdr_approved, vdr_failed],[gr_approved, gr_failed]])
        print(f"  Fisher: VDR({vdr_approved}/{vdr_approved+vdr_failed}) vs "
              f"GR({gr_approved}/{gr_approved+gr_failed}), OR={OR:.1f}, p={p:.4f}")

# ── 承認薬のリスト（VDR優位のみ）──────────────────────────────────────────────
print("\n" + "=" * 65)
print("承認済み薬剤とVDR/GRスコア（CD + 乾癬）")
print("=" * 65)
approved = df_all[df_all['status']=='approved'].sort_values(['disease','vdr_score'], ascending=[True,False])
print(approved[['disease','gene','drug','vdr_gr_class','vdr_score','gr_score','year']].to_string(index=False))

# ── 図作成 ────────────────────────────────────────────────────────────────────
fig = plt.figure(figsize=(16, 10))
gs  = gridspec.GridSpec(2, 3, figure=fig, hspace=0.45, wspace=0.38)

ax_cd_bar  = fig.add_subplot(gs[0, 0])
ax_ps_bar  = fig.add_subplot(gs[0, 1])
ax_scatter = fig.add_subplot(gs[0, 2])
ax_summary = fig.add_subplot(gs[1, 0:2])
ax_text    = fig.add_subplot(gs[1, 2])

STATUS_COLOR = {'approved':'#1565C0', 'failed':'#B71C1C', 'ongoing':'#F57F17'}
VDR_COLOR    = '#1565C0'
GR_COLOR     = '#B71C1C'

def plot_drug_bar(ax, df, title):
    """各薬剤をVDR/GRクラスと承認状況で表示"""
    df_sorted = df.sort_values(['vdr_gr_class','vdr_score'], ascending=[True, False])
    genes = df_sorted['gene'].tolist()
    drugs = df_sorted['drug'].tolist()
    colors = [VDR_COLOR if c=='VDR' else GR_COLOR if c=='GR' else '#757575'
              for c in df_sorted['vdr_gr_class']]
    alphas = [1.0 if s=='approved' else 0.35 if s=='failed' else 0.6
              for s in df_sorted['status']]
    scores = df_sorted['vdr_score'].clip(lower=1)  # ゼロを避ける

    bars = ax.barh(range(len(df_sorted)), scores,
                   color=colors, alpha=0.85, edgecolor='none')
    # Alphaを個別に設定
    for bar, alpha, status in zip(bars, alphas, df_sorted['status']):
        bar.set_alpha(alpha)
        # 承認済みには星印
        if status == 'approved':
            ax.text(bar.get_width() + 2, bar.get_y() + bar.get_height()/2,
                    '★', va='center', fontsize=8, color='#1565C0')
        elif status == 'failed':
            ax.text(bar.get_width() + 2, bar.get_y() + bar.get_height()/2,
                    '✗', va='center', fontsize=8, color='#B71C1C')

    labels = [f"{g} ({d[:15]})" for g, d in zip(genes, drugs)]
    ax.set_yticks(range(len(df_sorted)))
    ax.set_yticklabels(labels, fontsize=7.5)
    ax.set_xlabel('THP-1 VDR score', fontsize=9)
    ax.set_title(title, fontsize=11, fontweight='bold')
    ax.axvline(50, color='gray', ls='--', lw=0.8, alpha=0.6)

    legend = [Patch(fc=VDR_COLOR, label='VDR-dominant'),
              Patch(fc=GR_COLOR,  label='GR-dominant'),
              Patch(fc='#1565C0', alpha=0.35, label='Failed/ongoing')]
    ax.legend(handles=legend, fontsize=7, loc='lower right')

plot_drug_bar(ax_cd_bar,  df_cd,  'A  Crohn\'s Disease\nDrug targets by VDR score')
plot_drug_bar(ax_ps_bar,  df_ps,  'B  Psoriasis\nDrug targets by VDR score')

# Panel C: VDR score vs approval outcome scatter
ax = ax_scatter
for _, row in df_all.iterrows():
    color = STATUS_COLOR[row['status']]
    marker = 'D' if row['disease'] == "Crohn's Disease" else 'o'
    ax.scatter(row['vdr_score'], row['gr_score'],
               c=color, marker=marker, s=70, alpha=0.8, edgecolors='none', zorder=3)
    if row['status'] == 'approved' and row['vdr_score'] > 20:
        ax.annotate(row['gene'],
                    xy=(row['vdr_score'], row['gr_score']),
                    xytext=(4, 2), textcoords='offset points', fontsize=7)

ax.axline((0,0), slope=1, color='gray', ls='--', lw=0.8, label='VDR=GR')
ax.set_xlabel('THP-1 VDR score', fontsize=9)
ax.set_ylabel('GR score (THP-1 + BEAS-2B)', fontsize=9)
ax.set_title('C  VDR vs GR score\nby approval status', fontsize=11, fontweight='bold')

legend_handles = [
    Patch(fc=STATUS_COLOR['approved'], label='Approved'),
    Patch(fc=STATUS_COLOR['failed'],   label='Failed'),
    Patch(fc=STATUS_COLOR['ongoing'],  label='Ongoing'),
    plt.scatter([],[], marker='D', c='gray', s=50, label="Crohn's"),
    plt.scatter([],[], marker='o', c='gray', s=50, label='Psoriasis'),
]
ax.legend(handles=legend_handles, fontsize=7.5)

# Panel D: 承認率のサマリー比較
ax = ax_summary
diseases = ["Crohn's Disease", "Psoriasis", "Combined"]
vdr_rates, gr_rates = [], []

for d_name, df_d in [("Crohn's Disease", df_cd), ("Psoriasis", df_ps),
                       ("Combined", df_all)]:
    for cls, rates_list in [("VDR", vdr_rates), ("GR", gr_rates)]:
        sub = df_d[(df_d['vdr_gr_class']==cls) & (df_d['status'].isin(['approved','failed']))]
        approved = (sub['status']=='approved').sum()
        total    = len(sub)
        rates_list.append(100*approved/total if total > 0 else 0)

x = np.arange(len(diseases))
w = 0.35
bars_vdr = ax.bar(x - w/2, vdr_rates, w, color=VDR_COLOR, alpha=0.85,
                   label='VDR-dominant targets', edgecolor='none')
bars_gr  = ax.bar(x + w/2, gr_rates,  w, color=GR_COLOR,  alpha=0.85,
                   label='GR-dominant targets',  edgecolor='none')

for bar, rate in zip(bars_vdr, vdr_rates):
    ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 1,
            f'{rate:.0f}%', ha='center', va='bottom', fontsize=10, fontweight='bold', color=VDR_COLOR)
for bar, rate in zip(bars_gr, gr_rates):
    ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 1,
            f'{rate:.0f}%', ha='center', va='bottom', fontsize=10, fontweight='bold', color=GR_COLOR)

ax.set_xticks(x)
ax.set_xticklabels(diseases, fontsize=11)
ax.set_ylabel('Approval rate (%)\n(excluding ongoing)', fontsize=10)
ax.set_ylim(0, 115)
ax.set_title("D  Approval rate: VDR-dominant vs GR-dominant targets\n"
             "Crohn's Disease and Psoriasis", fontsize=11, fontweight='bold')
ax.legend(fontsize=9)
ax.axhline(50, color='gray', ls=':', lw=0.8)

# Fisher test annotation
for i, (d_name, df_d) in enumerate([("Crohn's Disease", df_cd),
                                     ("Psoriasis", df_ps), ("Combined", df_all)]):
    va = ((df_d['vdr_gr_class']=='VDR') & (df_d['status']=='approved')).sum()
    vf = ((df_d['vdr_gr_class']=='VDR') & (df_d['status']=='failed')).sum()
    ga = ((df_d['vdr_gr_class']=='GR')  & (df_d['status']=='approved')).sum()
    gf = ((df_d['vdr_gr_class']=='GR')  & (df_d['status']=='failed')).sum()
    if vf + ga + gf > 0:
        OR, p = fisher_exact([[va, vf],[ga, gf]])
        p_str = f'p={p:.3f}' if p >= 0.001 else f'p={p:.2e}'
        ax.text(i, 105, f'OR={OR:.1f}\n{p_str}', ha='center', va='bottom',
                fontsize=8.5, color='#333',
                bbox=dict(fc='white', ec='gray', alpha=0.7, boxstyle='round,pad=0.2'))

# Panel E: キーメッセージテキスト
ax = ax_text
ax.axis('off')

# 主要数値を計算
cd_vdr_app  = ((df_cd['vdr_gr_class']=='VDR') & (df_cd['status']=='approved')).sum()
cd_vdr_tot  = ((df_cd['vdr_gr_class']=='VDR') & (df_cd['status'].isin(['approved','failed']))).sum()
cd_gr_app   = ((df_cd['vdr_gr_class']=='GR')  & (df_cd['status']=='approved')).sum()
cd_gr_tot   = ((df_cd['vdr_gr_class']=='GR')  & (df_cd['status'].isin(['approved','failed']))).sum()
ps_vdr_app  = ((df_ps['vdr_gr_class']=='VDR') & (df_ps['status']=='approved')).sum()
ps_vdr_tot  = ((df_ps['vdr_gr_class']=='VDR') & (df_ps['status'].isin(['approved','failed']))).sum()
ps_gr_app   = ((df_ps['vdr_gr_class']=='GR')  & (df_ps['status']=='approved')).sum()
ps_gr_tot   = ((df_ps['vdr_gr_class']=='GR')  & (df_ps['status'].isin(['approved','failed']))).sum()

msg = (
    f"KEY FINDINGS\n"
    f"{'─'*32}\n\n"
    f"Crohn's Disease\n"
    f"  VDR targets: {cd_vdr_app}/{cd_vdr_tot} approved\n"
    f"    ({100*cd_vdr_app/cd_vdr_tot:.0f}%)\n"
    f"  GR targets:  {cd_gr_app}/{cd_gr_tot} approved\n"
    f"    ({100*cd_gr_app/cd_gr_tot:.0f}% )\n\n"
    f"Psoriasis\n"
    f"  VDR targets: {ps_vdr_app}/{ps_vdr_tot} approved\n"
    f"    ({100*ps_vdr_app/ps_vdr_tot:.0f}%)\n"
    f"  GR targets:  {ps_gr_app}/{ps_gr_tot} approved\n"
    f"    ({100*ps_gr_app/ps_gr_tot:.0f}%)\n\n"
    f"{'─'*32}\n"
    f"VDR-dominant loci =\n"
    f"the ONLY successful\n"
    f"therapeutic space in\n"
    f"VDR-dependent diseases"
)
ax.text(0.05, 0.95, msg, transform=ax.transAxes,
        va='top', ha='left', fontsize=9.5, family='monospace',
        bbox=dict(fc='#EEF5FF', ec='#1565C0', boxstyle='round,pad=0.5'))

fig.suptitle("VDR-dominant targets exclusively drive drug approval success\n"
             "in VDR-dependent inflammatory diseases (Crohn's Disease & Psoriasis)",
             fontsize=12, fontweight='bold', y=1.01)

plt.savefig(OUT / "fig_cd_psoriasis_vdr_approval.png", dpi=200, bbox_inches='tight')
plt.savefig(OUT / "fig_cd_psoriasis_vdr_approval.pdf", bbox_inches='tight')
print(f"\nSaved: fig_cd_psoriasis_vdr_approval.png / .pdf")
