"""
Multi-Disease VDR/GR Pharmacology Landscape
============================================
拡張解析: CD / UC / Ps / RA / Asthma / Rhinitis / AD
候補疾患: irAE / IgG4-RD / SLE / IgAN
2026-05-23

スコア規則 (cd_psoriasis_honest.py と統一):
  VDR = thp1_vdr  (THP-1 MACS2)
  GR  = gr_beas2b (BEAS-2B) — THP-1内にある場合; なければ ReMap GR
  両方ゼロ → "unknown" (Fisher除外)
"""

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['font.family'] = 'DejaVu Sans'
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.patches import Patch, FancyBboxPatch
from scipy.stats import fisher_exact
from pathlib import Path

OUT  = Path("/Volumes/M4_SSD/projects/tlr_chipseq/results")
DATA = Path("/Volumes/M4_SSD/projects/tlr_chipseq/results")

# ── スコア読み込み ─────────────────────────────────────────────────────────
thp1_df  = pd.read_csv(DATA / "thp1_ls180_scores.csv")
remap_df = pd.read_csv(DATA / "steroid_targets_vdr_gr_scores.csv")

# THP-1優先。なければReMap
GENE_SCORES = {}
for _, r in remap_df.iterrows():
    GENE_SCORES[r['gene']] = (r['VDR'], r['GR'], "ReMap")
for _, r in thp1_df.iterrows():
    GENE_SCORES[r['gene']] = (r['thp1_vdr'], r['gr_beas2b'], "THP-1")

# 手動補完・上書き (cd_psoriasis_honest.py と同じ)
GENE_SCORES.update({
    "IL23A":   (211.05,  70.34, "THP-1"),
    "ITGB7":   ( 62.27,  17.23, "ReMap"),
    "TNF":     ( 36.86,  19.49, "ReMap"),   # ReMap使用(THP-1 VDR=0)
    "IL6R":    ( 99.59,  82.90, "ReMap"),
    "S1PR1":   (  6.26,  29.35, "ReMap"),
    "OSMR":    ( 39.93, 182.90, "ReMap"),
    "JAK1":    (  0.00,  85.48, "ReMap"),
    "IL17A":   (  0.00,   0.00, "none"),
    "IL22":    (  0.00,   0.00, "none"),
    "IL5RA":   ( 20.35,  61.24, "ReMap"),
    "CXCR2":   ( 69.42,   8.12, "ReMap"),
    "TNFSF13B":(  0.00,   0.00, "none"),   # BAFF: データなし
    "MS4A1":   (  0.00,   8.92, "ReMap"),  # CD20: GR優位
})

def vdr_class(gene):
    if gene not in GENE_SCORES:
        return "unknown", 0, 0
    v, g, _ = GENE_SCORES[gene]
    if v == 0 and g == 0:
        return "unknown", 0, 0
    return ("VDR" if v > g else "GR"), v, g

# ─────────────────────────────────────────────────────────────────────────────
# 薬剤リスト定義
# ─────────────────────────────────────────────────────────────────────────────
CD_DRUGS = [
    ("IL23A","approved","ustekinumab","anti-IL23"),
    ("IL23A","approved","risankizumab","anti-IL23"),
    ("IL23A","approved","mirikizumab","anti-IL23"),
    ("ITGB7","approved","vedolizumab","anti-α4β7"),
    ("TNF",  "approved","infliximab","anti-TNF"),
    ("TNF",  "approved","adalimumab","anti-TNF"),
    ("TNF",  "approved","certolizumab","anti-TNF"),
    ("JAK1", "approved","upadacitinib","JAK1i"),
    ("S1PR1","approved","ozanimod","S1PR1"),
    ("SMAD7","failed","mongersen","Ph3失敗"),
    ("MMP9", "failed","andecaliximab","Ph2/3失敗"),
    ("CCR9", "failed","vercirnon","Ph3失敗"),
    ("IFNG", "failed","fontolizumab","Ph2中止"),
    ("IL6R", "failed","tocilizumab","CD適応なし"),
    ("IL1B", "failed","canakinumab","CD適応なし"),
    ("OSMR", "failed","vixarelimab","Ph2失敗"),
]

UC_DRUGS = [
    ("IL23A","approved","ustekinumab","anti-IL23"),
    ("IL23A","approved","risankizumab","anti-IL23"),
    ("IL23A","approved","mirikizumab","anti-IL23"),
    ("ITGB7","approved","vedolizumab","anti-α4β7"),
    ("TNF",  "approved","infliximab","anti-TNF"),
    ("TNF",  "approved","adalimumab","anti-TNF"),
    ("TNF",  "approved","golimumab","anti-TNF"),
    ("JAK1", "approved","tofacitinib","JAK"),
    ("JAK1", "approved","upadacitinib","JAK1i"),
    ("S1PR1","approved","ozanimod","S1PR1"),
    ("S1PR1","approved","etrasimod","S1PR1"),
    ("SMAD7","failed","mongersen","Ph3失敗"),
    ("MMP9", "failed","andecaliximab","Ph2/3失敗"),
    ("CCR9", "failed","vercirnon","Ph3失敗"),
    ("OSMR", "failed","vixarelimab","Ph2失敗"),
    ("IL6R", "failed","tocilizumab","Ph2失敗"),
    ("IL1B", "failed","canakinumab","UC適応なし"),
]

PS_DRUGS = [
    ("IL23A","approved","ustekinumab","anti-IL23"),
    ("IL23A","approved","guselkumab","anti-IL23"),
    ("IL23A","approved","risankizumab","anti-IL23"),
    ("IL23A","approved","tildrakizumab","anti-IL23"),
    ("TNF",  "approved","etanercept","anti-TNF"),
    ("TNF",  "approved","infliximab","anti-TNF"),
    ("TNF",  "approved","adalimumab","anti-TNF"),
    ("IL17A","approved","secukinumab","anti-IL17A"),
    ("IL17A","approved","ixekizumab","anti-IL17A"),
    ("IL17A","approved","bimekizumab","anti-IL17A/F"),
    ("TYK2", "approved","deucravacitinib","TYK2i"),
    ("IL22", "failed","fezakinumab","Ph2 limited"),
    ("MAPK14","failed","p38i","多数失敗"),
    ("ICAM1","failed","alicaforsen","Ph3失敗"),
    ("MMP9", "failed","MMPi","失敗"),
]

# ── RA: TNF/IL6R/CD86全承認、GR(MAPK14/MMP9)失敗 ─────────────────────────
RA_DRUGS = [
    ("TNF",   "approved","infliximab","anti-TNF"),
    ("TNF",   "approved","adalimumab","anti-TNF"),
    ("TNF",   "approved","etanercept","anti-TNF"),
    ("TNF",   "approved","certolizumab","anti-TNF"),
    ("TNF",   "approved","golimumab","anti-TNF"),
    ("IL6R",  "approved","tocilizumab","anti-IL6R"),
    ("IL6R",  "approved","sarilumab","anti-IL6R"),
    ("CD86",  "approved","abatacept","CTLA4-Ig/CD86"),
    ("JAK1",  "approved","tofacitinib","JAKi"),
    ("JAK1",  "approved","baricitinib","JAK1/2i"),
    ("JAK1",  "approved","upadacitinib","JAK1i"),
    ("JAK1",  "approved","filgotinib","JAK1i"),
    ("IL1B",  "approved","anakinra","anti-IL1Ra"),      # GR: largely superseded
    ("MAPK14","failed","p38i (losmapimod等)","Ph3失敗(RA)"),
    ("MMP9",  "failed","MMPi (marimastat等)","Ph3失敗"),
    ("IL17A", "failed","secukinumab","RA Ph2/3否定"),  # no ChIP data
]

# ── Asthma ──────────────────────────────────────────────────────────────────
ASTHMA_DRUGS = [
    ("IGHE",  "approved","omalizumab","anti-IgE"),
    ("IL5",   "approved","mepolizumab","anti-IL5"),
    ("IL5",   "approved","reslizumab","anti-IL5"),
    ("IL5RA", "approved","benralizumab","anti-IL5Rα"),
    ("IL4R",  "approved","dupilumab","anti-IL4Rα"),
    ("TSLP",  "approved","tezepelumab","anti-TSLP"),
    ("IL33",  "approved","itepekimab","anti-IL33"),
    ("TNF",   "failed","golimumab","severe asthma Ph3失敗"),
    ("MAPK14","failed","p38i (various)","Ph2/3失敗"),
    ("IL17A", "failed","secukinumab","好中球性喘息 無効"),
    ("IL1B",  "failed","canakinumab","非好酸球性喘息 無効"),
    ("CXCR2", "failed","navarixin","好中球性喘息 Ph2/3失敗"),
]

# ── Allergic Rhinitis / CRSwNP ───────────────────────────────────────────────
# CCR3はTHP-1ではVDR=0,GR=24.77 → GR判定(Paper Gの225スコアと異なる)
RHINITIS_DRUGS = [
    ("IGHE",  "approved","omalizumab","anti-IgE (アレルギー性鼻炎)"),
    ("IL4R",  "approved","dupilumab","CRSwNP承認"),
    ("IL5",   "approved","mepolizumab","CRSwNP承認"),
    ("TSLP",  "approved","tezepelumab","CRSwNP承認(2024)"),
    ("CCR3",  "failed","bertilimumab","CRS Ph2失敗"),
    ("MAPK14","failed","p38i","気道炎症 Ph2/3失敗"),
    ("IL1B",  "failed","canakinumab","CRSwNP適応なし"),
]

# ── Atopic Dermatitis ─────────────────────────────────────────────────────────
# 注: IL4R/IL13が全てBEAS-2B GR優位 → GR判定。VDR承認薬なし。
AD_DRUGS = [
    ("IL4R",  "approved","dupilumab","anti-IL4Rα"),
    ("IL13",  "approved","tralokinumab","anti-IL13"),
    ("IL13",  "approved","lebrikizumab","anti-IL13"),
    ("JAK1",  "approved","abrocitinib","JAK1i"),
    ("JAK1",  "approved","upadacitinib","JAK1i"),
    ("JAK1",  "approved","baricitinib","JAK1/2i"),
    ("MAPK14","failed","p38i","AD Ph2/3失敗"),
    ("CCR3",  "failed","anti-CCR3","好酸球性AD Ph2失敗"),
    ("TNFAIP3","failed","NF-κBi","失敗"),
]

# ─────────────────────────────────────────────────────────────────────────────
# 候補疾患: 承認薬プロファイルのみ（Fisher検定困難/失敗定義不明）
# ─────────────────────────────────────────────────────────────────────────────
CANDIDATE_DISEASES = {
    "irAE\n(ICI targets)": {
        "desc": "ICIが遮断する標的のVDR/GRプロファイル\n→ 治療ギャップ軸の同定",
        "genes": [
            ("CD274/PD-L1","approved","atezolizumab","VDR=186 > GR=80"),
            ("LAG3",       "approved","relatlimab",  "VDR=46  > GR=19"),
            ("HAVCR2/TIM3","approved","cobolimab",   "VDR=55  ≈ GR=56 (両軸)"),
            ("PDCD1/PD-1", "approved","pembrolizumab","VDR=0   < GR=76 (GR優位)"),
            ("CTLA4",      "approved","ipilimumab",  "VDR=0   < GR=27 (GR優位)"),
        ],
        "note": "PD-L1/LAG3=VDR優位の免疫ブレーキ\n管理(ステロイド)=GR軸 → 軸ミスマッチ仮説",
    },
    "IgG4-RD\n(Pipeline)": {
        "desc": "現行治療=GR(ステロイド)+CD20(rituximab)\nパイプライン候補がVDR優位",
        "genes": [
            ("CD19",     "pipeline","inebilizumab","VDR=84  >> GR=0  (VDR-only)"),
            ("CD38",     "pipeline","daratumumab", "VDR=95  > GR=19 (VDR優位)"),
            ("MS4A1/CD20","approved","rituximab",   "VDR=0   < GR=9  (GR優位・現行標準)"),
        ],
        "note": "次世代薬(CD19/CD38)がVDR優位\n→ 承認されれば IgG4 も VDR-positive に移行予測",
    },
    "SLE": {
        "desc": "IFNAR1(anifrolumab承認)がGR優位\n→ 本フレームワーク適用限界あり",
        "genes": [
            ("IFNAR1",    "approved","anifrolumab", "VDR=11  < GR=44 (GR優位・例外)"),
            ("CD38",      "pipeline","daratumumab", "VDR=95  > GR=19 (VDR優位・試験中)"),
            ("TNFSF13B",  "approved","belimumab",   "VDR=0   GR=0   (スコアなし)"),
        ],
        "note": "pDC特異的IFN軸が主役 → THP-1/BEAS-2Bスコアの限界\n但しVitD欠乏と疾患活動性の相関は強い",
    },
    "IgA腎症\n(IgAN)": {
        "desc": "腸内細菌-腸管-腎臓軸疾患\n→ VDR依存の根拠はあるが承認薬マッピング困難",
        "genes": [
            ("C5AR1",    "approved","iptacopan",  "VDR=0/33(THP1/ReMap) < GR=56 (GR優位)"),
            ("IFNAR1",   "pipeline","iptacopan等","VDR=11 < GR=44"),
            ("TNFSF13B", "pipeline","atacicept",  "スコアなし"),
        ],
        "note": "腸管VDR-microbiome-IgA軸は疾患根本に関与\n但し現行承認薬は補体/免疫抑制標的 → 将来VDR研究が必要",
    },
}

# ─────────────────────────────────────────────────────────────────────────────
# データフレーム構築・Fisher検定
# ─────────────────────────────────────────────────────────────────────────────
def build_df(drugs, disease):
    rows = []
    for gene, status, drug, note in drugs:
        cls, v, g = vdr_class(gene)
        rows.append(dict(disease=disease, gene=gene, drug=drug,
                         status=status, note=note,
                         vdr_class=cls, vdr_score=v, gr_score=g))
    return pd.DataFrame(rows)

TESTED = [
    ("Crohn's Disease",    CD_DRUGS),
    ("Ulcerative Colitis", UC_DRUGS),
    ("Psoriasis",          PS_DRUGS),
    ("RA",                 RA_DRUGS),
    ("Asthma",             ASTHMA_DRUGS),
    ("Rhinitis/CRSwNP",    RHINITIS_DRUGS),
    ("Atopic Dermatitis",  AD_DRUGS),
]

dfs    = {}
result = {}

print("=" * 72)
print("Multi-Disease VDR/GR Pharmacology Landscape")
print("スコア: THP-1 VDR vs BEAS-2B GR (統一判定基準)")
print("=" * 72)

for name, drugs in TESTED:
    df = build_df(drugs, name)
    dfs[name] = df

    va = int(((df['vdr_class']=='VDR')&(df['status']=='approved')).sum())
    vf = int(((df['vdr_class']=='VDR')&(df['status']=='failed')).sum())
    ga = int(((df['vdr_class']=='GR') &(df['status']=='approved')).sum())
    gf = int(((df['vdr_class']=='GR') &(df['status']=='failed')).sum())
    unk_a = int(((df['vdr_class']=='unknown')&(df['status']=='approved')).sum())
    unk_f = int(((df['vdr_class']=='unknown')&(df['status']=='failed')).sum())

    vdr_rate = 100*va/(va+vf) if (va+vf)>0 else float('nan')
    gr_rate  = 100*ga/(ga+gf) if (ga+gf)>0 else float('nan')

    # Fisher: VDRの承認優位性を検定 (one-sided, VDR > GR)
    if (va+vf)>0 and (ga+gf)>0:
        OR, p = fisher_exact([[va,vf],[ga,gf]], alternative='greater')
    else:
        OR, p = float('nan'), float('nan')

    result[name] = dict(va=va, vf=vf, ga=ga, gf=gf,
                        unk_a=unk_a, unk_f=unk_f,
                        vdr_rate=vdr_rate, gr_rate=gr_rate, OR=OR, p=p)

    sig = ("★★★" if p<0.001 else "★★" if p<0.01 else "★" if p<0.05 else "NS") if not np.isnan(p) else "--"
    or_s = "∞" if (not np.isnan(OR) and np.isinf(OR)) else (f"{OR:.1f}" if not np.isnan(OR) else "--")
    p_s  = f"{p:.4f}" if not np.isnan(p) else "--"
    print(f"\n── {name} ──")
    print(f"  VDR: 承認={va} 失敗={vf}  →  {vdr_rate:.0f}%")
    print(f"  GR:  承認={ga} 失敗={gf}  →  {gr_rate:.0f}%")
    print(f"  Unk: 承認={unk_a} 失敗={unk_f} (Fisher除外)")
    print(f"  OR={or_s}, p={p_s}  {sig}")

# 保存
res_df = pd.DataFrame(result).T.reset_index().rename(columns={'index':'disease'})
res_df.to_csv(OUT/"multi_disease_gradient.csv", index=False)

all_df = pd.concat([dfs[n] for n,_ in TESTED], ignore_index=True)
all_df.to_csv(OUT/"multi_disease_drug_list.csv", index=False)

# ─────────────────────────────────────────────────────────────────────────────
# 図の作成
# ─────────────────────────────────────────────────────────────────────────────
COLOR_VDR = '#1565C0'
COLOR_GR  = '#B71C1C'
COLOR_NS  = '#546E7A'
COLOR_CAND= '#78909C'

fig = plt.figure(figsize=(24, 18))
gs  = gridspec.GridSpec(3, 4, figure=fig, hspace=0.58, wspace=0.42,
                        height_ratios=[1.4, 1.0, 1.0])

# Panel A: Disease landscape volcano (メイン)
ax_land = fig.add_subplot(gs[0, :3])
# Panel B: 凡例・候補疾患サマリー
ax_cand = fig.add_subplot(gs[0, 3])
# Panels C-I: waterfall (7疾患)
ax_list = [fig.add_subplot(gs[1+(i//4), i%4]) for i in range(7)]

# ── Panel A: Disease Landscape ─────────────────────────────────────────────
ax = ax_land
disease_names = list(result.keys())

for dname in disease_names:
    r = result[dname]
    vr   = r['vdr_rate']
    pval = r['p']
    OR   = r['OR']
    va, vf, ga, gf = r['va'], r['vf'], r['ga'], r['gf']

    if np.isnan(vr) or np.isnan(pval):
        continue

    y = -np.log10(pval) if pval > 0 else 5.0
    x = vr

    # 色決定
    if pval < 0.05 and vr > 50:
        color = COLOR_VDR
        lw = 2.5
    elif pval < 0.05 and vr < 50:
        color = COLOR_GR
        lw = 2.5
    elif vr >= 75:
        color = '#90CAF9'   # 方向一致だがNS
        lw = 1.5
    else:
        color = COLOR_NS
        lw = 1.0

    n_scored = va+vf+ga+gf
    size = 150 + n_scored * 25

    ax.scatter(x, y, s=size, c=color, edgecolors='white', lw=lw,
               zorder=4, alpha=0.85)

    # ラベル位置調整
    offx, offy = 2.5, 0.05
    ha = 'left'
    if dname in ["Atopic Dermatitis"]:
        offy = -0.18; ha = 'center'
    if dname in ["Rhinitis/CRSwNP"]:
        offx = -3; ha = 'right'
    ax.annotate(dname.replace("'s"," "),
                xy=(x, y), xytext=(offx, offy),
                textcoords='offset points',
                fontsize=10.5,
                fontweight='bold' if pval<0.05 else 'normal',
                ha=ha, va='center',
                color='black')

# 有意水準ライン
ax.axhline(-np.log10(0.05), color='gray', ls='--', lw=1.0, alpha=0.7,
           label='p = 0.05')
ax.axvline(50, color='gray', ls=':', lw=0.8, alpha=0.5)

# 領域ラベル
ax.fill_betweenx([0, -np.log10(0.05)], 50, 105,
                 color='#B71C1C', alpha=0.04)
ax.fill_betweenx([-np.log10(0.05), 6], 50, 105,
                 color='#1565C0', alpha=0.06)
ax.text(100, 0.15, 'GR-dominant\nNS zone', ha='right', va='bottom',
        fontsize=8.5, color=COLOR_GR, alpha=0.7)
ax.text(100, -np.log10(0.05)+0.1, 'VDR-positive\nsignificant',
        ha='right', va='bottom', fontsize=8.5, color=COLOR_VDR)
ax.text(5, -np.log10(0.05)+0.1, 'VDR-positive\nNS (underpowered?)',
        ha='left', va='bottom', fontsize=8.5, color='#90CAF9')

ax.set_xlabel('VDR-dominant target approval rate (%)', fontsize=12)
ax.set_ylabel('−log₁₀(p)  [Fisher exact, one-sided]', fontsize=12)
ax.set_title('A  VDR/GR Pharmacology Landscape\n'
             'Tested diseases: CD / UC / Psoriasis / RA / Asthma / Rhinitis / AD',
             fontsize=12, fontweight='bold')
ax.set_xlim(-5, 110)
ax.set_ylim(-0.1, 4.0)

# bubble size legend
for n, lab in [(5,"N=5"),(15,"N=15")]:
    ax.scatter([], [], s=150+n*25, c='gray', alpha=0.6, label=f'n scored={n}')
ax.legend(fontsize=8.5, loc='upper left', framealpha=0.9)

# ── Panel B: 候補疾患サマリー ─────────────────────────────────────────────
ax = ax_cand
ax.axis('off')

y_pos = 0.97
for cname, cinfo in CANDIDATE_DISEASES.items():
    ax.text(0.02, y_pos, f"◆ {cname}", transform=ax.transAxes,
            fontsize=9.5, fontweight='bold', va='top')
    y_pos -= 0.04
    # gene profile
    for gene, status, drug, note in cinfo['genes']:
        cls, _, _ = vdr_class(gene.split('/')[0])
        color = COLOR_VDR if cls == 'VDR' else (COLOR_GR if cls == 'GR' else 'gray')
        marker = '▲' if cls=='VDR' else '▼' if cls=='GR' else '●'
        ax.text(0.04, y_pos,
                f"  {marker} {gene[:10]:10s} {note[:28]}",
                transform=ax.transAxes, fontsize=7.5, va='top', color=color)
        y_pos -= 0.032
    # note
    ax.text(0.03, y_pos, cinfo['note'],
            transform=ax.transAxes, fontsize=7.0, va='top',
            color='#333', style='italic')
    y_pos -= 0.06
    ax.plot([0, 1], [y_pos+0.02, y_pos+0.02], color='#ccc', lw=0.5,
            transform=ax.transAxes)

ax.set_title('B  Candidate Diseases\n(失敗定義困難/データ不足)',
             fontsize=10, fontweight='bold', loc='left')

# 凡例箱
ax.text(0.02, 0.03,
        "▲ VDR-dominant target\n▼ GR-dominant target\n● Unknown/no data",
        transform=ax.transAxes, fontsize=8, va='bottom',
        bbox=dict(fc='#f5f5f5', ec='gray', boxstyle='round,pad=0.4'))

# ── Panels C-I: Waterfall ────────────────────────────────────────────────────
def waterfall(ax, df, title, result_dict):
    df = df.copy()
    df['sv'] = df['vdr_class'].map({'VDR':0,'unknown':1,'GR':2})
    df_s = df.sort_values(['status','sv','vdr_score'],
                          ascending=[True,True,False])

    labels = [f"{r['gene']}\n({r['drug'][:9]})" for _,r in df_s.iterrows()]
    scores = df_s['vdr_score'].values
    status = df_s['status'].values
    cls    = df_s['vdr_class'].values

    cols = []
    for s, c in zip(status, cls):
        if s == 'approved':
            cols.append(COLOR_VDR if c=='VDR' else '#FF6F00' if c=='GR' else '#43A047')
        else:
            cols.append(COLOR_GR  if c=='GR'  else '#EF9A9A' if c=='VDR' else '#9E9E9E')

    ax.barh(range(len(df_s)), scores+0.5, color=cols, alpha=0.82, edgecolor='none')
    ax.set_yticks(range(len(df_s)))
    ax.set_yticklabels(labels, fontsize=5.5)
    ax.set_xlabel('VDR score', fontsize=7.5)

    r = result_dict
    p = r['p']; OR = r['OR']
    or_s = "∞" if (not np.isnan(OR) and np.isinf(OR)) else (f"OR={OR:.1f}" if not np.isnan(OR) else "")
    p_s  = f"p={p:.3f}" if not np.isnan(p) else ""
    sig  = "★" if (not np.isnan(p) and p<0.05) else ""
    vr   = r['vdr_rate']; gr = r['gr_rate']
    vr_s = f"{vr:.0f}%" if not np.isnan(vr) else "NA"
    gr_s = f"{gr:.0f}%" if not np.isnan(gr) else "NA"

    ax.set_title(f"{title}\nVDR:{vr_s} GR:{gr_s}  {or_s} {p_s}{sig}",
                 fontsize=8, fontweight='bold')

    for i, (s, sc) in enumerate(zip(status, scores)):
        ax.text(sc+0.8, i, '★' if s=='approved' else '✗',
                va='center', fontsize=7,
                color=COLOR_VDR if s=='approved' else COLOR_GR)

panel_labels = list("CDEFGHI")
disease_titles = ["C  CD", "D  UC", "E  Psoriasis",
                  "F  RA", "G  Asthma", "H  Rhinitis/CRSwNP", "I  Atopic Derm."]
for i, (dname, _) in enumerate(TESTED):
    waterfall(ax_list[i], dfs[dname], disease_titles[i], result[dname])

fig.suptitle(
    "VDR/GR Pharmacology Landscape: From Proven (IBD/RA) to Candidate (irAE/IgG4/SLE/IgAN)\n"
    "VDR-dominant target enrichment among approved biologics — cell-type context determines framework applicability",
    fontsize=13, fontweight='bold', y=1.01)

plt.savefig(OUT/"fig_multi_disease_landscape.png", dpi=200, bbox_inches='tight')
plt.savefig(OUT/"fig_multi_disease_landscape.pdf",            bbox_inches='tight')
print("\nSaved: fig_multi_disease_landscape.png / .pdf")

# ── 最終サマリー ──────────────────────────────────────────────────────────
print("\n" + "=" * 72)
print("FINAL SUMMARY")
print("=" * 72)
print(f"\n{'Disease':25s} {'VDR%':>6} {'GR%':>6} {'OR':>8} {'p':>8} {'sig':>5}")
print("-" * 72)
for name in [n for n,_ in TESTED]:
    r = result[name]
    vr_s = f"{r['vdr_rate']:.0f}" if not np.isnan(r['vdr_rate']) else " --"
    gr_s = f"{r['gr_rate']:.0f}"  if not np.isnan(r['gr_rate'])  else " --"
    or_s = ("∞" if np.isinf(r['OR']) else f"{r['OR']:.1f}") if not np.isnan(r['OR']) else "--"
    p_s  = f"{r['p']:.4f}" if not np.isnan(r['p']) else "  --"
    sig  = ("★★" if r['p']<0.01 else "★" if r['p']<0.05 else "NS") if not np.isnan(r['p']) else "--"
    print(f"  {name:25s} {vr_s:>5}%  {gr_s:>5}%  {or_s:>7}  {p_s:>8}  {sig}")

print("\n候補疾患 (Fisher検定なし):")
for cname in CANDIDATE_DISEASES:
    print(f"  {cname.replace(chr(10),' '):30s} → VDR profile analysis only")

print("\nSaved: multi_disease_gradient.csv")
print("       multi_disease_drug_list.csv")
print("       fig_multi_disease_landscape.png/pdf")
