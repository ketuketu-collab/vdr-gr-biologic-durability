#!/usr/bin/env python3
"""
誘導期 vs 維持期: VDR/GR 二相治療モデル
=========================================
コンセプト:
  急性期（誘導）: GR軸も含め炎症を制圧（どちらも効く）
  慢性期（維持）: VDR軸で維持（GR→CYP24A1→VitD枯渇を避ける）

現行治療が既にこのモデルに沿っているか、乖離しているかを可視化。
喘息/ADも含めた全疾患横断。

2026-05-23
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyArrowPatch
from pathlib import Path

OUT = Path("/Volumes/M4_SSD/projects/tlr_chipseq/results")

# ─── Tier色 ──────────────────────────────────────────────────────────────
T_COL  = {1: "#1565C0", 2: "#42A5F5", 3: "#E53935", 4: "#757575"}
T_LABEL= {1: "T1 VDR直接", 2: "T2 VDR下流", 3: "T3 GR回路", 4: "T4 回路外"}

# ─── 疾患 × 治療フェーズ データ ──────────────────────────────────────────
# (drug, gene, tier, lt_status, phase_note)
# phase: "ind"=誘導, "maint"=維持, "both"=両方, "bridge"=ステロイド橋渡し
# lt_status: "ok"=安全, "box"=FDA Black Box, "fail"=失敗, "emerg"=新興/VDR仮説

DISEASES = {

    "Crohn's\nDisease": {
        "axis": "VDR",
        "steroid": True,  # ステロイド橋渡しあり
        "induction": [
            ("infliximab",    "TNF",   1, "ok"),
            ("adalimumab",    "TNF",   1, "ok"),
            ("vedolizumab",   "ITGB7", 1, "ok"),
            ("ustekinumab",   "IL23A", 1, "ok"),
            ("risankizumab",  "IL23A", 1, "ok"),
            ("upadacitinib",  "JAK1",  3, "box"),
        ],
        "maintenance": [
            ("infliximab",    "TNF",   1, "ok"),
            ("adalimumab",    "TNF",   1, "ok"),
            ("vedolizumab",   "ITGB7", 1, "ok"),
            ("ustekinumab",   "IL23A", 1, "ok"),
            ("risankizumab",  "IL23A", 1, "ok"),
            ("upadacitinib",  "JAK1",  3, "box"),
            ("ozanimod",      "S1PR1", 4, "ok"),
        ],
        "vitd_note": "VitD欠乏がCD根本病因\n→VitD+Hypoacyl LPS = 根治戦略",
    },

    "Ulcerative\nColitis": {
        "axis": "VDR",
        "steroid": True,
        "induction": [
            ("infliximab",    "TNF",   1, "ok"),
            ("golimumab",     "TNF",   1, "ok"),
            ("vedolizumab",   "ITGB7", 1, "ok"),
            ("ustekinumab",   "IL23A", 1, "ok"),
            ("mirikizumab",   "IL23A", 1, "ok"),
            ("tofacitinib",   "JAK1",  3, "box"),
            ("upadacitinib",  "JAK1",  3, "box"),
        ],
        "maintenance": [
            ("infliximab",    "TNF",   1, "ok"),
            ("golimumab",     "TNF",   1, "ok"),
            ("vedolizumab",   "ITGB7", 1, "ok"),
            ("ustekinumab",   "IL23A", 1, "ok"),
            ("mirikizumab",   "IL23A", 1, "ok"),
            ("tofacitinib",   "JAK1",  3, "box"),
            ("upadacitinib",  "JAK1",  3, "box"),
            ("ozanimod",      "S1PR1", 4, "ok"),
            ("etrasimod",     "S1PR1", 4, "ok"),
        ],
        "vitd_note": "IBD全般にVitD欠乏関連\nVitD補充で再燃リスク低下(観察研究)",
    },

    "Psoriasis": {
        "axis": "VDR",
        "steroid": False,  # 全身ステロイドは乾癬では基本使わない（リバウンド）
        "induction": [
            ("ustekinumab",    "IL23A", 1, "ok"),
            ("guselkumab",     "IL23A", 1, "ok"),
            ("risankizumab",   "IL23A", 1, "ok"),
            ("tildrakizumab",  "IL23A", 1, "ok"),
            ("infliximab",     "TNF",   1, "ok"),
            ("adalimumab",     "TNF",   1, "ok"),
            ("etanercept",     "TNF",   1, "ok"),
            ("secukinumab",    "IL17A", 2, "ok"),
            ("ixekizumab",     "IL17A", 2, "ok"),
            ("bimekizumab",    "IL17A", 2, "ok"),
            ("deucravacitinib","TYK2",  3, "ok"),
        ],
        "maintenance": [
            ("risankizumab",   "IL23A", 1, "ok"),
            ("guselkumab",     "IL23A", 1, "ok"),
            ("ustekinumab",    "IL23A", 1, "ok"),
            ("adalimumab",     "TNF",   1, "ok"),
            ("secukinumab",    "IL17A", 2, "ok"),
            ("ixekizumab",     "IL17A", 2, "ok"),
            ("deucravacitinib","TYK2",  3, "ok"),
        ],
        "vitd_note": "局所VitD製剤（カルシポトリオール）\n= 乾癬治療の柱。全身VDR維持が根本",
    },

    "RA": {
        "axis": "VDR",
        "steroid": True,
        "induction": [
            ("infliximab",    "TNF",   1, "ok"),
            ("adalimumab",    "TNF",   1, "ok"),
            ("etanercept",    "TNF",   1, "ok"),
            ("golimumab",     "TNF",   1, "ok"),
            ("tocilizumab",   "IL6R",  1, "ok"),
            ("sarilumab",     "IL6R",  1, "ok"),
            ("abatacept",     "CD86",  1, "ok"),
            ("tofacitinib",   "JAK1",  3, "box"),
            ("baricitinib",   "JAK1",  3, "box"),
            ("upadacitinib",  "JAK1",  3, "box"),
        ],
        "maintenance": [
            ("infliximab",    "TNF",   1, "ok"),
            ("adalimumab",    "TNF",   1, "ok"),
            ("etanercept",    "TNF",   1, "ok"),
            ("tocilizumab",   "IL6R",  1, "ok"),
            ("abatacept",     "CD86",  1, "ok"),
            ("tofacitinib",   "JAK1",  3, "box"),
            ("baricitinib",   "JAK1",  3, "box"),
        ],
        "vitd_note": "ORAL Surveillance: 抗TNF(T1)が\nJAK(T3)より長期安全性優位",
    },

    "Asthma": {
        "axis": "GR",   # GR軸疾患だが…
        "steroid": True,  # 急性増悪に全身ステロイド
        "induction": [
            ("mepolizumab",   "IL5",   3, "ok"),
            ("benralizumab",  "IL5RA", 3, "ok"),
            ("dupilumab",     "IL4R",  3, "ok"),
            ("tezepelumab",   "TSLP",  1, "ok"),   # ← VDR! 上気道VDR軸
            ("omalizumab",    "IGHE",  4, "ok"),
        ],
        "maintenance": [
            ("mepolizumab",   "IL5",   3, "ok"),
            ("benralizumab",  "IL5RA", 3, "ok"),
            ("dupilumab",     "IL4R",  3, "ok"),
            ("tezepelumab",   "TSLP",  1, "ok"),   # ← VDR維持の候補
            ("omalizumab",    "IGHE",  4, "ok"),
        ],
        "vitd_note": "VitD補充→喘息増悪↓\n(Jolliffe et al. Lancet Resp Med 2021)\ntezepelumab(T1/TSLP)が全表現型で優位",
    },

    "Atopic\nDermatitis": {
        "axis": "GR",
        "steroid": True,
        "induction": [
            ("dupilumab",    "IL4R",  3, "ok"),
            ("tralokinumab", "IL13",  3, "ok"),
            ("lebrikizumab", "IL13",  3, "ok"),
            ("abrocitinib",  "JAK1",  3, "box"),
            ("upadacitinib", "JAK1",  3, "box"),
            ("baricitinib",  "JAK1",  3, "box"),
        ],
        "maintenance": [
            ("dupilumab",    "IL4R",  3, "ok"),    # ◎ 安全・長期
            ("tralokinumab", "IL13",  3, "ok"),
            ("lebrikizumab", "IL13",  3, "ok"),
            ("abrocitinib",  "JAK1",  3, "box"),   # ⚠️ Black Box
            ("upadacitinib", "JAK1",  3, "box"),   # ⚠️ Black Box
        ],
        "vitd_note": "VDR直接標的薬なし(現時点)\nVitD欠乏=AD重症化関連\nVDR維持薬は未開拓",
    },

    "Rhinitis/\nCRSwNP": {
        "axis": "GR",
        "steroid": True,
        "induction": [
            ("dupilumab",    "IL4R",  3, "ok"),
            ("mepolizumab",  "IL5",   3, "ok"),
            ("omalizumab",   "IGHE",  4, "ok"),
            ("tezepelumab",  "TSLP",  1, "ok"),   # ← VDR 2024承認
        ],
        "maintenance": [
            ("dupilumab",    "IL4R",  3, "ok"),
            ("omalizumab",   "IGHE",  4, "ok"),
            ("tezepelumab",  "TSLP",  1, "ok"),   # ← VDR 最新
        ],
        "vitd_note": "tezepelumab(T1/TSLP=VDR)が\n2024年CRSwNP承認 → VDR維持の先行例",
    },
}

# ─── VDR率計算 ────────────────────────────────────────────────────────────
def vdr_rate(drug_list):
    t1 = sum(1 for d in drug_list if d[2] == 1)
    total = len(drug_list)
    return t1 / total * 100 if total > 0 else 0

def box_rate(drug_list):
    box = sum(1 for d in drug_list if d[3] == "box")
    total = len(drug_list)
    return box / total * 100 if total > 0 else 0

print("=" * 60)
print("誘導期 vs 維持期 VDR率 サマリー")
print("=" * 60)
for dis, data in DISEASES.items():
    vi = vdr_rate(data["induction"])
    vm = vdr_rate(data["maintenance"])
    bi = box_rate(data["induction"])
    bm = box_rate(data["maintenance"])
    name = dis.replace("\n", " ")
    print(f"\n{name} ({data['axis']}軸疾患):")
    print(f"  誘導期  VDR率={vi:.0f}%  Black Box率={bi:.0f}%")
    print(f"  維持期  VDR率={vm:.0f}%  Black Box率={bm:.0f}%")

# ─── Figure ──────────────────────────────────────────────────────────────
ndis = len(DISEASES)
fig, axes = plt.subplots(2, ndis, figsize=(22, 11),
                          gridspec_kw={"height_ratios": [1.8, 1]})
fig.suptitle(
    "Induction vs Maintenance Phase — VDR/GR Circuit Classification\n"
    "Concept: GR for acute control  →  VDR for long-term maintenance",
    fontsize=13, fontweight="bold", y=0.99)

PHASE_LABEL = {"induction": "Induction\n(急性期制圧)", "maintenance": "Maintenance\n(長期維持)"}

for col, (dis_name, data) in enumerate(DISEASES.items()):

    ax_drug = axes[0, col]
    ax_stat = axes[1, col]

    # 背景: VDR軸(青) vs GR軸(赤)疾患
    bg = "#EEF4FD" if data["axis"] == "VDR" else "#FFF3E0"
    ax_drug.set_facecolor(bg)

    # ── 薬剤バー ──
    for pi, phase in enumerate(["induction", "maintenance"]):
        drugs = data[phase]
        n = len(drugs)
        bar_h = 0.7 / n if n > 0 else 0.7
        for di, (drug, gene, tier, lt) in enumerate(drugs):
            y = pi + (di - (n-1)/2) * bar_h
            col_val = T_COL[tier]
            alpha = 0.85
            lw = 2.5 if lt == "box" else 0.5
            ec = "#B71C1C" if lt == "box" else "white"
            ax_drug.barh(y, 1, left=0, height=bar_h*0.85,
                         color=col_val, alpha=alpha,
                         edgecolor=ec, linewidth=lw)
            fontsize = max(5.5, min(7, 60/n))
            ax_drug.text(0.03, y, f"{drug[:12]}", va="center",
                         fontsize=fontsize, color="white", fontweight="bold")
            if lt == "box":
                ax_drug.text(0.97, y, "⚠", va="center", ha="right",
                             fontsize=6.5, color="#FFD54F")

    # 区切り線
    ax_drug.axhline(1.0, color="#999", lw=1, ls="--")

    ax_drug.set_xlim(0, 1); ax_drug.set_ylim(-0.7, 1.7)
    ax_drug.set_yticks([0, 1])
    ax_drug.set_yticklabels(["Induction", "Maintenance"], fontsize=7.5)
    ax_drug.set_xticks([])
    ax_drug.set_title(dis_name, fontsize=9, fontweight="bold",
                      color="#1B5E20" if data["axis"]=="VDR" else "#E65100")

    # ステロイド橋渡し注記
    if data["steroid"]:
        ax_drug.text(0.5, 0.02, "← Steroid bridge (GR)", ha="center",
                     fontsize=6, color="#B71C1C", transform=ax_drug.transAxes,
                     style="italic")

    # ── 統計パネル ──
    vi = vdr_rate(data["induction"])
    vm = vdr_rate(data["maintenance"])
    bi = box_rate(data["induction"])
    bm = box_rate(data["maintenance"])

    x = [0, 1]
    ax_stat.bar([0], [vi], color=T_COL[1], alpha=0.75, width=0.35, label="VDR(T1)%")
    ax_stat.bar([1], [vm], color=T_COL[1], alpha=0.75, width=0.35)
    ax_stat.bar([0], [bi], bottom=[vi], color="#EF9A9A", alpha=0.8,
                width=0.35, label="BlackBox(T3)%")
    ax_stat.bar([1], [bm], bottom=[vm], color="#EF9A9A", alpha=0.8, width=0.35)

    ax_stat.set_xlim(-0.5, 1.5)
    ax_stat.set_ylim(0, 110)
    ax_stat.set_xticks([0, 1])
    ax_stat.set_xticklabels(["Ind.", "Maint."], fontsize=7.5)
    ax_stat.set_ylabel("% drugs" if col == 0 else "", fontsize=7)
    ax_stat.axhline(50, color="#ccc", lw=0.7, ls=":")

    # VDR変化方向矢印
    arrow_col = "#1565C0" if vm >= vi else "#E53935"
    delta = vm - vi
    if abs(delta) > 3:
        ax_stat.annotate("", xy=(1.2, vm+2), xytext=(1.2, vi+2),
                         arrowprops=dict(arrowstyle="->", color=arrow_col, lw=1.5))

    # VitD note
    ax_stat.text(0.5, -0.35, data["vitd_note"],
                 transform=ax_stat.transAxes,
                 ha="center", va="top", fontsize=5.8,
                 color="#2E7D32" if "VitD" in data["vitd_note"] else "#555",
                 style="italic")

# 共通凡例
legend_patches = [
    mpatches.Patch(color=T_COL[1], label="T1 VDR直接標的"),
    mpatches.Patch(color=T_COL[2], label="T2 VDR下流標的"),
    mpatches.Patch(color=T_COL[3], label="T3 GR回路標的"),
    mpatches.Patch(color=T_COL[4], label="T4 回路外"),
    mpatches.Patch(color="#EF9A9A", label="FDA Black Box⚠️"),
    mpatches.Patch(color="#EEF4FD", label="VDR軸疾患 (background)"),
    mpatches.Patch(color="#FFF3E0", label="GR軸疾患 (background)"),
]
fig.legend(handles=legend_patches, loc="lower center", ncol=7,
           fontsize=8, framealpha=0.9,
           bbox_to_anchor=(0.5, -0.01))

plt.subplots_adjust(left=0.04, right=0.99, top=0.94, bottom=0.18,
                    hspace=0.55, wspace=0.35)

for ext in ["png", "pdf"]:
    fig.savefig(OUT / f"fig_induction_maintenance_phase.{ext}",
                dpi=200, bbox_inches="tight")

print("\nSaved: fig_induction_maintenance_phase.png / .pdf")

# ─── Key insight summary ─────────────────────────────────────────────────
print("\n" + "=" * 60)
print("KEY INSIGHTS")
print("=" * 60)
print("""
VDR軸疾患（CD/UC/Ps/RA）:
  → 誘導期: VDR+GR薬 混在（ステロイド橋渡し含む）
  → 維持期: VDR(T1)薬が主役。JAK(T3)はBlack Boxで制限
  → 理想パターン「GR誘導→VDR維持」は既に現行標準治療に近い

GR軸疾患（Asthma/AD）:
  → 現行維持: GR軸生物製剤（dupilumab/mepolizumab）が安全に使える
    → 疾患軸に合致しているため安全（T3でもOK in GR-axis disease）
  → 新展開: tezepelumab(TSLP=T1/VDR)が喘息全表現型で最強維持
  → VitD補充: 喘息(Jolliffe 2021)・AD ともに補助的エビデンス

「GR誘導→VDR維持」の臨床的根拠:
  GR長期 → CYP24A1誘導 → VitD慢性枯渇 → VDR回路崩壊
  VDR維持薬 → GR-CYP24A1サイクルに入らない → 持続的寛解

JAKi(T3) in VDR疾患（CD/RA）:
  誘導は◎。しかし維持でFDA Black Box → ORAL Surveillanceで実証
  → 「誘導のみ使用、VDR薬に橋渡し」が最も合理的
""")
