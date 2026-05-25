#!/usr/bin/env python3
"""
治療タイムライン: GR誘導 → VDR維持 二相モデル
=================================================
コンセプト:
  急性期: GR軸薬（ステロイド/JAKi）で迅速に制圧
  移行期: GR漸減しながらVDR薬が効果を発揮
  維持期: VDR薬のみで長期寛解

機序:
  GR長期 → CYP24A1誘導 → VitD慢性枯渇 → VDR回路崩壊
  → 早期にGRを離脱してVDR薬に移行することが合理的

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

# ─── 疾患データ ───────────────────────────────────────────────────────────
# GR薬とVDR薬の使用期間（週）と、理想・現行の違い
# gr_end: GR薬を離脱すべき週（理想）
# vdr_start: VDR薬開始週
# vdr_drug: 維持期のVDR薬
# gr_drug: 誘導期のGR薬
# status: "ideal"=理想パターン実現済, "partial"=一部, "gap"=VDR薬なし

DISEASES = [
    {
        "name": "Crohn's\nDisease",
        "axis": "VDR",
        "gr_drug": "Steroid (bridge)",
        "vdr_drug": "anti-IL23 / anti-TNF\nvedolizumab",
        "gr_bridge": (0, 10),      # ステロイド使用期間(週)
        "jak_induction": (0, 12),  # JAKi 誘導使用（任意）
        "vdr_range": (2, 104),     # VDR薬使用期間
        "gr_maint_risk": True,     # GR長期使用のリスクあり（JAKi維持がある）
        "status": "ideal",
        "note": "現行標準治療に近い\nJAKi維持→Black Box⚠️",
        "vitd_opport": True,
    },
    {
        "name": "Ulcerative\nColitis",
        "axis": "VDR",
        "gr_drug": "Steroid / JAKi",
        "vdr_drug": "anti-IL23 / anti-TNF\nvedolizumab",
        "gr_bridge": (0, 10),
        "jak_induction": (0, 12),
        "vdr_range": (2, 104),
        "gr_maint_risk": True,
        "status": "ideal",
        "note": "S1PR1薬(T4)も維持で使用",
        "vitd_opport": True,
    },
    {
        "name": "Psoriasis",
        "axis": "VDR",
        "gr_drug": "(全身ステロイド\n基本使わない)",
        "vdr_drug": "anti-IL23\n(anti-TNF/anti-IL17A)",
        "gr_bridge": None,         # ステロイド橋渡しなし
        "jak_induction": None,
        "vdr_range": (0, 104),
        "gr_maint_risk": False,
        "status": "ideal",
        "note": "最もクリーンなVDR維持\n局所VitD(カルシポトリオール)も使用",
        "vitd_opport": True,
    },
    {
        "name": "RA",
        "axis": "VDR",
        "gr_drug": "Low-dose steroid\n+ MTX",
        "vdr_drug": "anti-TNF / anti-IL6R\nabatacept",
        "gr_bridge": (0, 12),
        "jak_induction": (0, 12),
        "vdr_range": (0, 104),
        "gr_maint_risk": True,
        "status": "partial",
        "note": "JAKi維持→ORAL Surveillance⚠️\n抗TNF長期維持が最も安全",
        "vitd_opport": True,
    },
    {
        "name": "Asthma",
        "axis": "GR",
        "gr_drug": "全身ステロイド(急性増悪)\n吸入ステロイド(維持)",
        "vdr_drug": "tezepelumab(TSLP/T1)\n← VDR軸・新登場",
        "gr_bridge": (0, 2),       # 急性増悪時のみ全身ステロイド
        "jak_induction": None,
        "vdr_range": (0, 104),     # tezepelumabは持続使用
        "gr_maint_risk": False,    # 吸入ステロイドは局所→全身影響小
        "status": "partial",
        "note": "吸入GRは維持でも使用\ntezepelumab(T1)が最強維持へ移行中\nVitD補充→増悪↓(Jolliffe 2021)",
        "vitd_opport": True,
    },
    {
        "name": "Atopic\nDermatitis",
        "axis": "GR",
        "gr_drug": "外用ステロイド\n(+全身ステロイド)",
        "vdr_drug": "dupilumab(IL4R/T3)\n← VDR薬は現時点でなし",
        "gr_bridge": (0, 4),
        "jak_induction": (0, 12),  # JAKi 誘導
        "vdr_range": None,         # VDR直接標的薬なし
        "gr_maint_risk": True,
        "status": "gap",
        "note": "VDR直接維持薬なし→空白領域\ndupilumab(T3)が最良維持だが\nJAKi(T3)維持はBlack Box⚠️",
        "vitd_opport": True,
    },
    {
        "name": "Rhinitis/\nCRSwNP",
        "axis": "GR",
        "gr_drug": "鼻腔内ステロイド\n(局所GR)",
        "vdr_drug": "tezepelumab(TSLP/T1)\n2024年CRSwNP承認",
        "gr_bridge": (0, 104),     # 局所ステロイドは継続（全身影響小）
        "jak_induction": None,
        "vdr_range": (0, 104),
        "gr_maint_risk": False,
        "status": "partial",
        "note": "局所GR+tezepelumab(T1)の\n併用維持が新標準へ",
        "vitd_opport": True,
    },
]

# ─── Figure ──────────────────────────────────────────────────────────────
fig = plt.figure(figsize=(22, 14))
fig.suptitle(
    "Treatment Timeline: GR-axis (Induction)  →  VDR-axis (Maintenance)\n"
    "Mechanistic rationale: GR long-term → CYP24A1 → VitD depletion → VDR circuit failure",
    fontsize=13, fontweight="bold", y=0.99)

# グリッド: 上段=タイムライン、下段=説明
from matplotlib.gridspec import GridSpec
ndis = len(DISEASES)
gs = GridSpec(3, ndis, figure=fig,
              height_ratios=[3.5, 1, 0.8],
              hspace=0.6, wspace=0.3,
              left=0.04, right=0.99, top=0.93, bottom=0.04)

# 時間軸
T_MAX = 104  # 週
T_TICKS = [0, 8, 12, 26, 52, 104]
T_LABELS = ["0", "8w", "12w", "26w", "1yr", "2yr"]

C_GR    = "#E53935"   # GR薬（赤）
C_VDR   = "#1565C0"   # VDR薬（青）
C_JAK   = "#FF7043"   # JAKi（橙赤）
C_CS_I  = "#EF9A9A"   # 吸入/局所ステロイド（薄赤）
C_GAP   = "#BDBDBD"   # VDR薬なし（グレー）
C_VITD  = "#2E7D32"   # VitD（緑）

STATUS_COLOR = {"ideal": "#E8F5E9", "partial": "#FFF8E1", "gap": "#FFEBEE"}
STATUS_LABEL = {"ideal": "現行≈理想", "partial": "移行中", "gap": "VDR薬 未開拓"}

for col, dis in enumerate(DISEASES):
    ax = fig.add_subplot(gs[0, col])
    ax_note = fig.add_subplot(gs[1, col])
    ax_mech = fig.add_subplot(gs[2, col])

    # 背景
    ax.set_facecolor(STATUS_COLOR[dis["status"]])

    # ── タイムライン描画 ──────────────────────────────────────────────────
    y_positions = {"gr": 2.5, "jak": 1.5, "vdr": 0.5}

    # ステロイド橋渡し
    if dis["gr_bridge"]:
        s, e = dis["gr_bridge"]
        alpha_gr = 0.9 if e < 20 else 0.5  # 長期なら薄く
        color_gr = C_GR if e < 20 else C_CS_I
        ax.barh(y_positions["gr"], e-s, left=s, height=0.6,
                color=color_gr, alpha=alpha_gr, label="GR(steroid)")
        label_gr = "Steroid\n(bridge)" if e < 20 else "Local GR\n(inhaled/nasal)"
        ax.text(s+1, y_positions["gr"], label_gr, va="center",
                fontsize=6, color="white", fontweight="bold")
        # GR離脱矢印（短期橋渡しの場合）
        if e < 20:
            ax.annotate("", xy=(e+2, y_positions["gr"]),
                        xytext=(e, y_positions["gr"]),
                        arrowprops=dict(arrowstyle="->", color="#B71C1C", lw=1.2))
    else:
        ax.text(2, y_positions["gr"], "(ステロイド橋渡しなし)",
                va="center", fontsize=6, color="#999", style="italic")

    # JAKi（誘導期のみ推奨）
    if dis["jak_induction"]:
        s, e = dis["jak_induction"]
        ax.barh(y_positions["jak"], e-s, left=s, height=0.6,
                color=C_JAK, alpha=0.85, label="JAKi (T3)")
        ax.text(s+1, y_positions["jak"], "JAKi\n(T3, ind.)", va="center",
                fontsize=6, color="white", fontweight="bold")
        # 理想: 誘導後离脱
        ax.annotate("", xy=(e+2, y_positions["jak"]),
                    xytext=(e, y_positions["jak"]),
                    arrowprops=dict(arrowstyle="->", color="#BF360C", lw=1.2))
        ax.text(e+3, y_positions["jak"], "taper→stop\n⚠️BBox if continued",
                va="center", fontsize=5.5, color="#BF360C", style="italic")

    # VDR薬（維持）
    if dis["vdr_range"]:
        s, e = dis["vdr_range"]
        ax.barh(y_positions["vdr"], e-s, left=s, height=0.6,
                color=C_VDR, alpha=0.85, label="VDR(T1)")
        label_vdr = dis["vdr_drug"].split("\n")[0][:14]
        ax.text(s+1, y_positions["vdr"], label_vdr, va="center",
                fontsize=6, color="white", fontweight="bold")
        # 「維持」ラベル
        ax.text(T_MAX*0.7, y_positions["vdr"], "◎ 維持",
                va="center", fontsize=6.5, color="#90CAF9", fontweight="bold")
    else:
        # VDR薬なし
        ax.barh(y_positions["vdr"], T_MAX, left=0, height=0.6,
                color=C_GAP, alpha=0.4, hatch="///")
        ax.text(T_MAX/2, y_positions["vdr"], "VDR維持薬\n未開拓",
                va="center", ha="center", fontsize=7, color="#757575",
                fontweight="bold")

    # VitD補充ライン（常時・緑点線）
    if dis["vitd_opport"]:
        ax.plot([0, T_MAX], [y_positions["vdr"]-0.5, y_positions["vdr"]-0.5],
                color=C_VITD, lw=1.5, ls="--", alpha=0.7)
        ax.text(T_MAX*0.5, y_positions["vdr"]-0.5,
                "VitD補充（根本修復）", va="center", ha="center",
                fontsize=5.5, color=C_VITD, style="italic")

    # 軸
    ax.set_xlim(0, T_MAX)
    ax.set_ylim(y_positions["vdr"]-0.8, y_positions["gr"]+0.7)
    ax.set_xticks(T_TICKS)
    ax.set_xticklabels(T_LABELS, fontsize=7)
    ax.set_yticks([])
    ax.set_xlabel("Time (weeks)", fontsize=7)

    # y軸ラベル（最左列のみ）
    if col == 0:
        ax.set_yticks([y_positions["gr"], y_positions["jak"], y_positions["vdr"]])
        ax.set_yticklabels(["GR\n(steroid)", "JAKi\n(bridge)", "VDR\n(maint)"],
                            fontsize=6.5)

    # タイトル
    status_c = {"ideal":"#2E7D32","partial":"#E65100","gap":"#B71C1C"}[dis["status"]]
    ax.set_title(
        f"{dis['name']}\n[{dis['axis']}軸] {STATUS_LABEL[dis['status']]}",
        fontsize=8.5, fontweight="bold",
        color=status_c)

    # 転換ゾーン（誘導→維持）
    ax.axvspan(8, 14, alpha=0.07, color="gold", zorder=0)
    ax.text(11, y_positions["gr"]+0.5, "移行\n期", ha="center",
            fontsize=6, color="#795548")

    # ── 注記パネル ─────────────────────────────────────────────────────
    ax_note.axis("off")
    ax_note.text(0.5, 0.9, dis["note"],
                 transform=ax_note.transAxes,
                 ha="center", va="top", fontsize=7,
                 color="#333",
                 bbox=dict(fc=STATUS_COLOR[dis["status"]],
                           ec=status_c, boxstyle="round,pad=0.3",
                           alpha=0.8))

    # ── 機序パネル ─────────────────────────────────────────────────────
    ax_mech.axis("off")
    if dis["gr_maint_risk"]:
        mech_text = "GR長期→CYP24A1↑\n→VitD枯渇→VDR崩壊"
        mech_col  = "#B71C1C"
    else:
        mech_text = "GR早期離脱\n→VitD温存"
        mech_col  = "#1565C0"
    ax_mech.text(0.5, 0.7, mech_text,
                 transform=ax_mech.transAxes,
                 ha="center", va="center", fontsize=6.5,
                 color=mech_col, fontweight="bold",
                 style="italic")

# ── 共通凡例 ──────────────────────────────────────────────────────────────
legend_patches = [
    mpatches.Patch(color=C_GR,   label="GR薬（全身ステロイド）: 急性制圧"),
    mpatches.Patch(color=C_CS_I, label="局所GR（吸入/鼻腔内）: 全身影響小"),
    mpatches.Patch(color=C_JAK,  label="JAKi T3: 誘導のみ推奨（維持→Black Box⚠️）"),
    mpatches.Patch(color=C_VDR,  label="VDR薬 T1: 長期維持◎（Black Boxなし）"),
    mpatches.Patch(color=C_GAP,  label="VDR維持薬なし（創薬空白）", alpha=0.5, hatch="///"),
    mpatches.Patch(color=C_VITD, label="VitD補充（根本修復・全疾患共通）", alpha=0.7),
    mpatches.Patch(color="#FFF8E1",ec="#E65100",label="移行中（現行≠理想）"),
    mpatches.Patch(color="#FFEBEE",ec="#B71C1C",label="VDR維持未開拓（創薬機会）"),
]
fig.legend(handles=legend_patches, loc="lower center", ncol=4,
           fontsize=8, framealpha=0.9,
           bbox_to_anchor=(0.5, -0.01))

for ext in ["png", "pdf"]:
    fig.savefig(OUT / f"fig_treatment_timeline_concept.{ext}",
                dpi=200, bbox_inches="tight")

print("Saved: fig_treatment_timeline_concept.png / .pdf")

print("""
============================================================
二相治療モデル サマリー
============================================================

【理想パターン（現行標準治療に近い）】
  CD / UC / RA:
    Week 0-10:  ステロイド橋渡し（GR）+ VDR薬開始
    Week 10+:   ステロイド離脱、VDR薬（抗IL-23/抗TNF）で維持
    JAKiは誘導期のみ使用が合理的（維持継続→Black Box累積）

  Psoriasis:
    全身ステロイドなし → 最初からVDR薬で誘導＋維持
    最もクリーンなVDR治療

【移行中】
  Asthma:
    全身GR(急性増悪) → 吸入GR(局所) + tezepelumab(T1/VDR)
    VDR軸への移行が始まっている

  Rhinitis/CRSwNP:
    鼻腔内GR(局所) + tezepelumab(T1/VDR)で維持
    局所GRは全身影響小 → 許容範囲

【VDR維持薬が未開拓】
  Atopic Dermatitis:
    dupilumab(T3)が最良だが VDR直接標的薬なし
    → 創薬機会。VitD補充が補助的選択肢

【全疾患共通】
  VitD補充 = GR→CYP24A1→VitD枯渇を補う根本介入
  VDR維持薬の効果を底支えする役割
""")
