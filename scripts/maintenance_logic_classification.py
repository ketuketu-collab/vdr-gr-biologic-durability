#!/usr/bin/env python3
"""
維持療法予測: VDR/GR 作用機序ロジック分類
==========================================
ChIP-seqスコアを使わず「生物学的回路への帰属」でTier分類。

コアロジック:
  VDR欠乏（VitD欠乏 + dysbiosis）が慢性炎症の根本。
  「どの回路」の薬かによって長期維持効果が変わる。

Tier定義:
  T1: VDR回路・直接標的  — VDRが通常抑制する遺伝子を遮断
      (IL23A/TNF/ITGB7/IL6R/CD86/LAG3/PD-L1)
      → GR→CYP24A1→VitD枯渇なし → 長期維持◎

  T2: VDR回路・下流標的  — VDR標的の下流エフェクター
      (IL17A: IL23A→Th17→IL17A)
      → 末端遮断。根本VDR欠乏は残存 → 長期は疾患依存

  T3: GR回路・標的       — JAK-STAT/TYK2/S1PR経路（GR優位シグナル）
      JAK阻害薬 / TYK2i / MAPK14
      → GR軸経由で効く → 急性◎ / 長期⚠️（FDA Black Box）

  T4: 回路外             — VDR/GR双方に非依存
      IGHE(IgE中和) / S1PR1(リンパ球移動)
      → 機序中立 → 病態特異的

Short-term note: 急性期はGR軸でも効く（炎症を抑える経路は複数）
Long-term note: 維持でVDR軸かGR軸かが決定的に重要

2026-05-23
"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from pathlib import Path

OUT = Path("/Volumes/M4_SSD/projects/tlr_chipseq/results")

# ─── Tier定義 ─────────────────────────────────────────────────────────────
TIERS = {
    1: dict(label="T1\nVDR回路\n直接標的",  color="#1565C0", pred="◎予測",
            logic="VDRが抑制する遺伝子を遮断\n→ GR→CYP24A1→VitD枯渇なし"),
    2: dict(label="T2\nVDR回路\n下流標的",  color="#42A5F5", pred="△予測(疾患依存)",
            logic="IL23A→Th17→IL17A など\n→ 末端遮断・根本VDR欠乏残存"),
    3: dict(label="T3\nGR回路\n標的",       color="#E53935", pred="⚠️予測",
            logic="JAK-STAT/TYK2 GR優位経路\n→ CYP24A1誘導リスク"),
    4: dict(label="T4\n回路外",             color="#9E9E9E", pred="△予測(中立)",
            logic="S1PR1/IGHE: VDR/GR回路外\n→ 機序中立"),
}

# ─── LT(長期)ステータス定義 ────────────────────────────────────────────────
# ◎=3: 長期維持◎, 安全シグナルなし
# △=2: 長期使用可, 要モニタリング
# ⚠️=1: FDA Black Box Warning (感染/癌/CV/血栓)
# ❌=0: 臨床試験失敗 / 長期使用不可
LT_SCORE   = {"◎": 3, "△": 2, "⚠️": 1, "❌": 0}
LT_MARKER  = {"◎": "◎", "△": "△", "⚠️": "⚠", "❌": "✕"}
LT_COLOR   = {"◎": "#1B5E20", "△": "#F57F17", "⚠️": "#B71C1C", "❌": "#424242"}

# ─── 薬剤データ ─────────────────────────────────────────────────────────────
# (drug, gene, tier, diseases, approved, lt_status, lt_evidence, disease_axis)
# disease_axis: "VDR" = IBD/RA/Ps, "GR" = Asthma/AD (フレームワーク対象外)
DRUGS = [
    # ── T1: VDR回路・直接標的 ───────────────────────────────────────────────
    # IL23A (VDR=211): CD/UC/Ps の核心。VDRが抑制する最上流炎症遺伝子
    ("ustekinumab",   "IL23A", 1, "CD/UC/Ps", True,  "◎", "UNIFI長期・PSUMMIT長期◎",          "VDR"),
    ("risankizumab",  "IL23A", 1, "CD/UC/Ps", True,  "◎", "SEQUENCE 1年優位・維持承認◎",       "VDR"),
    ("mirikizumab",   "IL23A", 1, "UC/CD",    True,  "◎", "UC/CD維持承認2023-24",               "VDR"),
    ("guselkumab",    "IL23A", 1, "Ps",       True,  "◎", "Ps長期5年◎",                         "VDR"),
    ("tildrakizumab", "IL23A", 1, "Ps",       True,  "◎", "Ps長期◎",                            "VDR"),
    # ITGB7 (VDR=62): 腸管リンパ球ホーミング。VDR直接標的
    ("vedolizumab",   "ITGB7", 1, "CD/UC",    True,  "◎", "最良長期安全性・GEMINI長期",          "VDR"),
    # TNF (VDR=157): VDR直接抑制標的の代表
    ("infliximab",    "TNF",   1, "CD/UC/RA/Ps",True,"◎", "20年以上実績・ACCENT I長期維持",     "VDR"),
    ("adalimumab",    "TNF",   1, "CD/UC/RA/Ps",True,"◎", "最多長期実績・10年コホート◎",        "VDR"),
    ("certolizumab",  "TNF",   1, "CD/RA",    True,  "◎", "長期実績あり",                       "VDR"),
    ("golimumab",     "TNF",   1, "UC/RA",    True,  "◎", "UC/RA維持承認・長期◎",               "VDR"),
    ("etanercept",    "TNF",   1, "RA/Ps",    True,  "◎", "RA/Ps長期実績◎",                     "VDR"),
    # IL6R (VDR=100, GR=83): VDR優位だがGRも存在 → Tier1の中では慎重
    ("tocilizumab",   "IL6R",  1, "RA",       True,  "△", "RA長期使用可・感染+脂質注意",        "VDR"),
    ("sarilumab",     "IL6R",  1, "RA",       True,  "△", "RA承認・同クラス",                   "VDR"),
    # CD86 (VDR=27, GR=0): 完全VDR専管
    ("abatacept",     "CD86",  1, "RA",       True,  "◎", "RA長期維持◎・安全性良好",            "VDR"),
    # TSLP (VDR=24, GR=0): VDR専管
    ("tezepelumab",   "TSLP",  1, "Asthma/CRS",True, "◎", "TSLP長期◎・NAVIGATOR試験",          "GR"),

    # ── T2: VDR回路・下流標的 (IL23A→Th17→IL17A) ──────────────────────────
    # IL17A: ChIPデータなし(VDR=0,GR=0)。IL23A下流エフェクター
    # 重要: Ps◎ (皮膚IL17Aは主エフェクター) vs CD❌(腸管粘膜バリア破壊)
    ("secukinumab",   "IL17A", 2, "Ps/RA",    True,  "◎", "Ps長期◎・CD禁忌(腸管バリア破壊)",  "VDR"),
    ("ixekizumab",    "IL17A", 2, "Ps/RA",    True,  "◎", "Ps長期◎・CD禁忌",                   "VDR"),
    ("bimekizumab",   "IL17A/F",2,"Ps",       True,  "◎", "Ps長期◎(IL17A+F二重阻害)",          "VDR"),
    ("fezakinumab",   "IL22",  2, "Ps",       False, "❌", "Ps Ph2 limited・中止",               "VDR"),

    # ── T3: GR回路・標的 ───────────────────────────────────────────────────
    # JAK1 (GR=85, VDR=0): GR優位経路。FDA Black Box全JAK阻害薬
    ("upadacitinib",  "JAK1",  3, "CD/UC/RA/AD",True,"⚠️","FDA Black Box(感染/癌/CV/VTE)",    "VDR"),
    ("tofacitinib",   "JAK1",  3, "UC/RA",    True,  "⚠️", "FDA Black Box・ORAL Surveillance:CV+癌↑","VDR"),
    ("baricitinib",   "JAK1",  3, "RA/AD",    True,  "⚠️", "FDA Black Box・同クラス",           "VDR"),
    ("filgotinib",    "JAK1",  3, "RA",       True,  "⚠️", "FDA Black Box・同クラス",           "VDR"),
    ("abrocitinib",   "JAK1",  3, "AD",       True,  "⚠️", "FDA Black Box・同クラス",           "GR"),
    # TYK2 (GR=73, VDR=20): GR優位。JAK阻害薬よりは選択的
    ("deucravacitinib","TYK2", 3, "Ps",       True,  "△", "TYK2i・JAK非依存・長期安全性蓄積中","VDR"),
    # MAPK14 p38i: GR優位・多数失敗
    ("p38i_various",  "MAPK14",3, "CD/RA/Ps", False, "❌", "多数Ph3失敗・GR優位標的",            "VDR"),
    # MMP9 (GR=72): GR優位
    ("andecaliximab", "MMP9",  3, "CD/UC",    False, "❌", "Ph2/3失敗・GR優位",                  "VDR"),
    ("MMPi_various",  "MMP9",  3, "RA/Ps",    False, "❌", "Ph3失敗・GR優位",                    "VDR"),
    # CCR9, OSMR, IFNG: GR優位・失敗
    ("vercirnon",     "CCR9",  3, "CD",       False, "❌", "Ph3失敗",                            "VDR"),
    ("fontolizumab",  "IFNG",  3, "CD",       False, "❌", "Ph2中止・GR優位",                    "VDR"),
    ("vixarelimab",   "OSMR",  3, "CD/UC",    False, "❌", "Ph2失敗・OSMR GR=183(最高値)",       "VDR"),
    ("mongersen",     "SMAD7", 3, "CD",       False, "❌", "Ph3失敗",                            "VDR"),
    # IL1B (GR=45): GR優位・Still病のみ例外
    ("anakinra",      "IL1B",  3, "Still/RA", True,  "△", "Still病承認・RA限定的使用",           "VDR"),
    ("canakinumab",   "IL1B",  3, "Still/痛風",True, "△", "Still病/痛風承認・腸管VDR疾患外",    "VDR"),
    # GR軸サイトカイン生物製剤: 喘息/ADで適切（GR軸疾患）→ safety問題なし
    ("mepolizumab",   "IL5",   3, "Asthma",   True,  "◎", "IL5・好酸球性喘息長期◎",             "GR"),
    ("reslizumab",    "IL5",   3, "Asthma",   True,  "◎", "IL5・喘息承認",                      "GR"),
    ("benralizumab",  "IL5RA", 3, "Asthma",   True,  "◎", "IL5Rα・喘息長期◎",                  "GR"),
    ("dupilumab",     "IL4R",  3, "AD/Asthma",True,  "◎", "IL4Rα・AD/喘息長期安全性◎",         "GR"),
    # GR軸薬・VDR軸疾患で失敗
    ("golimumab_asth","TNF",   1, "Asthma",   False, "❌", "TNF→severe asthma Ph3失敗",          "GR"),

    # ── T4: 回路外 (VDR/GR非依存) ─────────────────────────────────────────
    ("omalizumab",    "IGHE",  4, "Asthma/CRS",True, "◎", "IgE中和・長期安全性◎",               "GR"),
    ("ozanimod",      "S1PR1", 4, "UC/MS",    True,  "△", "リンパ球移動制御・心電図監視要",      "VDR"),
    ("etrasimod",     "S1PR1", 4, "UC",       True,  "△", "S1PR1・同クラス",                    "VDR"),
    # CCR3 (GR=25, VDR=0): GR優位・20年失敗
    ("bertilimumab",  "CCR3",  3, "CRS",      False, "❌", "Ph2失敗・CCR3 GR優位（20年失敗予測通り）","GR"),
]

df = pd.DataFrame(DRUGS, columns=["drug","gene","tier","diseases","approved",
                                   "lt_status","lt_evidence","disease_axis"])
df["lt_score"]    = df["lt_status"].map(LT_SCORE)
df["tier_label"]  = df["tier"].map(lambda t: TIERS[t]["label"])
df["tier_color"]  = df["tier"].map(lambda t: TIERS[t]["color"])
df["lt_color"]    = df["lt_status"].map(LT_COLOR)

df.to_csv(OUT / "maintenance_logic_classification.csv", index=False)

# VDR軸疾患のみ抽出（フレームワーク検証対象）
vdr_df = df[df["disease_axis"] == "VDR"].copy()
gr_df  = df[df["disease_axis"] == "GR"].copy()

# ─── Print summary ─────────────────────────────────────────────────────────
print("=" * 60)
print("維持療法ロジック分類 — VDR軸疾患（CD/UC/Ps/RA）")
print("=" * 60)
for t in sorted(vdr_df["tier"].unique()):
    sub = vdr_df[vdr_df["tier"] == t]
    cnts = sub["lt_status"].value_counts()
    print(f"\nTier {t} ({TIERS[t]['pred']})  n={len(sub)}")
    for st in ["◎","△","⚠️","❌"]:
        n = cnts.get(st, 0)
        if n: print(f"  {st}: {n}")

print("\n" + "=" * 60)
print("GR軸疾患（Asthma/AD）— 参考（フレームワーク外）")
print("=" * 60)
for t in sorted(gr_df["tier"].unique()):
    sub = gr_df[gr_df["tier"] == t]
    cnts = sub["lt_status"].value_counts()
    print(f"\nTier {t}  n={len(sub)}")
    for st in ["◎","△","⚠️","❌"]:
        n = cnts.get(st, 0)
        if n: print(f"  {st}: {n}")

# ─── Figure ────────────────────────────────────────────────────────────────
fig = plt.figure(figsize=(18, 14))
fig.suptitle(
    "VDR/GR Mechanism Logic Classification for Long-term Maintenance\n"
    "Circuit membership predicts maintenance safety — not short-term efficacy",
    fontsize=13, fontweight="bold", y=0.98)

gs = fig.add_gridspec(3, 3, hspace=0.55, wspace=0.4,
                      left=0.06, right=0.97, top=0.93, bottom=0.05)

ax_logic = fig.add_subplot(gs[0, :])   # Panel A: Logic diagram
ax_main  = fig.add_subplot(gs[1:, :2]) # Panel B: Drug scatter (VDR-axis diseases)
ax_bar   = fig.add_subplot(gs[1, 2])   # Panel C: Summary bar
ax_il17  = fig.add_subplot(gs[2, 2])   # Panel D: IL17A special case

# ── Panel A: Logic tier 説明 ──────────────────────────────────────────────
ax_logic.set_xlim(0, 1); ax_logic.set_ylim(0, 1)
ax_logic.axis("off")
ax_logic.set_title("A  Circuit Logic — Maintenance Prediction by Tier",
                   fontsize=10, fontweight="bold", loc="left")

tier_box_props = dict(boxstyle="round,pad=0.5", alpha=0.15)
xpos = [0.05, 0.28, 0.52, 0.76]
for i, (t, tinfo) in enumerate(TIERS.items()):
    if t > 4: continue
    ax_logic.add_patch(mpatches.FancyBboxPatch(
        (xpos[i]-0.01, 0.05), 0.21, 0.87,
        boxstyle="round,pad=0.02", facecolor=tinfo["color"],
        alpha=0.12, edgecolor=tinfo["color"], lw=2, transform=ax_logic.transAxes))
    ax_logic.text(xpos[i]+0.1, 0.85, tinfo["label"],
                  ha="center", va="top", fontsize=9, fontweight="bold",
                  color=tinfo["color"], transform=ax_logic.transAxes)
    ax_logic.text(xpos[i]+0.1, 0.60, tinfo["pred"],
                  ha="center", va="top", fontsize=9,
                  color=LT_COLOR.get(tinfo["pred"].split("予")[0], "#333"),
                  transform=ax_logic.transAxes, fontweight="bold")
    ax_logic.text(xpos[i]+0.1, 0.42, tinfo["logic"],
                  ha="center", va="top", fontsize=7.5,
                  color="#333", transform=ax_logic.transAxes)

# 矢印: GR→CYP24A1→VitD枯渇のリスク
ax_logic.annotate(
    "GR→CYP24A1→VitD depletion risk increases →",
    xy=(0.97, 0.04), xytext=(0.05, 0.04),
    xycoords="axes fraction", textcoords="axes fraction",
    arrowprops=dict(arrowstyle="->", color="#B71C1C", lw=1.5),
    fontsize=8.5, color="#B71C1C", va="center")

# ── Panel B: Drug scatter — VDR-axis diseases ────────────────────────────
np.random.seed(42)
ax_main.set_xlim(0.5, 4.5); ax_main.set_ylim(-0.4, 3.5)
ax_main.set_xticks([1,2,3,4])
ax_main.set_xticklabels(["T1\nVDR直接", "T2\n下流標的", "T3\nGR回路", "T4\n回路外"],
                         fontsize=9)
ax_main.set_yticks([0,1,2,3])
ax_main.set_yticklabels(["❌ 失敗", "⚠️ Black Box", "△ 要注意", "◎ 長期維持◎"],
                         fontsize=9)
ax_main.set_xlabel("Circuit Tier (Logic-based classification)", fontsize=10)
ax_main.set_ylabel("Long-term Maintenance Status", fontsize=10)
ax_main.set_title("B  VDR-axis Diseases (CD / UC / Ps / RA)\n"
                  "Logic tier predicts long-term maintenance status",
                  fontsize=10, fontweight="bold", loc="left")

# 背景帯
for t, col in [(1,"#E3F2FD"), (3,"#FFEBEE")]:
    ax_main.axvspan(t-0.45, t+0.45, alpha=0.3, color=col, zorder=0)

ax_main.axhline(0.5, color="#ccc", lw=0.5, ls="--", zorder=0)
ax_main.axhline(1.5, color="#ccc", lw=0.5, ls="--", zorder=0)
ax_main.axhline(2.5, color="#ccc", lw=0.5, ls="--", zorder=0)

for _, row in vdr_df.iterrows():
    jx = row["tier"] + np.random.uniform(-0.3, 0.3)
    jy = row["lt_score"] + np.random.uniform(-0.12, 0.12)
    size = 180 if row["approved"] else 90
    marker = "o" if row["approved"] else "x"
    ax_main.scatter(jx, jy, c=row["lt_color"], s=size, marker=marker,
                    edgecolors="white" if row["approved"] else row["lt_color"],
                    linewidths=0.8, zorder=3, alpha=0.85)
    # 薬名ラベル（主要薬のみ）
    key_drugs = {"infliximab","tofacitinib","upadacitinib","ustekinumab",
                 "vedolizumab","abatacept","secukinumab","deucravacitinib",
                 "mongersen","andecaliximab","ozanimod","tocilizumab"}
    if row["drug"] in key_drugs:
        ax_main.annotate(row["drug"], (jx, jy),
                         textcoords="offset points", xytext=(5, 2),
                         fontsize=6.5, color="#333", zorder=4)

# FDA Black Box ゾーン注記
ax_main.text(3.0, 1.0, "FDA Black Box\nWarning zone",
             ha="center", va="center", fontsize=8, color="#B71C1C",
             style="italic", transform=ax_main.transData,
             bbox=dict(fc="#FFEBEE", ec="#E53935", boxstyle="round,pad=0.3", alpha=0.7))

# 凡例
from matplotlib.lines import Line2D
legend_els = [
    Line2D([0],[0], marker="o", color="w", markerfacecolor=LT_COLOR["◎"],   markersize=9, label="◎ 長期維持◎"),
    Line2D([0],[0], marker="o", color="w", markerfacecolor=LT_COLOR["△"],   markersize=9, label="△ 要注意"),
    Line2D([0],[0], marker="o", color="w", markerfacecolor=LT_COLOR["⚠️"],  markersize=9, label="⚠️ FDA Black Box"),
    Line2D([0],[0], marker="o", color="w", markerfacecolor=LT_COLOR["❌"],  markersize=9, label="❌ 失敗"),
    Line2D([0],[0], marker="o", color="w", markerfacecolor="#555", markersize=9, label="承認済"),
    Line2D([0],[0], marker="x", color="#555", markersize=8,                      label="試験失敗"),
]
ax_main.legend(handles=legend_els, fontsize=7.5, loc="upper right",
               ncol=2, framealpha=0.85)

# ── Panel C: Tier別LT比率 (VDR軸疾患) ────────────────────────────────────
tier_order = [1, 2, 3, 4]
statuses   = ["◎", "△", "⚠️", "❌"]
st_colors  = [LT_COLOR[s] for s in statuses]

counts = []
for t in tier_order:
    sub = vdr_df[vdr_df["tier"] == t]
    row = [sub[sub["lt_status"] == s].shape[0] for s in statuses]
    counts.append(row)

counts = np.array(counts)
totals = counts.sum(axis=1, keepdims=True)
totals[totals == 0] = 1
pcts = counts / totals * 100

xlabels = ["T1\nVDR直接", "T2\n下流", "T3\nGR回路", "T4\n外"]
bottom = np.zeros(len(tier_order))
for j, (st, col) in enumerate(zip(statuses, st_colors)):
    vals = pcts[:, j]
    bars = ax_bar.bar(range(len(tier_order)), vals, bottom=bottom,
                      color=col, label=st, width=0.6, alpha=0.88)
    for k, (v, b) in enumerate(zip(vals, bottom)):
        if v > 8:
            ax_bar.text(k, b + v/2, f"{v:.0f}%", ha="center", va="center",
                        fontsize=7.5, color="white", fontweight="bold")
    bottom += vals

ax_bar.set_xticks(range(len(tier_order))); ax_bar.set_xticklabels(xlabels, fontsize=8)
ax_bar.set_ylabel("Long-term status (%)", fontsize=8)
ax_bar.set_ylim(0, 110)
ax_bar.set_title("C  VDR-axis Diseases:\nLT Status by Tier", fontsize=9,
                 fontweight="bold", loc="left")
ax_bar.legend(fontsize=7, loc="upper right")

# n数
for k, t in enumerate(tier_order):
    n = vdr_df[vdr_df["tier"] == t].shape[0]
    ax_bar.text(k, 103, f"n={n}", ha="center", fontsize=7.5, color="#555")

# ── Panel D: IL17A 特殊ケース (Tier2の疾患依存性) ────────────────────────
ax_il17.axis("off")
ax_il17.set_title("D  Tier2 (downstream) — Disease context matters",
                  fontsize=9, fontweight="bold", loc="left")

il17_text = (
    "IL17A = Downstream of IL23A (VDR-direct)\n"
    "  IL23A(VDR=211) → Th17 → IL17A\n\n"
    "Psoriasis (skin):  ◎ 長期維持\n"
    "  → 皮膚IL17Aは主エフェクター\n"
    "  → 遮断で皮疹消失・長期維持\n\n"
    "Crohn's Disease (gut):  ❌ 禁忌\n"
    "  → 腸管IL17Aは粘膜バリア維持に必須\n"
    "  → 遮断でCD増悪(IBD発症報告あり)\n\n"
    "→ VDR回路下流を遮断すると\n"
    "  残った機能まで失われる\n"
    "  (根本=VDR欠乏)は未修正のまま\n\n"
    "比較: 上流(IL23A)遮断は\n"
    "  CD・Ps両方で長期◎ [T1]"
)
ax_il17.text(0.05, 0.97, il17_text, transform=ax_il17.transAxes,
             fontsize=7.8, va="top", family="monospace",
             bbox=dict(fc="#E8F5E9", ec="#388E3C", boxstyle="round,pad=0.5"))

# ─── Save ─────────────────────────────────────────────────────────────────
for ext in ["png", "pdf"]:
    fig.savefig(OUT / f"fig_maintenance_logic.{ext}",
                dpi=200, bbox_inches="tight")

print("\nSaved: fig_maintenance_logic.png / .pdf")
print(f"       maintenance_logic_classification.csv  (n={len(df)})")

# ─── 最終サマリー出力 ─────────────────────────────────────────────────────
print("\n" + "=" * 60)
print("KEY FINDING: Logic tier vs Long-term outcome (VDR軸疾患)")
print("=" * 60)
for t in [1,2,3,4]:
    sub = vdr_df[vdr_df["tier"] == t]
    if sub.empty: continue
    good = sub[sub["lt_status"].isin(["◎"])].shape[0]
    warn = sub[sub["lt_status"].isin(["⚠️","❌"])].shape[0]
    tot  = sub.shape[0]
    print(f"  T{t} {TIERS[t]['pred']:15s}  ◎={good}/{tot}"
          f"  ⚠️+❌={warn}/{tot}")

print("\n→ T1(VDR直接): ◎が最多、FDA Black Box なし")
print("→ T3(GR回路):  JAK阻害薬は全薬FDA Black Box Warning")
print("→ T2(下流):    疾患依存（Ps◎ vs CD❌ for IL17A）")
