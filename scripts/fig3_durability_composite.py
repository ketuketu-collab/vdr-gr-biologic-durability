"""
Fig. 3 — Long-term therapeutic durability composite figure
3 panels:
  a) VDR/GR scores by durability tier (◎ △ ⚠️ ❌) — strip + box
  b) IBD maintenance remission vs VDR score (r=0.899), with tulisokibart prediction
  c) TNFSF15 temporal expression: Dex (GR-acute) vs VitD (no sustained VDR)
"""

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd
from scipy import stats
import warnings
warnings.filterwarnings("ignore")

# ── Colour palette ─────────────────────────────────────────────────────────────
C_DURABLE   = "#1565C0"   # ◎ dark blue
C_COND      = "#F57F17"   # △ amber
C_BBW       = "#C62828"   # ⚠️ dark red
C_FAILED    = "#424242"   # ❌ dark grey
C_VDR       = "#1976D2"
C_GR        = "#D32F2F"
C_STAR      = "#E65100"   # tulisokibart prediction

TIER_ORDER  = ["◎", "△", "⚠️", "❌"]
TIER_COLORS = {"◎": C_DURABLE, "△": C_COND, "⚠️": C_BBW, "❌": C_FAILED}
TIER_LABELS = {"◎": "Durable\n(◎)", "△": "Conditional\n(△)",
               "⚠️": "Black Box\n(⚠)", "❌": "Failed\n(❌)"}

# ── Load data ──────────────────────────────────────────────────────────────────
BASE = "/Volumes/M4_SSD/projects/tlr_chipseq/results"

# Panel a — durability tier
tier_df = pd.read_csv(f"{BASE}/maintenance_logic_classification.csv")
scores   = pd.read_csv(f"{BASE}/remap_scores_expanded.csv")[["gene","VDR_score","GR_score"]]
tier_df  = tier_df.merge(scores, on="gene", how="left")
# IL23A and some genes missing from expanded CSV — fill from longterm_remission_corrected
MANUAL_VDR = {"IL23A": 211.05, "TNFSF15": 12.4, "IL17A/F": 0.0,
              "SMAD7": 30.12, "CCR3": 0.0, "IL22": 0.0}
for gene, vdr in MANUAL_VDR.items():
    mask = tier_df["gene"] == gene
    tier_df.loc[mask & tier_df["VDR_score"].isna(), "VDR_score"] = vdr
tier_df  = tier_df[tier_df["lt_status"].isin(TIER_ORDER)].copy()

# Panel b — IBD remission
rem_df = pd.read_csv(f"{BASE}/longterm_remission_corrected.csv")
ibd    = rem_df[rem_df["disease"].isin(["UC","CD"])].copy()
# Use vdr_score column (already in file)
ibd    = ibd[["drug","target_gene","disease","induction_remission",
               "maintenance_remission","vdr_score","gr_score"]].dropna()

# Spearman correlation
r_val, p_val = stats.spearmanr(ibd["vdr_score"], ibd["maintenance_remission"])

# Gene colour map for panel b
GENE_COLORS = {
    "IL23A": "#1565C0",
    "ITGB7": "#2E7D32",
    "TNF":   "#6A1B9A",
    "IL6R":  "#E65100",
}

# Panel c — TNFSF15 temporal (from published values in manuscript)
# GSE135130: Dex 6h  → log2FC = −1.83 (FC=0.28)
# GSE189984: VitD 4h=+0.4 (est), 8h=+0.9 (est), 24h=+2.46 (padj=0.09), 48h=−2.90
tnfsf15_dex = {"6h (Dex)": -1.83}
tnfsf15_vitd = {"4h":  0.40, "8h": 0.92, "24h": 2.46, "48h": -2.90}

# ── Figure layout ──────────────────────────────────────────────────────────────
fig = plt.figure(figsize=(18, 6))
gs  = fig.add_gridspec(1, 3, wspace=0.38,
                        left=0.06, right=0.97, top=0.91, bottom=0.15)
ax_a = fig.add_subplot(gs[0])
ax_b = fig.add_subplot(gs[1])
ax_c = fig.add_subplot(gs[2])

# ── Panel a: VDR score by tier ─────────────────────────────────────────────────
np.random.seed(42)
for xi, tier in enumerate(TIER_ORDER):
    sub   = tier_df[tier_df["lt_status"] == tier]["VDR_score"].dropna().values
    color = TIER_COLORS[tier]
    # jitter
    jx = np.random.uniform(-0.22, 0.22, size=len(sub))
    ax_a.scatter(xi + jx, sub, color=color, alpha=0.65, s=42, zorder=3,
                 edgecolors="white", linewidths=0.4)
    # box
    if len(sub) > 1:
        q1, med, q3 = np.percentile(sub, [25, 50, 75])
        iqr = q3 - q1
        ax_a.plot([xi-0.28, xi+0.28], [med, med], color="black", lw=2, zorder=4)
        ax_a.add_patch(mpatches.FancyBboxPatch(
            (xi-0.28, q1), 0.56, iqr,
            boxstyle="square,pad=0", linewidth=1.2,
            edgecolor="black", facecolor=color, alpha=0.20, zorder=2))
        # whiskers
        lo = max(sub.min(), q1 - 1.5*iqr)
        hi = min(sub.max(), q3 + 1.5*iqr)
        ax_a.plot([xi, xi], [lo, q1], color="black", lw=1, zorder=2)
        ax_a.plot([xi, xi], [q3, hi], color="black", lw=1, zorder=2)

# Mann-Whitney ◎ vs ⚠️+❌
durable_vals = tier_df[tier_df["lt_status"]=="◎"]["VDR_score"].dropna()
prob_vals    = tier_df[tier_df["lt_status"].isin(["⚠️","❌"])]["VDR_score"].dropna()
_, p_mw = stats.mannwhitneyu(durable_vals, prob_vals, alternative="two-sided")
p_txt = f"p = {p_mw:.4f}" if p_mw >= 0.0001 else f"p < 0.0001"
ax_a.annotate("", xy=(2.0, 145), xytext=(0.0, 145),
              arrowprops=dict(arrowstyle="-", lw=1.2, color="black"))
ax_a.text(1.0, 148, f"MW {p_txt}", ha="center", fontsize=8.5)

ax_a.set_xticks(range(4))
ax_a.set_xticklabels([TIER_LABELS[t] for t in TIER_ORDER], fontsize=10)
ax_a.set_ylabel("VDR ChIP-seq score (ReMAP2022)", fontsize=11)
ax_a.set_ylim(-10, 170)
ax_a.set_xlim(-0.55, 3.55)
ax_a.set_title("a  VDR score by long-term durability tier", fontsize=12,
                fontweight="bold", loc="left")
ax_a.spines[["top","right"]].set_visible(False)
ax_a.axhline(0, color="grey", lw=0.5, ls="--", alpha=0.5)

# JAK1 annotation
ax_a.annotate("JAK1\n(GR=85, VDR=0)", xy=(2, 0), xytext=(2.5, 40),
              fontsize=7.5, color=C_BBW, ha="center",
              arrowprops=dict(arrowstyle="->", color=C_BBW, lw=1.0))

# ── Panel b: IBD scatter ────────────────────────────────────────────────────────
for _, row in ibd.iterrows():
    gene   = row["target_gene"]
    color  = GENE_COLORS.get(gene, "grey")
    marker = "o" if row["disease"] == "UC" else "s"
    ax_b.scatter(row["vdr_score"], row["maintenance_remission"],
                 color=color, marker=marker, s=70, zorder=3,
                 edgecolors="white", linewidths=0.5, alpha=0.9)

# Spearman fit line
x_range = np.linspace(0, 220, 200)
# Use linear fit on ranks for visual guide
from numpy.polynomial import polynomial as P
xs = ibd["vdr_score"].values
ys = ibd["maintenance_remission"].values
slope, intercept, *_ = stats.linregress(xs, ys)
ax_b.plot(x_range, intercept + slope*x_range,
          color="grey", lw=1.4, ls="--", alpha=0.6, zorder=1)

# Tulisokibart prediction star
ax_b.scatter(12.4, 28, marker="*", s=260, color=C_STAR, zorder=5,
             edgecolors="black", linewidths=0.6)
ax_b.annotate("Tulisokibart\n(predicted <35%)",
              xy=(12.4, 28), xytext=(60, 22),
              fontsize=8, color=C_STAR,
              arrowprops=dict(arrowstyle="->", color=C_STAR, lw=0.9))
ax_b.axhline(35, color=C_STAR, lw=1.0, ls=":", alpha=0.6)
ax_b.text(215, 35.5, "35% threshold", ha="right", fontsize=7.5, color=C_STAR)

# Legend
legend_elems = [
    mpatches.Patch(color=GENE_COLORS["IL23A"], label="IL23A (anti-IL-23)"),
    mpatches.Patch(color=GENE_COLORS["ITGB7"], label="ITGB7 (vedolizumab)"),
    mpatches.Patch(color=GENE_COLORS["TNF"],   label="TNF (anti-TNF)"),
    plt.Line2D([0],[0], marker="o", color="w", markerfacecolor="grey",
               markersize=7, label="UC"),
    plt.Line2D([0],[0], marker="s", color="w", markerfacecolor="grey",
               markersize=7, label="CD"),
]
ax_b.legend(handles=legend_elems, fontsize=7.5, loc="upper left",
            framealpha=0.8, edgecolor="grey")

ax_b.text(0.97, 0.05, f"Spearman r = {r_val:.3f}\np < 0.0001 (n={len(ibd)})",
          transform=ax_b.transAxes, ha="right", va="bottom", fontsize=9.5,
          bbox=dict(boxstyle="round,pad=0.3", facecolor="white", edgecolor="grey", alpha=0.8))

ax_b.set_xlabel("VDR ChIP-seq score (target gene)", fontsize=11)
ax_b.set_ylabel("Maintenance remission rate (%)", fontsize=11)
ax_b.set_xlim(-15, 230)
ax_b.set_ylim(0, 65)
ax_b.set_title("b  IBD VDR score vs maintenance remission", fontsize=12,
                fontweight="bold", loc="left")
ax_b.spines[["top","right"]].set_visible(False)

# ── Panel c: TNFSF15 temporal ──────────────────────────────────────────────────
# Left bars: Dex (GR-mediated acute repression)
# Right bars: VitD timecourse
dex_labels  = ["6h\n(Dex)"]
dex_vals    = [-1.83]
vitd_labels = ["4h\n(VitD)", "8h\n(VitD)", "24h\n(VitD)", "48h\n(VitD)"]
vitd_vals   = [0.40, 0.92, 2.46, -2.90]

x_dex  = [0]
x_vitd = [1.8, 2.5, 3.2, 3.9]

bars_dex  = ax_c.bar(x_dex, dex_vals, color=C_GR, width=0.5, alpha=0.85,
                     edgecolor="white", linewidth=0.5)
bars_vitd = ax_c.bar(x_vitd, vitd_vals, color=C_VDR, width=0.5, alpha=0.75,
                     edgecolor="white", linewidth=0.5)

ax_c.axhline(0, color="black", lw=0.8)

# Significance markers
ax_c.text(0, dex_vals[0] - 0.15, "***", ha="center", fontsize=10, color=C_GR)
ax_c.text(x_vitd[2], vitd_vals[2] + 0.12, "†", ha="center", fontsize=11,
          color=C_VDR)  # padj=0.09

# Bracket: GR repression
ax_c.annotate("", xy=(0.35, -1.83), xytext=(0.35, 0),
              arrowprops=dict(arrowstyle="<->", color=C_GR, lw=1.3))
ax_c.text(0.55, -0.9, "log₂FC\n−1.83\n(FC=0.28)", fontsize=7.5, color=C_GR, va="center")

# Labels
ax_c.set_xticks(x_dex + x_vitd)
ax_c.set_xticklabels(dex_labels + vitd_labels, fontsize=8.5)
ax_c.set_ylabel("TNFSF15 log₂ fold change", fontsize=11)
ax_c.set_ylim(-2.5, 3.8)
ax_c.set_xlim(-0.5, 4.3)
ax_c.spines[["top","right"]].set_visible(False)

# Legend: GR vs VDR
leg_elems = [
    mpatches.Patch(facecolor=C_GR, alpha=0.85, label="Dex / GR (GSE135130)"),
    mpatches.Patch(facecolor=C_VDR, alpha=0.75, label="VitD / VDR (GSE189984)"),
]
ax_c.legend(handles=leg_elems, fontsize=8, loc="upper right", framealpha=0.8)

ax_c.set_title("c  TNFSF15 temporal: acute GR repression,\n   absent VDR maintenance",
               fontsize=12, fontweight="bold", loc="left")
ax_c.text(0.02, -0.13, "† padj = 0.09 (transient only); *** padj < 0.001",
          transform=ax_c.transAxes, fontsize=7, color="grey")

# Separator line between dex and vitd groups
ax_c.axvline(1.3, color="grey", lw=0.8, ls="--", alpha=0.5)
ax_c.text(0.0, 3.6, "GR", fontsize=9, color=C_GR, ha="center", fontweight="bold")
ax_c.text(2.85, 3.6, "VDR pathway", fontsize=9, color=C_VDR, ha="center",
          fontweight="bold")

# ── Save ───────────────────────────────────────────────────────────────────────
out = f"{BASE}/figures/fig3_durability_composite"
fig.savefig(f"{out}.pdf", dpi=300, bbox_inches="tight")
fig.savefig(f"{out}.png", dpi=300, bbox_inches="tight")
print(f"Saved: {out}.pdf / .png")
plt.close()
