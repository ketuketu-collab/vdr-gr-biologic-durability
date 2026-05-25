"""
Fig. 5 — Mechanistic model and prospective prediction
3 panels:
  a) GR → CYP24A1 → VitD depletion → VDR circuit failure schematic
  b) Pipeline prediction scatter (VDR/GR ratio; tulisokibart highlighted)
  c) Induction-to-maintenance trajectory: VDR-dominant vs GR-dominant ± VitD
"""

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.patheffects as pe
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch
import numpy as np
import pandas as pd
import warnings
warnings.filterwarnings("ignore")

C_VDR   = "#1565C0"
C_GR    = "#C62828"
C_STAR  = "#E65100"
C_GREEN = "#2E7D32"
C_AMBER = "#F57F17"
C_GREY  = "#757575"

BASE = "/Volumes/M4_SSD/projects/tlr_chipseq/results"

# ── Load pipeline data (tulisokibart already fixed to TNFSF15 in CSV) ─────────
pipe = pd.read_csv(f"{BASE}/pipeline_prediction.csv")

# Reference approved drugs for panel b
reference = pd.DataFrame([
    dict(drug="anti-IL-23", gene="IL23A", vdr=211.05, gr=70.34,  group="ref_approved"),
    dict(drug="vedolizumab", gene="ITGB7", vdr=62.27, gr=17.23,  group="ref_approved"),
    dict(drug="anti-TNF",   gene="TNF",   vdr=36.86,  gr=19.49,  group="ref_approved"),
    dict(drug="JAK-i (⚠ BBW)", gene="JAK1", vdr=12.81, gr=85.48, group="ref_bbw"),
])

# ── Figure layout ──────────────────────────────────────────────────────────────
fig = plt.figure(figsize=(18, 6.5))
gs  = fig.add_gridspec(1, 3, wspace=0.40,
                        left=0.05, right=0.97, top=0.90, bottom=0.13)
ax_a = fig.add_subplot(gs[0])
ax_b = fig.add_subplot(gs[1])
ax_c = fig.add_subplot(gs[2])

# ══════════════════════════════════════════════════════════════════════════════
# Panel a — Mechanistic schematic: GR→CYP24A1→VitD depletion→VDR circuit failure
# ══════════════════════════════════════════════════════════════════════════════
ax_a.set_xlim(0, 10)
ax_a.set_ylim(0, 10)
ax_a.axis("off")
ax_a.set_title("a  GR→CYP24A1→VitD depletion mechanism", fontsize=12,
                fontweight="bold", loc="left")

def box(ax, x, y, w, h, label, sublabel=None, color="#1565C0", fontsize=10):
    ax.add_patch(FancyBboxPatch((x, y), w, h,
                  boxstyle="round,pad=0.15", facecolor=color, alpha=0.15,
                  edgecolor=color, linewidth=2))
    ax.text(x + w/2, y + h/2 + (0.15 if sublabel else 0),
            label, ha="center", va="center", fontsize=fontsize,
            fontweight="bold", color=color)
    if sublabel:
        ax.text(x + w/2, y + h/2 - 0.3, sublabel,
                ha="center", va="center", fontsize=7.5, color=C_GREY)

def arrow(ax, x1, y1, x2, y2, color="black", label="", lw=2, style="->"):
    ax.annotate("", xy=(x2, y2), xytext=(x1, y1),
                arrowprops=dict(arrowstyle=style, color=color, lw=lw,
                                connectionstyle="arc3,rad=0.0"))
    if label:
        mx, my = (x1+x2)/2, (y1+y2)/2
        ax.text(mx+0.12, my, label, fontsize=8, color=color,
                ha="left", va="center", style="italic")

# Row 1: Glucocorticoid stimulus → GR activation
box(ax_a, 0.3, 8.2, 2.8, 1.0, "Steroid / GC\nstimulus", color=C_GR, fontsize=9)
arrow(ax_a, 3.1, 8.7, 4.2, 8.7, color=C_GR, label="")
box(ax_a, 4.2, 8.2, 2.8, 1.0, "GR activation\n(NR3C1)", color=C_GR, fontsize=9)

# GR → CYP24A1 (downward arrow)
arrow(ax_a, 5.6, 8.2, 5.6, 7.1, color=C_GR, label="Direct\nbinding\n(fold=14.7×)", lw=2)
box(ax_a, 4.2, 6.0, 2.8, 1.0, "CYP24A1↑↑\n(VitD catabolism)", color=C_GR, fontsize=9)

# CYP24A1 → VitD depletion (left)
arrow(ax_a, 4.2, 6.5, 3.1, 6.5, color=C_GR, label="")
box(ax_a, 0.3, 6.0, 2.8, 1.0, "1,25(OH)₂D₃\n↓ depleted", color=C_AMBER, fontsize=9)

# VitD depleted → VDR ligand pool gone (downward)
arrow(ax_a, 1.7, 6.0, 1.7, 4.9, color=C_AMBER, label="")
box(ax_a, 0.3, 3.8, 2.8, 1.0, "VDR ligand\npool ↓↓", color=C_AMBER, fontsize=9)

# VDR ligand pool → VDR circuit failure (right)
arrow(ax_a, 3.1, 4.3, 4.2, 4.3, color=C_AMBER, label="")
box(ax_a, 4.2, 3.8, 2.8, 1.0, "VDR circuit\nfailure", color=C_VDR, fontsize=9)

# VDR circuit failure → chronic remission loss
arrow(ax_a, 5.6, 3.8, 5.6, 2.7, color=C_VDR, label="")
box(ax_a, 4.2, 1.6, 2.8, 1.0, "Maintenance\nremission ↓", color="#424242", fontsize=9)

# Evidence labels
ax_a.text(7.3, 8.7, "① ORAL Surveillance\n   tofacitinib ↑CV/cancer",
          fontsize=7.5, color=C_GREY, va="center")
ax_a.text(7.3, 6.5, "② ENCODE ChIP-seq\n   BEAS-2B + Dex",
          fontsize=7.5, color=C_GREY, va="center")
ax_a.text(7.3, 4.3, "③ NHANES n=31,799\n   VitD↓ → ER visit↑ p=0.005",
          fontsize=7.5, color=C_GREY, va="center")

# "Proposed intervention" — VitD supplement arrow
ax_a.add_patch(FancyArrowPatch((3.5, 4.3), (3.5, 5.0),
               arrowstyle="->", color=C_GREEN, lw=2,
               mutation_scale=15,
               connectionstyle="arc3,rad=0"))
ax_a.text(3.5, 5.2, "VitD Rx\n(proposed)", ha="center", fontsize=8,
          color=C_GREEN, fontweight="bold")

# ══════════════════════════════════════════════════════════════════════════════
# Panel b — Predictive validation: pipeline biologics (VDR vs GR scatter)
# ══════════════════════════════════════════════════════════════════════════════
ax_b.set_title("b  Predictive validation of pipeline biologics\n"
               "   using VDR/GR ChIP-seq rule",
               fontsize=11, fontweight="bold", loc="left")

# Shaded zones
lim_b = 320
xs = np.linspace(0, lim_b, 300)
ax_b.fill_between(xs, 0, xs,    alpha=0.05, color=C_VDR, zorder=0)
ax_b.fill_between(xs, xs, lim_b, alpha=0.05, color=C_GR,  zorder=0)
ax_b.plot(xs, xs, "--", color="grey", lw=0.8, alpha=0.45, zorder=1)
ax_b.text(155, 140, "VDR=GR", fontsize=7, color="grey", rotation=38, alpha=0.6)
ax_b.text(145, 25,  "VDR-dominant\n(durable)", fontsize=8, color=C_VDR, alpha=0.75)
ax_b.text(15,  250, "GR-dominant\n(failure risk)", fontsize=8, color=C_GR,  alpha=0.75)

# Reference approved drugs
ref_labels = {
    "anti-IL-23": (10, -14),
    "vedolizumab": (6,   4),
    "anti-TNF":    (6,   4),
    "JAK-i (⚠ BBW)": (6, 4),
}
for _, row in reference.iterrows():
    col = C_VDR if row["group"] == "ref_approved" else C_GR
    ax_b.scatter(row["vdr"], row["gr"], marker="D", s=110, color=col,
                 edgecolors="black", linewidths=0.6, zorder=4, alpha=0.9)
    ox, oy = ref_labels.get(row["drug"], (6, 4))
    ax_b.annotate(row["drug"], xy=(row["vdr"], row["gr"]),
                  xytext=(ox, oy), textcoords="offset points",
                  fontsize=7.5, color=col, fontweight="bold")

# Pipeline agents (excluding tulisokibart)
for _, row in pipe.iterrows():
    if row["drug"] == "tulisokibart":
        continue
    if row["prediction"] == "predicted_failure":
        col, mk, ms = C_GR,  "X", 80
    elif row["prediction"] == "approved_recent":
        col, mk, ms = "#2E7D32", "o", 65
    else:
        col, mk, ms = C_VDR, "o", 65
    ax_b.scatter(row["vdr_score"], row["gr_score"],
                 marker=mk, s=ms, color=col, alpha=0.75, zorder=3,
                 edgecolors="white", linewidths=0.4)
    # Label agents above GR>50
    if row["gr_score"] > 50 or row["vdr_score"] > 20:
        ax_b.annotate(f"{row['drug']}\n({row['gene']})",
                      xy=(row["vdr_score"], row["gr_score"]),
                      xytext=(5, 3), textcoords="offset points",
                      fontsize=6.5, color=col)

# ── Tulisokibart ★ ──
ax_b.scatter(12.4, 286.0, marker="*", s=450, color=C_STAR, zorder=6,
             edgecolors="black", linewidths=0.8)
ax_b.annotate(
    "Tulisokibart (TNFSF15)\n"
    "VDR=12.4  GR=286.0\n"
    "Ratio = 0.043\n"
    "★ Pre-registered: FAILURE\n"
    "   osf.io/tnp63  (2026-05-25)\n"
    "   Validation: ATLAS-UC H1 2027",
    xy=(12.4, 286.0), xytext=(65, 195),
    fontsize=8, color=C_STAR, fontweight="bold",
    arrowprops=dict(arrowstyle="->", color=C_STAR, lw=1.5,
                    connectionstyle="arc3,rad=-0.15"),
    bbox=dict(boxstyle="round,pad=0.35", facecolor="white",
              edgecolor=C_STAR, alpha=0.92))

# Legend
leg_b = [
    plt.Line2D([0],[0], marker="D", color="w", markerfacecolor=C_VDR,
               markersize=8, label="Approved ◎ (reference)"),
    plt.Line2D([0],[0], marker="D", color="w", markerfacecolor=C_GR,
               markersize=8, label="Approved ⚠ Black Box"),
    plt.Line2D([0],[0], marker="o", color="w", markerfacecolor="#2E7D32",
               markersize=7, label="Recently approved"),
    plt.Line2D([0],[0], marker="X", color="w", markerfacecolor=C_GR,
               markersize=8, label="Pipeline: predicted failure"),
    plt.Line2D([0],[0], marker="*", color="w", markerfacecolor=C_STAR,
               markersize=11, label="Tulisokibart (pre-registered ★)"),
]
ax_b.legend(handles=leg_b, fontsize=7.5, loc="lower right",
            framealpha=0.88, edgecolor="grey")

ax_b.set_xlabel("VDR (NR1I1) ChIP-seq score  [ReMAP2022]", fontsize=10)
ax_b.set_ylabel("GR (NR3C1) ChIP-seq score  [ReMAP2022]", fontsize=10)
ax_b.set_xlim(-15, 250)
ax_b.set_ylim(-10, 320)
ax_b.spines[["top","right"]].set_visible(False)

# ══════════════════════════════════════════════════════════════════════════════
# Panel c — Sequential GR→VDR therapeutic paradigm
# ══════════════════════════════════════════════════════════════════════════════
ax_c.set_title("c  Sequential GR→VDR therapeutic paradigm",
               fontsize=12, fontweight="bold", loc="left")

weeks     = [0, 8, 12, 24, 44, 52, 60, 104]

# ① GR-dominant chronic use: decent induction, rapid decline (current failure mode)
gr_chronic   = [0, 28, 30, 28, 22, 18, 14, 10]

# ② VDR-dominant from start (e.g., anti-IL-23): modest induction, sustained maintenance
vdr_only     = [0, 22, 25, 35, 43, 44, 44, 45]

# ③ Sequential: GR induction (week 0–12) → switch to VDR maintenance (proposed optimal)
#    Induction matches GR arm; after switch rises to VDR-dominant level + VitD benefit
seq_switch   = [0, 28, 30, 36, 46, 50, 52, 53]

# ④ Sequential + VitD at transition (fully optimised)
seq_vitd     = [0, 28, 31, 40, 50, 54, 56, 57]

ax_c.plot(weeks, gr_chronic,  "--x",  color=C_GR,    lw=2.0, ms=7,
          label="GR-dominant chronic\n  (current failure mode)",
          zorder=3, markeredgewidth=2)
ax_c.plot(weeks, vdr_only,    "-o",   color=C_VDR,   lw=2.0, ms=6,
          label="VDR-dominant from induction\n  (slow onset, durable)",
          zorder=3)
ax_c.plot(weeks, seq_switch,  "-s",   color="#6A1B9A", lw=2.5, ms=6,
          label="Sequential: GR induction\n  → switch to VDR at wk 12 ★",
          zorder=4)
ax_c.plot(weeks, seq_vitd,    "-^",   color=C_GREEN,  lw=2.0, ms=6,
          ls=(0,(5,2)),
          label="Sequential + VitD at switch\n  (proposed optimal)",
          zorder=4)

# TNF inhibitor reference band
ax_c.axhline(35, color=C_AMBER, lw=1.0, ls=":", alpha=0.7)
ax_c.text(106, 35.5, "TNF-i\nbenchmark\n(35%)", ha="right", fontsize=7.5,
          color=C_AMBER)

# Induction vs maintenance phase shading
ax_c.axvspan(0,  12,  alpha=0.07, color=C_GR,  zorder=0)
ax_c.axvspan(12, 104, alpha=0.04, color=C_VDR, zorder=0)
ax_c.text(6,  62, "Induction\n(GR phase)",   ha="center", fontsize=8, color=C_GR,
          fontweight="bold", va="top")
ax_c.text(58, 62, "Maintenance\n(VDR phase)", ha="center", fontsize=8, color=C_VDR,
          fontweight="bold", va="top")
ax_c.axvline(12, color="grey", lw=1.2, ls="--", alpha=0.6)

# Switch point annotation
ax_c.annotate("Switch to\nVDR-dominant\n+ VitD Rx",
              xy=(12, 31), xytext=(18, 18),
              fontsize=8, color="#6A1B9A", fontweight="bold",
              arrowprops=dict(arrowstyle="->", color="#6A1B9A", lw=1.2,
                              connectionstyle="arc3,rad=0.2"))

# GR chronic failure annotation
ax_c.annotate("GR→CYP24A1→VitD↓\n→ VDR circuit fails",
              xy=(44, 20), xytext=(55, 9),
              fontsize=7.5, color=C_GR,
              arrowprops=dict(arrowstyle="->", color=C_GR, lw=0.9))

ax_c.set_xlabel("Weeks from induction", fontsize=11)
ax_c.set_ylabel("Estimated remission rate (%)", fontsize=11)
ax_c.set_xlim(-5, 110)
ax_c.set_ylim(0, 68)
ax_c.set_xticks([0, 12, 24, 44, 52, 104])
ax_c.legend(fontsize=7.5, loc="lower right", framealpha=0.88, edgecolor="grey")
ax_c.spines[["top","right"]].set_visible(False)
ax_c.text(0.02, -0.14,
          "Trajectories estimated from published Phase III data; sequential and VitD arms are\n"
          "hypothetical proposals based on the GR→CYP24A1→VitD depletion mechanism.",
          transform=ax_c.transAxes, fontsize=6.5, color="grey")

# ── Save ───────────────────────────────────────────────────────────────────────
out = f"{BASE}/figures/fig5_mechanistic_pipeline"
fig.savefig(f"{out}.pdf", dpi=300, bbox_inches="tight")
fig.savefig(f"{out}.png", dpi=300, bbox_inches="tight")
print(f"Saved: {out}.pdf / .png")
plt.close()
