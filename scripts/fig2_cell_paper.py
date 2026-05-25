#!/usr/bin/env python3
"""
Figure 2 for Cell paper:
VDR/GR dual nuclear receptor axis establishes complementary TLR regulation
Panel A: Human THP-1 ChIP-seq — VDR vs GR at TLR loci
Panel B: Cross-species conservation (human vs mouse)
Panel C: 2-brake model schematic
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.patches as FancyArrowPatch
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch
from matplotlib.lines import Line2D
import warnings
warnings.filterwarnings("ignore")

OUT = "/Volumes/M4_SSD/projects/tlr_chipseq/results/figures/fig2_cell_2brake.pdf"
OUT_PNG = "/Volumes/M4_SSD/projects/tlr_chipseq/results/figures/fig2_cell_2brake.png"

# ── Color palette ──────────────────────────────────────────────
VDR_COL  = "#2166AC"   # deep blue
GR_COL   = "#D6604D"   # warm red
VDR_LITE = "#92C5DE"
GR_LITE  = "#F4A582"
GREY     = "#BDBDBD"
BG_ZONE  = "#EFF7FF"   # GR-free zone highlight

# ── Data ───────────────────────────────────────────────────────
# Human THP-1 (immune-cell filtered, best score across windows)
human_loci = ["TLR10", "TLR1", "TLR6", "TLR2"]
human_vdr  = [3.98,    3.98,   8.46,   11.37]
human_gr   = [0.0,     0.0,    12.25,  56.14]

# Mouse (macrophage GSE124722 / BMDM)
mouse_loci = ["Tlr1", "Tlr6", "Tlr2"]
mouse_vdr  = [9.05,   14.64,  10.49]   # macrophage
mouse_gr   = [0.0,    4.48,   16.74]   # BMDM

# ── Figure layout ──────────────────────────────────────────────
fig = plt.figure(figsize=(18, 6.5))
fig.patch.set_facecolor("white")

gs = fig.add_gridspec(1, 3, wspace=0.38,
                      left=0.06, right=0.97, top=0.88, bottom=0.16)
ax_a = fig.add_subplot(gs[0])
ax_b = fig.add_subplot(gs[1])
ax_c = fig.add_subplot(gs[2])

# ══════════════════════════════════════════════════════════════
# PANEL A  — Human THP-1: VDR vs GR across TLR loci
# ══════════════════════════════════════════════════════════════
x = np.arange(len(human_loci))
w = 0.32

# GR-free zone background
ax_a.axvspan(-0.5, 1.5, color=BG_ZONE, alpha=0.6, zorder=0, label="_nolegend_")
ax_a.text(0.5, 62, "GR-free zone", ha="center", va="bottom",
          fontsize=8.5, color=VDR_COL, fontstyle="italic", fontweight="bold")

bars_vdr = ax_a.bar(x - w/2, human_vdr, w, color=VDR_COL, label="VDR",
                    edgecolor="white", linewidth=0.5, zorder=3)
bars_gr  = ax_a.bar(x + w/2, human_gr,  w, color=GR_COL,  label="GR/NR3C1",
                    edgecolor="white", linewidth=0.5, zorder=3)

# Value labels
for bar, val in zip(bars_vdr, human_vdr):
    if val > 0:
        ax_a.text(bar.get_x() + bar.get_width()/2, val + 0.8,
                  f"{val:.1f}", ha="center", va="bottom", fontsize=7.5,
                  color=VDR_COL, fontweight="bold")
for bar, val in zip(bars_gr, human_gr):
    if val > 0:
        ax_a.text(bar.get_x() + bar.get_width()/2, val + 0.8,
                  f"{val:.1f}", ha="center", va="bottom", fontsize=7.5,
                  color=GR_COL, fontweight="bold")
    else:
        ax_a.text(bar.get_x() + bar.get_width()/2, 0.8,
                  "N.D.", ha="center", va="bottom", fontsize=6.5,
                  color=GREY)

ax_a.set_xticks(x)
ax_a.set_xticklabels(human_loci, fontsize=11, fontweight="bold")
ax_a.set_ylabel("ChIP-seq binding score\n(ReMap2022, THP-1)", fontsize=10)
ax_a.set_ylim(0, 70)
ax_a.spines[["top", "right"]].set_visible(False)
ax_a.legend(fontsize=9, frameon=False, loc="upper left")
ax_a.set_title("Human THP-1 (monocyte)", fontsize=11, pad=8, fontweight="bold")

# Bracket for GR-free zone
ax_a.annotate("", xy=(-0.48, 66), xytext=(1.48, 66),
              arrowprops=dict(arrowstyle="|-|", color=VDR_COL, lw=1.5))

ax_a.text(-0.02, 72, "A", transform=ax_a.transAxes,
          fontsize=15, fontweight="bold", va="top")

# ══════════════════════════════════════════════════════════════
# PANEL B  — Cross-species conservation heatmap
# ══════════════════════════════════════════════════════════════
# Normalize: max across all values for color scale
all_vals = np.array([
    [3.98,  0.0,    9.05,  0.0  ],   # TLR1/Tlr1
    [8.46,  12.25,  14.64, 4.48 ],   # TLR6/Tlr6
    [11.37, 56.14,  10.49, 16.74],   # TLR2/Tlr2
])
labels_row = ["TLR1 / Tlr1", "TLR6 / Tlr6", "TLR2 / Tlr2"]
labels_col = ["VDR\n(Human)", "GR\n(Human)", "VDR\n(Mouse)", "GR\n(Mouse)"]

# Use separate colormaps: blue for VDR, red for GR
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable

# Max for each TF across both species
vdr_max = max(3.98, 9.05, 8.46, 14.64, 11.37, 10.49)  # ~14.64
gr_max  = max(0, 0, 12.25, 4.48, 56.14, 16.74)         # 56.14

norm_vdr = Normalize(0, vdr_max)
norm_gr  = Normalize(0, gr_max)

import matplotlib.cm as cm
cmap_vdr = cm.Blues
cmap_gr  = cm.Reds

cell_w, cell_h = 0.18, 0.22
x0, y0 = 0.12, 0.68

for ri, (row, lbl) in enumerate(zip(all_vals, labels_row)):
    y_pos = y0 - ri * (cell_h + 0.04)
    # Row label
    ax_b.text(0.02, y_pos + cell_h/2, lbl, ha="left", va="center",
              fontsize=9.5, transform=ax_b.transAxes, fontweight="bold")
    for ci, val in enumerate(row):
        x_pos = x0 + ci * (cell_w + 0.025)
        is_vdr = (ci % 2 == 0)
        if is_vdr:
            color = cmap_vdr(norm_vdr(val))
        else:
            color = cmap_gr(norm_gr(val)) if val > 0 else "#F5F5F5"

        rect = FancyBboxPatch((x_pos, y_pos), cell_w, cell_h,
                               boxstyle="round,pad=0.01",
                               facecolor=color, edgecolor="white",
                               linewidth=1.5, transform=ax_b.transAxes,
                               zorder=3)
        ax_b.add_patch(rect)

        # Text
        txt_color = "white" if (val > 20 or (is_vdr and val > 10)) else "#333333"
        label_txt = f"{val:.1f}" if val > 0 else "N.D."
        ax_b.text(x_pos + cell_w/2, y_pos + cell_h/2, label_txt,
                  ha="center", va="center", fontsize=9, fontweight="bold",
                  color=txt_color, transform=ax_b.transAxes)

# Column headers
col_colors = [VDR_COL, GR_COL, VDR_COL, GR_COL]
for ci, (lbl, cc) in enumerate(zip(labels_col, col_colors)):
    x_pos = x0 + ci * (cell_w + 0.025) + cell_w/2
    ax_b.text(x_pos, y0 + cell_h + 0.04, lbl, ha="center", va="bottom",
              fontsize=9, color=cc, fontweight="bold",
              transform=ax_b.transAxes)

# Species divider
ax_b.plot([0.57, 0.57], [0.12, 0.88], color=GREY, lw=1, ls="--",
          transform=ax_b.transAxes)
ax_b.text(0.35, 0.94, "Human", ha="center", va="bottom",
          fontsize=9.5, transform=ax_b.transAxes, color="#333333",
          fontweight="bold")
ax_b.text(0.74, 0.94, "Mouse", ha="center", va="bottom",
          fontsize=9.5, transform=ax_b.transAxes, color="#333333",
          fontweight="bold")

ax_b.set_xlim(0, 1); ax_b.set_ylim(0, 1)
ax_b.axis("off")
ax_b.set_title("Cross-species conservation", fontsize=11, pad=8, fontweight="bold")
ax_b.text(-0.05, 1.02, "B", transform=ax_b.transAxes,
          fontsize=15, fontweight="bold", va="top")

# Data source annotation
ax_b.text(0.5, 0.03, "Human: ReMap2022 THP-1 (GSE89431, GSE99887)\nMouse: ReMap2022 macrophage (GSE124722) / BMDM",
          ha="center", va="bottom", fontsize=7, color=GREY,
          transform=ax_b.transAxes, style="italic")

# ══════════════════════════════════════════════════════════════
# PANEL C  — 2-Brake model schematic
# ══════════════════════════════════════════════════════════════
ax_c.set_xlim(0, 10); ax_c.set_ylim(0, 10)
ax_c.axis("off")
ax_c.set_title("Dual nuclear receptor model", fontsize=11, pad=8, fontweight="bold")
ax_c.text(-0.05, 1.02, "C", transform=ax_c.transAxes,
          fontsize=15, fontweight="bold", va="top")

def box(ax, x, y, w, h, color, label, sublabel="", fontsize=9.5, alpha=0.92):
    rect = FancyBboxPatch((x - w/2, y - h/2), w, h,
                           boxstyle="round,pad=0.15",
                           facecolor=color, edgecolor="white",
                           linewidth=1.5, alpha=alpha, zorder=3)
    ax.add_patch(rect)
    ax.text(x, y + (0.15 if sublabel else 0), label,
            ha="center", va="center", fontsize=fontsize,
            fontweight="bold", color="white", zorder=4)
    if sublabel:
        ax.text(x, y - 0.35, sublabel, ha="center", va="center",
                fontsize=7.5, color="white", alpha=0.9, zorder=4)

def arrow(ax, x1, y1, x2, y2, color, lw=2):
    ax.annotate("", xy=(x2, y2), xytext=(x1, y1),
                arrowprops=dict(arrowstyle="-|>", color=color,
                                lw=lw, mutation_scale=14))

# ── Left side: VDR / TLR10 axis ──
box(ax_c, 2.5, 8.8, 2.8, 0.9, VDR_COL, "Vitamin D  →  VDR", fontsize=8.5)
box(ax_c, 2.5, 7.2, 2.8, 0.9, VDR_LITE, "TLR1 / TLR10", "GR-free zone", fontsize=9)
box(ax_c, 2.5, 5.6, 2.8, 0.9, "#1A6B3A", "Bacteroides", "hypoacyl LPS", fontsize=8.5)

arrow(ax_c, 2.5, 8.35, 2.5, 7.65, VDR_COL)
arrow(ax_c, 2.5, 6.75, 2.5, 6.05, "#1A6B3A")

# TLR10 + Bacteroides → TLR2/TLR10 complex
ax_c.annotate("", xy=(2.5, 4.75), xytext=(2.5, 5.15),
              arrowprops=dict(arrowstyle="-|>", color="#555555", lw=1.5))
box(ax_c, 2.5, 4.2, 2.8, 0.9, "#2CA25F", "TLR2 / TLR10 complex", fontsize=8)

arrow(ax_c, 2.5, 3.75, 2.5, 3.1, "#2CA25F", lw=2)

# Thermostat label
ax_c.text(0.6, 7.2, "🌡 Thermostat\n(gentle)", ha="center", va="center",
          fontsize=8, color=VDR_COL, style="italic",
          bbox=dict(boxstyle="round,pad=0.3", fc=BG_ZONE, ec=VDR_LITE, lw=1))

# ── Right side: GR / TLR6 axis ──
box(ax_c, 7.5, 8.8, 2.8, 0.9, GR_COL, "Steroids  →  GR", fontsize=8.5)
box(ax_c, 7.5, 7.2, 2.8, 0.9, GR_LITE, "TLR6 / TLR2", "strong brake", fontsize=9)
box(ax_c, 7.5, 5.6, 2.8, 0.9, "#B35806", "Pathogen", "acyl lipoproteins", fontsize=8.5)

arrow(ax_c, 7.5, 8.35, 7.5, 7.65, GR_COL)
arrow(ax_c, 7.5, 5.15, 7.5, 4.4, GR_COL)

ax_c.annotate("", xy=(7.5, 4.75), xytext=(7.5, 5.15),
              arrowprops=dict(arrowstyle="-|>", color=GR_COL, lw=1.5))
box(ax_c, 7.5, 4.2, 2.8, 0.9, GR_COL, "NF-κB suppression", "via GR direct", fontsize=8)

# Fire truck label
ax_c.text(9.4, 7.2, "🚒 Fire truck\n(strong)", ha="center", va="center",
          fontsize=8, color=GR_COL, style="italic",
          bbox=dict(boxstyle="round,pad=0.3", fc="#FFF0EE", ec=GR_LITE, lw=1))

# ── Center: A20 convergence ──
box(ax_c, 5.0, 2.7, 3.2, 0.9, "#525252", "A20 (TNFAIP3)", "NF-κB suppression", fontsize=9)

arrow(ax_c, 2.5, 3.75, 3.5, 3.05, "#2CA25F", lw=2)
arrow(ax_c, 7.5, 3.75, 6.5, 3.05, GR_COL, lw=2)

# ── Output: homeostasis ──
box(ax_c, 5.0, 1.3, 4.0, 0.9, "#252525", "Mucosal homeostasis", fontsize=9.5)
arrow(ax_c, 5.0, 2.25, 5.0, 1.75, "#252525", lw=2)

# Divider
ax_c.axvline(x=5.0, ymin=0.38, ymax=0.98, color=GREY, lw=1, ls=":", alpha=0.5)

# ── Figure title ──────────────────────────────────────────────
fig.suptitle(
    "VDR and GR establish a dual-brake transcriptional architecture at TLR loci",
    fontsize=13, fontweight="bold", y=0.98
)

fig.savefig(OUT, dpi=300, facecolor="white")
fig.savefig(OUT_PNG, dpi=150, facecolor="white")
plt.close()
print(f"Saved: {OUT}")
print(f"Saved: {OUT_PNG}")
