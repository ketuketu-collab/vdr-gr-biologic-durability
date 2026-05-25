#!/usr/bin/env python3
"""
Publication-quality figure: VDR occupies active enhancers upstream of NOD2 in human monocytes.

5-panel genomic track:
  Panel 1: H3K27ac signal (CD14+ monocyte, ENCODE ENCFF523ZCA) — active enhancer mark
  Panel 2: DNase-seq signal (THP-1, ENCODE ENCFF686JSS) — open chromatin
  Panel 3: VDR ChIP-seq peaks (ReMap2022; THP-1 highlighted)
  Panel 4: NOD2 gene model (minus strand)
  Panel 5: Schematic — VitD → VDR → enhancer → NOD2 suppression

Region: chr16:50,630,000-50,780,000  (≈150 kb, NOD2 locus + upstream enhancers)
NOD2 TSS: 50,756,049 (minus strand → transcription goes LEFT; upstream = higher coord)
"""

import gzip, re
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
import matplotlib.gridspec as gridspec
import pyBigWig

# ── File paths ────────────────────────────────────────────────────────────────
H3K27AC_BW  = "/Volumes/M4_SSD/projects/tlr10_chipseq_vdr_gr/data/encode_ENCFF523ZCA.bigWig"
DNASE_BW    = "/Volumes/M4_SSD/projects/tlr10_chipseq_vdr_gr/data/encode_ENCFF686JSS.bigWig"
VDR_BED_GZ  = "/Volumes/M4_SSD/projects/tlr10_chipseq_vdr_gr/data/remap2022_VDR_all_macs2_hg38.bed.gz"

OUT_PDF = "/Volumes/M4_SSD/projects/tlr_chipseq/results/figures/NOD2_VDR_enhancer_final_v2.pdf"
OUT_PNG = "/Volumes/M4_SSD/projects/tlr_chipseq/results/figures/NOD2_VDR_enhancer_final_v2.png"

# ── Locus ─────────────────────────────────────────────────────────────────────
CHROM      = "chr16"
LOC_START  = 50_630_000
LOC_END    = 50_780_000
NOD2_TSS   = 50_756_049   # minus-strand: transcribes rightward→leftward
NOD2_TES   = 50_700_000   # approximate gene end (minus strand)

# Enhancer regions (from previous analysis; E2a+E2b merged into E2)
ENHANCERS = [
    {"name": "E1", "start": 50_689_000, "end": 50_700_000,
     "score": 126, "label": "E1\n(−60 kb)"},
    {"name": "E2", "start": 50_706_000, "end": 50_720_000,
     "score": 45,  "label": "E2\n(−43–47 kb)"},
]

# ── Colors ────────────────────────────────────────────────────────────────────
COL_H3K27AC = "#E08028"   # orange
COL_DNASE   = "#4A90D9"   # blue
COL_VDR_THP = "#C0392B"   # red — THP-1 peaks
COL_VDR_OTH = "#BBBBBB"   # gray — other cell types
COL_ENH     = "#FDAE6B"   # enhancer highlight
COL_GENE    = "#2C3E50"   # gene body
COL_ARROW   = "#27AE60"   # chromatin loop arrow


# ── Utility: read bigWig signal ───────────────────────────────────────────────
def get_signal(bw_path, chrom, start, end, nbins=1500):
    bw = pyBigWig.open(bw_path)
    vals = bw.stats(chrom, start, end, type="mean", nBins=nbins)
    bw.close()
    vals = np.array([v if v is not None else 0.0 for v in vals])
    return vals


# ── Utility: read VDR peaks ───────────────────────────────────────────────────
def get_vdr_peaks(bed_gz, chrom, start, end):
    peaks = []
    with gzip.open(bed_gz, "rt") as fh:
        for line in fh:
            if line.startswith("#"): continue
            if not line.startswith(chrom): continue
            cols = line.rstrip().split("\t")
            if len(cols) < 5: continue
            if cols[0] != chrom: continue
            s, e = int(cols[1]), int(cols[2])
            if e < start or s > end: continue
            name = cols[3]
            try:
                score = float(cols[4])
            except:
                score = 0.0
            summit = int(cols[6]) if len(cols) > 6 else (s + e) // 2
            is_thp1 = "THP" in name or "THP-1" in name
            peaks.append({"start": s, "end": e, "score": score,
                          "summit": summit, "name": name, "is_thp1": is_thp1})
    return sorted(peaks, key=lambda x: x["score"])  # lowest first → highest on top


# ── Load data ─────────────────────────────────────────────────────────────────
print("Loading bigWig signals …")
x_bins = np.linspace(LOC_START, LOC_END, 1500)
h3k27ac = get_signal(H3K27AC_BW, CHROM, LOC_START, LOC_END)
dnase    = get_signal(DNASE_BW,   CHROM, LOC_START, LOC_END)

print("Loading VDR peaks …")
vdr_peaks = get_vdr_peaks(VDR_BED_GZ, CHROM, LOC_START, LOC_END)
print(f"  Found {len(vdr_peaks)} VDR peaks in region; "
      f"{sum(p['is_thp1'] for p in vdr_peaks)} THP-1")

# ── Figure layout ─────────────────────────────────────────────────────────────
fig = plt.figure(figsize=(14, 10))
fig.patch.set_facecolor("white")

gs = gridspec.GridSpec(
    5, 1,
    height_ratios=[2.5, 1.8, 2.0, 1.0, 1.4],
    hspace=0.06,
    left=0.12, right=0.96, top=0.93, bottom=0.07
)

axes = [fig.add_subplot(gs[i]) for i in range(5)]

# Shared x-axis
for ax in axes[:-1]:
    ax.set_xlim(LOC_START, LOC_END)
    ax.tick_params(labelbottom=False)
axes[-1].set_xlim(LOC_START, LOC_END)

# ── Enhancer shading (all panels) ─────────────────────────────────────────────
for ax in axes[:4]:
    for enh in ENHANCERS:
        ax.axvspan(enh["start"], enh["end"],
                   color=COL_ENH, alpha=0.25, zorder=0)

# ── Panel 0: H3K27ac ─────────────────────────────────────────────────────────
ax0 = axes[0]
ax0.fill_between(x_bins, h3k27ac, color=COL_H3K27AC, alpha=0.85, linewidth=0)
ax0.set_ylim(0, max(h3k27ac.max() * 1.15, 1))
ax0.set_ylabel("Signal\n(rpm)", fontsize=8.5, labelpad=4)
ax0.spines[["top","right","bottom"]].set_visible(False)
ax0.tick_params(left=True, labelsize=8)

# Track label
ax0.text(0.005, 0.88, "H3K27ac  (CD14⁺ monocyte)", transform=ax0.transAxes,
         fontsize=9.5, color=COL_H3K27AC, fontweight="bold", va="top")
ax0.text(0.005, 0.70, "Active enhancer mark  ·  ENCODE ENCFF523ZCA",
         transform=ax0.transAxes, fontsize=7.5, color="#666666", va="top",
         style="italic")

# ── Panel 1: DNase-seq ───────────────────────────────────────────────────────
ax1 = axes[1]
ax1.fill_between(x_bins, dnase, color=COL_DNASE, alpha=0.85, linewidth=0)
ax1.set_ylim(0, max(dnase.max() * 1.15, 1))
ax1.set_ylabel("Signal\n(rpm)", fontsize=8.5, labelpad=4)
ax1.spines[["top","right","bottom"]].set_visible(False)
ax1.tick_params(left=True, labelsize=8)

ax1.text(0.005, 0.88, "DNase-seq  (THP-1)", transform=ax1.transAxes,
         fontsize=9.5, color=COL_DNASE, fontweight="bold", va="top")
ax1.text(0.005, 0.70, "Open chromatin  ·  ENCODE ENCFF686JSS",
         transform=ax1.transAxes, fontsize=7.5, color="#666666", va="top",
         style="italic")

# ── Panel 2: VDR ChIP-seq peaks ──────────────────────────────────────────────
ax2 = axes[2]
ax2.set_ylim(0, 160)
ax2.set_ylabel("Score", fontsize=8.5, labelpad=4)
ax2.spines[["top","right","bottom"]].set_visible(False)
ax2.tick_params(left=True, labelsize=8)
ax2.axhline(0, color="#CCCCCC", lw=0.5)

ax2.text(0.005, 0.95, "VDR ChIP-seq peaks  (ReMap2022)",
         transform=ax2.transAxes, fontsize=9.5, fontweight="bold",
         color=COL_VDR_THP, va="top")

# Score threshold line
ax2.axhline(20, color="#DDDDDD", lw=0.8, ls="--", zorder=0)
ax2.text(LOC_START + 2000, 21, "score=20", fontsize=6.5, color="#AAAAAA")

# Plot peaks as rectangles; stack THP-1 on top of others
PEAK_H = 12  # peak rectangle height in score units
labeled_thp1 = set()
labeled_scores = []  # (x_center, score, label_txt) to avoid overlap

for pk in vdr_peaks:
    s = max(pk["start"], LOC_START)
    e = min(pk["end"],   LOC_END)
    w = e - s
    y_bot = 2  # all peaks at y=2 baseline
    color = COL_VDR_THP if pk["is_thp1"] else COL_VDR_OTH
    height = max(pk["score"] * 0.85, 3)
    rect = mpatches.Rectangle(
        (s, y_bot), w, height,
        facecolor=color, edgecolor="none",
        alpha=0.9 if pk["is_thp1"] else 0.5, zorder=3
    )
    ax2.add_patch(rect)

    # Label high-scoring THP-1 peaks (avoid overlap)
    if pk["is_thp1"] and pk["score"] >= 30:
        xc = (s + e) / 2
        # Check distance from existing labels
        too_close = any(abs(xc - lx) < 8000 and abs(pk["score"] - ls) < 20
                        for lx, ls, _ in labeled_scores)
        if not too_close:
            ax2.text(xc, y_bot + height + 2,
                     f"{pk['score']:.0f}",
                     ha="center", va="bottom", fontsize=7.5,
                     color=COL_VDR_THP, fontweight="bold", zorder=5)
            labeled_scores.append((xc, pk["score"], str(pk["score"])))

# Legend
leg_patches = [
    mpatches.Patch(color=COL_VDR_THP, alpha=0.9, label="VDR peak — THP-1"),
    mpatches.Patch(color=COL_VDR_OTH, alpha=0.7, label="VDR peak — other"),
]
ax2.legend(handles=leg_patches, fontsize=7.5, frameon=False,
           loc="upper right", bbox_to_anchor=(0.99, 0.96))

# ── Enhancer labels with vertical lines ──────────────────────────────────────
for enh in ENHANCERS:
    xm = (enh["start"] + enh["end"]) / 2
    ax2.axvline(xm, ymin=0.0, ymax=0.95, color=COL_ENH, lw=1.2,
                ls=":", alpha=0.8, zorder=1)
    ax2.text(xm, 150, enh["label"], ha="center", va="top",
             fontsize=8, color="#C05000", fontweight="bold",
             bbox=dict(boxstyle="round,pad=0.2", fc="white", ec=COL_ENH, lw=1))

# ── Panel 3: NOD2 gene model ──────────────────────────────────────────────────
ax3 = axes[3]
ax3.set_ylim(-0.5, 1.5)
ax3.axis("off")

# Gene body line
gene_y = 0.5
ax3.plot([NOD2_TES, NOD2_TSS], [gene_y, gene_y],
         color=COL_GENE, lw=2.5, solid_capstyle="round", zorder=2)

# Exon blocks (approximate positions)
exons = [
    (50_752_000, 50_756_049),  # Exon 1 (near TSS, minus strand)
    (50_744_000, 50_745_500),
    (50_741_000, 50_742_000),
    (50_737_000, 50_738_500),
    (50_733_500, 50_735_000),
    (50_729_000, 50_731_000),
    (50_725_000, 50_727_000),
    (50_721_000, 50_723_000),
    (50_717_000, 50_719_000),
    (50_713_000, 50_715_000),
    (50_709_000, 50_711_500),
    (50_703_000, 50_706_000),
]
exon_h = 0.35
for exs, exe in exons:
    exs = max(exs, LOC_START)
    exe = min(exe, LOC_END)
    rect = mpatches.Rectangle(
        (exs, gene_y - exon_h/2), exe - exs, exon_h,
        facecolor=COL_GENE, edgecolor="none", zorder=3
    )
    ax3.add_patch(rect)

# Transcription direction arrow (minus strand → rightward gene, transcription leftward)
ax3.annotate("", xy=(NOD2_TES - 2000, gene_y),
             xytext=(NOD2_TSS + 2000, gene_y),
             arrowprops=dict(arrowstyle="-|>", color=COL_GENE, lw=1.5,
                             mutation_scale=14))

# Gene name — place above the gene body, offset from loop arcs
ax3.text(NOD2_TSS - 28000, gene_y - 0.38,
         "NOD2",
         ha="center", va="top", fontsize=11, fontweight="bold",
         color=COL_GENE)
ax3.text(NOD2_TSS - 28000, gene_y - 0.62,
         "(minus strand, chr16)",
         ha="center", va="top", fontsize=7.5, color="#666666")

# TSS marker
ax3.axvline(NOD2_TSS, ymin=0, ymax=1, color="#888888", lw=1, ls="--", alpha=0.6)
ax3.text(NOD2_TSS + 1500, 1.32, "TSS", ha="left", va="top",
         fontsize=7.5, color="#888888")

# Chromatin loop arcs from each enhancer to TSS
from matplotlib.patches import Arc
for i, enh in enumerate(ENHANCERS):
    xm = (enh["start"] + enh["end"]) / 2
    mid_x = (xm + NOD2_TSS) / 2
    width = abs(NOD2_TSS - xm)
    arc_h = 0.55 + i * 0.18  # stagger arc heights
    arc = Arc((mid_x, gene_y), width, arc_h,
              angle=0, theta1=0, theta2=180,
              color=COL_ARROW, lw=1.6, ls="-", alpha=0.75, zorder=4)
    ax3.add_patch(arc)

# Single arrowhead near TSS
ax3.annotate("", xy=(NOD2_TSS - 2000, gene_y + 0.04),
             xytext=(NOD2_TSS - 7000, gene_y + 0.32),
             arrowprops=dict(arrowstyle="-|>", color=COL_ARROW,
                             lw=1.4, mutation_scale=11))

# Loop label — placed in clear space
ax3.text(NOD2_TSS - 52000, gene_y + 0.70,
         "Chromatin loop\n(enhancer → promoter)",
         ha="center", va="bottom", fontsize=7.5, color=COL_ARROW,
         style="italic",
         bbox=dict(boxstyle="round,pad=0.25", fc="white",
                   ec=COL_ARROW, lw=0.9, alpha=0.85))

# Enhancer labels below gene track
for enh in ENHANCERS:
    xm = (enh["start"] + enh["end"]) / 2
    ax3.text(xm, -0.38, enh["label"].split("\n")[1],
             ha="center", va="top", fontsize=7.5,
             color="#C05000", fontweight="bold")

ax3.set_xlim(LOC_START, LOC_END)

# ── Panel 4: Genomic ruler / x-axis ──────────────────────────────────────────
ax4 = axes[4]
ax4.set_xlim(LOC_START, LOC_END)
ax4.set_ylim(0, 1)
ax4.axis("off")

# Scale bar 10 kb
sb_start = LOC_START + 5000
sb_end   = sb_start + 10_000
ax4.plot([sb_start, sb_end], [0.85, 0.85], color="black", lw=2,
         solid_capstyle="butt")
ax4.plot([sb_start, sb_start], [0.75, 0.95], color="black", lw=1.5)
ax4.plot([sb_end,   sb_end],   [0.75, 0.95], color="black", lw=1.5)
ax4.text((sb_start + sb_end)/2, 0.55, "10 kb",
         ha="center", va="top", fontsize=9, color="black")

# Chromosome label
ax4.text(LOC_END - 2000, 0.85,
         f"{CHROM}:{LOC_START:,}–{LOC_END:,}  (hg38)",
         ha="right", va="center", fontsize=8, color="#555555", style="italic")

# ── Enhancer shading labels (top of Panel 0) ─────────────────────────────────
for enh in ENHANCERS:
    xm = (enh["start"] + enh["end"]) / 2
    axes[0].axvline(xm, color=COL_ENH, lw=1.5, ls=":", alpha=0.5, zorder=0)
    axes[0].text(xm, axes[0].get_ylim()[1] * 0.97,
                 enh["label"].split("\n")[0],
                 ha="center", va="top", fontsize=8,
                 color="#C05000", fontweight="bold")

# ── NOD2 TSS marker on all panels ────────────────────────────────────────────
for ax in axes[:3]:
    ax.axvline(NOD2_TSS, color="#888888", lw=1, ls="--", alpha=0.5)

# ── Figure title ─────────────────────────────────────────────────────────────
fig.suptitle(
    "VDR occupies active enhancers upstream of NOD2 in human monocytes",
    fontsize=13, fontweight="bold", y=0.97
)

# Subtitle
fig.text(0.5, 0.945,
         "H3K27ac-marked open chromatin (ENCODE) co-localizes with VDR peaks (ReMap2022, THP-1 + 1,25(OH)₂D₃)",
         ha="center", fontsize=9, color="#555555", style="italic")

# ── Save ─────────────────────────────────────────────────────────────────────
print("Saving …")
fig.savefig(OUT_PDF, bbox_inches="tight", facecolor="white")
fig.savefig(OUT_PNG, dpi=200, bbox_inches="tight", facecolor="white")
plt.close()
print(f"PDF → {OUT_PDF}")
print(f"PNG → {OUT_PNG}")
