#!/usr/bin/env python3
"""
TLR10 locus: VDR occupies active enhancers (parallel to NOD2 figure)
Region: chr4:38,720,000-38,780,000
TLR10 TSS: 38,767,648 (minus strand)
Key enhancers:
  E1: -23.5 kb (summit 38,744,097)  VDR score=28.2, H3K27ac=16
  E2: -5.9 kb  (summit 38,761,739)  VDR score=13.2, H3K27ac=151
"""
import gzip
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.gridspec as gridspec
from matplotlib.patches import Arc
import pyBigWig

H3K27AC_BW = "/Volumes/M4_SSD/projects/tlr10_chipseq_vdr_gr/data/encode_ENCFF523ZCA.bigWig"
DNASE_BW   = "/Volumes/M4_SSD/projects/tlr10_chipseq_vdr_gr/data/encode_ENCFF686JSS.bigWig"
VDR_BED_GZ = "/Volumes/M4_SSD/projects/tlr10_chipseq_vdr_gr/data/remap2022_VDR_all_macs2_hg38.bed.gz"

OUT_PDF = "/Volumes/M4_SSD/projects/tlr_chipseq/results/figures/TLR10_VDR_enhancer_final.pdf"
OUT_PNG = "/Volumes/M4_SSD/projects/tlr_chipseq/results/figures/TLR10_VDR_enhancer_final.png"

CHROM     = "chr4"
LOC_START = 38_718_000
LOC_END   = 38_780_000
TLR10_TSS = 38_767_648  # minus strand
TLR10_TES = 38_728_000  # approximate

ENHANCERS = [
    {"name": "E1", "start": 38_743_200, "end": 38_745_000,
     "score": 28.2, "label": "E1\n(−24 kb)"},
    {"name": "E2", "start": 38_760_800, "end": 38_762_200,
     "score": 13.2, "label": "E2\n(−6 kb)"},
]

COL_H3K27AC = "#E08028"
COL_DNASE   = "#4A90D9"
COL_VDR_THP = "#C0392B"
COL_VDR_OTH = "#BBBBBB"
COL_ENH     = "#FDAE6B"
COL_GENE    = "#2C3E50"
COL_ARROW   = "#27AE60"

def get_signal(bw_path, chrom, start, end, nbins=1000):
    bw = pyBigWig.open(bw_path)
    vals = bw.stats(chrom, start, end, type="mean", nBins=nbins)
    bw.close()
    return np.array([v if v is not None else 0.0 for v in vals])

def get_vdr_peaks(bed_gz, chrom, start, end):
    peaks = []
    with gzip.open(bed_gz, "rt") as fh:
        for line in fh:
            if not line.startswith(chrom[:4]): continue
            cols = line.rstrip().split("\t")
            if cols[0] != chrom: continue
            s, e = int(cols[1]), int(cols[2])
            if e < start or s > end: continue
            name = cols[3]
            try: score = float(cols[4])
            except: score = 0.0
            summit = int(cols[6]) if len(cols) > 6 else (s+e)//2
            is_thp1 = "THP" in name
            peaks.append({"start": s, "end": e, "score": score,
                          "summit": summit, "name": name, "is_thp1": is_thp1})
    return sorted(peaks, key=lambda x: x["score"])

print("Loading signals …")
x_bins   = np.linspace(LOC_START, LOC_END, 1000)
h3k27ac  = get_signal(H3K27AC_BW, CHROM, LOC_START, LOC_END)
dnase    = get_signal(DNASE_BW,   CHROM, LOC_START, LOC_END)
vdr_peaks = get_vdr_peaks(VDR_BED_GZ, CHROM, LOC_START, LOC_END)
print(f"  VDR peaks: {len(vdr_peaks)}, THP-1: {sum(p['is_thp1'] for p in vdr_peaks)}")

# ── Figure ────────────────────────────────────────────────────
fig = plt.figure(figsize=(13, 9))
fig.patch.set_facecolor("white")
gs = gridspec.GridSpec(4, 1, height_ratios=[2.5, 1.8, 1.8, 1.0],
                       hspace=0.06, left=0.12, right=0.96, top=0.93, bottom=0.07)
axes = [fig.add_subplot(gs[i]) for i in range(4)]

for ax in axes:
    ax.set_xlim(LOC_START, LOC_END)
    ax.tick_params(labelbottom=False)

# Enhancer shading
for ax in axes[:3]:
    for enh in ENHANCERS:
        ax.axvspan(enh["start"], enh["end"], color=COL_ENH, alpha=0.28, zorder=0)

# ── Panel 0: H3K27ac ────────────────────────────────
ax0 = axes[0]
ax0.fill_between(x_bins, h3k27ac, color=COL_H3K27AC, alpha=0.85, lw=0)
ax0.set_ylim(0, max(h3k27ac.max()*1.15, 1))
ax0.set_ylabel("Signal\n(rpm)", fontsize=8.5)
ax0.spines[["top","right","bottom"]].set_visible(False)
ax0.tick_params(left=True, labelsize=8)
ax0.text(0.005, 0.90, "H3K27ac  (CD14⁺ monocyte)", transform=ax0.transAxes,
         fontsize=9.5, color=COL_H3K27AC, fontweight="bold", va="top")
ax0.text(0.005, 0.72, "Active enhancer mark  ·  ENCODE ENCFF523ZCA",
         transform=ax0.transAxes, fontsize=7.5, color="#666666", va="top", style="italic")
# Enhancer labels
for enh in ENHANCERS:
    xm = (enh["start"]+enh["end"])/2
    ax0.text(xm, ax0.get_ylim()[1]*0.96, enh["label"].split("\n")[0],
             ha="center", va="top", fontsize=8.5, color="#C05000", fontweight="bold")
    ax0.axvline(xm, color=COL_ENH, lw=1.5, ls=":", alpha=0.6, zorder=1)

# ── Panel 1: DNase ────────────────────────────────
ax1 = axes[1]
ax1.fill_between(x_bins, dnase, color=COL_DNASE, alpha=0.85, lw=0)
ax1.set_ylim(0, max(dnase.max()*1.15, 1))
ax1.set_ylabel("Signal\n(rpm)", fontsize=8.5)
ax1.spines[["top","right","bottom"]].set_visible(False)
ax1.tick_params(left=True, labelsize=8)
ax1.text(0.005, 0.90, "DNase-seq  (THP-1)", transform=ax1.transAxes,
         fontsize=9.5, color=COL_DNASE, fontweight="bold", va="top")
ax1.text(0.005, 0.72, "Open chromatin  ·  ENCODE ENCFF686JSS",
         transform=ax1.transAxes, fontsize=7.5, color="#666666", va="top", style="italic")
for enh in ENHANCERS:
    ax1.axvline((enh["start"]+enh["end"])/2, color=COL_ENH, lw=1.5, ls=":", alpha=0.5)

# ── Panel 2: VDR peaks ─────────────────────────────
ax2 = axes[2]
ax2.set_ylim(0, 50)
ax2.set_ylabel("Score", fontsize=8.5)
ax2.spines[["top","right","bottom"]].set_visible(False)
ax2.tick_params(left=True, labelsize=8)
ax2.axhline(0, color="#CCCCCC", lw=0.5)
ax2.axhline(10, color="#EEEEEE", lw=0.8, ls="--", zorder=0)
ax2.text(LOC_START+1000, 10.5, "score=10", fontsize=6.5, color="#AAAAAA")
ax2.text(0.005, 0.95, "VDR ChIP-seq peaks  (ReMap2022)",
         transform=ax2.transAxes, fontsize=9.5, fontweight="bold", color=COL_VDR_THP, va="top")

labeled = []
for pk in vdr_peaks:
    s = max(pk["start"], LOC_START)
    e = min(pk["end"],   LOC_END)
    h = max(pk["score"]*0.9, 2)
    color = COL_VDR_THP if pk["is_thp1"] else COL_VDR_OTH
    rect = mpatches.Rectangle((s, 2), e-s, h,
                               facecolor=color, edgecolor="none",
                               alpha=0.9 if pk["is_thp1"] else 0.5, zorder=3)
    ax2.add_patch(rect)
    if pk["is_thp1"] and pk["score"] >= 10:
        xc = (s+e)/2
        if not any(abs(xc-lx)<5000 for lx,_ in labeled):
            ax2.text(xc, 2+h+1, f"{pk['score']:.1f}",
                     ha="center", va="bottom", fontsize=8, color=COL_VDR_THP,
                     fontweight="bold", zorder=5)
            labeled.append((xc, pk["score"]))

# Enhancer labels + vertical lines
for enh in ENHANCERS:
    xm = (enh["start"]+enh["end"])/2
    ax2.axvline(xm, color=COL_ENH, lw=1.2, ls=":", alpha=0.8)
    ax2.text(xm, 46, enh["label"],
             ha="center", va="top", fontsize=8, color="#C05000", fontweight="bold",
             bbox=dict(boxstyle="round,pad=0.2", fc="white", ec=COL_ENH, lw=1))

leg = [mpatches.Patch(color=COL_VDR_THP, alpha=0.9, label="VDR — THP-1"),
       mpatches.Patch(color=COL_VDR_OTH, alpha=0.7, label="VDR — other")]
ax2.legend(handles=leg, fontsize=7.5, frameon=False, loc="upper right")

# TSS marker on all panels
for ax in axes[:3]:
    ax.axvline(TLR10_TSS, color="#888888", lw=1, ls="--", alpha=0.5)

# ── Panel 3: Gene model ───────────────────────────
ax3 = axes[3]
ax3.set_xlim(LOC_START, LOC_END)
ax3.set_ylim(-0.6, 1.6)
ax3.axis("off")
gene_y = 0.5

# Gene body
ax3.plot([TLR10_TES, TLR10_TSS], [gene_y, gene_y],
         color=COL_GENE, lw=2.5, solid_capstyle="round", zorder=2)

# Exons (approximate, minus strand)
exons = [
    (38_766_200, 38_767_648),  # Exon 1 (TSS)
    (38_759_500, 38_760_700),
    (38_756_000, 38_757_000),
    (38_752_000, 38_753_500),
    (38_748_000, 38_749_500),
    (38_744_000, 38_745_200),
    (38_740_000, 38_741_500),
    (38_736_000, 38_737_500),
    (38_732_000, 38_733_500),
    (38_728_000, 38_730_000),
]
for exs, exe in exons:
    exs = max(exs, LOC_START); exe = min(exe, LOC_END)
    ax3.add_patch(mpatches.Rectangle((exs, gene_y-0.18), exe-exs, 0.36,
                                      facecolor=COL_GENE, edgecolor="none", zorder=3))

# Direction arrow (minus strand → transcription leftward)
ax3.annotate("", xy=(TLR10_TES-2000, gene_y),
             xytext=(TLR10_TSS+2000, gene_y),
             arrowprops=dict(arrowstyle="-|>", color=COL_GENE, lw=1.5, mutation_scale=13))

# Chromatin loop arcs
for i, enh in enumerate(ENHANCERS):
    xm = (enh["start"]+enh["end"])/2
    mid_x = (xm+TLR10_TSS)/2
    width = abs(TLR10_TSS-xm)
    arc_h = 0.50 + i*0.22
    ax3.add_patch(Arc((mid_x, gene_y), width, arc_h,
                      angle=0, theta1=0, theta2=180,
                      color=COL_ARROW, lw=1.6, alpha=0.75, zorder=4))

ax3.annotate("", xy=(TLR10_TSS-2000, gene_y+0.05),
             xytext=(TLR10_TSS-6000, gene_y+0.28),
             arrowprops=dict(arrowstyle="-|>", color=COL_ARROW, lw=1.4, mutation_scale=11))

# Gene label
ax3.text((TLR10_TSS+TLR10_TES)/2, gene_y-0.38,
         "TLR10  (minus strand, chr4)",
         ha="center", va="top", fontsize=10, fontweight="bold", color=COL_GENE)

# TSS label
ax3.text(TLR10_TSS+500, 1.35, "TSS", ha="left", va="top", fontsize=7.5, color="#888888")

# Loop label
ax3.text(TLR10_TSS-30000, gene_y+0.72,
         "Chromatin loop\n(enhancer → promoter)",
         ha="center", va="bottom", fontsize=7.5, color=COL_ARROW, style="italic",
         bbox=dict(boxstyle="round,pad=0.25", fc="white", ec=COL_ARROW, lw=0.9, alpha=0.85))

# Scale bar
sb_s = LOC_START + 3000
sb_e = sb_s + 10_000
ax3.plot([sb_s, sb_e], [0.85, 0.85], color="black", lw=2, solid_capstyle="butt",
         transform=ax3.get_xaxis_transform(), clip_on=False)
ax3.text((sb_s+sb_e)/2, 0.55, "10 kb", ha="center", va="top", fontsize=9, color="black",
         transform=ax3.get_xaxis_transform(), clip_on=False)
ax3.text(LOC_END-1000, 0.85,
         f"{CHROM}:{LOC_START:,}–{LOC_END:,}  (hg38)",
         ha="right", va="center", fontsize=7.5, color="#555555", style="italic",
         transform=ax3.get_xaxis_transform(), clip_on=False)

# ── Title ─────────────────────────────────────────
fig.suptitle(
    "VDR occupies active enhancers upstream of TLR10 in human monocytes",
    fontsize=13, fontweight="bold", y=0.97)
fig.text(0.5, 0.945,
         "H3K27ac-marked open chromatin (ENCODE) co-localizes with VDR peaks (ReMap2022, THP-1 + 1,25(OH)₂D₃)",
         ha="center", fontsize=9, color="#555555", style="italic")

fig.savefig(OUT_PDF, bbox_inches="tight", facecolor="white")
fig.savefig(OUT_PNG, dpi=200, bbox_inches="tight", facecolor="white")
plt.close()
print(f"PDF → {OUT_PDF}")
print(f"PNG → {OUT_PNG}")
