#!/usr/bin/env python3
"""
BigWig pile-up visualization: VDR ChIP-seq signal at TLR loci (hg38).
Run after 06_align_vdr.sh completes.
Also overlays ReMap2022 peak annotations.
"""
import gzip
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import pyBigWig

PROJ    = Path("/Volumes/M4_SSD/projects/tlr_chipseq")
FIG_DIR = PROJ / "results" / "figures"
FIG_DIR.mkdir(parents=True, exist_ok=True)

BIGWIGS = {
    "VDR +VitD (THP-1)": PROJ/"results/coverage/vitd_merged.bw",
    "VDR −VitD (THP-1)": PROJ/"results/coverage/vehicle_merged.bw",
}

LOCI = {
    "TLR10":   ("chr4",  38767648,  "-", "#E74C3C"),
    "TLR1":    ("chr4",  38770059,  "+", "#3498DB"),
    "TLR6":    ("chr4",  38806672,  "+", "#2ECC71"),
    "TLR2":    ("chr4",  153684080, "+", "#E67E22"),
    "TLR4":    ("chr9",  117704402, "+", "#9B59B6"),
    "TNFAIP3": ("chr6",  137866317, "+", "#1ABC9C"),
    "DUSP1":   ("chr5",  172768090, "+", "#34495E"),
}
FLANK = 15000
NBINS = 1000

VDR_BED = Path("/Volumes/M4_SSD/projects/tlr10_chipseq_vdr_gr/data/remap2022_VDR_all_macs2_hg38.bed.gz")

def get_vdr_thp1_peaks(chrom, start, end):
    peaks = []
    with gzip.open(VDR_BED, "rt") as f:
        for line in f:
            cols = line.rstrip().split("\t")
            if cols[0] != chrom: continue
            s, e = int(cols[1]), int(cols[2])
            if e < start or s > end: continue
            if "THP" not in cols[3]: continue
            sc = float(cols[4]) if cols[4]!='.' else 0.0
            summit = int(cols[6]) if len(cols)>6 else (s+e)//2
            peaks.append({"start":s,"end":e,"score":sc,"summit":summit,"name":cols[3]})
    return peaks


fig, axes = plt.subplots(len(LOCI), 1, figsize=(14, 3.5*len(LOCI)),
                          gridspec_kw={"hspace":0.7})

bw_handles = {}
for label, bw_path in BIGWIGS.items():
    if bw_path.exists():
        bw_handles[label] = pyBigWig.open(str(bw_path))

if not bw_handles:
    print("No bigWig files found — run 06_align_vdr.sh first")
    print("(Creating placeholder figure)")

colors = {"VDR +VitD (THP-1)": "#E67E22", "VDR −VitD (THP-1)": "#AAAAAA"}

for ax, (gene, (chrom, tss, strand, gcolor)) in zip(axes, LOCI.items()):
    reg_start = max(0, tss - FLANK)
    reg_end   = tss + FLANK
    x_pos = np.linspace(reg_start, reg_end, NBINS)

    for label, bw in bw_handles.items():
        try:
            vals = bw.stats(chrom, reg_start, reg_end, nBins=NBINS)
            vals = np.array([v if v is not None else 0.0 for v in vals])
            ax.fill_between(x_pos, vals, alpha=0.5, color=colors[label], label=label)
            ax.plot(x_pos, vals, lw=0.5, color=colors[label], alpha=0.7)
        except Exception:
            pass

    # Overlay THP-1 VDR peaks as triangles
    peaks = get_vdr_thp1_peaks(chrom, reg_start, reg_end)
    ymax = ax.get_ylim()[1] if ax.get_ylim()[1] > 0 else 1.0
    for pk in peaks:
        ax.axvspan(pk["start"], pk["end"], alpha=0.2, color="#E67E22", zorder=1)
        ax.plot(pk["summit"], ymax*0.9, "v", color="#E67E22", ms=8, zorder=5)
        ax.text(pk["summit"], ymax*0.92, f"VDR\n{pk['score']:.1f}",
                ha="center", va="bottom", fontsize=6.5, color="darkorange", zorder=6)

    # TSS line
    ax.axvline(tss, color=gcolor, lw=1.5, linestyle="--", alpha=0.8, zorder=3)
    ax.axvspan(tss-2000, tss+2000, alpha=0.05, color=gcolor, zorder=1)
    # Gene body bar
    gene_end = tss + (1000 if strand=="+" else -1000)
    ybot = ax.get_ylim()[0] - 0.05 * abs(ax.get_ylim()[1]-ax.get_ylim()[0])

    ax.set_xlim(reg_start, reg_end)
    ax.set_ylabel("CPM", fontsize=8)
    ax.set_title(f"{gene}  [{chrom}:{tss:,}  {strand} strand]", fontsize=9,
                 color=gcolor, fontweight="bold")
    ax.set_xlabel(f"Position ({chrom})", fontsize=7)
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)

    xt = np.linspace(reg_start, reg_end, 7).astype(int)
    ax.set_xticks(xt)
    ax.set_xticklabels([f"{t:,}" for t in xt], fontsize=6.5)

    if ax == axes[0]:
        ax.legend(fontsize=8, loc="upper right")

for bw in bw_handles.values():
    bw.close()

fig.suptitle("VDR ChIP-seq signal tracks: THP-1 cells ±VitD (hg38, GSE89431)\n"
             "▼ = VDR peak (ReMap2022, THP-1 specific)",
             fontsize=10, y=1.01)

plt.savefig(FIG_DIR / "vdr_bigwig_pileup.pdf", bbox_inches="tight")
plt.savefig(FIG_DIR / "vdr_bigwig_pileup.png", dpi=300, bbox_inches="tight")
plt.close()
print(f"Figure → {FIG_DIR/'vdr_bigwig_pileup.pdf'}")
print("Done.")
