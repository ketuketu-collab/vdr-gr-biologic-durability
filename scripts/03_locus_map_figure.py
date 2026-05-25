#!/usr/bin/env python3
"""
Locus map figure: VDR and GR/NR3C1 peak positions at TLR loci.
Uses ReMap2022 aggregated peaks (hg38). Shows peak score, study, and cell type.
Highlights THP-1 / macrophage / monocyte / B-cell peaks.
"""
import gzip, re
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D

BASE_OLD = Path("/Volumes/M4_SSD/projects/tlr10_chipseq_vdr_gr/data")
FIG_DIR  = Path("/Volumes/M4_SSD/projects/tlr_chipseq/results/figures")
FIG_DIR.mkdir(parents=True, exist_ok=True)

VDR_BED   = BASE_OLD / "remap2022_VDR_all_macs2_hg38.bed.gz"
NR3C1_BED = BASE_OLD / "remap2022_NR3C1_all_macs2_hg38.bed.gz"

# hg38 gene bodies and TSS (strand-aware)
GENES = {
    "TLR10":   {"chrom": "chr4",  "tss": 38767648,  "end": 38739275,  "strand": "-", "color": "#E74C3C"},
    "TLR1":    {"chrom": "chr4",  "tss": 38770059,  "end": 38795755,  "strand": "+", "color": "#3498DB"},
    "TLR6":    {"chrom": "chr4",  "tss": 38806672,  "end": 38831036,  "strand": "+", "color": "#2ECC71"},
    "TLR2":    {"chrom": "chr4",  "tss": 153684080, "end": 153705699, "strand": "+", "color": "#E67E22"},
    "TLR4":    {"chrom": "chr9",  "tss": 117704402, "end": 117724735, "strand": "+", "color": "#9B59B6"},
    "TNFAIP3": {"chrom": "chr6",  "tss": 137866317, "end": 137883313, "strand": "+", "color": "#1ABC9C"},
    "DUSP1":   {"chrom": "chr5",  "tss": 172768090, "end": 172771195, "strand": "+", "color": "#34495E"},
}

WINDOW = 15000  # ±15kb for display

IMMUNE_COLORS = {
    "THP-1": "#E74C3C",
    "macrophage": "#C0392B",
    "monocyte":   "#E67E22",
    "GM12878":    "#8E44AD",
    "ALL":        "#9B59B6",
    "SUP-B15":    "#9B59B6",
    "NALM":       "#9B59B6",
    "A-549":      "#BDC3C7",
    "U2OS":       "#BDC3C7",
    "default":    "#95A5A6",
}

def cell_color(name):
    for k, c in IMMUNE_COLORS.items():
        if k.lower() in name.lower():
            return c
    return IMMUNE_COLORS["default"]


def get_peaks(bed_path, chrom, start, end):
    rows = []
    with gzip.open(bed_path, "rt") as f:
        for line in f:
            cols = line.rstrip().split("\t")
            if cols[0] != chrom: continue
            s, e = int(cols[1]), int(cols[2])
            if e < start or s > end: continue
            score = float(cols[4]) if cols[4] != '.' else 0.0
            summit = int(cols[6]) if len(cols) > 6 else (s + e) // 2
            rows.append({"start": s, "end": e, "name": cols[3],
                         "score": score, "summit": summit})
    return pd.DataFrame(rows) if rows else pd.DataFrame(
        columns=["start","end","name","score","summit"])


# ── Build figures for TLR cluster (chr4:38.7-38.9Mb) and individual genes ────
fig_w, fig_h = 16, 18
fig, axes = plt.subplots(len(GENES), 2, figsize=(fig_w, fig_h),
                          gridspec_kw={"hspace": 0.6, "wspace": 0.05})

genes_list = list(GENES.keys())

for row_idx, gene in enumerate(genes_list):
    g = GENES[gene]
    chrom, tss, strand = g["chrom"], g["tss"], g["strand"]
    gene_end = g["end"]
    gene_color = g["color"]

    reg_start = tss - WINDOW
    reg_end   = tss + WINDOW

    vdr_df  = get_peaks(VDR_BED,   chrom, reg_start, reg_end)
    gr_df   = get_peaks(NR3C1_BED, chrom, reg_start, reg_end)

    for col_idx, (df, tf, tf_color) in enumerate([
        (vdr_df,  "VDR",    "#E67E22"),
        (gr_df,   "NR3C1",  "#8E44AD"),
    ]):
        ax = axes[row_idx, col_idx]

        # Gene body
        gstart = min(tss, gene_end) - reg_start
        gwidth = abs(tss - gene_end)
        ax.barh(y=-0.5, width=gwidth, left=gstart, height=0.3,
                color=gene_color, alpha=0.4, zorder=2)
        # TSS arrow
        arrow_x = tss - reg_start
        arrowdir = 1 if strand == "+" else -1
        ax.annotate("", xy=(arrow_x + arrowdir * 400, -0.35),
                    xytext=(arrow_x, -0.35),
                    arrowprops=dict(arrowstyle="->", color=gene_color, lw=1.5))

        # Plot peaks as lollipops
        y_max = 0
        for _, pk in df.iterrows():
            summit_x = pk["summit"] - reg_start
            score = max(pk["score"], 0.05)
            cc = cell_color(pk["name"])
            ax.vlines(summit_x, 0, score, color=cc, linewidth=1.2, alpha=0.7, zorder=3)
            ax.plot(summit_x, score, "o", color=cc, markersize=4, alpha=0.9, zorder=4)
            y_max = max(y_max, score)

        # ±2kb promoter shading
        promo_w = 2000
        ax.axvspan(WINDOW - promo_w, WINDOW + promo_w, alpha=0.07, color="gray", zorder=1)
        ax.axvline(WINDOW, color="gray", lw=0.7, linestyle="--", alpha=0.5, zorder=1)

        # THP-1 + macrophage peak count
        if not df.empty:
            immune_n = sum(any(k in r["name"] for k in ["THP","macrophage","monocyte","GM12878","ALL","SUP-B15"])
                           for _, r in df.iterrows())
        else:
            immune_n = 0

        ax.set_xlim(0, 2 * WINDOW)
        ax.set_ylim(-0.8, max(y_max * 1.15, 1))
        ax.set_ylabel("Score", fontsize=7)
        if row_idx == 0:
            ax.set_title(f"{tf} ChIP-seq peaks\n(ReMap2022, ±15 kb from TSS)", fontsize=9)

        # Gene label
        ax.text(0.02, 0.97, f"{gene} ({chrom}:{tss:,})", transform=ax.transAxes,
                fontsize=7.5, va="top", ha="left", fontweight="bold", color=gene_color)
        ax.text(0.02, 0.82, f"n={len(df)} peaks total\n{immune_n} in immune/B cells",
                transform=ax.transAxes, fontsize=6.5, va="top", ha="left")

        # X axis: kb from TSS
        xticks = [0, WINDOW//2, WINDOW, WINDOW + WINDOW//2, 2*WINDOW]
        ax.set_xticks(xticks)
        ax.set_xticklabels([f"{(t - WINDOW)/1000:.0f} kb" for t in xticks], fontsize=6)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

# Legend
legend_elems = [
    Line2D([0],[0], marker="o", color="w", markerfacecolor="#E74C3C",  ms=6, label="THP-1"),
    Line2D([0],[0], marker="o", color="w", markerfacecolor="#C0392B",  ms=6, label="Macrophage"),
    Line2D([0],[0], marker="o", color="w", markerfacecolor="#E67E22",  ms=6, label="Monocyte"),
    Line2D([0],[0], marker="o", color="w", markerfacecolor="#8E44AD",  ms=6, label="GM12878/B-cell/ALL"),
    Line2D([0],[0], marker="o", color="w", markerfacecolor="#BDC3C7",  ms=6, label="Other"),
    mpatches.Patch(facecolor="gray", alpha=0.15, label="±2 kb promoter"),
]
fig.legend(handles=legend_elems, loc="lower center", ncol=3, fontsize=8,
           bbox_to_anchor=(0.5, -0.01), frameon=True)

fig.suptitle("VDR and GR/NR3C1 ChIP-seq peaks at TLR/A20/DUSP1 loci (hg38)\n"
             "ReMap2022 aggregated peaks — lollipop = peak summit (height = fold-change score)",
             fontsize=10, y=1.01)

plt.savefig(FIG_DIR / "locus_map_all_genes.pdf", bbox_inches="tight")
plt.savefig(FIG_DIR / "locus_map_all_genes.png", dpi=200, bbox_inches="tight")
plt.close()

print(f"Figure → {FIG_DIR/'locus_map_all_genes.pdf'}")


# ── TLR cluster zoomed: TLR10/TLR1/TLR6 on one axis (chr4:38.72-38.85Mb) ────
CLUSTER_START = 38720000
CLUSTER_END   = 38870000
cluster_genes = {
    "TLR10": {"tss": 38767648, "end": 38739275, "strand": "-", "color": "#E74C3C"},
    "TLR1":  {"tss": 38770059, "end": 38795755, "strand": "+", "color": "#3498DB"},
    "TLR6":  {"tss": 38806672, "end": 38831036, "strand": "+", "color": "#2ECC71"},
}

fig2, axes2 = plt.subplots(2, 1, figsize=(14, 7),
                            gridspec_kw={"hspace": 0.5})
for ax2, (bed, tf, color) in zip(axes2, [
    (VDR_BED,   "VDR (+VitD/1,25-OH-2D3)",  "#E67E22"),
    (NR3C1_BED, "GR/NR3C1 (+Dex)",           "#8E44AD"),
]):
    df = get_peaks(bed, "chr4", CLUSTER_START, CLUSTER_END)
    y_max = 0
    for _, pk in df.iterrows():
        summit_x = pk["summit"]
        sc = max(pk["score"], 0.05)
        cc = cell_color(pk["name"])
        ax2.vlines(summit_x, 0, sc, color=cc, lw=1.3, alpha=0.75)
        ax2.plot(summit_x, sc, "o", color=cc, ms=5, alpha=0.9)
        y_max = max(y_max, sc)

    # Gene bodies
    for gname, gi in cluster_genes.items():
        gstart = min(gi["tss"], gi["end"])
        gend   = max(gi["tss"], gi["end"])
        yp = -2.0 - list(cluster_genes.keys()).index(gname) * 0.8
        ax2.barh(y=yp, width=gend-gstart, left=gstart, height=0.5,
                 color=gi["color"], alpha=0.5)
        ax2.annotate(gname, xy=((gstart+gend)/2, yp + 0.25), fontsize=8,
                     ha="center", va="bottom", color=gi["color"], fontweight="bold")
        ax2.axvspan(gi["tss"]-2000, gi["tss"]+2000, alpha=0.06, color=gi["color"])

    ax2.set_xlim(CLUSTER_START, CLUSTER_END)
    ax2.set_ylim(-5, max(y_max*1.1, 5))
    ax2.set_title(f"{tf} ChIP-seq — chr4:38.72-38.85 Mb (TLR10/TLR1/TLR6 cluster)", fontsize=10)
    ax2.set_ylabel("Score")
    xt = np.arange(CLUSTER_START, CLUSTER_END+1, 20000)
    ax2.set_xticks(xt)
    ax2.set_xticklabels([f"{x/1e6:.3f} Mb" for x in xt], fontsize=7, rotation=30, ha="right")
    ax2.spines["top"].set_visible(False)
    ax2.spines["right"].set_visible(False)

fig2.legend(handles=legend_elems, loc="lower center", ncol=3, fontsize=8,
            bbox_to_anchor=(0.5, -0.02), frameon=True)
fig2.suptitle("VDR vs GR/NR3C1 ChIP-seq: TLR10/TLR1/TLR6 cluster (chr4:38Mb, hg38)\n"
              "Shaded: ±2 kb from TSS (promoter region)", fontsize=10)
plt.savefig(FIG_DIR / "tlr_cluster_chr4_pileup.pdf", bbox_inches="tight")
plt.savefig(FIG_DIR / "tlr_cluster_chr4_pileup.png", dpi=200, bbox_inches="tight")
plt.close()
print(f"Figure → {FIG_DIR/'tlr_cluster_chr4_pileup.pdf'}")

print("\nAll figures done.")
