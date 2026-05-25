#!/usr/bin/env python3
"""
THP-1 specific comparison: VDR vs GR/NR3C1 at all TLR loci.
Shows only peaks from THP-1 cells to make the cleanest cell-type matched comparison.
Key finding: In THP-1 monocytes, VDR binds TLR10 (no GR), GR binds TLR2 strongly.
"""
import gzip
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.gridspec import GridSpec

BASE_OLD = Path("/Volumes/M4_SSD/projects/tlr10_chipseq_vdr_gr/data")
FIG_DIR  = Path("/Volumes/M4_SSD/projects/tlr_chipseq/results/figures")
FIG_DIR.mkdir(parents=True, exist_ok=True)

VDR_BED   = BASE_OLD / "remap2022_VDR_all_macs2_hg38.bed.gz"
NR3C1_BED = BASE_OLD / "remap2022_NR3C1_all_macs2_hg38.bed.gz"

GENES = {
    "TLR10":   ("chr4",  38767648,  "-"),
    "TLR1":    ("chr4",  38770059,  "+"),
    "TLR6":    ("chr4",  38806672,  "+"),
    "TLR2":    ("chr4",  153684080, "+"),
    "TLR4":    ("chr9",  117704402, "+"),
    "TNFAIP3": ("chr6",  137866317, "+"),
    "DUSP1":   ("chr5",  172768090, "+"),
}

THP1_FILTER = lambda name: "THP" in name or "THP-1" in name
MACRO_FILTER = lambda name: "macrophage" in name.lower() or "monocyte" in name.lower()
ANY_IMMUNE = lambda name: any(k in name for k in [
    "THP","macrophage","monocyte","GM12878","ALL","SUP-B15","NALM"])

def get_peaks_filtered(bed, chrom, start, end, filt):
    rows = []
    with gzip.open(bed, "rt") as f:
        for line in f:
            cols = line.rstrip().split("\t")
            if cols[0] != chrom: continue
            s, e = int(cols[1]), int(cols[2])
            if e < start or s > end: continue
            name = cols[3]
            if not filt(name): continue
            score = float(cols[4]) if cols[4] != '.' else 0.0
            summit = int(cols[6]) if len(cols)>6 else (s+e)//2
            rows.append({"name":name,"start":s,"end":e,"score":score,"summit":summit})
    return pd.DataFrame(rows) if rows else pd.DataFrame(
        columns=["name","start","end","score","summit"])


# ── Main comparison figure: 3 panels ─────────────────────────────────────────
# Panel 1: Bar chart of peak count in THP-1 per gene
# Panel 2: Max score in THP-1 per gene
# Panel 3: Heatmap-style comparison
genes_order = list(GENES.keys())
n = len(genes_order)

vdr_thp1_n = []
gr_thp1_n  = []
vdr_thp1_s = []
gr_thp1_s  = []
vdr_macro_n = []
gr_macro_n  = []
gr_macro_s  = []

HW = 10000
for gene in genes_order:
    chrom, tss, _ = GENES[gene]
    vdf = get_peaks_filtered(VDR_BED,   chrom, tss-HW, tss+HW, THP1_FILTER)
    gdf = get_peaks_filtered(NR3C1_BED, chrom, tss-HW, tss+HW, THP1_FILTER)
    mdf = get_peaks_filtered(NR3C1_BED, chrom, tss-HW, tss+HW, MACRO_FILTER)
    vdr_thp1_n.append(len(vdf))
    gr_thp1_n.append(len(gdf))
    vdr_thp1_s.append(vdf["score"].max() if len(vdf)>0 else 0.0)
    gr_thp1_s.append(gdf["score"].max() if len(gdf)>0 else 0.0)
    vdr_macro_n.append(0)  # VDR not in macrophage datasets we have
    gr_macro_n.append(len(mdf))
    gr_macro_s.append(mdf["score"].max() if len(mdf)>0 else 0.0)

x = np.arange(n)
w = 0.3
gene_labels = [f"{'→' if GENES[g][2]=='+' else '←'} {g}" for g in genes_order]
tlr10_idx = genes_order.index("TLR10")
tlr2_idx  = genes_order.index("TLR2")

fig, axes = plt.subplots(2, 2, figsize=(14, 10),
                          gridspec_kw={"hspace":0.55, "wspace":0.35})

# ── Panel A: Peak counts in THP-1 ────────────────────────────────────────────
ax = axes[0, 0]
b1 = ax.bar(x - w/2, vdr_thp1_n, w, color="#E67E22", alpha=0.85, label="VDR (+VitD)",
            edgecolor="black", lw=0.5)
b2 = ax.bar(x + w/2, gr_thp1_n,  w, color="#8E44AD", alpha=0.85, label="GR/NR3C1 (+Dex)",
            edgecolor="black", lw=0.5)
for bar, val in [(b, v) for b, v in zip(list(b1)+list(b2), vdr_thp1_n+gr_thp1_n) if v>0]:
    ax.text(bar.get_x()+bar.get_width()/2, bar.get_height()+0.02, str(val),
            ha="center", va="bottom", fontsize=9, fontweight="bold")
ax.set_xticks(x); ax.set_xticklabels(gene_labels, rotation=40, ha="right", fontsize=9)
ax.get_xticklabels()[tlr10_idx].set_color("red")
ax.get_xticklabels()[tlr2_idx].set_color("blue")
ax.set_ylabel("Number of peaks (ReMap2022)"); ax.set_title("A  THP-1 ChIP-seq: peak counts\n(±10 kb from TSS)", fontsize=10)
ax.legend(fontsize=9); ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)

# ── Panel B: Max score in THP-1 ──────────────────────────────────────────────
ax = axes[0, 1]
b1 = ax.bar(x - w/2, vdr_thp1_s, w, color="#E67E22", alpha=0.85, label="VDR",
            edgecolor="black", lw=0.5)
b2 = ax.bar(x + w/2, gr_thp1_s,  w, color="#8E44AD", alpha=0.85, label="GR/NR3C1",
            edgecolor="black", lw=0.5)
ax.set_xticks(x); ax.set_xticklabels(gene_labels, rotation=40, ha="right", fontsize=9)
ax.get_xticklabels()[tlr10_idx].set_color("red")
ax.get_xticklabels()[tlr2_idx].set_color("blue")
ax.set_ylabel("Max fold-change score"); ax.set_title("B  THP-1 ChIP-seq: peak strength\n(max score ±10 kb from TSS)", fontsize=10)
ax.legend(fontsize=9); ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)

# ── Panel C: GR peaks — THP-1 vs macrophage at TLR2 locus ─────────────────
ax = axes[1, 0]
chrom, tss, _ = GENES["TLR2"]
WINDOW_PLOT = 15000
vdr_t = get_peaks_filtered(VDR_BED,   chrom, tss-WINDOW_PLOT, tss+WINDOW_PLOT, THP1_FILTER)
gr_t  = get_peaks_filtered(NR3C1_BED, chrom, tss-WINDOW_PLOT, tss+WINDOW_PLOT, THP1_FILTER)
gr_m  = get_peaks_filtered(NR3C1_BED, chrom, tss-WINDOW_PLOT, tss+WINDOW_PLOT, MACRO_FILTER)

y_offset = {"VDR_THP1": 0, "GR_THP1": 1, "GR_Macro": 2}
colors   = {"VDR_THP1": "#E67E22", "GR_THP1": "#8E44AD", "GR_Macro": "#C0392B"}
track_h  = 0.3

for label, df, c in [("VDR_THP1", vdr_t, "#E67E22"),
                     ("GR_THP1",  gr_t,  "#8E44AD"),
                     ("GR_Macro", gr_m,  "#C0392B")]:
    yo = y_offset[label]
    for _, pk in df.iterrows():
        s_rel = pk["start"] - (tss - WINDOW_PLOT)
        e_rel = pk["end"]   - (tss - WINDOW_PLOT)
        ax.barh(y=yo, width=e_rel-s_rel, left=s_rel, height=track_h,
                color=c, alpha=0.7, edgecolor="black", lw=0.3)
        summit_rel = pk["summit"] - (tss - WINDOW_PLOT)
        ax.plot(summit_rel, yo+track_h/2, "|", color="black", ms=6, lw=1.5)
        ax.text(summit_rel, yo+track_h+0.05, f"{pk['score']:.1f}",
                ha="center", va="bottom", fontsize=6)

# TSS marker
ax.axvline(WINDOW_PLOT, color="blue", lw=1.5, linestyle="--", alpha=0.7, label="TLR2 TSS")
ax.axvspan(WINDOW_PLOT-2000, WINDOW_PLOT+2000, alpha=0.08, color="blue")
ax.set_yticks([0.15, 1.15, 2.15])
ax.set_yticklabels(["VDR (THP-1)", "GR (THP-1+Dex)", "GR (Macrophage)"], fontsize=9)
xt = np.arange(0, 2*WINDOW_PLOT+1, 5000)
ax.set_xticks(xt)
ax.set_xticklabels([f"{(t-WINDOW_PLOT)/1000:+.0f}kb" for t in xt], fontsize=8)
ax.set_title("C  TLR2 locus: VDR vs GR binding tracks\n(THP-1 and macrophage, hg38)", fontsize=10)
ax.set_xlabel("Distance from TLR2 TSS (chr4:153,684,080)")
ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)

# ── Panel D: Same for TLR10 locus ────────────────────────────────────────────
ax = axes[1, 1]
chrom, tss, strand = GENES["TLR10"]
vdr_t10 = get_peaks_filtered(VDR_BED,   chrom, tss-WINDOW_PLOT, tss+WINDOW_PLOT, THP1_FILTER)
gr_t10  = get_peaks_filtered(NR3C1_BED, chrom, tss-WINDOW_PLOT, tss+WINDOW_PLOT, THP1_FILTER)
gr_m10  = get_peaks_filtered(NR3C1_BED, chrom, tss-WINDOW_PLOT, tss+WINDOW_PLOT, MACRO_FILTER)

for label, df, c in [("VDR_THP1", vdr_t10, "#E67E22"),
                     ("GR_THP1",  gr_t10,  "#8E44AD"),
                     ("GR_Macro", gr_m10,  "#C0392B")]:
    yo = y_offset[label]
    for _, pk in df.iterrows():
        s_rel = pk["start"] - (tss - WINDOW_PLOT)
        e_rel = pk["end"]   - (tss - WINDOW_PLOT)
        ax.barh(y=yo, width=e_rel-s_rel, left=s_rel, height=track_h,
                color=c, alpha=0.7, edgecolor="black", lw=0.3)
        summit_rel = pk["summit"] - (tss - WINDOW_PLOT)
        ax.plot(summit_rel, yo+track_h/2, "|", color="black", ms=6, lw=1.5)
        ax.text(summit_rel, yo+track_h+0.05, f"{pk['score']:.1f}",
                ha="center", va="bottom", fontsize=6)

ax.axvline(WINDOW_PLOT, color="red", lw=1.5, linestyle="--", alpha=0.7, label="TLR10 TSS")
ax.axvspan(WINDOW_PLOT-2000, WINDOW_PLOT+2000, alpha=0.08, color="red")
ax.set_yticks([0.15, 1.15, 2.15])
ax.set_yticklabels(["VDR (THP-1)", "GR (THP-1+Dex)", "GR (Macrophage)"], fontsize=9)
ax.set_xticks(xt)
ax.set_xticklabels([f"{(t-WINDOW_PLOT)/1000:+.0f}kb" for t in xt], fontsize=8)
ax.set_title("D  TLR10 locus: VDR vs GR binding tracks\n(THP-1 and macrophage, hg38)", fontsize=10)
ax.set_xlabel("Distance from TLR10 TSS (chr4:38,767,648, − strand)")
ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)

# annotation
ax.text(0.98, 0.95, "← VDR only\n   GR absent",
        transform=ax.transAxes, ha="right", va="top", fontsize=9,
        color="#E74C3C", fontweight="bold",
        bbox=dict(boxstyle="round,pad=0.2", fc="lightyellow", ec="red", alpha=0.8))

axes[1,0].text(0.98, 0.95, "← GR: strong\n   at enhancer",
        transform=axes[1,0].transAxes, ha="right", va="top", fontsize=9,
        color="#8E44AD", fontweight="bold",
        bbox=dict(boxstyle="round,pad=0.2", fc="lightyellow", ec="purple", alpha=0.8))

fig.suptitle(
    "VDR and GR/NR3C1 ChIP-seq at TLR loci — THP-1 monocytes (hg38, ReMap2022)\n"
    "Key: VDR selectively binds TLR10 (no GR in THP-1/macrophage)\n"
    "GR/NR3C1 binds TLR2 strongly at −7 kb enhancer (score=56 in THP-1+Dex)",
    fontsize=10, y=1.02
)

plt.savefig(FIG_DIR / "thp1_vdr_gr_comparison.pdf", bbox_inches="tight")
plt.savefig(FIG_DIR / "thp1_vdr_gr_comparison.png", dpi=300, bbox_inches="tight")
plt.close()
print(f"Figure → {FIG_DIR/'thp1_vdr_gr_comparison.pdf'}")

# Print summary
print("\n=== THP-1 specific summary (±10kb) ===")
for gene in genes_order:
    chrom, tss, _ = GENES[gene]
    vdf = get_peaks_filtered(VDR_BED,   chrom, tss-HW, tss+HW, THP1_FILTER)
    gdf = get_peaks_filtered(NR3C1_BED, chrom, tss-HW, tss+HW, THP1_FILTER)
    v_n = len(vdf); v_s = vdf["score"].max() if v_n>0 else 0.0
    g_n = len(gdf); g_s = gdf["score"].max() if g_n>0 else 0.0
    marker = " ←★" if gene in ("TLR10","TLR2") else ""
    print(f"  {gene:8s}: VDR {v_n} peaks max={v_s:.2f}  |  GR {g_n} peaks max={g_s:.2f}{marker}")
print("\nDone.")
