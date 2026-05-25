#!/usr/bin/env python3
"""
Extended ReMap2022 analysis: VDR + NR3C1 peaks at all TLR loci + A20 + DUSP1 (hg38)
"""
import gzip, sys, warnings
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

warnings.filterwarnings("ignore")

BASE_OLD  = Path("/Volumes/M4_SSD/projects/tlr10_chipseq_vdr_gr/data")
BASE_NEW  = Path("/Volumes/M4_SSD/projects/tlr_chipseq")
RES_DIR   = BASE_NEW / "results" / "peaks"
FIG_DIR   = BASE_NEW / "results" / "figures"
RES_DIR.mkdir(parents=True, exist_ok=True)
FIG_DIR.mkdir(parents=True, exist_ok=True)

VDR_BED   = BASE_OLD / "remap2022_VDR_all_macs2_hg38.bed.gz"
NR3C1_BED = BASE_OLD / "remap2022_NR3C1_all_macs2_hg38.bed.gz"

# hg38 TSS coordinates (strand-aware)
GENES_HG38 = {
    "TLR10":  ("chr4",  38767648,  "-"),  # TSS = end of gene body
    "TLR1":   ("chr4",  38770059,  "+"),
    "TLR6":   ("chr4",  38806672,  "+"),
    "TLR2":   ("chr4",  153684080, "+"),
    "TLR4":   ("chr9",  117704402, "+"),
    "TNFAIP3":("chr6",  137866317, "+"),  # A20
    "DUSP1":  ("chr5",  172768090, "+"),  # MKP-1
}
PROMOTER_HW = 5000  # ±5 kb for peaks


def parse_bed_gz_region(bed_path: Path, chrom: str, start: int, end: int) -> pd.DataFrame:
    records = []
    with gzip.open(bed_path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 3:
                continue
            c, s, e = cols[0], int(cols[1]), int(cols[2])
            if c != chrom:
                continue
            if e < start or s > end:
                continue
            name   = cols[3] if len(cols) > 3 else "."
            signal = float(cols[6]) if len(cols) > 6 else 0.0
            summit = int(cols[9]) if len(cols) > 9 else (s + e) // 2
            records.append({
                "chrom": c, "start": s, "end": e, "name": name,
                "signal": signal, "summit_abs": s + summit,
            })
    return pd.DataFrame(records) if records else pd.DataFrame(
        columns=["chrom", "start", "end", "name", "signal", "summit_abs"]
    )


def analyse_tf_at_loci(tf_name: str, bed_path: Path, genes: dict, hw: int) -> pd.DataFrame:
    rows = []
    for gene, (chrom, tss, strand) in genes.items():
        p_start = tss - hw
        p_end   = tss + hw
        df = parse_bed_gz_region(bed_path, chrom, p_start, p_end)
        n_peaks = len(df)
        max_sig = df["signal"].max() if n_peaks > 0 else 0.0
        studies = "; ".join(sorted(set(
            r.split(".")[0] for r in df["name"].tolist()
        ))) if n_peaks > 0 else ""
        rows.append({
            "TF": tf_name, "gene": gene, "chrom": chrom, "tss": tss, "strand": strand,
            "promoter_start": p_start, "promoter_end": p_end,
            "n_peaks": n_peaks, "max_signal": round(max_sig, 2),
            "studies": studies,
        })
    return pd.DataFrame(rows)


print("=== VDR peaks at TLR/A20/DUSP1 loci (ReMap2022, hg38) ===")
vdr_df  = analyse_tf_at_loci("VDR",   VDR_BED,   GENES_HG38, PROMOTER_HW)
gr_df   = analyse_tf_at_loci("NR3C1", NR3C1_BED, GENES_HG38, PROMOTER_HW)

combined = pd.merge(
    vdr_df[["gene","chrom","tss","strand","n_peaks","max_signal","studies"]],
    gr_df[["gene","n_peaks","max_signal","studies"]],
    on="gene", suffixes=("_VDR","_NR3C1")
)
combined.to_csv(RES_DIR / "remap_all_loci_summary.csv", index=False)
print(combined.to_string())

# ── Figure: bar chart of peak counts at each locus ──────────────────────────
fig, axes = plt.subplots(1, 2, figsize=(12, 5), sharey=False)
genes_order = list(GENES_HG38.keys())
x = np.arange(len(genes_order))
width = 0.35

for ax, (col, color, label) in zip(axes, [
    ("n_peaks_VDR",   "#E67E22", "VDR (+VitD)"),
    ("n_peaks_NR3C1", "#8E44AD", "NR3C1/GR (+Dex)"),
]):
    vals = [combined.loc[combined["gene"] == g, col].values[0] for g in genes_order]
    bars = ax.bar(x, vals, color=color, alpha=0.8, edgecolor="black", linewidth=0.5)
    for bar, v in zip(bars, vals):
        if v > 0:
            ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.05,
                    str(v), ha="center", va="bottom", fontsize=9, fontweight="bold")
    ax.set_xticks(x)
    ax.set_xticklabels(genes_order, rotation=45, ha="right", fontsize=10)
    ax.set_ylabel("Number of ChIP-seq peaks (ReMap2022)")
    ax.set_title(f"{label}\nReMap2022 aggregated peaks (±5 kb from TSS)")
    ax.set_ylim(0, max(vals) + 2 if max(vals) > 0 else 3)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

# Highlight TLR10 and TLR2
for ax in axes:
    for i, g in enumerate(genes_order):
        if g in ("TLR10", "TLR2"):
            ax.get_xticklabels()[i].set_color("red" if g == "TLR10" else "blue")
            ax.get_xticklabels()[i].set_fontweight("bold")

plt.tight_layout()
out_pdf = FIG_DIR / "remap_peak_counts_all_loci.pdf"
out_png = FIG_DIR / "remap_peak_counts_all_loci.png"
plt.savefig(out_pdf, bbox_inches="tight")
plt.savefig(out_png, dpi=300, bbox_inches="tight")
plt.close()
print(f"\nFigure → {out_pdf}")

# ── Detailed peak table ───────────────────────────────────────────────────────
print("\n=== Detailed peak calls ===")
for tf_name, df in [("VDR", vdr_df), ("NR3C1", gr_df)]:
    print(f"\n{tf_name}:")
    sub = df[df["n_peaks"] > 0][["gene","n_peaks","max_signal","studies"]]
    if sub.empty:
        print("  (no peaks at any locus)")
    else:
        print(sub.to_string(index=False))

print("\n=== DONE: remap_all_loci_summary.csv ===")
