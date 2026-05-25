#!/usr/bin/env python3
"""
Detailed ReMap2022 peak analysis with correct score parsing and cell-type filtering.
Key insight from ReMap2022 BED9 format:
  col4 = score (fold change / -log10 pval depending on dataset)
  col6 = summit coordinate (NOT signal value)
"""
import gzip, re
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

BASE_OLD  = Path("/Volumes/M4_SSD/projects/tlr10_chipseq_vdr_gr/data")
BASE_NEW  = Path("/Volumes/M4_SSD/projects/tlr_chipseq")
RES_DIR   = BASE_NEW / "results" / "peaks"
FIG_DIR   = BASE_NEW / "results" / "figures"
for d in (RES_DIR, FIG_DIR): d.mkdir(parents=True, exist_ok=True)

VDR_BED   = BASE_OLD / "remap2022_VDR_all_macs2_hg38.bed.gz"
NR3C1_BED = BASE_OLD / "remap2022_NR3C1_all_macs2_hg38.bed.gz"

# hg38 TSS (strand-aware)
GENES = {
    "TLR10":   ("chr4",  38767648,  "-"),
    "TLR1":    ("chr4",  38770059,  "+"),
    "TLR6":    ("chr4",  38806672,  "+"),
    "TLR2":    ("chr4",  153684080, "+"),
    "TLR4":    ("chr9",  117704402, "+"),
    "TNFAIP3": ("chr6",  137866317, "+"),
    "DUSP1":   ("chr5",  172768090, "+"),
    "NOD2":    ("chr16", 50756049,  "-"),  # hg38, NM_022162
    "VDR":     ("chr12", 47844418,  "+"),  # VDR自身のプロモーター（自己調節確認用）
}

# Cell types relevant to innate immunity / the paper
IMMUNE_KEYWORDS = [
    "THP-1", "THP1", "macrophage", "monocyte", "BMDM",
    "GM12878",  # B lymphocyte - especially relevant for TLR10
    "ALL",      # Acute lymphoblastic leukemia (B-cell)
    "SUP-B15",  # B-cell leukemia
    "NALM",     # B-cell leukemia
]

def is_immune(name: str) -> bool:
    return any(kw.lower() in name.lower() for kw in IMMUNE_KEYWORDS)


def parse_region(bed_path: Path, chrom: str, start: int, end: int) -> pd.DataFrame:
    records = []
    with gzip.open(bed_path, "rt") as fh:
        for line in fh:
            if line.startswith("#"): continue
            c = line[:6]
            if not c.startswith(chrom[:4]): continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 5: continue
            if cols[0] != chrom: continue
            s, e = int(cols[1]), int(cols[2])
            if e < start or s > end: continue
            name   = cols[3]
            try:
                score = float(cols[4])
            except ValueError:
                score = 0.0
            summit = int(cols[6]) if len(cols) > 6 else (s + e) // 2
            # parse GSE/ENCODE id from name
            m = re.match(r"([^.]+)\.[^.]+\.(.+)", name)
            gse = m.group(1) if m else name
            celltype = m.group(2) if m else ""
            records.append({
                "chrom": chrom, "start": s, "end": e,
                "name": name, "score": score, "summit": summit,
                "gse": gse, "celltype": celltype,
                "immune": is_immune(name),
                "dist_to_win_start": s - start,
            })
    return pd.DataFrame(records) if records else pd.DataFrame(
        columns=["chrom","start","end","name","score","summit","gse","celltype","immune","dist_to_win_start"])


def summarise_tf(tf: str, bed: Path, hw: int = 5000) -> pd.DataFrame:
    rows = []
    for gene, (chrom, tss, strand) in GENES.items():
        df = parse_region(bed, chrom, tss - hw, tss + hw)
        n_all    = len(df)
        n_immune = len(df[df["immune"]]) if not df.empty else 0
        max_score_all    = df["score"].max() if n_all > 0 else 0.0
        max_score_immune = df.loc[df["immune"],"score"].max() if n_immune > 0 else 0.0
        imm_studies = "; ".join(sorted(set(
            df.loc[df["immune"],"gse"].tolist()))) if n_immune > 0 else ""
        rows.append({
            "TF": tf, "gene": gene,
            "n_peaks_total": n_all, "n_peaks_immune": n_immune,
            "max_score_total": round(max_score_all, 2),
            "max_score_immune": round(max_score_immune, 2),
            "immune_studies": imm_studies,
        })
    return pd.DataFrame(rows)


print("Running detailed peak analysis …")
vdr_sum  = summarise_tf("VDR",   VDR_BED)
gr_sum   = summarise_tf("NR3C1", NR3C1_BED)

# Merge
merged = pd.merge(
    vdr_sum[["gene","n_peaks_total","n_peaks_immune","max_score_total","max_score_immune","immune_studies"]],
    gr_sum [["gene","n_peaks_total","n_peaks_immune","max_score_total","max_score_immune","immune_studies"]],
    on="gene", suffixes=("_VDR","_GR")
)
merged.to_csv(RES_DIR / "peak_summary_detailed.csv", index=False)

print("\n=== PEAK SUMMARY (±5 kb from TSS) ===")
print(merged[["gene","n_peaks_immune_VDR","max_score_immune_VDR",
              "n_peaks_immune_GR","max_score_immune_GR"]].to_string(index=False))

# ── Figure: grouped bar chart (immune-cell peaks only) ───────────────────────
genes_order = list(GENES.keys())
n = len(genes_order)
x = np.arange(n)
w = 0.35

fig, axes = plt.subplots(1, 2, figsize=(14, 5))

for ax, (col_n, col_s, tf, color, title) in zip(axes, [
    ("n_peaks_immune_VDR",  "max_score_immune_VDR",  "VDR",   "#E67E22",
     "VDR ChIP-seq peaks\n(immune/B-cell datasets, ±5 kb from TSS)"),
    ("n_peaks_immune_GR",   "max_score_immune_GR",   "NR3C1", "#8E44AD",
     "GR/NR3C1 ChIP-seq peaks\n(immune/B-cell datasets, ±5 kb from TSS)"),
]):
    counts = [merged.loc[merged["gene"]==g, col_n].values[0] for g in genes_order]
    scores = [merged.loc[merged["gene"]==g, col_s].values[0] for g in genes_order]
    bars = ax.bar(x, counts, color=color, alpha=0.8, edgecolor="black", lw=0.7)
    for bar, cnt, sc in zip(bars, counts, scores):
        if cnt > 0:
            ax.text(bar.get_x()+bar.get_width()/2, bar.get_height()+0.05,
                    f"{cnt}", ha="center", va="bottom", fontsize=9, fontweight="bold")
    ax2 = ax.twinx()
    ax2.plot(x, scores, "o--", color="darkred", alpha=0.6, ms=6, label="max score")
    ax2.set_ylabel("Max ChIP score (fold change)", color="darkred", fontsize=9)
    ax2.tick_params(axis="y", labelcolor="darkred")
    ax.set_xticks(x)
    lbls = ax.set_xticklabels(genes_order, rotation=45, ha="right", fontsize=10)
    for lbl, g in zip(lbls, genes_order):
        if g == "TLR10": lbl.set_color("red"); lbl.set_fontweight("bold")
        if g == "TLR2":  lbl.set_color("blue"); lbl.set_fontweight("bold")
    ax.set_ylabel("Number of peaks")
    ax.set_title(title, fontsize=10)
    ax.spines["top"].set_visible(False)

plt.tight_layout()
plt.savefig(FIG_DIR / "peak_summary_immune_only.pdf", bbox_inches="tight")
plt.savefig(FIG_DIR / "peak_summary_immune_only.png", dpi=300, bbox_inches="tight")
plt.close()

# ── THP-1 specific table ──────────────────────────────────────────────────────
print("\n=== THP-1 SPECIFIC PEAKS ===")
thp1_rows = []
for gene, (chrom, tss, strand) in GENES.items():
    vdf = parse_region(VDR_BED,   chrom, tss-5000, tss+5000)
    gdf = parse_region(NR3C1_BED, chrom, tss-5000, tss+5000)
    vdr_thp1 = vdf[vdf["celltype"].str.contains("THP", case=False, na=False)]
    gr_thp1  = gdf[gdf["celltype"].str.contains("THP", case=False, na=False)]
    thp1_rows.append({
        "gene": gene,
        "VDR_THP1_n": len(vdr_thp1),
        "VDR_THP1_maxscore": round(vdr_thp1["score"].max(), 2) if len(vdr_thp1) else 0,
        "VDR_THP1_peaks": "; ".join(vdr_thp1["name"].tolist()),
        "GR_THP1_n": len(gr_thp1),
        "GR_THP1_maxscore": round(gr_thp1["score"].max(), 2) if len(gr_thp1) else 0,
        "GR_THP1_peaks": "; ".join(gr_thp1["name"].tolist()),
    })
thp1_df = pd.DataFrame(thp1_rows)
thp1_df.to_csv(RES_DIR / "thp1_specific_peaks.csv", index=False)
print(thp1_df[["gene","VDR_THP1_n","VDR_THP1_maxscore","GR_THP1_n","GR_THP1_maxscore"]].to_string(index=False))

print("\n=== MACROPHAGE/MONOCYTE SPECIFIC NR3C1 AT TLR2 ===")
gdf = parse_region(NR3C1_BED, "chr4", 153679080, 153689080)
mac = gdf[gdf["celltype"].str.contains("macrophage|monocyte|THP", case=False, na=False)]
print(mac[["name","score","summit"]].to_string(index=False))

print(f"\nFigure → {FIG_DIR/'peak_summary_immune_only.pdf'}")
print("DONE.")
