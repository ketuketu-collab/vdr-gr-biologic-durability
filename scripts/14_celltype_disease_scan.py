#!/usr/bin/env python3
"""
Cell-type specific VDR/GR score by disease area
LS180 (colon epithelial) → IBD/UC
LX2   (hepatic stellate)  → ILD/liver fibrosis
LCLGM10861 (B-cell-like) → autoimmune/MS
kidney-cortex             → SLE/nephritis
GM12878 (B cell, GR only) → autoimmune/MS
macrophage (GR only)      → RA/IBD
BEAS-2B (GR only)         → asthma
"""
import gzip, os
import pandas as pd
from collections import defaultdict

VDR_BED = "/Volumes/M4_SSD/projects/tlr10_chipseq_vdr_gr/data/remap2022_VDR_all_macs2_hg38.bed.gz"
GR_BED  = "/Volumes/M4_SSD/projects/tlr10_chipseq_vdr_gr/data/remap2022_NR3C1_all_macs2_hg38.bed.gz"
OUT_DIR = "/Volumes/M4_SSD/projects/tlr_chipseq/results"

# ── 細胞種定義 ──────────────────────────────────────────────────────────
CELL_TYPES = {
    # (BED file, keyword, label)
    "VDR_THP1":       (VDR_BED, "THP-1",          "VDR"),
    "VDR_LS180":      (VDR_BED, "LS180",           "VDR"),
    "VDR_LX2":        (VDR_BED, "LX2",             "VDR"),
    "VDR_LCL":        (VDR_BED, "LCLGM10861",      "VDR"),
    "VDR_kidney":     (VDR_BED, "kidney-cortex",   "VDR"),
    "GR_THP1":        (GR_BED,  "THP-1",           "GR"),
    "GR_macrophage":  (GR_BED,  "macrophage",      "GR"),
    "GR_GM12878":     (GR_BED,  "GM12878",         "GR"),
    "GR_BEAS2B":      (GR_BED,  "BEAS-2B",         "GR"),
    "GR_IMR90":       (GR_BED,  "IMR-90",          "GR"),
}

# ── 疾患別遺伝子リスト ───────────────────────────────────────────────────
DISEASE_GENES = {
    "IBD/UC": {
        "IL23A":  ("chr5",  40390734,  40410734),
        "TNF":    ("chr6",  31575567,  31595567),
        "TYK2":   ("chr19", 10351457,  10371457),
        "ITGB7":  ("chr12", 54699526,  54719526),
        "NOD2":   ("chr16", 50713000,  50733000),
        "IL10":   ("chr1",  206941650, 206961650),
        "MADCAM1":("chr19", 5786452,   5806452),
        "JAK1":   ("chr1",  64833531,  64853531),
        "TLR10":  ("chr4",  38757000,  38777000),
        "STAT3":  ("chr17", 42313327,  42333327),
    },
    "RA": {
        "TNF":    ("chr6",  31575567,  31595567),
        "IL6R":   ("chr1",  154385583, 154405583),
        "CD86":   ("chr3",  122055493, 122075493),
        "MS4A1":  ("chr11", 60455330,  60475330),
        "IL1B":   ("chr2",  112829751, 112849751),
        "PTPRC":  ("chr1",  198602162, 198622162),
    },
    "ILD/fibrosis": {
        "VEGFA":  ("chr6",  43770209,  43790209),
        "PDGFRB": ("chr5",  150110823, 150130823),
        "TGFB1":  ("chr19", 41836941,  41856941),
        "FGFR1":  ("chr8",  38395133,  38415133),
        "IL6R":   ("chr1",  154385583, 154405583),
    },
    "SLE/nephritis": {
        "CD38":   ("chr4",  15779893,  15799893),
        "IFNAR1": ("chr21", 34694908,  34714908),
        "IL6R":   ("chr1",  154385583, 154405583),
        "C5AR1":  ("chr19", 55066012,  55086012),
        "TNFSF13B":("chr13",108472287, 108492287),
    },
    "asthma/rhinitis": {
        "IL13":   ("chr5",  132673986, 132693986),
        "TSLP":   ("chr5",  110401135, 110421135),
        "IL4R":   ("chr16", 27338042,  27358042),
        "IL5":    ("chr5",  132563299, 132583299),
        "IL33":   ("chr9",  6215786,   6235786),
        "CCR3":   ("chr3",  46207342,  46227342),
        "SIGLEC8":("chr19", 51874012,  51894012),
    },
    "MS/autoimmune": {
        "IL23A":  ("chr5",  40390734,  40410734),
        "CD52":   ("chr1",  26686832,  26706832),
        "ITGA4":  ("chr2",  182388872, 182408872),
        "ITGAL":  ("chr16", 30480897,  30500897),
        "LAG3":   ("chr12", 6772090,   6792090),
        "PDCD1":  ("chr2",  241849694, 241869694),
        "MS4A1":  ("chr11", 60455330,  60475330),
        "CD19":   ("chr16", 28932905,  28952905),
    },
    "irAE/checkpoint": {
        "CD274":  ("chr9",  5450503,   5470503),
        "PDCD1":  ("chr2",  241849694, 241869694),
        "LAG3":   ("chr12", 6772090,   6792090),
        "CTLA4":  ("chr2",  203867771, 203887771),
        "HAVCR2": ("chr5",  156631630, 156651630),
    },
}

# ── BED読み込み（細胞種フィルタ付き）───────────────────────────────────
def load_bed(bed_gz, keyword):
    peaks = defaultdict(list)
    with gzip.open(bed_gz, 'rt') as f:
        for line in f:
            if not line.startswith('chr'):
                continue
            parts = line.rstrip().split('\t')
            if len(parts) < 5:
                continue
            if keyword.lower() not in parts[3].lower():
                continue
            try:
                chrom = parts[0]
                start = int(parts[1])
                end   = int(parts[2])
                score = float(parts[4])
                peaks[chrom].append((start, end, score))
            except:
                pass
    return peaks

def score_locus(peaks, chrom, wstart, wend):
    hits = [s for (s, e, sc) in peaks.get(chrom, [])
            if s < wend and e > wstart]
    scores = [sc for (s, e, sc) in peaks.get(chrom, [])
              if s < wend and e > wstart]
    return len(hits), max(scores) if scores else 0.0

# ── 全細胞種のpeakを読み込み ─────────────────────────────────────────────
print("Loading BED files...")
cell_peaks = {}
for ct_name, (bed_file, keyword, _) in CELL_TYPES.items():
    print(f"  {ct_name} ({keyword})...", end=" ", flush=True)
    cell_peaks[ct_name] = load_bed(bed_file, keyword)
    total = sum(len(v) for v in cell_peaks[ct_name].values())
    print(f"{total} peaks")

# ── スコア計算 ────────────────────────────────────────────────────────────
all_rows = []
for disease, genes in DISEASE_GENES.items():
    for gene, (chrom, wstart, wend) in genes.items():
        row = {"Disease": disease, "Gene": gene}
        for ct_name in CELL_TYPES:
            n, s = score_locus(cell_peaks[ct_name], chrom, wstart, wend)
            row[ct_name] = round(s, 1)
        all_rows.append(row)

df = pd.DataFrame(all_rows)

# ── 保存 ──────────────────────────────────────────────────────────────────
out_csv = f"{OUT_DIR}/celltype_disease_vdr_gr.csv"
df.to_csv(out_csv, index=False)
print(f"\nSaved: {out_csv}")
print(f"Shape: {df.shape}")

# ── サマリー表示 ──────────────────────────────────────────────────────────
vdr_cols = [c for c in df.columns if c.startswith("VDR")]
gr_cols  = [c for c in df.columns if c.startswith("GR")]

print("\n" + "="*90)
for disease in DISEASE_GENES:
    sub = df[df["Disease"] == disease].copy()
    print(f"\n【{disease}】")
    print(f"  {'Gene':<12}", end="")
    for c in vdr_cols + gr_cols:
        print(f"  {c:<14}", end="")
    print()
    print("  " + "-"*80)
    for _, row in sub.iterrows():
        print(f"  {row['Gene']:<12}", end="")
        for c in vdr_cols + gr_cols:
            v = row[c]
            marker = "★" if v > 50 else ("◆" if v > 10 else " ")
            print(f"  {v:>6.1f}{marker:<7}", end="")
        print()
