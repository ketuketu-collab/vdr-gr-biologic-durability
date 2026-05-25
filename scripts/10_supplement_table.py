#!/usr/bin/env python3
"""
Supplement Table: VDR and GR ChIP-seq binding scores
THP-1 specific (VDR: +1,25(OH)2D3; GR/NR3C1: +Dexamethasone)
ReMap2022, hg38
"""
import gzip, os
import pandas as pd
from collections import defaultdict

VDR_BED = "/Volumes/M4_SSD/projects/tlr10_chipseq_vdr_gr/data/remap2022_VDR_all_macs2_hg38.bed.gz"
GR_BED  = "/Volumes/M4_SSD/projects/tlr10_chipseq_vdr_gr/data/remap2022_NR3C1_all_macs2_hg38.bed.gz"

def load_bed_thp1(bed_gz, cell_keyword="THP-1"):
    """ファイルを一度読み込み、THP-1ピークをchromごとに格納"""
    peaks = defaultdict(list)  # chrom -> [(start, end, score, study)]
    with gzip.open(bed_gz, 'rt') as f:
        for line in f:
            if not line.startswith('chr'):
                continue
            parts = line.rstrip().split('\t')
            if len(parts) < 5:
                continue
            name = parts[3]
            if cell_keyword.lower() not in name.lower():
                continue
            try:
                chrom = parts[0]
                start = int(parts[1])
                end   = int(parts[2])
                score = float(parts[4])
                study = name.split('.')[0]
                peaks[chrom].append((start, end, score, study))
            except:
                pass
    return peaks

print("VDRデータ読み込み中...")
vdr_peaks = load_bed_thp1(VDR_BED, "THP-1")
print(f"  THP-1 VDRピーク: {sum(len(v) for v in vdr_peaks.values())}")

print("GRデータ読み込み中...")
gr_peaks = load_bed_thp1(GR_BED, "THP-1")
print(f"  THP-1 GRピーク: {sum(len(v) for v in gr_peaks.values())}")

# TSS ± 10kb window, hg38
loci = {
    # TLR family (chr4 cluster)
    "TLR10": ("chr4",  38757000,  38777000),
    "TLR1":  ("chr4",  38760000,  38780000),
    "TLR6":  ("chr4",  38796000,  38816000),
    "TLR2":  ("chr4",  153674000, 153694000),
    "TLR4":  ("chr9",  117695000, 117715000),
    "TLR5":  ("chr1",  223247000, 223267000),
    "TLR3":  ("chr4",  185306000, 185326000),
    "TLR7":  ("chrX",  12874000,  12894000),
    "TLR8":  ("chrX",  12921000,  12941000),
    "TLR9":  ("chr3",  52256000,  52276000),
    # NLR family
    "NOD1":  ("chr7",  30466000,  30486000),
    "NOD2":  ("chr16", 50713000,  50733000),
    "NLRP3": ("chr1",  247368000, 247388000),
    # Checkpoint molecules
    "CD274": ("chr9",  5440000,   5480000),    # PD-L1
    "PDCD1": ("chr2",  241840000, 241860000),  # PD-1
    "CTLA4": ("chr2",  203858000, 203878000),
    "HAVCR2":("chr5",  156990000, 157020000),  # TIM-3
    "LAG3":  ("chr12", 6762000,   6782000),
    "TIGIT": ("chr3",  114352000, 114372000),
    # Downstream signaling
    "TNFAIP3":("chr6", 137856000, 137876000),  # A20 (positive control GR)
    "DUSP1": ("chr5",  172758000, 172778000),  # MKP-1 (positive control GR)
    "NFKB1": ("chr4",  102491000, 102511000),
    "NFKBIA":("chr14", 35358000,  35378000),
    # VDR target genes (positive controls VDR)
    "CAMP":  ("chr3",  48599000,  48619000),   # cathelicidin
    "CYP24A1":("chr20",52740000,  52760000),
    "S100A8":("chr1",  153330000, 153350000),
}

def get_best_hit(peaks_dict, chrom, start, end):
    """事前読み込み済みピークから指定領域のベストヒットを返す"""
    if chrom not in peaks_dict:
        return 0.0, "", 0
    best_score = 0.0
    best_study = ""
    n_peaks = 0
    for (ps, pe, score, study) in peaks_dict[chrom]:
        if pe < start or ps > end:
            continue
        n_peaks += 1
        if score > best_score:
            best_score = score
            best_study = study
    return round(best_score, 2), best_study, n_peaks

rows = []
for gene, (chrom, start, end) in loci.items():
    vdr_score, vdr_study, vdr_n = get_best_hit(vdr_peaks, chrom, start, end)
    gr_score,  gr_study,  gr_n  = get_best_hit(gr_peaks,  chrom, start, end)

    # binding pattern classification
    if vdr_score >= 3.0 and gr_score < 1.0:
        pattern = "VDR only"
    elif gr_score >= 3.0 and vdr_score < 1.0:
        pattern = "GR only"
    elif vdr_score >= 3.0 and gr_score >= 3.0:
        pattern = "Both"
    else:
        pattern = "Neither"

    rows.append({
        "Gene":              gene,
        "Chr":               chrom,
        "Category":          ("TLR" if gene.startswith("TLR") else
                              "NLR" if gene in ("NOD1","NOD2","NLRP3") else
                              "Checkpoint" if gene in ("CD274","PDCD1","CTLA4","HAVCR2","LAG3","TIGIT") else
                              "Downstream/Control"),
        "VDR peaks (THP-1)": vdr_n,
        "VDR max score":     vdr_score,
        "VDR dataset":       vdr_study,
        "GR peaks (THP-1)":  gr_n,
        "GR max score":      gr_score,
        "GR dataset":        gr_study,
        "Binding pattern":   pattern,
    })

df = pd.DataFrame(rows)
df = df.sort_values(["Category", "Gene"]).reset_index(drop=True)

# ── 出力 ─────────────────────────────────────────────────────────────────────
out_dir = "/Volumes/M4_SSD/projects/tlr_chipseq/results"
os.makedirs(out_dir, exist_ok=True)

# CSV (full)
df.to_csv(f"{out_dir}/supplement_table_vdr_gr_chipseq.csv", index=False)

# TSV (for Word/Excel paste)
df.to_csv(f"{out_dir}/supplement_table_vdr_gr_chipseq.tsv", sep="\t", index=False)

# Pretty print
print("\n=== Supplement Table: VDR and GR ChIP-seq (THP-1) ===\n")
print(df.to_string(index=False))

# Pattern summary
print("\n=== Pattern Summary ===")
print(df.groupby(["Category", "Binding pattern"]).size().to_string())

print(f"\n保存: {out_dir}/supplement_table_vdr_gr_chipseq.csv")
