#!/usr/bin/env python3
"""
Full innate immune receptor VDR/GR ChIP-seq scan
TLR, NLR, RLR, CLR, cGAS-STING, downstream signaling
THP-1 cells, ReMap2022 hg38
"""
import gzip, os
import pandas as pd
from collections import defaultdict

VDR_BED = "/Volumes/M4_SSD/projects/tlr10_chipseq_vdr_gr/data/remap2022_VDR_all_macs2_hg38.bed.gz"
GR_BED  = "/Volumes/M4_SSD/projects/tlr10_chipseq_vdr_gr/data/remap2022_NR3C1_all_macs2_hg38.bed.gz"

def load_bed_thp1(bed_gz, cell_keyword="THP-1"):
    peaks = defaultdict(list)
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
                chrom = parts[0]; start = int(parts[1]); end = int(parts[2])
                score = float(parts[4]); study = name.split('.')[0]
                peaks[chrom].append((start, end, score, study))
            except:
                pass
    return peaks

print("Loading VDR peaks..."); vdr_peaks = load_bed_thp1(VDR_BED)
print(f"  THP-1 VDR: {sum(len(v) for v in vdr_peaks.values())} peaks")
print("Loading GR peaks...");  gr_peaks  = load_bed_thp1(GR_BED)
print(f"  THP-1 GR:  {sum(len(v) for v in gr_peaks.values())} peaks")

# ── Gene windows (TSS ± 10 kb, hg38) ──────────────────────────────────
loci = {
    # ── TLR family ──────────────────────────────────────────────────
    "TLR1":   ("chr4",  38760000,  38780000),   # cluster (overlaps TLR10)
    "TLR10":  ("chr4",  38757000,  38777000),   # cluster
    "TLR6":   ("chr4",  38796000,  38816000),   # cluster
    "TLR2":   ("chr4",  153674000, 153694000),
    "TLR3":   ("chr4",  185306000, 185326000),
    "TLR4":   ("chr9",  117695000, 117715000),
    "TLR5":   ("chr1",  223247000, 223267000),
    "TLR7":   ("chrX",  12874000,  12894000),
    "TLR8":   ("chrX",  12921000,  12941000),
    "TLR9":   ("chr3",  52256000,  52276000),
    # ── NLR family ──────────────────────────────────────────────────
    "NOD1":   ("chr7",  30466000,  30486000),
    "NOD2":   ("chr16", 50713000,  50733000),
    "NLRP1":  ("chr17", 5549419,   5569419),
    "NLRP3":  ("chr1",  247368000, 247388000),
    "NLRP6":  ("chr11", 271897,    291897),
    "NLRC4":  ("chr2",  32235111,  32255111),
    "NLRC5":  ("chr16", 57026861,  57046861),
    "CIITA":  ("chr16", 10894621,  10914621),
    # ── RLR (RIG-I-like receptors) ───────────────────────────────────
    "RIGI":   ("chr9",  32480825,  32500825),   # DDX58
    "IFIH1":  ("chr2",  162282879, 162302879),  # MDA5
    "DHX58":  ("chr17", 42097068,  42117068),   # LGP2
    # ── cGAS-STING pathway ───────────────────────────────────────────
    "CGAS":   ("chr6",  73427620,  73447620),
    "STING1": ("chr5",  139469231, 139489231),
    "IRF3":   ("chr19", 49652715,  49672715),
    "IRF7":   ("chr11", 604267,    624267),
    # ── CLR (C-type lectin receptors) ────────────────────────────────
    "CLEC7A": ("chr12", 10113517,  10133517),   # Dectin-1
    "CLEC6A": ("chr12", 8457146,   8477146),    # Dectin-2
    "CLEC4E": ("chr12", 8527105,   8547105),    # Mincle
    "CD209":  ("chr19", 7733778,   7753778),    # DC-SIGN
    "MRC1":   ("chr10", 17850089,  17870089),   # Mannose receptor
    # ── Adaptor molecules ────────────────────────────────────────────
    "MYD88":  ("chr3",  38130788,  38150788),
    "TICAM1": ("chr19", 4813822,   4833822),    # TRIF
    "IRAK4":  ("chr12", 43768625,  43788625),
    "IRAK1":  ("chrX",  154005215, 154025215),
    "TRAF6":  ("chr11", 36487030,  36507030),
    # ── Transcription factors / effectors ────────────────────────────
    "NFKB1":  ("chr4",  102491000, 102511000),
    "NFKBIA": ("chr14", 35358000,  35378000),
    "IFNB1":  ("chr9",  21067523,  21087523),
    "IL6":    ("chr7",  22718943,  22738943),
    "IL1B":   ("chr2",  112823283, 112843283),
    "TNF":    ("chr6",  31566950,  31586950),
    # ── Negative regulators ──────────────────────────────────────────
    "TNFAIP3":("chr6",  137856000, 137876000),  # A20
    "TOLLIP": ("chr11", 1288000,   1308000),
    "SIGIRR": ("chr11", 64039000,  64059000),
    "DUSP1":  ("chr5",  172758000, 172778000),  # MKP-1
    # ── Checkpoint (irAE targets, for comparison) ────────────────────
    "CD274":  ("chr9",  5440000,   5480000),    # PD-L1
    "PDCD1":  ("chr2",  241840000, 241860000),  # PD-1
    "LAG3":   ("chr12", 6762000,   6782000),
    "HAVCR2": ("chr5",  156990000, 157020000),  # TIM-3
    "CTLA4":  ("chr2",  203858000, 203878000),
    "TIGIT":  ("chr3",  114352000, 114372000),
    # ── VDR/GR positive controls ─────────────────────────────────────
    "CAMP":   ("chr3",  48599000,  48619000),
    "CYP24A1":("chr20", 52740000,  52760000),
    "S100A8": ("chr1",  153330000, 153350000),
}

# Category mapping
cat_map = {
    "TLR1":"TLR","TLR2":"TLR","TLR3":"TLR","TLR4":"TLR","TLR5":"TLR",
    "TLR6":"TLR","TLR7":"TLR","TLR8":"TLR","TLR9":"TLR","TLR10":"TLR",
    "NOD1":"NLR","NOD2":"NLR","NLRP1":"NLR","NLRP3":"NLR","NLRP6":"NLR",
    "NLRC4":"NLR","NLRC5":"NLR","CIITA":"NLR",
    "RIGI":"RLR","IFIH1":"RLR","DHX58":"RLR",
    "CGAS":"cGAS-STING","STING1":"cGAS-STING","IRF3":"cGAS-STING","IRF7":"cGAS-STING",
    "CLEC7A":"CLR","CLEC6A":"CLR","CLEC4E":"CLR","CD209":"CLR","MRC1":"CLR",
    "MYD88":"Adaptor","TICAM1":"Adaptor","IRAK4":"Adaptor","IRAK1":"Adaptor","TRAF6":"Adaptor",
    "NFKB1":"Effector","NFKBIA":"Effector","IFNB1":"Effector",
    "IL6":"Effector","IL1B":"Effector","TNF":"Effector",
    "TNFAIP3":"NegReg","TOLLIP":"NegReg","SIGIRR":"NegReg","DUSP1":"NegReg",
    "CD274":"Checkpoint","PDCD1":"Checkpoint","LAG3":"Checkpoint",
    "HAVCR2":"Checkpoint","CTLA4":"Checkpoint","TIGIT":"Checkpoint",
    "CAMP":"VDR_ctrl","CYP24A1":"VDR_ctrl","S100A8":"VDR_ctrl",
}

def get_best_hit(peaks_dict, chrom, start, end):
    if chrom not in peaks_dict:
        return 0.0, "", 0
    best_score = 0.0; best_study = ""; n_peaks = 0
    for (ps, pe, score, study) in peaks_dict[chrom]:
        if pe < start or ps > end:
            continue
        n_peaks += 1
        if score > best_score:
            best_score = score; best_study = study
    return round(best_score, 2), best_study, n_peaks

rows = []
for gene, (chrom, start, end) in loci.items():
    vdr_score, _, vdr_n = get_best_hit(vdr_peaks, chrom, start, end)
    gr_score,  _, gr_n  = get_best_hit(gr_peaks,  chrom, start, end)
    if vdr_score >= 3.0 and gr_score < 1.0:   pattern = "VDR only"
    elif gr_score >= 3.0 and vdr_score < 1.0: pattern = "GR only"
    elif vdr_score >= 3.0 and gr_score >= 3.0: pattern = "Both"
    else:                                       pattern = "Neither"
    rows.append({
        "Gene": gene, "Category": cat_map.get(gene, "Other"),
        "Chr": chrom,
        "VDR_n": vdr_n, "VDR_score": vdr_score,
        "GR_n":  gr_n,  "GR_score":  gr_score,
        "Pattern": pattern,
    })

df = pd.DataFrame(rows)
cat_order = ["TLR","NLR","RLR","cGAS-STING","CLR","Adaptor","Effector","NegReg","Checkpoint","VDR_ctrl"]
df["cat_idx"] = df["Category"].map({c: i for i, c in enumerate(cat_order)})
df = df.sort_values(["cat_idx","Gene"]).drop("cat_idx", axis=1).reset_index(drop=True)

out_dir = "/Volumes/M4_SSD/projects/tlr_chipseq/results"
df.to_csv(f"{out_dir}/innate_immune_full_vdr_gr.csv", index=False)
df.to_csv(f"{out_dir}/innate_immune_full_vdr_gr.tsv", sep="\t", index=False)

print("\n=== Full Innate Immune VDR/GR ChIP-seq Table ===\n")
print(df.to_string(index=False))

print("\n=== Pattern Summary ===")
print(df.groupby(["Category","Pattern"]).size().unstack(fill_value=0).to_string())
print(f"\nSaved: {out_dir}/innate_immune_full_vdr_gr.csv")
