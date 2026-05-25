#!/usr/bin/env python3
"""
VDR and GR ChIP-seq binding at all major PRR loci (TLR, NLR, RLR, cGAS-STING)
Human THP-1 hg38
"""
import subprocess, os
import pandas as pd
import numpy as np

VDR_BED = "/Volumes/M4_SSD/projects/tlr10_chipseq_vdr_gr/data/remap2022_VDR_all_macs2_hg38.bed.gz"
GR_BED  = "/Volumes/M4_SSD/projects/tlr_chipseq/data/remap2022_NR3C1_THP1_hg38.bed.gz"

# ── PRR gene coordinates (hg38, TSS ± 10kb) ──────────────────────────────────
prr_genes = {
    # TLR family
    "TLR1":  ("chr4",  38742000,  38762000),
    "TLR2":  ("chr4",  153686000, 153706000),
    "TLR6":  ("chr4",  38750000,  38770000),
    "TLR10": ("chr4",  38760000,  38780000),
    "TLR4":  ("chr9",  117705000, 117725000),
    "TLR5":  ("chr1",  223257000, 223277000),
    "TLR3":  ("chr4",  185316000, 185336000),
    "TLR7":  ("chrX",  12884000,  12904000),
    "TLR8":  ("chrX",  12931000,  12951000),
    "TLR9":  ("chr3",  52266000,  52286000),
    # NLR family
    "NOD1":  ("chr7",  30476000,  30496000),
    "NOD2":  ("chr16", 50723000,  50743000),
    "NLRP3": ("chr1",  247378000, 247398000),
    "NLRP1": ("chr17", 5481000,   5501000),
    "NLRC4": ("chr2",  32416000,  32436000),
    "CARD8": ("chr19", 8540000,   8560000),
    # RLR family
    "DDX58": ("chr9",  32449000,  32469000),   # RIG-I
    "IFIH1": ("chr2",  163034000, 163054000),  # MDA5
    "DHX58": ("chr17", 80067000,  80087000),   # LGP2
    # cGAS-STING pathway
    "CGAS":  ("chr6",  77161000,  77181000),
    "STING1":("chr5",  138848000, 138868000),  # TMEM173
    "IRF3":  ("chr19", 49659000,  49679000),
    "IRF7":  ("chr11", 614000,    634000),
    # Key downstream
    "TNFAIP3":("chr6", 137866000, 137886000),  # A20
    "NFKB1": ("chr4",  102501000, 102521000),
    "NFKBIA":("chr14", 35368000,  35388000),   # IkBa
}

def get_max_score(bed_gz, chrom, start, end, cell_filter=None):
    cmd = f"tabix {bed_gz} {chrom}:{start}-{end} 2>/dev/null"
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if not result.stdout.strip():
        return 0.0
    best = 0.0
    for line in result.stdout.strip().split('\n'):
        parts = line.split('\t')
        if len(parts) < 5:
            continue
        if cell_filter and cell_filter not in parts[3]:
            continue
        try:
            score = float(parts[4])
            best = max(best, score)
        except:
            pass
    return best

rows = []
for gene, (chrom, start, end) in prr_genes.items():
    vdr = get_max_score(VDR_BED, chrom, start, end, "THP-1")
    gr  = get_max_score(GR_BED,  chrom, start, end)
    rows.append({"gene": gene, "VDR": round(vdr,2), "GR": round(gr,2)})

df = pd.DataFrame(rows)
df["ratio_VDR_GR"] = (df["VDR"] / (df["GR"] + 0.1)).round(2)
df["pattern"] = df.apply(lambda r:
    "VDR-only" if r["VDR"] > 2 and r["GR"] < 1 else
    "GR-only"  if r["GR"]  > 2 and r["VDR"] < 1 else
    "Both"     if r["VDR"] > 2 and r["GR"]  > 2 else
    "Neither", axis=1)

df = df.sort_values("VDR", ascending=False)
print(df.to_string(index=False))
print(f"\nVDR-only: {(df.pattern=='VDR-only').sum()}")
print(f"GR-only:  {(df.pattern=='GR-only').sum()}")
print(f"Both:     {(df.pattern=='Both').sum()}")
print(f"Neither:  {(df.pattern=='Neither').sum()}")
