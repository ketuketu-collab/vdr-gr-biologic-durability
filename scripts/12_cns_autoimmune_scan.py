#!/usr/bin/env python3
"""
CNS autoimmune disease (MS/NMOSD/MOGAD) VDR/GR ChIP-seq scan
THP-1 cells, ReMap2022 hg38
Targets: approved DMTs + failed trials + autoantigens
"""
import gzip
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

print("Loading VDR peaks (THP-1)..."); vdr_peaks = load_bed_thp1(VDR_BED)
print(f"  VDR: {sum(len(v) for v in vdr_peaks.values())} peaks")
print("Loading GR peaks (THP-1)..."); gr_peaks = load_bed_thp1(GR_BED)
print(f"  GR : {sum(len(v) for v in gr_peaks.values())} peaks")

# ── CNS autoimmune targets (TSS ± 10 kb, hg38) ───────────────────────
loci = {
    # === MS approved DMTs ===
    "MS4A1":    ("chr11", 60547000,  60577000),  # CD20 — ocrelizumab, ofatumumab
    "ITGA4":    ("chr2",  181447000, 181477000), # α4 — natalizumab (VLA-4)
    "ITGB1":    ("chr10", 32900000,  32930000),  # β1 — pairs w/ITGA4
    "CD52":     ("chr1",  26312000,  26332000),  # alemtuzumab
    "S1PR1":    ("chr1",  101230000, 101260000), # fingolimod, ozanimod
    "S1PR5":    ("chr19", 10330000,  10360000),  # siponimod
    "DHODH":    ("chr16", 72040000,  72070000),  # teriflunomide
    "NFE2L2":   ("chr2",  177220000, 177250000), # Nrf2 (DMF pathway)
    "KEAP1":    ("chr19", 10486000,  10516000),  # DMF target
    # === MS/MS-related Bruton's tyrosine kinase (recent) ===
    "BTK":      ("chrX",  101349000, 101379000), # evobrutinib, tolebrutinib
    "TEC":      ("chr4",  47795000,  47825000),  # related kinase
    # === NMOSD approved ===
    "IL6R":     ("chr1",  154395000, 154425000), # satralizumab
    "C5":       ("chr9",  121399000, 121429000), # eculizumab, ravulizumab
    "CD19":     ("chr16", 28931000,  28961000),  # inebilizumab
    # === MOGAD / NMOSD autoantigens ===
    "AQP4":     ("chr18", 26851000,  26881000),  # NMOSD autoantigen
    "MOG":      ("chr6",  29650000,  29680000),  # MOGAD autoantigen
    "MBP":      ("chr18", 76989000,  77019000),  # myelin basic protein
    "PLP1":     ("chrX",  103770000, 103800000), # myelin proteolipid
    # === IFN-β axis (IFN-β1a/1b drug class) ===
    "IFNAR1":   ("chr21", 33324000,  33380000),  # type I IFN receptor
    "IFNAR2":   ("chr21", 33252000,  33282000),
    "IFNB1":    ("chr9",  21067000,  21097000),  # interferon-β
    # === Th17 axis / failed in MS but relevant ===
    "IL17A":    ("chr6",  52186000,  52216000),  # secukinumab (failed MS)
    "IL17RA":   ("chr22", 17084000,  17114000),
    "IL12B":    ("chr5",  159314000, 159344000), # ustekinumab (failed MS)
    "IL23A":    ("chr12", 56330000,  56360000),  # guselkumab (failed MS)
    "RORC":     ("chr1",  151806000, 151836000), # RORγt master TF
    # === B-cell axis (broader than CD20) ===
    "TNFSF13B": ("chr13", 108250000, 108280000), # BAFF — belimumab (NMOSD trials)
    "TNFRSF13B":("chr17", 16940000,  16970000),  # TACI — atacicept (failed MS)
    "TNFRSF17": ("chr16", 11960000,  11990000),  # BCMA
    "CXCR5":    ("chr11", 118752000, 118782000), # Tfh migration
    # === T-cell costimulation / migration ===
    "CTLA4":    ("chr2",  203858000, 203888000), # abatacept trials
    "CD40LG":   ("chrX",  136648000, 136678000), # failed (thrombosis)
    "ITGAL":    ("chr16", 30472000,  30502000),  # LFA-1 α
    "ITGB2":    ("chr21", 44885000,  44915000),  # LFA-1/Mac-1 β
    "CCR7":     ("chr17", 40540000,  40570000),
    "CXCR3":    ("chrX",  71608000,  71638000),
    # === Complement axis (broader) ===
    "C3":       ("chr19", 6677000,   6707000),
    "C1QA":     ("chr1",  22636000,  22666000),
    "CFB":      ("chr6",  31935000,  31965000),
    # === Microglia / TREM ===
    "TREM2":    ("chr6",  41155000,  41185000),
    "CSF1R":    ("chr5",  150053000, 150083000), # PLX5622-like microglia depletion
    # === EBV-related (MS etiology, Lancet 2022 Kuhle/Bjornevik) ===
    "EBNA1":    None,  # viral; skip
    # === VDR/GR positive controls ===
    "CAMP":     ("chr3",  48599000,  48619000),
    "CYP24A1":  ("chr20", 52740000,  52770000),
}
loci = {k: v for k, v in loci.items() if v is not None}

# Drug annotation
drug = {
    "MS4A1":    "ocrelizumab/ofatumumab (anti-CD20)",
    "ITGA4":    "natalizumab (anti-VLA-4)",
    "ITGB1":    "natalizumab pair",
    "CD52":     "alemtuzumab",
    "S1PR1":    "fingolimod/ozanimod",
    "S1PR5":    "siponimod",
    "DHODH":    "teriflunomide",
    "NFE2L2":   "DMF/diroximel (Nrf2 axis)",
    "KEAP1":    "DMF target",
    "BTK":      "evobrutinib/tolebrutinib (trial)",
    "TEC":      "BTKi off-target",
    "IL6R":     "satralizumab (NMOSD)",
    "C5":       "eculizumab/ravulizumab (NMOSD)",
    "CD19":     "inebilizumab (NMOSD)",
    "AQP4":     "NMOSD autoantigen",
    "MOG":      "MOGAD autoantigen",
    "MBP":      "myelin autoantigen",
    "PLP1":     "myelin autoantigen",
    "IFNAR1":   "IFN-β class target",
    "IFNAR2":   "IFN-β class",
    "IFNB1":    "interferon-β (drug itself)",
    "IL17A":    "secukinumab — FAILED MS",
    "IL17RA":   "brodalumab",
    "IL12B":    "ustekinumab — FAILED MS",
    "IL23A":    "guselkumab — FAILED MS",
    "RORC":     "RORγt inhibitors",
    "TNFSF13B": "belimumab (NMOSD trial)",
    "TNFRSF13B":"atacicept — WORSENED MS",
    "TNFRSF17": "anti-BCMA CAR-T",
    "CXCR5":    "Tfh migration",
    "CTLA4":    "abatacept (small trial)",
    "CD40LG":   "FAILED (thrombosis)",
    "ITGAL":    "efalizumab (withdrawn-PML)",
    "ITGB2":    "LFA-1 β",
    "CCR7":     "lymph node homing",
    "CXCR3":    "Th1 migration",
    "C3":       "pegcetacoplan-class",
    "C1QA":     "ANX005-class",
    "CFB":      "iptacopan",
    "TREM2":    "microglia receptor",
    "CSF1R":    "microglia depletion",
    "CAMP":     "VDR positive ctrl",
    "CYP24A1":  "VDR positive ctrl",
}

cat = {
    "MS4A1":"MS-B-cell", "CD19":"MS-B-cell", "TNFSF13B":"MS-B-cell",
    "TNFRSF13B":"MS-B-cell", "TNFRSF17":"MS-B-cell", "CXCR5":"MS-B-cell",
    "ITGA4":"MS-migration", "ITGB1":"MS-migration", "ITGAL":"MS-migration",
    "ITGB2":"MS-migration", "CCR7":"MS-migration", "CXCR3":"MS-migration",
    "CD52":"MS-broad-immune",
    "S1PR1":"MS-S1P", "S1PR5":"MS-S1P",
    "DHODH":"MS-metabolic", "NFE2L2":"MS-metabolic", "KEAP1":"MS-metabolic",
    "BTK":"MS-kinase", "TEC":"MS-kinase",
    "IL6R":"NMOSD", "C5":"NMOSD", "AQP4":"NMOSD-Ag",
    "MOG":"MOGAD-Ag", "MBP":"Myelin-Ag", "PLP1":"Myelin-Ag",
    "IFNAR1":"IFN-axis","IFNAR2":"IFN-axis","IFNB1":"IFN-axis",
    "IL17A":"Th17-failed","IL17RA":"Th17-failed","IL12B":"Th17-failed",
    "IL23A":"Th17-failed","RORC":"Th17-failed",
    "CTLA4":"T-costim","CD40LG":"T-costim",
    "C3":"Complement","C1QA":"Complement","CFB":"Complement",
    "TREM2":"Microglia","CSF1R":"Microglia",
    "CAMP":"VDR-ctrl","CYP24A1":"VDR-ctrl",
}

def get_best(peaks_dict, chrom, start, end):
    if chrom not in peaks_dict:
        return 0.0, 0
    best = 0.0; n = 0
    for ps, pe, sc, _ in peaks_dict[chrom]:
        if pe < start or ps > end:
            continue
        n += 1
        if sc > best:
            best = sc
    return round(best, 2), n

rows = []
for gene, (chrom, s, e) in loci.items():
    v, vn = get_best(vdr_peaks, chrom, s, e)
    g, gn = get_best(gr_peaks,  chrom, s, e)
    if   v >= 3 and g <  1: pat = "VDR only"
    elif g >= 3 and v <  1: pat = "GR only"
    elif v >= 3 and g >= 3: pat = "Both"
    else:                    pat = "Neither"
    rows.append({"Gene":gene, "Category":cat.get(gene,"Other"),
                 "VDR_score":v, "VDR_n":vn,
                 "GR_score":g,  "GR_n":gn,
                 "Pattern":pat, "Drug":drug.get(gene,"")})

df = pd.DataFrame(rows).sort_values(["Category","VDR_score"], ascending=[True, False]).reset_index(drop=True)
out = "/Volumes/M4_SSD/projects/tlr_chipseq/results/cns_autoimmune_vdr_gr.csv"
df.to_csv(out, index=False)

print("\n=== CNS Autoimmune (MS/NMOSD/MOGAD) VDR/GR ChIP-seq ===\n")
print(df.to_string(index=False))
print("\n=== Pattern summary ===")
print(df["Pattern"].value_counts().to_string())
print("\n=== By category ===")
print(df.groupby(["Category","Pattern"]).size().unstack(fill_value=0).to_string())
print(f"\nSaved: {out}")
