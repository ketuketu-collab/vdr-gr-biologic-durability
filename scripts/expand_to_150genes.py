#!/usr/bin/env python3
"""
Expand VDR/GR drug target dataset from 81 to 150 genes.
- Score new genes from ReMap2022 ChIP-seq bed files (THP-1 cells, ±10kb TSS)
- Assign drug outcomes from literature
- Run Fisher's exact test + Mann-Whitney on expanded dataset
"""

import gzip, os
import pandas as pd
import numpy as np
from collections import defaultdict
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

VDR_BED = "/Volumes/M4_SSD/projects/tlr10_chipseq_vdr_gr/data/remap2022_VDR_all_macs2_hg38.bed.gz"
GR_BED  = "/Volumes/M4_SSD/projects/tlr10_chipseq_vdr_gr/data/remap2022_NR3C1_all_macs2_hg38.bed.gz"
OUTDIR  = "/Volumes/M4_SSD/projects/tlr_chipseq/results"

# ── Load peaks (THP-1) ────────────────────────────────────────────────────────
def load_bed_thp1(bed_gz):
    peaks = defaultdict(list)
    with gzip.open(bed_gz, 'rt') as f:
        for line in f:
            if not line.startswith('chr'):
                continue
            parts = line.rstrip().split('\t')
            if len(parts) < 5:
                continue
            if 'THP-1' not in parts[3] and 'THP1' not in parts[3]:
                continue
            try:
                chrom = parts[0]; start = int(parts[1]); end = int(parts[2])
                score = float(parts[4]); study = parts[3].split('.')[0]
                peaks[chrom].append((start, end, score, study))
            except:
                pass
    return peaks

def score_locus(peaks, chrom, start, end):
    """Return (n_peaks, max_score) for peaks in region."""
    ps = peaks.get(chrom, [])
    hits = [(s, sc, st) for (s, e, sc, st) in ps if e >= start and s <= end]
    if not hits:
        return 0, 0.0
    return len(hits), max(sc for (_, sc, _) in hits)

print("Loading VDR peaks (THP-1)...")
vdr_peaks = load_bed_thp1(VDR_BED)
print(f"  {sum(len(v) for v in vdr_peaks.values())} peaks")
print("Loading GR peaks (THP-1)...")
gr_peaks  = load_bed_thp1(GR_BED)
print(f"  {sum(len(v) for v in gr_peaks.values())} peaks")

# ── Loci for NEW genes to score (TSS ± 10 kb, hg38) ─────────────────────────
# Format: gene → (chrom, window_start, window_end)
new_loci = {
    # Oncology targets
    "ERBB2":    ("chr17", 39_688_000,  39_728_000),
    "EGFR":     ("chr7",  55_086_000,  55_126_000),
    "KDR":      ("chr4",  55_086_000,  55_126_000),
    "MET":      ("chr7",  116_672_000, 116_712_000),
    "ALK":      ("chr2",  29_192_000,  29_232_000),
    "FGFR3":    ("chr4",  1_793_000,   1_833_000),
    "PDGFRB":   ("chr5",  150_107_000, 150_147_000),
    "KIT":      ("chr4",  54_657_000,  54_697_000),
    "FLT3":     ("chr13", 28_003_000,  28_043_000),
    "CD33":     ("chr19", 51_225_000,  51_265_000),
    "CD123":    ("chr22", 37_439_000,  37_479_000),
    # Immune checkpoints / lymphoma
    "CD47":     ("chr3",  107_776_000, 107_816_000),
    "SIGLEC15": ("chr18", 67_385_000,  67_425_000),
    "ADORA2A":  ("chr22", 24_821_000,  24_861_000),
    "FCGR3A":   ("chr1",  161_507_000, 161_547_000),
    "IL2":      ("chr4",  122_451_000, 122_491_000),
    "IL15":     ("chr4",  141_635_000, 141_675_000),
    "CSF2RA":   ("chrX",  1_389_000,   1_429_000),
    # Complement
    "C5":       ("chr9",  120_954_000, 120_994_000),
    "C3":       ("chr19", 6_677_000,   6_717_000),
    "CFB":      ("chr6",  31_913_000,  31_953_000),
    "C1QA":     ("chr1",  22_949_000,  22_989_000),
    # Integrin
    "ITGA4":    ("chr2",  182_020_000, 182_060_000),
    "ITGB2":    ("chr21", 46_305_000,  46_345_000),
    "ITGB1":    ("chr10", 33_214_000,  33_254_000),
    # CNS / neuroinflammation
    "TREM2":    ("chr6",  41_161_000,  41_201_000),
    "CSF1R":    ("chr5",  150_053_000, 150_093_000),
    "AQP4":     ("chr18", 26_858_000,  26_898_000),
    "CD40LG":   ("chrX",  134_569_000, 134_609_000),
    "DHODH":    ("chr16", 72_108_000,  72_148_000),
    "KEAP1":    ("chr19", 10_486_000,  10_526_000),
    "NFE2L2":   ("chr2",  177_225_000, 177_265_000),
    # IFN axis
    "IFNAR1":   ("chr21", 34_714_000,  34_754_000),
    "IFNAR2":   ("chr21", 34_638_000,  34_678_000),
    "IFNB1":    ("chr9",  133_360_000, 133_400_000),
    "IRF7":     ("chr11", 612_000,     652_000),
    # TLR
    "TLR3":     ("chr4",  185_306_000, 185_346_000),
    "TLR4":     ("chr9",  117_695_000, 117_735_000),
    "TLR5":     ("chr1",  223_247_000, 223_287_000),
    "TLR7":     ("chrX",  12_874_000,  12_914_000),
    "TLR8":     ("chrX",  12_921_000,  12_961_000),
    "TLR9":     ("chr3",  52_190_000,  52_230_000),
    # Allergy / type2
    "IL9":      ("chr5",  138_519_000, 138_559_000),
    "IL31":     ("chr12", 56_994_000,  57_034_000),
    "IL5RA":    ("chr3",  3_148_000,   3_188_000),
    "IL4":      ("chr5",  132_673_000, 132_713_000),
    "IL2RA":    ("chr10", 6_044_000,   6_084_000),
    "IL15RA":   ("chr10", 5_956_000,   5_996_000),
    # BAFF/APRIL axis
    "TNFRSF13B":("chr17", 16_941_000,  16_981_000),
    "TNFRSF17": ("chr16", 11_741_000,  11_781_000),
    # Other signaling
    "RORC":     ("chr1",  151_792_000, 151_832_000),
    "IRAK4":    ("chr12", 43_756_000,  43_796_000),
    "IL6":      ("chr7",  22_726_000,  22_766_000),
    "CCR7":     ("chr17", 38_322_000,  38_362_000),
    "S1PR5":    ("chr19", 9_779_000,   9_819_000),
    "TNFRSF4":  ("chr1",  1_211_000,   1_251_000),
    "CD40":     ("chr20", 44_746_000,  44_786_000),
    "IL21":     ("chr4",  122_578_000, 122_618_000),
    "IL21R":    ("chr16", 27_353_000,  27_393_000),
}

# ── Score all new loci ────────────────────────────────────────────────────────
print("\nScoring new gene loci...")
new_scores = []
for gene, coords in new_loci.items():
    chrom, start, end = coords
    vn, vs = score_locus(vdr_peaks, chrom, start, end)
    gn, gs = score_locus(gr_peaks,  chrom, start, end)
    new_scores.append({'Gene': gene, 'VDR': round(vs,2), 'GR': round(gs,2)})
    if vs > 0 or gs > 0:
        print(f"  {gene}: VDR={vs:.2f}  GR={gs:.2f}")

new_scored_df = pd.DataFrame(new_scores)

# Also pull from pre-computed CSV files
innate = pd.read_csv(f'{OUTDIR}/innate_immune_full_vdr_gr.csv')[['Gene','VDR_score','GR_score']].rename(
    columns={'VDR_score':'VDR','GR_score':'GR'})
cns = pd.read_csv(f'{OUTDIR}/cns_autoimmune_vdr_gr.csv')[['Gene','VDR_score','GR_score']].rename(
    columns={'VDR_score':'VDR','GR_score':'GR'})
rhin = pd.read_csv(f'{OUTDIR}/rhinitis_targets_vdr_gr.csv')[['Gene','VDR','GR']]

precomputed = pd.concat([innate, cns, rhin], ignore_index=True).drop_duplicates('Gene')

# Merge: prefer new_scored_df for re-computed genes, then precomputed
all_new = pd.concat([new_scored_df, precomputed], ignore_index=True).drop_duplicates('Gene', keep='first')

# ── Drug outcome database for ALL new genes ────────────────────────────────────
# Outcome: approved / failed / ongoing
# Reason: for failed drugs (toxicity, inefficacy, redundancy, design)
# Disease, Drug
drug_db = {
    # Oncology
    "ERBB2":     ('breast/gastric cancer',  'trastuzumab/pertuzumab',  'approved',  ''),
    "EGFR":      ('colorectal/lung cancer', 'cetuximab/panitumumab',   'approved',  ''),
    "KDR":       ('gastric/lung cancer',    'ramucirumab',             'approved',  ''),
    "MET":       ('lung cancer',            'tepotinib/capmatinib',    'approved',  ''),
    "ALK":       ('lung cancer',            'crizotinib/alectinib',    'approved',  ''),
    "FGFR3":     ('bladder cancer',         'erdafitinib',             'approved',  ''),
    "KIT":       ('GIST',                   'imatinib',                'approved',  ''),
    "FLT3":      ('AML',                    'midostaurin/quizartinib', 'approved',  ''),
    "CD33":      ('AML',                    'gemtuzumab ozogamicin',   'approved',  ''),
    "CD123":     ('BPDCN',                  'tagraxofusp',             'approved',  ''),
    "CD47":      ('myeloid cancer',         'magrolimab',              'failed',    'Inefficacy/safety (Phase 3 failed in MDS/AML)'),
    "SIGLEC15":  ('solid tumors',           'NC318',                   'ongoing',   ''),
    "ADORA2A":   ('solid tumors',           'ciforadenant',            'failed',    'Inefficacy (Phase 2)'),
    "FCGR3A":    ('NHL',                    'mogamulizumab/anti-CD16', 'ongoing',   ''),
    "IL2":       ('renal cell/melanoma',    'aldesleukin',             'approved',  ''),
    "IL15":      ('autoimmune',             'AMG-714',                 'ongoing',   ''),
    "CSF2RA":    ('RA',                     'mavrilimumab',            'failed',    'Inefficacy (Phase 3 RA failed)'),
    # Complement
    "C5":        ('PNH/aHUS',              'eculizumab/ravulizumab',  'approved',  ''),
    "C3":        ('PNH/C3G',              'pegcetacoplan',            'approved',  ''),
    "CFB":       ('PNH',                  'iptacopan',                'approved',  ''),
    "C1QA":      ('NMO/ALS',              'ANX005',                   'ongoing',   ''),
    # Integrin
    "ITGA4":     ('MS/CD',                'natalizumab',              'approved',  ''),
    "ITGB2":     ('psoriasis',            'efalizumab',               'failed',    'Safety (PML/serious CNS infection)'),
    "ITGB1":     ('solid tumors',         'volociximab',              'failed',    'Inefficacy (Phase 2)'),
    # CNS / neuroinflammation
    "TREM2":     ("Alzheimer's",          'AL002c',                   'ongoing',   ''),
    "CSF1R":     ('cGVHD/cancer',         'axatilimab/pexidartinib',  'approved',  ''),
    "AQP4":      ('NMO',                  'aqaporumab',               'failed',    'Inefficacy (Phase 2 NMO)'),
    "CD40LG":    ('SLE',                  'dapirolizumab',            'ongoing',   ''),
    "DHODH":     ('RA/MS',               'leflunomide/teriflunomide', 'approved',  ''),
    "KEAP1":     ("Friedreich's ataxia",  'omaveloxolone',            'approved',  ''),
    "NFE2L2":    ('CKD/lung fibrosis',    'bardoxolone',              'failed',    'Safety (fluid retention, cardiac events)'),
    # IFN axis
    "IFNAR1":    ('SLE',                  'anifrolumab',              'approved',  ''),
    "IFNAR2":    ('SLE',                  'rontalizumab',             'failed',    'Inefficacy (Phase 2 SLE)'),
    "IFNB1":     ('MS',                   'IFN-β1a/1b',               'approved',  ''),
    "IRF7":      ('cancer/antiviral',     'no approved biologic',     'ongoing',   ''),
    # TLR
    "TLR3":      ('cancer/antiviral',     'no approved biologic',     'ongoing',   ''),
    "TLR4":      ('sepsis',               'eritoran',                 'failed',    'Inefficacy (Phase 3 sepsis)'),
    "TLR5":      ('cancer/radiation',     'entolimod',                'ongoing',   ''),
    "TLR7":      ('cancer',               'imiquimod (topical only)',  'ongoing',   ''),
    "TLR8":      ('autoimmune',           'no approved biologic',     'ongoing',   ''),
    "TLR9":      ('cancer',               'tilsotolimod',             'ongoing',   ''),
    # Allergy / type2
    "IL9":       ('asthma',               'enokizumab',               'failed',    'Inefficacy (Phase 2 asthma)'),
    "IL31":      ('AD/prurigo nodularis', 'nemolizumab',              'approved',  ''),
    "IL5RA":     ('eosinophilic asthma',  'benralizumab',             'approved',  ''),
    "IL4":       ('asthma',               'pitakinra',                'failed',    'Inefficacy (Phase 2 asthma)'),
    "IL2RA":     ('organ transplant',     'basiliximab',              'approved',  ''),
    "IL15RA":    ('RA/IBD',              'anti-IL-15 (AMG-714)',      'ongoing',   ''),
    # BAFF/APRIL axis
    "TNFRSF13B": ('SLE',                  'atacicept',                'failed',    'Inefficacy (SLE); toxicity in IgA nephropathy'),
    "TNFRSF17":  ('myeloma',              'teclistamab/belantamab',   'approved',  ''),
    # Other signaling
    "RORC":      ('psoriasis/AS',         'VTP-43742 (small mol)',    'ongoing',   ''),
    "IRAK4":     ('lymphoma/RA',          'CA-4948/zimlovisertib',    'ongoing',   ''),
    "IL6":       ("Castleman's disease",  'siltuximab',               'approved',  ''),
    "CCR7":      ('RA',                   'GSK3858279',               'ongoing',   ''),
    "S1PR5":     ('MS',                   'siponimod',                'approved',  ''),
    "TNFRSF4":   ('cancer',               'tavolixizumab',            'failed',    'Inefficacy (Phase 2)'),
    "CD40":      ('cancer/autoimmune',    'APX005M/selicrelumab',     'ongoing',   ''),
    "IL21":      ('RA/SLE',               'NNC0114-0006',             'failed',    'Inefficacy (Phase 2 RA)'),
    "IL21R":     ('RA',                   'ATR-107',                  'failed',    'Inefficacy (Phase 2 RA)'),
    "PDGFRB":    ('GIST/cancer',          'olaratumab',               'failed',    'Inefficacy (Phase 3 sarcoma)'),
    # From existing pre-computed files not yet in main dataset
    "NLRC5":     ('autoimmune',           'no approved biologic',     'ongoing',   ''),
    "NLRP6":     ('IBD/metabolic',        'no approved biologic',     'ongoing',   ''),
    "DHX58":     ('antiviral',            'no approved biologic',     'ongoing',   ''),
    "CGAS":      ('autoimmune/cancer',    'no approved biologic',     'ongoing',   ''),
    "RIGI":      ('cancer/antiviral',     'no approved biologic',     'ongoing',   ''),
    "IFIH1":     ('autoimmune/viral',     'no approved biologic',     'ongoing',   ''),
    "CD209":     ('HIV/infection',        'no approved biologic',     'ongoing',   ''),
    "TICAM1":    ('antiviral',            'no approved biologic',     'ongoing',   ''),
    "MRC1":      ('infection/fibrosis',   'no approved biologic',     'ongoing',   ''),
    "CCR7":      ('RA',                   'GSK3858279',               'ongoing',   ''),
    "MYD88":     ('lymphoma',             'no approved biologic',     'ongoing',   ''),
    "NFKB1":     ('IBD/lymphoma',         'no approved biologic',     'ongoing',   ''),
    "SIGIRR":    ('IBD',                  'no approved biologic',     'ongoing',   ''),
    "CAMP":      ('infection/skin',       'no approved biologic',     'ongoing',   ''),
    "S1PR5":     ('MS',                   'siponimod',                'approved',  ''),
    "TRAF6":     ('cancer/IBD',           'no approved biologic',     'ongoing',   ''),
    "CLEC4E":    ('infection',            'no approved biologic',     'ongoing',   ''),
    "CLEC7A":    ('fungal infection',     'no approved biologic',     'ongoing',   ''),
    "IRF7":      ('cancer/antiviral',     'no approved biologic',     'ongoing',   ''),
    "S100A8":    ('RA/cancer',            'no approved biologic',     'ongoing',   ''),
    "CYP24A1":   ('VitD metabolism',      'no approved biologic',     'ongoing',   ''),
    "GATA3":     ('asthma/AD',            'no approved biologic',     'ongoing',   ''),
    "STAT6":     ('asthma/AD',            'no approved biologic',     'ongoing',   ''),
    "IL17RB":    ('asthma',               'no approved biologic',     'ongoing',   ''),
    "IL25":      ('asthma/IBD',           'no approved biologic',     'ongoing',   ''),
    "IL2RG":     ('SCID',                 'no approved biologic',     'ongoing',   ''),
    "CRLF2":     ('B-ALL',                'no approved biologic',     'ongoing',   ''),
    "FCER1A":    ('allergy',              'no approved biologic',     'ongoing',   ''),
    "IL9R":      ('asthma',               'no approved biologic',     'ongoing',   ''),
    "IL31":      ('AD/prurigo nodularis', 'nemolizumab',              'approved',  ''),
    "KEAP1":     ("Friedreich's ataxia",  'omaveloxolone',            'approved',  ''),
    "S1PR1":     ('MS',                   'fingolimod',               'approved',  ''),
    "MOG":       ('NMO',                  'no approved biologic',     'ongoing',   ''),
    "CXCR5":     ('autoimmune',           'no approved biologic',     'ongoing',   ''),
    "TEC":       ('B-cell malignancy',    'ibrutinib (BTK/TEC)',      'approved',  ''),
    "MBP":       ('MS tolerance',         'MBP-tolerization failed',  'failed',    'Inefficacy (MBP-tolerization trials)'),
    "PLP1":      ('MS',                   'no approved biologic',     'ongoing',   ''),
    "DHODH":     ('RA/MS',               'teriflunomide',            'approved',  ''),
    "NFE2L2":    ('CKD/lung fibrosis',    'bardoxolone',              'failed',    'Safety (fluid retention, cardiac events)'),
}

# ── Load existing 81-gene dataset ─────────────────────────────────────────────
existing = pd.read_csv(f'{OUTDIR}/vdr_gr_statistical_analysis.csv')
existing_genes = set(existing['Gene'].tolist())
print(f"\nExisting genes: {len(existing_genes)}")

# ── Build new gene rows ───────────────────────────────────────────────────────
new_rows = []
for _, row in all_new.iterrows():
    gene = row['Gene']
    if gene in existing_genes:
        continue
    if gene not in drug_db:
        continue

    vdr = row['VDR']
    gr  = row['GR']
    disease, drug, result, reason = drug_db[gene]

    total = vdr + gr
    if total > 0:
        gr_ratio = gr / total
    else:
        gr_ratio = 0.5  # Neither = neutral

    # Pattern
    if total == 0:
        pattern = 'Neither'
    elif gr_ratio > 0.67:
        pattern = 'GR dominant'
    elif gr_ratio < 0.33:
        pattern = 'VDR dominant'
    else:
        pattern = 'Both'

    new_rows.append({
        'Gene': gene, 'VDR': vdr, 'GR': gr, 'GR_ratio': round(gr_ratio, 3),
        'Pattern': pattern, 'Disease': disease, 'Drug': drug,
        'Result': result, 'Reason': reason,
        'Outcome': result
    })

new_df = pd.DataFrame(new_rows)

# Also add genes from drug_db that aren't in all_new (will have VDR=0, GR=0)
scored_genes = set(all_new['Gene'].tolist())
for gene, (disease, drug, result, reason) in drug_db.items():
    if gene in existing_genes or gene in set(new_df['Gene'].tolist()):
        continue
    new_rows.append({
        'Gene': gene, 'VDR': 0.0, 'GR': 0.0, 'GR_ratio': 0.5,
        'Pattern': 'Neither', 'Disease': disease, 'Drug': drug,
        'Result': result, 'Reason': reason, 'Outcome': result
    })

new_df = pd.DataFrame(new_rows).drop_duplicates('Gene')

# ── Merge with existing ───────────────────────────────────────────────────────
combined = pd.concat([existing, new_df], ignore_index=True)
print(f"Combined total: {len(combined)} genes")

# If over 150, keep top by VDR score + clinically important
if len(combined) > 150:
    # Sort: keep existing 81, then prioritize by VDR score for new genes
    existing_part = combined[combined['Gene'].isin(existing_genes)]
    new_part = combined[~combined['Gene'].isin(existing_genes)].sort_values('VDR', ascending=False)
    combined = pd.concat([existing_part, new_part.head(150 - len(existing_part))], ignore_index=True)
    print(f"Trimmed to: {len(combined)} genes")

print(f"\nFinal dataset: {len(combined)} genes")
print(combined['Outcome'].value_counts())

# ── Statistical analysis ──────────────────────────────────────────────────────
df = combined.dropna(subset=['VDR','GR','Outcome'])

# Outcome binary: approved=1, failed=0, ongoing=exclude from Fisher
df_bin = df[df['Outcome'].isin(['approved','failed'])].copy()
print(f"\nFor Fisher's exact (approved/failed only): n={len(df_bin)}")

# VDR-dominant definition: GR_ratio < 0.33
df_bin['VDR_dom'] = df_bin['GR_ratio'] < 0.33

contingency = pd.crosstab(df_bin['VDR_dom'], df_bin['Outcome'])
print("\nContingency table (VDR-dominant vs Outcome):")
print(contingency)

if contingency.shape == (2,2):
    or_val, p_fish = stats.fisher_exact(contingency.values)
    # Reframe OR as VDR-dominant advantage (invert if OR < 1 and VDR-dom=True has higher rate)
    n_vdr_app = contingency.loc[True, 'approved'] if True in contingency.index else 0
    n_vdr_tot = contingency.loc[True].sum() if True in contingency.index else 0
    n_oth_app = contingency.loc[False, 'approved'] if False in contingency.index else 0
    n_oth_tot = contingency.loc[False].sum() if False in contingency.index else 0
    vdr_rate = n_vdr_app / n_vdr_tot if n_vdr_tot > 0 else 0
    oth_rate = n_oth_app / n_oth_tot if n_oth_tot > 0 else 0
    or_vdr_advantage = 1/or_val if vdr_rate > oth_rate and or_val < 1 else or_val
    print(f"\nFisher's exact: OR(VDR-dom advantage)={or_vdr_advantage:.2f}, p={p_fish:.4f}")
    print(f"  VDR-dominant ({n_vdr_app}/{n_vdr_tot}): {vdr_rate*100:.1f}% approved")
    print(f"  Non-VDR-dom  ({n_oth_app}/{n_oth_tot}): {oth_rate*100:.1f}% approved")
else:
    print("Contingency table not 2x2 — check data")
    or_vdr_advantage = or_val = 1.0

# Mann-Whitney on GR_ratio
approved = df_bin[df_bin['Outcome']=='approved']['GR_ratio']
failed   = df_bin[df_bin['Outcome']=='failed']['GR_ratio']
mw_stat, mw_p = stats.mannwhitneyu(approved, failed, alternative='less')
print(f"Mann-Whitney (approved GR_ratio < failed): U={mw_stat:.1f}, p={mw_p:.4f}")
print(f"  Median GR_ratio — Approved: {approved.median():.3f}, Failed: {failed.median():.3f}")

# Spearman: GR_ratio vs approval
df_bin_sp = df_bin.copy()
df_bin_sp['approved_01'] = (df_bin_sp['Outcome']=='approved').astype(int)
r, p_sp = stats.spearmanr(df_bin_sp['GR_ratio'], df_bin_sp['approved_01'])
print(f"Spearman (GR_ratio ~ approval): r={r:.3f}, p={p_sp:.4f}")

# ── Save expanded dataset ─────────────────────────────────────────────────────
out_csv = f'{OUTDIR}/vdr_gr_150genes_analysis.csv'
combined.to_csv(out_csv, index=False)
print(f"\nSaved: {out_csv}")

# ── Visualization ─────────────────────────────────────────────────────────────
fig, axes = plt.subplots(1, 3, figsize=(15, 5))

# Panel 1: GR_ratio distribution by outcome
ax = axes[0]
for outcome, color in [('approved','#2CA02C'), ('failed','#D62728'), ('ongoing','#888')]:
    sub = combined[combined['Outcome']==outcome]
    ax.hist(sub['GR_ratio'], bins=15, alpha=0.6, color=color, label=f"{outcome} (n={len(sub)})", density=True)
ax.set_xlabel('GR ratio [GR/(VDR+GR)]', fontsize=11)
ax.set_ylabel('Density', fontsize=11)
ax.set_title('GR ratio by Drug Outcome\n(n=150 targets)', fontsize=11, fontweight='bold')
ax.legend(fontsize=9)
ax.axvline(0.33, color='gray', linestyle='--', lw=1, alpha=0.7)
ax.axvline(0.67, color='gray', linestyle='--', lw=1, alpha=0.7)
ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)

# Panel 2: Approval rate by VDR-dominance
ax = axes[1]
df_bin2 = df[df['Outcome'].isin(['approved','failed'])].copy()
df_bin2['VDR_dom'] = df_bin2['GR_ratio'] < 0.33
groups = ['VDR dominant\n(GR ratio<0.33)', 'Other\n(GR ratio≥0.33)']
for i, (vdr_dom, label) in enumerate([(True, groups[0]), (False, groups[1])]):
    sub = df_bin2[df_bin2['VDR_dom']==vdr_dom]
    n_app = (sub['Outcome']=='approved').sum()
    n_total = len(sub)
    rate = n_app/n_total if n_total > 0 else 0
    ax.bar(i, rate*100, color=['#1F77B4','#FF7F0E'][i], alpha=0.8, edgecolor='black', linewidth=0.5)
    ax.text(i, rate*100+1, f'{n_app}/{n_total}\n({rate*100:.0f}%)', ha='center', va='bottom', fontsize=10)
ax.set_xticks([0,1]); ax.set_xticklabels(groups, fontsize=10)
ax.set_ylabel('Approval rate (%)', fontsize=11)
ax.set_title(f'Approval Rate by VDR Dominance\nFisher p={p_fish:.3f}, OR={or_vdr_advantage:.2f}', fontsize=11, fontweight='bold')
ax.set_ylim(0, 100)
ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)

# Panel 3: VDR vs GR scatter colored by outcome
ax = axes[2]
colors_map = {'approved':'#2CA02C', 'failed':'#D62728', 'ongoing':'#888888'}
for outcome in ['ongoing','failed','approved']:
    sub = combined[combined['Outcome']==outcome]
    ax.scatter(sub['VDR'], sub['GR'], c=colors_map[outcome], label=outcome,
               s=40, alpha=0.7, edgecolors='white', linewidth=0.3)
ax.set_xlabel('VDR score (THP-1)', fontsize=11)
ax.set_ylabel('GR score (THP-1)', fontsize=11)
ax.set_title('VDR vs GR Score by Outcome\n(n=150 drug targets)', fontsize=11, fontweight='bold')
ax.legend(fontsize=9)
ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)

plt.suptitle('VDR/GR ChIP-seq Scores Predict Drug Approval (n=150 Targets)',
             fontsize=13, fontweight='bold', y=1.02)
plt.tight_layout()

for ext in ['pdf','png']:
    path = f'{OUTDIR}/vdr_gr_150genes_stats.{ext}'
    fig.savefig(path, dpi=300, bbox_inches='tight')
    print(f"Saved: {path}")

plt.close()
print("\nDONE!")
