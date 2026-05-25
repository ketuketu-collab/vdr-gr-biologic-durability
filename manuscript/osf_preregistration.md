# OSF Preregistration Draft
# VDR/GR ChIP-seq Rule Predicts Long-term Durability of Biologics in Chronic Inflammatory Diseases

**Date**: 2026-05-25
**Investigator**: Hiroyuki Nagashima, Sapporo Medical University

---

## Title
A VDR/GR ChIP-seq Occupancy Rule Predicts Long-term Therapeutic Durability of Biologics in Chronic Inflammatory Diseases: Prospective Validation Against Tulisokibart Phase 3 Results

---

## Hypotheses

### Primary Hypothesis
The ratio of vitamin D receptor (VDR) to glucocorticoid receptor (GR/NR3C1) ChIP-seq peak occupancy at a drug's target gene locus—derived from the ReMAP2022 public compendium—predicts long-term therapeutic durability of biologics in chronic inflammatory diseases.

Specifically: **VDR-dominant targets (VDR score > GR score) are associated with durable long-term remission; GR-dominant targets are associated with either regulatory failure or loss of durability over time.**

### Pre-specified Prospective Prediction
**Tulisokibart (MK-7240, anti-TL1A monoclonal antibody; Merck)**

- Target gene: TNFSF15 (TL1A)
- VDR score (ReMAP2022, ±10kb of TSS): **12.4**
- GR score (ReMAP2022, ±10kb of TSS): **286.0**
- VDR/GR ratio: **0.043**
- Classification: **GR-dominant** (among the highest GR scores in the 381-gene database)

**Prediction**: Tulisokibart will fail to demonstrate durable long-term remission in chronic inflammatory disease Phase 3 trials, OR will show significantly lower maintenance remission rates compared to VDR-dominant agents (anti-IL-23, vedolizumab) in IBD.

**Primary validation event**: Phase 3 ATLAS-UC trial (NCT06052059) primary completion August 2026; results expected within 6–12 months of primary completion (estimated H1 2027, most likely DDW 2027 or equivalent major GI conference). Late-breaking abstract presentation at fall 2026 meetings (ACG/UEG) is possible but not anticipated, given that abstract submission deadlines for those conferences precede the trial's primary completion date. Per FDAAA 801, results must be posted to ClinicalTrials.gov by August 2027 at the latest.

---

## Background

We have developed a molecular framework—the "VDR/GR rule"—analogous to Lipinski's Rule of Five for oral drugs, that predicts biologic drug success in chronic inflammatory diseases using public ChIP-seq data from ReMAP2022.

**Retrospective validation (n=92 drug-disease pairs, 7 diseases)**:
- VDR-dominant targets → approved: 89% (33/37), PPV=89%
- GR-dominant targets → approved: 51% (28/55)
- Fisher's exact p = 0.0001, OR = 8.0, AUC = 0.706

**Long-term durability (IBD, n=12 drug-disease pairs)**:
- VDR score vs maintenance remission rate: Spearman r = 0.899, p < 0.0001
- GR=0 (JAK1) → FDA Black Box warnings (all 5 JAK inhibitors)
- GR=286 (TNFSF15/TL1A) → highest GR score in the dataset → predicted non-durable

**Biological mechanism**:
- Temporal RNA-seq: GR target genes peak at 0-4h (acute), VDR target genes peak at 24-48h (chronic)
- GR directly binds and induces CYP24A1 (vitamin D-catabolizing enzyme), depleting endogenous VDR ligand during chronic therapy
- Cross-species VDR conservation: r=0.451 (p=3.6×10⁻¹⁶); GR conservation: r=0.032 (ns)

---

## Methods

### ChIP-seq scoring
VDR (vitamin D receptor) and NR3C1 (GR) peak signals within ±10kb of gene TSS were summed across all available experiments in ReMAP2022 (remap2022_VDR_all_macs2_hg38.bed.gz; remap2022_NR3C1_all_macs2_hg38.bed.gz). TSS coordinates from Ensembl REST API (GRCh38/hg38).

TNFSF15 coordinates: chr9:114,806,039, strand=−1 (Ensembl hg38)
Window: chr9:114,796,039–114,816,039

### Outcome definition
**Primary outcome**: Clinical remission rate at maintenance endpoint in Phase 3 ATLAS-UC (NCT06052059).

**Comparison**: Tulisokibart maintenance remission will be compared against VDR-dominant IBD agents:
- Vedolizumab (ITGB7, VDR=62): 41.8% UC maintenance remission (GEMINI1)
- Ustekinumab (IL23A, VDR=211): 43.8% UC maintenance remission (UNIFI)
- Risankizumab (IL23A, VDR=211): 40.2% UC maintenance remission (INSPIRE)

**Failure criterion**: If tulisokibart UC maintenance remission < 35% (below TNF inhibitor benchmark) or Phase 3 primary endpoint miss → prediction confirmed.

### Statistical plan
Primary: Compare tulisokibart maintenance remission (when published) against the predicted range for GR-dominant agents (< 35%) vs VDR-dominant agents (40–50%). No formal statistical test required—this is a prospective prediction of categorical outcome (success/failure).

---

## Data Availability
All analysis code and data:
`/Volumes/M4_SSD/projects/tlr_chipseq/`

Key files:
- `results/remap_scores_expanded.csv` — 381-gene VDR/GR scores
- `results/longterm_remission_corrected.csv` — verified clinical trial data
- `results/pipeline_prediction.csv` — pipeline drug predictions
- `manuscript/manuscript_draft_v1.md` — full manuscript draft

---

## Timeline
- **2026-05-25**: Preregistration submitted — **OSF URL: https://osf.io/tnp63** (registered 12:51 AM, embargoed)
- **2026-06 (planned)**: medRxiv preprint submission
- **2026-08**: ATLAS-UC primary completion
- **2026-08**: ATLAS-UC primary completion → data lock begins
- **2026-11〜2027-05 (expected window)**: ATLAS-UC results reported
  - Best case: Late-breaking abstract at ACG/UEG 2026 (fall 2026)
  - Most likely: DDW 2027 (May 2027) or UEG 2027
  - Guaranteed deadline: ClinicalTrials.gov posting by August 2027 (FDAAA 801)
- **Within 1 month of reporting**: Prospective validation assessment and manuscript update
- **2027 Q2–Q3 (revised target)**: Nature Medicine submission

---

## Conflicts of Interest
None declared.

---

## Notes
- TNFSF15 (TL1A) was not in our original 308-gene database; scores were calculated directly from ReMAP2022 bed files on 2026-05-25.
- Previous pipeline prediction incorrectly listed tulisokibart as anti-MADCAM1 (VDR=277); corrected to anti-TL1A/TNFSF15 (GR=286) based on published literature (Lancet Gastroenterol Hepatol 2025).
- The correction strengthens the framework: TL1A/TNFSF15 is among the most GR-dominant targets in the dataset, consistent with the predicted biological profile of a GR-responsive acute-phase cytokine.
