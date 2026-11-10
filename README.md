# VDR/GR ChIP-seq Biologic Durability Rule

[![OSF Preregistration](https://img.shields.io/badge/OSF-Preregistered-blue)](https://osf.io/tnp63)

## Overview

Code and data for: **"VDR/GR ChIP-seq Occupancy Predicts Long-term Durability of Biologics in Chronic Inflammatory Diseases"**

Target journal: *Nature Medicine*  
OSF Preregistration: https://osf.io/tnp63 (registered 2026-05-25, embargoed)

## Core Finding

The ratio of chromatin occupancy by vitamin D receptor (VDR) vs glucocorticoid receptor (GR/NR3C1) at a drug's target gene—computed from public ChIP-seq data (ReMAP2022)—predicts:
- Biologic approval rate: PPV=89%, Fisher p=0.0001, AUC=0.706
- IBD maintenance remission: Spearman r=0.899, p<0.0001
- Long-term durability tier (durable vs problematic VDR score: 67.3 vs 15.4, p=0.046)

## Prospective Prediction (Pre-registered)

Tulisokibart (MK-7240, anti-TL1A/TNFSF15): VDR=12.4, GR=286.0 (highest GR in 381-gene dataset)  
**Predicted: ATLAS-UC Phase III failure (maintenance remission <35%)**  
Validation: NCT06052059, primary completion August 2026, results expected H1 2027

## Repository Structure

```
scripts/          Analysis and figure generation scripts
  fig3_durability_composite.py    Fig. 3 (durability tiers + IBD scatter + temporal)
  fig5_mechanistic_pipeline.py    Fig. 5 (mechanism + pipeline + trajectory)
  fig_pipeline_validation.py      Suppl. Fig. 3 (pipeline ranking)
results/          Output CSV tables and figures
  remap_scores_expanded.csv       381-gene VDR/GR scores (ReMAP2022)
  pipeline_prediction.csv         Pipeline drug predictions
  longterm_remission_corrected.csv  Clinical remission data (n=19)
  figures/                        Generated figures (PDF + PNG)
manuscript/       Manuscript and preregistration
  manuscript_draft_v2.md          Nature Medicine draft
  osf_preregistration.md          OSF preregistration document
```

## Stratified ATAC motif analysis (GRE vs DR3)

`run_stratified.py` stratifies THP-1 ATAC-seq peaks by LPS response (gained / lost / shared)
and by distance to the nearest TSS (promoter ≤1 kb / proximal ≤10 kb / distal), then tests
GRE (JASPAR MA0113.3, NR3C1) and DR3 (MA0074.1, RXRA::VDR) motif enrichment in gained and lost
peaks against shared peaks of the same stratum (Fisher's exact test, plus a GC-stratified
Cochran–Mantel–Haenszel test; BH-FDR).

```bash
python3 run_stratified.py \
  --unstim THP1_ATAC_untreated_peaks.narrowPeak \
  --stim   THP1_ATAC_LPS_peaks.narrowPeak \
  --genome hg38.fa \
  --gre MA0113.3.meme --dr3 MA0074.1.meme \
  --tss gencode_tss.bed \
  --out results_stratified
```

Outputs: `peaks_annotated.tsv.gz`, `stratum_counts.tsv`, `enrichment.tsv`, `cooccurrence.tsv`,
`run_params.json`. Thresholds (`--half-width`, `--promoter-bp`, `--proximal-bp`, `--motif-p`,
`--min-overlap-frac`, `--gc-bins`) are listed in `python3 run_stratified.py -h`.

## Data Sources

- ChIP-seq: [ReMAP2022](https://remap.univ-amu.fr/) (hg38)
- Temporal RNA-seq: GSE189984, GSE93735, GSE135130
- CYP24A1 ChIP-seq: ENCODE ENCSR000AKV (BEAS-2B)
- NHANES: CDC public data (2001–2018)

## Requirements

```bash
pip install pandas numpy scipy scikit-learn matplotlib pybedtools
```

## Author

Hiroyuki Nagashima, MD  
Department of Gastroenterology and Hepatology  
Sapporo Medical University School of Medicine  
ketuketu@mac.com

## License

CC0 1.0 Universal — see [LICENSE](LICENSE)
