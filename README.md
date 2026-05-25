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
