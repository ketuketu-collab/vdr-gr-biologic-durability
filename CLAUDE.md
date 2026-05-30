# CLAUDE.md

Guidance for AI assistants (Claude Code) working in this repository.

## What this repository is

This is a **computational biology / bioinformatics research project**, not a
software product. It contains the analysis code, derived data tables, figures,
and manuscript for a single scientific paper:

> **"VDR/GR ChIP-seq Occupancy Predicts Long-term Therapeutic Durability of
> Biologics in Chronic Inflammatory Diseases"**
> Target journal: *Nature Medicine* · OSF preregistration: https://osf.io/tnp63
> Author: Hiroyuki Nagashima, MD (Sapporo Medical University)

**Core scientific claim ("the VDR/GR rule"):** the ratio of chromatin occupancy
by the vitamin D receptor (VDR) versus the glucocorticoid receptor (GR / NR3C1)
at a drug's target gene — computed from public ChIP-seq data (ReMAP2022, hg38) —
predicts biologic approval, durable remission, and Phase III outcomes for chronic
inflammatory diseases. The headline prospective prediction is that tulisokibart
(anti-TL1A/TNFSF15, GR-dominant) will fail ATLAS-UC Phase III maintenance.

Because the work is **preregistered and embargoed**, treat the scientific
content as a fixed record. See "Scientific integrity" below.

## Repository layout

```
scripts/      Python (and a few bash) analysis + figure-generation scripts
results/      Derived output tables (*.csv) — tracked in git
results/figures/  Generated figures (*.pdf + *.png) — tracked in git
manuscript/   manuscript_draft_v2.md  (Nature Medicine draft)
              osf_preregistration.md  (preregistration document)
README.md     Public-facing project summary, data sources, requirements
.gitignore    Excludes large/raw bioinformatics data; whitelists tracked outputs
```

There is **no build system, package manifest, test suite, or CI**. The repo is a
collection of standalone scripts plus their committed outputs.

### What is and isn't tracked

`.gitignore` deliberately excludes raw/large genomics data (`*.bed`, `*.bam`,
`*.bigwig`, `*.fastq`, `*.pkl`, `data/raw/`, `logs/`) and then explicitly
re-includes the durable artifacts: `results/*.csv`, `scripts/*.py`,
`manuscript/*.md`, and `results/figures/*.{pdf,png}`. When adding new outputs,
keep to this pattern — commit small derived tables and figures, never raw
sequencing data.

## How the scripts work

- **Language:** Python 3 (`#!/usr/bin/env python3`), plus a couple of bash
  scripts (`06_align_vdr.sh`, `run_vdr_pipeline.sh`) for the upstream ChIP-seq
  alignment that is normally run once on raw FASTQ.
- **Core libraries:** `pandas`, `numpy`, `scipy`, `scikit-learn`, `matplotlib`,
  and for genomics `pybedtools` / `pyBigWig`. Plotting scripts always set the
  non-interactive backend first:
  ```python
  import matplotlib
  matplotlib.use("Agg")
  ```
- **Numbered scripts (`01_`–`14_`)** form a rough upstream→downstream pipeline
  (ReMap peak extraction → peak analysis → locus/comparison figures →
  supplement tables → disease scans). **Named scripts** are later, more targeted
  analyses (ChEMBL validation, disease maps, logistic models, specific figures).
  `run_vdr_pipeline.sh` chains the alignment + pileup steps.
- **Scripts are run individually**, e.g. `python scripts/fig3_durability_composite.py`.
  Each script reads CSVs from `results/`, writes new CSVs and/or figures back
  into `results/` and `results/figures/`.

### ⚠️ Hardcoded absolute paths — the most important gotcha

Essentially **every script hardcodes the original author's local macOS path**:

```python
BASE = "/Volumes/M4_SSD/projects/tlr_chipseq/results"
```

These paths (`/Volumes/M4_SSD/...`, `~/miniforge3`, `~/genomes`) **do not exist
in this checkout or in cloud environments.** Scripts will fail with
`FileNotFoundError` if run as-is. The committed `results/` directory mirrors that
`results/` folder, so to actually run a script you must repoint `BASE` (and any
`RES_DIR` / `FIG_DIR` / raw-data paths) at the repo's own directories.

When asked to make scripts runnable/portable, the right fix is to replace these
hardcoded paths with repo-relative paths derived from the script location, e.g.:

```python
from pathlib import Path
BASE = Path(__file__).resolve().parent.parent / "results"
```

Do **not** do this silently across all 55 scripts as a drive-by; confirm scope
with the user first, since it changes a large number of files at once.

### External data dependencies

- **Local raw data** that is gitignored and absent here: ReMap2022 BED files
  (`remap2022_VDR_all_macs2_hg38.bed.gz`, `..._NR3C1_...`), BAM/bigWig coverage.
- **Network APIs** that some scripts call directly: ChEMBL
  (`https://www.ebi.ac.uk/chembl/api/data`), ClinicalTrials.gov v2, NCBI E-utils,
  and GEO FTP. These require outbound network access and may be rate-limited.

## Coding conventions (match the existing style)

- Each script opens with a module docstring summarizing what it produces (often
  naming the manuscript figure/table it generates, e.g. "Fig. 3 — ...").
- Section separators use box-drawing comment banners:
  `# ── Load data ─────────────────────────`.
- Constants (color palettes, gene coordinates, tier symbols) are defined in
  `UPPER_CASE` near the top. Figures use a consistent palette (VDR blue
  `#1976D2`, GR red `#D32F2F`) and durability-tier glyphs (`◎ △ ⚠️ ❌`).
- `warnings.filterwarnings("ignore")` is used throughout to keep output clean.
- Gene/disease/score data flows through pandas DataFrames keyed by `gene`;
  VDR/GR values live in columns like `VDR_score`, `GR_score`, `vdr_score`,
  `gr_score` (casing varies between files — check the specific CSV).
- Code comments frequently cite data provenance inline (GSE accession numbers,
  published log2FC values, hg38 coordinates). Preserve these citations when
  editing.

## Scientific integrity (important)

- This is preregistered, embargoed research tied to a real manuscript and a real
  prospective clinical prediction. **Do not invent, "clean up," or alter
  numerical results, statistics, scores, or claims** to make things look better.
  Reported numbers (p-values, AUC, correlation coefficients, VDR/GR scores) must
  trace back to the data and code.
- If you find a genuine bug that changes a result, **surface it to the user**
  rather than quietly changing committed outputs or the manuscript.
- Manuscript and preregistration text (`manuscript/*.md`) is a record of what was
  registered. Edit only when explicitly asked, and flag any change that affects a
  registered hypothesis or reported finding.

## Git workflow

- Active development branch for this work: **`claude/claude-md-docs-5OshZ`**.
  Develop, commit, and push there; create it locally if needed. Do not push to
  `main` without explicit permission.
- Push with `git push -u origin <branch>`; retry transient network failures with
  exponential backoff.
- Do **not** open a pull request unless the user explicitly asks.
- Commit messages are short and descriptive (see history, e.g. "Add sequential
  GR→VDR paradigm to Discussion and Fig. 5c"). Keep that style.

## Environment / requirements

```bash
pip install pandas numpy scipy scikit-learn matplotlib pybedtools pyBigWig openpyxl
```

The upstream alignment steps (`*.sh`) additionally assume a conda env named
`chipseq` plus bowtie2/MACS2 and an hg38 index — these only matter if
reprocessing raw ChIP-seq, which is not part of normal work in this checkout.
