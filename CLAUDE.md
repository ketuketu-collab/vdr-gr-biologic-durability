# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this repository is

This is the analysis code and results snapshot backing a scientific manuscript:
**"VDR/GR ChIP-seq Occupancy Predicts Long-term Durability of Biologics in Chronic Inflammatory Diseases"** (target: *Nature Medicine*; OSF preregistration https://osf.io/tnp63).

It is **not an application** — there is no build system, no test suite, no linter, no package manifest, and no CLI/argument parsing. It is a collection of ~50 standalone Python analysis/figure scripts plus two bash ChIP-seq pipeline scripts. Each script runs top-to-bottom and writes CSV tables and/or PDF+PNG figures.

The central quantity is the **VDR/GR score**: for a drug's target gene, the ChIP-seq occupancy (from ReMap2022, hg38) by the vitamin D receptor (VDR) vs the glucocorticoid receptor (GR/NR3C1) at the gene's TSS ± 5 kb. The thesis is that a high VDR-to-GR ratio predicts durable biologic efficacy (approval rate, IBD maintenance remission, durability tier), while GR-dominant targets predict acute-but-not-durable response.

## Critical convention: scripts use hardcoded local macOS paths

**Every script in `scripts/` reads from and writes to absolute paths on the original author's machine, not repo-relative paths.** You will see, in essentially all 55 scripts:

- `/Volumes/M4_SSD/projects/tlr_chipseq/results/` — primary results dir (76 references)
- `/Volumes/M4_SSD/projects/tlr10_chipseq_vdr_gr/data/` — raw ChIP-seq / ReMap2022 BED inputs (28 references)
- `/Volumes/M4_SSD/ref/gencode.v43.primary_assembly.annotation.gtf.gz` — gene annotation

The git repo's `results/` directory is a **curated snapshot of outputs**, copied back from that local `results/` dir. It is not where the scripts write. **A script will not run in this checkout as-is** — you must repoint its `BASE` / `RESULTS` / `DATA` / `CSV` path constants (defined near the top of each file) to the local `results/` and to wherever the large inputs live. When asked to run or modify a script, expect to update these path constants first.

## Inputs that are NOT in the repo

Large genomic files are gitignored (see `.gitignore`: `*.bed`, `*.bam`, `*.bigwig`, `*.fastq`, `*.pkl`, `data/raw/`). Downstream scripts depend on inputs that exist only on the author's machine:

- **ReMap2022 BED files** — e.g. `remap2022_VDR_all_macs2_hg38.bed.gz`, `remap2022_NR3C1_all_macs2_hg38.bed.gz` (the raw binding-site source for all scoring).
- **`tableS4_phase1_vdr_gr_scores.xlsx`** — drug-approval status, read by several scripts. Note it has a **Japanese sheet name (`全体 VDR降順`)** and Japanese status strings parsed in code: `承認済` → approved, `失敗`/`中止` → failed, `Phase…` → ongoing. Preserve these literals when editing.
- ChIP-seq FASTQ/BAM/bigWig for the alignment pipeline.

What IS tracked (and what `.gitignore` explicitly re-includes): `results/*.csv`, `results/figures/*.pdf|*.png`, `scripts/*.py`, `manuscript/*.md`.

## Two distinct pipelines

**1. ChIP-seq read processing (bash, conda).** `scripts/run_vdr_pipeline.sh` orchestrates `06_align_vdr.sh` (bowtie2 align → samtools sort/index → `bamCoverage` CPM-normalized bigWig → MACS2 peak calling on VitD-vs-vehicle) and `07_bigwig_pileup.py`. Requires a conda env named `chipseq` (`bowtie2`, `samtools`, `macs2`, deepTools) and a bowtie2 hg38 index at `~/genomes/GRCh38_noalt_as/`. Scripts are idempotent via `step_done`/`step_skip` marker files and "skip if output exists" guards. This stage produces the VDR ChIP-seq tracks; most users won't rerun it.

**2. Downstream analysis (Python).** The bulk of the repo. Each script loads CSVs from `results/` (and sometimes the ReMap BEDs / Excel), computes statistics, and emits an updated CSV and/or a manuscript figure. These are independent — there is no shared library or import graph between them; conventions are copied script-to-script. Representative flow:
- `remap_breadth_analysis.py` — fetches TSS coords (Ensembl REST, batched by 50), overlaps TSS±5kb against the ReMap BEDs, computes per-gene VDR/GR cell-type breadth and reproducibility → `remap_breadth_per_gene.csv`.
- `logistic_vdr_approval.py` — logistic regression (VDR/GR scores + reproducibility → approval), bootstrap CIs, CV ROC-AUC → `vdr_logistic_scored_genes.csv` + figure.
- `fig3_durability_composite.py`, `fig5_mechanistic_pipeline.py`, `fig2_cell_paper.py`, etc. — final manuscript figures.

Script naming: numbered prefixes (`01_`…`14_`) are an early ordered ChIP-seq/locus exploration series; later analyses use descriptive names (`fig*`, `logistic_*`, `*_vdr_gr*`, `disease_vdr_map_v2/v3`). Versioned names (`_v2`, `_v3`, `_expanded`, `_corrected`, `_honest`) supersede earlier ones — prefer the highest version / "corrected" variant when one exists.

## Running things

There is no entry point or runner. Run a single Python analysis directly (after fixing its path constants):

```bash
python scripts/logistic_vdr_approval.py
```

Dependencies (no lockfile; install manually):

```bash
pip install pandas numpy scipy scikit-learn matplotlib pybedtools openpyxl
# ChIP-seq pipeline only (via conda): bowtie2 samtools macs2 deeptools
```

## Figure conventions (follow these when adding/editing figures)

- Headless backend: `import matplotlib; matplotlib.use("Agg")` before `pyplot`, and `warnings.filterwarnings("ignore")`.
- Always save **both** `.pdf` and `.png` at `dpi=300`, `bbox_inches="tight"`, into the local `results/figures/`.
- Color semantics are consistent: **VDR = blue** (`#1976D2`/`#2166AC`), **GR = red** (`#D32F2F`/`#D6604D`), approved = green, failed = dark red, with orange for the prospective tulisokibart prediction.
- Panel titles use a left-aligned letter prefix, e.g. `ax.set_title("a  …", loc="left", fontweight="bold")`; top/right spines hidden.
- Durability tiers are encoded with the literal glyphs `◎ △ ⚠️ ❌` (durable / conditional / black-box / failed) — keep these exact characters.

## Key results files

- `results/remap_scores_expanded.csv` — per-gene VDR/GR scores (the 381-gene scored set).
- `results/remap_breadth_per_gene.csv` — adds cell-type breadth + reproducibility (input to the logistic model).
- `results/longterm_remission_corrected.csv` — clinical IBD remission data (use this over `longterm_remission_analysis.csv`).
- `results/maintenance_logic_classification.csv` — durability-tier labels (`lt_status` column).
- `results/pipeline_prediction.csv` / `prospective_predictions.csv` — forward predictions, incl. the preregistered tulisokibart (anti-TL1A) ATLAS-UC failure call.

## Manuscript

`manuscript/manuscript_draft_v2.md` is the live draft; `manuscript/osf_preregistration.md` is the registered preregistration. Statistics quoted in the manuscript (PPV, Fisher/Spearman/Mann-Whitney p-values, AUCs, specific log2FC values) are produced by the scripts and are sometimes also hardcoded into figure annotations — if you change an analysis, check both the figure script's inline numbers and the manuscript text for consistency.

## Git workflow

Active development branch: `claude/claude-md-docs-FaIYb`. Develop, commit, and push there; do not push to `main` without explicit permission. Use `git push -u origin <branch>`, retrying on network errors. Do not open a pull request unless explicitly asked.
