# Robustness findings — internal memo

**Date:** 2026-05-30
**Scope:** Re-analysis of the three headline claims using only committed data.
**Reproduce:** `python3 scripts/robustness_check.py` → `results/robustness_summary.csv`

These checks were run before external review to find the weak points first.
Two design choices, if challenged by a reviewer, materially change the
conclusions. They should be addressed before any preprint/journal submission.

---

## (A) ChIP-seq depth confound

VDR/GR occupancy is summed across all ReMAP2022 experiments **without
normalising for the number of experiments**. GR (NR3C1) is far more heavily
ChIP'd than VDR:

| | VDR | GR |
|---|---|---|
| #experiments (median) | 2 | 7 |
| #experiments (max) | 7 | 28 |

Raw score correlates strongly with experiment count: **VDR r=0.46, GR r=0.65**
(both p<1e-14). So a raw "GR-dominant" label partly encodes *study popularity*,
not biology.

Dividing each score by its experiment count (mean signal per experiment)
**flips the VDR/GR dominance direction for 98/381 genes (26%)**. The paradigm
example collapses:

| gene | raw VDR/GR | per-experiment VDR/GR |
|---|---|---|
| **JAK1** (the "GR-dominant black-box" archetype) | 12.8 / 85.5 | **6.4 / 7.8 (≈ parity)** |
| TNF | 36.9 / 19.5 | 36.9 / 4.9 |
| ITGB7 | 62.3 / 17.2 | 31.1 / 2.9 |

## (1) Approval prediction (n=92)

| metric | published (raw, n=92) | raw, matched n=66 | depth-corrected, n=66 |
|---|---|---|---|
| PPV (VDR-dominant) | 0.89 | 0.85 | 0.72 |
| Odds ratio | 8.0 | 3.59 | **1.06** |
| Fisher p | 0.0001 | **0.053** | **1.00** |
| AUC | 0.706 | 0.549 | 0.578 |

Two independent problems:
1. The published p=0.0001 leans on the **26 genes that have no experiment-count
   data**; on the 66 genes with measurable depth the raw signal is already only
   marginal (p=0.053) and AUC ≈ 0.55 (near chance).
2. Depth correction removes the signal entirely (OR→1.06, p→1.0).

## (2) IBD maintenance remission (the r=0.899 result)

- Per drug-disease pair (as published): n=12, Spearman **r=0.899, p=0.00007** —
  reproduces.
- But there are only **3 unique target genes / 3 unique VDR values**
  (TNF=37, ITGB7=62, IL23A=211). n=12 is pseudoreplication; the true number of
  independent x-values is 3.
- Collapsed to one point per gene: n=3, r=1.0 (a monotonic 3-point trend; the
  p-value is not meaningful).
- Leave-one-gene-out: dropping **TNF** → r=0.632, **p=0.13 (non-significant)**.
  A single gene controls the result.

## (A2) Scoring provenance — the table is a patchwork of ≥2 incompatible methods

Tracing how `remap_scores_expanded.csv` is built (`expand_steroid_targets.py`)
revealed the scores are **not** computed by the published Methods:

- The committed generator scores genes as **`10·(n_celltypes) + n_experiments`**
  within a **±5 kb** window — a pure count, with **no peak-height information**.
  (This is why "score ∝ experiment count" in section A is near-tautological.)
- But that integer formula reproduces only **122/381 (VDR)** and **65/381 (GR)**
  committed rows. The other **274/381 rows carry decimal scores** from an
  earlier signal-based method that is **not committed** to the repo.
- The Methods text states yet a **third** definition: "±10 kb … sum of peak
  signal intensities … pybedtools/BEDTools". This matches neither code path.
- TNFSF15's pre-registered scores (VDR=12.4 / GR=286.0) are decimals that fit
  none of the above, so they came from a fourth, ad-hoc computation.

**Consequence:** the 381-gene VDR/GR scores are not on a common scale, no
committed script regenerates the table, and TNFSF15 cannot be "merged in"
consistently. A single canonical pipeline must be defined and run over all
genes at once.

## (3) "TNFSF15 GR=286, highest in the 381-gene dataset"

- **TNFSF15 is not present** in `remap_scores_expanded.csv` (cited as
  Supplementary Table 1, 381 rows). No row has GR≈286.
- Max GR in the committed table is **EGFR=370**; **5 genes exceed GR=286**
  (IL1R1, EGFR, PMEL, MAP2K1, MUC1).
- The pre-registration states TNFSF15 was scored separately on 2026-05-25, so
  this may be a versioning gap — but the manuscript wording ("highest GR score
  in our 381-gene dataset") is not supported by the file it cites, and will be
  flagged on OSF/Supplementary cross-check.

---

## Recommended fixes (priority order)

1. **Re-define the score** to remove the depth confound (per-experiment mean, or
   peak-presence rank), then re-run every downstream statistic. Decide whether
   anything survives before doing more.
2. **Report the IBD relationship at the gene level** (3–6 points) and drop the
   "n=12, p<0.0001" framing; or model drug-disease pairs with a mixed-effects
   model treating target gene as a random effect.
3. **Add a TF negative control** (e.g. CTCF) to show the prediction is specific
   to VDR/GR and not a generic occupancy/expression proxy. *(Requires raw ReMAP
   bed data — not in this repo yet.)*
4. **Reconcile TNFSF15**: add it to the 381-gene table from a single dated
   pipeline run, and correct the "highest GR" wording to "highest among the
   chronic-inflammatory targets" (if true).
5. Narrow scope to the core claim; move NHANES / cross-species / CYP24A1 to a
   separate paper.

---

## Canonical re-scoring pipeline (addresses fix #1 and #4)

`scripts/score_remap_occupancy.py` defines ONE explicit rule and re-scores every
gene + TNFSF15 in a single run. For each gene it reports, within ±W of the TSS
(default ±10 kb, matching the Methods text), both:
- `*_signal_sum` — sum of peak signalValue (the Methods-stated definition), and
- `*_signal_per_exp` — depth-normalised (per-experiment) score,

plus `n_exp` / `n_celltype`, and a factual GR rank for TNFSF15 (replacing the
unverified "highest GR" wording). Output goes to
`results/remap_scores_canonical.csv` and does **not** overwrite the legacy table.

```bash
# Requires the two ReMAP2022 hg38 bed.gz files (not in repo; from remap.univ-amu.fr)
python3 scripts/score_remap_occupancy.py \
    --vdr-bed /path/remap2022_VDR_all_macs2_hg38.bed.gz \
    --gr-bed  /path/remap2022_NR3C1_all_macs2_hg38.bed.gz
```

TNFSF15's TSS (chr9:114,806,039, GRCh38) is seeded in
`results/gene_tss_hg38.csv` so it is always scored identically to the rest;
other genes' TSS are fetched from Ensembl REST and cached on first run.
**Status:** logic validated on synthetic peaks; real numbers require running
where the ReMAP bed files and network are available (both unavailable in the
web sandbox — ReMAP/Ensembl are network-blocked here).
