# NEXT SESSION — resume instructions

Goal: produce the **real** canonical VDR/GR score table (re-scoring every gene +
TNFSF15 under one consistent rule) and the final verdict table for the three
headline claims. The scripts are already written, committed, and validated on
synthetic data. Only the actual data run remains — it needs network access that
was blocked in the session where these scripts were created.

Work branch: `claude/directory-file-listing-5aQNt`

---

## Prerequisite (one-time, done by the human before the session)

In this repo's Claude-Code-on-the-web **Environment → Network policy**, switch to
a custom allowlist and add these hosts, then start a NEW session:

| host | why |
|---|---|
| `remap.univ-amu.fr` | download ReMAP2022 bed files |
| `remapdata.univ-amu.fr` | ReMAP storage mirror (whichever serves the files) |
| `rest.ensembl.org` | per-gene TSS coordinates |

Sanity check at session start:
```bash
curl -sI --max-time 8 https://rest.ensembl.org/info/ping?content-type=application/json | head -1
curl -sI --max-time 8 https://remap.univ-amu.fr/ | head -1
```
Both should return a non-403 status. If still 403, the policy did not take
effect — re-check the environment config and start another new session.

---

## Run (what Claude should execute)

```bash
# 0) deps
pip install -q pandas numpy scipy scikit-learn

# 1) fetch the two ReMAP2022 per-TF bed files (hg38, "all" MACS2 peaks).
#    Get the exact current URLs from https://remap.univ-amu.fr/download_page
#    (TF = VDR and NR3C1). Save them as below:
mkdir -p data
#   VDR  -> data/remap2022_VDR_all_macs2_hg38.bed.gz
#   NR3C1-> data/remap2022_NR3C1_all_macs2_hg38.bed.gz
# Example (verify the path on the download page first):
# curl -L -o data/remap2022_VDR_all_macs2_hg38.bed.gz \
#   "https://remap.univ-amu.fr/storage/remap2022/hg38/MACS2/TF/remap2022_VDR_all_macs2_hg38_v1_0.bed.gz"
# curl -L -o data/remap2022_NR3C1_all_macs2_hg38.bed.gz \
#   "https://remap.univ-amu.fr/storage/remap2022/hg38/MACS2/TF/remap2022_NR3C1_all_macs2_hg38_v1_0.bed.gz"

# 2) canonical re-scoring (every gene + TNFSF15, one rule, +-10kb signal sum
#    plus depth-normalised per-experiment score). TSS auto-fetched from Ensembl
#    and cached to results/gene_tss_hg38.csv (TNFSF15 already seeded).
python3 scripts/score_remap_occupancy.py \
    --vdr-bed data/remap2022_VDR_all_macs2_hg38.bed.gz \
    --gr-bed  data/remap2022_NR3C1_all_macs2_hg38.bed.gz
# -> results/remap_scores_canonical.csv

# 3) re-judge the three headline claims from the canonical table
python3 scripts/reanalyze_canonical.py
# -> results/canonical_reanalysis_summary.csv

# 4) commit the REAL outputs (do NOT commit anything from a synthetic run)
git add results/remap_scores_canonical.csv results/canonical_reanalysis_summary.csv \
        results/gene_tss_hg38.csv
git commit -m "Add canonical VDR/GR scores and re-analysis (real ReMAP run)"
git push -u origin claude/directory-file-listing-5aQNt
```

Do NOT add `data/*.bed.gz` to git (large; keep them out — add to .gitignore if
needed).

---

## What to report back

From `results/canonical_reanalysis_summary.csv`, state for each claim how it
holds up under `signal_sum` vs `depth_normalised`:

1. **Approval** — PPV / odds ratio / Fisher p / AUC. Does the signal survive
   depth-normalisation, or collapse (as the crude check suggested: OR→~1, p→~1)?
2. **IBD maintenance remission** — report BOTH per-pair and per-gene Spearman.
   The per-gene n is only ~3, so treat per-pair p-values as pseudoreplicated.
3. **TNFSF15 GR rank** — the script prints TNFSF15's actual GR rank. Use it to
   replace the manuscript wording "highest GR score in our 381-gene dataset"
   with the true rank (the committed legacy table already shows ≥5 genes above
   GR=286, so "highest" is almost certainly wrong).

Then advise on manuscript edits (Methods must describe the ACTUAL scoring rule:
window size, signal-sum vs count, depth handling).

---

## Background (why this work exists)

See `manuscript/robustness_findings.md` for the full audit. Key points:
- The legacy `results/remap_scores_expanded.csv` mixes ≥2 incompatible scoring
  methods; no committed script regenerates it; the Methods text matches neither.
- Raw VDR/GR scores correlate strongly with ChIP-seq experiment count
  (study-depth confound); depth-normalisation flipped dominance for 26% of genes
  and collapsed the approval signal in a crude check.
- The IBD r=0.899 (n=12) rests on only 3 unique target genes.
- TNFSF15 (the pre-registered prediction target) was absent from the table.

Scripts involved:
- `scripts/score_remap_occupancy.py` — canonical single-rule scorer (THIS run)
- `scripts/reanalyze_canonical.py` — re-judges headline claims from canonical table
- `scripts/robustness_check.py` — audit of the legacy table (already run/committed)
