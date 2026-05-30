"""
Robustness / sensitivity analysis for the VDR/GR rule
=====================================================

Re-evaluates the three headline claims using ONLY the committed data, and
quantifies how much they depend on two design choices that a reviewer will
challenge:

  (A) ChIP-seq depth confound — VDR/GR occupancy scores are summed across all
      ReMAP2022 experiments WITHOUT normalising for the number of experiments.
      GR (NR3C1) is far more heavily ChIP'd than VDR, so a raw "GR-dominant"
      label may partly reflect study popularity rather than biology.

  (B) Pseudoreplication — the IBD maintenance-remission correlation (r=0.899,
      n=12) is built from drug-disease pairs that share identical per-gene VDR
      scores. The number of INDEPENDENT x-values is the number of target genes
      (3), not 12.

For each claim we report the original ("raw") result and a depth-corrected /
gene-level result side by side. Output: results/robustness_summary.csv.

Run from the repository root:  python3 scripts/robustness_check.py
"""

import os
import numpy as np
import pandas as pd
from scipy import stats

# Resolve paths relative to this file so the script runs inside the repo.
HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.normpath(os.path.join(HERE, "..", "results"))

rows = []   # collected for the summary table


def rec(section, metric, raw, corrected, note=""):
    rows.append(dict(section=section, metric=metric,
                     raw=raw, corrected=corrected, note=note))


# ──────────────────────────────────────────────────────────────────────────────
# Load scores + per-factor experiment counts
# ──────────────────────────────────────────────────────────────────────────────
scores = pd.read_csv(f"{RESULTS}/remap_scores_expanded.csv")

# Depth-corrected score = mean peak signal per experiment.
for tf in ["VDR", "GR"]:
    scores[f"{tf}_norm"] = scores[f"{tf}_score"] / scores[f"{tf}_experiments"].replace(0, np.nan)

print("=" * 72)
print("(A) ChIP-seq depth confound")
print("=" * 72)
for tf in ["VDR", "GR"]:
    m = scores[f"{tf}_experiments"] > 0
    r, p = stats.spearmanr(scores.loc[m, f"{tf}_score"], scores.loc[m, f"{tf}_experiments"])
    print(f"  {tf} score vs #experiments : Spearman r={r:.3f}  p={p:.2e}  (n={int(m.sum())})")
    rec("A_depth_confound", f"{tf}_score_vs_nexp_spearman_r", round(r, 3), "", f"p={p:.2e}")
print(f"  #experiments  median: VDR={scores.VDR_experiments.median():.0f} "
      f"GR={scores.GR_experiments.median():.0f}   "
      f"max: VDR={scores.VDR_experiments.max():.0f} GR={scores.GR_experiments.max():.0f}")

# How many genes flip VDR/GR dominance after depth correction?
raw_dom = scores.VDR_score > scores.GR_score
norm_dom = scores.VDR_norm.fillna(0) > scores.GR_norm.fillna(0)
flip = (raw_dom != norm_dom)
print(f"  Dominance direction flips after normalisation: {int(flip.sum())}/{len(scores)} "
      f"({100 * flip.mean():.0f}%)")
rec("A_depth_confound", "dominance_flip_fraction", "", f"{100*flip.mean():.0f}%",
    f"{int(flip.sum())}/{len(scores)} genes flip VDR/GR direction")

for g in ["JAK1", "IL6R", "TNF", "ITGB7", "IL23A"]:
    r = scores[scores.gene == g]
    if len(r):
        r = r.iloc[0]
        print(f"    {g:6s} raw VDR/GR={r.VDR_score:6.1f}/{r.GR_score:6.1f}  "
              f"-> norm={r.VDR_norm:5.1f}/{r.GR_norm:5.1f}")


# ──────────────────────────────────────────────────────────────────────────────
# (1) Approval analysis — Fisher + AUC, raw vs depth-corrected
# ──────────────────────────────────────────────────────────────────────────────
print("\n" + "=" * 72)
print("(1) Approval (n=92): raw vs depth-corrected VDR/GR")
print("=" * 72)
drugs = pd.read_csv(f"{RESULTS}/multi_disease_drug_list.csv")
drugs = drugs.merge(
    scores[["gene", "VDR_norm", "GR_norm"]], on="gene", how="left")
drugs["approved"] = (drugs.status == "approved").astype(int)

# Apples-to-apples: restrict BOTH analyses to the rows where a depth-corrected
# score exists, so the raw-vs-corrected difference is not confounded by n.
matched = drugs.dropna(subset=["VDR_norm", "GR_norm"]).copy()
print(f"  (matched subset with experiment counts: n={len(matched)} of {len(drugs)})")


def approval_block(df, vcol, gcol, label):
    d = df.dropna(subset=[vcol, gcol]).copy()
    d["vdr_dom"] = d[vcol] > d[gcol]
    a = d[d.vdr_dom]
    b = d[~d.vdr_dom]
    # Fisher 2x2: dominance x approved
    tab = [[int(a.approved.sum()), int((1 - a.approved).sum())],
           [int(b.approved.sum()), int((1 - b.approved).sum())]]
    orr, fp = stats.fisher_exact(tab)
    ppv = tab[0][0] / max(sum(tab[0]), 1)
    # AUC of continuous ratio
    ratio = d[vcol] / d[gcol].replace(0, np.nan)
    mask = ratio.notna()
    try:
        from sklearn.metrics import roc_auc_score
        auc = roc_auc_score(d.approved[mask], ratio[mask])
    except Exception:
        auc = float("nan")
    print(f"  [{label}] n={len(d)}  PPV(VDR-dom)={ppv:.2f}  OR={orr:.2f}  "
          f"Fisher p={fp:.4f}  AUC={auc:.3f}")
    return ppv, orr, fp, auc, len(d)


ppv_r, or_r, fp_r, auc_r, n_r = approval_block(matched, "vdr_score", "gr_score", "raw (matched n)")
ppv_n, or_n, fp_n, auc_n, n_n = approval_block(matched, "VDR_norm", "GR_norm", "depth-corrected")
rec("1_approval", "PPV_vdr_dominant", round(ppv_r, 2), round(ppv_n, 2), f"matched subset n={n_r}")
rec("1_approval", "fisher_p", round(fp_r, 4), round(fp_n, 4), "")
rec("1_approval", "AUC_vdr_gr_ratio", round(auc_r, 3), round(auc_n, 3), "")


# ──────────────────────────────────────────────────────────────────────────────
# (2) IBD durability — raw n=12 vs gene-level, with leave-one-gene-out
# ──────────────────────────────────────────────────────────────────────────────
print("\n" + "=" * 72)
print("(2) IBD maintenance remission vs VDR score")
print("=" * 72)
rem = pd.read_csv(f"{RESULTS}/longterm_remission_corrected.csv")
ibd = rem[rem.disease.isin(["UC", "CD"])].dropna(
    subset=["vdr_score", "maintenance_remission"]).copy()

r_pair, p_pair = stats.spearmanr(ibd.vdr_score, ibd.maintenance_remission)
n_genes = ibd.target_gene.nunique()
print(f"  Per drug-disease PAIR (as published): n={len(ibd)}  "
      f"Spearman r={r_pair:.3f}  p={p_pair:.5f}")
print(f"  ... but only {n_genes} unique target genes / {ibd.vdr_score.nunique()} "
      f"unique VDR values -> pseudoreplication")

# Collapse to one point per gene (true independent units)
gene = ibd.groupby("target_gene").agg(
    vdr=("vdr_score", "first"),
    maint=("maintenance_remission", "mean")).reset_index()
r_gene, p_gene = stats.spearmanr(gene.vdr, gene.maint)
print(f"  Per GENE (collapsed): n={len(gene)}  Spearman r={r_gene:.3f}  p={p_gene:.4f}")
rec("2_IBD", "spearman_r", round(r_pair, 3), round(r_gene, 3),
    f"raw=per-pair n={len(ibd)} (p={p_pair:.5f}); corrected=per-gene n={len(gene)} (p={p_gene:.4f})")

print("  Leave-one-gene-out (per-pair correlation):")
for g in ibd.target_gene.unique():
    sub = ibd[ibd.target_gene != g]
    if sub.vdr_score.nunique() < 2:
        print(f"    drop {g:8s}: <2 genes remain, undefined")
        continue
    r, p = stats.spearmanr(sub.vdr_score, sub.maintenance_remission)
    print(f"    drop {g:8s} (n={len(sub)}): r={r:.3f}  p={p:.4f}")


# ──────────────────────────────────────────────────────────────────────────────
# (3) TNFSF15 provenance check
# ──────────────────────────────────────────────────────────────────────────────
print("\n" + "=" * 72)
print("(3) TNFSF15 / 'highest GR in 381-gene dataset' claim")
print("=" * 72)
in_table = "TNFSF15" in set(scores.gene)
gmax = scores.GR_score.max()
gmax_gene = scores.loc[scores.GR_score.idxmax(), "gene"]
n_above_286 = int((scores.GR_score >= 286).sum())
print(f"  TNFSF15 present in remap_scores_expanded.csv (Suppl. Table 1)? {in_table}")
print(f"  Max GR in committed table: {gmax:.1f} ({gmax_gene})")
print(f"  Genes with GR>=286 in the table: {n_above_286}  "
      f"-> claim 'GR=286 is the highest' is not supported by this file")
rec("3_TNFSF15", "in_381_table", in_table, "", "")
rec("3_TNFSF15", "n_genes_GR_ge_286", n_above_286, "",
    f"max GR={gmax:.1f} ({gmax_gene}); claimed TNFSF15 GR=286 absent from file")


# ──────────────────────────────────────────────────────────────────────────────
# Write summary table
# ──────────────────────────────────────────────────────────────────────────────
out = pd.DataFrame(rows)
out_path = f"{RESULTS}/robustness_summary.csv"
out.to_csv(out_path, index=False)
print("\n" + "=" * 72)
print(f"Summary written to {out_path}")
print("=" * 72)
print(out.to_string(index=False))
