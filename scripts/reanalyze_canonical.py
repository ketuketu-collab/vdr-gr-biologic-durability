#!/usr/bin/env python3
"""
Re-run the headline statistics on the CANONICAL score table
===========================================================

Closes the loop after `score_remap_occupancy.py`. Given
`results/remap_scores_canonical.csv` (one consistent scoring rule for every gene
+ TNFSF15), this regenerates the three headline claims under BOTH score
definitions emitted by the canonical pipeline:

  * signal_sum     — the Methods-stated "sum of peak signal intensities"
  * signal_per_exp — depth-normalised (per-experiment) score

and writes a single before/after verdict table so the effect of removing the
ChIP-seq-depth confound is explicit.

Re-uses the SAME clinical inputs as the published analysis:
  results/multi_disease_drug_list.csv        (n=92 approval pairs)
  results/longterm_remission_corrected.csv   (IBD maintenance remission)

Genes are matched to canonical scores by symbol; drug-disease pairs whose target
gene is missing from the canonical table are dropped (and reported).

USAGE
-----
  # 1) produce the canonical table (needs ReMAP bed files):
  python3 scripts/score_remap_occupancy.py --vdr-bed ... --gr-bed ...
  # 2) re-judge every headline claim from it:
  python3 scripts/reanalyze_canonical.py

Output: results/canonical_reanalysis_summary.csv
"""

import os
import sys
import numpy as np
import pandas as pd
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.normpath(os.path.join(HERE, "..", "results"))


def main(canon_path=None):
    canon_path = canon_path or f"{RESULTS}/remap_scores_canonical.csv"
    if not os.path.exists(canon_path):
        sys.exit(f"ERROR: {canon_path} not found.\n"
                 f"Run scripts/score_remap_occupancy.py first (needs ReMAP bed files).")
    canon = pd.read_csv(canon_path)
    canon["gene"] = canon["gene"].astype(str).str.upper()
    need = {"gene", "VDR_signal_sum", "GR_signal_sum",
            "VDR_signal_per_exp", "GR_signal_per_exp"}
    missing = need - set(canon.columns)
    if missing:
        sys.exit(f"ERROR: canonical table missing columns: {sorted(missing)}")

    # two score definitions: (VDR col, GR col, label)
    DEFS = [("VDR_signal_sum", "GR_signal_sum", "signal_sum"),
            ("VDR_signal_per_exp", "GR_signal_per_exp", "depth_normalised")]
    rows = []

    def g2s(col):
        return dict(zip(canon.gene, canon[col]))

    # ── (1) Approval: Fisher + AUC ────────────────────────────────────────────
    drugs = pd.read_csv(f"{RESULTS}/multi_disease_drug_list.csv")
    drugs["gene"] = drugs.gene.astype(str).str.upper()
    drugs["approved"] = (drugs.status == "approved").astype(int)
    matched = drugs[drugs.gene.isin(set(canon.gene))].copy()
    n_drop = len(drugs) - len(matched)
    print(f"(1) Approval: {len(matched)}/{len(drugs)} pairs matched to canonical "
          f"scores ({n_drop} dropped: target gene absent)")
    try:
        from sklearn.metrics import roc_auc_score
    except Exception:
        roc_auc_score = None
    for vcol, gcol, label in DEFS:
        d = matched.copy()
        d["v"] = d.gene.map(g2s(vcol)); d["g"] = d.gene.map(g2s(gcol))
        d = d.dropna(subset=["v", "g"])
        d["vdr_dom"] = d.v > d.g
        a, b = d[d.vdr_dom], d[~d.vdr_dom]
        tab = [[int(a.approved.sum()), int((1 - a.approved).sum())],
               [int(b.approved.sum()), int((1 - b.approved).sum())]]
        orr, fp = stats.fisher_exact(tab)
        ppv = tab[0][0] / max(sum(tab[0]), 1)
        ratio = d.v / d.g.replace(0, np.nan)
        m = ratio.notna()
        auc = (roc_auc_score(d.approved[m], ratio[m])
               if roc_auc_score and m.any() and d.approved[m].nunique() == 2 else float("nan"))
        print(f"    [{label:16s}] n={len(d)}  PPV={ppv:.2f}  OR={orr:.2f}  "
              f"Fisher p={fp:.4f}  AUC={auc:.3f}")
        rows += [
            dict(claim="approval", definition=label, metric="PPV_vdr_dom", value=round(ppv, 3)),
            dict(claim="approval", definition=label, metric="odds_ratio", value=round(orr, 3)),
            dict(claim="approval", definition=label, metric="fisher_p", value=round(fp, 4)),
            dict(claim="approval", definition=label, metric="AUC", value=round(auc, 3)),
        ]

    # ── (2) IBD maintenance remission ─────────────────────────────────────────
    rem = pd.read_csv(f"{RESULTS}/longterm_remission_corrected.csv")
    ibd = rem[rem.disease.isin(["UC", "CD"])].copy()
    ibd["gene"] = ibd.target_gene.astype(str).str.upper()
    ibd = ibd[ibd.gene.isin(set(canon.gene))].dropna(subset=["maintenance_remission"])
    print(f"(2) IBD: {len(ibd)} pairs / {ibd.gene.nunique()} genes matched")
    for vcol, _, label in DEFS:
        d = ibd.copy()
        d["v"] = d.gene.map(g2s(vcol))
        d = d.dropna(subset=["v"])
        r_pair, p_pair = stats.spearmanr(d.v, d.maintenance_remission)
        gene = d.groupby("gene").agg(v=("v", "first"),
                                     maint=("maintenance_remission", "mean"))
        if len(gene) >= 2:
            r_gene, p_gene = stats.spearmanr(gene.v, gene.maint)
        else:
            r_gene, p_gene = float("nan"), float("nan")
        print(f"    [{label:16s}] per-pair r={r_pair:.3f} (p={p_pair:.4f}, n={len(d)})  "
              f"| per-gene r={r_gene:.3f} (p={p_gene:.4f}, n={len(gene)})")
        rows += [
            dict(claim="IBD", definition=label, metric="spearman_r_per_pair", value=round(r_pair, 3)),
            dict(claim="IBD", definition=label, metric="spearman_r_per_gene", value=round(r_gene, 3)),
        ]

    # ── (3) TNFSF15 GR-rank claim ─────────────────────────────────────────────
    print("(3) TNFSF15 'highest GR' claim:")
    if "TNFSF15" in set(canon.gene):
        t = canon[canon.gene == "TNFSF15"].iloc[0]
        for _, gcol, label in DEFS:
            rank = int((canon[gcol] > t[gcol]).sum()) + 1
            print(f"    [{label:16s}] TNFSF15 GR={t[gcol]:.2f}  rank #{rank}/{len(canon)}"
                  + ("  (= highest)" if rank == 1 else "  (NOT highest)"))
            rows.append(dict(claim="TNFSF15", definition=label,
                             metric="GR_rank", value=rank))
    else:
        print("    TNFSF15 absent from canonical table (run scorer with seed TSS).")

    out = pd.DataFrame(rows)
    out_path = f"{RESULTS}/canonical_reanalysis_summary.csv"
    out.to_csv(out_path, index=False)
    print(f"\nWrote verdict table -> {out_path}")
    print(out.pivot_table(index=["claim", "metric"], columns="definition",
                          values="value", sort=False).to_string())


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else None)
