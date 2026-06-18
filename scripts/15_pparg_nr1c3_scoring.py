#!/usr/bin/env python3
"""
15_pparg_nr1c3_scoring.py
=========================
Add PPARγ (gene PPARG / nuclear receptor NR1C3) ChIP-seq occupancy as a THIRD
regulator axis alongside VDR and GR, then test the hypothesis that PPARγ sits
even further toward the chronic-homeostatic / durable pole than VDR.

Scoring convention is IDENTICAL to the VDR/GR pipeline in
`expand_steroid_targets.py` so the new PPARG_score is directly comparable to the
existing VDR_score / GR_score columns in results/remap_scores_expanded.csv:

    score = (#distinct base cell types within ±HALF of TSS) * 10
            + (#distinct ChIP-seq experiments within ±HALF of TSS)

    HALF = 5000  (±5 kb of TSS; matches the existing expanded scores)

TSS coordinates are fetched from the Ensembl REST API (same batch endpoint as
the original pipeline). The PPARG ReMAP2022 BED is supplied via --pparg-bed or
the REMAP_PPARG_BED env var; download it from https://remap.univ-amu.fr/
(per-TF "All peaks" track for PPARG, hg38) e.g.:
    remap2022_PPARG_all_macs2_hg38.bed.gz

Outputs (under results/):
    remap_scores_pparg.csv        gene, VDR_score, GR_score, PPARG_score, ...
    pparg_durability_stats.txt    Mann-Whitney + Spearman + AUC comparison

Usage:
    python scripts/15_pparg_nr1c3_scoring.py \
        --pparg-bed /path/to/remap2022_PPARG_all_macs2_hg38.bed.gz

If the BED is not found the script still runs the *analysis* portion on any
pre-computed remap_scores_pparg.csv, and otherwise prints clear instructions.

Hypothesis under test (HN, 2026-06):
    GR  →  VDR  →  PPARγ  is a single temporal durability axis.
    Predictions:
      (P1) durable (◎) targets score HIGHER on PPARγ than problematic (⚠️/❌),
           with a separation at least as strong as VDR (Mann-Whitney).
      (P2) PPARγ score correlates with IBD maintenance remission at least as
           strongly as VDR (Spearman; VDR benchmark r = 0.899).
      (P3) adding PPARγ does not collapse the cancer-negative-control specificity.
"""
import argparse
import gzip
import json
import os
import re
import sys
import time
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

try:
    from sklearn.metrics import roc_auc_score
    HAVE_SKLEARN = True
except Exception:
    HAVE_SKLEARN = False

# ── paths (override via env) ────────────────────────────────────────────────
REPO = Path(__file__).resolve().parent.parent
RESULTS = Path(os.environ.get("VDRGR_RESULTS", REPO / "results"))
HALF = int(os.environ.get("REMAP_HALF_WINDOW", "5000"))  # ±5 kb, matches pipeline
ENSEMBL = "https://rest.ensembl.org/lookup/symbol/homo_sapiens"


# ── TSS lookup (Ensembl batch) ──────────────────────────────────────────────
def fetch_tss_batch(genes):
    data = json.dumps({"symbols": genes}).encode()
    req = urllib.request.Request(
        ENSEMBL, data=data,
        headers={"Content-Type": "application/json", "Accept": "application/json"},
    )
    try:
        with urllib.request.urlopen(req, timeout=60) as r:
            d = json.loads(r.read())
    except Exception as e:
        print(f"  Ensembl error: {e}", file=sys.stderr)
        return {}
    tss = {}
    for gene, info in d.items():
        if not info:
            continue
        chrom = str(info.get("seq_region_name", ""))
        if re.match(r"^\d+$|^X$|^Y$", chrom):
            chrom = "chr" + chrom
            strand = info.get("strand", 1)
            tss_pos = info["start"] if strand == 1 else info["end"]
            tss[gene.upper()] = (chrom, tss_pos)
    return tss


def fetch_all_tss(genes):
    tss_map = {}
    for i in range(0, len(genes), 50):
        chunk = genes[i:i + 50]
        tss_map.update(fetch_tss_batch(chunk))
        time.sleep(0.5)
    return tss_map


# ── ReMAP BED indexing + scoring (identical convention to VDR/GR) ───────────
def build_index(bed_gz):
    idx = {}
    with gzip.open(bed_gz, "rt") as f:
        for line in f:
            p = line.rstrip().split("\t")
            if len(p) < 4:
                continue
            chrom, start, end, name = p[0], int(p[1]), int(p[2]), p[3]
            parts = name.split(".")
            gse_id = parts[0] if parts else "unknown"
            celltype = parts[2] if len(parts) > 2 else "unknown"
            base_cell = re.sub(r"_[A-Z0-9_]+$", "", celltype, flags=re.IGNORECASE)
            idx.setdefault(chrom, []).append((start, end, base_cell, gse_id))
    return idx


def score_gene(chrom, tss_pos, idx, n_total):
    lo, hi = tss_pos - HALF, tss_pos + HALF
    hits = [e for e in idx.get(chrom, []) if e[1] >= lo and e[0] <= hi]
    cells = set(h[2] for h in hits)
    exps = set(h[3] for h in hits)
    score = len(cells) * 10 + len(exps)
    repro = len(exps) / n_total if n_total else 0
    return score, len(cells), len(exps), round(repro, 4)


def compute_pparg_scores(bed_path, genes):
    print(f"Building PPARG ReMAP index from {bed_path} ...")
    idx = build_index(bed_path)
    all_exp = {gse for chrom_data in idx.values() for *_, gse in chrom_data}
    n_total = len(all_exp)
    print(f"  PPARG experiments (reproducibility denominator): {n_total}")

    print(f"Fetching TSS for {len(genes)} genes from Ensembl ...")
    tss_map = fetch_all_tss(genes)
    print(f"  TSS resolved: {len(tss_map)}/{len(genes)}")

    rows = []
    for gene, (chrom, tss) in tss_map.items():
        s, c, e, r = score_gene(chrom, tss, idx, n_total)
        rows.append({
            "gene": gene, "PPARG_score": round(s, 2),
            "PPARG_celltypes": c, "PPARG_experiments": e, "PPARG_repro": r,
        })
    return pd.DataFrame(rows)


# ── base gene table (augment expanded set with clinical anchor genes) ───────
def build_base(scores_path):
    """Union of the expanded-score gene list with clinically essential genes
    (IL23A, TNFSF15, NOD2, CYP24A1, ...) that are absent from the 381-gene
    expanded file but carry VDR/GR scores in the longterm remission table on
    the *same* scoring scale. Without this, IL23A et al. are silently dropped
    from the durability analyses."""
    base = pd.read_csv(scores_path)
    base["gene"] = base["gene"].astype(str).str.upper()
    add_frames = []
    # explicit anchors used in the manuscript but absent from the expanded set,
    # so PPARγ gets scored for them even if they carry no VDR/GR here
    anchors = pd.DataFrame({"gene": ["IL23A", "TNFSF15", "NOD2", "CYP24A1"]})
    add_frames.append(anchors)
    lt = RESULTS / "longterm_remission_corrected.csv"
    if lt.exists():
        d = pd.read_csv(lt)
        gcol = "target_gene" if "target_gene" in d else "gene"
        a = d[[gcol, "vdr_score", "gr_score"]].rename(
            columns={gcol: "gene", "vdr_score": "VDR_score", "gr_score": "GR_score"})
        a["gene"] = a["gene"].astype(str).str.upper()
        add_frames.append(a)
    if add_frames:
        extra = pd.concat(add_frames, ignore_index=True).drop_duplicates("gene")
        missing = extra[~extra["gene"].isin(base["gene"])]
        if len(missing):
            print(f"  augmenting with {len(missing)} clinical genes absent from "
                  f"{Path(scores_path).name}: {sorted(missing['gene'])}")
            base = pd.concat([base, missing], ignore_index=True, sort=False)
    return base.drop_duplicates("gene", keep="first")


def partial_spearman(x, y, z):
    """Partial Spearman of x,y controlling for z: rank-transform, residualize x
    and y on z linearly, correlate the residuals. Returns (r, p, n)."""
    d = pd.DataFrame({"x": x, "y": y, "z": z}).dropna()
    if len(d) < 5:
        return float("nan"), float("nan"), len(d)
    rx, ry, rz = (stats.rankdata(d[c]) for c in ("x", "y", "z"))

    def resid(a, b):
        return a - np.polyval(np.polyfit(b, a, 1), b)

    r, p = stats.pearsonr(resid(rx, rz), resid(ry, rz))
    return r, p, len(d)


def add_within_tf_z(df, cols):
    """Within-regulator z-score so absolute magnitudes are comparable across
    VDR/GR/PPARγ despite very different ChIP-seq experiment counts (PPARγ is
    sparse: ~9 experiments vs VDR/GR). Rank-based tests are unaffected; this
    only makes mean-level (P1) comparisons scale-fair."""
    for c in cols:
        if c in df:
            v = df[c].astype(float)
            sd = v.std(ddof=0)
            df[c.replace("_score", "_z")] = (v - v.mean()) / sd if sd else 0.0
    return df


# ── analysis ────────────────────────────────────────────────────────────────
def run_analysis(scored: pd.DataFrame, out_txt: Path):
    lines = []

    def log(s=""):
        print(s)
        lines.append(s)

    add_within_tf_z(scored, ["VDR_score", "GR_score", "PPARG_score"])

    log("=" * 70)
    log("PPARγ (NR1C3) durability-axis analysis")
    log("=" * 70)
    log(f"genes scored: {len(scored)}   (PPARG_score non-null: "
        f"{scored['PPARG_score'].notna().sum()})")
    for col in ("VDR_score", "GR_score", "PPARG_score"):
        if col in scored:
            v = scored[col].dropna()
            log(f"  {col:12s} median={v.median():7.2f}  mean={v.mean():7.2f}  max={v.max():7.2f}")
    n_pp = int(scored["PPARG_experiments"].dropna().sum()) if "PPARG_experiments" in scored else None
    log("  NOTE: PPARγ ChIP-seq is sparse/adipocyte-biased; absolute magnitudes")
    log("  are NOT comparable across regulators. Use _z (within-TF) for P1 and")
    log("  rank-based Spearman (P2) for the durability claim.")

    # (P1) durability tiers: durable (◎) vs problematic (⚠️/❌)
    tier_path = RESULTS / "maintenance_logic_classification.csv"
    if tier_path.exists():
        tiers = pd.read_csv(tier_path)
        gcol = "gene" if "gene" in tiers else "target_gene"
        m = tiers.merge(scored, left_on=gcol, right_on="gene", how="left")
        durable = m[m["lt_status"] == "◎"]
        problem = m[m["lt_status"].isin(["⚠️", "❌"])]
        log("\n--- (P1) durable (◎) vs problematic (⚠️/❌)  [within-TF z-score] ---")
        log(f"  n durable={len(durable)}  n problematic={len(problem)}")
        for raw, z in (("VDR_score", "VDR_z"), ("GR_score", "GR_z"),
                       ("PPARG_score", "PPARG_z")):
            a, b = durable[raw].dropna(), problem[raw].dropna()
            az, bz = durable[z].dropna(), problem[z].dropna()
            if len(a) >= 2 and len(b) >= 2:
                _, p = stats.mannwhitneyu(a, b, alternative="two-sided")
                log(f"  {raw:12s} durable_z={az.mean():+5.2f}  problematic_z={bz.mean():+5.2f}  "
                    f"(raw {a.mean():6.2f} vs {b.mean():6.2f})  MWU p={p:.4f}")
    else:
        log(f"\n(P1) skipped: {tier_path} not found")

    # (P2) PPARγ vs IBD maintenance remission (VDR benchmark r=0.899)
    rem_path = RESULTS / "longterm_remission_corrected.csv"
    if rem_path.exists():
        rem = pd.read_csv(rem_path)
        gcol = "target_gene" if "target_gene" in rem else "gene"
        m = rem.merge(scored[["gene", "PPARG_score"]], left_on=gcol,
                      right_on="gene", how="left")
        log("\n--- (P2) score vs IBD maintenance remission (Spearman) ---")
        if "vdr_score" in m:
            r, p = stats.spearmanr(m["vdr_score"], m["maintenance_remission"],
                                   nan_policy="omit")
            log(f"  VDR   r={r:.3f}  p={p:.4g}  (benchmark)")
        sub = m.dropna(subset=["PPARG_score", "maintenance_remission"])
        if len(sub) >= 4:
            r, p = stats.spearmanr(sub["PPARG_score"], sub["maintenance_remission"])
            log(f"  PPARγ r={r:.3f}  p={p:.4g}  (n={len(sub)})")
        else:
            log("  PPARγ: insufficient overlap for correlation")

        # (P2-IBD) UC+CD only — the apples-to-apples manuscript r=0.899 setting
        ibd = rem[rem["disease"].isin(["UC", "CD"])].copy()
        ibd = ibd.merge(scored[["gene", "PPARG_score"]], left_on=gcol,
                        right_on="gene", how="left")
        log(f"\n--- (P2-IBD) UC+CD only, n={len(ibd)}  (manuscript r=0.899 setting) ---")
        for label, col in (("VDR", "vdr_score"), ("GR", "gr_score"),
                           ("PPARγ", "PPARG_score")):
            sub = ibd.dropna(subset=[col, "maintenance_remission"])
            if len(sub) >= 4:
                r, p = stats.spearmanr(sub[col], sub["maintenance_remission"])
                log(f"  {label:6s} r={r:+.3f}  p={p:.4g}  (n={len(sub)})")
            else:
                log(f"  {label:6s} insufficient overlap (n={len(sub)})")
        miss = ibd[ibd["PPARG_score"].isna()][gcol].unique()
        if len(miss):
            log(f"  (PPARG missing for: {sorted(miss)} — verify these were scored)")
    else:
        log(f"\n(P2) skipped: {rem_path} not found")

    # (P3) approval AUC comparison (VDR vs PPARγ) if approval labels exist
    appr_path = RESULTS / "gene_disease_phase_expanded.csv"
    if appr_path.exists() and HAVE_SKLEARN:
        appr = pd.read_csv(appr_path)
        statuscol = next((c for c in ("best_status", "mol_status", "status")
                          if c in appr), None)
        gcol = next((c for c in ("gene", "target_gene") if c in appr), None)
        log("\n--- (P3) approval discrimination (AUC) ---")
        if statuscol and gcol:
            a = appr[appr[statuscol].isin(["approved", "failed"])].copy()
            a["y"] = (a[statuscol] == "approved").astype(int)
            # approval table already carries VDR_score; bring only PPARG_score
            a = a.merge(scored[["gene", "PPARG_score"]],
                        left_on=gcol, right_on="gene", how="left",
                        suffixes=("", "_s"))
            for col in ("VDR_score", "PPARG_score"):
                s = a.dropna(subset=[col, "y"])
                if s["y"].nunique() == 2 and len(s) >= 10:
                    auc = roc_auc_score(s["y"], s[col])
                    log(f"  {col:12s} AUC={auc:.3f}  (n={len(s)})")
                else:
                    log(f"  {col:12s} insufficient data (n={len(s)})")
        else:
            log(f"  skipped: no usable status/gene column in {appr_path.name}")
    else:
        log("\n(P3) skipped: approval table or sklearn unavailable")

    # (P4) regulator collinearity — tests the "bridge" hypothesis:
    # if PPARγ is intermediate between GR (induction) and VDR (maintenance),
    # it should share variance with BOTH, whereas VDR and GR are near-orthogonal.
    log("\n--- (P4) regulator collinearity (Spearman across genes) ---")

    def _coll(df, label):
        d = df.dropna(subset=["VDR_score", "GR_score", "PPARG_score"])
        pg, pv = (stats.spearmanr(d["PPARG_score"], d[c])[0] for c in ("GR_score", "VDR_score"))
        vg = stats.spearmanr(d["VDR_score"], d["GR_score"])[0]
        log(f"  {label} (n={len(d)}):  PPARγ–GR r={pg:+.3f}   "
            f"PPARγ–VDR r={pv:+.3f}   VDR–GR r={vg:+.3f}")

    _coll(scored, "all genes")
    appr_path = RESULTS / "gene_disease_phase_expanded.csv"
    if appr_path.exists():
        tg = set(pd.read_csv(appr_path)["gene"].astype(str).str.upper())
        _coll(scored[scored["gene"].isin(tg)], "immune drug-targets")
    log("  bridge prediction: PPARγ correlates with BOTH GR and VDR; VDR–GR weak.")

    # (P5) does the PPARγ durability signal survive controlling for GR?
    rem_path = RESULTS / "longterm_remission_corrected.csv"
    if rem_path.exists():
        rem = pd.read_csv(rem_path)
        gcol = "target_gene" if "target_gene" in rem else "gene"
        ibd = rem[rem["disease"].isin(["UC", "CD"])].merge(
            scored[["gene", "PPARG_score"]], left_on=gcol, right_on="gene", how="left")
        log("\n--- (P5) partial Spearman vs IBD maintenance, control = GR (UC+CD) ---")
        for label, col in (("VDR", "vdr_score"), ("PPARγ", "PPARG_score")):
            r, p, n = partial_spearman(ibd[col], ibd["maintenance_remission"], ibd["gr_score"])
            log(f"  {label:6s} r(score, maint | GR) = {r:+.3f}  p={p:.4g}  (n={n})")

    log("\nInterpretation:")
    log("  - PPARγ predicts DURABILITY (P1/P2) but not APPROVAL (P3): maintenance-")
    log("    specific, approval-decoupled.")
    log("  - If P4 shows PPARγ correlated with BOTH GR and VDR, PPARγ is a BRIDGE")
    log("    between induction (GR) and maintenance (VDR) — the resolution phase.")
    out_txt.write_text("\n".join(lines) + "\n")
    print(f"\nStats written → {out_txt}")


# ── main ────────────────────────────────────────────────────────────────────
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pparg-bed", default=os.environ.get("REMAP_PPARG_BED"),
                    help="ReMAP2022 PPARG all-peaks BED(.gz), hg38")
    ap.add_argument("--scores", default=str(RESULTS / "remap_scores_expanded.csv"),
                    help="existing VDR/GR scores (gene list source)")
    ap.add_argument("--out", default=str(RESULTS / "remap_scores_pparg.csv"))
    ap.add_argument("--analysis-only", action="store_true",
                    help="skip scoring; run analysis on existing --out file")
    args = ap.parse_args()

    base = build_base(args.scores)
    genes = sorted(base["gene"].astype(str).str.upper().unique())

    out_path = Path(args.out)
    if args.analysis_only or not args.pparg_bed:
        if out_path.exists():
            scored = pd.read_csv(out_path)
        elif not args.pparg_bed:
            print("ERROR: no --pparg-bed supplied and no pre-computed scores at",
                  out_path, file=sys.stderr)
            print("\nDownload the PPARG track from https://remap.univ-amu.fr/ and run:",
                  file=sys.stderr)
            print("  python scripts/15_pparg_nr1c3_scoring.py "
                  "--pparg-bed remap2022_PPARG_all_macs2_hg38.bed.gz", file=sys.stderr)
            sys.exit(2)
        else:
            scored = pd.read_csv(out_path)
    else:
        bed = Path(args.pparg_bed)
        if not bed.exists():
            print(f"ERROR: PPARG BED not found: {bed}", file=sys.stderr)
            sys.exit(2)
        pparg = compute_pparg_scores(bed, genes)
        scored = base.merge(pparg, on="gene", how="left")
        scored.to_csv(out_path, index=False)
        print(f"Scores written → {out_path}  ({scored['PPARG_score'].notna().sum()} genes scored)")

    run_analysis(scored, RESULTS / "pparg_durability_stats.txt")


if __name__ == "__main__":
    main()
