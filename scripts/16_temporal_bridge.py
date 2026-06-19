#!/usr/bin/env python3
"""
16_temporal_bridge.py
=====================
Decisive, occupancy-independent test of the induction→maintenance BRIDGE
hypothesis for PPARγ:

    GR targets peak EARLY (acute, 0-4 h)
    PPARγ targets peak in the MIDDLE (resolution / bridge, ~8-12 h)
    VDR targets peak LATE (maintenance, 24-48 h)

If PPARγ-module target genes reach their expression peak at times strictly
between the GR module (early) and the VDR module (late), the temporal bridge is
supported independently of ChIP-seq occupancy.

INPUT — a gene × timepoint expression matrix (CSV), e.g. log2 fold-change vs
vehicle/0 h, or normalized counts. One column holds gene symbols; the remaining
numeric columns are timepoints. Hours are taken from --timepoints or parsed from
the column names (first number found, e.g. "t8", "8h", "hr_8" → 8).

    gene,0h,2h,4h,8h,12h,24h,48h
    TSC22D3,0,3.5,3.1,1.2,0.4,0.1,0.0
    CD36,0,0.2,0.8,2.4,2.1,1.0,0.6
    CYP24A1,0,0.1,0.3,1.0,2.0,6.0,10.3
    ...

Ideal dataset: a single inflammation→resolution time course in macrophages
(e.g. LPS or LPS→efferocytosis/IL-4, 0-48 h). With stimulus-specific datasets
(Dex for GR, VitD for VDR) the cross-stimulus caveat must be stated — running
each module on its own dataset still establishes the relative peak ordering but
not within one stimulus.

Usage:
    python scripts/16_temporal_bridge.py --matrix expr_timecourse.csv
    python scripts/16_temporal_bridge.py --matrix m.csv --timepoints 0,2,4,8,12,24,48

Outputs (results/):
    temporal_bridge_peaks.csv     per-gene peak hour + module
    temporal_bridge_stats.txt     module peak-hour summary + ordering test
    temporal_bridge.pdf/.png      peak hour by module (strip + median)
"""
import argparse
import os
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

REPO = Path(__file__).resolve().parent.parent
RESULTS = Path(os.environ.get("VDRGR_RESULTS", REPO / "results"))

# ── gene modules (resolution-focused) ───────────────────────────────────────
# GR: classic immediate-early glucocorticoid targets (acute)
GR_EARLY = ["TSC22D3", "FKBP5", "DUSP1", "ZBTB16", "PER1", "KLF13", "TXNIP", "DDIT4"]
# PPARγ: M2 / efferocytosis / resolution (bridge) — incl. ALOX15 (makes the
# endogenous PPARγ ligand 15-HETE → feed-forward) and MERTK (efferocytosis)
PPARG_BRIDGE = ["CD36", "MRC1", "ANGPTL4", "FABP4", "MERTK", "ALOX15",
                "CD163", "PPARG", "LPL", "PLIN2", "GDF15", "VDR"]
# VDR: canonical late vitamin-D targets (maintenance)
VDR_LATE = ["CYP24A1", "CAMP", "NOD2", "TLR10", "DEFB4A", "IL37", "CD14"]

MODULES = {"GR (acute)": GR_EARLY,
           "PPARγ (bridge)": PPARG_BRIDGE,
           "VDR (maintenance)": VDR_LATE}


def parse_hours(cols, override):
    if override:
        hrs = [float(x) for x in override.split(",")]
        if len(hrs) != len(cols):
            sys.exit(f"--timepoints has {len(hrs)} values but {len(cols)} timepoint columns")
        return hrs
    hrs = []
    for c in cols:
        m = re.search(r"[-+]?\d*\.?\d+", str(c))
        if not m:
            sys.exit(f"cannot parse hour from column '{c}'; pass --timepoints")
        hrs.append(float(m.group()))
    return hrs


def load_matrix(path, gene_col, timepoints):
    df = pd.read_csv(path)
    if gene_col is None:
        gene_col = df.columns[0]
    df[gene_col] = df[gene_col].astype(str).str.upper()
    tcols = [c for c in df.columns if c != gene_col]
    # keep only numeric timepoint columns
    tcols = [c for c in tcols if pd.api.types.is_numeric_dtype(df[c])]
    if len(tcols) < 3:
        sys.exit("need at least 3 numeric timepoint columns")
    hours = parse_hours(tcols, timepoints)
    order = np.argsort(hours)
    tcols = [tcols[i] for i in order]
    hours = [hours[i] for i in order]
    return df, gene_col, tcols, hours


def peak_hour(row, tcols, hours):
    vals = row[tcols].astype(float).values
    if np.all(np.isnan(vals)):
        return np.nan
    return hours[int(np.nanargmax(vals))]


def write_demo_matrix(path):
    """Synthetic gene × timepoint matrix where GR peaks early, PPARγ mid, VDR late.
    Serves as a pipeline sanity-check and a CSV format template for real data."""
    hours = [0, 2, 4, 8, 12, 24, 48]
    peaks = {  # module gene → (peak hour, amplitude)
        "TSC22D3": 2, "FKBP5": 3, "DUSP1": 4, "ZBTB16": 2,        # GR early
        "CD36": 8, "MRC1": 12, "ANGPTL4": 10, "MERTK": 8,
        "ALOX15": 12, "PPARG": 10,                                 # PPARγ mid
        "CYP24A1": 48, "CAMP": 24, "NOD2": 24, "TLR10": 48,        # VDR late
    }
    rng = np.random.default_rng(0)
    rows = []
    for g, pk in peaks.items():
        vals = [round(float(np.exp(-((h - pk) ** 2) / (2 * 6.0 ** 2))
                            + rng.normal(0, 0.03)), 3) for h in hours]
        rows.append({"gene": g, **dict(zip([f"{h}h" for h in hours], vals))})
    df = pd.DataFrame(rows)
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, index=False)
    return str(path)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--matrix", help="gene × timepoint expression CSV")
    ap.add_argument("--gene-col", default=None, help="gene-symbol column (default: first)")
    ap.add_argument("--timepoints", default=None, help="comma hours, e.g. 0,2,4,8,12,24,48")
    ap.add_argument("--out", default=str(RESULTS / "temporal_bridge_peaks.csv"))
    ap.add_argument("--demo", action="store_true",
                    help="generate a synthetic matrix (GR early/PPARγ mid/VDR late) and run "
                         "— sanity-checks the pipeline and writes a CSV format template")
    args = ap.parse_args()

    if args.demo:
        args.matrix = write_demo_matrix(RESULTS / "temporal_bridge_demo_matrix.csv")
        print(f"[demo] synthetic matrix → {args.matrix}\n")
    elif not args.matrix:
        sys.exit("provide --matrix <file.csv> (no angle brackets), or --demo to see a worked example")

    df, gcol, tcols, hours = load_matrix(args.matrix, args.gene_col, args.timepoints)
    print(f"timepoints (h): {hours}")

    g2mod = {g: m for m, genes in MODULES.items() for g in genes}
    sub = df[df[gcol].isin(g2mod)].copy()
    sub["module"] = sub[gcol].map(g2mod)
    sub["peak_h"] = sub.apply(lambda r: peak_hour(r, tcols, hours), axis=1)
    sub = sub.dropna(subset=["peak_h"])
    sub[[gcol, "module", "peak_h"] + tcols].to_csv(args.out, index=False)

    lines = ["=" * 64, "Temporal bridge test: peak hour by regulator module", "=" * 64]
    found = {}
    for m in MODULES:
        ph = sub.loc[sub["module"] == m, "peak_h"]
        found[m] = ph
        genes = ", ".join(sorted(sub.loc[sub["module"] == m, gcol]))
        if len(ph):
            lines.append(f"  {m:20s} n={len(ph):2d}  median peak = {ph.median():5.1f} h  "
                         f"mean = {ph.mean():5.1f} h")
            lines.append(f"      genes: {genes}")
        else:
            lines.append(f"  {m:20s} n=0  (no module genes found in matrix)")

    gr, pg, vd = (found[k] for k in ("GR (acute)", "PPARγ (bridge)", "VDR (maintenance)"))
    lines.append("\nBridge ordering (median peak hours):")
    if len(gr) and len(pg) and len(vd):
        mgr, mpg, mvd = gr.median(), pg.median(), vd.median()
        ok = mgr <= mpg <= mvd
        lines.append(f"  GR {mgr:.1f}  ≤?  PPARγ {mpg:.1f}  ≤?  VDR {mvd:.1f}   "
                     f"→ {'SUPPORTED' if ok else 'NOT supported'}")
        if len(gr) >= 3 and len(pg) >= 3:
            lines.append(f"  PPARγ vs GR  later? Mann-Whitney p="
                         f"{stats.mannwhitneyu(pg, gr, alternative='greater').pvalue:.4f}")
        if len(pg) >= 3 and len(vd) >= 3:
            lines.append(f"  PPARγ vs VDR earlier? Mann-Whitney p="
                         f"{stats.mannwhitneyu(pg, vd, alternative='less').pvalue:.4f}")
        if len(gr) >= 3 and len(pg) >= 3 and len(vd) >= 3:
            kw = stats.kruskal(gr, pg, vd)
            lines.append(f"  Kruskal-Wallis across modules: H={kw.statistic:.2f}, p={kw.pvalue:.4f}")
    else:
        lines.append("  insufficient module coverage for ordering test")

    lines.append("\nCaveat: if GR/PPARγ/VDR modules come from different stimuli or")
    lines.append("cell types, peak ordering is suggestive, not within-stimulus proof.")
    txt = "\n".join(lines)
    print("\n" + txt)
    (RESULTS / "temporal_bridge_stats.txt").write_text(txt + "\n")

    # figure
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(7, 4.5))
        colors = {"GR (acute)": "#8E44AD", "PPARγ (bridge)": "#2E8B57",
                  "VDR (maintenance)": "#E67E22"}
        for i, m in enumerate(MODULES):
            ph = found[m]
            if not len(ph):
                continue
            x = np.random.normal(i, 0.06, len(ph))
            ax.scatter(x, ph, color=colors[m], s=45, alpha=0.8, edgecolor="black", lw=0.4)
            ax.hlines(ph.median(), i - 0.25, i + 0.25, color=colors[m], lw=3)
        ax.set_xticks(range(len(MODULES)))
        ax.set_xticklabels(list(MODULES.keys()))
        ax.set_ylabel("Target-gene expression peak (hours)")
        ax.set_title("Temporal bridge: GR → PPARγ → VDR peak ordering")
        ax.spines[["top", "right"]].set_visible(False)
        plt.tight_layout()
        plt.savefig(RESULTS / "temporal_bridge.pdf")
        plt.savefig(RESULTS / "temporal_bridge.png", dpi=300)
        print(f"\nFigure → {RESULTS / 'temporal_bridge.pdf'}")
    except Exception as e:
        print(f"(figure skipped: {e})", file=sys.stderr)

    print(f"Per-gene peaks → {args.out}")


if __name__ == "__main__":
    main()
