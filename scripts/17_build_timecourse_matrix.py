#!/usr/bin/env python3
"""
17_build_timecourse_matrix.py
=============================
Assemble a gene × timepoint expression matrix (the input for
scripts/16_temporal_bridge.py) from common processed-RNA-seq layouts, e.g.
GSE189984 (1,25-D3 PBMC, 4/8/24/48 h).

Three input modes:

(A) --deseq : one DESeq2 result file PER TIMEPOINT (log2FC vs vehicle), the
    layout matching the manuscript methods. Map each file to its hour:
        --deseq vd_4h.csv=4 vd_8h.csv=8 vd_24h.csv=24 vd_48h.csv=48
    Gene + log2FoldChange columns are auto-detected; baseline 0 h is set to 0.

(B) --counts + --design : one wide normalized-count matrix (genes × samples)
    plus a design table mapping sample→hour. Per-gene log2 fold-change is
    computed vs the earliest hour:  log2((mean_t + 1) / (mean_0 + 1)).
        --counts counts.csv --design design.csv
    design.csv: columns  sample,hour   (sample names must match count columns)

(C) --long : a long table with columns gene,hour,value (value = log2FC or expr).
        --long table.csv

Output: a gene × hour CSV (columns "0h","4h",...) ready for scripts/16:
    python scripts/17_build_timecourse_matrix.py --deseq ... --out m.csv
    python scripts/16_temporal_bridge.py --matrix m.csv

Use --genes-only to keep just the GR/PPARγ/VDR module genes (smaller file).
"""
import argparse
import os
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parent.parent
RESULTS = Path(os.environ.get("VDRGR_RESULTS", REPO / "results"))

GENE_CANDIDATES = ["gene", "Gene", "symbol", "Symbol", "gene_name", "gene_symbol",
                   "GeneSymbol", "hgnc_symbol", "external_gene_name", "geneid", "Geneid"]
LFC_CANDIDATES = ["log2FoldChange", "log2FC", "log2fc", "logFC", "lfc",
                  "log2_fold_change", "Log2FoldChange"]


def _pick(cols, candidates):
    for c in candidates:
        if c in cols:
            return c
    return None


def _read_table(path, **kw):
    """Read CSV/TSV/Excel by extension (xlsx needs openpyxl)."""
    p = str(path).lower()
    if p.endswith((".xlsx", ".xls")):
        return pd.read_excel(path, **kw)
    if p.endswith((".tsv", ".txt")):
        return pd.read_csv(path, sep="\t", **kw)
    return pd.read_csv(path, **kw)


def _module_genes():
    import importlib.util
    spec = importlib.util.spec_from_file_location("m16", REPO / "scripts" / "16_temporal_bridge.py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return {g for genes in m.MODULES.values() for g in genes}


def from_deseq(specs):
    frames = {}
    gene_col = None
    for spec in specs:
        if "=" not in spec:
            sys.exit(f"--deseq expects file=hour, got '{spec}'")
        path, hr = spec.rsplit("=", 1)
        hr = float(hr)
        d = pd.read_csv(path)
        # gene column may be the index (unnamed first col)
        gc = _pick(d.columns, GENE_CANDIDATES)
        if gc is None:
            gc = d.columns[0]  # assume first column is the gene id
        lc = _pick(d.columns, LFC_CANDIDATES)
        if lc is None:
            sys.exit(f"{path}: no log2FC column found (looked for {LFC_CANDIDATES})")
        s = d[[gc, lc]].copy()
        s[gc] = s[gc].astype(str).str.upper()
        s = s.dropna(subset=[lc]).drop_duplicates(gc).set_index(gc)[lc]
        frames[f"{hr:g}h"] = s
        gene_col = "gene"
    mat = pd.DataFrame(frames)
    mat.insert(0, "0h", 0.0)  # vehicle baseline = 0 log2FC
    mat = mat.reset_index().rename(columns={"index": "gene"})
    return mat


def from_counts(counts_path, design_path, treatment=None, control=None):
    counts = _read_table(counts_path, index_col=0)
    design = _read_table(design_path)
    if not {"sample", "hour"}.issubset(design.columns):
        sys.exit("--design needs columns: sample,hour (+ optional condition)")
    design["sample"] = design["sample"].astype(str)
    missing = [s for s in design["sample"] if s not in counts.columns]
    if missing:
        sys.exit(f"design samples not in counts columns: {missing[:5]}...")
    hours = sorted(design["hour"].unique())

    def mean_cols(sel):
        cols = sel["sample"].tolist()
        return counts[cols].mean(axis=1) if cols else None

    out = {}
    if treatment and control:
        # per-hour log2FC(treatment / control) — the stimulus-induced trajectory
        if "condition" not in design.columns:
            sys.exit("--treatment/--control require a 'condition' column in --design")
        for h in hours:
            t = mean_cols(design[(design.hour == h) & (design.condition == treatment)])
            c = mean_cols(design[(design.hour == h) & (design.condition == control)])
            if t is None or c is None:
                continue
            out[f"{h:g}h"] = np.log2((t + 1) / (c + 1))
    else:
        # log2FC vs earliest hour (single-arm trajectory)
        base = mean_cols(design[design.hour == hours[0]])
        out[f"{hours[0]:g}h"] = np.zeros(len(counts))
        for h in hours:
            out[f"{h:g}h"] = np.log2((mean_cols(design[design.hour == h]) + 1) / (base + 1))

    mat = pd.DataFrame(out, index=counts.index)
    mat.index = mat.index.astype(str).str.upper()
    mat = mat.reset_index()
    mat.columns = ["gene"] + list(mat.columns[1:])
    return mat


def from_long(path):
    d = pd.read_csv(path)
    need = {"gene", "hour", "value"}
    if not need.issubset(d.columns):
        sys.exit(f"--long needs columns: {need}")
    d["gene"] = d["gene"].astype(str).str.upper()
    mat = d.pivot_table(index="gene", columns="hour", values="value", aggfunc="mean")
    mat.columns = [f"{float(c):g}h" for c in mat.columns]
    return mat.reset_index()


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--deseq", nargs="+", metavar="FILE=HOUR",
                    help="per-timepoint DESeq2 result files mapped to hours")
    ap.add_argument("--counts", help="wide normalized-count matrix, CSV/TSV/XLSX (genes × samples)")
    ap.add_argument("--design", help="sample,hour(,condition) table for --counts")
    ap.add_argument("--treatment", help="condition label for stimulus (log2FC vs --control per hour)")
    ap.add_argument("--control", help="condition label for vehicle/control")
    ap.add_argument("--long", help="long table gene,hour,value")
    ap.add_argument("--genes-only", action="store_true",
                    help="keep only GR/PPARγ/VDR module genes")
    ap.add_argument("--out", default=str(RESULTS / "timecourse_matrix.csv"))
    args = ap.parse_args()

    if args.deseq:
        mat = from_deseq(args.deseq)
    elif args.counts and args.design:
        mat = from_counts(args.counts, args.design, args.treatment, args.control)
    elif args.long:
        mat = from_long(args.long)
    else:
        sys.exit("choose one input: --deseq ... | --counts c --design d | --long f")

    # order timepoint columns by hour
    gcol = mat.columns[0]
    tcols = [c for c in mat.columns if c != gcol]
    tcols = sorted(tcols, key=lambda c: float(re.search(r"[-+]?\d*\.?\d+", c).group()))
    mat = mat[[gcol] + tcols]

    if args.genes_only:
        keep = _module_genes()
        mat = mat[mat[gcol].isin(keep)]

    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    mat.to_csv(args.out, index=False)
    print(f"matrix: {mat.shape[0]} genes × {len(tcols)} timepoints {tcols}")
    print(f"written → {args.out}")
    # quick peek at module genes present
    present = mat[gcol].isin(_module_genes())
    print(f"module genes present: {present.sum()}")
    print(f"\nnext:  python scripts/16_temporal_bridge.py --matrix {args.out}")


if __name__ == "__main__":
    main()
