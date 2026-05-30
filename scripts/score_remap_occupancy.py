#!/usr/bin/env python3
"""
Canonical ReMAP2022 VDR / GR occupancy scoring — single, reproducible pipeline
==============================================================================

WHY THIS SCRIPT EXISTS
----------------------
The committed `results/remap_scores_expanded.csv` (Supplementary Table 1, 381
genes) was produced by a patchwork of scoring definitions:
  * `expand_steroid_targets.py` uses a COUNT score  =  10*n_celltypes + n_exp
    (no peak-height information), within a +-5 kb window.
  * ~270/381 rows carry decimal scores inherited from an earlier signal-based
    method that is not committed.
  * The Methods text states a THIRD definition (+-10 kb, "sum of peak signal
    intensities", pybedtools/BEDTools).
  * TNFSF15 (tulisokibart target, VDR=12.4 / GR=286.0 in the OSF
    pre-registration) is absent from the table and matches none of the above.

Because these definitions are mutually inconsistent, TNFSF15 cannot be "merged
in" meaningfully. This script defines ONE explicit scoring rule and re-scores
EVERY gene (the existing universe + TNFSF15) in a SINGLE run, emitting both the
signal-sum score and a depth-normalised (per-experiment) score so the
study-depth confound is transparent.

DEFINITION (defaults; all parameterisable)
------------------------------------------
For a gene with TSS at (chrom, pos), consider every ReMAP2022 peak overlapping
[pos - W, pos + W] (default W = 10 kb, strand-aware TSS):
  signal_sum        = sum of peak signalValue (BED col 7)   <- Methods-stated
  n_experiments     = # distinct ReMAP datasets (GSE id, name field part 0)
  n_celltypes       = # distinct cell types (name field part 2, suffix-stripped)
  signal_per_exp    = signal_sum / n_experiments            <- depth-corrected
ReMAP name field format: "<GSEid>.<TF>.<celltype>" (dot-delimited).

USAGE
-----
  python3 scripts/score_remap_occupancy.py \
      --vdr-bed /path/remap2022_VDR_all_macs2_hg38.bed.gz \
      --gr-bed  /path/remap2022_NR3C1_all_macs2_hg38.bed.gz

Runs offline if every gene's TSS is in the cache (results/gene_tss_hg38.csv);
otherwise missing TSS are fetched from Ensembl REST (rest.ensembl.org) and
written back to the cache. The two bed.gz files are NOT in this repo; obtain
them from https://remap.univ-amu.fr/ (hg38, "all" MACS2 peaks).

OUTPUT
------
  results/remap_scores_canonical.csv   (does NOT overwrite the legacy table)
"""

import argparse
import gzip
import json
import os
import re
import sys
import time
import urllib.request

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.normpath(os.path.join(HERE, "..", "results"))

# TNFSF15 TSS from the OSF pre-registration (osf.io/tnp63), Ensembl GRCh38:
# chr9:114,806,039, strand -1.  Seeded so the prospective-prediction target is
# always scored by the same rule as everything else.
TNFSF15_SEED = dict(gene="TNFSF15", chrom="chr9", tss=114806039, strand=-1)


# ──────────────────────────────────────────────────────────────────────────────
# TSS resolution
# ──────────────────────────────────────────────────────────────────────────────
def load_tss_cache(path):
    if os.path.exists(path):
        df = pd.read_csv(path)
        return {r.gene.upper(): (str(r.chrom), int(r.tss), int(r.strand))
                for r in df.itertuples()}
    return {}


def save_tss_cache(path, tss):
    rows = [dict(gene=g, chrom=c, tss=p, strand=s) for g, (c, p, s) in sorted(tss.items())]
    pd.DataFrame(rows).to_csv(path, index=False)


def fetch_tss_ensembl(genes, tss, cache_path, batch=45):
    """Fill missing TSS from Ensembl REST; update cache incrementally."""
    missing = [g for g in genes if g.upper() not in tss]
    if not missing:
        return tss
    url = "https://rest.ensembl.org/lookup/symbol/homo_sapiens"
    print(f"  resolving {len(missing)} TSS via Ensembl REST ...")
    for i in range(0, len(missing), batch):
        chunk = missing[i:i + batch]
        data = json.dumps({"symbols": chunk}).encode()
        req = urllib.request.Request(
            url, data=data,
            headers={"Content-Type": "application/json", "Accept": "application/json"})
        try:
            with urllib.request.urlopen(req, timeout=60) as r:
                d = json.loads(r.read())
        except Exception as e:
            print(f"    Ensembl error ({e}); {len(chunk)} genes left unresolved")
            continue
        for gene, info in (d or {}).items():
            if not info:
                continue
            chrom = str(info.get("seq_region_name", ""))
            if not re.match(r"^(\d+|X|Y)$", chrom):
                continue
            strand = int(info.get("strand", 1))
            pos = int(info["start"]) if strand == 1 else int(info["end"])
            tss[gene.upper()] = ("chr" + chrom, pos, strand)
        save_tss_cache(cache_path, tss)   # persist as we go
        time.sleep(0.4)
    return tss


# ──────────────────────────────────────────────────────────────────────────────
# ReMAP bed indexing + scoring
# ──────────────────────────────────────────────────────────────────────────────
def build_index(bed_gz):
    """chrom -> list of (start, end, celltype, gse_id, signal)."""
    idx = {}
    n_exp = set()
    with gzip.open(bed_gz, "rt") as f:
        for line in f:
            if line.startswith(("#", "track", "browser")):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 4:
                continue
            chrom, start, end, name = p[0], int(p[1]), int(p[2]), p[3]
            signal = float(p[6]) if len(p) > 6 and p[6] not in (".", "") else 0.0
            parts = name.split(".")
            gse = parts[0] if parts else "unknown"
            celltype = parts[2] if len(parts) > 2 else "unknown"
            base_cell = re.sub(r"_[A-Z0-9_]+$", "", celltype, flags=re.IGNORECASE)
            idx.setdefault(chrom, []).append((start, end, base_cell, gse, signal))
            n_exp.add(gse)
    # sort each chrom by start for a quick scan
    for c in idx:
        idx[c].sort(key=lambda x: x[0])
    return idx, len(n_exp)


def score_gene(chrom, tss, idx, window):
    lo, hi = tss - window, tss + window
    sig_sum = 0.0
    cells, exps = set(), set()
    for s, e, cell, gse, signal in idx.get(chrom, []):
        if s > hi:
            break
        if e < lo:
            continue
        sig_sum += signal
        cells.add(cell)
        exps.add(gse)
    n_exp = len(exps)
    per_exp = sig_sum / n_exp if n_exp else 0.0
    return round(sig_sum, 3), n_exp, len(cells), round(per_exp, 3)


# ──────────────────────────────────────────────────────────────────────────────
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--vdr-bed", required=True, help="remap2022_VDR_all_macs2_hg38.bed.gz")
    ap.add_argument("--gr-bed", required=True, help="remap2022_NR3C1_all_macs2_hg38.bed.gz")
    ap.add_argument("--window", type=int, default=10000, help="half-window around TSS (bp); default 10000")
    ap.add_argument("--gene-universe", default=f"{RESULTS}/remap_scores_expanded.csv",
                    help="CSV with a 'gene' column defining which genes to score")
    ap.add_argument("--tss-cache", default=f"{RESULTS}/gene_tss_hg38.csv")
    ap.add_argument("--out", default=f"{RESULTS}/remap_scores_canonical.csv")
    ap.add_argument("--no-network", action="store_true", help="never call Ensembl; skip uncached TSS")
    args = ap.parse_args()

    # Gene universe = existing table genes + TNFSF15 (deduplicated, upper-cased)
    universe = pd.read_csv(args.gene_universe)["gene"].astype(str).str.upper().tolist()
    genes = sorted(set(universe) | {"TNFSF15"})
    print(f"Scoring {len(genes)} genes (incl. TNFSF15), window=+-{args.window} bp")

    # TSS: cache (seeded with TNFSF15) -> Ensembl
    tss = load_tss_cache(args.tss_cache)
    tss.setdefault("TNFSF15", (TNFSF15_SEED["chrom"], TNFSF15_SEED["tss"], TNFSF15_SEED["strand"]))
    if not args.no_network:
        tss = fetch_tss_ensembl(genes, tss, args.tss_cache)
    save_tss_cache(args.tss_cache, tss)
    resolved = [g for g in genes if g in tss]
    print(f"TSS resolved: {len(resolved)}/{len(genes)}"
          + ("" if len(resolved) == len(genes)
             else f"  (unresolved skipped: {sorted(set(genes)-set(resolved))[:8]}...)"))

    print("Indexing ReMAP bed files ...")
    vdr_idx, n_vdr_total = build_index(args.vdr_bed)
    gr_idx, n_gr_total = build_index(args.gr_bed)
    print(f"  VDR datasets total={n_vdr_total}   GR datasets total={n_gr_total}")

    rows = []
    for g in resolved:
        chrom, pos, strand = tss[g]
        vs, ve, vc, vpe = score_gene(chrom, pos, vdr_idx, args.window)
        gs, ge, gc, gpe = score_gene(chrom, pos, gr_idx, args.window)
        rows.append(dict(
            gene=g, chrom=chrom, tss=pos, strand=strand,
            VDR_signal_sum=vs, VDR_n_exp=ve, VDR_n_celltype=vc, VDR_signal_per_exp=vpe,
            GR_signal_sum=gs, GR_n_exp=ge, GR_n_celltype=gc, GR_signal_per_exp=gpe,
            VDR_dominant_signal=vs > gs,
            VDR_dominant_perexp=vpe > gpe,
        ))
    out = pd.DataFrame(rows).sort_values("gene").reset_index(drop=True)
    out.to_csv(args.out, index=False)

    print(f"\nWrote {len(out)} genes -> {args.out}")
    t = out[out.gene == "TNFSF15"]
    if len(t):
        r = t.iloc[0]
        print("\nTNFSF15 (re-scored by the canonical rule):")
        print(f"  VDR: signal_sum={r.VDR_signal_sum} (n_exp={r.VDR_n_exp}) "
              f"per_exp={r.VDR_signal_per_exp}")
        print(f"  GR : signal_sum={r.GR_signal_sum} (n_exp={r.GR_n_exp}) "
              f"per_exp={r.GR_signal_per_exp}")
        rank_signal = int((out.GR_signal_sum > r.GR_signal_sum).sum()) + 1
        rank_perexp = int((out.GR_signal_per_exp > r.GR_signal_per_exp).sum()) + 1
        print(f"  GR rank by signal_sum: #{rank_signal}/{len(out)}   "
              f"by per_exp: #{rank_perexp}/{len(out)}")
        print("  -> use these ranks to state the 'highest GR' claim accurately.")


if __name__ == "__main__":
    main()
