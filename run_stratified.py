#!/usr/bin/env python3
"""
Stratified ATAC-seq motif analysis: GRE (NR3C1) vs DR3 (VDR::RXR) in
LPS-responsive chromatin.

Peaks from an unstimulated and a stimulated ATAC-seq sample are stratified
along two axes:

  1. Response class
       gained : stim peak with no overlapping unstim peak
       lost   : unstim peak with no overlapping stim peak
       shared : stim peak that overlaps an unstim peak (reference set)
  2. Distance from summit to the nearest TSS
       promoter : <= --promoter-bp           (default 1 kb)
       proximal : <= --proximal-bp           (default 10 kb)
       distal   :  > --proximal-bp

For every peak a summit-centred window (+/- --half-width bp) is scanned on
both strands with each motif. A peak is a motif hit when its best site
reaches p <= --motif-p (FIMO-style p-value from the exact score
distribution under the background base composition).

Enrichment of each motif in gained / lost peaks is tested against shared
peaks of the same distance stratum:
  - Fisher's exact test (unadjusted)
  - Cochran-Mantel-Haenszel test stratified by GC-content bins
    (GC-adjusted common odds ratio), since ATAC peak sets differ in GC and
    both motifs are GC-sensitive.
BH-FDR is applied across all tests in a run.

Usage:
  python3 run_stratified.py \\
    --unstim THP1_ATAC_untreated_peaks.narrowPeak \\
    --stim   THP1_ATAC_LPS_peaks.narrowPeak \\
    --genome hg38.fa \\
    --gre MA0113.3.meme --dr3 MA0074.1.meme \\
    --tss gencode_tss.bed \\
    --out results_stratified

Outputs (in --out):
  peaks_annotated.tsv.gz   one row per peak with class, stratum, GC, motif scores
  stratum_counts.tsv       hit counts / fractions per class x stratum x motif
  enrichment.tsv           Fisher + GC-stratified CMH tests vs shared peaks
  cooccurrence.tsv         GRE/DR3 joint hit patterns per class x stratum
  run_params.json          arguments, thresholds, input sizes

Requirements: numpy, pandas, scipy
"""

import argparse
import gzip
import json
import math
import mmap
import os
import sys

import numpy as np
import pandas as pd
from scipy import stats

NARROWPEAK_COLS = ["chrom", "start", "end", "name", "score", "strand",
                   "signal", "pval", "qval", "summit"]
DIST_STRATA = ["promoter", "proximal", "distal"]
CLASSES = ["gained", "lost", "shared"]


def log(msg):
    print(f"[run_stratified] {msg}", file=sys.stderr, flush=True)


# ─── Input parsing ──────────────────────────────────────────────────────────

def open_text(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def read_narrowpeak(path):
    rows = []
    with open_text(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 3:
                continue
            f = f + ["."] * (10 - len(f)) if len(f) < 10 else f[:10]
            rows.append(f)
    df = pd.DataFrame(rows, columns=NARROWPEAK_COLS)
    df["start"] = df["start"].astype(np.int64)
    df["end"] = df["end"].astype(np.int64)
    for c in ("signal", "pval", "qval"):
        df[c] = pd.to_numeric(df[c], errors="coerce")
    summit = pd.to_numeric(df["summit"], errors="coerce")
    # narrowPeak summit is an offset from start; -1 / missing -> peak midpoint
    mid = (df["end"] - df["start"]) // 2
    summit = summit.where(summit >= 0, mid).fillna(mid).astype(np.int64)
    df["summit_pos"] = df["start"] + summit
    # MACS2 emits one line per summit for multi-summit peaks; keep the
    # strongest summit per interval so each region is counted once.
    n0 = len(df)
    df = (df.sort_values("signal", ascending=False, na_position="last")
            .drop_duplicates(["chrom", "start", "end"])
            .sort_values(["chrom", "start"]).reset_index(drop=True))
    if len(df) < n0:
        log(f"  {path}: collapsed {n0 - len(df)} extra summits of multi-summit peaks")
    return df


def read_tss(path):
    rows = []
    with open_text(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            f = line.rstrip("\n").split("\t")
            chrom, s, e = f[0], int(f[1]), int(f[2])
            name = f[3] if len(f) > 3 else "."
            strand = f[5] if len(f) > 5 else "."
            # 1-bp TSS records are used as is; for wider records the TSS is
            # the 5' end according to strand.
            if e - s <= 1:
                pos = s
            else:
                pos = e - 1 if strand == "-" else s
            rows.append((chrom, pos, name))
    df = pd.DataFrame(rows, columns=["chrom", "pos", "gene"])
    return df.drop_duplicates(["chrom", "pos"]).sort_values(["chrom", "pos"])


def read_meme(path):
    """Return (name, prob matrix [w x 4, ACGT], background [4]) of the first motif."""
    with open(path) as fh:
        lines = [ln.rstrip("\n") for ln in fh]
    bg = np.full(4, 0.25)
    name, mat, i = None, None, 0
    while i < len(lines):
        ln = lines[i].strip()
        if ln.startswith("Background letter frequencies"):
            toks = lines[i + 1].split()
            vals = dict(zip(toks[0::2], map(float, toks[1::2])))
            if all(k in vals for k in "ACGT"):
                bg = np.array([vals[k] for k in "ACGT"])
                bg = bg / bg.sum()
        elif ln.startswith("MOTIF") and name is None:
            toks = ln.split()
            name = "_".join(toks[1:3]) if len(toks) > 2 else toks[1]
        elif ln.startswith("letter-probability matrix") and mat is None:
            w = None
            for t, v in zip(ln.split()[2::2], ln.split()[3::2]):
                if t.rstrip("=") == "w":
                    w = int(v)
            rows, j = [], i + 1
            while j < len(lines) and (w is None or len(rows) < w):
                toks = lines[j].split()
                if len(toks) == 4:
                    try:
                        rows.append([float(x) for x in toks])
                    except ValueError:
                        break
                elif rows:
                    break
                j += 1
            mat = np.array(rows)
            i = j
            continue
        i += 1
    if mat is None or mat.ndim != 2 or mat.shape[1] != 4 or len(mat) == 0:
        raise ValueError(f"{path}: no letter-probability matrix found")
    mat = mat / mat.sum(axis=1, keepdims=True)
    return name or os.path.basename(path), mat, bg


# ─── FASTA access (faidx-compatible index, mmap) ────────────────────────────

class Fasta:
    def __init__(self, path, index_dir):
        if path.endswith(".gz"):
            raise SystemExit(f"{path}: gzipped FASTA is not supported; "
                             "decompress it first (gunzip -k).")
        self.path = path
        fai = path + ".fai"
        if not os.path.exists(fai):
            fai = os.path.join(index_dir, os.path.basename(path) + ".fai")
            if not os.path.exists(fai):
                log(f"building FASTA index {fai} (one-time scan)")
                self._build_index(path, fai)
        self.index = {}
        with open(fai) as fh:
            for line in fh:
                n, length, off, lb, lw = line.split("\t")[:5]
                self.index[n] = (int(length), int(off), int(lb), int(lw))
        self._fh = open(path, "rb")
        self._mm = mmap.mmap(self._fh.fileno(), 0, access=mmap.ACCESS_READ)

    @staticmethod
    def _build_index(path, fai):
        entries = []
        name = None
        with open(path, "rb") as fh:
            offset = 0
            for raw in fh:
                if raw.startswith(b">"):
                    if name is not None:
                        entries.append((name, length, seq_off, lb, lw))
                    name = raw[1:].split()[0].decode()
                    length, seq_off, lb, lw, short_seen = 0, offset + len(raw), None, None, False
                else:
                    bases = len(raw.rstrip(b"\r\n"))
                    if lb is None:
                        lb, lw = bases, len(raw)
                    elif short_seen and bases:
                        raise SystemExit(f"{path}: non-uniform line length in {name}; "
                                         "reformat the FASTA (e.g. samtools faidx).")
                    elif bases != lb:
                        short_seen = True
                    length += bases
                offset += len(raw)
            if name is not None:
                entries.append((name, length, seq_off, lb or 0, lw or 0))
        with open(fai, "w") as out:
            for e in entries:
                out.write("\t".join(map(str, e)) + "\n")

    def fetch(self, chrom, start, end):
        """0-based half-open; positions off the chromosome are returned as N."""
        length, off, lb, lw = self.index[chrom]
        s, e = max(0, start), min(length, end)
        if e <= s:
            return b"N" * (end - start)
        b0 = off + (s // lb) * lw + s % lb
        b1 = off + ((e - 1) // lb) * lw + (e - 1) % lb + 1
        seq = self._mm[b0:b1].replace(b"\n", b"").replace(b"\r", b"")
        return b"N" * (s - start) + seq.upper() + b"N" * (end - e)


# ─── Peak stratification ────────────────────────────────────────────────────

def overlaps_any(query, ref, min_frac):
    """Boolean per query peak: overlaps >= min_frac of its own length with any ref peak."""
    hit = np.zeros(len(query), dtype=bool)
    for chrom, q in query.groupby("chrom"):
        r = ref[ref["chrom"] == chrom]
        if r.empty:
            continue
        rs = r["start"].to_numpy()
        re_ = r["end"].to_numpy()
        order = np.argsort(rs)
        rs, re_ = rs[order], re_[order]
        max_len = int((re_ - rs).max())
        qs, qe = q["start"].to_numpy(), q["end"].to_numpy()
        lo = np.searchsorted(rs, qs - max_len, side="left")
        hi = np.searchsorted(rs, qe, side="left")
        res = np.zeros(len(q), dtype=bool)
        for k in range(len(q)):
            if hi[k] <= lo[k]:
                continue
            ov = np.minimum(qe[k], re_[lo[k]:hi[k]]) - np.maximum(qs[k], rs[lo[k]:hi[k]])
            if ov.size and ov.max() >= max(1, min_frac * (qe[k] - qs[k])):
                res[k] = True
        hit[q.index.to_numpy()] = res
    return hit


def nearest_tss(peaks, tss):
    dist = np.full(len(peaks), np.nan)
    gene = np.full(len(peaks), ".", dtype=object)
    for chrom, p in peaks.groupby("chrom"):
        t = tss[tss["chrom"] == chrom]
        if t.empty:
            continue
        tp = t["pos"].to_numpy()
        tg = t["gene"].to_numpy()
        x = p["summit_pos"].to_numpy()
        i = np.clip(np.searchsorted(tp, x), 1, len(tp)) - 1
        j = np.clip(i + 1, 0, len(tp) - 1)
        di, dj = np.abs(x - tp[i]), np.abs(x - tp[j])
        best = np.where(dj < di, j, i)
        dist[p.index.to_numpy()] = np.minimum(di, dj)
        gene[p.index.to_numpy()] = tg[best]
    return dist, gene


def distance_stratum(dist, promoter_bp, proximal_bp):
    out = np.full(len(dist), "no_tss", dtype=object)
    ok = ~np.isnan(dist)
    out[ok & (dist <= promoter_bp)] = "promoter"
    out[ok & (dist > promoter_bp) & (dist <= proximal_bp)] = "proximal"
    out[ok & (dist > proximal_bp)] = "distal"
    return out


# ─── Motif scanning ─────────────────────────────────────────────────────────

ENCODE = np.full(256, 4, dtype=np.uint8)
for _i, _c in enumerate(b"ACGT"):
    ENCODE[_c] = _i


class Motif:
    SCALE = 100  # integer score resolution for the exact p-value DP

    def __init__(self, name, prob, bg, pseudocount=0.01):
        self.name = name
        self.width = len(prob)
        p = (prob + pseudocount * bg) / (1 + pseudocount)
        lo = np.log2(p / bg)
        self.int_pssm = np.round(lo * self.SCALE).astype(np.int64)        # w x 4
        rc = self.int_pssm[::-1, ::-1]                                    # reverse complement
        # column 4 = N: scored at the column minimum so N never creates a site
        self.fwd = np.hstack([self.int_pssm, self.int_pssm.min(1, keepdims=True)])
        self.rev = np.hstack([rc, rc.min(1, keepdims=True)])
        self._build_pvalues(bg)

    def _build_pvalues(self, bg):
        mins = self.int_pssm.min(axis=1)
        offset = int(mins.sum())
        span = int((self.int_pssm.max(axis=1) - mins).sum())
        dist = np.zeros(span + 1)
        dist[0] = 1.0
        cur = 0
        for row, m in zip(self.int_pssm, mins):
            nxt = np.zeros(span + 1)
            width = int(row.max() - m)
            for b in range(4):
                sh = int(row[b] - m)
                nxt[sh:sh + cur + 1] += dist[:cur + 1] * bg[b]
            dist, cur = nxt, cur + width
        sf = np.cumsum(dist[::-1])[::-1]                                   # P(S >= s)
        self._offset = offset
        self._sf = np.clip(sf, 0, 1)

    def pvalue(self, score):
        idx = np.clip(np.asarray(score) - self._offset, 0, len(self._sf) - 1)
        return self._sf[idx]

    def threshold(self, p):
        ok = np.nonzero(self._sf <= p)[0]
        return int(ok[0] + self._offset) if ok.size else None

    def best_scores(self, codes):
        """codes: n x L uint8 (0-4). Returns best score over both strands per row."""
        n, L = codes.shape
        npos = L - self.width + 1
        best = np.full(n, np.iinfo(np.int64).min, dtype=np.int64)
        for mat in (self.fwd, self.rev):
            acc = np.zeros((n, npos), dtype=np.int64)
            for j in range(self.width):
                acc += mat[j][codes[:, j:j + npos]]
            best = np.maximum(best, acc.max(axis=1))
        return best


def extract_windows(fasta, peaks, half_width):
    L = 2 * half_width
    codes = np.full((len(peaks), L), 4, dtype=np.uint8)
    gc = np.full(len(peaks), np.nan)
    in_genome = np.zeros(len(peaks), dtype=bool)
    for k, (chrom, summit) in enumerate(zip(peaks["chrom"].to_numpy(),
                                            peaks["summit_pos"].to_numpy())):
        if chrom not in fasta.index:
            continue
        seq = fasta.fetch(chrom, int(summit) - half_width, int(summit) + half_width)
        arr = ENCODE[np.frombuffer(seq, dtype=np.uint8)]
        codes[k] = arr
        acgt = (arr < 4).sum()
        if acgt:
            gc[k] = ((arr == 1) | (arr == 2)).sum() / acgt
            in_genome[k] = True
    return codes, gc, in_genome


# ─── Statistics ─────────────────────────────────────────────────────────────

def cmh_test(tables):
    """Cochran-Mantel-Haenszel over 2x2 tables [[a,b],[c,d]]; returns (OR_MH, p)."""
    num = den = sum_a = sum_e = sum_v = 0.0
    for a, b, c, d in tables:
        n = a + b + c + d
        if n < 2 or (a + b) == 0 or (c + d) == 0 or (a + c) == 0 or (b + d) == 0:
            continue
        num += a * d / n
        den += b * c / n
        sum_a += a
        sum_e += (a + b) * (a + c) / n
        sum_v += (a + b) * (c + d) * (a + c) * (b + d) / (n * n * (n - 1))
    if sum_v == 0:
        return np.nan, np.nan
    chi2 = (abs(sum_a - sum_e) - 0.5) ** 2 / sum_v
    or_mh = num / den if den > 0 else np.inf
    return or_mh, float(stats.chi2.sf(max(chi2, 0.0), 1))


def bh_fdr(p):
    p = np.asarray(p, dtype=float)
    out = np.full(len(p), np.nan)
    ok = ~np.isnan(p)
    if ok.sum():
        pv = p[ok]
        order = np.argsort(pv)
        ranked = pv[order] * ok.sum() / (np.arange(ok.sum()) + 1)
        ranked = np.minimum.accumulate(ranked[::-1])[::-1]
        res = np.empty_like(pv)
        res[order] = np.minimum(ranked, 1)
        out[ok] = res
    return out


def gc_bins(gc_ref, n_bins):
    qs = np.nanquantile(gc_ref, np.linspace(0, 1, n_bins + 1))
    qs = np.unique(qs)
    qs[0], qs[-1] = -np.inf, np.inf
    return qs


# ─── Main ───────────────────────────────────────────────────────────────────

def parse_args(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--unstim", required=True, help="unstimulated narrowPeak (.gz ok)")
    ap.add_argument("--stim", required=True, help="stimulated narrowPeak (.gz ok)")
    ap.add_argument("--genome", required=True, help="genome FASTA (uncompressed; .fai used if present)")
    ap.add_argument("--gre", required=True, help="GRE motif, MEME format (e.g. JASPAR MA0113.3 NR3C1)")
    ap.add_argument("--dr3", required=True, help="DR3 motif, MEME format (e.g. JASPAR MA0074.1 RXRA::VDR)")
    ap.add_argument("--tss", required=True, help="TSS BED (chrom, start, end[, gene, score, strand])")
    ap.add_argument("--out", required=True, help="output directory")
    ap.add_argument("--half-width", type=int, default=100,
                    help="scan window = summit +/- this many bp (default 100)")
    ap.add_argument("--min-overlap-frac", type=float, default=0.0,
                    help="fraction of a peak that must overlap the other condition to count as "
                         "shared; 0 = any 1-bp overlap (default 0)")
    ap.add_argument("--promoter-bp", type=int, default=1000, help="promoter cut-off (default 1000)")
    ap.add_argument("--proximal-bp", type=int, default=10000, help="proximal cut-off (default 10000)")
    ap.add_argument("--motif-p", type=float, default=1e-4,
                    help="per-site p-value threshold for a motif hit (default 1e-4)")
    ap.add_argument("--gc-bins", type=int, default=5,
                    help="GC quantile bins for the CMH test (default 5)")
    ap.add_argument("--chunk", type=int, default=20000, help="peaks per scanning chunk")
    return ap.parse_args(argv)


def main(argv=None):
    a = parse_args(argv)
    if a.promoter_bp >= a.proximal_bp:
        raise SystemExit("--promoter-bp must be smaller than --proximal-bp")
    for p in (a.unstim, a.stim, a.genome, a.gre, a.dr3, a.tss):
        if not os.path.exists(p):
            raise SystemExit(f"input not found: {p}")
    os.makedirs(a.out, exist_ok=True)

    log("reading peaks")
    unstim, stim = read_narrowpeak(a.unstim), read_narrowpeak(a.stim)
    log(f"  unstim {len(unstim):,} peaks | stim {len(stim):,} peaks")

    stim_shared = overlaps_any(stim, unstim, a.min_overlap_frac)
    unstim_shared = overlaps_any(unstim, stim, a.min_overlap_frac)
    stim["peak_class"] = np.where(stim_shared, "shared", "gained")
    lost = unstim[~unstim_shared].copy()
    lost["peak_class"] = "lost"
    stim["source"], lost["source"] = "stim", "unstim"
    peaks = pd.concat([stim, lost], ignore_index=True)
    log(f"  gained {int((~stim_shared).sum()):,} | shared {int(stim_shared.sum()):,} | "
        f"lost {len(lost):,}")

    log("annotating nearest TSS")
    tss = read_tss(a.tss)
    dist, gene = nearest_tss(peaks, tss)
    peaks["tss_dist"] = dist
    peaks["nearest_gene"] = gene
    peaks["dist_stratum"] = distance_stratum(dist, a.promoter_bp, a.proximal_bp)

    log("loading motifs")
    motifs = {}
    for key, path in (("GRE", a.gre), ("DR3", a.dr3)):
        name, prob, bg = read_meme(path)
        motifs[key] = Motif(name, prob, bg)
        thr = motifs[key].threshold(a.motif_p)
        log(f"  {key}: {name} (w={motifs[key].width}), score threshold {thr} at p<={a.motif_p:g}")
        if 2 * a.half_width < motifs[key].width:
            raise SystemExit("--half-width too small for motif width")

    log("extracting sequences and scanning")
    fasta = Fasta(a.genome, a.out)
    n = len(peaks)
    gc = np.full(n, np.nan)
    in_genome = np.zeros(n, dtype=bool)
    best = {k: np.zeros(n, dtype=np.int64) for k in motifs}
    for s in range(0, n, a.chunk):
        sub = peaks.iloc[s:s + a.chunk]
        codes, gc_s, ok_s = extract_windows(fasta, sub, a.half_width)
        gc[s:s + len(sub)] = gc_s
        in_genome[s:s + len(sub)] = ok_s
        for k, m in motifs.items():
            best[k][s:s + len(sub)] = m.best_scores(codes)
        log(f"  {min(s + a.chunk, n):,}/{n:,}")
    missing = int((~in_genome).sum())
    if missing:
        chroms = sorted(set(peaks.loc[~in_genome, "chrom"]))[:5]
        log(f"  WARNING: {missing:,} peaks not in genome FASTA (e.g. {chroms}); excluded")
    if missing == n:
        raise SystemExit("no peak chromosome matched the FASTA (chr1 vs 1 naming?)")

    peaks["gc"] = gc
    for k, m in motifs.items():
        peaks[f"{k}_best_score"] = best[k] / Motif.SCALE
        peaks[f"{k}_best_p"] = m.pvalue(best[k])
        peaks[f"{k}_hit"] = peaks[f"{k}_best_p"] <= a.motif_p
    peaks = peaks[in_genome].reset_index(drop=True)

    cols = ["chrom", "start", "end", "summit_pos", "name", "signal", "qval", "source",
            "peak_class", "tss_dist", "nearest_gene", "dist_stratum", "gc"]
    for k in motifs:
        cols += [f"{k}_best_score", f"{k}_best_p", f"{k}_hit"]
    peaks[cols].to_csv(os.path.join(a.out, "peaks_annotated.tsv.gz"),
                       sep="\t", index=False, compression="gzip")

    log("computing stratified statistics")
    strata = DIST_STRATA + ["all"]
    counts, enrich, cooc = [], [], []
    for st in strata:
        sel = peaks if st == "all" else peaks[peaks["dist_stratum"] == st]
        for cls in CLASSES:
            g = sel[sel["peak_class"] == cls]
            row = {"dist_stratum": st, "peak_class": cls, "n_peaks": len(g),
                   "median_gc": float(g["gc"].median()) if len(g) else np.nan}
            for k in motifs:
                h = int(g[f"{k}_hit"].sum())
                row[f"{k}_hits"] = h
                row[f"{k}_frac"] = h / len(g) if len(g) else np.nan
            counts.append(row)
            if len(g):
                gh, dh = g["GRE_hit"], g["DR3_hit"]
                cooc.append({"dist_stratum": st, "peak_class": cls, "n_peaks": len(g),
                             "GRE_and_DR3": int((gh & dh).sum()),
                             "GRE_only": int((gh & ~dh).sum()),
                             "DR3_only": int((~gh & dh).sum()),
                             "neither": int((~gh & ~dh).sum())})

        ref = sel[sel["peak_class"] == "shared"]
        if len(ref) == 0:
            continue
        edges = gc_bins(ref["gc"].to_numpy(), a.gc_bins)
        for cls in ("gained", "lost"):
            test = sel[sel["peak_class"] == cls]
            if len(test) == 0:
                continue
            for k in motifs:
                a1 = int(test[f"{k}_hit"].sum()); b1 = len(test) - a1
                c1 = int(ref[f"{k}_hit"].sum()); d1 = len(ref) - c1
                or_f, p_f = stats.fisher_exact([[a1, b1], [c1, d1]])
                tb = []
                tbin = np.digitize(test["gc"], edges[1:-1])
                rbin = np.digitize(ref["gc"], edges[1:-1])
                for b in range(len(edges) - 1):
                    th = test[f"{k}_hit"].to_numpy()[tbin == b]
                    rh = ref[f"{k}_hit"].to_numpy()[rbin == b]
                    tb.append((int(th.sum()), int(len(th) - th.sum()),
                               int(rh.sum()), int(len(rh) - rh.sum())))
                or_mh, p_mh = cmh_test(tb)
                enrich.append({
                    "dist_stratum": st, "motif": k, "motif_id": motifs[k].name,
                    "test_class": cls, "ref_class": "shared",
                    "test_hits": a1, "test_n": len(test),
                    "ref_hits": c1, "ref_n": len(ref),
                    "test_frac": a1 / len(test), "ref_frac": c1 / len(ref),
                    "fisher_OR": or_f, "fisher_p": p_f,
                    "cmh_gc_OR": or_mh, "cmh_gc_p": p_mh,
                    "test_median_gc": float(test["gc"].median()),
                    "ref_median_gc": float(ref["gc"].median()),
                })

    counts = pd.DataFrame(counts)
    enrich = pd.DataFrame(enrich)
    if len(enrich):
        # FDR within the three distance strata; the pooled "all" rows are a
        # summary of the same peaks and are corrected separately.
        for mask in (enrich["dist_stratum"] != "all", enrich["dist_stratum"] == "all"):
            enrich.loc[mask, "fisher_fdr"] = bh_fdr(enrich.loc[mask, "fisher_p"])
            enrich.loc[mask, "cmh_gc_fdr"] = bh_fdr(enrich.loc[mask, "cmh_gc_p"])
    counts.to_csv(os.path.join(a.out, "stratum_counts.tsv"), sep="\t", index=False)
    enrich.to_csv(os.path.join(a.out, "enrichment.tsv"), sep="\t", index=False)
    pd.DataFrame(cooc).to_csv(os.path.join(a.out, "cooccurrence.tsv"), sep="\t", index=False)

    params = vars(a).copy()
    params.update({
        "n_unstim_peaks": len(unstim), "n_stim_peaks": len(stim),
        "n_peaks_not_in_genome": missing,
        "n_peaks_no_tss_chrom": int((peaks["dist_stratum"] == "no_tss").sum()),
        "class_counts": peaks["peak_class"].value_counts().to_dict(),
        "motifs": {k: {"id": m.name, "width": m.width,
                       "score_threshold": (m.threshold(a.motif_p) or 0) / Motif.SCALE}
                   for k, m in motifs.items()},
    })
    with open(os.path.join(a.out, "run_params.json"), "w") as fh:
        json.dump(params, fh, indent=2, default=str)

    show = enrich[enrich["dist_stratum"] != "all"] if len(enrich) else enrich
    if len(show):
        with pd.option_context("display.width", 160, "display.max_columns", 20):
            print(show[["dist_stratum", "motif", "test_class", "test_frac", "ref_frac",
                        "fisher_OR", "fisher_fdr", "cmh_gc_OR", "cmh_gc_fdr"]]
                  .to_string(index=False, float_format=lambda x: f"{x:.3g}"))
    log(f"done -> {a.out}")


if __name__ == "__main__":
    main()
