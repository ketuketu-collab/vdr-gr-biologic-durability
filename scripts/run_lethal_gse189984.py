#!/usr/bin/env python3
"""
GSE189984（ヒト PBMC, 10 nM 1,25(OH)2D3 vs 0.1% EtOH, 4/8/24/48 h, n=3）で
「VitD で下がった遺伝子は KO 致死遺伝子に濃縮しているか」を検定する。

手順
----
1. GEO 補足ファイル GSE189984_Normalized_counts.xlsx（DESeq2 正規化カウント）
   から生カウントを復元する。正規化値 = 生カウント / size factor なので、
   各サンプルの最小非ゼロ値を 1 カウントとみなして割り戻す。整数に戻ることを
   検証する（誤差が 1e-6 を超えたら停止）。
2. pyDESeq2 で時間点ごとに 125D vs EtOH（3 vs 3, 対応なし）。
   主結果は apeglm 縮小 LFC、感度解析として MLE LFC。
3. ヒト→マウスのオルソログは MGI HOM_MouseHumanSequence.rpt を
   Entrez Gene ID で結合（シンボルの改名に影響されない）。
4. IMPC viability.csv の遺伝子シンボルは MGI ID で MGI の現行シンボルに揃える。
5. lethal_enrichment.py を発現量（baseMean）で層別して実行する。

入力ファイルの取得
------------------
  https://ftp.ncbi.nlm.nih.gov/geo/series/GSE189nnn/GSE189984/suppl/GSE189984_Normalized_counts.xlsx
  https://ftp.ebi.ac.uk/pub/databases/impc/all-data-releases/latest/results/viability.csv.gz
  https://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt

使い方
------
python3 scripts/run_lethal_gse189984.py \
    --counts-xlsx GSE189984_Normalized_counts.xlsx \
    --impc viability.csv.gz --mgi HOM_MouseHumanSequence.rpt \
    --out results/lethal_enrichment
"""

import argparse
import os
import subprocess
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
TIMES = ["04h", "08h", "24h", "48h"]


def recover_counts(xlsx):
    d = pd.read_excel(xlsx)
    samples = [c for c in d.columns if c.startswith(("EtOH_", "125D_"))]
    d = d.dropna(subset=["entrezgene_id"])
    d["entrez"] = d["entrezgene_id"].astype(int).astype(str)
    # 同じ Entrez に複数 Ensembl がある場合は合計カウントが最大の行を残す
    d = (d.assign(_tot=d[samples].sum(axis=1))
          .sort_values("_tot", ascending=False)
          .drop_duplicates("entrez"))
    raw = {}
    for s in samples:
        x = d[s].to_numpy(float)
        sf_inv = x[x > 0].min()
        r = x / sf_inv
        err = np.abs(r - np.round(r)).max()
        if err > 1e-6:
            sys.exit(f"{s}: 生カウントへの復元に失敗（最大誤差 {err:.2e}）")
        raw[s] = np.round(r).astype(int)
    counts = pd.DataFrame(raw, index=d["entrez"].to_numpy())
    info = d.set_index("entrez")[["hgnc_symbol", "ensembl_gene_id"]]
    print(f"[復元] {counts.shape[0]} 遺伝子 × {counts.shape[1]} サンプル（整数誤差 OK）")
    return counts, info


def run_deseq(counts, info, t, outdir):
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats

    cols = [c for c in counts.columns if f"_{t}_" in c]
    sub = counts[cols]
    sub = sub[sub.sum(axis=1) > 0]
    meta = pd.DataFrame({"condition": ["VD" if c.startswith("125D") else "EtOH" for c in cols]},
                        index=cols)
    dds = DeseqDataSet(counts=sub.T, metadata=meta, design="~condition", quiet=True)
    dds.deseq2()
    ds = DeseqStats(dds, contrast=["condition", "VD", "EtOH"], quiet=True)
    ds.summary()
    res = ds.results_df.copy()
    res = res.rename(columns={"log2FoldChange": "lfc_mle"})
    ds.lfc_shrink(coeff="condition[T.VD]")
    res["lfc_shrunk"] = ds.results_df["log2FoldChange"]
    res.index.name = "entrez"
    res = res.join(info)
    path = os.path.join(outdir, f"deg_GSE189984_{t}.csv")
    res.reset_index().to_csv(path, index=False)
    n_sig = int((res["padj"] < 0.05).sum())
    print(f"[DESeq2 {t}] {len(res)} 遺伝子, padj<0.05: {n_sig} "
          f"(down {int(((res['padj'] < 0.05) & (res['lfc_mle'] < 0)).sum())})")
    return path


def build_ortho(mgi, path):
    m = pd.read_csv(mgi, sep="\t", dtype=str)
    key = "DB Class Key"
    hu = m[m["NCBI Taxon ID"] == "9606"][[key, "EntrezGene ID"]].rename(
        columns={"EntrezGene ID": "human_entrez"})
    mo = m[m["NCBI Taxon ID"] == "10090"][[key, "Symbol", "Mouse MGI ID"]].rename(
        columns={"Symbol": "mouse_symbol", "Mouse MGI ID": "mgi_id"})
    pairs = hu.merge(mo, on=key).dropna(subset=["human_entrez", "mouse_symbol"])
    pairs[["human_entrez", "mouse_symbol"]].drop_duplicates().to_csv(path, index=False)
    return mo.dropna().drop_duplicates("mgi_id").set_index("mgi_id")["mouse_symbol"]


def harmonize_impc(impc, mgi_symbol, path):
    v = pd.read_csv(impc)
    cur = v["Gene Accession Id"].map(mgi_symbol)
    n_changed = int((cur.notna() & (cur != v["Gene Symbol"])).sum())
    v["Gene Symbol"] = cur.fillna(v["Gene Symbol"])
    v.to_csv(path, index=False)
    print(f"[IMPC] {v['Gene Symbol'].nunique()} 遺伝子, MGI 現行シンボルへ更新 {n_changed} 行")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--counts-xlsx", required=True)
    ap.add_argument("--impc", required=True)
    ap.add_argument("--mgi", required=True)
    ap.add_argument("--out", default="results/lethal_enrichment")
    ap.add_argument("--n-perm", type=int, default=10000)
    args = ap.parse_args()

    work = os.path.join(args.out, "inputs")
    os.makedirs(work, exist_ok=True)

    counts, info = recover_counts(args.counts_xlsx)
    ortho = os.path.join(work, "ortho_human_entrez_to_mouse.csv")
    mgi_symbol = build_ortho(args.mgi, ortho)
    impc = os.path.join(work, "impc_viability_mgi_symbols.csv")
    harmonize_impc(args.impc, mgi_symbol, impc)

    for t in TIMES:
        deg = run_deseq(counts, info, t, work)
        for lfc_col in ("lfc_shrunk", "lfc_mle"):
            out = os.path.join(args.out, f"GSE189984_{t}_{lfc_col}")
            cmd = [sys.executable, os.path.join(HERE, "lethal_enrichment.py"),
                   "--deg", deg, "--deg-gene", "entrez", "--deg-lfc", lfc_col,
                   "--deg-expr", "baseMean",
                   "--impc", impc, "--impc-gene", "Gene Symbol",
                   "--impc-call", "Viability Phenotype HOMs/HEMIs",
                   "--ortho", ortho, "--ortho-human", "human_entrez",
                   "--ortho-mouse", "mouse_symbol",
                   "--n-perm", str(args.n_perm), "--out", out]
            with open(os.path.join(out + ".log"), "w") as log:
                subprocess.run(cmd, check=True, stdout=log, stderr=subprocess.STDOUT)
            print(f"[lethal] {t} {lfc_col} → {out}/")


if __name__ == "__main__":
    main()
