#!/usr/bin/env python3
"""
標的遺伝子リスト（ヒト）について、マウス KO で致死になるかを注釈する。

情報源（2 つを統合）
--------------------
1. IMPC viability.csv
   統一プロトコルでの離乳前生存判定（lethal / subviable / viable）。
2. MGI MGI_GenePheno.rpt + MGI_PhenotypicAllele.rpt + MPheno_OBO.ontology
   文献由来の KO 表現型。1 遺伝子の遺伝子型のうち、全アレルが
   "Null/knockout" 型の homozygous / hemizygous / compound het を
   「KO」とみなす（条件付き KO・点変異・ノックインは除く）。
   致死 = MP:0010770 preweaning lethality（prenatal / perinatal を含む）、
          MP:0011400 lethality complete penetrance、
          MP:0008569 lethality at weaning とその下位語。
   MP:0002083 premature death（成体での早期死亡）は致死に含めず別欄。

オルソログ: MGI HOM_MouseHumanSequence.rpt（ヒトシンボル→マウス）。
1:1 でない場合も全マウスオルソログを使い、ortholog_type に明記する。

統合判定 (ko_lethal_call)
-------------------------
  lethal            : IMPC lethal、または MGI KO で完全浸透の致死
  lethal_partial    : IMPC subviable、または MGI KO で不完全浸透・浸透度不明の致死のみ
  viable            : KO の表現型データがあり、致死の報告なし
  no_ko_data        : KO データなし
  no_mouse_ortholog : マウスオルソログなし

使い方
------
python3 scripts/ko_lethality_annotation.py \
    --genes results/remap_scores_expanded.csv --gene-col gene \
    --impc viability.csv.gz --mgi-dir <MGI レポートのディレクトリ> \
    --out results/ko_lethality
"""

import argparse
import collections
import os
import re

import pandas as pd

LETHAL_ROOTS = {"MP:0010770", "MP:0011400", "MP:0008569"}
PRENATAL_ROOT = "MP:0002080"
PERINATAL_ROOT = "MP:0002081"
PREMATURE_DEATH_ROOT = "MP:0002083"


def load_mp(path):
    names, parents = {}, collections.defaultdict(set)
    cur = None
    alt = {}
    for line in open(path, encoding="utf-8"):
        line = line.rstrip("\n")
        if line == "[Term]":
            cur = None
        elif line.startswith("id: "):
            cur = line[4:]
        elif cur and line.startswith("name: "):
            names[cur] = line[6:]
        elif cur and line.startswith("is_a: "):
            parents[cur].add(line[6:16])
        elif cur and line.startswith("alt_id: "):
            alt[line[8:]] = cur
    children = collections.defaultdict(set)
    for c, ps in parents.items():
        for p in ps:
            children[p].add(c)

    def desc(root):
        out, stack = set(), [root]
        while stack:
            t = stack.pop()
            if t not in out:
                out.add(t)
                stack.extend(children[t])
        return out

    return names, desc, alt


def split_genotype(comp):
    """'A<tm1(cre/ERT2)>/A<+>' を '/' で分割（<> 内の / は無視）。"""
    parts, depth, buf = [], 0, ""
    for ch in comp:
        if ch == "<":
            depth += 1
        elif ch == ">":
            depth -= 1
        if ch == "/" and depth == 0:
            parts.append(buf)
            buf = ""
        else:
            buf += ch
    parts.append(buf)
    return [p.strip() for p in parts]


def load_mgi_ko(mgi_dir, names, desc, alt):
    lethal = set().union(*(desc(r) for r in LETHAL_ROOTS))
    prenatal = desc(PRENATAL_ROOT)
    perinatal = desc(PERINATAL_ROOT)
    premature = desc(PREMATURE_DEATH_ROOT)

    alle = pd.read_csv(os.path.join(mgi_dir, "MGI_PhenotypicAllele.rpt"), sep="\t",
                       comment="#", header=None, dtype=str)
    null_alleles = set(alle.loc[alle[4].fillna("").str.contains("Null/knockout"), 1])

    gp = pd.read_csv(os.path.join(mgi_dir, "MGI_GenePheno.rpt"), sep="\t", header=None, dtype=str,
                     names=["genotype", "alleles", "allele_ids", "background", "mp",
                            "pmid", "marker_id", "genotype_id"])
    gp["mp"] = gp["mp"].map(lambda t: alt.get(t, t))

    def ko_zygosity(comp):
        parts = split_genotype(comp)
        if len(parts) != 2:
            return None
        kinds = []
        for p in parts:
            if p in null_alleles:
                kinds.append("null")
            elif p.endswith("<+>"):
                kinds.append("wt")
            elif p in ("Y", "0"):
                kinds.append("hemi")
            else:
                return None          # null 以外の変異アレルを含む → KO とみなさない
        nn = kinds.count("null")
        if nn == 2:
            return "hom"
        if nn == 1 and "hemi" in kinds:
            return "hemi"
        if nn == 1 and "wt" in kinds:
            return "het"
        return None

    gp["zyg"] = gp["genotype"].map(ko_zygosity)
    ko = gp[gp["zyg"].notna()].copy()

    def penetrance(t):
        n = names.get(t, "")
        if "incomplete penetrance" in n:
            return "incomplete"
        if "complete penetrance" in n:
            return "complete"
        return "unspecified"

    ko["is_lethal"] = ko["mp"].isin(lethal)
    ko["pen"] = ko["mp"].map(penetrance)

    rows = {}
    for mid, g in ko.groupby("marker_id"):
        homs = g[g["zyg"].isin(["hom", "hemi"])]
        hets = g[g["zyg"] == "het"]
        hl = homs[homs["is_lethal"]]
        stages = []
        if hl["mp"].isin(prenatal).any():
            stages.append("prenatal")
        if hl["mp"].isin(perinatal).any():
            stages.append("perinatal")
        if len(hl) and not stages:
            stages.append("postnatal/preweaning")
        rows[mid] = {
            "mgi_ko_hom_genotypes": homs["genotype_id"].nunique(),
            "mgi_ko_hom_lethal": bool(len(hl)),
            "mgi_ko_hom_lethal_complete": bool((hl["pen"] == "complete").any()),
            "mgi_lethal_stage": ";".join(stages),
            "mgi_lethal_terms": ";".join(sorted({names.get(t, t) for t in hl["mp"]})),
            "mgi_lethal_pmids": ";".join(sorted({p for p in hl["pmid"].dropna()
                                                 for p in str(p).split(",") if p})),
            "mgi_ko_het_lethal": bool(hets["is_lethal"].any()),
            "mgi_ko_premature_death": bool(homs["mp"].isin(premature).any()),
            # 文献（PubMed 付き）注釈だけでの判定。IMPC 由来の注釈には PubMed が無い
            "lit_ko_hom_genotypes": homs.loc[homs["pmid"].notna(), "genotype_id"].nunique(),
            "lit_ko_hom_lethal": bool(hl["pmid"].notna().any()),
        }
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--genes", required=True)
    ap.add_argument("--gene-col", default="gene")
    ap.add_argument("--impc", required=True)
    ap.add_argument("--mgi-dir", required=True)
    ap.add_argument("--out", default="results/ko_lethality")
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    genes = pd.read_csv(args.genes)
    human = genes[args.gene_col].astype(str).str.strip()

    # オルソログ
    hom = pd.read_csv(os.path.join(args.mgi_dir, "HOM_MouseHumanSequence.rpt"), sep="\t", dtype=str)
    key = "DB Class Key"
    hu = hom[hom["NCBI Taxon ID"] == "9606"][[key, "Symbol"]].rename(columns={"Symbol": "human"})
    mo = hom[hom["NCBI Taxon ID"] == "10090"][[key, "Symbol", "Mouse MGI ID"]].rename(
        columns={"Symbol": "mouse", "Mouse MGI ID": "mgi_id"})
    pairs = hu.merge(mo, on=key)[["human", "mouse", "mgi_id"]].drop_duplicates()
    h_n = pairs.groupby("human")["mouse"].transform("nunique")
    m_n = pairs.groupby("mouse")["human"].transform("nunique")
    pairs["ortholog_type"] = [
        "1:1" if a == 1 and b == 1 else ("1:many" if a > 1 and b == 1 else
                                         ("many:1" if a == 1 else "many:many"))
        for a, b in zip(h_n, m_n)]

    # IMPC（MGI ID で引く）
    v = pd.read_csv(args.impc)
    v["call"] = v["Viability Phenotype HOMs/HEMIs"].str.lower()
    rank = {"lethal": 2, "subviable": 1, "viable": 0}
    impc = {}
    for mid, g in v.groupby("Gene Accession Id"):
        calls = set(g["call"])
        best = max(calls, key=lambda c: rank.get(c, -1))
        impc[mid] = {"impc_call": best,
                     "impc_conflicting": len(calls) > 1 or g["Comment"].notna().any()}

    names, desc, alt = load_mp(os.path.join(args.mgi_dir, "MPheno_OBO.ontology"))
    mgi = load_mgi_ko(args.mgi_dir, names, desc, alt)

    out = []
    for h in human:
        p = pairs[pairs["human"] == h]
        if p.empty:
            out.append({"gene": h, "mouse_ortholog": "", "ortholog_type": "none",
                        "ko_lethal_call": "no_mouse_ortholog"})
            continue
        for _, r in p.iterrows():
            row = {"gene": h, "mouse_ortholog": r["mouse"], "mgi_id": r["mgi_id"],
                   "ortholog_type": r["ortholog_type"]}
            row.update(impc.get(r["mgi_id"], {"impc_call": "", "impc_conflicting": False}))
            row.update(mgi.get(r["mgi_id"], {"mgi_ko_hom_genotypes": 0}))
            ic = row["impc_call"]
            if ic == "lethal" or row.get("mgi_ko_hom_lethal_complete"):
                call = "lethal"
            elif ic == "subviable" or row.get("mgi_ko_hom_lethal"):
                call = "lethal_partial"
            elif ic == "viable" or row["mgi_ko_hom_genotypes"] > 0:
                call = "viable"
            else:
                call = "no_ko_data"
            row["ko_lethal_call"] = call
            # IMPC と文献 KO の食い違い（背景系統・アレル設計の違いで起こりうる）
            lit_n = row.get("lit_ko_hom_genotypes", 0)
            lit_l = row.get("lit_ko_hom_lethal", False)
            conflict = ""
            if ic in ("lethal", "subviable") and lit_n > 0 and not lit_l:
                conflict = "IMPC致死_文献KOは致死報告なし"
            elif ic == "viable" and lit_l:
                conflict = "IMPC生存_文献KOは致死"
            row["impc_vs_literature"] = conflict
            out.append(row)

    df = pd.DataFrame(out)
    extra = [c for c in genes.columns if c != args.gene_col]
    df = df.merge(genes.rename(columns={args.gene_col: "gene"})[["gene"] + extra],
                  on="gene", how="left")
    df.to_csv(os.path.join(args.out, "ko_lethality_per_ortholog.csv"), index=False)

    # 遺伝子単位にまとめる（複数オルソログは最も重い判定）
    sev = {"lethal": 4, "lethal_partial": 3, "viable": 2, "no_ko_data": 1, "no_mouse_ortholog": 0}
    gl = (df.assign(_s=df["ko_lethal_call"].map(sev))
            .sort_values("_s", ascending=False).drop_duplicates("gene").drop(columns="_s"))
    gl = gl.set_index("gene").loc[pd.unique(human)].reset_index()
    gl.to_csv(os.path.join(args.out, "ko_lethality_per_gene.csv"), index=False)

    print(f"{len(gl)} 遺伝子")
    print(gl["ko_lethal_call"].value_counts().to_string())
    return gl


if __name__ == "__main__":
    main()
