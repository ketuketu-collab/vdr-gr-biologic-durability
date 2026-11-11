#!/usr/bin/env python3
"""
標的遺伝子の「細胞種特異的（条件付き）KO」で致死・早期死亡になるかを MGI から抽出する。

T 細胞全体（Cd4-cre）では平気だが Treg（Foxp3-cre）だけ消すと死ぬ、のような
「一分画だけで致死」を眺めるための表を作る。

判定
----
- MGI_PhenoGenoMP.rpt の遺伝子型のうち、
    * Cre アレル（Tg(X-cre) や X<tm(cre)>、icre / cre/ERT2 を含む）を 1 つ持ち、
    * 標的遺伝子の座位に野生型アレルを含まない（flox/flox, flox/null など）
  ものを条件付き KO とする。Cre の駆動遺伝子 X を「ドライバー」とする。
- 全身・生殖系列で働くドライバー（EIIa, Sox2, CMV, ROSA26-CreERT2, UBC-CreERT2,
  CAG-CreER など）は「一分画」ではないので is_global_driver=True として区別する。
- 致死 = preweaning lethality / lethality complete penetrance / lethality at weaning
  とその下位語。早期死亡 = MP:0002083 premature death とその下位語
  （Treg 欠損のように離乳後に自己免疫で死ぬ型はこちらに入る）。

出力
----
conditional_ko_by_driver.csv : 遺伝子 × ドライバーごとの判定
conditional_ko_summary.csv   : 遺伝子ごとに「死ぬドライバー / 死なないドライバー」を並べた表
"""

import argparse
import os
import re
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ko_lethality_annotation import (LETHAL_ROOTS, PREMATURE_DEATH_ROOT,  # noqa: E402
                                     load_mp, split_genotype)

CRE_RE = re.compile(r"[(/\-]i?cre(?:/ert2?|ert2?|/esr1\*?|er)?\b", re.I)
GLOBAL_DRIVERS = {
    "EIIA", "ELLA", "SOX2", "CMV", "ACTB", "HPRT", "HPRT1", "ZP3", "DDX4", "STRA8",
    "GT(ROSA)26SOR", "UBC", "CAG", "UBIQUITIN", "PGK1", "MEOX2", "NANOG", "POU5F1",
    "POLR2A", "ACTB", "EEF1A1", "H2AZ2",
}
# よく使われるドライバーの対象細胞（未登録のものはドライバー名だけ表示）
DRIVER_CELL = {
    "FOXP3": "Treg", "CD4": "T細胞(CD4+CD8+期以降)", "LCK": "T細胞(胸腺)",
    "CD8A": "CD8 T細胞", "E8I": "CD8 T細胞", "CD2": "T/NK細胞", "TNFRSF4": "活性化T/Treg",
    "RORC": "RORγt+ (Th17/ILC3/DP胸腺細胞)", "GZMB": "細胞傷害性T/NK", "NCR1": "NK細胞",
    "CD19": "B細胞", "CD79A": "B細胞", "MB1": "B細胞", "CR2": "成熟B細胞", "AICDA": "胚中心B細胞",
    "CD21": "成熟B細胞", "PRDM1": "形質細胞", "LYZ2": "骨髄系(マクロファージ/好中球)",
    "LYSM": "骨髄系(マクロファージ/好中球)", "CSF1R": "マクロファージ/単球", "CX3CR1": "マクロファージ/ミクログリア",
    "ITGAX": "樹状細胞", "CD11C": "樹状細胞", "MRP8": "好中球", "S100A8": "好中球", "MCPT8": "好塩基球",
    "CPA3": "マスト細胞", "VAV1": "全血球系", "VAV": "全血球系", "MX1": "血球系/肝(誘導)",
    "PTPRC": "全血球系", "TEK": "内皮+血球系", "TIE2": "内皮+血球系", "CDH5": "内皮", "PDGFB": "内皮",
    "LYVE1": "リンパ管内皮/一部マクロファージ", "PROX1": "リンパ管内皮", "VIL1": "腸上皮", "VIL": "腸上皮",
    "ALB": "肝細胞", "KRT14": "表皮基底", "KRT5": "表皮基底", "NES": "神経前駆", "SYN1": "神経細胞",
    "CAMK2A": "前脳神経細胞", "GFAP": "アストロサイト", "MYH6": "心筋", "MYL7": "心房筋",
    "INS2": "膵β細胞", "PDX1": "膵", "COL1A1": "線維芽/骨芽", "COL1A2": "線維芽/骨芽",
    "SFTPC": "肺胞II型上皮", "NKX2-1": "肺/甲状腺上皮", "PDGFRB": "周皮細胞/間質", "TAGLN": "平滑筋",
    "MYH11": "平滑筋", "ACTA2": "平滑筋", "ADIPOQ": "脂肪細胞", "LEPR": "骨髄間質", "PF4": "巨核球/血小板",
    "KIT": "c-Kit+", "ALB1": "肝細胞", "CKMM": "骨格筋/心筋", "MCK": "骨格筋/心筋",
    "MMTV": "乳腺など(MMTV)", "WAP": "乳腺", "COL2A1": "軟骨", "BGLAP": "骨芽細胞", "SLC6A3": "ドパミン神経",
    "MYL2": "心室筋", "TYR": "メラノサイト", "NKX2-5": "心臓前駆", "PRRX1": "四肢間葉", "FABP4": "脂肪/マクロファージ",
    "SCGB1A1": "気道クラブ細胞", "MPZ": "シュワン細胞", "SOX9": "Sox9+", "WNT1": "神経堤", "SHH": "Shh+", "LGR5": "腸幹細胞", "OLIG2": "オリゴデンドロサイト",
}


def allele_gene(sym):
    if sym.startswith("Tg("):
        inner = sym[3:sym.index(")")] if ")" in sym else sym[3:]
        return inner.split("-")[0].split("/")[0]
    return sym.split("<")[0]


def driver_of(sym):
    # 別の座位に挿入されたトランスジーン（例: Speer6-ps1<Tg(Alb-cre)21Mgn>）は中の Tg(...) を使う
    m = re.search(r"Tg\(([^)]*)\)", sym)
    if m:
        return m.group(1).split("-")[0].split("/")[0].upper()
    return allele_gene(sym).upper()


def parse_genotype(comp):
    loci = []
    depth, buf = 0, ""
    for ch in comp:
        if ch == "<":
            depth += 1
        elif ch == ">":
            depth -= 1
        if ch == "," and depth == 0:
            loci.append(buf)
            buf = ""
        else:
            buf += ch
    loci.append(buf)
    return [split_genotype(l.strip()) for l in loci if l.strip()]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--orthologs", default="results/ko_lethality/ko_lethality_per_ortholog.csv")
    ap.add_argument("--mgi-dir", required=True)
    ap.add_argument("--out", default="results/ko_lethality")
    args = ap.parse_args()

    names, desc, alt = load_mp(os.path.join(args.mgi_dir, "MPheno_OBO.ontology"))
    lethal = set().union(*(desc(r) for r in LETHAL_ROOTS))
    # 「腫瘍の出現が早い」は致死ではないので早期死亡から外す
    premature = desc(PREMATURE_DEATH_ROOT) - desc("MP:0014237")

    orth = pd.read_csv(args.orthologs).dropna(subset=["mgi_id"])
    target = {r.mouse_ortholog: (r.gene, r.mgi_id, r.ko_lethal_call) for r in orth.itertuples()}

    gp = pd.read_csv(os.path.join(args.mgi_dir, "MGI_PhenoGenoMP.rpt"), sep="\t", header=None,
                     dtype=str, names=["genotype", "alleles", "background", "mp", "pmid",
                                       "marker_ids", "genotype_id"])
    gp = gp[gp["genotype"].str.contains(CRE_RE, na=False)]
    gp["mp"] = gp["mp"].map(lambda t: alt.get(t, t))

    # 標的側は flox（Conditional ready）か null のアレルに限る。
    # 活性化型・ヒト化・ノックイン（Cre で発現 ON になる gain-of-function など）は除く。
    attr = pd.read_csv(os.path.join(args.mgi_dir, "MGI_PhenotypicAllele.rpt"), sep="\t",
                       comment="#", header=None, dtype=str).set_index(1)[4].fillna("")
    bad = ("Humanized sequence", "Inserted expressed sequence", "Activated",
           "Dominant negative", "Modified isoform", "Modified regulatory region", "Hypomorph")

    # MGI の属性に Activated が付いていないが、実際は Cre で活性化する gain-of-function のアレル
    known_gof = {"Braf<tm1Mmcm>"}   # BRAF V600E (Braf CA; Dankort 2007)

    def ko_allele(a):
        if a in known_gof:
            return False
        t = attr.get(a, "")
        return (("Conditional ready" in t or "Null/knockout" in t)
                and not any(b in t for b in bad))

    rows = []
    for gid, g in gp.groupby("genotype_id"):
        loci = parse_genotype(g["genotype"].iloc[0])
        cre = [a for l in loci for a in l if CRE_RE.search(a)]
        if len({driver_of(a) for a in cre}) != 1:
            continue
        drv = driver_of(cre[0])
        inducible = bool(re.search(r"ert|esr1|cre/?er\b", cre[0], re.I))
        def benign(locus):
            # Cre 以外の座位で許すのは野生型・トランスジーン片側・レポーターのみ
            for a in locus:
                if a.endswith("<+>") or a in ("0", "Y"):
                    continue
                t = attr.get(a, "")
                if "Reporter" in t and not any(x in t for x in ("Null/knockout", "Activated",
                                                                 "Dominant negative")):
                    continue
                return False
            return True

        for l in loci:
            if any(CRE_RE.search(a) for a in l) or len(l) != 2:
                continue
            others = [o for o in loci if o is not l and not any(CRE_RE.search(a) for a in o)]
            if not all(benign(o) for o in others):
                continue          # 他の遺伝子の変異を含む（二重 KO など）
            gene = allele_gene(l[0])
            if gene not in target or allele_gene(l[1]) != gene:
                continue
            if any(a.endswith("<+>") for a in l):
                continue          # 片側が野生型 → ヘテロなので除外
            if not all(ko_allele(a) for a in l):
                continue
            if not any("Conditional ready" in attr.get(a, "") for a in l):
                continue          # flox を含まない＝全身 KO に Cre を掛けただけ
            mps = set(g["mp"])
            hg, mid, call = target[gene]
            rows.append({
                "gene": hg, "mouse": gene, "global_ko_call": call,
                "driver": drv, "cell_type": DRIVER_CELL.get(drv, ""),
                "inducible": inducible, "is_global_driver": drv in GLOBAL_DRIVERS,
                "lethal": bool(mps & lethal), "premature_death": bool(mps & premature),
                "lethal_terms": ";".join(sorted({names.get(t, t) for t in mps & (lethal | premature)})),
                "genotype": g["genotype"].iloc[0], "genotype_id": gid,
                "pmids": ";".join(sorted({p for p in g["pmid"].dropna() for p in p.split(",")})),
            })

    d = pd.DataFrame(rows)
    by = (d.groupby(["gene", "mouse", "global_ko_call", "driver", "cell_type",
                     "is_global_driver"], as_index=False)
            .agg(n_genotypes=("genotype_id", "nunique"), inducible=("inducible", "any"),
                 lethal=("lethal", "any"), premature_death=("premature_death", "any"),
                 lethal_terms=("lethal_terms", lambda s: ";".join(sorted({t for x in s for t in x.split(";") if t}))),
                 pmids=("pmids", lambda s: ";".join(sorted({t for x in s for t in x.split(";") if t})))))
    by["dies"] = by["lethal"] | by["premature_death"]
    by.to_csv(os.path.join(args.out, "conditional_ko_by_driver.csv"), index=False)

    def lab(r):
        s = r.driver.capitalize() if r.driver.isupper() else r.driver
        s += f"({r.cell_type})" if r.cell_type else ""
        s += "[誘導]" if r.inducible else ""
        if r.lethal:
            s += "[致死(一部)]" if "incomplete" in str(r.lethal_terms) and "complete penetrance" not in str(r.lethal_terms).replace("incomplete penetrance", "") else "[致死]"
        elif r.premature_death:
            s += "[早期死亡]"
        return s

    sub = by[~by["is_global_driver"]].copy()
    sub["label"] = sub.apply(lab, axis=1)
    summ = (sub.groupby(["gene", "mouse", "global_ko_call"])
               .apply(lambda g: pd.Series({
                   "dies_drivers": ", ".join(g.loc[g.dies, "label"]),
                   "survives_drivers": ", ".join(g.loc[~g.dies, "label"]),
                   "n_dies": int(g.dies.sum()), "n_drivers": len(g)}), include_groups=False)
               .reset_index())
    summ = summ[summ["n_dies"] > 0].sort_values(["global_ko_call", "gene"])
    summ.to_csv(os.path.join(args.out, "conditional_ko_summary.csv"), index=False,
                encoding="utf-8-sig")
    print(f"条件付き KO データのある標的遺伝子: {sub['gene'].nunique()}")
    print(f"一分画 KO で致死/早期死亡があるもの: {len(summ)}")
    print(summ.groupby("global_ko_call").size().to_string())
    write_list(args.out, summ, by)
    return summ


IMMUNE_DRIVERS = {"FOXP3", "CD4", "LCK", "CD8A", "E8I", "CD2", "TNFRSF4", "RORC", "CD19", "CD79A",
                  "CR2", "AICDA", "LYZ2", "CSF1R", "CX3CR1", "ITGAX", "VAV1", "NCR1", "S100A8",
                  "PTPRC", "GZMB", "MCPT8", "CPA3", "PRDM1", "MX1", "TEK"}


def write_list(outdir, summ, by):
    """眺める用の一覧（conditional_ko_list.md）と、日本語列名の summary を書く。"""
    pg = pd.read_csv(os.path.join(outdir, "ko_lethality_per_gene.csv")).set_index("gene")
    lab = {"lethal": "致死", "lethal_partial": "一部致死", "viable": "生存", "no_ko_data": "データなし"}
    s = summ.fillna("").copy()
    s["全身KO"] = [lab.get(c, c) + ("+離乳後早期死亡" if pg.loc[g, "mgi_ko_premature_death"] is True
                                   or pg.loc[g, "mgi_ko_premature_death"] == True else "")  # noqa: E712
                 for g, c in zip(s["gene"], s["global_ko_call"])]
    s["文献PMID"] = [";".join(sorted({p for x in by[(by.gene == g) & by.dies].pmids.dropna()
                                       for p in str(x).split(";") if p})) for g in s["gene"]]
    s = s.rename(columns={"gene": "遺伝子", "mouse": "マウス", "dies_drivers": "死ぬ Cre（対象細胞）",
                          "survives_drivers": "死なない Cre"})
    cols = ["遺伝子", "マウス", "全身KO", "死ぬ Cre（対象細胞）", "死なない Cre", "文献PMID"]
    s[cols].to_csv(os.path.join(outdir, "conditional_ko_summary.csv"), index=False, encoding="utf-8-sig")
    A = s[s["global_ko_call"] != "lethal"][cols]
    C = s[s["global_ko_call"] == "lethal"][cols]

    i = by[by.driver.isin(IMMUNE_DRIVERS) & ~by.is_global_driver]
    i = i.groupby("gene").filter(lambda g: g.dies.any() and (~g.dies).any())

    def fmt(g, m):
        return ", ".join(f"{d.capitalize()}({c})" if c else d.capitalize()
                         for d, c in zip(g.loc[m, "driver"], g.loc[m, "cell_type"].fillna("")))
    rows = [{"遺伝子": gene, "死ぬ免疫系 Cre": fmt(g, g.dies), "死なない免疫系 Cre": fmt(g, ~g.dies)}
            for gene, g in i.groupby("gene")]
    B = pd.DataFrame(rows, columns=["遺伝子", "死ぬ免疫系 Cre", "死なない免疫系 Cre"])

    with open(os.path.join(outdir, "conditional_ko_list.md"), "w", encoding="utf-8") as f:
        f.write("# 細胞種特異的 KO（Cre/flox）で死ぬ標的遺伝子\n\n"
                "MGI の条件付き KO。標的遺伝子が flox/flox か flox/null、Cre が 1 種類、他の遺伝子の"
                "変異を含まないものに限る。全身・生殖系列ドライバー（EIIa, Sox2, CMV, ROSA26-CreERT2 など）"
                "は除外。[致死] は離乳前致死、[早期死亡] は離乳後の早期死亡（MGI の premature death。"
                "腫瘍の早期発生は含めない）。\n\n")
        f.write(f"## A. 全身 KO では致死でないのに、一分画の KO で死ぬ（{len(A)}）\n\n"
                + A.to_markdown(index=False) + "\n\n")
        f.write(f"## B. 免疫系の中で、分画によって死ぬ・死なないが分かれる（{len(B)}）\n\n"
                + B.to_markdown(index=False) + "\n\n")
        f.write(f"## C. 全身 KO でも致死で、どの分画の KO で死ぬかが分かっているもの（{len(C)}）\n\n"
                + C.to_markdown(index=False) + "\n\n")
        f.write("## 注意\n\n"
                "- MGI の注釈には実験条件が付いていない。例えば IL4R の Lyz2-cre の早期死亡"
                "（PMID 15142530）は住血吸虫感染下での死亡で、無処置での死亡ではない。\n"
                "- 「死なない Cre」は、その組合せで致死・早期死亡の注釈がないという意味で、"
                "生存を確認したとは限らない。\n"
                "- Cre ラインには異所性発現があるもの（例: Gfap-cre 8Gtm は神経前駆にも入る）があり、"
                "対象細胞の欄は代表的な発現部位。\n"
                "- MGI で Activated と付いていない gain-of-function アレルは自動では除けない。"
                "Braf<tm1Mmcm>（BRAF V600E）は手で除外した。\n")


if __name__ == "__main__":
    main()
