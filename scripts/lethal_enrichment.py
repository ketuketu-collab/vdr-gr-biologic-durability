#!/usr/bin/env python3
"""
「発現が下がった遺伝子は KO 致死遺伝子に濃縮しているか」を検定する。

設計の要点
----------
1. 主結果はカットオフなしの順位検定（Mann-Whitney + GSEA 風の running score）。
   小さい fold change を捨てないため。
2. hypergeometric は感度解析として複数カットオフで回す。
3. 背景集合 = IMPC で viability 判定済み ∩ 自データで検出された遺伝子。
4. down 側と up 側の濃縮を別々に出し、非対称性を明示する。

必要な入力
----------
--deg   : 自分の発現解析結果 (CSV/TSV)。遺伝子名列と logFC 列が必要。
--impc  : IMPC viability レポート (CSV)。
          https://www.mousephenotype.org/data/release の
          "Viability" / "genotype-phenotype" 系ファイル。
          遺伝子シンボル列と viability 判定列（Lethal / Subviable / Viable）が必要。
--ortho : (任意) ヒトデータの場合のヒト→マウス オルソログ対応表。

使い方
------
python3 lethal_enrichment.py \
    --deg my_deg.csv --deg-gene gene --deg-lfc log2FoldChange \
    --impc impc_viability.csv --impc-gene gene_symbol --impc-call viability \
    --out results/
"""

import argparse
import os
import sys

import numpy as np
import pandas as pd
from scipy import stats


# ---------------------------------------------------------------- 入力の読み込み

def read_table(path):
    """拡張子から区切り文字を推定して読む。"""
    sep = "\t" if path.endswith((".tsv", ".txt")) else ","
    return pd.read_csv(path, sep=sep)


def normalize_symbol(s):
    """遺伝子シンボルの表記揺れを吸収する（大文字化・空白除去）。"""
    return s.astype(str).str.strip().str.upper()


def load_lethal_sets(impc_path, gene_col, call_col):
    """IMPC の viability 判定から lethal / subviable / viable の3集合を作る。

    背景集合は「判定が付いている遺伝子」全体。ここを全遺伝子にすると
    濃縮が過大評価されるので、必ず判定済みのみに限定する。
    """
    impc = read_table(impc_path)
    for col in (gene_col, call_col):
        if col not in impc.columns:
            sys.exit(f"IMPC ファイルに列 '{col}' がありません。実際の列: {list(impc.columns)}")

    impc = impc[[gene_col, call_col]].dropna()
    impc[gene_col] = normalize_symbol(impc[gene_col])
    calls = impc[call_col].astype(str).str.lower()

    # IMPC の表記は "lethal" / "subviable" / "viable" を含む形で来る
    lethal = set(impc.loc[calls.str.contains("lethal"), gene_col])
    subviable = set(impc.loc[calls.str.contains("subviable"), gene_col])
    # "viable" は "subviable" にも部分一致するので差集合で取る
    viable = set(impc.loc[calls.str.contains("viable"), gene_col]) - subviable - lethal
    tested = lethal | subviable | viable

    return {"lethal": lethal, "subviable": subviable,
            "viable": viable, "tested": tested}


def map_to_mouse(df, gene_col, ortho_path, human_col, mouse_col):
    """ヒト遺伝子をマウスオルソログに変換する。1対多は全て残す。"""
    ortho = read_table(ortho_path)
    ortho[human_col] = normalize_symbol(ortho[human_col])
    ortho[mouse_col] = normalize_symbol(ortho[mouse_col])
    merged = df.merge(ortho[[human_col, mouse_col]],
                      left_on=gene_col, right_on=human_col, how="inner")
    merged[gene_col] = merged[mouse_col]
    return merged.drop(columns=[c for c in (human_col, mouse_col) if c in merged.columns])


# ---------------------------------------------------------------- 検定

def rank_test(lfc, is_lethal):
    """Mann-Whitney U: lethal 遺伝子の logFC が非 lethal より低いか。

    カットオフを使わないので、小さい fold change の情報を捨てない。
    これが主結果。
    """
    a = lfc[is_lethal]
    b = lfc[~is_lethal]
    if len(a) < 3 or len(b) < 3:
        return None

    u, p_less = stats.mannwhitneyu(a, b, alternative="less")      # lethal が下側
    _, p_greater = stats.mannwhitneyu(a, b, alternative="greater")  # lethal が上側

    # 効果量: rank-biserial correlation（-1 なら完全に下側）
    n1, n2 = len(a), len(b)
    rbc = 2 * u / (n1 * n2) - 1

    return {
        "n_lethal": n1,
        "n_other": n2,
        "median_lfc_lethal": float(np.median(a)),
        "median_lfc_other": float(np.median(b)),
        "rank_biserial_r": float(rbc),
        "p_lethal_down": float(p_less),
        "p_lethal_up": float(p_greater),
    }


def running_score(lfc, is_lethal, n_perm=10000, seed=0):
    """GSEA 風の running enrichment score と置換検定による p 値。

    logFC の昇順（最も下がったものが先頭）に並べ、lethal 遺伝子が
    先頭側に固まっているほど ES が大きくなる。
    """
    order = np.argsort(lfc)            # 昇順 = 下がった順
    hits = is_lethal[order]
    weights = np.abs(lfc[order])

    def es(h):
        hit_w = weights * h
        total = hit_w.sum()
        if total == 0:
            return 0.0
        inc = np.cumsum(hit_w) / total
        dec = np.cumsum(~h) / max((~h).sum(), 1)
        curve = inc - dec
        return float(curve[np.argmax(np.abs(curve))])

    observed = es(hits)

    rng = np.random.default_rng(seed)
    null = np.empty(n_perm)
    for i in range(n_perm):
        null[i] = es(rng.permutation(hits))

    # 片側（下がった側への濃縮）。+1 補正で p=0 を避ける
    p = (np.sum(null >= observed) + 1) / (n_perm + 1)
    return {"ES": observed, "p_perm": float(p), "n_perm": n_perm}


def hypergeom_at(lfc, is_lethal, cut, direction):
    """指定カットオフでの hypergeometric 検定（感度解析用）。"""
    if direction == "down":
        selected = lfc <= -abs(cut)
    else:
        selected = lfc >= abs(cut)

    N = len(lfc)                      # 背景（判定済み ∩ 検出済み）
    K = int(is_lethal.sum())          # 背景中の lethal
    n = int(selected.sum())           # 選ばれた遺伝子数
    k = int((selected & is_lethal).sum())  # そのうち lethal

    if n == 0:
        return None

    # P(X >= k): 過剰濃縮の片側検定
    p = stats.hypergeom.sf(k - 1, N, K, n)
    expected = n * K / N
    return {
        "direction": direction,
        "cutoff": abs(cut),
        "n_selected": n,
        "n_lethal_in_selected": k,
        "expected_lethal": round(expected, 2),
        "fold_enrichment": round(k / expected, 3) if expected > 0 else np.nan,
        "p_hypergeom": float(p),
    }


def asymmetry_test(lfc, is_lethal, cut):
    """down 側と up 側で lethal 比率が違うかを Fisher で直接比較する。

    「上がった側ではなく下がった側に本質がある」の直接の証拠になる。
    """
    down = lfc <= -abs(cut)
    up = lfc >= abs(cut)
    if down.sum() == 0 or up.sum() == 0:
        return None

    table = [[int((down & is_lethal).sum()), int((down & ~is_lethal).sum())],
             [int((up & is_lethal).sum()), int((up & ~is_lethal).sum())]]
    odds, p = stats.fisher_exact(table, alternative="greater")
    return {
        "cutoff": abs(cut),
        "down_lethal": table[0][0], "down_other": table[0][1],
        "up_lethal": table[1][0], "up_other": table[1][1],
        "odds_ratio_down_vs_up": float(odds),
        "p_fisher": float(p),
    }


# ---------------------------------------------------------------- 実行

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--deg", required=True, help="発現解析結果 (CSV/TSV)")
    ap.add_argument("--deg-gene", default="gene", help="遺伝子シンボルの列名")
    ap.add_argument("--deg-lfc", default="log2FoldChange", help="logFC の列名")
    ap.add_argument("--impc", required=True, help="IMPC viability レポート")
    ap.add_argument("--impc-gene", default="gene_symbol")
    ap.add_argument("--impc-call", default="viability")
    ap.add_argument("--ortho", help="ヒト→マウス オルソログ表 (任意)")
    ap.add_argument("--ortho-human", default="human_symbol")
    ap.add_argument("--ortho-mouse", default="mouse_symbol")
    ap.add_argument("--cutoffs", default="0.1,0.25,0.5,1.0",
                    help="感度解析に使う |logFC| カットオフ")
    ap.add_argument("--out", default="results")
    args = ap.parse_args()

    os.makedirs(args.out, exist_ok=True)

    # --- DEG 読み込み
    deg = read_table(args.deg)
    for col in (args.deg_gene, args.deg_lfc):
        if col not in deg.columns:
            sys.exit(f"DEG ファイルに列 '{col}' がありません。実際の列: {list(deg.columns)}")

    deg = deg[[args.deg_gene, args.deg_lfc]].dropna()
    deg.columns = ["gene", "lfc"]
    deg["gene"] = normalize_symbol(deg["gene"])

    if args.ortho:
        deg = map_to_mouse(deg, "gene", args.ortho,
                           args.ortho_human, args.ortho_mouse)

    # 1遺伝子が複数行ある場合は |logFC| が最大の行を代表にする
    deg = (deg.assign(_abs=deg["lfc"].abs())
              .sort_values("_abs", ascending=False)
              .drop_duplicates("gene")
              .drop(columns="_abs"))

    # --- IMPC 読み込みと背景集合の確定
    sets = load_lethal_sets(args.impc, args.impc_gene, args.impc_call)

    # ここが要点: 背景は「判定済み ∩ 検出済み」
    bg = deg[deg["gene"].isin(sets["tested"])].copy()
    if len(bg) < 100:
        print(f"[警告] 背景集合が {len(bg)} 遺伝子しかありません。"
              "シンボルの表記系（マウス/ヒト、Ensembl ID など）が一致しているか確認してください。",
              file=sys.stderr)

    lfc = bg["lfc"].to_numpy(dtype=float)
    is_lethal = bg["gene"].isin(sets["lethal"]).to_numpy()

    print(f"背景集合: {len(bg)} 遺伝子 "
          f"(IMPC 判定済み {len(sets['tested'])} との積集合)")
    print(f"  うち lethal: {int(is_lethal.sum())}")
    print(f"  うち subviable: {int(bg['gene'].isin(sets['subviable']).sum())}")

    # --- 主結果: 順位検定
    rank = rank_test(lfc, is_lethal)
    rs = running_score(lfc, is_lethal)

    # --- 感度解析: 複数カットオフでの hypergeometric
    cutoffs = [float(c) for c in args.cutoffs.split(",")]
    hyper_rows = []
    for c in cutoffs:
        for d in ("down", "up"):
            r = hypergeom_at(lfc, is_lethal, c, d)
            if r:
                hyper_rows.append(r)

    # --- 非対称性
    asym_rows = [r for r in (asymmetry_test(lfc, is_lethal, c) for c in cutoffs) if r]

    # --- 出力
    bg.assign(is_lethal=is_lethal).to_csv(
        os.path.join(args.out, "background_annotated.csv"), index=False)
    pd.DataFrame(hyper_rows).to_csv(
        os.path.join(args.out, "hypergeometric_sensitivity.csv"), index=False)
    pd.DataFrame(asym_rows).to_csv(
        os.path.join(args.out, "down_vs_up_asymmetry.csv"), index=False)

    print("\n=== 主結果: 順位検定（カットオフなし） ===")
    if rank:
        for k, v in rank.items():
            print(f"  {k}: {v}")
    print(f"  running ES: {rs['ES']:.4f}  p(置換, {rs['n_perm']}回): {rs['p_perm']:.5f}")

    print("\n=== 感度解析: hypergeometric ===")
    if hyper_rows:
        print(pd.DataFrame(hyper_rows).to_string(index=False))

    print("\n=== down vs up の非対称性 (Fisher) ===")
    if asym_rows:
        print(pd.DataFrame(asym_rows).to_string(index=False))

    print(f"\n出力先: {args.out}/")


if __name__ == "__main__":
    main()
