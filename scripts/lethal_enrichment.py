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
5. 発現量の交絡を制御する。致死遺伝子は高発現のハウスキーピング遺伝子に
   偏っており、logFC の分布（LFC 縮小の強さ・増殖低下による一斉低下）は
   発現量に依存する。--deg-expr を与えると、発現量で層別した置換検定
   （ラベルを同じ発現量ビンの中だけで入れ替える）で全統計量を再評価する。
6. subviable の扱いを明示する（既定は解析から除外）。3 通りの扱いで
   順位検定を回した感度表も出す。
7. オルソログは既定で 1:1 のみ使う（1 対多で同じ値が複製されるのを防ぐ）。
   重複行は |logFC| 最大ではなく、最も発現の高い行（発現量なしなら中央値）
   を代表にする（極端値への偏りを防ぐ）。

必要な入力
----------
--deg      : 自分の発現解析結果 (CSV/TSV)。遺伝子名列と logFC 列が必要。
--deg-expr : (推奨) 発現量の列（DESeq2 の baseMean、edgeR の logCPM など）。
             順位で層別するので、対数・線形どちらのスケールでもよい。
--impc     : IMPC viability レポート (CSV/TSV)。
             https://www.mousephenotype.org/data/release の
             "Viability" / "genotype-phenotype" 系ファイル。
             遺伝子シンボル列と viability 判定列（Lethal / Subviable / Viable を
             含む文字列。"Homozygous - Lethal" のような形でもよい）が必要。
             指定した列名が無い場合は、よくある列名（marker_symbol など）を
             大文字小文字を無視して探す。
--ortho    : (任意) ヒトデータの場合のヒト→マウス オルソログ対応表。

使い方
------
python3 lethal_enrichment.py \
    --deg my_deg.csv --deg-gene gene --deg-lfc log2FoldChange --deg-expr baseMean \
    --impc impc_viability.csv --impc-gene gene_symbol --impc-call viability \
    --out results/
"""

import argparse
import os
import re
import sys

import numpy as np
import pandas as pd
from scipy import stats


IMPC_GENE_CANDIDATES = ["gene_symbol", "marker_symbol", "gene symbol",
                        "marker symbol", "symbol", "gene"]
IMPC_CALL_CANDIDATES = ["viability", "viability phenotype homs/hemis",
                        "viability_call", "outcome", "call"]

SUBVIABLE_MODES = ("exclude", "other", "lethal")


# ---------------------------------------------------------------- 入力の読み込み

def read_table(path):
    """拡張子から区切り文字を推定して読む（.gz も可）。"""
    base = path[:-3] if path.endswith(".gz") else path
    sep = "\t" if base.endswith((".tsv", ".txt")) else ","
    return pd.read_csv(path, sep=sep)


def normalize_symbol(s):
    """遺伝子シンボルの表記揺れを吸収する（大文字化・空白除去）。"""
    return s.astype(str).str.strip().str.upper()


def _norm_colname(c):
    return re.sub(r"[\s_]+", " ", str(c)).strip().lower()


def resolve_col(df, requested, candidates, label):
    """列名を解決する。指定名 → 大文字小文字・空白/_ を無視した一致 → 候補名の順。"""
    if requested in df.columns:
        return requested
    lookup = {_norm_colname(c): c for c in df.columns}
    for name in [requested] + list(candidates):
        hit = lookup.get(_norm_colname(name))
        if hit is not None:
            print(f"[情報] {label} の列として '{hit}' を使います（指定: '{requested}'）。",
                  file=sys.stderr)
            return hit
    sys.exit(f"{label} の列 '{requested}' が見つかりません。実際の列: {list(df.columns)}")


def load_lethal_sets(impc_path, gene_col, call_col):
    """IMPC の viability 判定から lethal / subviable / viable の3集合を作る。

    背景集合は「判定が付いている遺伝子」全体。ここを全遺伝子にすると
    濃縮が過大評価されるので、必ず判定済みのみに限定する。

    1遺伝子に複数の判定（アレル・センター・性別違い）がある場合は
    lethal > subviable > viable の優先順位で1つにまとめ、件数を報告する。
    """
    impc = read_table(impc_path)
    gene_col = resolve_col(impc, gene_col, IMPC_GENE_CANDIDATES, "IMPC 遺伝子")
    call_col = resolve_col(impc, call_col, IMPC_CALL_CANDIDATES, "IMPC 判定")

    impc = impc[[gene_col, call_col]].dropna().copy()
    impc["gene"] = normalize_symbol(impc[gene_col])
    calls = impc[call_col].astype(str).str.lower()

    # "viable" は "subviable" にも部分一致するので、語として判定する
    is_sub = calls.str.contains("subviable")
    is_let = calls.str.contains("lethal")
    is_via = calls.str.contains(r"(?<!sub)viable", regex=True)

    lethal_raw = set(impc.loc[is_let, "gene"])
    sub_raw = set(impc.loc[is_sub, "gene"])
    via_raw = set(impc.loc[is_via, "gene"])

    lethal = lethal_raw
    subviable = sub_raw - lethal
    viable = via_raw - lethal - subviable
    tested = lethal | subviable | viable

    n_conflict = len((lethal_raw & sub_raw) | (lethal_raw & via_raw) | (sub_raw & via_raw))
    if n_conflict:
        print(f"[情報] 判定が複数ある遺伝子 {n_conflict} 個を "
              "lethal > subviable > viable の優先順位でまとめました。", file=sys.stderr)

    return {"lethal": lethal, "subviable": subviable,
            "viable": viable, "tested": tested}


def map_to_mouse(df, gene_col, ortho_path, human_col, mouse_col, mode="one2one"):
    """ヒト遺伝子をマウスオルソログに変換する。

    mode="one2one": ヒト・マウスともに相手が1つだけの対応のみ使う（既定）。
                    1対多だと同じ logFC が複数のマウス遺伝子に複製され、
                    独立な観測として数えられてしまうため。
    mode="all"    : 全対応を残す（旧挙動）。
    """
    ortho = read_table(ortho_path)
    human_col = resolve_col(ortho, human_col, ["human_symbol", "human gene name"], "オルソログ(ヒト)")
    mouse_col = resolve_col(ortho, mouse_col, ["mouse_symbol", "mouse gene name"], "オルソログ(マウス)")
    ortho = pd.DataFrame({"_h": normalize_symbol(ortho[human_col]),
                          "_m": normalize_symbol(ortho[mouse_col])}).drop_duplicates()

    if mode == "one2one":
        h_n = ortho.groupby("_h")["_m"].transform("nunique")
        m_n = ortho.groupby("_m")["_h"].transform("nunique")
        n_before = len(ortho)
        ortho = ortho[(h_n == 1) & (m_n == 1)]
        print(f"[情報] オルソログ 1:1 のみ使用: {len(ortho)}/{n_before} 対応", file=sys.stderr)

    merged = df.merge(ortho, left_on=gene_col, right_on="_h", how="inner")
    merged[gene_col] = merged["_m"]
    return merged.drop(columns=["_h", "_m"])


def collapse_duplicates(deg):
    """1遺伝子が複数行ある場合に代表値を決める。

    |logFC| 最大を選ぶと極端値に偏るので、発現量があれば最も発現の高い行、
    無ければ logFC の中央値を使う。
    """
    n_dup = int(deg["gene"].duplicated().sum())
    if n_dup == 0:
        return deg
    if "expr" in deg.columns:
        out = (deg.sort_values("expr", ascending=False)
                  .drop_duplicates("gene"))
        how = "最も発現の高い行"
    else:
        out = deg.groupby("gene", as_index=False)["lfc"].median()
        how = "logFC の中央値"
    print(f"[情報] 重複行 {n_dup} 行を {how} でまとめました。", file=sys.stderr)
    return out


def assign_labels(bg, sets, subviable_mode):
    """subviable の扱いに応じて解析対象と lethal ラベルを決める。

    exclude: subviable を背景から除く（lethal vs viable のきれいな比較。既定）
    other  : subviable を非 lethal 側に入れる（効果が薄まる）
    lethal : subviable を lethal 側に入れる
    """
    is_sub = bg["gene"].isin(sets["subviable"])
    if subviable_mode == "exclude":
        bg = bg[~is_sub]
        is_lethal = bg["gene"].isin(sets["lethal"])
    elif subviable_mode == "lethal":
        is_lethal = bg["gene"].isin(sets["lethal"]) | is_sub
    else:
        is_lethal = bg["gene"].isin(sets["lethal"])
    return bg, is_lethal.to_numpy()


def make_strata(expr, n_bins):
    """発現量の順位で等頻度ビンを作る。スケール（log/線形）に依存しない。"""
    if expr is None or n_bins <= 1:
        return np.zeros(len(expr) if expr is not None else 0, dtype=int)
    r = stats.rankdata(expr, method="ordinal")
    return np.minimum(((r - 1) * n_bins // len(expr)).astype(int), n_bins - 1)


# ---------------------------------------------------------------- 検定（解析解）

def rank_test(lfc, is_lethal):
    """Mann-Whitney U: lethal 遺伝子の logFC が非 lethal より低いか。

    カットオフを使わないので、小さい fold change の情報を捨てない。
    発現量の交絡は制御しない（層別版は permutation_tests を参照）。
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


def hypergeom_at(lfc, is_lethal, cut, direction):
    """指定カットオフでの hypergeometric 検定（感度解析用）。"""
    selected = lfc <= -abs(cut) if direction == "down" else lfc >= abs(cut)

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


# ---------------------------------------------------------------- 置換検定

def _es(h, weights):
    """GSEA 風 running enrichment score（h は logFC 昇順に並べた lethal ラベル）。

    正なら下がった側、負なら上がった側に lethal が固まっている。
    """
    hit_w = weights * h
    total = hit_w.sum()
    if total == 0:
        return 0.0
    inc = np.cumsum(hit_w) / total
    miss = ~h
    dec = np.cumsum(miss) / max(miss.sum(), 1)
    curve = inc - dec
    return float(curve[np.argmax(np.abs(curve))])


def _perm_p(null, obs, side):
    """+1 補正つきの片側置換 p。"""
    if side == "greater":
        return float((np.sum(null >= obs) + 1) / (len(null) + 1))
    return float((np.sum(null <= obs) + 1) / (len(null) + 1))


def permutation_tests(lfc, is_lethal, strata, cutoffs, n_perm=10000, seed=0):
    """lethal ラベルを strata 内だけで入れ替える置換検定で全統計量を評価する。

    strata が全て同じ値なら通常の（層別なしの）置換検定になる。
    発現量ビンを strata に与えると「発現量が同程度の遺伝子と比べて
    lethal が下がっているか」を問うことになり、発現量の交絡を除ける。

    返り値
    ------
    rank : rank-biserial r（観測値・置換 null の平均・両側それぞれの片側 p）
    es   : running ES（観測値・null 平均・p_down / p_up）
    hyper: カットオフ×方向ごとの lethal 数（観測・null 期待値・p）
    asym : カットオフごとの down−up の lethal 比率差（観測・null 平均・p）
    """
    n = len(lfc)
    n1 = int(is_lethal.sum())
    n2 = n - n1
    ranks = stats.rankdata(lfc)
    order = np.argsort(lfc, kind="stable")
    weights = np.abs(lfc[order])

    sel = {}
    for c in cutoffs:
        sel[(c, "down")] = lfc <= -abs(c)
        sel[(c, "up")] = lfc >= abs(c)

    def rbc_of(h):
        u = ranks[h].sum() - n1 * (n1 + 1) / 2
        return 2 * u / (n1 * n2) - 1

    def asym_of(h, c):
        d, u = sel[(c, "down")], sel[(c, "up")]
        nd, nu = d.sum(), u.sum()
        if nd == 0 or nu == 0:
            return np.nan
        return (h & d).sum() / nd - (h & u).sum() / nu

    def all_stats(h):
        return (rbc_of(h), _es(h[order], weights),
                [int((h & sel[k]).sum()) for k in sel],
                [asym_of(h, c) for c in cutoffs])

    obs_r, obs_es, obs_k, obs_a = all_stats(is_lethal)

    strata_idx = [np.flatnonzero(strata == s) for s in np.unique(strata)]
    rng = np.random.default_rng(seed)
    null_r = np.empty(n_perm)
    null_es = np.empty(n_perm)
    null_k = np.empty((n_perm, len(sel)))
    null_a = np.empty((n_perm, len(cutoffs)))
    perm = is_lethal.copy()
    for i in range(n_perm):
        for idx in strata_idx:
            perm[idx] = rng.permutation(is_lethal[idx])
        r, e, k, a = all_stats(perm)
        null_r[i], null_es[i], null_k[i], null_a[i] = r, e, k, a

    rank = {"rank_biserial_r": float(obs_r),
            "null_mean_r": float(null_r.mean()),
            "r_minus_null": float(obs_r - null_r.mean()),
            "p_lethal_down": _perm_p(null_r, obs_r, "less"),
            "p_lethal_up": _perm_p(null_r, obs_r, "greater")}
    es = {"ES": obs_es, "null_mean_ES": float(null_es.mean()),
          "p_down": _perm_p(null_es, obs_es, "greater"),
          "p_up": _perm_p(null_es, obs_es, "less")}

    hyper = {}
    for j, key in enumerate(sel):
        exp = float(null_k[:, j].mean())
        hyper[key] = {"expected_lethal_perm": round(exp, 2),
                      "fold_enrichment_perm": round(obs_k[j] / exp, 3) if exp > 0 else np.nan,
                      "p_perm": _perm_p(null_k[:, j], obs_k[j], "greater")}

    asym = {}
    for j, c in enumerate(cutoffs):
        col = null_a[:, j]
        if np.isnan(obs_a[j]):
            continue
        asym[c] = {"diff_frac_down_minus_up": float(obs_a[j]),
                   "null_mean_diff": float(np.nanmean(col)),
                   "p_perm": _perm_p(col[~np.isnan(col)], obs_a[j], "greater")}

    return {"rank": rank, "es": es, "hyper": hyper, "asym": asym}


# ---------------------------------------------------------------- 実行

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--deg", required=True, help="発現解析結果 (CSV/TSV)")
    ap.add_argument("--deg-gene", default="gene", help="遺伝子シンボルの列名")
    ap.add_argument("--deg-lfc", default="log2FoldChange", help="logFC の列名")
    ap.add_argument("--deg-expr", default=None,
                    help="発現量の列名（baseMean / logCPM など）。与えると発現量で層別した置換検定を行う")
    ap.add_argument("--expr-bins", type=int, default=20,
                    help="発現量の層の数（模擬データで 10 層は偽陽性率 6%%、20 層で 5%%）")
    ap.add_argument("--impc", required=True, help="IMPC viability レポート")
    ap.add_argument("--impc-gene", default="gene_symbol")
    ap.add_argument("--impc-call", default="viability")
    ap.add_argument("--subviable", choices=SUBVIABLE_MODES, default="exclude",
                    help="subviable の扱い: exclude=除外(既定) / other=非lethal側 / lethal=lethal側")
    ap.add_argument("--ortho", help="ヒト→マウス オルソログ表 (任意)")
    ap.add_argument("--ortho-human", default="human_symbol")
    ap.add_argument("--ortho-mouse", default="mouse_symbol")
    ap.add_argument("--ortho-mode", choices=("one2one", "all"), default="one2one",
                    help="one2one=1:1 オルソログのみ(既定) / all=全対応")
    ap.add_argument("--cutoffs", default="0.1,0.25,0.5,1.0",
                    help="感度解析に使う |logFC| カットオフ")
    ap.add_argument("--n-perm", type=int, default=10000)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out", default="results")
    args = ap.parse_args()

    os.makedirs(args.out, exist_ok=True)

    # --- DEG 読み込み
    deg = read_table(args.deg)
    cols = {"gene": args.deg_gene, "lfc": args.deg_lfc}
    if args.deg_expr:
        cols["expr"] = args.deg_expr
    for col in cols.values():
        if col not in deg.columns:
            sys.exit(f"DEG ファイルに列 '{col}' がありません。実際の列: {list(deg.columns)}")

    deg = deg[list(cols.values())].dropna()
    deg.columns = list(cols.keys())
    deg["gene"] = normalize_symbol(deg["gene"])

    if args.ortho:
        deg = map_to_mouse(deg, "gene", args.ortho,
                           args.ortho_human, args.ortho_mouse, args.ortho_mode)

    deg = collapse_duplicates(deg)

    # --- IMPC 読み込みと背景集合の確定
    sets = load_lethal_sets(args.impc, args.impc_gene, args.impc_call)

    # 背景は「判定済み ∩ 検出済み」
    bg_all = deg[deg["gene"].isin(sets["tested"])].copy()
    n_sub = int(bg_all["gene"].isin(sets["subviable"]).sum())
    bg, is_lethal = assign_labels(bg_all, sets, args.subviable)

    if len(bg) < 100:
        print(f"[警告] 背景集合が {len(bg)} 遺伝子しかありません。"
              "シンボルの表記系（マウス/ヒト、Ensembl ID など）が一致しているか確認してください。",
              file=sys.stderr)
    if is_lethal.sum() < 3 or (~is_lethal).sum() < 3:
        sys.exit(f"lethal {int(is_lethal.sum())} / 非 lethal {int((~is_lethal).sum())} 遺伝子では"
                 "検定できません（各 3 以上必要）。")

    lfc = bg["lfc"].to_numpy(dtype=float)
    cutoffs = [float(c) for c in args.cutoffs.split(",")]

    if args.deg_expr:
        strata = make_strata(bg["expr"].to_numpy(dtype=float), args.expr_bins)
        strat_label = f"発現量 {len(np.unique(strata))} 層で層別"
        # 発現量の偏りそのものを表示する（交絡の大きさの目安）
        ex = bg["expr"].to_numpy(dtype=float)
        _, p_expr = stats.mannwhitneyu(ex[is_lethal], ex[~is_lethal])
        expr_note = (f"  発現量の中央値: lethal {np.median(ex[is_lethal]):.3g} / "
                     f"その他 {np.median(ex[~is_lethal]):.3g} (MW p={p_expr:.2e})")
    else:
        strata = np.zeros(len(bg), dtype=int)
        strat_label = "層別なし"
        expr_note = None
        print("[警告] --deg-expr が無いので発現量の交絡を制御できません。"
              "致死遺伝子は高発現に偏るため、結果が発現量の効果である可能性が残ります。",
              file=sys.stderr)

    print(f"背景集合: {len(bg)} 遺伝子 "
          f"(IMPC 判定済み {len(sets['tested'])} との積集合; subviable={args.subviable})")
    print(f"  うち lethal ラベル: {int(is_lethal.sum())}")
    print(f"  subviable（背景中）: {n_sub}")
    if expr_note:
        print(expr_note)

    # --- 主結果: 順位検定（解析解）と置換検定
    rank = rank_test(lfc, is_lethal)
    perm = permutation_tests(lfc, is_lethal, strata, cutoffs,
                             n_perm=args.n_perm, seed=args.seed)

    # --- subviable の扱いの感度解析（順位検定のみ）
    sub_rows = []
    for mode in SUBVIABLE_MODES:
        b, y = assign_labels(bg_all, sets, mode)
        r = rank_test(b["lfc"].to_numpy(dtype=float), y)
        if r:
            sub_rows.append({"subviable_mode": mode, **r})

    # --- 感度解析: 複数カットオフでの hypergeometric
    hyper_rows = []
    for c in cutoffs:
        for d in ("down", "up"):
            r = hypergeom_at(lfc, is_lethal, c, d)
            if r:
                hyper_rows.append({**r, **perm["hyper"][(c, d)]})

    # --- 非対称性
    asym_rows = []
    for c in cutoffs:
        r = asymmetry_test(lfc, is_lethal, c)
        if r:
            asym_rows.append({**r, **perm["asym"].get(c, {})})

    # --- 出力
    bg.assign(is_lethal=is_lethal, stratum=strata).to_csv(
        os.path.join(args.out, "background_annotated.csv"), index=False)
    main_rows = []
    if rank:
        main_rows.append({"test": "mann_whitney_analytic", "stratification": "none", **rank})
    main_rows.append({"test": "rank_biserial_permutation", "stratification": strat_label,
                      **perm["rank"]})
    main_rows.append({"test": "running_ES_permutation", "stratification": strat_label,
                      **perm["es"]})
    pd.DataFrame(main_rows).to_csv(os.path.join(args.out, "rank_tests.csv"), index=False)
    pd.DataFrame(sub_rows).to_csv(
        os.path.join(args.out, "subviable_sensitivity.csv"), index=False)
    pd.DataFrame(hyper_rows).to_csv(
        os.path.join(args.out, "hypergeometric_sensitivity.csv"), index=False)
    pd.DataFrame(asym_rows).to_csv(
        os.path.join(args.out, "down_vs_up_asymmetry.csv"), index=False)

    print("\n=== 主結果: 順位検定（カットオフなし） ===")
    if rank:
        print("  [Mann-Whitney 解析解・層別なし]")
        for k, v in rank.items():
            print(f"    {k}: {v}")
    pr, pe = perm["rank"], perm["es"]
    print(f"  [置換検定 {args.n_perm} 回・{strat_label}]")
    print(f"    rank-biserial r: {pr['rank_biserial_r']:.4f}  null 平均: {pr['null_mean_r']:.4f}  "
          f"差: {pr['r_minus_null']:.4f}")
    print(f"    p(lethal が下側): {pr['p_lethal_down']:.5f}  p(上側): {pr['p_lethal_up']:.5f}")
    print(f"    running ES: {pe['ES']:.4f}  null 平均: {pe['null_mean_ES']:.4f}  "
          f"p_down: {pe['p_down']:.5f}  p_up: {pe['p_up']:.5f}")

    print("\n=== subviable の扱いによる感度（Mann-Whitney） ===")
    if sub_rows:
        print(pd.DataFrame(sub_rows)[["subviable_mode", "n_lethal", "n_other",
                                      "rank_biserial_r", "p_lethal_down",
                                      "p_lethal_up"]].to_string(index=False))

    print("\n=== 感度解析: hypergeometric（p_perm は層別置換） ===")
    if hyper_rows:
        print(pd.DataFrame(hyper_rows).to_string(index=False))

    print("\n=== down vs up の非対称性（Fisher / p_perm は層別置換） ===")
    if asym_rows:
        print(pd.DataFrame(asym_rows).to_string(index=False))

    print(f"\n出力先: {args.out}/")


if __name__ == "__main__":
    main()
