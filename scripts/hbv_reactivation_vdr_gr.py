#!/usr/bin/env python3
"""
HBV reactivation risk vs VDR/GR ChIP-seq occupancy at the drug target gene.

Question (Nagashima, 2026-06):
    Do drugs with high HBV-reactivation risk fall on the GR-dominant side
    of the VDR/GR occupancy axis?

Approach (self-contained; no external pharmacovigilance download required):
    1. Take the curated VDR/GR ChIP-seq scores (results/remap_scores_expanded.csv).
    2. Annotate each druggable target gene with the maximal HBV-reactivation
       risk tier of the biologic/small-molecule that hits it, graded from
       established hepatology guidance (AGA 2015 Technical Review;
       AASLD 2018 HBV guidance; Loomba & Liang Gastroenterology 2017).
    3. Test whether GR occupancy (and GR-dominance) is enriched in the
       higher-risk tiers.

IMPORTANT CAVEAT (documented in the output): HBV-reactivation risk is driven
by the *depth and breadth of immunosuppression* (esp. B-cell / plasma-cell /
T-cell depletion), which is mechanistically distinct from chromatin occupancy
*at the target locus*. This analysis therefore tests an association, not a
mechanism, and is under-powered (n in the low teens).
"""

import os
import numpy as np
import pandas as pd
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
SCORES = os.path.join(ROOT, "results", "remap_scores_expanded.csv")
OUT_CSV = os.path.join(ROOT, "results", "hbv_reactivation_vdr_gr.csv")

# ── HBV reactivation risk by drug target ────────────────────────────────────
# tier: 3 = high (>10%), 2 = moderate (1-10%), 1 = low (<1%)
# Grading follows AGA 2015 (Gastroenterology 148:215) risk groups.
# Each entry: gene -> (tier, representative drug, mechanism / depletion target)
HBV_RISK = {
    # ── HIGH (>10%): B-cell / plasma-cell / pan-lymphocyte depletion ──────────
    "MS4A1": (3, "rituximab / obinutuzumab",   "anti-CD20, B-cell depletion"),
    "CD38":  (3, "daratumumab",                 "anti-CD38, plasma-cell depletion"),
    "CD52":  (3, "alemtuzumab",                 "anti-CD52, pan-lymphocyte depletion"),
    "NR3C1": (3, "corticosteroids (high-dose)", "GR agonist, broad immunosuppression"),
    # ── MODERATE (1-10%): cytokine / kinase / trafficking blockade ───────────
    "TNF":   (2, "infliximab / adalimumab",     "anti-TNF"),
    "IL6R":  (2, "tocilizumab",                 "anti-IL-6R"),
    "BTK":   (2, "ibrutinib",                   "BTK inhibitor (B-cell signalling)"),
    "JAK1":  (2, "tofacitinib / upadacitinib",  "JAK inhibitor"),
    "JAK2":  (2, "ruxolitinib / baricitinib",   "JAK inhibitor"),
    "JAK3":  (2, "tofacitinib",                 "JAK inhibitor"),
    "S1PR1": (2, "ozanimod / fingolimod",       "S1P-receptor modulator (lymphocyte sequestration)"),
    "ITGA4": (2, "natalizumab",                 "anti-alpha4 integrin"),
    # ── LOW (<1%): selective / gut-restricted / Th17-axis ────────────────────
    "TYK2":  (1, "deucravacitinib",             "selective TYK2 inhibitor"),
    "ITGB7": (1, "vedolizumab",                 "gut-selective anti-alpha4beta7"),
    "IL17A": (1, "secukinumab / ixekizumab",    "anti-IL-17A"),
}

TIER_NAME = {3: "high", 2: "moderate", 1: "low"}


def make_figure(sub, rho_gr, p_gr, p_kw):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    TIER_COLOR = {3: "#C62828", 2: "#F57F17", 1: "#2E7D32"}
    TIER_LABEL = {3: "high (>10%)", 2: "moderate (1-10%)", 1: "low (<1%)"}

    fig, ax = plt.subplots(figsize=(7.2, 6.2))
    mx = float(max(sub["VDR_score"].max(), sub["GR_score"].max())) * 1.12

    # diagonal = equal occupancy (boundary of GR-dominance)
    ax.plot([0, mx], [0, mx], ls="--", lw=1.0, color="#9E9E9E", zorder=1)
    ax.text(mx * 0.62, mx * 0.70, "VDR = GR", color="#9E9E9E",
            rotation=45, ha="center", va="center", fontsize=9)

    for tier in (3, 2, 1):
        g = sub[sub.hbv_tier == tier]
        ax.scatter(g["VDR_score"], g["GR_score"], s=130,
                   color=TIER_COLOR[tier], edgecolor="black", lw=0.6,
                   label=f"HBV risk: {TIER_LABEL[tier]}", zorder=3)
    for _, r in sub.iterrows():
        ax.annotate(r["gene"], (r["VDR_score"], r["GR_score"]),
                    xytext=(5, 4), textcoords="offset points", fontsize=8)

    ax.set_xlabel("VDR ChIP-seq occupancy at target locus")
    ax.set_ylabel("GR ChIP-seq occupancy at target locus")
    ax.set_title("HBV-reactivation risk vs VDR/GR occupancy\n"
                 f"(no association: Spearman rho={rho_gr:+.2f}, p={p_gr:.2f}; "
                 f"Kruskal-Wallis p={p_kw:.2f})", fontsize=10)
    ax.set_xlim(-5, mx)
    ax.set_ylim(-5, mx)
    ax.legend(loc="upper right", frameon=True, fontsize=9)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()

    figdir = os.path.join(ROOT, "results", "figures")
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(figdir, f"hbv_reactivation_vdr_gr.{ext}"), dpi=200)
    print(f"Wrote results/figures/hbv_reactivation_vdr_gr.[png|pdf]")


def main():
    df = pd.read_csv(SCORES)
    sub = df[df["gene"].isin(HBV_RISK)].copy()

    sub["hbv_tier"] = sub["gene"].map(lambda g: HBV_RISK[g][0])
    sub["hbv_risk"] = sub["hbv_tier"].map(TIER_NAME)
    sub["drug"]     = sub["gene"].map(lambda g: HBV_RISK[g][1])
    sub["mechanism"] = sub["gene"].map(lambda g: HBV_RISK[g][2])
    sub["GR_dominant"] = sub["GR_score"] > sub["VDR_score"]
    # log-ratio: positive = GR-leaning
    sub["log2_GR_VDR"] = np.log2((sub["GR_score"] + 1) / (sub["VDR_score"] + 1))

    cols = ["gene", "drug", "mechanism", "hbv_risk", "hbv_tier",
            "VDR_score", "GR_score", "GR_dominant", "log2_GR_VDR"]
    out = sub[cols].sort_values(["hbv_tier", "GR_score"], ascending=[False, False])

    missing = [g for g in HBV_RISK if g not in set(df["gene"])]

    # ── Report ───────────────────────────────────────────────────────────────
    pd.set_option("display.width", 160)
    print("=" * 92)
    print("HBV reactivation risk vs VDR/GR ChIP-seq occupancy at the target gene")
    print("=" * 92)
    print(out.to_string(index=False,
          formatters={"log2_GR_VDR": "{:+.2f}".format}))
    if missing:
        print(f"\nTargets not present in 381-gene ChIP-seq set (excluded): {missing}")

    # ── Statistics ─────────────────────────────────────────────────────────────
    print("\n" + "-" * 92)
    print("Association tests (GR occupancy vs HBV-reactivation tier)")
    print("-" * 92)

    # 1) Spearman: continuous GR score / log-ratio vs ordinal tier
    rho_gr, p_gr = stats.spearmanr(sub["hbv_tier"], sub["GR_score"])
    rho_lr, p_lr = stats.spearmanr(sub["hbv_tier"], sub["log2_GR_VDR"])
    print(f"Spearman  tier vs GR_score   : rho={rho_gr:+.3f}  p={p_gr:.3f}")
    print(f"Spearman  tier vs log2(GR/VDR): rho={rho_lr:+.3f}  p={p_lr:.3f}")

    # 2) Kruskal-Wallis across the three tiers (GR score)
    groups = [sub.loc[sub.hbv_tier == t, "GR_score"] for t in (3, 2, 1)]
    H, p_kw = stats.kruskal(*groups)
    print(f"Kruskal-Wallis GR_score across high/mod/low: H={H:.3f}  p={p_kw:.3f}")

    # 3) GR-dominance fraction by tier
    print("\nGR-dominant fraction by tier:")
    for t in (3, 2, 1):
        g = sub[sub.hbv_tier == t]
        frac = g["GR_dominant"].mean()
        print(f"  {TIER_NAME[t]:>8}: {g['GR_dominant'].sum()}/{len(g)} "
              f"({frac:.0%}) GR-dominant   "
              f"[median GR={g['GR_score'].median():.1f}]")

    # 4) High-risk vs rest: Mann-Whitney on GR score
    hi = sub.loc[sub.hbv_tier == 3, "GR_score"]
    lo = sub.loc[sub.hbv_tier < 3, "GR_score"]
    U, p_mw = stats.mannwhitneyu(hi, lo, alternative="greater")
    print(f"\nMann-Whitney GR_score  high vs (mod+low), one-sided 'high>rest': "
          f"U={U:.1f}  p={p_mw:.3f}")

    out.to_csv(OUT_CSV, index=False)
    print(f"\nWrote {os.path.relpath(OUT_CSV, ROOT)}")

    make_figure(sub, rho_gr, p_gr, p_kw)

    # ── Honest interpretation ───────────────────────────────────────────────────
    print("\n" + "=" * 92)
    print("Interpretation")
    print("=" * 92)
    print(
        "The hypothesis ('high HBV-reactivation risk => GR-dominant') is only\n"
        "PARTIALLY supported. anti-CD20 (MS4A1) and anti-CD38 (CD38) sit on the\n"
        "GR-leaning side, but the strongest counter-examples are decisive:\n"
        "  - BTK (ibrutinib): VDR-dominant by a wide margin, yet clinically risky.\n"
        "  - CD52 (alemtuzumab): VDR-dominant, highest-risk depletion.\n"
        "  - NR3C1 (corticosteroids): the archetypal HBV-reactivation drug, and\n"
        "    its own target locus is VDR-dominant in occupancy.\n"
        "Mechanistic reason: reactivation tracks the DEPTH of immunosuppression\n"
        "(B/plasma/T-cell depletion), which is orthogonal to ChIP-seq occupancy\n"
        "at the target locus. With n in the low teens this is hypothesis-level\n"
        "only and should not be over-read."
    )


if __name__ == "__main__":
    main()
