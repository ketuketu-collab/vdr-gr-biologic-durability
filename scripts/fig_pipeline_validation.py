"""
Predictive validation of pipeline biologics using VDR/GR ChIP-seq rule

Two-panel figure:
  a) 2D scatter: VDR score vs GR score for all agents
     - approved reference drugs (diamonds)
     - pipeline agents (circles = predicted success, × = predicted failure)
     - tulisokibart (★ pre-registered prediction, TNFSF15)
  b) Ranked horizontal bar: VDR/GR ratio per pipeline agent
     - colour-coded by prediction
     - threshold line at ratio=1.0 (VDR=GR)
"""

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
import warnings
warnings.filterwarnings("ignore")

# ── Colours ───────────────────────────────────────────────────────────────────
C_SUCCESS  = "#1565C0"   # predicted durable
C_FAILURE  = "#B71C1C"   # predicted failure
C_APPROVED = "#2E7D32"   # recently approved (reference)
C_STAR     = "#E65100"   # tulisokibart
C_GREY     = "#616161"

BASE = "/Volumes/M4_SSD/projects/tlr_chipseq/results"

# ── Reference approved drugs (for scatter background) ─────────────────────────
reference = pd.DataFrame([
    dict(drug="anti-IL-23\n(ustekinumab/risa)", gene="IL23A", vdr=211.05, gr=70.34,  group="ref_approved"),
    dict(drug="vedolizumab",    gene="ITGB7",  vdr=62.27,  gr=17.23,  group="ref_approved"),
    dict(drug="anti-TNF",       gene="TNF",    vdr=36.86,  gr=19.49,  group="ref_approved"),
    dict(drug="tocilizumab",    gene="IL6R",   vdr=99.59,  gr=82.90,  group="ref_approved"),
    dict(drug="JAK inhibitors\n(BBW ⚠)",gene="JAK1", vdr=12.81, gr=85.48, group="ref_bbw"),
])

# ── Pipeline agents (fixed CSV) ───────────────────────────────────────────────
pipe = pd.read_csv(f"{BASE}/pipeline_prediction.csv")

# Friendly labels for plot
NAME_MAP = {
    "tulisokibart":   "tulisokibart\n(TNFSF15)",
    "amlitelimab":    "amlitelimab\n(OX40L/TNFSF4)",
    "nemolizumab":    "nemolizumab\n(IL31RA)",
    "spesolimab":     "spesolimab\n(IL36R)",
    "mavrilimumab":   "mavrilimumab\n(CSF2RA)",
    "otilimab":       "otilimab\n(CSF2)",
    "ianalumab":      "ianalumab\n(BAFF-R)",
    "dapirolizumab":  "dapirolizumab\n(CD40LG)",
    "iscalimab":      "iscalimab\n(CD40)",
    "telitacicept":   "telitacicept\n(BAFF)",
    "tezepelumab":    "tezepelumab\n(TSLP) ✓approved",
    "itepekimab":     "itepekimab\n(IL33) ✓approved",
    "anifrolumab":    "anifrolumab\n(IFNAR1) ✓approved",
    "bimekizumab":    "bimekizumab\n(IL17F) ✓approved",
}

pipe["label"] = pipe["drug"].map(NAME_MAP).fillna(pipe["drug"])
pipe["ratio_safe"] = pipe["vdr_score"] / (pipe["gr_score"].replace(0, 0.1))

# ── Figure layout ─────────────────────────────────────────────────────────────
fig, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(18, 8),
                                  gridspec_kw=dict(wspace=0.35))
fig.suptitle("Predictive validation of pipeline biologics using VDR/GR ChIP-seq rule",
             fontsize=14, fontweight="bold", y=0.98)

# ════════════════════════════════════════════════════════════════════════════
# Panel a — 2D scatter: VDR vs GR score
# ════════════════════════════════════════════════════════════════════════════

# Shaded zones
lim = 320
xs = np.linspace(0, lim, 300)
ax_a.fill_between(xs, 0, xs, alpha=0.05, color=C_SUCCESS, zorder=0)   # VDR > GR
ax_a.fill_between(xs, xs, lim, alpha=0.05, color=C_FAILURE, zorder=0) # GR > VDR
ax_a.plot(xs, xs, "--", color="grey", lw=0.9, alpha=0.5, zorder=1)
ax_a.text(220, 200, "VDR = GR", fontsize=8, color="grey", rotation=38, alpha=0.6)
ax_a.text(180, 30,  "VDR-dominant\n(predicted durable)",
          fontsize=9, color=C_SUCCESS, alpha=0.7, ha="center")
ax_a.text(60, 240,  "GR-dominant\n(predicted failure)",
          fontsize=9, color=C_FAILURE, alpha=0.7, ha="center")

# Reference approved drugs
for _, row in reference.iterrows():
    col = C_APPROVED if row["group"] == "ref_approved" else "#C62828"
    ax_a.scatter(row["vdr"], row["gr"], marker="D", s=130, color=col,
                 edgecolors="black", linewidths=0.7, zorder=4, alpha=0.9)
    offset = (6, 4)
    if row["gene"] == "IL23A": offset = (6, -12)
    if row["gene"] == "JAK1":  offset = (6, 4)
    ax_a.annotate(row["drug"], xy=(row["vdr"], row["gr"]),
                  xytext=offset, textcoords="offset points",
                  fontsize=8, color=col, fontweight="bold")

# Pipeline agents
for _, row in pipe.iterrows():
    if row["drug"] == "tulisokibart":
        continue  # Plot separately
    if row["prediction"] == "predicted_failure":
        col, marker, ms, zorder = C_FAILURE, "X", 100, 3
    elif row["prediction"] == "approved_recent":
        col, marker, ms, zorder = C_APPROVED, "o", 70, 3
    else:
        col, marker, ms, zorder = C_SUCCESS, "o", 80, 3

    ax_a.scatter(row["vdr_score"], row["gr_score"],
                 marker=marker, s=ms, color=col,
                 edgecolors="white", linewidths=0.4, alpha=0.8, zorder=zorder)

    # Label high-GR or high-VDR agents
    if row["gr_score"] > 50 or row["vdr_score"] > 20:
        short = row["drug"].replace("\n"," ")
        ax_a.annotate(f"{short}\n({row['gene']})",
                      xy=(row["vdr_score"], row["gr_score"]),
                      xytext=(7, 3), textcoords="offset points",
                      fontsize=7, color=col)

# Tulisokibart — star
ax_a.scatter(12.4, 286.0, marker="*", s=500, color=C_STAR, zorder=6,
             edgecolors="black", linewidths=0.8)
ax_a.annotate(
    "Tulisokibart (TNFSF15)\nVDR=12.4  GR=286.0\nRatio=0.043\n"
    "★ Pre-registered failure\n   osf.io/tnp63 (2026-05-25)\n"
    "   ATLAS-UC results: Nov 2026",
    xy=(12.4, 286.0), xytext=(70, 220),
    fontsize=8.5, color=C_STAR, fontweight="bold",
    arrowprops=dict(arrowstyle="->", color=C_STAR, lw=1.5,
                    connectionstyle="arc3,rad=-0.2"),
    bbox=dict(boxstyle="round,pad=0.4", facecolor="white",
              edgecolor=C_STAR, alpha=0.9))

ax_a.set_xlabel("VDR (NR1I1) ChIP-seq score  [ReMAP2022]", fontsize=12)
ax_a.set_ylabel("GR (NR3C1) ChIP-seq score  [ReMAP2022]", fontsize=12)
ax_a.set_xlim(-15, 260)
ax_a.set_ylim(-10, 320)
ax_a.spines[["top","right"]].set_visible(False)
ax_a.set_title("a", fontsize=13, fontweight="bold", loc="left")

# Legend
leg_a = [
    Line2D([0],[0], marker="D", color="w", markerfacecolor=C_APPROVED,
           markersize=10, label="Approved reference (◎ durable)"),
    Line2D([0],[0], marker="D", color="w", markerfacecolor="#C62828",
           markersize=10, label="Approved (⚠ Black Box warning)"),
    Line2D([0],[0], marker="o", color="w", markerfacecolor=C_SUCCESS,
           markersize=9, label="Pipeline: predicted durable"),
    Line2D([0],[0], marker="X", color="w", markerfacecolor=C_FAILURE,
           markersize=9, label="Pipeline: predicted failure"),
    Line2D([0],[0], marker="o", color="w", markerfacecolor=C_APPROVED,
           markersize=9, label="Pipeline: recently approved"),
    Line2D([0],[0], marker="*", color="w", markerfacecolor=C_STAR,
           markersize=13, label="Tulisokibart (pre-registered)"),
]
ax_a.legend(handles=leg_a, fontsize=8, loc="lower right",
            framealpha=0.9, edgecolor="grey")

# ════════════════════════════════════════════════════════════════════════════
# Panel b — Ranked horizontal bar: VDR/GR ratio
# ════════════════════════════════════════════════════════════════════════════

# Build complete set: pipeline + reference approved
all_agents = []

for _, row in reference.iterrows():
    ratio = row["vdr"] / (row["gr"] + 0.1)
    all_agents.append(dict(
        label=row["drug"],
        gene=row["gene"],
        ratio=ratio,
        vdr=row["vdr"],
        gr=row["gr"],
        category="approved_ref" if row["group"]=="ref_approved" else "approved_bbw",
    ))

for _, row in pipe.iterrows():
    if row["prediction"] in ("approved_recent",):
        cat = "approved_recent"
    elif row["prediction"] == "predicted_failure":
        cat = "failure"
    else:
        cat = "success"
    all_agents.append(dict(
        label=row["label"],
        gene=row["gene"],
        ratio=float(row["ratio_safe"]),
        vdr=row["vdr_score"],
        gr=row["gr_score"],
        category=cat,
    ))

df_bar = pd.DataFrame(all_agents).sort_values("ratio", ascending=True).reset_index(drop=True)

CAT_COL = {
    "approved_ref":    C_APPROVED,
    "approved_bbw":    "#C62828",
    "approved_recent": "#66BB6A",
    "success":         C_SUCCESS,
    "failure":         C_FAILURE,
}
CAT_ALPHA = {
    "approved_ref": 0.90, "approved_bbw": 0.90,
    "approved_recent": 0.75, "success": 0.75, "failure": 0.75,
}

y_pos = np.arange(len(df_bar))

for yi, (_, row) in enumerate(df_bar.iterrows()):
    is_tuli = "tulisokibart" in row["label"].lower()
    col   = C_STAR if is_tuli else CAT_COL[row["category"]]
    alpha = 1.0   if is_tuli else CAT_ALPHA[row["category"]]
    lw    = 2.0   if is_tuli else 0.3

    ax_b.barh(yi, min(row["ratio"], 8), color=col, alpha=alpha,
              edgecolor="black" if is_tuli else "white", linewidth=lw,
              zorder=3 if is_tuli else 2)

    # Score annotation on bar
    ratio_disp = f"{row['ratio']:.2f}" if row["ratio"] < 8 else ">8"
    ax_b.text(min(row["ratio"], 8) + 0.08, yi, ratio_disp,
              va="center", fontsize=7.5,
              color=C_STAR if is_tuli else "black",
              fontweight="bold" if is_tuli else "normal")

    # Drug label
    gene_str = f"  ({row['gene']})"
    label_str = row["label"].split("\n")[0]  # first line only
    ax_b.text(-0.12, yi, label_str + gene_str,
              va="center", ha="right", fontsize=8.0,
              color=C_STAR if is_tuli else "black",
              fontweight="bold" if is_tuli else "normal")

# Threshold line: VDR = GR → ratio = 1.0
ax_b.axvline(1.0, color="grey", lw=1.5, ls="--", alpha=0.7, zorder=1)
ax_b.text(1.05, len(df_bar)-0.5, "VDR = GR\n(ratio=1)", fontsize=8,
          color="grey", va="top")

# Tulisokibart annotation
tuli_yi = df_bar[df_bar["label"].str.contains("tulisokibart", case=False)].index
if len(tuli_yi) > 0:
    tyi = tuli_yi[0]
    ax_b.annotate("← Pre-registered FAILURE\n   ATLAS-UC results: Nov 2026\n   osf.io/tnp63",
                  xy=(0.05, tyi), xytext=(3.5, tyi + 1.5),
                  fontsize=8.5, color=C_STAR, fontweight="bold",
                  arrowprops=dict(arrowstyle="->", color=C_STAR, lw=1.3))

ax_b.set_yticks([])
ax_b.set_xlabel("VDR/GR ChIP-seq ratio  (higher = more VDR-dominant)", fontsize=11)
ax_b.set_xlim(-0.1, 9.0)
ax_b.set_ylim(-0.8, len(df_bar) - 0.2)
ax_b.spines[["top","right","left"]].set_visible(False)
ax_b.set_title("b  VDR/GR ratio ranked — pipeline agents vs approved references",
               fontsize=12, fontweight="bold", loc="left")

# Shading
ax_b.axvspan(0, 1.0, alpha=0.04, color=C_FAILURE)
ax_b.axvspan(1.0, 9.0, alpha=0.04, color=C_SUCCESS)

# Legend
leg_b = [
    mpatches.Patch(facecolor=C_APPROVED, alpha=0.9, label="Approved ◎ durable"),
    mpatches.Patch(facecolor="#C62828",  alpha=0.9, label="Approved ⚠ Black Box"),
    mpatches.Patch(facecolor="#66BB6A",  alpha=0.8, label="Recently approved"),
    mpatches.Patch(facecolor=C_SUCCESS,  alpha=0.75, label="Pipeline: predicted durable"),
    mpatches.Patch(facecolor=C_FAILURE,  alpha=0.75, label="Pipeline: predicted failure"),
    mpatches.Patch(facecolor=C_STAR,     alpha=1.0,  label="Tulisokibart (pre-registered)"),
]
ax_b.legend(handles=leg_b, fontsize=8, loc="lower right",
            framealpha=0.9, edgecolor="grey")

# ── Save ──────────────────────────────────────────────────────────────────────
out = f"{BASE}/figures/fig_pipeline_validation"
fig.savefig(f"{out}.pdf", dpi=300, bbox_inches="tight")
fig.savefig(f"{out}.png", dpi=300, bbox_inches="tight")
print(f"Saved: {out}.pdf / .png")
plt.close()
