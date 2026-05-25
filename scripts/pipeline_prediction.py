#!/usr/bin/env python3
"""
Pipeline drug prediction: VDR/GR rule applied to Phase II/III biologics.

Scoring: same as remap_scores_expanded.csv (TSS ± 5 kb, n_celltypes*10 + n_experiments)
Classifier: VDR_score > GR_score  →  VDR-dominant  →  predicted success
            VDR_score ≤ GR_score  →  GR-dominant   →  predicted failure

Validation on existing data:
  VDR-dominant approval rate: 89.2%  (33/37)
  GR-dominant  approval rate: 50.9%  (28/55)
  Fisher OR = 7.96, p = 0.0001

Outputs:
  results/pipeline_prediction.csv
  results/figures/pipeline_prediction.png  (300 dpi, Arial)
  results/figures/pipeline_prediction.pdf
"""
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats
from sklearn.metrics import roc_auc_score

warnings.filterwarnings("ignore")

# ── Paths ─────────────────────────────────────────────────────────────────────
RES_DIR  = Path("/Volumes/M4_SSD/projects/tlr_chipseq/results")
FIG_DIR  = RES_DIR / "figures"
FIG_DIR.mkdir(parents=True, exist_ok=True)

EXISTING_CSV = RES_DIR / "multi_disease_drug_list.csv"
REMAP_CSV    = RES_DIR / "remap_scores_expanded.csv"

# ── Pipeline drugs (Phase II/III, 2024-2026) ──────────────────────────────────
PIPELINE = [
    # disease, gene, drug, phase_label
    ("AD/Asthma",    "TSLP",      "tezepelumab",   "approved_recent"),
    ("AD/Asthma",    "IL33",      "itepekimab",    "approved_recent"),
    ("AD",           "TNFSF4",    "amlitelimab",   "phase3"),
    ("AD",           "IL31RA",    "nemolizumab",   "phase3"),
    ("AD",           "IL36R",     "spesolimab",    "phase3_AD"),
    ("IBD",          "MADCAM1",   "tulisokibart",  "phase2"),
    ("IBD",          "IL36R",     "spesolimab",    "phase2_CD"),
    ("RA",           "CSF2RA",    "mavrilimumab",  "phase3"),
    ("RA",           "CSF2",      "otilimab",      "phase3"),
    ("RA",           "TNFRSF13C", "ianalumab",     "phase3"),
    ("SLE",          "IFNAR1",    "anifrolumab",   "approved_recent"),
    ("SLE",          "CD40LG",    "dapirolizumab", "phase3"),
    ("SLE",          "CD40",      "iscalimab",     "phase3"),
    ("SLE",          "TNFRSF13B", "telitacicept",  "phase3"),
    ("Psoriasis/RA", "IL17F",     "bimekizumab",   "approved_recent"),
]

# Manually computed scores for genes absent from remap_scores_expanded.csv
# (TSS±5kb, n_celltypes*10 + n_experiments, same formula)
MANUAL_SCORES = {
    "IL36R":     (0.0,  12.0),
    "CSF2RA":    (13.0, 25.0),
    "TNFRSF13B": (0.0,  22.0),
}

# ── Step 1: load existing ReMAP scores ────────────────────────────────────────
remap = pd.read_csv(REMAP_CSV)[["gene", "VDR_score", "GR_score"]]
gene_scores = {}
for _, row in remap.iterrows():
    gene_scores[row["gene"]] = (row["VDR_score"], row["GR_score"])
for gene, (v, g) in MANUAL_SCORES.items():
    gene_scores[gene] = (v, g)

# ── Step 2: build pipeline dataframe ──────────────────────────────────────────
rows = []
for disease, gene, drug, phase in PIPELINE:
    v, g = gene_scores.get(gene, (0.0, 0.0))
    vdr_dom = v > g
    if phase == "approved_recent":
        prediction = "approved_recent"
    else:
        prediction = "predicted_success" if vdr_dom else "predicted_failure"
    rows.append({
        "disease":    disease,
        "gene":       gene,
        "drug":       drug,
        "phase":      phase,
        "vdr_score":  v,
        "gr_score":   g,
        "ratio":      v / (g + 1),
        "vdr_class":  "VDR" if vdr_dom else "GR",
        "prediction": prediction,
    })

pipe_df = pd.DataFrame(rows)

# ── Step 3: statistics from existing data ─────────────────────────────────────
existing = pd.read_csv(EXISTING_CSV)
sub      = existing[existing["status"].isin(["approved", "failed"])].copy()
sub["vdr_dom"] = sub["vdr_score"] > sub["gr_score"]

ct = pd.crosstab(sub["vdr_dom"], sub["status"])
vdr_app  = ct.loc[True, "approved"]
vdr_fail = ct.loc[True, "failed"]
gr_app   = ct.loc[False, "approved"]
gr_fail  = ct.loc[False, "failed"]

odds_ratio, p_fisher = stats.fisher_exact(
    [[vdr_app, vdr_fail], [gr_app, gr_fail]]
)
vdr_rate = vdr_app / (vdr_app + vdr_fail)
gr_rate  = gr_app  / (gr_app  + gr_fail)

# Mann-Whitney
app_v  = sub.loc[sub["status"] == "approved", "vdr_score"]
fail_v = sub.loc[sub["status"] == "failed",   "vdr_score"]
_, p_mw = stats.mannwhitneyu(app_v, fail_v, alternative="greater")

# AUC from VDR-dominant binary rule
sub["pred_bin"] = sub["vdr_dom"].astype(int)
sub["label"]    = (sub["status"] == "approved").astype(int)
auc_binary = roc_auc_score(sub["label"], sub["pred_bin"])

print("=== Existing data statistics ===")
print(f"  VDR-dominant approval rate: {vdr_rate:.1%}  ({vdr_app}/{vdr_app+vdr_fail})")
print(f"  GR-dominant  approval rate: {gr_rate:.1%}  ({gr_app}/{gr_app+gr_fail})")
print(f"  Fisher OR = {odds_ratio:.2f},  p = {p_fisher:.4f}")
print(f"  Mann-Whitney (VDR approved > failed): p = {p_mw:.4f}")
print(f"  Binary AUC (VDR-dominant rule): {auc_binary:.3f}")

# ── Print pipeline summary ─────────────────────────────────────────────────────
print("\n=== Pipeline Drug Predictions ===")
print(f"{'Drug':26s} {'Gene':12s} {'VDR':>7s} {'GR':>7s} {'Class':>5s}  {'Prediction'}")
print("-" * 78)
for _, row in pipe_df.iterrows():
    print(f"  {row['drug']:24s} {row['gene']:12s} {row['vdr_score']:7.1f} "
          f"{row['gr_score']:7.1f} {row['vdr_class']:>5s}  {row['prediction']}")

phase_only = pipe_df[~pipe_df["phase"].str.startswith("approved")]
n_pipeline = len(phase_only)
n_success  = (phase_only["prediction"] == "predicted_success").sum()
n_failure  = (phase_only["prediction"] == "predicted_failure").sum()
n_rec_app  = (pipe_df["prediction"] == "approved_recent").sum()

print(f"\n--- Pipeline summary ---")
print(f"  Recently approved (in list):      {n_rec_app}")
print(f"  Phase II/III total (unique drugs): {n_pipeline}")
print(f"  Predicted SUCCESS (VDR > GR):     {n_success}  ({100*n_success/max(n_pipeline,1):.1f}%)")
print(f"  Predicted FAILURE (GR ≥ VDR):     {n_failure}  ({100*n_failure/max(n_pipeline,1):.1f}%)")

# ── Save CSV ───────────────────────────────────────────────────────────────────
out_csv = RES_DIR / "pipeline_prediction.csv"
pipe_df.to_csv(out_csv, index=False)
print(f"\nSaved: {out_csv}")

# ── Step 4: Figure ─────────────────────────────────────────────────────────────
print("\nGenerating figure …")

plt.rcParams.update({
    "font.family": "Arial",
    "font.size":    11,
    "axes.linewidth": 1.2,
    "pdf.fonttype":  42,
    "ps.fonttype":   42,
})

fig, ax = plt.subplots(figsize=(10.5, 8.5))

# ─── Existing drugs ────────────────────────────────────────────────────────────
app  = existing[existing["status"] == "approved"]
fail = existing[existing["status"] == "failed"]

ax.scatter(
    np.log1p(app["vdr_score"]),  np.log1p(app["gr_score"]),
    color="#E67E22", s=72, alpha=0.70, zorder=3,
    label=f"Approved biologic (existing, n={len(app)})"
)
ax.scatter(
    np.log1p(fail["vdr_score"]), np.log1p(fail["gr_score"]),
    color="#2980B9", s=72, alpha=0.70, zorder=3,
    label=f"Failed biologic (existing, n={len(fail)})"
)

# ─── Pipeline drugs ────────────────────────────────────────────────────────────
COLOR_MAP  = {"approved_recent":   "#E67E22",
              "predicted_success": "#27AE60",
              "predicted_failure": "#C0392B"}
MARKER_MAP = {"approved_recent":   "o",
              "predicted_success": "*",
              "predicted_failure": "*"}
SIZE_MAP   = {"approved_recent":   150,
              "predicted_success": 360,
              "predicted_failure": 360}
EDGE_MAP   = {"approved_recent":   "#922B21",
              "predicted_success": "#1E8449",
              "predicted_failure": "#7B241C"}
LABEL_MAP  = {"approved_recent":   "Pipeline: approved (2023-2024)",
              "predicted_success": "Pipeline: predicted SUCCESS ★",
              "predicted_failure": "Pipeline: predicted FAILURE ★"}

seen_labels = set()

for _, row in pipe_df.iterrows():
    pred = row["prediction"]
    xv = float(np.log1p(row["vdr_score"]))
    yv = float(np.log1p(row["gr_score"]))

    leg = LABEL_MAP[pred] if LABEL_MAP[pred] not in seen_labels else "_nolegend_"
    seen_labels.add(LABEL_MAP[pred])

    ax.scatter(xv, yv,
               color=COLOR_MAP[pred], marker=MARKER_MAP[pred],
               s=SIZE_MAP[pred], edgecolors=EDGE_MAP[pred], linewidths=1.0,
               zorder=6, label=leg)

    # Drug name annotation
    ax.annotate(
        row["drug"],
        (xv, yv),
        xytext=(5, 4),
        textcoords="offset points",
        fontsize=7.0, color="#222222", zorder=7
    )

# ─── Diagonal y = x  (VDR = GR boundary) ─────────────────────────────────────
xlim = ax.get_xlim()
ylim = ax.get_ylim()
mn = min(xlim[0], ylim[0]) - 0.3
mx = max(xlim[1], ylim[1]) + 0.3
ax.plot([mn, mx], [mn, mx],
        color="#95A5A6", ls="-", lw=1.4, zorder=2, label="VDR = GR (decision boundary)")
ax.set_xlim(mn, mx)
ax.set_ylim(mn, mx)

# ─── Shaded region annotations ────────────────────────────────────────────────
y_bot = ax.get_ylim()[0] + 0.10
x_bot = ax.get_xlim()[0] + 0.08

ax.text(ax.get_xlim()[1] - 0.1, y_bot,
        "VDR-dominant\n→ success zone",
        fontsize=8.5, color="#1E8449", va="bottom",
        ha="right", style="italic", alpha=0.85)
ax.text(x_bot, ax.get_ylim()[1] - 0.1,
        "GR-dominant\n→ failure zone",
        fontsize=8.5, color="#922B21", va="top",
        ha="left", style="italic", alpha=0.85)

# ─── Labels & title ───────────────────────────────────────────────────────────
ax.set_xlabel("VDR peak score  [log(1 + score)]", fontsize=13, fontweight="bold")
ax.set_ylabel("GR (NR3C1) peak score  [log(1 + score)]", fontsize=13, fontweight="bold")
ax.set_title(
    "Predictive validation of pipeline biologics using VDR/GR ChIP-seq rule\n"
    f"(ReMAP2022, TSS ± 5 kb  |  existing data: OR = {odds_ratio:.1f}, "
    f"Fisher p = {p_fisher:.4g}  |  VDR-dom approval rate {vdr_rate:.0%} vs GR-dom {gr_rate:.0%})",
    fontsize=11, fontweight="bold", pad=11
)

ax.legend(loc="upper left", fontsize=8.5, framealpha=0.92,
          edgecolor="#CCCCCC", handlelength=1.5)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.tick_params(direction="out", length=4)

plt.tight_layout(pad=1.5)

out_png = FIG_DIR / "pipeline_prediction.png"
out_pdf = FIG_DIR / "pipeline_prediction.pdf"
fig.savefig(out_png, dpi=300, bbox_inches="tight")
fig.savefig(out_pdf, dpi=300, bbox_inches="tight")
print(f"Saved: {out_png}")
print(f"Saved: {out_pdf}")
plt.close(fig)

print("\nDone.")
