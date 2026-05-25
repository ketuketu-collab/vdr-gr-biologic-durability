#!/usr/bin/env python3
"""CNS autoimmune VDR/GR ChIP-seq figure"""
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

df = pd.read_csv("/Volumes/M4_SSD/projects/tlr_chipseq/results/cns_autoimmune_vdr_gr.csv")

# Drop autoantigens (not drug targets) and positive controls for main fig
exclude_cat = {"NMOSD-Ag","MOGAD-Ag","Myelin-Ag","VDR-ctrl"}
df_drug = df[~df["Category"].isin(exclude_cat)].copy()
df_drug = df_drug[(df_drug["VDR_score"] > 0) | (df_drug["GR_score"] > 0)]
df_drug = df_drug.sort_values("VDR_score", ascending=True)

fig, ax = plt.subplots(figsize=(10, 11))
y = np.arange(len(df_drug))
ax.barh(y - 0.2, df_drug["VDR_score"], height=0.4, color="#1f77b4", label="VDR (THP-1)")
ax.barh(y + 0.2, df_drug["GR_score"],  height=0.4, color="#d62728", label="GR (THP-1)")

# annotate drug names
for i, (_, row) in enumerate(df_drug.iterrows()):
    drug = str(row["Drug"]).split(" — ")[0].split(" (")[0]
    ax.text(max(row["VDR_score"], row["GR_score"]) + 3, i, drug,
            va="center", fontsize=8, color="#444")

ax.set_yticks(y)
ax.set_yticklabels(df_drug["Gene"], fontsize=9)
ax.set_xlabel("ChIP-seq score (MACS2, ReMap2022, THP-1)")
ax.set_title("CNS autoimmune (MS / NMOSD) drug targets — VDR vs GR\n"
             "TSS ± 10 kb, hg38, ReMap2022", fontsize=11)
ax.legend(loc="lower right")
ax.grid(axis="x", alpha=0.3)
ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)

plt.tight_layout()
out_png = "/Volumes/M4_SSD/projects/tlr_chipseq/results/cns_autoimmune_vdr_gr.png"
out_pdf = "/Volumes/M4_SSD/projects/tlr_chipseq/results/figures/cns_autoimmune_vdr_gr.pdf"
plt.savefig(out_png, dpi=150, bbox_inches="tight")
plt.savefig(out_pdf, bbox_inches="tight")
print(f"Saved: {out_png}")
print(f"Saved: {out_pdf}")

# === Pattern pie ===
fig2, ax2 = plt.subplots(figsize=(6, 6))
pat_counts = df_drug["Pattern"].value_counts()
colors = {"VDR only": "#1f77b4", "Both": "#9467bd", "GR only": "#d62728", "Neither": "#999"}
ax2.pie(pat_counts.values, labels=pat_counts.index,
        colors=[colors.get(p, "#aaa") for p in pat_counts.index],
        autopct=lambda p: f"{p:.0f}%\n(n={int(round(p*pat_counts.sum()/100))})",
        startangle=90, wedgeprops=dict(edgecolor="white", linewidth=2))
ax2.set_title(f"CNS autoimmune drug targets (n={pat_counts.sum()})\n"
              "VDR/GR pattern distribution", fontsize=11)
out2 = "/Volumes/M4_SSD/projects/tlr_chipseq/results/cns_autoimmune_pie.png"
plt.savefig(out2, dpi=150, bbox_inches="tight")
print(f"Saved: {out2}")
