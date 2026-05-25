#!/usr/bin/env python3
"""
Logistic regression on all 308 genes using Excel drug phase data
- Model A: 承認済>=1 (202) vs 未承認 (106)  — balanced labels
- Model B: P1→II通過率 as continuous outcome (linear regression)
- Model C: P1失敗のみ(21) vs 承認済(202) — strict failure signal
"""

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from scipy import stats
from sklearn.linear_model import LogisticRegression, LinearRegression
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_curve, auc
from sklearn.model_selection import StratifiedKFold, KFold, cross_val_predict
import warnings
warnings.filterwarnings('ignore')

RESULTS = '/Volumes/M4_SSD/projects/tlr_chipseq/results/'

# ── load & merge ───────────────────────────────────────────────────────
df_xl  = pd.read_excel(RESULTS + 'steroid_309genes_analysis.xlsx',
                       sheet_name='全体(308遺伝子)VDR-GR降順')
df_csv = pd.read_csv(RESULTS + 'remap_breadth_per_gene.csv')
df = df_xl.merge(df_csv[['gene','VDR_repro','GR_repro','VDR_dominant']],
                 left_on='遺伝子', right_on='gene', how='left')
df = df.rename(columns={'VDR':'VDR_score','GR':'GR_score','遺伝子':'gene_name'})
df['log_VDR']   = np.log1p(df['VDR_score'])
df['log_GR']    = np.log1p(df['GR_score'])
df['log_ratio'] = df['log_VDR'] - df['log_GR']
print(f"Dataset: {len(df)} genes")

COL_VDR = '#2166AC'
COL_GR  = '#D6604D'
COL_APP = '#1a9641'
COL_FAIL= '#d7191c'

FEATURES      = ['log_VDR', 'log_GR', 'VDR_repro', 'GR_repro']
FEATURE_LABELS= ['log(VDR score)', 'log(GR score)', 'VDR repro', 'GR repro']
N_BOOT = 2000
rng = np.random.default_rng(42)

def bootstrap_coefs(X_sc, y, is_binary=True, n_boot=N_BOOT):
    coefs = np.zeros((n_boot, X_sc.shape[1]))
    for i in range(n_boot):
        idx = rng.integers(0, len(y), len(y))
        if is_binary and len(np.unique(y[idx])) < 2:
            coefs[i] = np.nan; continue
        if is_binary:
            m = LogisticRegression(C=1.0, max_iter=1000, random_state=0)
        else:
            m = LinearRegression()
        m.fit(X_sc[idx], y[idx])
        coefs[i] = m.coef_[0] if is_binary else m.coef_
    return coefs

# ══════════════════════════════════════════════════════════════════════
# MODEL A: 承認済>=1 vs 未承認 (n=308)
# ══════════════════════════════════════════════════════════════════════
df['y_approved'] = (df['承認済'] >= 1).astype(int)
print(f"\nModel A: approved={df['y_approved'].sum()}  not_approved={(df['y_approved']==0).sum()}")

X   = df[FEATURES].values
y_a = df['y_approved'].values
sc_a = StandardScaler(); Xa = sc_a.fit_transform(X)

lr_a = LogisticRegression(C=1.0, max_iter=1000, random_state=42)
lr_a.fit(Xa, y_a)

boot_a = bootstrap_coefs(Xa, y_a)
ci_lo_a = np.nanpercentile(boot_a, 2.5, axis=0)
ci_hi_a = np.nanpercentile(boot_a, 97.5, axis=0)

print("Model A coefficients:")
for f, coef, lo, hi in zip(FEATURE_LABELS, lr_a.coef_[0], ci_lo_a, ci_hi_a):
    sig = "***" if (lo>0 or hi<0) else ""
    print(f"  {f:20s}: {coef:+.3f}  95%CI [{lo:+.3f}, {hi:+.3f}]{sig}")

cv_a = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
prob_a = cross_val_predict(LogisticRegression(C=1.0, max_iter=1000, random_state=42),
                           Xa, y_a, cv=cv_a, method='predict_proba')[:,1]
fpr_a, tpr_a, _ = roc_curve(y_a, prob_a)
auc_a = auc(fpr_a, tpr_a)
print(f"CV AUC: {auc_a:.3f}")

# ══════════════════════════════════════════════════════════════════════
# MODEL B: P1→II通過率 as continuous outcome (linear regression)
# ══════════════════════════════════════════════════════════════════════
y_b = df['P1→II通過率'].values.astype(float)
sc_b = StandardScaler(); Xb = sc_b.fit_transform(X)

lr_b = LinearRegression()
lr_b.fit(Xb, y_b)

boot_b = bootstrap_coefs(Xb, y_b, is_binary=False)
ci_lo_b = np.nanpercentile(boot_b, 2.5, axis=0)
ci_hi_b = np.nanpercentile(boot_b, 97.5, axis=0)

# Spearman between VDR score and pass rate
r_vdr, p_vdr = stats.spearmanr(df['log_VDR'], y_b)
r_gr,  p_gr  = stats.spearmanr(df['log_GR'],  y_b)
print(f"\nModel B (P1→II通過率):")
print(f"  Spearman VDR: r={r_vdr:.3f}  p={p_vdr:.3e}")
print(f"  Spearman GR:  r={r_gr:.3f}  p={p_gr:.3e}")
print("  Linear regression coefficients:")
for f, coef, lo, hi in zip(FEATURE_LABELS, lr_b.coef_, ci_lo_b, ci_hi_b):
    sig = "***" if (lo>0 or hi<0) else ""
    print(f"  {f:20s}: {coef:+.3f}  95%CI [{lo:+.3f}, {hi:+.3f}]{sig}")

cv_b = KFold(n_splits=5, shuffle=True, random_state=42)
pred_b = cross_val_predict(LinearRegression(), Xb, y_b, cv=cv_b)
r_cv, p_cv = stats.pearsonr(y_b, pred_b)
print(f"  CV r={r_cv:.3f}  p={p_cv:.3e}")

# ══════════════════════════════════════════════════════════════════════
# MODEL C: strict — P1失敗のみ(0) vs 承認済(1)
# ══════════════════════════════════════════════════════════════════════
mask_c = ((df['承認済'] >= 1) | ((df['承認済']==0) & (df['P1止まり'] >= 1)))
df_c   = df[mask_c].copy()
df_c['y_strict'] = (df_c['承認済'] >= 1).astype(int)
print(f"\nModel C (strict): approved={df_c['y_strict'].sum()}  p1_fail={(df_c['y_strict']==0).sum()}")

Xc_raw = df_c[FEATURES].values
y_c    = df_c['y_strict'].values
sc_c   = StandardScaler(); Xc = sc_c.fit_transform(Xc_raw)

lr_c = LogisticRegression(C=1.0, max_iter=1000, random_state=42)
lr_c.fit(Xc, y_c)

boot_c = bootstrap_coefs(Xc, y_c)
ci_lo_c = np.nanpercentile(boot_c, 2.5, axis=0)
ci_hi_c = np.nanpercentile(boot_c, 97.5, axis=0)

print("Model C coefficients:")
for f, coef, lo, hi in zip(FEATURE_LABELS, lr_c.coef_[0], ci_lo_c, ci_hi_c):
    sig = "***" if (lo>0 or hi<0) else ""
    print(f"  {f:20s}: {coef:+.3f}  95%CI [{lo:+.3f}, {hi:+.3f}]{sig}")

cv_c = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
prob_c = cross_val_predict(LogisticRegression(C=1.0, max_iter=1000, random_state=42),
                           Xc, y_c, cv=cv_c, method='predict_proba')[:,1]
fpr_c, tpr_c, _ = roc_curve(y_c, prob_c)
auc_c = auc(fpr_c, tpr_c)
print(f"CV AUC: {auc_c:.3f}")

# ── predicted probability for all 308 genes ───────────────────────────
df['pred_prob_A'] = lr_a.predict_proba(Xa)[:,1]
df['pred_passrate'] = lr_b.predict(Xb)

top_novel = (df[df['承認済']==0]
             .sort_values('pred_prob_A', ascending=False)
             .head(15))
print("\nTop 15 novel candidates (未承認, high predicted prob):")
print(top_novel[['gene_name','VDR_score','GR_score','P1止まり','VDR_dominant','pred_prob_A']].to_string())

# ══════════════════════════════════════════════════════════════════════
# FIGURE
# ══════════════════════════════════════════════════════════════════════
fig = plt.figure(figsize=(18, 14))
gs = GridSpec(2, 3, figure=fig, hspace=0.45, wspace=0.42)
ax_a = fig.add_subplot(gs[0, 0])
ax_b = fig.add_subplot(gs[0, 1])
ax_c = fig.add_subplot(gs[0, 2])
ax_d = fig.add_subplot(gs[1, 0:2])
ax_e = fig.add_subplot(gs[1, 2])

colors_feat = [COL_VDR, COL_GR, COL_VDR, COL_GR]

# ── Panel A: Model A coefficient forest ───────────────────────────────
for i,(coef,lo,hi,col) in enumerate(zip(lr_a.coef_[0], ci_lo_a, ci_hi_a, colors_feat)):
    ax_a.plot([lo,hi],[i,i], color=col, lw=2.5)
    ax_a.plot(coef, i, 'o', color=col, ms=10, zorder=5)
ax_a.axvline(0, color='gray', ls='--', lw=1)
ax_a.set_yticks(range(4)); ax_a.set_yticklabels(FEATURE_LABELS, fontsize=10)
ax_a.set_xlabel('Standardized coefficient', fontsize=10)
ax_a.set_title(f'A   Logistic: 承認済≥1 vs 未承認\n(n=308, AUC={auc_a:.2f})',
               fontsize=11, fontweight='bold', loc='left')

# ── Panel B: Model B linear forest ────────────────────────────────────
for i,(coef,lo,hi,col) in enumerate(zip(lr_b.coef_, ci_lo_b, ci_hi_b, colors_feat)):
    ax_b.plot([lo,hi],[i,i], color=col, lw=2.5)
    ax_b.plot(coef, i, 's', color=col, ms=10, zorder=5)
ax_b.axvline(0, color='gray', ls='--', lw=1)
ax_b.set_yticks(range(4)); ax_b.set_yticklabels(FEATURE_LABELS, fontsize=10)
ax_b.set_xlabel('Standardized coefficient', fontsize=10)
ax_b.set_title(f'B   Linear: P1→II通過率 (continuous)\n(n=308, CV r={r_cv:.2f})',
               fontsize=11, fontweight='bold', loc='left')

# ── Panel C: ROC curves ────────────────────────────────────────────────
ax_c.plot(fpr_a, tpr_a, color=COL_VDR, lw=2,
          label=f'Model A: 承認≥1 (AUC={auc_a:.2f})')
ax_c.plot(fpr_c, tpr_c, color=COL_GR,  lw=2, ls='--',
          label=f'Model C: strict (AUC={auc_c:.2f})')
ax_c.plot([0,1],[0,1],'k--',lw=0.8,alpha=0.4)
ax_c.set_xlabel('False positive rate', fontsize=10)
ax_c.set_ylabel('True positive rate', fontsize=10)
ax_c.legend(fontsize=9)
ax_c.set_title('C   CV ROC curves', fontsize=11, fontweight='bold', loc='left')

# ── Panel D: scatter log_VDR vs pred_prob_A, colored by approval ──────
approved_mask = df['y_approved'] == 1
ax_d.scatter(df.loc[~approved_mask,'log_VDR'], df.loc[~approved_mask,'pred_prob_A'],
             c='#cccccc', s=30, alpha=0.6, label='未承認', zorder=1)
ax_d.scatter(df.loc[approved_mask,'log_VDR'],  df.loc[approved_mask,'pred_prob_A'],
             c=COL_APP, s=50, alpha=0.8, label='承認済≥1', zorder=3,
             edgecolors='white', lw=0.5)

# highlight top novel
for _, row in top_novel.head(8).iterrows():
    ax_d.annotate(row['gene_name'],
                  (row['log_VDR'], row['pred_prob_A']),
                  fontsize=7, color='#333333',
                  xytext=(4,2), textcoords='offset points')

# smoothed trend line
from scipy.ndimage import uniform_filter1d
sorted_idx = np.argsort(df['log_VDR'].values)
xs_sorted  = df['log_VDR'].values[sorted_idx]
ys_sorted  = df['pred_prob_A'].values[sorted_idx]
smoothed   = uniform_filter1d(ys_sorted, size=20)
ax_d.plot(xs_sorted, smoothed, color=COL_VDR, lw=2, ls='-', alpha=0.7, label='Trend')

ax_d.set_xlabel('log(1 + VDR score)', fontsize=11)
ax_d.set_ylabel('Predicted approval probability', fontsize=11)
ax_d.set_title('D   Predicted approval probability vs VDR score (n=308)',
               fontsize=11, fontweight='bold', loc='left')
ax_d.legend(fontsize=9)

# ── Panel E: top novel candidates bar ─────────────────────────────────
top15 = top_novel.head(15)
bars = ax_e.barh(range(len(top15)), top15['pred_prob_A'].values[::-1],
                 color=[COL_VDR if v else COL_GR for v in top15['VDR_dominant'].values[::-1]],
                 alpha=0.85, height=0.7)
ax_e.set_yticks(range(len(top15)))
ax_e.set_yticklabels(top15['gene_name'].values[::-1], fontsize=9)
ax_e.axvline(df[df['y_approved']==1]['pred_prob_A'].median(),
             color=COL_APP, ls='--', lw=1.2, label='Approved median')
ax_e.set_xlabel('Predicted approval probability', fontsize=10)
ax_e.set_title('E   Top novel candidates\n(未承認遺伝子)', fontsize=11, fontweight='bold', loc='left')
ax_e.legend(fontsize=8)
for bar, val, dom in zip(bars, top15['pred_prob_A'].values[::-1],
                          top15['VDR_dominant'].values[::-1]):
    label = f"{val:.2f}" + (" ▲" if dom else "")
    ax_e.text(val+0.005, bar.get_y()+bar.get_height()/2,
              label, va='center', fontsize=8)
ax_e.set_xlim(0, 1.0)

for ax in [ax_a,ax_b,ax_c,ax_d,ax_e]:
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

fig.suptitle('VDR score predicts drug approval across all 308 chronic inflammatory targets\n'
             '(▲ = VDR-dominant gene)',
             fontsize=13, fontweight='bold', y=0.99)

out = RESULTS + 'fig_logistic_308genes'
fig.savefig(out + '.pdf', dpi=300, bbox_inches='tight', facecolor='white')
fig.savefig(out + '.png', dpi=200, bbox_inches='tight', facecolor='white')
print(f"\nSaved: {out}.pdf / .png")

# ── save scored table ──────────────────────────────────────────────────
out_df = df[['gene_name','VDR_score','GR_score','VDR_dominant',
             '承認済','P1止まり','P1→II通過率',
             'pred_prob_A','pred_passrate']].sort_values('pred_prob_A', ascending=False)
out_df.to_csv(RESULTS + 'logistic_308genes_scored.csv', index=False)
print("Saved: logistic_308genes_scored.csv")

print("\n=== SUMMARY ===")
print(f"Model A (承認≥1 vs 未承認, n=308):    AUC={auc_a:.3f}")
print(f"Model B (P1→II通過率, linear, n=308): CV r={r_cv:.3f}  p={p_cv:.2e}")
print(f"Model C (strict approved vs P1fail):  AUC={auc_c:.3f}")
print(f"\nSpearman VDR vs P1→II通過率: r={r_vdr:.3f}  p={p_vdr:.2e}")
print(f"Spearman GR  vs P1→II通過率: r={r_gr:.3f}  p={p_gr:.2e}")
print(f"\nApproved median pred_prob: {df[df['y_approved']==1]['pred_prob_A'].median():.3f}")
print(f"Not-approved median:       {df[df['y_approved']==0]['pred_prob_A'].median():.3f}")
print("\nDONE")
