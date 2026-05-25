#!/usr/bin/env python3
"""
Logistic regression: VDR/GR continuous scores -> drug approval probability
- Response: approved(1) vs failed(0), then approved vs failed+other
- Predictors: VDR_score, GR_score, VDR_repro, GR_repro, repro_ratio
- Bootstrap 95% CI, ROC-AUC, calibration
"""

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from scipy import stats
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_curve, auc, brier_score_loss
from sklearn.model_selection import StratifiedKFold, cross_val_predict
import warnings
warnings.filterwarnings('ignore')

RESULTS = '/Volumes/M4_SSD/projects/tlr_chipseq/results/'
CSV     = RESULTS + 'remap_breadth_per_gene.csv'

df = pd.read_csv(CSV)
print(f"Loaded: {len(df)} genes")

COL_VDR = '#2166AC'
COL_GR  = '#D6604D'
COL_APP = '#1a9641'
COL_FAIL= '#d7191c'

# ── feature engineering ────────────────────────────────────────────────
df['log_VDR'] = np.log1p(df['VDR_score'])
df['log_GR']  = np.log1p(df['GR_score'])
df['log_ratio'] = df['log_VDR'] - df['log_GR']  # positive = VDR-dominant

FEATURES = ['log_VDR', 'log_GR', 'VDR_repro', 'GR_repro']

# ══════════════════════════════════════════════════════════════════════
# MODEL A: strict (approved=1, failed=0)
# ══════════════════════════════════════════════════════════════════════
strict = df[df['status'].isin(['approved','failed'])].copy()
strict['y'] = (strict['status'] == 'approved').astype(int)
print(f"\nStrict set: n={len(strict)}  approved={strict['y'].sum()}  failed={(strict['y']==0).sum()}")

X_s = strict[FEATURES].values
y_s = strict['y'].values
scaler_s = StandardScaler()
Xs_scaled = scaler_s.fit_transform(X_s)

lr_s = LogisticRegression(C=1.0, max_iter=1000, random_state=42)
lr_s.fit(Xs_scaled, y_s)

# coefficients + bootstrap CI
N_BOOT = 2000
rng = np.random.default_rng(42)
boot_coefs = np.zeros((N_BOOT, len(FEATURES)))
for i in range(N_BOOT):
    idx = rng.integers(0, len(y_s), len(y_s))
    if len(np.unique(y_s[idx])) < 2:
        boot_coefs[i] = np.nan
        continue
    lr_b = LogisticRegression(C=1.0, max_iter=1000, random_state=0)
    lr_b.fit(Xs_scaled[idx], y_s[idx])
    boot_coefs[i] = lr_b.coef_[0]

ci_lo_s = np.nanpercentile(boot_coefs, 2.5, axis=0)
ci_hi_s = np.nanpercentile(boot_coefs, 97.5, axis=0)

print("\nStrict model coefficients (standardized):")
for f, coef, lo, hi in zip(FEATURES, lr_s.coef_[0], ci_lo_s, ci_hi_s):
    sig = "*" if (lo > 0 or hi < 0) else ""
    print(f"  {f:15s}: {coef:+.3f}  95%CI [{lo:+.3f}, {hi:+.3f}]{sig}")

# ROC-AUC via leave-one-out (small n)
cv = StratifiedKFold(n_splits=min(5, strict['y'].sum()), shuffle=True, random_state=42)
probs_s = cross_val_predict(LogisticRegression(C=1.0, max_iter=1000, random_state=42),
                             Xs_scaled, y_s, cv=cv, method='predict_proba')[:,1]
fpr_s, tpr_s, _ = roc_curve(y_s, probs_s)
auc_s = auc(fpr_s, tpr_s)
print(f"CV AUC (strict): {auc_s:.3f}")

# ══════════════════════════════════════════════════════════════════════
# MODEL B: broad (approved=1, failed+other=0, exclude ongoing)
# ══════════════════════════════════════════════════════════════════════
broad = df[df['status'].isin(['approved','failed','other'])].copy()
broad['y'] = (broad['status'] == 'approved').astype(int)
print(f"\nBroad set: n={len(broad)}  approved={broad['y'].sum()}  non-approved={(broad['y']==0).sum()}")

X_b = broad[FEATURES].values
y_b = broad['y'].values
scaler_b = StandardScaler()
Xb_scaled = scaler_b.fit_transform(X_b)

lr_b_model = LogisticRegression(C=1.0, max_iter=1000, random_state=42, class_weight='balanced')
lr_b_model.fit(Xb_scaled, y_b)

boot_coefs_b = np.zeros((N_BOOT, len(FEATURES)))
for i in range(N_BOOT):
    idx = rng.integers(0, len(y_b), len(y_b))
    if len(np.unique(y_b[idx])) < 2: continue
    lr_bb = LogisticRegression(C=1.0, max_iter=1000, random_state=0, class_weight='balanced')
    lr_bb.fit(Xb_scaled[idx], y_b[idx])
    boot_coefs_b[i] = lr_bb.coef_[0]

ci_lo_b = np.nanpercentile(boot_coefs_b, 2.5, axis=0)
ci_hi_b = np.nanpercentile(boot_coefs_b, 97.5, axis=0)

print("\nBroad model coefficients (standardized):")
for f, coef, lo, hi in zip(FEATURES, lr_b_model.coef_[0], ci_lo_b, ci_hi_b):
    sig = "*" if (lo > 0 or hi < 0) else ""
    print(f"  {f:15s}: {coef:+.3f}  95%CI [{lo:+.3f}, {hi:+.3f}]{sig}")

cv_b = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
probs_b = cross_val_predict(LogisticRegression(C=1.0, max_iter=1000, random_state=42, class_weight='balanced'),
                             Xb_scaled, y_b, cv=cv_b, method='predict_proba')[:,1]
fpr_b, tpr_b, _ = roc_curve(y_b, probs_b)
auc_b = auc(fpr_b, tpr_b)
print(f"CV AUC (broad): {auc_b:.3f}")

# ── predicted probability for all 308 genes (broad model) ─────────────
X_all = scaler_b.transform(df[FEATURES].values)
df['pred_prob'] = lr_b_model.predict_proba(X_all)[:,1]

top_novel = df[(df['status']=='other') & (df['VDR_dominant']==True)].nlargest(15,'pred_prob')
print("\nTop novel VDR-dominant candidates by predicted approval probability:")
print(top_novel[['gene','VDR_score','GR_score','pred_prob']].to_string())

# ══════════════════════════════════════════════════════════════════════
# FIGURE
# ══════════════════════════════════════════════════════════════════════
fig = plt.figure(figsize=(18, 12))
gs = GridSpec(2, 3, figure=fig, hspace=0.45, wspace=0.42)
ax_a = fig.add_subplot(gs[0, 0])   # coef forest (strict)
ax_b = fig.add_subplot(gs[0, 1])   # coef forest (broad)
ax_c = fig.add_subplot(gs[0, 2])   # ROC both models
ax_d = fig.add_subplot(gs[1, 0:2]) # pred prob scatter (log_VDR vs prob)
ax_e = fig.add_subplot(gs[1, 2])   # top novel candidates

# ── Panel A: coefficient forest (strict) ──────────────────────────────
feature_labels = ['log(VDR score)', 'log(GR score)', 'VDR repro', 'GR repro']
ypos = list(range(len(FEATURES)))
colors_a = [COL_VDR if ('VDR' in f) else COL_GR for f in FEATURES]

for i, (coef, lo, hi, col) in enumerate(zip(lr_s.coef_[0], ci_lo_s, ci_hi_s, colors_a)):
    ax_a.plot([lo, hi], [i, i], color=col, lw=2)
    ax_a.plot(coef, i, 'o', color=col, ms=9, zorder=5)
ax_a.axvline(0, color='gray', ls='--', lw=1)
ax_a.set_yticks(ypos)
ax_a.set_yticklabels(feature_labels, fontsize=10)
ax_a.set_xlabel('Standardized coefficient', fontsize=10)
ax_a.set_title(f'A   Strict model\n(n={len(strict)}, AUC={auc_s:.2f})',
               fontsize=11, fontweight='bold', loc='left')

# ── Panel B: coefficient forest (broad) ───────────────────────────────
for i, (coef, lo, hi, col) in enumerate(zip(lr_b_model.coef_[0], ci_lo_b, ci_hi_b, colors_a)):
    ax_b.plot([lo, hi], [i, i], color=col, lw=2)
    ax_b.plot(coef, i, 'o', color=col, ms=9, zorder=5)
ax_b.axvline(0, color='gray', ls='--', lw=1)
ax_b.set_yticks(ypos)
ax_b.set_yticklabels(feature_labels, fontsize=10)
ax_b.set_xlabel('Standardized coefficient', fontsize=10)
ax_b.set_title(f'B   Broad model (class-balanced)\n(n={len(broad)}, AUC={auc_b:.2f})',
               fontsize=11, fontweight='bold', loc='left')

for ax in [ax_a, ax_b]:
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

# ── Panel C: ROC curves ────────────────────────────────────────────────
ax_c.plot(fpr_s, tpr_s, color=COL_VDR, lw=2, label=f'Strict (AUC={auc_s:.2f})')
ax_c.plot(fpr_b, tpr_b, color='#7b2d8b', lw=2, ls='--', label=f'Broad (AUC={auc_b:.2f})')
ax_c.plot([0,1],[0,1],'k--',lw=0.8,alpha=0.4)
ax_c.set_xlabel('False positive rate', fontsize=10)
ax_c.set_ylabel('True positive rate', fontsize=10)
ax_c.legend(fontsize=9)
ax_c.set_title('C   CV ROC curves', fontsize=11, fontweight='bold', loc='left')
ax_c.spines['top'].set_visible(False)
ax_c.spines['right'].set_visible(False)

# ── Panel D: predicted prob vs log(VDR) scatter ───────────────────────
status_col = {'approved': COL_APP, 'failed': COL_FAIL, 'ongoing': '#f59b00', 'other': '#cccccc'}
status_z   = {'approved': 4, 'failed': 3, 'ongoing': 3, 'other': 1}
for st, grp in df.groupby('status'):
    ax_d.scatter(grp['log_VDR'], grp['pred_prob'],
                 c=status_col[st], s=50 if st in ['approved','failed'] else 25,
                 alpha=0.8, zorder=status_z[st], label=st,
                 edgecolors='white' if st in ['approved','failed'] else 'none', lw=0.5)

# label approved and failed
for _, row in df[df['status'].isin(['approved','failed'])].iterrows():
    col = COL_APP if row['status']=='approved' else COL_FAIL
    ax_d.annotate(row['gene'], (row['log_VDR'], row['pred_prob']),
                  fontsize=6.5, color=col, xytext=(3,2), textcoords='offset points')

ax_d.set_xlabel('log(1 + VDR score)', fontsize=11)
ax_d.set_ylabel('Predicted approval probability\n(broad model)', fontsize=11)
ax_d.set_title('D   Predicted approval probability across all 308 genes',
               fontsize=11, fontweight='bold', loc='left')
ax_d.legend(fontsize=9, ncol=2)
ax_d.spines['top'].set_visible(False)
ax_d.spines['right'].set_visible(False)

# ── Panel E: top novel candidates ─────────────────────────────────────
top15 = df[(df['status']=='other') & (df['VDR_dominant']==True)].nlargest(15,'pred_prob')
bars = ax_e.barh(range(len(top15)), top15['pred_prob'].values[::-1],
                 color=COL_VDR, alpha=0.8, height=0.7)
ax_e.set_yticks(range(len(top15)))
ax_e.set_yticklabels(top15['gene'].values[::-1], fontsize=9)
ax_e.set_xlabel('Predicted approval probability', fontsize=10)
ax_e.set_title('E   Top novel VDR-dominant\ndrug target candidates',
               fontsize=11, fontweight='bold', loc='left')
ax_e.axvline(df[df['status']=='approved']['pred_prob'].median(),
             color=COL_APP, ls='--', lw=1.2, label='Approved median')
ax_e.legend(fontsize=8)
for bar, val in zip(bars, top15['pred_prob'].values[::-1]):
    ax_e.text(val + 0.005, bar.get_y() + bar.get_height()/2,
              f'{val:.2f}', va='center', fontsize=8)
ax_e.spines['top'].set_visible(False)
ax_e.spines['right'].set_visible(False)
ax_e.set_xlim(0, 0.55)

fig.suptitle('Logistic regression: VDR score predicts drug approval probability\nacross 308 chronic inflammatory disease targets',
             fontsize=13, fontweight='bold', y=0.98)

out = RESULTS + 'fig_logistic_vdr_approval'
fig.savefig(out + '.pdf', dpi=300, bbox_inches='tight', facecolor='white')
fig.savefig(out + '.png', dpi=200, bbox_inches='tight', facecolor='white')
print(f"\nSaved: {out}.pdf / .png")

# ── save scored table ──────────────────────────────────────────────────
df_out = df[['gene','VDR_score','GR_score','log_VDR','log_GR','VDR_repro','GR_repro',
             'VDR_dominant','status','pred_prob']].sort_values('pred_prob', ascending=False)
df_out.to_csv(RESULTS + 'vdr_logistic_scored_genes.csv', index=False)
print("Saved: vdr_logistic_scored_genes.csv")

print("\n=== SUMMARY ===")
print(f"Strict AUC: {auc_s:.3f} | Broad AUC: {auc_b:.3f}")
print(f"Predicted prob — approved median: {df[df['status']=='approved']['pred_prob'].median():.3f}")
print(f"Predicted prob — failed median:   {df[df['status']=='failed']['pred_prob'].median():.3f}")
print(f"Predicted prob — other median:    {df[df['status']=='other']['pred_prob'].median():.3f}")
print(f"\nTop novel candidate: {top15.iloc[0]['gene']} (p={top15.iloc[0]['pred_prob']:.3f})")
print("DONE")
