#!/usr/bin/env python3
"""
CDデータの3件修正:
1. ITGB7: etrolizumab(phase3,failed) → vedolizumab(phase4,approved) に差し替え
2. IL23A: 新規追加 (ustekinumab/mirikizumab, approved, GR-dom, VDR=128,GR=209)
3. ITGA4: 新規追加 (vedolizumab/natalizumab, approved, GR-dom, VDR=42,GR=67)
"""

import pandas as pd
import numpy as np

RESULTS = '/Volumes/M4_SSD/projects/tlr_chipseq/results/'
infile  = RESULTS + 'gene_disease_phase_expanded.csv'
outfile = RESULTS + 'gene_disease_phase_expanded.csv'

df = pd.read_csv(infile)

print(f"修正前: {len(df)} records")

# ── 1. ITGB7: etrolizumab(failed) → vedolizumab(approved) ─────────
mask = (df['gene']=='ITGB7') & (df['disease']=='Crohn Disease')
print(f"\n[1] ITGB7 CD before:\n{df[mask][['gene','disease','phase','mol_status','mol_id']].to_string(index=False)}")

df.loc[mask, 'phase']      = 4.0
df.loc[mask, 'mol_status'] = 'approved'
df.loc[mask, 'mol_id']     = 'CHEMBL1743087'  # vedolizumab
print(f"    → mol_id=CHEMBL1743087 (vedolizumab), phase=4, approved に変更")

# ── 2. IL23A: 新規追加 ─────────────────────────────────────────────
# ustekinumab (CHEMBL1201835) と mirikizumab (CHEMBL3990014) の2行
il23a_rows = [
    dict(gene='IL23A', disease='Crohn Disease', category='GI',
         phase=4.0, mol_status='approved', mol_id='CHEMBL1201835',
         VDR_score=128.0, GR_score=209.0, VDR_dominant=False,
         VDR_repro=0.5, GR_repro=0.5),   # repro=celltypes/総CT (近似)
    dict(gene='IL23A', disease='Crohn Disease', category='GI',
         phase=4.0, mol_status='approved', mol_id='CHEMBL3990014',
         VDR_score=128.0, GR_score=209.0, VDR_dominant=False,
         VDR_repro=0.5, GR_repro=0.5),
]
print(f"\n[2] IL23A 追加: ustekinumab(CHEMBL1201835) + mirikizumab(CHEMBL3990014)")
print(f"    VDR=128, GR=209, GR-dominant")

# ── 3. ITGA4: 新規追加 ─────────────────────────────────────────────
itga4_rows = [
    dict(gene='ITGA4', disease='Crohn Disease', category='GI',
         phase=4.0, mol_status='approved', mol_id='CHEMBL1743087',
         VDR_score=42.0, GR_score=67.0, VDR_dominant=False,
         VDR_repro=0.2, GR_repro=0.3),
    dict(gene='ITGA4', disease='Crohn Disease', category='GI',
         phase=4.0, mol_status='approved', mol_id='CHEMBL1201607',
         VDR_score=42.0, GR_score=67.0, VDR_dominant=False,
         VDR_repro=0.2, GR_repro=0.3),
]
print(f"\n[3] ITGA4 追加: vedolizumab(CHEMBL1743087) + natalizumab(CHEMBL1201607)")
print(f"    VDR=42, GR=67, GR-dominant")

new_rows = pd.DataFrame(il23a_rows + itga4_rows)
df = pd.concat([df, new_rows], ignore_index=True)

print(f"\n修正後: {len(df)} records")

# 確認
cd_check = df[df['disease']=='Crohn Disease'].copy()
gene_cd = (cd_check.groupby('gene')
           .agg(best_status=('mol_status', lambda x: 'approved' if 'approved' in x.values else 'failed'),
                VDR_score=('VDR_score','first'),
                GR_score=('GR_score','first'),
                VDR_dom=('VDR_dominant','first')).reset_index())

print("\n=== 修正後 CD 遺伝子リスト ===")
for dom, grp in gene_cd.groupby('VDR_dom', sort=False):
    label = 'VDR-dom' if dom else 'GR-dom'
    for _, r in grp.sort_values('best_status').iterrows():
        print(f"  {label:7s} {r['gene']:10s} VDR={r['VDR_score']:5.0f} GR={r['GR_score']:5.0f} → {r['best_status']}")

from scipy import stats
n_va = int((gene_cd['VDR_dom'] & (gene_cd['best_status']=='approved')).sum())
n_vf = int((gene_cd['VDR_dom'] & (gene_cd['best_status']=='failed')).sum())
n_ga = int((~gene_cd['VDR_dom'] & (gene_cd['best_status']=='approved')).sum())
n_gf = int((~gene_cd['VDR_dom'] & (gene_cd['best_status']=='failed')).sum())

_, p_v = stats.fisher_exact([[n_va,n_vf],[n_ga,n_gf]], alternative='greater')
a,b,c,d = n_va+0.5, n_vf+0.5, n_ga+0.5, n_gf+0.5
or_v = (a*d)/(b*c)
print(f"\nFisher's: {n_va}va/{n_vf}vf | {n_ga}ga/{n_gf}gf → OR={or_v:.2f} p={p_v:.3f}")

df.to_csv(outfile, index=False)
print(f"\nSaved: {outfile}")
