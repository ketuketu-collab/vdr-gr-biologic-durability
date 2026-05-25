#!/usr/bin/env python3
"""
独立検証2: Prospective prediction
現在 Phase2/3 の VDR-dominant 遺伝子を「将来承認候補」としてリストアップ
→ 論文投稿後に結果が出ることで前向き検証になる
"""

import pandas as pd
import numpy as np

RESULTS = '/Volumes/M4_SSD/projects/tlr_chipseq/results/'

df = pd.read_csv(RESULTS + 'gene_disease_phase_expanded.csv')
remap = pd.read_csv(RESULTS + 'remap_scores_expanded.csv')

# ongoing (phase2/3, not yet approved/failed) を取り出す
# gene_disease_phase_expanded には ongoing は除外済みなので
# steroid_disease_indications から再構築
indic = pd.read_csv(RESULTS + 'steroid_disease_indications.csv')

# 承認済み・失敗確定の mol_id
approved_mols = set(df[df['mol_status']=='approved']['mol_id'])
failed_mols   = set(df[df['mol_status']=='failed']['mol_id'])

# phase2以上で承認でも失敗でもない → ongoing候補
ongoing = indic[
    (indic['max_phase'] >= 2) &
    (~indic['mol_id'].isin(approved_mols)) &
    (~indic['mol_id'].isin(failed_mols))
].copy()
print(f"Ongoing phase2+ records: {len(ongoing)}, unique mols: {ongoing['mol_id'].nunique()}")

# gene×disease テーブルから VDR-dom遺伝子の phase 情報取得
# expand_steroid_targets の全mol(ongoing含む)を再読み込みするのが理想だが
# ここでは remap_scores の VDR-dom遺伝子に対して indic のphaseを参照

# VDR-dominant 遺伝子リスト
vdr_dom_genes = remap[remap['VDR_dominant']==True]['gene'].tolist()

# gene_disease_phase_expanded の全レコード（approved+failedのみ）で最高phaseを確認
gene_max = (df.groupby('gene')
            .agg(max_phase=('phase','max'),
                 best_status=('mol_status', lambda x: 'approved' if 'approved' in x.values else 'failed'),
                 VDR_score=('VDR_score','first'),
                 GR_score=('GR_score','first'),
                 VDR_dom=('VDR_dominant','first'),
                 diseases=('disease', lambda x: ', '.join(sorted(x.unique()))))
            .reset_index())

# VDR-dom で phase2/3 止まり（承認なし、失敗確定でもない）= promising candidates
# ここではphase3到達済みをhighlight
candidates = gene_max[
    (gene_max['VDR_dom']==True) &
    (gene_max['best_status']=='failed') &  # failedだがphase3まで行った
    (gene_max['max_phase'] >= 3)
].sort_values('VDR_score', ascending=False)

print(f"\n=== VDR-dominant Phase3到達遺伝子（未承認） ===")
print("→ 適応症・薬剤変更で再挑戦価値が高い候補")
print(candidates[['gene','VDR_score','GR_score','max_phase','diseases']].to_string(index=False))

# VDR-dom で承認済みの遺伝子（参考）
approved_vdr = gene_max[
    (gene_max['VDR_dom']==True) &
    (gene_max['best_status']=='approved')
].sort_values('VDR_score', ascending=False)

print(f"\n=== VDR-dominant 承認済み遺伝子（参考） ===")
print(approved_vdr[['gene','VDR_score','GR_score','diseases']].to_string(index=False))

# VDR scoreが高いがまだ試されていない（未試験）遺伝子
tested_genes = set(df['gene'].unique())
untested_vdr = remap[
    (remap['VDR_dominant']==True) &
    (~remap['gene'].isin(tested_genes))
].sort_values('VDR_score', ascending=False).head(20)

print(f"\n=== VDR-dominant 未試験遺伝子 Top20（将来の標的候補） ===")
print(untested_vdr[['gene','VDR_score','GR_score']].to_string(index=False))

# サマリー CSV 保存
summary = pd.concat([
    candidates.assign(prediction_type='Phase3_retry_candidate'),
    untested_vdr.assign(prediction_type='Novel_target_candidate',
                        max_phase=0, best_status='untested',
                        diseases='').rename(columns={'VDR_dominant':'VDR_dom'})
], ignore_index=True)

cols = ['gene','VDR_score','GR_score','VDR_dom','max_phase','best_status','diseases','prediction_type']
summary[[c for c in cols if c in summary.columns]].to_csv(
    RESULTS + 'prospective_predictions.csv', index=False)
print(f"\nSaved: {RESULTS}prospective_predictions.csv")
print("DONE")
