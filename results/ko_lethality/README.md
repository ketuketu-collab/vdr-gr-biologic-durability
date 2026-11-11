# 標的遺伝子のマウス KO 致死性

対象: `results/remap_scores_expanded.csv` の 381 遺伝子（ヒト）。
再現: `scripts/ko_lethality_annotation.py`（入力の取得元・判定基準はスクリプト冒頭）。

- `ko_lethality_per_gene.csv`: 1 遺伝子 1 行。オルソログが複数ある場合は最も重い判定の行。
- `ko_lethality_per_ortholog.csv`: マウスオルソログごとの行。

## 情報源と判定

- **IMPC**（viability.csv, latest release）: 統一プロトコルによる離乳前生存判定。
- **MGI**（MGI_GenePheno / MGI_PhenotypicAllele / MP オントロジー）: 文献由来の表現型。
  全アレルが "Null/knockout" 型の homozygous / hemizygous を KO とみなす。
  条件付き KO・hypomorph・点変異は除く。
- **致死の定義:** MP の preweaning lethality（prenatal / perinatal を含む）、
  lethality complete penetrance、lethality at weaning。
  離乳後の早期死亡（premature death）は `mgi_ko_premature_death` 列に別記し、
  致死には含めない。例: Ctla4 KO（3–4 週で死亡）は viable + premature death。

| 判定 | 遺伝子数 | 基準 |
|---|---|---|
| lethal | 102 | IMPC lethal、または MGI KO で完全浸透の致死 |
| lethal_partial | 35 | IMPC subviable、または MGI KO で不完全浸透・浸透度不明の致死のみ |
| viable | 233 | KO データあり、致死報告なし |
| no_ko_data | 9 | FDPS, TYRP1, IGHE, TYMS, EEF2, KEAP1, CLEC4C, KCNK10, HSPA8 |
| no_mouse_ortholog | 2 | CXCL8, ERVW-1 |

KEAP1 と EEF2 は、MGI で KO アレルが hypomorph に分類されているため KO から除外された。

## 臨床ステータス別

| status | lethal | lethal_partial | viable |
|---|---|---|---|
| approved | 4 (JAK1 周産期, S1PR1 胎生, PDGFRB 胎生/周産期, VEGFA 胎生/周産期) | 0 | 25 |
| failed | 2 (MAPK14 胎生, IFNG ※) | 2 (TNFRSF9, CSF2) | 6 |
| ongoing | 1 (RIPK1 周産期) | 0 | 6 |
| steroid_new | 7 | 6 | 34 |
| other | 88 | 27 | 162 |

※ IFNG は IMPC（KOMP アレル）では lethal だが、文献の Ifng KO（16 遺伝子型）では
致死の報告がない。IMPC と文献が食い違う遺伝子は 16 個あり、`impc_vs_literature` 列に記した。

## 探索的な観察（事前登録外）

- KO 致死（lethal + lethal_partial）の遺伝子は、GR スコアが高い
  （中央値 49.5 vs 22.0, Mann-Whitney p = 2e-4）。VDR スコアの差は有意でない（p = 0.11）。
- VDR_dominant 遺伝子のうち KO 致死は 25%（20/80）。それ以外では 40%（117/290）
  （Fisher OR = 0.49, p = 0.013）。
- 注意: ChIP-seq のピーク数は、発現が広く高い遺伝子ほど多くなりやすい。致死遺伝子は
  そのような遺伝子に偏るため、この関連は交絡で説明できる可能性がある。
