# VitD で下がった遺伝子は KO 致死遺伝子に濃縮しているか（GSE189984）

再現: `scripts/run_lethal_gse189984.py`（入力の取得元はスクリプト冒頭に記載）。
検定本体: `scripts/lethal_enrichment.py`。

## データと処理

- GSE189984: ヒト PBMC、10 nM 1,25(OH)2D3 vs 0.1% EtOH、4/8/24/48 h、各 n=3。
  GEO 補足の DESeq2 正規化カウント（12,304 遺伝子にフィルタ済み）から生カウントを
  復元した（全サンプルで整数誤差 < 1e-9）。
- pyDESeq2 0.5.4、時間点ごとに 3 vs 3（対応なし）。主結果は apeglm 縮小 LFC、
  感度解析は MLE LFC。既報値の再現: NOD2 8 h +1.78、TLR10 24 h +1.20、
  CYP24A1 48 h +10.49。
- IMPC viability.csv（latest release; 8,698 遺伝子、lethal / subviable / viable）。
  シンボルは MGI ID で MGI の現行シンボルに更新。
- オルソログ: MGI HOM_MouseHumanSequence.rpt、ヒト Entrez ID で結合、1:1 のみ。
- 背景集合: 約 5,046 遺伝子（うち lethal 1,777）。subviable（527）は除外。
- 発現量交絡の制御: baseMean の順位で 20 層、層内ラベル置換 10,000 回。
  背景中の lethal 遺伝子は高発現（baseMean 中央値 242 vs 145, MW p = 4e-36）。

## 結果（`summary_GSE189984.csv`）

| 時間 | rank-biserial r (MLE) | p(lethal 下側) 層別 | p(lethal 上側) 層別 | down 側 \|LFC\|≥0.5 の lethal 濃縮（層別 fold） |
|---|---|---|---|---|
| 4 h  | +0.033 | 0.96 | 0.039  | 0.62（n=85） |
| 8 h  | +0.011 | 0.98 | 0.022  | 0.47（n=280） |
| 24 h | +0.031 | 0.99 | 0.010  | 0.79（n=1,016） |
| 48 h | +0.049 | 1.00 | 0.0015 | 0.47（n=443） |

- **どの時間点でも、下がった遺伝子が致死遺伝子に濃縮している証拠はない。**
  カットオフなしの順位検定、running ES、全カットオフの hypergeometric、
  down vs up の非対称性のいずれでも、下側の p は 0.2–1.0。
- むしろ逆方向。致死遺伝子の logFC はわずかに上側へずれ（24 h・48 h で p ≤ 0.01）、
  強く下がった遺伝子群では致死遺伝子が期待より少ない（fold 0.3–0.8）。
  ただし効果量は r = 0.03–0.05 と非常に小さい。
- subviable の扱い（除外 / 非 lethal 側 / lethal 側）、縮小 LFC と MLE LFC の
  どちらでも結論は同じ。

## 注意

- 4 h の縮小 LFC はほぼ全遺伝子が 0 付近に縮むため、running ES（|LFC| 重み）は
  不安定（ES = −0.99）。4 h は MLE の結果を見るのがよい。
- PBMC は混合細胞集団であり、細胞組成の変化と遺伝子発現変化を区別できない。
- 1 データセットのみの結果。GR 側（GSE93735 マウス BMDM、GSE135130 THP-1
  マイクロアレイ）は未解析。
