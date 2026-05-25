"""
Type 2免疫 × VDR ChIP-seq 解析
3つのアプローチを並列実行:
  A) ReMap2022 epithelial cell types (LS180, prostate epithelial) でType2遺伝子スコア
  B) GSE326394 HaCaT VitD eRNA BigWig → Type2遺伝子プロモーターのシグナル変化
  C) ENCODE GATA3 ChIP-seq → Th2/lymphocyte でのType2遺伝子スコア

出力: /Volumes/M4_SSD/projects/tlr_chipseq/results/type2_vdr_chipseq_results.xlsx
"""
import gzip, re, os, time, json
import urllib.request
import numpy as np
import pandas as pd
import pyBigWig
from pathlib import Path

OUT_DIR = Path("/Volumes/M4_SSD/projects/tlr_chipseq/results")
DATA_DIR = Path("/Volumes/M4_SSD/projects/tlr10_chipseq_vdr_gr/data")
TMP = Path("/tmp/type2_chipseq")
TMP.mkdir(exist_ok=True)

GTF = Path("/Volumes/M4_SSD/ref/gencode.v43.primary_assembly.annotation.gtf.gz")
VDR_BED = DATA_DIR / "remap2022_VDR_all_macs2_hg38.bed.gz"
GR_BED  = DATA_DIR / "remap2022_NR3C1_all_macs2_hg38.bed.gz"
PROMOTER_HW = 5000  # ±5kb

# ── Type 2 target genes ────────────────────────────────────────────────────────
TYPE2_GENES = [
    "TSLP","IL33","IL25","IL4R","IL13","IL5","IL4","IL31RA",
    "IGHE","CCR3","SIGLEC8","IL1RL1","CPA3","HDC","FCER1A",
    "IL17RB","HHIP","AREG","EREG","S100A8","S100A9",
    # 比較用骨髄系
    "IL23A","TNF","IL6R","VEGFA","CD274","IL10","TYK2","ITGB7","C5AR1"
]

print("=== Type 2 VDR ChIP-seq Analysis ===")
print(f"Target genes: {len(TYPE2_GENES)}")

# ── Step 1: GTFから遺伝子座標を取得 ──────────────────────────────────────────
print("\n[Step 1] Extracting gene coordinates from GTF...")

gene_coords = {}
with gzip.open(GTF, 'rt') as f:
    for line in f:
        if line.startswith('#'): continue
        parts = line.split('\t')
        if len(parts) < 9 or parts[2] != 'gene': continue
        info = parts[8]
        m = re.search(r'gene_name "([^"]+)"', info)
        if not m: continue
        gene = m.group(1)
        if gene not in TYPE2_GENES: continue
        chrom, start, end, strand = parts[0], int(parts[3])-1, int(parts[4]), parts[6]
        tss = start if strand == '+' else end
        gene_coords[gene] = {
            'chrom': chrom, 'tss': tss, 'strand': strand,
            'prom_start': max(0, tss - PROMOTER_HW),
            'prom_end': tss + PROMOTER_HW
        }

print(f"  Found coordinates for {len(gene_coords)}/{len(TYPE2_GENES)} genes")

# ── Step 2A: ReMap2022 epithelial cell types スコア ──────────────────────────
print("\n[Step 2A] ReMap2022 epithelial cell scoring...")

EPITHELIAL_CELLS = [
    'LS180_125', 'LS180',
    'primary-prostate-epithelial-cell', 'primary-prostate-epithelial-cell_ethanol',
    'kidney-cortex',
    'LCLGM10861_CALCITRIOL',  # LCL (B cell line)
]

def score_bed_for_genes(bed_gz, cell_filter=None, genes=None):
    """BEDファイルから各遺伝子プロモーターのChIP-seqピーク数をカウント"""
    scores = {g: 0 for g in (genes or TYPE2_GENES)}
    if genes is None:
        genes = TYPE2_GENES
    with gzip.open(bed_gz, 'rt') as f:
        for line in f:
            parts = line.split('\t')
            if len(parts) < 4: continue
            name = parts[3]
            if cell_filter:
                ct = '.'.join(name.split('.')[2:])
                if not any(ct == c for c in cell_filter): continue
            chrom, start, end = parts[0], int(parts[1]), int(parts[2])
            for gene in genes:
                if gene not in gene_coords: continue
                gc = gene_coords[gene]
                if gc['chrom'] != chrom: continue
                if start < gc['prom_end'] and end > gc['prom_start']:
                    scores[gene] += 1
    return scores

# THP-1のみ
print("  Scoring THP-1 VDR...")
thp1_cells = ['THP-1_1d_1-25-OH-2D3','THP-1_2h_1-25-OH-2D3','THP-1_CAL','THP-1_EtOH_1d','THP-1_EtOH_2h','THP-1','THP-1_2H_125D','THP-1_1H_125D']
scores_thp1_vdr = score_bed_for_genes(VDR_BED, cell_filter=thp1_cells)

# 上皮系細胞のみ
print("  Scoring Epithelial VDR (LS180, prostate, kidney)...")
scores_epithelial_vdr = score_bed_for_genes(VDR_BED, cell_filter=EPITHELIAL_CELLS)

# 全細胞合計
print("  Scoring All cells VDR...")
scores_all_vdr = score_bed_for_genes(VDR_BED, cell_filter=None)

# GR全体
print("  Scoring All cells GR...")
scores_all_gr = score_bed_for_genes(GR_BED, cell_filter=None)

print("  Done with ReMap2022 scoring")
print("  VDR scores (THP-1) for Type 2 genes:")
for g in ["TSLP","IL33","IL13","IL4R","IL5","CCR3","IL31RA","IGHE","IL25"]:
    if g in scores_thp1_vdr:
        print(f"    {g}: THP1={scores_thp1_vdr[g]}, Epithelial={scores_epithelial_vdr.get(g,0)}, All={scores_all_vdr.get(g,0)}, GR={scores_all_gr.get(g,0)}")

# ── Step 2B: GSE326394 HaCaT BigWig ──────────────────────────────────────────
print("\n[Step 2B] Downloading HaCaT BigWig from GSE326394...")

bw_files = {
    'HaCaT_VitD_plus':  'ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM9630nnn/GSM9630400/suppl/GSM9630400_HaCaT.VD.plus.bw',
    'HaCaT_VitD_minus': 'ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM9630nnn/GSM9630400/suppl/GSM9630400_HaCaT.VD.minus.bw',
    'HaCaT_ctrl_plus':  'ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM9630nnn/GSM9630399/suppl/GSM9630399_HaCaT.EtOH.plus.bw',
    'HaCaT_ctrl_minus': 'ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM9630nnn/GSM9630399/suppl/GSM9630399_HaCaT.EtOH.minus.bw',
}

def download_file(url, dest):
    if dest.exists():
        print(f"    {dest.name}: already exists ({dest.stat().st_size//1024//1024}MB)")
        return True
    print(f"    Downloading {dest.name}...", end=' ', flush=True)
    try:
        req = urllib.request.Request(url)
        with urllib.request.urlopen(req, timeout=300) as resp:
            with open(dest, 'wb') as f:
                while True:
                    chunk = resp.read(1024*1024)
                    if not chunk: break
                    f.write(chunk)
        print(f"OK ({dest.stat().st_size//1024//1024}MB)")
        return True
    except Exception as e:
        print(f"FAILED: {e}")
        return False

bw_local = {}
for label, url in bw_files.items():
    fname = url.split('/')[-1]
    dest = TMP / fname
    if download_file(url, dest):
        bw_local[label] = dest

def score_bigwig_at_promoters(bw_path, gene_coords, strand_filter=None):
    """BigWigファイルで各遺伝子プロモーターの平均シグナルを計算"""
    scores = {}
    try:
        bw = pyBigWig.open(str(bw_path))
        for gene, gc in gene_coords.items():
            chrom = gc['chrom']
            start = gc['prom_start']
            end = gc['prom_end']
            if strand_filter and gc['strand'] != strand_filter:
                scores[gene] = 0.0
                continue
            try:
                stats = bw.stats(chrom, start, end, type='mean', nBins=1)
                scores[gene] = stats[0] if stats and stats[0] else 0.0
            except:
                scores[gene] = 0.0
        bw.close()
    except Exception as e:
        print(f"  BigWig error: {e}")
        scores = {g: 0.0 for g in gene_coords}
    return scores

bw_scores = {}
for label, path in bw_local.items():
    print(f"  Scoring {label}...")
    strand = '+' if 'plus' in label else '-'
    bw_scores[label] = score_bigwig_at_promoters(path, gene_coords, strand_filter=None)

# VitD誘導シグナル = VitD - control (stranded eRNA)
if len(bw_scores) >= 3:
    print("  Calculating VitD-induced eRNA signal...")
    hackat_vitd_induced = {}
    for gene in TYPE2_GENES:
        if gene not in gene_coords:
            hackat_vitd_induced[gene] = 0.0
            continue
        strand = gene_coords[gene]['strand']
        key_vitd = f'HaCaT_VitD_{"plus" if strand=="+" else "minus"}'
        key_ctrl = f'HaCaT_ctrl_{"plus" if strand=="+" else "minus"}'
        vitd_sig = bw_scores.get(key_vitd, {}).get(gene, 0.0) or 0.0
        ctrl_sig = bw_scores.get(key_ctrl, {}).get(gene, 0.0) or 0.0
        hackat_vitd_induced[gene] = vitd_sig - ctrl_sig
    print("  HaCaT VitD-induced eRNA (top genes):")
    top = sorted(hackat_vitd_induced.items(), key=lambda x: -x[1])[:15]
    for g, s in top:
        print(f"    {g}: {s:.4f}")
else:
    hackat_vitd_induced = {g: 0.0 for g in TYPE2_GENES}
    print("  BigWig download incomplete, skipping eRNA analysis")

# ── Step 2C: ENCODE GATA3 ChIP-seq ──────────────────────────────────────────
print("\n[Step 2C] ENCODE GATA3 ChIP-seq (Th2 lymphocyte)...")

GATA3_ENCODE_URLS = [
    # ENCODE GATA3 ChIP-seq in Th2/lymphocyte cells
    # GM12878 (B-LCL), K562 not ideal; looking for CD4/Th2-specific
    # These are the best available GATA3 ChIP-seq in lymphoid context
    "https://www.encodeproject.org/files/ENCFF421HCD/@@download/ENCFF421HCD.bed.gz",  # GATA3 GM12878
    "https://www.encodeproject.org/files/ENCFF879OKY/@@download/ENCFF879OKY.bed.gz",  # GATA3 K562
]

# より直接的なアプローチ: GEOからGATA3 ChIP-seq in Th2 cells
# GSE40463: GATA3 ChIP-seq in human Th2 cells
GATA3_GEO_QUERY = "GATA3+ChIP-seq+Th2"

def fetch_gata3_geo():
    url = f'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gds&term=GATA3+ChIP-seq+Th2+human&retmax=10&retmode=json'
    try:
        res = urllib.request.urlopen(url, timeout=15)
        data = json.loads(res.read())
        ids = data.get('esearchresult',{}).get('idlist',[])
        url2 = f'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=gds&id={",".join(ids[:5])}&retmode=json'
        res2 = urllib.request.urlopen(url2, timeout=15)
        d2 = json.loads(res2.read())
        results = []
        for gid in ids[:5]:
            r = d2.get('result',{}).get(gid,{})
            results.append({'acc': r.get('accession',''), 'title': r.get('title','')})
        return results
    except Exception as e:
        return []

gata3_datasets = fetch_gata3_geo()
print("  GATA3 ChIP-seq datasets found:")
for d in gata3_datasets:
    print(f"    {d['acc']}: {d['title'][:80]}")

# GSE67273 or similar - GATA3 ChIP-seq in CD4+ Th2
# Try to find and download processed peak file
GATA3_CANDIDATES = [
    ('GSE67273', 'GATA3 Th2 CD4+'),
    ('GSE40463', 'GATA3 Th2'),
    ('GSE14068', 'GATA3 Th2 primary'),
]

gata3_scores = {g: 0 for g in TYPE2_GENES}
gata3_found = False

for gse, desc in GATA3_CANDIDATES:
    try:
        url = f'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={gse}&targ=gsm&form=text&view=quick'
        res = urllib.request.urlopen(url, timeout=10)
        txt = res.read().decode()
        # サンプルファイルを探す
        bed_urls = re.findall(r'ftp://[^\s]+\.bed\.gz', txt)
        narrowpeak_urls = re.findall(r'ftp://[^\s]+\.narrowPeak\.gz', txt)
        all_urls = bed_urls + narrowpeak_urls
        if all_urls:
            print(f"  {gse} ({desc}): found {len(all_urls)} BED files")
            for u in all_urls[:3]:
                print(f"    {u}")
            # 最初のファイルをダウンロード
            fname = all_urls[0].split('/')[-1]
            dest = TMP / fname
            if download_file(all_urls[0], dest):
                # スコア計算
                print(f"  Scoring {gse} GATA3 peaks...")
                with gzip.open(dest, 'rt') as f:
                    for line in f:
                        if line.startswith('#'): continue
                        parts = line.split('\t')
                        if len(parts) < 3: continue
                        try:
                            chrom, start, end = parts[0], int(parts[1]), int(parts[2])
                        except: continue
                        for gene in TYPE2_GENES:
                            if gene not in gene_coords: continue
                            gc = gene_coords[gene]
                            if gc['chrom'] != chrom: continue
                            if start < gc['prom_end'] and end > gc['prom_start']:
                                gata3_scores[gene] += 1
                gata3_found = True
                print(f"  GATA3 scores (top genes):")
                top_g3 = sorted([(g, gata3_scores[g]) for g in TYPE2_GENES if g in gata3_scores], key=lambda x: -x[1])[:10]
                for g, s in top_g3:
                    print(f"    {g}: {s}")
                break
    except Exception as e:
        print(f"  {gse}: {e}")
    if gata3_found: break

if not gata3_found:
    print("  GATA3 GEO data not found via direct URL. Trying ENCODE...")
    # ENCODE GATA3 peak files
    encode_gata3_beds = [
        "https://www.encodeproject.org/files/ENCFF421HCD/@@download/ENCFF421HCD.bed.gz",
    ]
    for url in encode_gata3_beds:
        fname = url.split('/')[-1]
        dest = TMP / fname
        try:
            if download_file(url, dest):
                with gzip.open(dest,'rt') as f:
                    for line in f:
                        if line.startswith('#'): continue
                        parts = line.split('\t')
                        if len(parts) < 3: continue
                        try: chrom,start,end = parts[0],int(parts[1]),int(parts[2])
                        except: continue
                        for gene in TYPE2_GENES:
                            if gene not in gene_coords: continue
                            gc = gene_coords[gene]
                            if gc['chrom'] != chrom: continue
                            if start < gc['prom_end'] and end > gc['prom_start']:
                                gata3_scores[gene] += 1
                gata3_found = True
                break
        except Exception as e:
            print(f"  ENCODE GATA3: {e}")

# ── Step 3: 結果統合 ──────────────────────────────────────────────────────────
print("\n[Step 3] Integrating results...")

# Type 2ラベル
TYPE2_LABELS = {
    "TSLP": ("Type2-Alarmin", "上皮→ILC2"),
    "IL33": ("Type2-Alarmin", "上皮→ILC2/肥満細胞"),
    "IL25": ("Type2-Alarmin", "上皮→ILC2"),
    "IL4R": ("Type2-Effector", "Th2/上皮"),
    "IL13": ("Type2-Effector", "Th2/ILC2"),
    "IL5": ("Type2-Effector", "Th2/ILC2"),
    "IL4": ("Type2-Effector", "Th2"),
    "IL31RA": ("Type2-Effector", "皮膚/Th2"),
    "IGHE": ("Type2-IgE", "B細胞"),
    "CCR3": ("Type2-Eosinophil", "好酸球"),
    "SIGLEC8": ("Type2-Eosinophil", "好酸球"),
    "IL1RL1": ("Type2-Receptor", "Th2/ILC2"),
    "CPA3": ("Type2-MastCell", "肥満細胞"),
    "FCER1A": ("Type2-IgE-Rec", "肥満細胞/好酸球"),
    "IL17RB": ("Type2-ILC2", "ILC2"),
    "HHIP": ("Type2-Other", "上皮"),
    "AREG": ("Type2-Other", "Th2/ILC2"),
    "S100A8": ("Type2-Alarmin2", "上皮/好中球"),
    "S100A9": ("Type2-Alarmin2", "上皮/好中球"),
    "IL23A": ("Myeloid", "マクロファージ"),
    "TNF": ("Myeloid", "マクロファージ"),
    "IL6R": ("Myeloid", "マクロファージ"),
    "VEGFA": ("Myeloid", "マクロファージ"),
    "CD274": ("Myeloid", "マクロファージ"),
    "IL10": ("Myeloid", "マクロファージ"),
    "TYK2": ("Myeloid", "マクロファージ"),
    "ITGB7": ("Myeloid", "リンパ球"),
    "C5AR1": ("Myeloid", "マクロファージ"),
}

# 既知の薬剤承認状況
DRUG_STATUS = {
    "TSLP": "承認済(テゼペルマブ)", "IL33": "承認済(イテペキマブ)",
    "IL25": "Ph2中", "IL4R": "承認済(デュピクセント)",
    "IL13": "承認済(デュピクセント/トラロキヌマブ)", "IL5": "承認済(メポリズマブ)",
    "IL4": "承認済(IL4R経由)", "IL31RA": "承認済(ネモリズマブ)",
    "IGHE": "承認済(オマリズマブ)", "CCR3": "Ph2失敗(20年間)",
    "SIGLEC8": "承認済(lirentelimab)", "IL1RL1": "Ph2中",
    "CPA3": "未開発", "FCER1A": "未開発", "IL17RB": "Ph2中",
    "HHIP": "未開発", "AREG": "未開発", "S100A8": "未開発", "S100A9": "未開発",
    "IL23A": "承認済(リサンキズマブ)", "TNF": "承認済(インフリキシマブ)",
    "IL6R": "承認済(トシリズマブ)", "VEGFA": "承認済(ベバシズマブ)",
    "CD274": "承認済(アテゾリズマブ)", "IL10": "未承認",
    "TYK2": "承認済(デュークラバシチニブ)", "ITGB7": "承認済(ベドリズマブ)",
    "C5AR1": "承認済(アバコパン)",
}

rows = []
for gene in TYPE2_GENES:
    if gene not in gene_coords:
        continue
    cat, cell = TYPE2_LABELS.get(gene, ("Unknown", "?"))
    rows.append({
        '遺伝子': gene,
        'カテゴリ': cat,
        '関連細胞': cell,
        'VDR_THP1': scores_thp1_vdr.get(gene, 0),
        'VDR_Epithelial': scores_epithelial_vdr.get(gene, 0),
        'VDR_AllCells': scores_all_vdr.get(gene, 0),
        'GR_AllCells': scores_all_gr.get(gene, 0),
        'HaCaT_VitD_eRNA': round(hackat_vitd_induced.get(gene, 0.0), 4),
        'GATA3_score': gata3_scores.get(gene, 0),
        '薬剤状況': DRUG_STATUS.get(gene, '?'),
    })

df = pd.DataFrame(rows)
df['VDR_THP1_pattern'] = df.apply(lambda r:
    'VDR dominant' if r['VDR_THP1'] > r['GR_AllCells'] else
    'GR dominant' if r['GR_AllCells'] > r['VDR_THP1'] else 'Equal', axis=1)

# ── Step 4: Excel出力 ─────────────────────────────────────────────────────────
print("\n[Step 4] Saving results...")

import openpyxl
from openpyxl.styles import PatternFill, Font, Alignment, Border, Side
from openpyxl.utils import get_column_letter
from openpyxl.formatting.rule import ColorScaleRule

wb = openpyxl.Workbook()

# シート1: 詳細データ
ws1 = wb.active
ws1.title = "Type2 VDR Analysis"

headers = list(df.columns)
thin = Side(border_style='thin', color='CCCCCC')
hdr_fill = PatternFill('solid', fgColor='1B3A5C')
hdr_font = Font(color='FFFFFF', bold=True, size=9)
for j, h in enumerate(headers, 1):
    c = ws1.cell(1, j, h)
    c.fill = hdr_fill; c.font = hdr_font
    c.alignment = Alignment(horizontal='center', wrap_text=True)
    c.border = Border(left=thin,right=thin,top=thin,bottom=thin)

cat_fills = {
    'Type2-Alarmin': 'FFF3E0',   # 上皮アラーミン
    'Type2-Effector': 'E3F2FD',  # エフェクターサイトカイン
    'Type2-IgE': 'FCE4EC',
    'Type2-Eosinophil': 'F3E5F5',
    'Type2-MastCell': 'FFF9C4',
    'Type2-IgE-Rec': 'FCE4EC',
    'Type2-ILC2': 'E8F5E9',
    'Type2-Receptor': 'E0F7FA',
    'Type2-Alarmin2': 'FFF3E0',
    'Type2-Other': 'FAFAFA',
    'Myeloid': 'D4EDDA',
    'Unknown': 'FFFFFF',
}

for i, row in df.iterrows():
    xl_row = i + 2
    fill_color = cat_fills.get(row['カテゴリ'], 'FFFFFF')
    fill = PatternFill('solid', fgColor=fill_color)
    for j, col in enumerate(headers, 1):
        val = row[col]
        c = ws1.cell(xl_row, j, val)
        c.fill = fill
        c.border = Border(left=thin,right=thin,top=thin,bottom=thin)
        c.alignment = Alignment(horizontal='center')
        if col == '遺伝子':
            c.font = Font(bold=True, size=10)
        elif 'VDR' in col or 'GR' in col or 'GATA3' in col:
            c.font = Font(size=10)
        else:
            c.font = Font(size=9)

# 列幅
col_widths = [10, 18, 14, 10, 14, 11, 11, 14, 12, 24, 16]
for j, w in enumerate(col_widths[:len(headers)], 1):
    ws1.column_dimensions[get_column_letter(j)].width = w
ws1.freeze_panes = 'A2'

# カラースケール: VDR_Epithelial列
if 'VDR_Epithelial' in headers:
    col_idx = headers.index('VDR_Epithelial') + 1
    col_letter = get_column_letter(col_idx)
    n_rows = len(df)
    ws1.conditional_formatting.add(
        f'{col_letter}2:{col_letter}{n_rows+1}',
        ColorScaleRule(start_type='min', start_color='FFFFFF',
                       end_type='max', end_color='1565C0')
    )

# シート2: 要約
ws2 = wb.create_sheet("Summary")
ws2['A1'] = "=== Type 2 vs Myeloid: VDR Binding Analysis ==="
ws2['A1'].font = Font(bold=True, size=12)
ws2['A3'] = "Category"
ws2['B3'] = "n genes"
ws2['C3'] = "VDR THP-1 mean"
ws2['D3'] = "VDR Epithelial mean"
ws2['E3'] = "GR All mean"
ws2['F3'] = "HaCaT VitD-eRNA mean"
ws2['G3'] = "GATA3 mean"
for c in ws2['A3:G3'][0]:
    c.fill = PatternFill('solid', fgColor='2C3E50')
    c.font = Font(color='FFFFFF', bold=True, size=10)
    c.alignment = Alignment(horizontal='center')

cats_order = ['Type2-Alarmin', 'Type2-Effector', 'Type2-Eosinophil',
              'Type2-IgE', 'Type2-MastCell', 'Type2-ILC2', 'Type2-Receptor', 'Myeloid']
row_idx = 4
for cat in cats_order:
    sub = df[df['カテゴリ'] == cat]
    if len(sub) == 0: continue
    ws2.cell(row_idx, 1, cat)
    ws2.cell(row_idx, 2, len(sub))
    ws2.cell(row_idx, 3, round(sub['VDR_THP1'].mean(), 1))
    ws2.cell(row_idx, 4, round(sub['VDR_Epithelial'].mean(), 1))
    ws2.cell(row_idx, 5, round(sub['GR_AllCells'].mean(), 1))
    ws2.cell(row_idx, 6, round(sub['HaCaT_VitD_eRNA'].mean(), 4))
    ws2.cell(row_idx, 7, round(sub['GATA3_score'].mean(), 1))
    fill = PatternFill('solid', fgColor='E3F2FD' if 'Type2' in cat else 'D4EDDA')
    for c in ws2[row_idx]:
        c.fill = fill
        c.alignment = Alignment(horizontal='center')
        c.font = Font(size=10)
    row_idx += 1

# Key findings
ws2.cell(row_idx+2, 1, "=== Key Findings ===").font = Font(bold=True, size=11)
row_idx += 3

findings = []
# TSLP: VDR優位?
tslp = df[df['遺伝子']=='TSLP']
if len(tslp):
    t = tslp.iloc[0]
    findings.append(f"TSLP: VDR_THP1={t['VDR_THP1']}, VDR_Epithelial={t['VDR_Epithelial']}, HaCaT_eRNA={t['HaCaT_VitD_eRNA']:.4f}")

il33 = df[df['遺伝子']=='IL33']
if len(il33):
    t = il33.iloc[0]
    findings.append(f"IL33: VDR_THP1={t['VDR_THP1']}, VDR_Epithelial={t['VDR_Epithelial']}, HaCaT_eRNA={t['HaCaT_VitD_eRNA']:.4f}")

ccr3 = df[df['遺伝子']=='CCR3']
if len(ccr3):
    t = ccr3.iloc[0]
    findings.append(f"CCR3 (20y failure): VDR_THP1={t['VDR_THP1']}, GATA3={t['GATA3_score']}")

for finding in findings:
    ws2.cell(row_idx, 1, finding)
    ws2.cell(row_idx, 1).font = Font(size=10)
    row_idx += 1

for j in range(1, 8):
    ws2.column_dimensions[get_column_letter(j)].width = 22

out_path = OUT_DIR / "type2_vdr_chipseq_results.xlsx"
wb.save(out_path)
print(f"\n✅ Saved: {out_path}")

# CSV保存
csv_path = OUT_DIR / "type2_vdr_chipseq_results.csv"
df.to_csv(csv_path, index=False)
print(f"✅ Saved: {csv_path}")

print("\n=== Final Summary ===")
print(df[['遺伝子','カテゴリ','VDR_THP1','VDR_Epithelial','GR_AllCells','HaCaT_VitD_eRNA','GATA3_score','薬剤状況']].to_string(index=False))
print("\nDone!")
