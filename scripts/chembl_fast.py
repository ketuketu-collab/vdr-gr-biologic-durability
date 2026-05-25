"""
ChEMBL fast validation - uses count-only API endpoints (no per-drug lookup).
Checks approved drugs for all 308 genes in ~5-10 min instead of 2h.
"""
import pandas as pd
import urllib.request
import json
import time

df_orig = pd.read_excel('/tmp/steroid_309genes.xlsx', sheet_name='全体(308遺伝子)VDR-GR降順')
genes = df_orig['遺伝子'].dropna().tolist()
print(f"Total genes: {len(genes)}")

BASE = "https://www.ebi.ac.uk/chembl/api/data"

def chembl_get(url, retries=3):
    for i in range(retries):
        try:
            req = urllib.request.Request(url, headers={'Accept': 'application/json'})
            res = urllib.request.urlopen(req, timeout=20)
            return json.loads(res.read())
        except Exception as e:
            if i < retries-1:
                time.sleep(3)
    return None

def get_target_id(gene):
    url = f"{BASE}/target/search.json?q={gene}&target_type=SINGLE+PROTEIN&limit=5"
    data = chembl_get(url)
    if not data:
        return None, None
    # Exact gene symbol match first
    for t in data.get('targets', []):
        for comp in t.get('target_components', []):
            for syn in comp.get('target_component_synonyms', []):
                if syn.get('syn_type') == 'GENE_SYMBOL' and syn.get('component_synonym', '').upper() == gene.upper():
                    return t['target_chembl_id'], t.get('pref_name', '')
    targets = data.get('targets', [])
    if targets:
        return targets[0]['target_chembl_id'], targets[0].get('pref_name', '')
    return None, None

def get_drug_counts_fast(target_id):
    """
    Use mechanism.json to get all drugs for a target.
    Then use molecule count endpoint (no per-molecule detail fetch).
    Returns counts by phase and approved drug names.
    """
    url = f"{BASE}/mechanism.json?target_chembl_id={target_id}&limit=500"
    data = chembl_get(url)
    if not data:
        return 0, 0, 0, 0, 0, [], []

    mechanisms = data.get('mechanisms', [])
    mol_ids = list(set(m['molecule_chembl_id'] for m in mechanisms if m.get('molecule_chembl_id')))

    if not mol_ids:
        return 0, 0, 0, 0, 0, [], []

    # Batch fetch molecules using filter endpoint
    id_list = ','.join(mol_ids[:100])  # ChEMBL allows comma-separated IDs
    mol_url = f"{BASE}/molecule.json?molecule_chembl_id__in={id_list}&limit=500"
    mol_data = chembl_get(mol_url)

    if not mol_data:
        return 0, 0, 0, 0, 0, [], []

    molecules = mol_data.get('molecules', [])

    counts = {1: 0, 2: 0, 3: 0, 4: 0, 0: 0}
    approved_names = []
    phase3_names = []

    for mol in molecules:
        phase = mol.get('max_phase', 0) or 0
        try:
            phase = int(phase)
        except:
            phase = 0
        counts[phase] = counts.get(phase, 0) + 1
        name = mol.get('pref_name', '') or mol.get('molecule_chembl_id', '')
        if phase == 4:
            approved_names.append(name)
        elif phase == 3:
            phase3_names.append(name)

    n_approved = counts.get(4, 0)
    n_phase3   = counts.get(3, 0)
    n_phase2   = counts.get(2, 0)
    n_phase1   = counts.get(1, 0)
    n_total    = sum(counts.values())

    return n_total, n_approved, n_phase3, n_phase2, n_phase1, approved_names, phase3_names

# Load existing partial results if available
existing = {}
import os
prev_file = '/Volumes/M4_SSD/projects/tlr_chipseq/results/chembl_drug_phases.csv'
if os.path.exists(prev_file):
    df_prev = pd.read_csv(prev_file)
    for _, row in df_prev.iterrows():
        existing[row['gene']] = row
    print(f"Loaded {len(existing)} existing results")

results = []
for i, gene in enumerate(genes):
    if (i+1) % 10 == 0:
        print(f"  {i+1}/{len(genes)} done...")

    # Use existing data if available
    if gene in existing:
        row = existing[gene]
        results.append({
            'gene': gene,
            'target_chembl_id': row.get('chembl_id', ''),
            'target_name': '',
            'n_total': int(row.get('total', 0)),
            'n_approved': int(row.get('approved', 0)),
            'n_phase3': int(row.get('p3', 0)),
            'n_phase2': int(row.get('p2', 0)),
            'n_phase1': int(row.get('p1', 0)),
            'n_withdrawn': 0,
            'approved_names': '',
            'phase3_names': '',
            'withdrawn_names': ''
        })
        continue

    tid, tname = get_target_id(gene)
    if not tid:
        results.append({
            'gene': gene, 'target_chembl_id': '', 'target_name': '',
            'n_total': 0, 'n_approved': 0, 'n_phase3': 0, 'n_phase2': 0, 'n_phase1': 0,
            'n_withdrawn': 0, 'approved_names': '', 'phase3_names': '', 'withdrawn_names': ''
        })
    else:
        n_total, n_approved, n_phase3, n_phase2, n_phase1, app_names, p3_names = get_drug_counts_fast(tid)
        results.append({
            'gene': gene, 'target_chembl_id': tid, 'target_name': tname,
            'n_total': n_total, 'n_approved': n_approved, 'n_phase3': n_phase3,
            'n_phase2': n_phase2, 'n_phase1': n_phase1, 'n_withdrawn': 0,
            'approved_names': '; '.join(app_names),
            'phase3_names': '; '.join(p3_names),
            'withdrawn_names': ''
        })
    time.sleep(0.5)

df_chembl = pd.DataFrame(results)
df_merged = df_orig.merge(df_chembl, left_on='遺伝子', right_on='gene', how='left')

out = '/Volumes/M4_SSD/projects/tlr_chipseq/results/chembl_validation_308genes.csv'
df_merged.to_csv(out, index=False)
print(f"\nSaved: {out}")
print(f"Total genes: {len(df_chembl)}")
print(f"Genes with approved drugs: {(df_chembl.n_approved > 0).sum()}/{len(df_chembl)}")
print(f"Genes with any drug activity: {(df_chembl.n_total > 0).sum()}/{len(df_chembl)}")

# VDR vs GR approval rate comparison
col_vdr = 'VDR_score' if 'VDR_score' in df_merged.columns else None
col_gr  = 'GR_score'  if 'GR_score'  in df_merged.columns else None
pat_col = next((c for c in df_merged.columns if 'パターン' in c or 'Pattern' in c), None)

if pat_col:
    vdr_mask = df_merged[pat_col].fillna('').str.contains('VDR') & ~df_merged[pat_col].fillna('').str.contains('GR')
    gr_mask  = df_merged[pat_col].fillna('').str.contains('GR')  & ~df_merged[pat_col].fillna('').str.contains('VDR')
    vdr_genes = df_merged[vdr_mask]
    gr_genes  = df_merged[gr_mask]
    vdr_app = (vdr_genes['n_approved'] > 0).sum()
    vdr_tot = len(vdr_genes)
    gr_app  = (gr_genes['n_approved'] > 0).sum()
    gr_tot  = len(gr_genes)
    print(f"\n=== ChEMBL Approval Rates ===")
    print(f"VDR-only: {vdr_app}/{vdr_tot} genes with approved drugs ({100*vdr_app/vdr_tot:.1f}%)" if vdr_tot else "VDR-only: N/A")
    print(f"GR-only:  {gr_app}/{gr_tot} genes with approved drugs ({100*gr_app/gr_tot:.1f}%)" if gr_tot else "GR-only: N/A")
else:
    print("\n(Pattern column not found for VDR/GR split)")

print("\nTop genes by approved drugs:")
print(df_chembl.nlargest(20, 'n_approved')[['gene','n_total','n_approved','n_phase3','approved_names']].to_string())
