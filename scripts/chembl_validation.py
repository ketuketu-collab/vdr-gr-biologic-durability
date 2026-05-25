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
            res = urllib.request.urlopen(req, timeout=15)
            return json.loads(res.read())
        except:
            if i < retries-1: time.sleep(2)
    return None

def get_target_id(gene):
    url = f"{BASE}/target/search.json?q={gene}&target_type=SINGLE+PROTEIN&limit=5"
    data = chembl_get(url)
    if not data: return None, None
    for t in data.get('targets', []):
        for comp in t.get('target_components', []):
            for syn in comp.get('target_component_synonyms', []):
                if syn.get('syn_type') == 'GENE_SYMBOL' and syn.get('component_synonym','').upper() == gene.upper():
                    return t['target_chembl_id'], t.get('pref_name','')
    targets = data.get('targets', [])
    if targets: return targets[0]['target_chembl_id'], targets[0].get('pref_name','')
    return None, None

def get_drugs(target_id):
    drugs = []
    url = f"{BASE}/mechanism.json?target_chembl_id={target_id}&limit=200"
    data = chembl_get(url)
    if not data: return drugs
    for m in data.get('mechanisms', []):
        mol_id = m.get('molecule_chembl_id')
        if not mol_id: continue
        mol = chembl_get(f"{BASE}/molecule/{mol_id}.json")
        if not mol: continue
        drugs.append({
            'mol_id': mol_id,
            'name': mol.get('pref_name',''),
            'max_phase': mol.get('max_phase', 0) or 0,
            'mol_type': mol.get('molecule_type',''),
            'first_approval': mol.get('first_approval',''),
            'withdrawn': mol.get('withdrawn_flag', False),
            'moa': m.get('mechanism_of_action',''),
            'action': m.get('action_type','')
        })
        time.sleep(0.1)
    return drugs

results = []
for i, gene in enumerate(genes):
    tid, tname = get_target_id(gene)
    if not tid:
        results.append({'gene': gene, 'target_chembl_id': '', 'target_name': '',
            'n_total': 0, 'n_approved': 0, 'n_phase3': 0, 'n_phase2': 0, 'n_phase1': 0,
            'n_withdrawn': 0, 'approved_names': '', 'phase3_names': '', 'withdrawn_names': ''})
    else:
        drugs = get_drugs(tid)
        app = [d for d in drugs if d['max_phase'] == 4]
        p3  = [d for d in drugs if d['max_phase'] == 3]
        p2  = [d for d in drugs if d['max_phase'] == 2]
        p1  = [d for d in drugs if d['max_phase'] == 1]
        wd  = [d for d in drugs if d['withdrawn']]
        results.append({'gene': gene, 'target_chembl_id': tid, 'target_name': tname,
            'n_total': len(drugs), 'n_approved': len(app), 'n_phase3': len(p3),
            'n_phase2': len(p2), 'n_phase1': len(p1), 'n_withdrawn': len(wd),
            'approved_names': '; '.join([d['name'] for d in app if d['name']]),
            'phase3_names': '; '.join([d['name'] for d in p3 if d['name']]),
            'withdrawn_names': '; '.join([d['name'] for d in wd if d['name']])})
    if (i+1) % 20 == 0:
        print(f"  {i+1}/{len(genes)} done...")
    time.sleep(0.3)

df_chembl = pd.DataFrame(results)
df_merged = df_orig.merge(df_chembl, left_on='遺伝子', right_on='gene', how='left')
out = '/Volumes/M4_SSD/projects/tlr_chipseq/results/chembl_validation_308genes.csv'
df_merged.to_csv(out, index=False)
print(f"\nSaved: {out}")
print(f"Genes with approved drugs: {(df_chembl.n_approved > 0).sum()}/{len(df_chembl)}")

# VDR vs GR approval rate comparison
vdr_pat = df_merged['パターン'].fillna('').str.contains('VDR')
gr_pat  = df_merged['パターン'].fillna('').str.contains('GR') & ~vdr_pat
vdr_tot = df_merged[vdr_pat]['n_total'].sum()
vdr_app = df_merged[vdr_pat]['n_approved'].sum()
gr_tot  = df_merged[gr_pat]['n_total'].sum()
gr_app  = df_merged[gr_pat]['n_approved'].sum()
print(f"\nVDR-dominant: {vdr_app}/{vdr_tot} approved ({100*vdr_app/vdr_tot:.1f}%)" if vdr_tot else "")
print(f"GR-dominant:  {gr_app}/{gr_tot} approved ({100*gr_app/gr_tot:.1f}%)" if gr_tot else "")
