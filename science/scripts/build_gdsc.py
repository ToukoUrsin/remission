"""Build gdsc_long.csv, gdsc_targets.json, famous_lines.csv from GDSC2 release 8.5 + Cell Model Passports.

Usage: python3 -I scripts/build_gdsc.py <science_dir>
Inputs (raw/, untrusted data, parsed only):
  raw/gdsc/GDSC2_fitted_dose_response_27Oct23.csv  (converted from the xlsx)
  raw/gdsc/Cell_Lines_Details.xlsx
  raw/cmp/model_list_latest.csv.gz
  raw/cmp/mutations_selected_genes.csv  (rows of mutations_all_latest.csv.gz for BRAF/KRAS/NRAS/HRAS/PIK3CA/PTEN/EGFR)
  raw/cmp/cnv_summary_latest.csv.gz
"""
import json, re, sys, os
import numpy as np
import pandas as pd

S = sys.argv[1]
R = os.path.join(S, 'raw')

DRUGS = {
    # class: [(drug_name, drug_id, role)]
    'RAF_inhibitor': [('Dabrafenib', 1373, 'requested'), ('PLX-4720', 1036, 'requested'), ('SB590885', 1061, 'extra')],
    'MEK_inhibitor': [('Trametinib', 1372, 'requested'), ('Selumetinib', 1062, 'requested'), ('Selumetinib', 1736, 'requested'),
                      ('PD0325901', 1060, 'extra'), ('Refametinib', 1014, 'extra')],
    'EGFR_inhibitor': [('Gefitinib', 1010, 'requested'), ('Erlotinib', 1168, 'requested'), ('Afatinib', 1032, 'requested'),
                       ('Osimertinib', 1919, 'extra'), ('Lapatinib', 1558, 'extra')],
    'PI3K_inhibitor': [('Alpelisib', 1560, 'requested'), ('Pictilisib', 1058, 'requested'), ('Taselisib', 1561, 'requested'),
                       ('Buparlisib', 1873, 'extra')],
    'BCL2_family_inhibitor': [('Navitoclax', 1011, 'requested'), ('Venetoclax', 1909, 'requested'), ('WEHI-539', 1997, 'requested_BCLxL_selective'),
                              ('ABT737', 1910, 'extra')],
}
NOT_IN_GDSC2 = ['Vemurafenib', 'Cetuximab', 'Encorafenib', 'Binimetinib', 'Panitumumab']
id2 = {did: (cls, name, role) for cls, L in DRUGS.items() for name, did, role in L}

d = pd.read_csv(os.path.join(R, 'gdsc/GDSC2_fitted_dose_response_27Oct23.csv'))
d = d[d.DRUG_ID.isin(id2)].copy()
for did, (cls, name, role) in id2.items():
    got = d.loc[d.DRUG_ID == did, 'DRUG_NAME'].unique().tolist()
    assert got == [name], (did, name, got)
d['drug_class'] = d.DRUG_ID.map(lambda x: id2[x][0])
d['drug_role'] = d.DRUG_ID.map(lambda x: id2[x][2])
d['drug_label'] = np.where(d.DRUG_NAME == 'Selumetinib', 'Selumetinib_' + d.DRUG_ID.astype(str), d.DRUG_NAME)

# --- annotations
ml = pd.read_csv(os.path.join(R, 'cmp/model_list_latest.csv.gz'), low_memory=False)
ml = ml.set_index('model_id')
cld = pd.read_excel(os.path.join(R, 'gdsc/Cell_Lines_Details.xlsx'), sheet_name='Cell line details')
cld.columns = [' '.join(str(c).split()) for c in cld.columns]
cld = cld.dropna(subset=['COSMIC identifier']).set_index(cld['COSMIC identifier'].dropna().astype(int))

mut = pd.read_csv(os.path.join(R, 'cmp/mutations_selected_genes.csv'), low_memory=False)
mut = mut[mut.coding == 't']
HOT_RAS = {12, 13, 59, 61, 117, 146}

def flags(g):
    out = {}
    br = g[g.gene_symbol == 'BRAF']
    out['BRAF_V600E'] = bool((br.protein_mutation == 'p.V600E').any())
    out['BRAF_V600_any'] = bool(br.protein_mutation.str.match(r'p\.V600[A-Z]$', na=False).any())
    out['BRAF_driver_any'] = bool((br.cancer_driver == 't').any())
    for gene in ['KRAS', 'NRAS', 'HRAS']:
        x = g[(g.gene_symbol == gene) & (g.effect == 'missense')]
        out[gene + '_mut'] = bool(x.aa_pos.isin(HOT_RAS).any())
    p = g[(g.gene_symbol == 'PIK3CA') & (g.cancer_driver == 't') & g.effect.isin(['missense', 'inframe'])]
    out['PIK3CA_mut'] = bool(len(p))
    pt = g[(g.gene_symbol == 'PTEN') & g.effect.isin(['nonsense', 'frameshift', 'ess_splice', 'start_lost', 'stop_lost'])]
    out['PTEN_truncating'] = bool(len(pt))
    eg = g[(g.gene_symbol == 'EGFR') & (g.cancer_driver == 't')]
    out['EGFR_driver_mut'] = bool(len(eg))
    def s(gene):
        x = g[(g.gene_symbol == gene) & (g.cancer_driver == 't')].protein_mutation.dropna().unique()
        return ';'.join(sorted(x))
    for gene in ['BRAF', 'KRAS', 'NRAS', 'HRAS', 'PIK3CA', 'PTEN', 'EGFR']:
        out[gene + '_driver_variants'] = s(gene)
    return pd.Series(out)

F = mut.groupby('model_id').apply(flags, include_groups=False)
has_mut = set(ml.index[ml.mutation_data == True]) | set(F.index)

cnv = pd.read_csv(os.path.join(R, 'cmp/cnv_summary_latest.csv.gz'))
pten_del = set(cnv[(cnv.symbol == 'PTEN') & (cnv.cn_category == 'Deletion')].model_id)
pten_cnv_known = set(cnv[cnv.symbol == 'PTEN'].model_id)

lines = d[['CELL_LINE_NAME', 'SANGER_MODEL_ID', 'COSMIC_ID', 'TCGA_DESC']].drop_duplicates('SANGER_MODEL_ID').set_index('SANGER_MODEL_ID')
A = lines.copy()
A['cmp_tissue'] = ml.reindex(A.index).tissue
A['cmp_cancer_type'] = ml.reindex(A.index).cancer_type
A['gdsc_tissue_1'] = [cld['GDSC Tissue descriptor 1'].get(c, np.nan) for c in A.COSMIC_ID]
A['gdsc_tissue_2'] = [cld['GDSC Tissue descriptor 2'].get(c, np.nan) for c in A.COSMIC_ID]
A['msi_status'] = [cld['Microsatellite instability Status (MSI)'].get(c, np.nan) for c in A.COSMIC_ID]
A['mutation_data_available'] = A.index.isin(has_mut)
fl = F.reindex(A.index)
boolcols = [c for c in F.columns if not c.endswith('_variants')]
for c in boolcols:
    fl[c] = fl[c].astype('boolean').fillna(False)
    fl.loc[~A.mutation_data_available, c] = pd.NA
for c in [c for c in F.columns if c.endswith('_variants')]:
    fl[c] = fl[c].fillna('')
A = A.join(fl)
A['PTEN_homdel'] = pd.array([(m in pten_del) if m in pten_cnv_known else pd.NA for m in A.index], dtype='boolean')
A['PTEN_null'] = (A.PTEN_truncating.fillna(False) | A.PTEN_homdel.fillna(False)).astype('boolean')
A.loc[A.PTEN_truncating.isna() & A.PTEN_homdel.isna(), 'PTEN_null'] = pd.NA
A['RAS_RAF_wt'] = ~(A.KRAS_mut | A.NRAS_mut | A.HRAS_mut | A.BRAF_driver_any)
A['solid'] = A.cmp_tissue.notna() & (A.cmp_tissue != 'Haematopoietic and Lymphoid') & (A.gdsc_tissue_1 != 'blood')


# --- annotation table for every CMP model (used by build_combo.py)
AL = ml[['model_name', 'COSMIC_ID', 'tissue', 'cancer_type', 'msi_status']].copy()
AL['mutation_data_available'] = AL.index.isin(has_mut)
fa = F.reindex(AL.index)
for c in boolcols:
    fa[c] = fa[c].astype('boolean').fillna(False)
    fa.loc[~AL.mutation_data_available, c] = pd.NA
for c in [c for c in F.columns if c.endswith('_variants')]:
    fa[c] = fa[c].fillna('')
AL = AL.join(fa)
AL['PTEN_homdel'] = pd.array([(m in pten_del) if m in pten_cnv_known else pd.NA for m in AL.index], dtype='boolean')
AL['PTEN_null'] = (AL.PTEN_truncating.fillna(False) | AL.PTEN_homdel.fillna(False)).astype('boolean')
AL.loc[AL.PTEN_truncating.isna() & AL.PTEN_homdel.isna(), 'PTEN_null'] = pd.NA
AL['RAS_RAF_wt'] = ~(AL.KRAS_mut | AL.NRAS_mut | AL.HRAS_mut | AL.BRAF_driver_any)
AL.index.name = 'model_id'
AL.to_csv(os.path.join(S, 'cell_line_annotations.csv'))

# --- long table
L = d.merge(A.drop(columns=['CELL_LINE_NAME', 'COSMIC_ID', 'TCGA_DESC']), left_on='SANGER_MODEL_ID', right_index=True, how='left')
L['ln_max_conc'] = np.log(L.MAX_CONC)
L['IC50_below_max'] = L.LN_IC50 < L.ln_max_conc
cols = ['CELL_LINE_NAME', 'COSMIC_ID', 'SANGER_MODEL_ID', 'TCGA_DESC', 'cmp_tissue', 'cmp_cancer_type', 'gdsc_tissue_1', 'gdsc_tissue_2', 'msi_status',
        'mutation_data_available', 'BRAF_V600E', 'KRAS_mut', 'NRAS_mut', 'PIK3CA_mut', 'PTEN_null',
        'BRAF_V600_any', 'HRAS_mut', 'EGFR_driver_mut', 'PTEN_truncating', 'PTEN_homdel', 'RAS_RAF_wt', 'solid',
        'BRAF_driver_variants', 'KRAS_driver_variants', 'NRAS_driver_variants', 'HRAS_driver_variants', 'PIK3CA_driver_variants', 'PTEN_driver_variants', 'EGFR_driver_variants',
        'drug_class', 'drug_label', 'DRUG_NAME', 'DRUG_ID', 'drug_role', 'PUTATIVE_TARGET',
        'LN_IC50', 'AUC', 'RMSE', 'Z_SCORE', 'MIN_CONC', 'MAX_CONC', 'IC50_below_max']
L = L[cols].sort_values(['drug_class', 'drug_label', 'CELL_LINE_NAME'])
L.to_csv(os.path.join(S, 'gdsc_long.csv'), index=False)

# --- groups
ct = A.cmp_cancer_type
crc = ct == 'Colorectal Carcinoma'
mel = ct == 'Melanoma'
lung = A.cmp_tissue == 'Lung'
nsclc = ct.isin(['Non-Small Cell Lung Carcinoma', 'Squamous Cell Lung Carcinoma'])
panc = ct == 'Pancreatic Carcinoma'
brca = ct == 'Breast Carcinoma'
mk = A.mutation_data_available
B = lambda c: A[c].fillna(False).astype(bool)
GROUPS = {
    'a_BRAF_V600E_melanoma': (mel & B('BRAF_V600E'), 'cancer_type Melanoma, BRAF p.V600E'),
    'b_BRAF_V600E_colorectal': (crc & B('BRAF_V600E'), 'cancer_type Colorectal Carcinoma, BRAF p.V600E'),
    'c_KRAS_mut_colorectal': (crc & B('KRAS_mut'), 'Colorectal, KRAS missense at codon 12/13/59/61/117/146'),
    'd_KRAS_mut_lung_or_pancreas': ((lung | panc) & B('KRAS_mut'), 'tissue Lung (any lung cancer type) or Pancreatic Carcinoma, KRAS hotspot'),
    'd1_KRAS_mut_NSCLC': (nsclc & B('KRAS_mut'), 'NSCLC incl. squamous, KRAS hotspot'),
    'd2_KRAS_mut_pancreas': (panc & B('KRAS_mut'), 'Pancreatic Carcinoma, KRAS hotspot'),
    'e_RAS_RAF_wt_colorectal': (crc & mk & B('RAS_RAF_wt'), 'Colorectal with mutation data, no KRAS/NRAS/HRAS hotspot and no BRAF driver mutation'),
    'f_PIK3CA_mut_or_PTEN_null_colorectal_or_breast': ((crc | brca) & (B('PIK3CA_mut') | B('PTEN_null')), 'Colorectal or Breast; PIK3CA driver missense/inframe or PTEN truncating/homozygous deletion'),
    'f1_PIK3CA_mut_or_PTEN_null_breast': (brca & (B('PIK3CA_mut') | B('PTEN_null')), 'Breast subset of f'),
    'f2_PIK3CA_mut_or_PTEN_null_colorectal': (crc & (B('PIK3CA_mut') | B('PTEN_null')), 'Colorectal subset of f'),
    'f3_breast_PIK3CA_PTEN_wt': (brca & mk & ~(B('PIK3CA_mut') | B('PTEN_null')), 'comparator: Breast with mutation data, PIK3CA wt and PTEN not null'),
    'g_all_solid': (A.solid, 'all solid-tumour lines (CMP tissue not Haematopoietic and Lymphoid)'),
    'h_BRAF_wt_NRAS_mut_melanoma': (mel & B('NRAS_mut') & ~B('BRAF_V600_any'), 'extra: NRAS-mutant, BRAF V600-wt melanoma'),
    'i_EGFR_mut_NSCLC': (nsclc & B('EGFR_driver_mut'), 'extra: NSCLC with EGFR driver mutation (CMP cancer_driver flag)'),
    'j_RAS_RAF_wt_solid': (A.solid & mk & B('RAS_RAF_wt'), 'extra: all solid lines with no RAS hotspot / BRAF driver'),
}
QS = [5, 25, 50, 75, 95]

def summ(x):
    o = {'n_lines': int(len(x))}
    if len(x) == 0:
        return o
    o['ln_IC50_uM_quantiles'] = {f'p{q}': round(float(np.percentile(x.LN_IC50, q)), 3) for q in QS}
    o['AUC_quantiles'] = {f'p{q}': round(float(np.percentile(x.AUC, q)), 4) for q in QS}
    o['screened_conc_uM'] = {'min': float(x.MIN_CONC.min()), 'max': float(x.MAX_CONC.max()),
                              'n_distinct_max': int(x.MAX_CONC.nunique())}
    o['frac_IC50_below_max_conc'] = round(float(x.IC50_below_max.mean()), 3)
    o['median_RMSE'] = round(float(x.RMSE.median()), 4)
    return o

out = {
    'description': 'GDSC2 (release 8.5, GDSC2_fitted_dose_response_27Oct23) per drug x genotype group summaries. LN_IC50 = natural log of IC50 in micromolar (extrapolated beyond the screened range by the GDSC NLME fit, so values above ln(max_conc) are extrapolations). AUC = area under the fitted viability curve over the screened range, normalised to 1 (1 = no effect).',
    'mutation_source': 'Cell Model Passports mutations_all_latest (coding variants); PTEN homozygous deletion from cnv_summary_latest (cn_category == Deletion). Lines without CMP mutation data are excluded from genotype groups but included in g_all_solid.',
    'flag_definitions': {
        'BRAF_V600E': 'protein_mutation == p.V600E', 'KRAS_mut/NRAS_mut/HRAS_mut': 'missense at codon 12, 13, 59, 61, 117 or 146',
        'PIK3CA_mut': 'CMP cancer_driver == t and missense/inframe', 'PTEN_null': 'truncating (nonsense/frameshift/essential splice) or homozygous deletion',
        'RAS_RAF_wt': 'no KRAS/NRAS/HRAS hotspot and no BRAF variant flagged cancer_driver'},
    'drugs_not_in_GDSC2': NOT_IN_GDSC2,
    'group_definitions': {k: v[1] for k, v in GROUPS.items()},
    'group_line_counts': {k: int(v[0].sum()) for k, v in GROUPS.items()},
    'drugs': {},
}
for cls, Ld in DRUGS.items():
    for name, did, role in Ld:
        lab = name if name != 'Selumetinib' else f'Selumetinib_{did}'
        x = L[L.DRUG_ID == did]
        dd = {'drug_class': cls, 'drug_id': did, 'role': role, 'putative_target': x.PUTATIVE_TARGET.iloc[0], 'groups': {}}
        for g, (mask, _) in GROUPS.items():
            ids = set(A.index[mask])
            dd['groups'][g] = summ(x[x.SANGER_MODEL_ID.isin(ids)])
        out['drugs'][lab] = dd
json.dump(out, open(os.path.join(S, 'gdsc_targets.json'), 'w'), indent=1)

# --- famous lines
FAM = ['A375', 'SK-MEL-28', 'COLO-205', 'HT-29', 'RKO', 'WiDr', 'HCT-116', 'SW620', 'LoVo', 'DiFi', 'NCI-H1975', 'PC-9', 'PC-14']
norm = lambda s: re.sub(r'[^a-z0-9]', '', str(s).lower())
rows = []
for f in FAM:
    x = L[L.CELL_LINE_NAME.map(norm) == norm(f)]
    if not len(x):
        rows.append({'requested_line': f, 'note': 'not in GDSC2'})
        continue
    r0 = x.iloc[0]
    base = {'requested_line': f, 'gdsc_name': r0.CELL_LINE_NAME, 'COSMIC_ID': r0.COSMIC_ID, 'cancer_type': r0.cmp_cancer_type,
            'BRAF': r0.BRAF_driver_variants, 'KRAS': r0.KRAS_driver_variants, 'NRAS': r0.NRAS_driver_variants,
            'PIK3CA': r0.PIK3CA_driver_variants, 'PTEN_null': r0.PTEN_null, 'EGFR': r0.EGFR_driver_variants,
            'note': 'GDSC/CMP PC-14; Cellosaurus lists PC-14 as misidentified (actually PC-9)' if f == 'PC-14' else ''}
    order = [n if n != 'Selumetinib' else f'Selumetinib_{i}' for Ld in DRUGS.values() for n, i, _ in Ld]
    x = x.set_index('drug_label').reindex([o for o in order if o in set(x.drug_label)]).reset_index()
    for _, r in x.iterrows():
        base[f'lnIC50_{r.drug_label}'] = round(r.LN_IC50, 3)
    for _, r in x.iterrows():
        base[f'AUC_{r.drug_label}'] = round(r.AUC, 3)
    rows.append(base)
pd.DataFrame(rows).to_csv(os.path.join(S, 'famous_lines.csv'), index=False)

# doubling times: map NCI-60 names to GDSC names
dt = pd.read_csv(os.path.join(S, 'doubling_times.csv'))
gmap = {norm(n): (n, c, m) for m, n, c in zip(A.index, A.CELL_LINE_NAME, A.COSMIC_ID)}
alias = {'a549': 'a549', 'hs578t': 'hs578t', 'mdamb231': 'mdamb231', 'colo205': 'colo205', 'sw620': 'sw620', 'hct116': 'hct116', 'ht29': 'ht29'}
dt['gdsc_name'] = [gmap.get(norm(n), (None,))[0] for n in dt.cell_line]
dt['COSMIC_ID'] = [gmap.get(norm(n), (None, None))[1] for n in dt.cell_line]
dt['SANGER_MODEL_ID'] = [gmap.get(norm(n), (None, None, None))[2] for n in dt.cell_line]
dt.to_csv(os.path.join(S, 'doubling_times.csv'), index=False)
print('lines', len(A), 'with mut data', int(A.mutation_data_available.sum()))
print(pd.Series(out['group_line_counts']))
print('doubling mapped', dt.gdsc_name.notna().sum(), '/', len(dt))
