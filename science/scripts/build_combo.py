"""Build gdsc_combo_targets.json from the GDSC anchor combination screen (Jaaks et al., Nature 2022).
Usage: python3 -I scripts/build_combo.py <science_dir>
Inputs: raw/gdsc_combo/{breast,colon,pancreas}_anchor_combo.csv.gz ; cell_line_annotations.csv (from build_gdsc.py)
"""
import json, os, sys, itertools
import numpy as np
import pandas as pd

S = sys.argv[1]
CLASS = {
    'Dabrafenib': 'RAF', 'Trametinib': 'MEK', 'SCH772984': 'ERK',
    'Afatinib': 'EGFR', 'Sapitinib': 'EGFR', 'Erlotinib': 'EGFR', 'Lapatinib': 'EGFR', 'Gefitinib': 'EGFR',
    'Pictilisib': 'PI3K', 'Taselisib': 'PI3K', 'Alpelisib': 'PI3K', 'Dactolisib': 'PI3K', 'MK-2206': 'AKT',
    'Navitoclax': 'BCL2fam', 'Venetoclax': 'BCL2fam',
    'Afatinib | Trametinib': 'EGFR+MEK',
}
frames = []
for t in ['breast', 'colon', 'pancreas']:
    d = pd.read_csv(os.path.join(S, f'raw/gdsc_combo/{t}_anchor_combo.csv.gz'))
    d.columns = [c.strip() for c in d.columns]
    d['screen_tissue'] = t
    frames.append(d)
d = pd.concat(frames, ignore_index=True)
d = d.rename(columns={'Cell Line name': 'cell_line', 'SDIM': 'SIDM', 'Anchor Name': 'anchor', 'Anchor Conc': 'anchor_conc_uM', 'Library Name': 'library',
                      'Maxc': 'library_maxc_uM', 'Library IC50': 'library_xmid_log2norm', 'Combo IC50': 'combo_xmid_log2norm', 'Library Emax': 'library_emax',
                      'Bliss Emax': 'bliss_expected_emax', 'Combo Emax': 'combo_emax', 'Delta Xmid': 'delta_xmid_log2', 'Delta Emax': 'delta_emax', 'Synergy?': 'synergy'})
# check the published synergy rule (Jaaks 2022: dIC50 >= 3 log2 units, i.e. 8-fold, or dEmax >= 0.2)
rule = ((d.delta_xmid_log2 >= 3) | (d.delta_emax >= 0.2)).astype(int)
agree = float((rule == d.synergy).mean())
d['a_cls'] = d.anchor.map(CLASS)
d['l_cls'] = d.library.map(CLASS)
sub = d[d.a_cls.notna() & d.l_cls.notna()].copy()
sub['pair'] = [' + '.join(sorted([a, b])) for a, b in zip(sub.a_cls, sub.l_cls)]
sub['drug_pair'] = [' + '.join(sorted([a, b])) for a, b in zip(sub.anchor, sub.library)]
sub = sub[sub.anchor != sub.library]

ann = pd.read_csv(os.path.join(S, 'cell_line_annotations.csv'), index_col='model_id', low_memory=False)
for c in ['BRAF_V600E', 'KRAS_mut', 'NRAS_mut', 'PIK3CA_mut', 'PTEN_null', 'RAS_RAF_wt', 'mutation_data_available']:
    sub[c] = sub.SIDM.map(ann[c])
B = lambda c: sub[c].astype('boolean').fillna(False).astype(bool)
GEN = {
    'all': np.ones(len(sub), bool),
    'BRAF_V600E': B('BRAF_V600E'),
    'KRAS_mut': B('KRAS_mut'),
    'RAS_RAF_wt': B('RAS_RAF_wt') & B('mutation_data_available'),
    'PIK3CA_mut_or_PTEN_null': B('PIK3CA_mut') | B('PTEN_null'),
}

def per_line(x):
    # one row per cell line: synergistic at any anchor concentration / orientation; max effect sizes
    g = x.groupby('SIDM')
    return pd.DataFrame({'syn': g.synergy.max(), 'dEmax_max': g.delta_emax.max(), 'dXmid_max': g.delta_xmid_log2.max(),
                         'dEmax_med': g.delta_emax.median(), 'dXmid_med': g.delta_xmid_log2.median()})

def summ(x):
    if len(x) == 0:
        return None
    p = per_line(x)
    return {'n_lines': int(len(p)), 'n_rows': int(len(x)),
            'frac_lines_synergistic': round(float(p.syn.mean()), 3),
            'median_over_lines_of_max_delta_emax': round(float(p.dEmax_max.median()), 3),
            'median_over_lines_of_max_delta_xmid_log2': round(float(p.dXmid_max.median()), 3),
            'median_delta_emax_all_rows': round(float(x.delta_emax.median()), 3),
            'median_delta_xmid_log2_all_rows': round(float(x.delta_xmid_log2.median()), 3),
            'anchor_concs_uM': sorted(set(round(float(c), 4) for c in x.anchor_conc_uM)),
            'synergistic_lines': sorted(set(x[x.synergy == 1].cell_line))[:40]}

out = {'source': 'GDSC Combinations anchor screen (Jaaks et al., Nature 2022, doi:10.1038/s41586-022-04437-2), processed files {breast,colon,pancreas}_anchor_combo.csv.gz from https://gdsc-combinations.depmap.sanger.ac.uk/downloads',
       'definitions': {
           'design': '2 anchor concentrations x 7-point library dose response; each row = cell line x anchor x anchor conc x library drug',
           'delta_emax': 'Bliss-expected Emax minus observed combo Emax (viability units 0-1; positive = more kill than Bliss)',
           'delta_xmid_log2': 'library IC50 shift caused by anchor (library Xmid - combo Xmid) on log2 scale; 3 = 8-fold',
           'synergy': 'file flag; Jaaks 2022 rule is dXmid >= 3 or dEmax >= 0.2 at either anchor concentration. Agreement of file flag with this rule over all rows: %.4f' % agree,
           'per_line': 'a line counts as synergistic if any anchor conc / orientation row is flagged',
           'drug_classes': CLASS,
           'genotype_groups': 'from Cell Model Passports mutations (see gdsc_targets.json flag_definitions)'},
       'pairs': {}}
GEN = {k: pd.Series(v, index=sub.index) for k, v in GEN.items()}
for pair, x in sub.groupby('pair'):
    P = {'drug_pairs_tested': sorted(x.drug_pair.unique()), 'by_tissue': {}, 'by_drug_pair': {}}
    for t, xt in x.groupby('screen_tissue'):
        P['by_tissue'][t] = {g: summ(xt[m.loc[xt.index]]) for g, m in GEN.items()}
        P['by_tissue'][t] = {k: v for k, v in P['by_tissue'][t].items() if v}
    for dp, xd in x.groupby('drug_pair'):
        P['by_drug_pair'][dp] = {t: {g: summ(xt[m.loc[xt.index]]) for g, m in GEN.items() if m.loc[xt.index].any()} for t, xt in xd.groupby('screen_tissue')}
    out['pairs'][pair] = P
json.dump(out, open(os.path.join(S, 'gdsc_combo_targets.json'), 'w'), indent=1)
print('rule agreement', agree)
for pair, P in out['pairs'].items():
    for t, G in P['by_tissue'].items():
        print(pair.ljust(22), t.ljust(9), ' | '.join(f"{g}:{v['n_lines']}/{v['frac_lines_synergistic']:.2f}/{v['median_over_lines_of_max_delta_emax']:.2f}" for g, v in G.items()))
