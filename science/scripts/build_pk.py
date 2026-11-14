"""Write pk_table.json. Numbers are transcribed from FDA labels (openFDA JSON, raw/fda_labels),
EMA SmPCs (raw/ema) and papers (raw/pubmed); each value carries its source. Molar conversions computed here."""
import json, re, sys, os

S = sys.argv[1]
AW = {'C': 12.011, 'H': 1.008, 'N': 14.007, 'O': 15.999, 'S': 32.06, 'F': 18.998, 'Cl': 35.45, 'I': 126.904, 'Br': 79.904}

def mw(formula):
    return round(sum(AW[e] * (int(n) if n else 1) for e, n in re.findall(r'([A-Z][a-z]?)(\d*)', formula)), 2)

def um(ng_per_ml, mwt):
    return round(ng_per_ml / mwt, 4)  # ng/mL / (g/mol) = µmol/L

FDA = 'FDA label (openFDA, raw/fda_labels/{}.json), section {}'
EMA = 'EMA SmPC (raw/ema/{}.pdf), section {}'

drugs = []

# --- RAF: dabrafenib
f = 'C23H20F3N5O2S2'; m = mw(f)
drugs.append({
    'class': 'RAF_inhibitor', 'drug': 'dabrafenib', 'brand': 'Tafinlar',
    'molecular_formula_free_base': f, 'mw_free_base_g_per_mol': m, 'mw_note': 'label MW 615.7 is the mesylate salt; plasma concentrations refer to free base',
    'label_dose': {'value': '150 mg orally twice daily (q12h), fasted', 'source': FDA.format('TAFINLAR', '2.2')},
    'terminal_half_life_h': {'value': 8, 'note': 'oral, absorption-limited; IV microdose t1/2 2.6 h; hydroxy-dabrafenib 10 h; carboxy/desmethyl 21-22 h', 'source': FDA.format('TAFINLAR', '12.3') + '; ' + EMA.format('tafinlar', '5.2')},
    'Cmax_ss': {'ng_per_ml': 1478, 'uM': um(1478, m), 'stat': 'geometric mean, 150 mg bid', 'source': EMA.format('tafinlar', '5.2')},
    'Ctrough_ss': {'ng_per_ml': 26, 'uM': um(26, m), 'stat': 'geometric mean pre-dose (Cτ), 150 mg bid', 'source': EMA.format('tafinlar', '5.2')},
    'AUC_tau_ss_ng_h_per_ml': {'value': 4341, 'tau_h': 12, 'source': EMA.format('tafinlar', '5.2')},
    'Cavg_ss': {'ng_per_ml': round(4341 / 12, 1), 'uM': um(4341 / 12, m), 'source': 'derived: AUC(0-tau)/12 h'},
    'plasma_protein_binding_pct': {'value': 99.7, 'fraction_unbound': 0.003, 'source': FDA.format('TAFINLAR', '12.3')},
    'time_to_steady_state': {'value': '~14 days (auto-induction: CL/F rises from 17 L/h single dose to 34 L/h after 2 weeks bid; accumulation ratio 0.7)', 'source': FDA.format('TAFINLAR', '12.3')},
    'notes': 'Active metabolites hydroxy- and desmethyl-dabrafenib contribute (metabolite/parent AUC ratios 0.9, 11 (carboxy), 0.7).',
})

# --- RAF alternative: encorafenib (partial)
f = 'C22H27ClFN7O4S'; m = mw(f)
drugs.append({
    'class': 'RAF_inhibitor', 'drug': 'encorafenib', 'brand': 'Braftovi', 'molecular_formula_free_base': f, 'mw_free_base_g_per_mol': m,
    'label_dose': {'value': '450 mg once daily (melanoma, with binimetinib 45 mg bid); 300 mg once daily (BRAF V600E CRC, with cetuximab)', 'source': FDA.format('BRAFTOVI', '2')},
    'terminal_half_life_h': {'value': 3.5, 'note': 'mean (CV 17%); SmPC: median 6.32 h (3.74-8.09)', 'source': FDA.format('BRAFTOVI', '12.3') + '; ' + EMA.format('braftovi', '5.2')},
    'Cmax_ss': {'value': 'unverified', 'note': 'not stated in FDA label or EMA SmPC; not collected'},
    'plasma_protein_binding_pct': {'value': 86, 'fraction_unbound': 0.14, 'source': FDA.format('BRAFTOVI', '12.3')},
    'time_to_steady_state': {'value': 'within 15 days; steady-state exposure 50% lower than day 1 (CYP3A4 auto-induction)', 'source': FDA.format('BRAFTOVI', '12.3')},
})

# --- MEK: trametinib
f = 'C26H23FIN5O4'; m = mw(f)
drugs.append({
    'class': 'MEK_inhibitor', 'drug': 'trametinib', 'brand': 'Mekinist',
    'molecular_formula_free_base': f, 'mw_free_base_g_per_mol': m, 'mw_note': 'label molecular mass 693.5 is the DMSO solvate',
    'label_dose': {'value': '2 mg orally once daily, fasted', 'source': FDA.format('MEKINIST', '2')},
    'terminal_half_life_h': {'value': 127, 'note': 'SmPC mean terminal t1/2 after single dose 127 h (5.3 days); FDA label: estimated elimination t1/2 3.9-4.8 days; Infante 2012 effective t1/2 ~4 days', 'source': EMA.format('mekinist', '5.2') + '; ' + FDA.format('MEKINIST', '12.3')},
    'Cmax_ss': {'ng_per_ml': 22.2, 'uM': um(22.2, m), 'stat': 'geometric mean, 2 mg qd', 'source': EMA.format('mekinist', '5.2') + ' (same values in Infante et al., Lancet Oncol 2012, PMID 22805291)'},
    'Ctrough_ss': {'ng_per_ml': 12.1, 'uM': um(12.1, m), 'stat': 'geometric mean pre-dose', 'source': EMA.format('mekinist', '5.2')},
    'AUC_tau_ss_ng_h_per_ml': {'value': 370, 'tau_h': 24, 'source': EMA.format('mekinist', '5.2')},
    'Cavg_ss': {'ng_per_ml': round(370 / 24, 2), 'uM': um(370 / 24, m), 'source': 'derived: AUC(0-tau)/24 h'},
    'peak_trough_ratio': {'value': 1.8, 'source': EMA.format('mekinist', '5.2')},
    'plasma_protein_binding_pct': {'value': 97.4, 'fraction_unbound': 0.026, 'source': FDA.format('MEKINIST', '12.3')},
    'time_to_steady_state': {'value': 'by day 15 (accumulation ratio 6.0)', 'source': EMA.format('mekinist', '5.2')},
})

# --- EGFR small molecule: erlotinib
f = 'C22H23N3O4'; m = mw(f)
drugs.append({
    'class': 'EGFR_inhibitor', 'drug': 'erlotinib', 'brand': 'Tarceva (generic label used for FDA text)',
    'molecular_formula_free_base': f, 'mw_free_base_g_per_mol': m, 'mw_note': 'label MW 429.9 is the HCl salt',
    'label_dose': {'value': '150 mg orally once daily (NSCLC), on an empty stomach; 100 mg qd with gemcitabine (pancreatic)', 'source': FDA.format('ERLOTINIB', '2')},
    'terminal_half_life_h': {'value': 36.2, 'note': 'median, population PK (n=591); Hidalgo 2001 phase I: 24.4 ± 14.6 h', 'source': FDA.format('ERLOTINIB', '12.3') + '; Lu et al., Clin Pharmacol Ther 2006, PMID 16890575'},
    'Cmax_ss': {'value': 'unverified', 'note': 'steady-state Cmax not stated in label/SmPC abstracts; single-dose 150 mg Cmax 1056 ng/mL in non-smoking healthy volunteers (EMA SmPC 5.2)', 'single_dose_Cmax_ng_per_ml': 1056, 'single_dose_Cmax_uM': um(1056, m)},
    'Ctrough_ss': {'ng_per_ml': 1200, 'uM': um(1200, m), 'stat': 'mean ± SD 1.20 ± 0.62 µg/mL at 150 mg/d', 'source': 'Hidalgo et al., J Clin Oncol 2001, PMID 11432895 (abstract)',
                   'smoking': 'steady-state trough 1.28 µg/mL in former/never smokers vs 0.65 µg/mL in current smokers (EMA SmPC 5.2)'},
    'Cavg_ss': {'ng_per_ml': round(150000 / 4.47 / 24, 0), 'uM': um(150000 / 4.47 / 24, m), 'source': 'derived: dose/(CL/F x 24 h), CL/F 4.47 L/h (EMA SmPC 5.2)'},
    'plasma_protein_binding_pct': {'value': 93, 'fraction_unbound': 0.07, 'note': 'SmPC says ~95%', 'source': FDA.format('ERLOTINIB', '12.3')},
    'time_to_steady_state': {'value': '7-8 days', 'source': FDA.format('ERLOTINIB', '12.3')},
})

# --- EGFR antibody: cetuximab
m = 152000
drugs.append({
    'class': 'EGFR_inhibitor', 'drug': 'cetuximab', 'brand': 'Erbitux', 'mw_g_per_mol': m, 'mw_note': 'approximately 152 kDa (label)',
    'label_dose': {'value': '400 mg/m2 IV loading dose then 250 mg/m2 weekly, or 500 mg/m2 every 2 weeks', 'source': FDA.format('ERBITUX', '2')},
    'terminal_half_life_h': {'value': 112, 'range': '63-230 h', 'source': FDA.format('ERBITUX', '12.3')},
    'Cmax_ss': {'ug_per_ml_range_of_means': [168, 235], 'uM_range': [round(168e3 / m, 3), round(235e3 / m, 3)], 'stat': 'mean peak across studies, 400/250 mg/m2 weekly', 'source': FDA.format('ERBITUX', '12.3')},
    'Ctrough_ss': {'ug_per_ml_range_of_means': [41, 85], 'uM_range': [round(41e3 / m, 3), round(85e3 / m, 3)], 'source': FDA.format('ERBITUX', '12.3')},
    'plasma_protein_binding_pct': {'value': 'not applicable (IgG1 antibody)'},
    'time_to_steady_state': {'value': 'by the third weekly infusion', 'source': FDA.format('ERBITUX', '12.3')},
    'notes': 'Nonlinear (target-mediated) clearance: CL fell from 0.08 to 0.02 L/h/m2 as dose rose from 20 to 200 mg/m2, plateau above 200 mg/m2.',
})

# --- PI3K: alpelisib
f = 'C19H22F3N5O2S'; m = mw(f)
drugs.append({
    'class': 'PI3K_inhibitor', 'drug': 'alpelisib', 'brand': 'Piqray', 'molecular_formula_free_base': f, 'mw_free_base_g_per_mol': m,
    'label_dose': {'value': '300 mg orally once daily with food', 'source': FDA.format('PIQRAY', '2')},
    'terminal_half_life_h': {'value': '8-9', 'source': FDA.format('PIQRAY', '12.3')},
    'Cmax_ss': {'ng_per_ml': 2480, 'uM': um(2480, m), 'stat': 'mean (CV 23%), SOLAR-1, 300 mg qd', 'source': FDA.format('PIQRAY', '12.3')},
    'AUC_tau_ss_ng_h_per_ml': {'value': 33224, 'tau_h': 24, 'source': FDA.format('PIQRAY', '12.3')},
    'Cavg_ss': {'ng_per_ml': round(33224 / 24, 0), 'uM': um(33224 / 24, m), 'source': 'derived: AUC/24 h'},
    'Ctrough_ss': {'value': 'unverified', 'note': 'not stated in label; one-compartment estimate with t1/2 8.5 h and Tmax ~3 h gives ~0.25 x Cmax ≈ 600 ng/mL (rough derivation, not a measured value)'},
    'plasma_protein_binding_pct': {'value': 89, 'fraction_unbound': 0.11, 'source': FDA.format('PIQRAY', '12.3')},
    'time_to_steady_state': {'value': 'within 3 days (accumulation 1.3-1.5)', 'source': FDA.format('PIQRAY', '12.3')},
})

# --- BCL-2 family: venetoclax and navitoclax
f = 'C45H50ClN7O7S'; m = mw(f)
drugs.append({
    'class': 'BCL2_family_inhibitor', 'drug': 'venetoclax', 'brand': 'Venclexta', 'molecular_formula_free_base': f, 'mw_free_base_g_per_mol': m,
    'label_dose': {'value': '400 mg orally once daily with a meal after a 5-week ramp-up (20/50/100/200/400 mg), CLL/SLL; AML dosing differs (label 2.3)', 'source': FDA.format('VENCLEXTA', '2')},
    'terminal_half_life_h': {'value': 26, 'source': FDA.format('VENCLEXTA', '12.3')},
    'Cmax_ss': {'ng_per_ml': 2100, 'uM': um(2100, m), 'stat': 'mean ± SD 2.1 ± 1.1 µg/mL, 400 mg qd, low-fat meal', 'source': FDA.format('VENCLEXTA', '12.3')},
    'AUC_tau_ss_ng_h_per_ml': {'value': 32800, 'tau_h': 24, 'source': FDA.format('VENCLEXTA', '12.3')},
    'Cavg_ss': {'ng_per_ml': round(32800 / 24, 0), 'uM': um(32800 / 24, m), 'source': 'derived: AUC(0-24)/24 h'},
    'Ctrough_ss': {'value': 'unverified', 'note': 'not stated in label'},
    'plasma_protein_binding_pct': {'value': '>99', 'fraction_unbound': '<0.01', 'source': FDA.format('VENCLEXTA', '12.3')},
    'time_to_steady_state': {'value': 'unverified; ~5 days expected from 26 h half-life (5 half-lives), label says PK does not change over time'},
    'notes': 'Selective BCL-2; solid-tumour relevance limited. Food raises exposure 3.4-5.3x.',
})
f = 'C47H55ClF3N5O6S3'; m = mw(f)
drugs.append({
    'class': 'BCL2_family_inhibitor', 'drug': 'navitoclax (ABT-263; not approved)', 'molecular_formula_free_base': f, 'mw_free_base_g_per_mol': m,
    'label_dose': {'value': 'no label. Phase I/II solid tumour + trametinib RP2D: navitoclax 250 mg daily d1-28 with trametinib 2 mg d1-14 (Corcoran et al., Clin Cancer Res 2024, PMID 38456660); solid-tumour monotherapy phase I used a 1-week 150 mg lead-in then continuous daily dosing (Gandhi et al., J Clin Oncol 2011, PMID 21282543, abstract); the 325 mg daily phase II dose is unverified'},
    'terminal_half_life_h': {'value': 17.1, 'sd': 4.3, 'source': 'Wilson et al., Lancet Oncol 2010, PMID 21094089, Table 3 (all doses, day -3)'},
    'Cmax_ss': {'ug_per_ml': 6.44, 'sd': 3.18, 'uM': um(6440, m), 'stat': 'mean, day 14, 315 mg daily, n=8 (lymphoid malignancies)', 'source': 'Wilson et al., Lancet Oncol 2010, PMID 21094089, Table 3 (PMC3025495)'},
    'AUC_0_24_day1_ug_h_per_ml_315mg': {'value': 58.5, 'source': 'Wilson 2010 Table 3'},
    'Cavg_day1_315mg': {'ug_per_ml': round(58.5 / 24, 2), 'uM': um(58500 / 24, m), 'source': 'derived: AUC0-24 day 1 / 24 h (underestimates steady state; day-14 Cmax was ~1.5x day 1)'},
    'Tmax_h': {'value': '~9', 'source': 'Wilson 2010'},
    'plasma_protein_binding_pct': {'value': 'unverified (>99.9% commonly cited)'},
    'time_to_steady_state': {'value': 'unverified; ~3-4 days expected from 17 h half-life'},
    'notes': 'Dose-limiting on-target thrombocytopenia (BCL-xL in platelets).',
})

out = {'description': 'Representative clinical PK at label dose, one or more drugs per class. uM values use free-base MW computed from the molecular formula. Total (bound + unbound) plasma concentrations unless stated; multiply by fraction_unbound for free drug.',
       'retrieved': '2026-11-14', 'drugs': drugs}
json.dump(out, open(os.path.join(S, 'pk_table.json'), 'w'), indent=1, ensure_ascii=False)
for d in drugs:
    print(d['drug'], d.get('mw_free_base_g_per_mol', d.get('mw_g_per_mol')), d.get('Cmax_ss', {}).get('uM'), d.get('Ctrough_ss', {}).get('uM'), d.get('Cavg_ss', {}).get('uM'))
