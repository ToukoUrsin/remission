"""Write clinical_trials.json. Numbers transcribed from abstracts / full texts saved in raw/pubmed and FDA labels in raw/fda_labels.
Every number not checked against the abstract, full text or label text is the string "unverified" (or carries verified=False).
Usage: python3 -I scripts/build_clinical.py <science_dir>"""
import json, os, sys

S = sys.argv[1]
U = 'unverified'

def ref(authors, journal, year, doi, pmid, pmcid=None):
    r = {'authors': authors, 'journal': journal, 'year': year, 'doi': doi, 'pmid': pmid}
    if pmcid:
        r['pmcid'] = pmcid
    return r

T = []

T.append({'id': 'BRIM-3', 'regimen': 'vemurafenib monotherapy', 'indication': 'BRAF V600E metastatic melanoma, first line', 'phase': 3,
          'dose_schedule': 'vemurafenib 960 mg orally twice daily (comparator dacarbazine 1000 mg/m2 IV q3w)',
          'arms': {
              'vemurafenib': {'n_randomized': 'unverified (675 total)', 'ORR': {'value_pct': 48, 'responders': 106, 'n_evaluable': 219, 'ci95_pct': [42, 55]},
                              'median_PFS_months': 5.3, 'grade_ge3_AE_pct': U,
                              'other_toxicity': '38% (129/336) needed dose modification for toxicity; cutaneous SCC or keratoacanthoma in 18%'},
              'dacarbazine': {'n_randomized': U, 'ORR': {'value_pct': 5, 'ci95_pct': U}, 'median_PFS_months': 1.6, 'grade_ge3_AE_pct': U}},
          'citation': ref('Chapman PB, Hauschild A, Robert C, et al.', 'N Engl J Med', 2011, '10.1056/NEJMoa1103782', 21639808, 'PMC3549296'),
          'notes': 'Interim analysis; PFS assessed in 549 patients. Overall grade >=3 AE rate not stated in the abstract/text checked.'})

T.append({'id': 'BREAK-3', 'regimen': 'dabrafenib monotherapy', 'indication': 'BRAF V600E metastatic melanoma, first line', 'phase': 3,
          'dose_schedule': 'dabrafenib 150 mg orally twice daily (comparator dacarbazine 1000 mg/m2 IV q3w)',
          'arms': {
              'dabrafenib': {'n_randomized': 187, 'ORR': {'value_pct': 52, 'n': 187, 'ci95_pct': [44, 59], 'source': 'FDA Tafinlar label section 14 (independent review)'},
                             'median_PFS_months': 5.1, 'grade_ge3_AE_pct': U,
                             'other_toxicity': 'abstract: grade 3-4 AEs uncommon; treatment-related grade >=2 AEs in 53% of dabrafenib patients'},
              'dacarbazine': {'n_randomized': 63, 'ORR': {'value_pct': 17, 'source': 'FDA Tafinlar label'}, 'median_PFS_months': 2.7, 'grade_ge3_AE_pct': U}},
          'citation': ref('Hauschild A, Grob JJ, Demidov LV, et al.', 'Lancet', 2012, '10.1016/S0140-6736(12)60868-X', 22735384),
          'notes': 'Hazard ratio for PFS 0.30.'})

T.append({'id': 'COMBI-d', 'regimen': 'dabrafenib + trametinib vs dabrafenib + placebo', 'indication': 'BRAF V600E/K metastatic melanoma, first line', 'phase': 3,
          'dose_schedule': 'dabrafenib 150 mg orally twice daily + trametinib 2 mg orally once daily',
          'arms': {
              'dabrafenib+trametinib': {'n_randomized': 211, 'ORR': {'value_pct': 67, 'label_value_pct': 66, 'label_ci95_pct': [60, 73]},
                                        'median_PFS_months': 9.3, 'median_PFS_months_updated_2015': 11.0, 'grade_ge3_AE_pct': 32,
                                        'grade_ge3_note': 'grade 3-4 AEs 32% (combination) vs 31% (monotherapy); pyrexia grade 3 in 6%'},
              'dabrafenib+placebo': {'n_randomized': 212, 'ORR': {'value_pct': 51}, 'median_PFS_months': 8.8, 'median_PFS_months_updated_2015': 8.8, 'grade_ge3_AE_pct': 31}},
          'citation': [ref('Long GV, Stroyakovskiy D, Gogas H, et al.', 'N Engl J Med', 2014, '10.1056/NEJMoa1406037', 25265492),
                       ref('Long GV, Stroyakovskiy D, Gogas H, et al. (updated OS analysis)', 'Lancet', 2015, '10.1016/S0140-6736(15)60898-4', 26037941)]})

T.append({'id': 'COMBI-v', 'regimen': 'dabrafenib + trametinib vs vemurafenib', 'indication': 'BRAF V600E/K metastatic melanoma, first line', 'phase': 3,
          'dose_schedule': 'dabrafenib 150 mg twice daily + trametinib 2 mg once daily vs vemurafenib 960 mg twice daily',
          'arms': {
              'dabrafenib+trametinib': {'n_randomized': U, 'ORR': {'value_pct': 64, 'ci95_pct': U}, 'median_PFS_months': 11.4, 'grade_ge3_AE_pct': U},
              'vemurafenib': {'n_randomized': U, 'ORR': {'value_pct': 51, 'ci95_pct': U}, 'median_PFS_months': 7.3, 'grade_ge3_AE_pct': U}},
          'citation': ref('Robert C, Karaszewska B, Schachter J, et al.', 'N Engl J Med', 2015, '10.1056/NEJMoa1412690', 25399551),
          'notes': 'n = 704 total (1:1; per-arm n not checked). Cutaneous SCC/keratoacanthoma 1% vs 18%.'})

T.append({'id': 'Kopetz-2015-vemurafenib-CRC', 'regimen': 'vemurafenib monotherapy', 'indication': 'BRAF V600E metastatic colorectal cancer, pretreated', 'phase': 2,
          'dose_schedule': 'vemurafenib 960 mg orally twice daily',
          'arms': {'vemurafenib': {'n': 21, 'ORR': {'value_pct': 5, 'responders': 1, 'n': 21, 'ci95_pct': [1, 24]}, 'median_PFS_months': 2.1, 'median_OS_months': 7.7,
                                   'grade_ge3_AE_pct': round(100 * 15 / 21, 1), 'grade_ge3_note': '15/21 patients had >=1 grade 3 AE; no grade 4 (Table 2)'}},
          'citation': ref('Kopetz S, Desai J, Chan E, et al.', 'J Clin Oncol', 2015, '10.1200/JCO.2015.63.2497', 26460303, 'PMC4669589')})

T.append({'id': 'BEACON-CRC', 'regimen': 'encorafenib + cetuximab (doublet); encorafenib + binimetinib + cetuximab (triplet); control = investigator choice of irinotecan or FOLFIRI + cetuximab',
          'indication': 'BRAF V600E metastatic colorectal cancer after 1-2 prior regimens', 'phase': 3,
          'dose_schedule': 'encorafenib 300 mg once daily; binimetinib 45 mg twice daily (triplet only); cetuximab 400 mg/m2 then 250 mg/m2 weekly',
          'primary_analysis_2019': {
              'triplet': {'ORR': {'value_pct': 26, 'ci95_pct': [18, 35], 'note': 'first 331 randomized patients'}, 'median_OS_months': 9.0, 'grade_ge3_AE_pct': 58},
              'doublet': {'ORR': {'value_pct': 20, 'ci95_pct': [13, 29], 'n': 113, 'source': 'FDA Braftovi label section 14'}, 'median_PFS_months': 4.2,
                          'median_PFS_source': 'FDA Braftovi label (4.2 vs 1.5 months)', 'median_OS_months': 8.4, 'grade_ge3_AE_pct': 50},
              'control': {'ORR': {'value_pct': 2}, 'median_PFS_months': 1.5, 'median_OS_months': 5.4, 'grade_ge3_AE_pct': 61},
              'citation': ref('Kopetz S, Grothey A, Yaeger R, et al.', 'N Engl J Med', 2019, '10.1056/NEJMoa1908075', 31566309)},
          'updated_analysis_2021': {
              'triplet': {'n': 224, 'ORR': {'value_pct': 26.8, 'ci95_pct': [21.1, 33.1]}, 'median_PFS_months': 4.5, 'PFS_ci95': [4.2, 5.4], 'median_OS_months': 9.3, 'grade_ge3_AE_pct': 65.8},
              'doublet': {'n': 220, 'ORR': {'value_pct': 19.5, 'ci95_pct': [14.5, 25.4]}, 'median_PFS_months': 4.3, 'PFS_ci95': [4.1, 5.4], 'median_OS_months': 9.3, 'grade_ge3_AE_pct': 57.4},
              'control': {'n': 221, 'ORR': {'value_pct': 1.8, 'ci95_pct': [0.5, 4.6]}, 'median_PFS_months': 1.5, 'median_OS_months': 5.9, 'grade_ge3_AE_pct': 64.2},
              'citation': ref('Tabernero J, Grothey A, Van Cutsem E, et al.', 'J Clin Oncol', 2021, '10.1200/JCO.20.02088', 33503393, 'PMC8078423')}})

T.append({'id': 'Corcoran-2015-DT-CRC', 'regimen': 'dabrafenib + trametinib', 'indication': 'BRAF V600E metastatic colorectal cancer', 'phase': '1/2',
          'dose_schedule': 'dabrafenib 150 mg twice daily + trametinib 2 mg once daily',
          'arms': {'dabrafenib+trametinib': {'n': 43, 'ORR': {'value_pct': 12, 'responders': 5, 'n': 43, 'ci95_pct': U, 'note': '1 CR; stable disease 56%. Corcoran 2018 Table 3 lists confirmed ORR 3/43 (7%)'},
                                             'median_PFS_months': 3.5, 'median_PFS_source': 'Corcoran 2018 Cancer Discov comparison table', 'grade_ge3_AE_pct': U,
                                             'pharmacodynamics': 'on-treatment pERK reduction 47 +/- 24% in biopsies (less than in melanoma)'}},
          'citation': ref('Corcoran RB, Atreya CE, Falchook GS, et al.', 'J Clin Oncol', 2015, '10.1200/JCO.2015.63.2471', 26392102)})

T.append({'id': 'Corcoran-2018-DTP-CRC', 'regimen': 'dabrafenib + panitumumab (D+P), trametinib + panitumumab (T+P), dabrafenib + trametinib + panitumumab (D+T+P)',
          'indication': 'BRAF V600E metastatic colorectal cancer', 'phase': '1/2',
          'dose_schedule': 'dabrafenib 150 mg twice daily; trametinib 2 mg once daily; panitumumab 6 mg/kg IV every 2 weeks (full doses at RP2D)',
          'arms': {
              'D+T+P': {'n': 91, 'ORR': {'value_pct': 21, 'ci95_pct': [13.1, 30.7]}, 'median_PFS_months': 4.2, 'PFS_ci95': [4.0, 5.6], 'median_OS_months': 9.1, 'grade_ge3_AE_pct': 70, 'grade_ge3_note': 'grade 3 or 4'},
              'D+P': {'n': 20, 'ORR': {'value_pct': 10, 'ci95_pct': [1.2, 31.7]}, 'median_PFS_months': 3.5, 'PFS_ci95': [2.8, 5.8], 'median_OS_months': 13.2, 'grade_ge3_AE_pct': 45},
              'T+P': {'n': 31, 'ORR': {'value_pct': 0, 'ci95_pct': [0, 11.2]}, 'median_PFS_months': 2.6, 'grade_ge3_AE_pct': U, 'grade_ge3_note': 'grade 3/4 dermatitis acneiform 18%'}},
          'citation': ref('Corcoran RB, André T, Atreya CE, et al.', 'Cancer Discov', 2018, '10.1158/2159-8290.CD-17-1226', 29431699, 'PMC5882509')})

T.append({'id': 'CO.17-KRAS', 'regimen': 'cetuximab monotherapy vs best supportive care', 'indication': 'chemorefractory metastatic colorectal cancer, by tumour KRAS status', 'phase': 3,
          'dose_schedule': 'cetuximab 400 mg/m2 IV loading then 250 mg/m2 weekly (standard label schedule; not re-checked in paper text: unverified)',
          'arms': {
              'KRAS_wild_type_cetuximab': {'n': U, 'ORR': {'value_pct': 12.8}, 'median_PFS_months': 3.7, 'median_OS_months': 9.5, 'PFS_HR_vs_BSC': 0.40, 'grade_ge3_AE_pct': U},
              'KRAS_wild_type_BSC': {'median_PFS_months': 1.9, 'median_OS_months': 4.8},
              'KRAS_mutant_cetuximab': {'ORR': {'value_pct': 1.2}, 'PFS_HR_vs_BSC': 0.99, 'grade_ge3_AE_pct': U}},
          'citation': ref('Karapetis CS, Khambata-Ford S, Jonker DJ, et al.', 'N Engl J Med', 2008, '10.1056/NEJMoa0804385', 18946061),
          'notes': 'ORR 12.8% vs 1.2% verified from NEJM article text via web search; arm sizes not checked.'})

T.append({'id': 'Amado-2008-panitumumab-KRAS', 'regimen': 'panitumumab monotherapy vs best supportive care', 'indication': 'chemorefractory metastatic colorectal cancer, by tumour KRAS status', 'phase': 3,
          'dose_schedule': 'panitumumab 6 mg/kg IV every 2 weeks (standard; not re-checked in text: unverified)',
          'arms': {
              'KRAS_wild_type_panitumumab': {'ORR': {'value_pct': 17}, 'median_PFS_weeks': 12.3, 'PFS_HR_vs_BSC': 0.45, 'grade_ge3_AE_pct': U},
              'KRAS_wild_type_BSC': {'median_PFS_weeks': 7.3},
              'KRAS_mutant_panitumumab': {'ORR': {'value_pct': 0}, 'PFS_HR_vs_BSC': 0.99, 'grade_ge3_AE_pct': U}},
          'citation': ref('Amado RG, Wolf M, Peeters M, et al.', 'J Clin Oncol', 2008, '10.1200/JCO.2007.14.7116', 18316791),
          'notes': 'KRAS status determined in 427/463 patients (92%); 43% mutant.'})

T.append({'id': 'MEKi-monotherapy-KRAS', 'regimen': 'MEK inhibitor monotherapy in KRAS-mutant tumours',
          'trials': [
              {'drug': 'trametinib vs docetaxel', 'indication': 'KRAS-mutant NSCLC, second line', 'phase': 2, 'dose_schedule': 'trametinib 2 mg once daily; docetaxel 75 mg/m2 q3w',
               'n_trametinib': 86, 'ORR': {'value_pct': 12, 'note': 'partial responses'}, 'median_PFS_weeks': 12, 'comparator_PFS_weeks': 11, 'PFS_HR': 1.14, 'grade_ge3_AE_pct': U,
               'citation': ref('Blumenschein GR Jr, Smit EF, Planchard D, et al.', 'Ann Oncol', 2015, '10.1093/annonc/mdv072', 25722381)},
              {'drug': 'selumetinib vs capecitabine', 'indication': 'metastatic colorectal cancer after 1-2 prior regimens (not KRAS-selected)', 'phase': 2, 'dose_schedule': 'selumetinib 100 mg twice daily',
               'n_selumetinib': 34, 'ORR': {'value_pct': 0, 'note': '0 responses; 10 stable disease'}, 'median_PFS_days': 81, 'comparator_PFS_days': 88, 'grade_ge3_AE_pct': U,
               'citation': ref('Bennouna J, Lang I, Valladares-Ayerbes M, et al.', 'Invest New Drugs', 2011, '10.1007/s10637-010-9392-8', 20127139)},
              {'drug': 'selumetinib vs capecitabine', 'indication': 'gemcitabine-refractory pancreatic cancer', 'phase': 2, 'dose_schedule': 'selumetinib 100 mg twice daily',
               'median_OS_months': 5.4, 'comparator_OS_months': 5.0, 'ORR': U, 'grade_ge3_AE_pct': U,
               'citation': ref('Bodoky G, Timcheva C, Spigel DR, et al.', 'Invest New Drugs', 2012, '10.1007/s10637-011-9687-4', 21594619)},
              {'drug': 'selumetinib vs pemetrexed', 'indication': 'advanced NSCLC, second/third line', 'phase': 2, 'dose_schedule': 'selumetinib 100 mg twice daily',
               'ORR': {'note': '2 partial responses on selumetinib (n unverified)'}, 'grade_ge3_AE_pct': U,
               'citation': ref('Hainsworth JD, Cebotaru CL, Kanarev V, et al.', 'J Thorac Oncol', 2010, '10.1097/JTO.0b013e3181e8b3a3', 20802351)},
              {'drug': 'trametinib phase I (all comers)', 'phase': 1, 'dose_schedule': '0.125-10 mg; RP2D 2 mg once daily', 'ORR': {'value_pct': 10, 'responders': 21, 'n': 206},
               'citation': ref('Infante JR, Fecher LA, Falchook GS, et al.', 'Lancet Oncol', 2012, '10.1016/S1470-2045(12)70269-3', 22805291)}],
          'notes': 'Jänne et al. Lancet Oncol 2013 (PMID 23200175, doi 10.1016/S1470-2045(12)70489-8) tested selumetinib 75 mg bid + docetaxel, a combination, so it is not monotherapy.'})

T.append({'id': 'MEKi+PI3K-pathway-toxicity', 'regimen': 'MEK inhibitor + PI3K/AKT/mTOR pathway inhibitor combinations',
          'trials': [
              {'study': 'Shimizu 2012: pooled phase I patients receiving PI3K-AKT-mTOR + RAS-RAF-MEK pathway inhibitor combinations (dual) vs single-pathway agents',
               'dose_schedule': 'various investigational combinations (not one regimen)',
               'drug_related_grade_ge3_AE_pct': {'dual_pathway': 53.9, 'single_pathway': 18.1}, 'DLT_rate_pct': {'dual_pathway': 18.4, 'single_pathway': 9.4},
               'efficacy': '5/5 patients with co-mutations treated on dual-pathway therapy had tumour shrinkage of 2-64%; overall ORR unverified',
               'citation': ref('Shimizu T, Tolcher AW, Papadopoulos KP, et al.', 'Clin Cancer Res', 2012, '10.1158/1078-0432.CCR-11-2381', 22261800)},
              {'study': 'Bedard 2015: buparlisib (pan-PI3K) + trametinib phase Ib, RAS/BRAF-mutant solid tumours', 'n': 113,
               'dose_schedule': 'MTD buparlisib 70 mg qd + trametinib 1.5 mg qd (DLT 5/15 = 33%); RP2D buparlisib 60 mg qd + trametinib 1.5 mg qd (DLT 1/10)',
               'treatment_related_grade_3_4_AE_pct': 65,
               'ORR': {'KRAS_mutant_ovarian': {'value_pct': 29, 'n': 21}, 'NSCLC': {'responders': 1, 'n': 17}, 'pancreatic': 'best response stable disease'},
               'citation': ref('Bedard PL, Tabernero J, Janku F, et al.', 'Clin Cancer Res', 2015, '10.1158/1078-0432.CCR-14-1814', 25500057)}]})

T.append({'id': 'Navitoclax+trametinib', 'regimen': 'navitoclax + trametinib', 'indication': 'KRAS- or NRAS-mutant advanced solid tumours', 'phase': '1/2',
          'dose_schedule': 'RP2D: navitoclax 150 mg daily lead-in days 1-7 of cycle 1, then 250 mg daily days 1-28; trametinib 2 mg daily days 1-14 of each 28-day cycle',
          'n_total': 91, 'n_dose_escalation': 38,
          'ORR': {'at_RP2D': {'value_pct': 16, 'responders': 8, 'n_evaluable': 49}, 'gynaecological': {'value_pct': 33, 'responders': 7, 'n': 21},
                  'colorectal_NSCLC_pancreatic': {'responders': 0, 'note': 'no partial responses'}},
          'median_PFS_months': {'gynaecological': 4.8, 'other_tumour_types': 1.87},
          'grade_ge3_AE_pct': U, 'grade_ge3_note': 'grade >=3 thrombocytopenia 2%; DLTs included grade 3 AST/ALT rise, grade 4 hypokalaemia, grade 3 hyperbilirubinaemia, grade 3 enterocolitis',
          'citation': ref('Corcoran RB, Do KT, Kim JE, et al.', 'Clin Cancer Res', 2024, '10.1158/1078-0432.CCR-23-3135', 38456660, 'PMC11061595')})

out = {'description': 'Published clinical outcomes for targeted-therapy regimens used to calibrate the Remission simulator. ORR = objective response rate (RECIST), PFS = progression-free survival. '
                      '"unverified" marks numbers not confirmed against the paper abstract, full text or FDA label text during collection.',
       'retrieved': '2026-11-14', 'trials': T}
json.dump(out, open(os.path.join(S, 'clinical_trials.json'), 'w'), indent=1, ensure_ascii=False)
print(len(T), 'entries')
