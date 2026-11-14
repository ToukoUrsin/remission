"""Write SOURCES.md: every file in raw/ with URL, date, size, licence; every paper with DOI.
Usage: python3 -I scripts/build_sources.py <science_dir>"""
import os, sys

S = sys.argv[1]
R = os.path.join(S, 'raw')
DATE = '2026-11-14'
GDSC = 'https://cog.sanger.ac.uk/cancerrxgene/GDSC_release8.5/'
CMP = 'https://cog.sanger.ac.uk/cmp/download/'
COMBO = 'https://gdsc-combinations.depmap.sanger.ac.uk'
EMA = 'https://www.ema.europa.eu/en/documents/product-information/{}-epar-product-information_en.pdf'
OFDA = 'https://api.fda.gov/drug/label.json?search=openfda.brand_name:"{}"'
EFETCH = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pmc&id={}'
LIC_SANGER = 'Wellcome Sanger Institute open data, free reuse with citation (exact licence not checked: unverified)'
LIC_PMC = 'PMC copy; publisher / author-manuscript terms (research use)'

def src(rel):
    d, f = rel.split('/', 1)
    if d == 'gdsc':
        if f.endswith('.csv'):
            return 'derived locally from the .xlsx (pandas)', LIC_SANGER
        return GDSC + f, LIC_SANGER
    if d == 'cmp':
        if f == 'mutations_selected_genes.csv':
            return 'derived locally: rows of mutations_all_latest.csv.gz for BRAF/KRAS/NRAS/HRAS/PIK3CA/PTEN/EGFR', LIC_SANGER
        return CMP + f, LIC_SANGER
    if d == 'gdsc_combo':
        if f.endswith('_anchor_combo.csv.gz'):
            return f'{COMBO}/downloads/{f.split("_")[0]}/anchor_combo', LIC_SANGER + '; paper CC BY 4.0'
        if f.endswith('_matrix_results.csv.gz'):
            return f'{COMBO}/downloads/{f.replace("_matrix_results.csv.gz", "")}/matrix_results', LIC_SANGER
        if f == 'anchor_glossary.pdf':
            return f'{COMBO}/downloads/glossary.pdf', LIC_SANGER
        if f == 'Matrix_glossary.pdf':
            return f'{COMBO} (matrix screen glossary link on downloads page)', LIC_SANGER
        return f'{COMBO} (portal page / Dash layout callbacks, used to find download links)', LIC_SANGER
    if d == 'ema':
        n = f.rsplit('.', 1)[0]
        if f.endswith('.txt'):
            return f'derived locally: pdftotext of {n}.pdf', 'EMA: reproduction authorised with source acknowledged'
        return EMA.format(n), 'EMA: reproduction authorised with source acknowledged'
    if d == 'fda_labels':
        n = f.rsplit('.', 1)[0]
        q = OFDA.format(n) if n != 'ERLOTINIB' else 'https://api.fda.gov/drug/label.json?search=openfda.generic_name:"erlotinib"'
        return q, 'openFDA: public domain / CC0'
    if d == 'nci60':
        return 'https://www.ebi.ac.uk/europepmc/webservices/rest/PMC5654365/fullTextXML', 'CC BY 4.0'
    if d == 'pubmed':
        if f.startswith('PMC'):
            pid = f[3:].split('.')[0]
            lic = {'3964313': 'CC BY-NC-SA 3.0', '4669589': 'CC BY-NC-ND 3.0 US', '6700527': 'CC BY 4.0', '8891012': 'CC BY 4.0'}.get(pid, LIC_PMC)
            via = EFETCH.format(pid) if pid not in ('4668282', '3025495') else f'https://www.ncbi.nlm.nih.gov/research/bionlp/RESTful/pmcoa.cgi/BioC_xml/PMC{pid}/unicode'
            return via, lic
        if f == 'pmc_navito.xml':
            return EFETCH.format('4668282,3025495'), LIC_PMC
        if f == 'esum1.json':
            return 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=pubmed&retmode=json&id=<PMIDs>', 'NLM PubMed metadata (public)'
        return 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&rettype=abstract&id=<PMIDs>', 'PubMed abstracts (publisher copyright; quoted for research)'
    return '?', '?'

rows = []
for d in sorted(os.listdir(R)):
    for f in sorted(os.listdir(os.path.join(R, d))):
        rel = f'{d}/{f}'
        url, lic = src(rel)
        rows.append(f'| raw/{rel} | {url} | {DATE} | {os.path.getsize(os.path.join(R, rel)):,} | {lic} |')

PAPERS = [
    ('Iorio F et al. A landscape of pharmacogenomic interactions in cancer. Cell 2016', '10.1016/j.cell.2016.06.017', 'GDSC'),
    ('van der Meer D et al. Cell Model Passports. Nucleic Acids Res 2019', '10.1093/nar/gky872', 'CMP'),
    ('Jaaks P et al. Effective drug combinations in breast, colon and pancreatic cancer cells. Nature 2022', '10.1038/s41586-022-04437-2', 'combo screen'),
    ('Marshall EA et al. Small non-coding RNA transcriptome of the NCI-60 cell line panel. Sci Data 2017 (Table 1 doubling times)', '10.1038/sdata.2017.157', 'doubling times'),
    ('Chapman PB et al. NEJM 2011 (BRIM-3)', '10.1056/NEJMoa1103782', 'trials'),
    ('Hauschild A et al. Lancet 2012 (BREAK-3)', '10.1016/S0140-6736(12)60868-X', 'trials'),
    ('Long GV et al. NEJM 2014 (COMBI-d)', '10.1056/NEJMoa1406037', 'trials'),
    ('Long GV et al. Lancet 2015 (COMBI-d OS)', '10.1016/S0140-6736(15)60898-4', 'trials'),
    ('Robert C et al. NEJM 2015 (COMBI-v)', '10.1056/NEJMoa1412690', 'trials'),
    ('Kopetz S et al. J Clin Oncol 2015 (vemurafenib CRC)', '10.1200/JCO.2015.63.2497', 'trials'),
    ('Kopetz S et al. NEJM 2019 (BEACON CRC)', '10.1056/NEJMoa1908075', 'trials'),
    ('Tabernero J et al. J Clin Oncol 2021 (BEACON update)', '10.1200/JCO.20.02088', 'trials'),
    ('Corcoran RB et al. J Clin Oncol 2015 (D+T CRC)', '10.1200/JCO.2015.63.2471', 'trials'),
    ('Corcoran RB et al. Cancer Discov 2018 (D+T+P CRC)', '10.1158/2159-8290.CD-17-1226', 'trials'),
    ('Karapetis CS et al. NEJM 2008 (CO.17 KRAS)', '10.1056/NEJMoa0804385', 'trials'),
    ('Amado RG et al. J Clin Oncol 2008 (panitumumab KRAS)', '10.1200/JCO.2007.14.7116', 'trials'),
    ('Blumenschein GR et al. Ann Oncol 2015 (trametinib KRAS NSCLC)', '10.1093/annonc/mdv072', 'trials'),
    ('Bennouna J et al. Invest New Drugs 2011 (selumetinib CRC), PMID 20127139', '10.1007/s10637-010-9392-8', 'trials'),
    ('Bodoky G et al. Invest New Drugs 2012 (selumetinib pancreas), PMID 21594619', '10.1007/s10637-011-9687-4', 'trials'),
    ('Hainsworth JD et al. J Thorac Oncol 2010 (selumetinib NSCLC), PMID 20802351', '10.1097/JTO.0b013e3181e8b3a3', 'trials'),
    ('Infante JR et al. Lancet Oncol 2012 (trametinib phase I)', '10.1016/S1470-2045(12)70269-3', 'trials, PK'),
    ('Jänne PA et al. Lancet Oncol 2013 (selumetinib + docetaxel)', '10.1016/S1470-2045(12)70489-8', 'trials (note)'),
    ('Shimizu T et al. Clin Cancer Res 2012 (dual-pathway phase I)', '10.1158/1078-0432.CCR-11-2381', 'trials'),
    ('Bedard PL et al. Clin Cancer Res 2015 (buparlisib + trametinib)', '10.1158/1078-0432.CCR-14-1814', 'trials'),
    ('Corcoran RB et al. Clin Cancer Res 2024 (navitoclax + trametinib)', '10.1158/1078-0432.CCR-23-3135', 'trials'),
    ('Hidalgo M et al. J Clin Oncol 2001 (erlotinib phase I), PMID 11432895', '10.1200/JCO.2001.19.13.3267', 'PK'),
    ('Lu JF et al. Clin Pharmacol Ther 2006 (erlotinib pop-PK), PMID 16890575', '10.1016/j.clpt.2006.04.007', 'PK'),
    ('Wilson WH et al. Lancet Oncol 2010 (navitoclax phase I), PMID 21094089, PMC3025495', '10.1016/S1470-2045(10)70261-8', 'PK'),
    ('Gandhi L et al. J Clin Oncol 2011 (navitoclax solid tumours), PMID 21282543', '10.1200/JCO.2010.31.6208', 'PK'),
    ('Lito P et al. Cancer Cell 2012', '10.1016/j.ccr.2012.10.009', 'anchors'),
    ('Corcoran RB et al. Cancer Discov 2012', '10.1158/2159-8290.CD-11-0341', 'anchors'),
    ('Prahallad A et al. Nature 2012', '10.1038/nature10868', 'anchors'),
    ('Poulikakos PI et al. Nature 2010', '10.1038/nature08902', 'anchors'),
    ('Hatzivassiliou G et al. Nature 2010', '10.1038/nature08833', 'anchors'),
    ('Turke AB et al. Cancer Res 2012', '10.1158/0008-5472.CAN-11-3747', 'anchors'),
    ('Misale S et al. Nature 2012', '10.1038/nature11156', 'anchors'),
    ('Diaz LA et al. Nature 2012', '10.1038/nature11219', 'anchors'),
    ('Spencer SL et al. Nature 2009', '10.1038/nature08012', 'anchors'),
    ('Niepel M et al. Cell Syst 2019', '10.1016/j.cels.2019.06.005', 'anchors'),
    ('Akbani R et al. Mol Cell Proteomics 2014 (RPPA workshop report)', '10.1074/mcp.O113.034918', 'anchors'),
    ('Spatial normalization of RPPA data. PLoS One 2014', '10.1371/journal.pone.0097213', 'anchors'),
    ('Hafner M et al. Nat Methods 2016', '10.1038/nmeth.3853', 'anchors'),
    ('Vlachogiannis G et al. Science 2018', '10.1126/science.aao2774', 'anchors'),
    ('Kornauth C et al. Cancer Discov 2022', '10.1158/2159-8290.CD-21-0538', 'anchors'),
    ('Klinger B et al. Mol Syst Biol 2013', '10.1038/msb.2013.29', 'anchors'),
    ('Eduati F et al. Cancer Res 2017', '10.1158/0008-5472.CAN-17-0078', 'anchors'),
    ('Eduati F et al. Mol Syst Biol 2020', '10.15252/msb.20188664', 'anchors'),
    ('Kirouac DC et al. npj Syst Biol Appl 2017', '10.1038/s41540-017-0016-1', 'anchors'),
    ('Fröhlich F et al. Cell Syst 2018', '10.1016/j.cels.2018.10.013', 'anchors'),
    ('Gerosa L et al. Cell Syst 2020', '10.1016/j.cels.2020.10.002', 'anchors'),
]
OTHER = [
    ('FDA labels (openFDA): Zelboraf, Tafinlar, Mekinist, Braftovi, Mektovi, Erbitux, Vectibix, erlotinib (generic), Iressa, Piqray, Venclexta', 'raw/fda_labels'),
    ('EMA SmPCs: Tafinlar, Mekinist, Braftovi, Tarceva, Piqray, Venclyxto, Erbitux', 'raw/ema'),
]
md = ['# Sources', '', f'All downloads on {DATE}. Raw files are untrusted data and were parsed only, never executed. Sizes are in bytes.',
      'Licence column: Sanger (GDSC / Cell Model Passports / GDSC Combinations) licences were not checked on the sites, so they are marked unverified.', '',
      '## Downloaded files', '', '| File | URL | Date | Size (B) | Licence |', '|---|---|---|---|---|'] + rows + [
      '', 'Deleted after download (not used): NCI DTP doubling-time HTML (page redirected), Tarceva openFDA query (404), EMA Iressa PDF (HTTP 429).', '',
      '## Papers', '', '| Paper | DOI | Used for |', '|---|---|---|'] + [f'| {p} | {d} | {u} |' for p, d, u in PAPERS] + [
      '', '## Regulatory documents', ''] + [f'- {a} ({b})' for a, b in OTHER] + ['']
open(os.path.join(S, 'SOURCES.md'), 'w').write('\n'.join(md))
print(len(rows), 'files,', len(PAPERS), 'papers')
