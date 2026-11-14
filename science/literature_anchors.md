# Literature anchors for the Remission simulator

Retrieved 2026-11-14. Quotes and numbers come from the PMC full text saved in `raw/pubmed/` (PMC id given), or from the PubMed abstract when no PMC copy exists. Anything not found in the text is marked **unverified**.

## (i) ERK rebound after RAF/MEK inhibition: time course

- **Lito et al., Cancer Cell 2012** (doi 10.1016/j.ccr.2012.10.009, PMID 23153539, PMC3713778).
  - In BRAF V600E melanoma, vemurafenib inhibits pERK maximally within about 2 h. pERK then rebounds, partly within about 8 h and more by 24–48 h.
  - The rebound was pronounced in 4 of 7 lines and marginal in 2 more. Dabrafenib produced it too.
  - Mechanism: loss of ERK-dependent feedback (SPRY2, DUSPs) raises Ras-GTP. Ras-GTP was not appreciably raised at 12 h but was induced by 24–48 h, and RAF-inhibitor-resistant RAF dimers follow.
  - Once the rebound has happened, re-adding vemurafenib inhibits pMEK much less well, but MEK inhibition still works.
  - In SkMel-28 xenografts, more than 50% of ERK stayed phosphorylated on single-agent RAF or MEK inhibition. RAF + MEK together gave near-complete inhibition.
  - Fold recovery as a number: **unverified** (shown only in Figure 1F blots).
- **Corcoran et al., Cancer Discov 2012** (doi 10.1158/2159-8290.CD-11-0341, PMID 22448344, PMC3308191).
  - In BRAF V600E colorectal (CRC) lines, vemurafenib suppresses pERK only transiently. pERK starts to recover at 3–6 h and reaches about 25–50% of baseline by 24 h.
  - In BRAF V600E melanoma lines, suppression is sustained over 24 h.
  - Reactivation runs through EGFR → RAS → CRAF.
- **Prahallad et al., Nature 2012** (doi 10.1038/nature10868, PMID 22281684; no PMC copy, abstract only).
  - BRAF V600E colon cancer cells show rapid feedback activation of EGFR after BRAF inhibition. Melanoma expresses little EGFR, so it lacks this feedback.
  - Hours and fold change: **unverified** (abstract only).

## (ii) EGFR / receptor signalling after BRAF inhibition: CRC vs melanoma

- **Corcoran 2012** (PMC3308191):
  - BRAF-mutant CRC lines have higher baseline EGFR and phospho-EGFR than BRAF-mutant melanoma lines.
  - In patient tumours, 60% of BRAF-mutant CRCs (n=10) had high P-EGFR (IHC score 2–3), against 18% of melanomas (n=11), p<0.05.
  - Vemurafenib did **not** markedly raise P-EGFR itself. Instead, existing EGFR signalling is de-repressed downstream, raising RAS-GTP and pERK.
- **Prahallad 2012** reports EGFR feedback activation in colon cells. Magnitude: **unverified**.
- **Simulator implication:** in CRC, model a higher baseline RTK input with ERK-dependent negative feedback rather than strong receptor induction. Melanoma gets low EGFR.

## (iii) Paradoxical ERK activation in BRAF wild-type cells

- **Poulikakos et al., Nature 2010** (doi 10.1038/nature08902, PMID 20179705, PMC3178447).
  - ATP-competitive RAF inhibitors (PLX4720 and others) rapidly increase pMEK and pERK in BRAF wild-type cells.
  - This requires RAS activity and drug binding to one protomer of a RAF dimer, which then transactivates the other.
  - Cells co-expressing BRAF V600E and mutant RAS are not inhibited.
- **Hatzivassiliou et al., Nature 2010** (doi 10.1038/nature08833, PMID 20130576).
  - RAF inhibitors prime wild-type RAF (CRAF) to activate MAPK and enhance growth of KRAS-mutant xenografts.
- Fold change: **unverified.** Neither the text nor the abstract states a number; it appears only as Figure 1 blots. If a figure is needed, digitise those blots.

## (iv) AKT activation after MEK inhibition

- **Turke et al., Cancer Res 2012** (doi 10.1158/0008-5472.CAN-11-3747, PMID 22552284, PMC3515079).
  - In EGFR-mutant and HER2-amplified lines, selumetinib (AZD6244), GSK1120212 (trametinib) and PD0325901 increased AKT phosphorylation at both S473 and T308. Phosphorylation of the AKT targets GSK3α/β, ATP citrate lyase and PRAS40 also rose. The effect was seen in an H1975 xenograft in vivo.
  - Mechanism: loss of ERK's inhibitory phosphorylation of the EGFR/HER2 juxtamembrane threonine increases ERBB3 phosphorylation, which raises PI3K activity and PIP3.
  - Fold change: **unverified** (blots only).
- **Klinger et al., Mol Syst Biol 2013** (PMC3964313): in CRC lines, MEK inhibition raises pAKT through EGFR-mediated feedback.

## (v) Pre-existing KRAS-mutant subclones under EGFR blockade

- **Misale et al., Nature 2012** (doi 10.1038/nature11156, PMID 22722830, PMC3927413).
  - In patients relapsing on cetuximab or panitumumab, 6 of 10 (60%) acquired secondary KRAS mutations; 1 more had KRAS amplification.
  - KRAS-mutant allele percentage in resistant tumours ranged from 0.4% to 17%.
  - In the parental Lim1215 line, about 0.2% of cells carried KRAS G13D, and another mutation was present at roughly 1:40,000.
  - KRAS-mutant alleles were detected in plasma up to 10 months before radiographic progression. Clinical progression typically comes at 9–18 months.
- **Diaz et al., Nature 2012** (doi 10.1038/nature11219, PMID 22722843, PMC3436069).
  - On panitumumab, 9 of 24 patients (38%) whose tumours started KRAS wild-type developed detectable KRAS mutations in serum; 3 of them developed several different mutations.
  - Mutations appeared consistently, generally 5–6 months after treatment began.
  - Modelling indicated the mutant subclones existed before therapy.

## (vi) Fractional killing and cell-to-cell variability in apoptosis

- **Spencer et al., Nature 2009** (doi 10.1038/nature08012, PMID 19363473, PMC2858974).
  - With TRAIL alone, a fraction of clonal HeLa cells always survives. With TRAIL plus cycloheximide, all cells die.
  - Time to death varies widely between cells and is set by pre-existing differences in protein levels.
  - Sister cells are correlated in time to death; the correlation decays with a half-life of about 1.5 h with protein synthesis and about 11 h with cycloheximide. Sisters lose memory of each other within about 50 h (2 generations).
  - The surviving fraction at a given TRAIL dose: **unverified** (figures only).

## (vii) Technical variability of the assays

- **Plate viability assays.** Niepel et al., Cell Syst 2019 (doi 10.1016/j.cels.2019.06.005, PMID 31302153, PMC6700527).
  - Technical replicates were "very similar".
  - Biological replicates at one centre (MCF 10A + trametinib) had log10 GR50 and GRmax within 1.4 SD across three scientists.
  - Between centres, initial GR50 measurements differed by up to 200-fold. Most of that came from context-dependent factors.
  - The commonly quoted ~10–20% CV for plate assays: **unverified** (no CV number in the text searched).
- **RPPA.**
  - Akbani-group papers use a slide QC threshold of mean biological-replicate CV ≤ 15% (Spatial normalization of RPPA data, PLoS One 2014, doi 10.1371/journal.pone.0097213, PMC4264691).
  - The RPPA Society workshop report (Akbani et al., Mol Cell Proteomics 2014, doi 10.1074/mcp.O113.034918, PMID 24777629, PMC4083105) says clinical assays generally need CV ≤ 15%.
  - These two figures come from a web-search summary and were not re-read in full text: **unverified**.
- **Suggested simulator noise:** about 10–15% multiplicative noise on viability and on phospho-readouts. This is a modelling choice, not a measured value.

## (viii) GR metrics: targeted drugs are often cytostatic

- **Hafner et al., Nat Methods 2016** (doi 10.1038/nmeth.3853, PMID 27135972, PMC4887336).
  - GR(c) = 2^(k(c)/k(0)) − 1. GR = 0 means complete cytostasis; values between 0 and −1 mean cell death.
  - IC50 and Emax depend strongly on division time and assay length; GR50 and GRmax do not.
  - In their simulations they used division times of 1.8, 2.4 and 3.9 d, the lower quartile, median and upper quartile of their cell-line panel.
  - Typical GRmax ranges by drug class: **unverified.** The text checked gives no class ranges; they are only in the figures and the LINCS GR browser data. To anchor numerically, take them from that data.

## (ix) Functional precision oncology predicting patient response

- **Vlachogiannis et al., Science 2018** (doi 10.1126/science.aao2774, PMID 29472484, PMC6112415).
  - Patient-derived organoids (PDOs) from metastatic GI cancers forecast response to targeted agents or chemotherapy with 100% sensitivity, 93% specificity, 88% PPV and 100% NPV (Fisher p<0.0001, table S7).
- **Kornauth et al., Cancer Discov 2022** (EXALT; doi 10.1158/2159-8290.CD-21-0538, PMID 34635570, PMC9762339).
  - Single-cell functional precision medicine (scFPM) in aggressive haematological cancers: 143 patients, of whom 56 were treated by scFPM.
  - 54% had clinical benefit (PFS ratio > 1.3); 12 of 56 (21%) had exceptional responses (PFS ≥ 3× expected).
  - Turnaround: median 5 days (1–33) from sampling to scFPM report, and 21 days (4–77) from sampling to treatment.
  - Note: haematological, not solid tumours.

## (x) Mechanistic and logic models of this network fitted to perturbation data

One-liners summarise the abstracts (full text read only for Klinger).


- **Klinger et al., Mol Syst Biol 2013** (doi 10.1038/msb.2013.29, PMID 23752269, PMC3964313). Modular response analysis model of EGFR/MAPK/PI3K in CRC lines, fitted to phospho-protein data after ligand and inhibitor perturbations. It predicted EGFR-mediated feedback to AKT under MEK inhibition and synergy of MEK + EGFR inhibition.
- **Eduati et al., Cancer Res 2017** (doi 10.1158/0008-5472.CAN-17-0078, PMID 28381545). CellNOpt logic ODE models fitted per line to 14 phosphoproteins × 43 perturbation conditions in 14 CRC lines. Model parameters predicted drug sensitivity and candidate combinations.
- **Eduati et al., Mol Syst Biol 2020** (doi 10.15252/msb.20188664, PMID 32073727). Microfluidic perturbation of pancreatic cancer biopsy cells, with caspase/apoptosis readouts fitted by patient-specific logic models to suggest personalised combinations.
- **Kirouac et al., npj Syst Biol Appl 2017** (doi 10.1038/s41540-017-0016-1, PMID 28649441). Quantitative systems pharmacology model of EGFR–RAF–MEK–ERK feedback plus tumour growth in BRAF V600E CRC, fitted to cell-line and xenograft data. It predicted an ERK-inhibitor ORR of about 17% (the trial saw 2/18).
- **Fröhlich et al., Cell Syst 2018** (doi 10.1016/j.cels.2018.10.013, PMID 30503647). Large mechanistic pan-cancer signalling model (>1,200 species, >2,600 reactions) parameterised with CCLE expression and fitted to drug-response data.
- **Gerosa et al., Cell Syst 2020** (doi 10.1016/j.cels.2020.10.002, PMID 33113355). Live-cell ERK activity in BRAF V600E melanoma under RAF/MEK inhibitors. Drug-adapted persister cells show sporadic ERK pulses that drive slow proliferation, modelled with a stochastic pulse model.
