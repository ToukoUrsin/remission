"""Validate the Remission cell model against public data. One command:

    science/.venv/bin/python science/validation/validate_model.py            # ~3-6 min on 16 procs

Writes docs/validation.json and docs/img/val_*.png. Reads only public inputs in
science/ and the model in remission/. Never changes model parameters.

Parts (see docs/MODEL.md):
  a  72 h dose-response per virtual tumour vs GDSC2 cell lines (IC50, AUC)
  b  five mechanistic checks against classic papers
  c  label-dose guideline regimens in the unselected virtual population vs trial ORRs
  d  model PK constants vs label PK
  e  benchmark cohort selection (summary of the development scan)
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time
from multiprocessing import Pool

import numpy as np
from scipy.optimize import curve_fit

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, ROOT)
from remission import cohort as C  # noqa: E402
from remission import model as M  # noqa: E402
from remission import patients as P  # noqa: E402

SCI = os.path.join(ROOT, "science")
DOCS = os.path.join(ROOT, "docs")
IMG = os.path.join(DOCS, "img")

# Each model drug class -> the GDSC2 compound it is compared with and that compound's screened range (uM).
GDSC_DRUG = {
    "RAFi": ("Dabrafenib", 0.010005, 10.0),
    "MEKi": ("Trametinib", 0.001001, 1.0),
    "EGFRi": ("Erlotinib", 0.002001, 2.0),
    "METi": ("Crizotinib", 0.001001, 10.0),
    "KRASi": ("KRAS (G12C) Inhibitor-12", 0.020011, 20.0),
    "PI3Ki": ("Alpelisib", 0.005003, 5.0),
    "BCLi": ("Navitoclax", 0.002001, 10.0),
    "MCLi": ("AZD5991", 0.010005, 10.0),
}
NCONC = 9
ASSAY_DAYS = 3.0


# ---------------------------------------------------------------------------
# GDSC-style curve fit: viability = 1 / (1 + exp((ln c - ln IC50) / s)), top 1, bottom 0
# (the two-parameter sigmoid of Vis et al. 2016 that GDSC2 fits; IC50 may be extrapolated).
# ---------------------------------------------------------------------------
def _sig(x, xmid, s):
    return 1.0 / (1.0 + np.exp((x - xmid) / s))


def fit_curve(conc, v):
    x = np.log(conc)
    lo, hi = x[0], x[-1]
    v = np.clip(v, 0.0, 1.2)
    try:
        # guess from the first crossing of 0.5
        below = np.where(v < 0.5)[0]
        x0 = x[below[0]] if len(below) else hi + 1.0
        popt, _ = curve_fit(_sig, x, v, p0=[x0, 0.6], bounds=([lo - 3.0, 0.05], [hi + 3.0, 6.0]), maxfev=4000)
        xmid, s = popt
    except Exception:
        xmid, s = hi + 3.0, 1.0
    xs = np.linspace(lo, hi, 200)
    auc = float(np.mean(_sig(xs, xmid, s)))
    return float(xmid), auc


# ---------------------------------------------------------------------------
# Per-patient work (runs in worker processes)
# ---------------------------------------------------------------------------
def dose_response(p):
    """72 h viability relative to untreated for the tumour clone mix, each single agent, 9 concentrations."""
    rows = [np.zeros(M.ND)]
    concs = {}
    for d in M.DRUGS:
        _, lo, hi = GDSC_DRUG[d]
        c = np.geomspace(lo, hi, NCONC)
        concs[d] = c
        for x in c:
            r = np.zeros(M.ND)
            r[M.IDX[d]] = x
            rows.append(r)
    out = M.simulate_wells(p["K"], p["tumour"], p["frac"], np.array(rows), ASSAY_DAYS, rng=None, read_days=[ASSAY_DAYS])
    N = out["N"][-1].sum(axis=1)
    v = N[1:] / N[0]
    res = {}
    for i, d in enumerate(M.DRUGS):
        vv = v[i * NCONC:(i + 1) * NCONC]
        ln_ic50, auc = fit_curve(concs[d], vv)
        res[d] = {"ln_ic50": ln_ic50, "auc": auc, "v": vv.round(4).tolist()}
    res["growth72"] = float(N[0] / M.SEED_CELLS)
    return res


def progression_day(T, T0):
    """First day the burden is >= 1.2x its running minimum and at least 2 % of baseline above it (RECIST-like PD)."""
    runmin = np.minimum.accumulate(T)
    hit = np.where((T >= 1.2 * runmin) & (T - runmin >= 0.02 * T0) & (np.arange(len(T)) > 0))[0]
    return int(hit[0]) if len(hit) else None


def clinical(p, lattice_reps):
    reg = C.constant(C.GUIDELINES[p["archetype"]][0])
    plan = M.expand_regimen(reg)
    mean = M.simulate_tissue_mean(p["K"], p["healthy"], p["tumour"], p["frac"], plan[None])
    Tm = mean["tumour"][0].sum(axis=1)
    out = {"mean_reward": float(mean["reward"][0]), "mean_cure": float(mean["cure"][0]), "harm": float(mean["harm"][0]),
           "mean_best": float(Tm[1:M.TREAT_DAYS + 1].min() / Tm[0] - 1.0), "mean_pd": progression_day(Tm, Tm[0])}
    if lattice_reps:
        lat = M.simulate_tissue(p["K"], p["healthy"], p["tumour"], p["frac"], reg,
                                seed=P.seed_from("validation-lattice", p["seed"]), reps=lattice_reps)
        T = lat["tumour"].sum(axis=2)          # [reps, days+1]
        best = (T[:, 1:M.TREAT_DAYS + 1].min(axis=1) / T[:, 0] - 1.0)
        out.update(lat_best=best.tolist(), lat_reward=lat["reward"].tolist(),
                   lat_elim=[bool(T[r, -1] == 0) for r in range(lattice_reps)],
                   lat_pd=[progression_day(T[r], T[r, 0]) for r in range(lattice_reps)],
                   lat_curve=(T.mean(axis=0)[::7] / T[:, 0].mean()).round(4).tolist())
    return out


def work(args):
    i, lattice_reps = args
    seed = P.seed_from("validation", i)
    p = P.make_patient(seed)
    major = p["tumour"][0]
    return {"i": i, "seed": seed, "arch": p["archetype"], "alts": list(major["alterations"]),
            "minors": [c["name"][6:] for c in p["tumour"][1:]], "frac": p["frac"].tolist(),
            "dr": dose_response(p), "clin": clinical(p, lattice_reps)}


# ---------------------------------------------------------------------------
# GDSC2 reference groups
# ---------------------------------------------------------------------------
def gdsc_groups():
    import pandas as pd
    g = pd.read_csv(os.path.join(SCI, "raw/gdsc/GDSC2_fitted_dose_response_27Oct23.csv"),
                    usecols=["SANGER_MODEL_ID", "DRUG_NAME", "LN_IC50", "AUC", "MAX_CONC"])
    a = pd.read_csv(os.path.join(SCI, "cell_line_annotations.csv"))
    a = a[a["mutation_data_available"] == True]  # noqa: E712
    a = a[a["tissue"] != "Haematopoietic and Lymphoid"]
    ct = a["cancer_type"]
    nsclc = ct.isin(["Non-Small Cell Lung Carcinoma", "Squamous Cell Lung Carcinoma"])
    t = lambda s: s.astype(str) == "True"  # noqa: E731
    groups = {
        "BRAF V600E melanoma": a[(ct == "Melanoma") & t(a["BRAF_V600E"])],
        "BRAF V600E colorectal": a[(ct == "Colorectal Carcinoma") & t(a["BRAF_V600E"])],
        "KRAS G12C (any solid)": a[a["KRAS_driver_variants"].astype(str).str.contains("G12C")],
        "KRAS-mutant colorectal": a[(ct == "Colorectal Carcinoma") & t(a["KRAS_mut"])],
        "EGFR-mutant NSCLC": a[nsclc & t(a["EGFR_driver_mut"])],
        "PIK3CA-mut or PTEN-null (any solid)": a[t(a["PIK3CA_mut"]) | t(a["PTEN_null"])],
        "melanoma + colorectal + NSCLC (all)": a[ct.isin(["Melanoma", "Colorectal Carcinoma"]) | nsclc],
    }
    ids = {k: set(v["model_id"]) for k, v in groups.items()}
    return g, ids


# Which simulated subgroup is compared with which GDSC group, per drug class (figure panels).
COMPARE = {
    "RAFi": [("BRAF_MEL", "BRAF V600E melanoma"), ("BRAF_CRC", "BRAF V600E colorectal")],
    "MEKi": [("BRAF_MEL", "BRAF V600E melanoma"), ("BRAF_CRC", "BRAF V600E colorectal"), ("KRAS", "KRAS G12C (any solid)"),
             ("KRAS", "KRAS-mutant colorectal")],
    "EGFRi": [("EGFR", "EGFR-mutant NSCLC"), ("EGFR", "EGFR-mutant NSCLC", "Osimertinib"), ("BRAF_CRC", "BRAF V600E colorectal")],
    "METi": [("all", "melanoma + colorectal + NSCLC (all)")],
    "KRASi": [("KRAS", "KRAS G12C (any solid)")],
    "PI3Ki": [("pi3k_alt", "PIK3CA-mut or PTEN-null (any solid)")],
    "BCLi": [("all", "melanoma + colorectal + NSCLC (all)")],
    "MCLi": [("all", "melanoma + colorectal + NSCLC (all)")],
}
SIM_LABEL = {"BRAF_MEL": "BRAF mel", "BRAF_CRC": "BRAF CRC", "KRAS": "KRAS G12C CRC", "EGFR": "EGFR lung",
             "all": "all", "pi3k_alt": "PIK3CA/PTEN"}


def sim_subset(res, key):
    if key == "all":
        return res
    if key == "pi3k_alt":
        return [r for r in res if ("pik3ca" in r["alts"] or "pten" in r["alts"])]
    return [r for r in res if r["arch"] == key]


def q(x):
    x = np.asarray(x, dtype=float)
    if len(x) == 0:
        return None
    return {"n": int(len(x)), "median": round(float(np.median(x)), 3),
            "q25": round(float(np.percentile(x, 25)), 3), "q75": round(float(np.percentile(x, 75)), 3)}


def part_a(res):
    g, ids = gdsc_groups()
    table = {}
    for d, pairs in COMPARE.items():
        rows = []
        for pair in pairs:
            skey, gkey = pair[:2]
            name, lo, hi = GDSC_DRUG[d]
            if len(pair) > 2:
                name = pair[2]
                lo, hi = {"Osimertinib": (0.001001, 1.0)}[name]
            gd = g[g["DRUG_NAME"] == name]
            sub = sim_subset(res, skey)
            s_ic = [r["dr"][d]["ln_ic50"] for r in sub]
            s_auc = [r["dr"][d]["auc"] for r in sub]
            gg = gd[gd["SANGER_MODEL_ID"].isin(ids[gkey])]
            d_ic, d_auc = gg["LN_IC50"].values, gg["AUC"].values
            row = {"sim_group": SIM_LABEL[skey], "gdsc_group": gkey, "gdsc_drug": name, "range_uM": [lo, hi],
                   "sim_ln_ic50": q(s_ic), "gdsc_ln_ic50": q(d_ic), "sim_auc": q(s_auc), "gdsc_auc": q(d_auc),
                   "sim_share_ic50_above_max": round(float(np.mean(np.array(s_ic) > np.log(hi))), 3) if len(s_ic) else None,
                   "gdsc_share_ic50_above_max": round(float(np.mean(d_ic > np.log(hi))), 3) if len(d_ic) else None}
            if len(s_ic) and len(d_ic):
                row["delta_log10_ic50_median"] = round((np.median(s_ic) - np.median(d_ic)) / np.log(10), 2)
                row["delta_auc_median"] = round(float(np.median(s_auc) - np.median(d_auc)), 3)
            row["_sim_ic"], row["_gdsc_ic"] = list(map(float, s_ic)), list(map(float, d_ic))
            rows.append(row)
        table[d] = rows
    return table


# ---------------------------------------------------------------------------
# Part b: mechanistic checks. Signalling time courses use the model's own functions.
# ---------------------------------------------------------------------------
def time_course(K, clone, conc, hours=48, dt_h=0.5):
    """E/E0, A/A0 and RAS-GTP/RAS-GTP0 of one clone at constant drug, from M.signalling + M.step_feedback
    (identical to the signalling inside M.simulate_wells; only cell numbers are left out)."""
    p = M.pack([clone])
    occ = M.occupancy(np.asarray(conc, float)[None, :], np.asarray(K, float), p["kshift"])
    E0, A0 = M.baseline(p)
    occ0 = np.zeros_like(occ)

    def ras(occ_, F_):
        Re = p["egfr"] * F_ * (1 - occ_[..., M.iEGFR])
        Ra = p["alt"] * F_ * (1 - occ_[..., M.iMET])
        return M.ARAS * (Re + Ra) * (1 + 0 * p["ras"])     # receptor-driven RAS-GTP (wild-type RAS pool)
    R0 = ras(occ0, np.ones(1))
    F = np.ones(1)
    dt = dt_h / 24.0
    n = int(round(hours / dt_h))
    t, Es, As, Rs = [], [], [], []
    for s in range(n + 1):
        E, A = M.signalling(p, occ, F)
        t.append(s * dt_h); Es.append(float(E[0] / E0[0])); As.append(float(A[0] / A0[0])); Rs.append(float(ras(occ, F)[0] / R0[0]))
        F = M.step_feedback(p, F, E, E0, dt)
    return np.array(t), np.array(Es), np.array(As), np.array(Rs)


def conc_vec(**kw):
    c = np.zeros(M.ND)
    for k, v in kw.items():
        c[M.IDX[k]] = v
    return c


CSS = {d: M.PK[d]["css"] for d in M.DRUGS}


def part_b(patients):
    out, curves = {}, {}
    by = lambda a: [p for p in patients if p["archetype"] == a]  # noqa: E731

    # (i) ERK rebound after RAFi in BRAF CRC vs melanoma; suppressed by adding EGFRi
    def band(pl, conc):
        X = np.array([time_course(p["K"], p["tumour"][0], conc)[1] for p in pl])
        return X
    t = time_course(patients[0]["K"], patients[0]["tumour"][0], conc_vec())[0]
    crc, mel = by("BRAF_CRC"), by("BRAF_MEL")
    E_crc = band(crc, conc_vec(RAFi=CSS["RAFi"]))
    E_mel = band(mel, conc_vec(RAFi=CSS["RAFi"]))
    E_crc_e = band(crc, conc_vec(RAFi=CSS["RAFi"], EGFRi=CSS["EGFRi"]))
    i24 = int(np.where(t == 24)[0][0])
    m = lambda X, j: float(np.median(X[:, j]))  # noqa: E731
    r = {"crc_nadir": m(E_crc, 0), "crc_24h": m(E_crc, i24), "mel_nadir": m(E_mel, 0), "mel_24h": m(E_mel, i24),
         "crc_rafi_egfri_24h": m(E_crc_e, i24), "literature": "Corcoran 2012: CRC pERK recovers to ~25-50 % of baseline by 24 h; melanoma suppression sustained over 24 h"}
    r["pass"] = bool(r["crc_24h"] >= 2 * r["crc_nadir"] and r["mel_24h"] <= 0.5 * r["crc_24h"] and r["crc_rafi_egfri_24h"] <= 0.5 * r["crc_24h"])
    r["crc_24h_in_25_50pct"] = bool(0.25 <= r["crc_24h"] <= 0.50)
    E1c, E1m = band(crc, conc_vec(RAFi=1.0)), band(mel, conc_vec(RAFi=1.0))
    E_crc_rm = band(crc, conc_vec(RAFi=CSS["RAFi"], MEKi=CSS["MEKi"]))
    r.update(crc_nadir_1uM=m(E1c, 0), crc_24h_1uM=m(E1c, i24), mel_nadir_1uM=m(E1m, 0), mel_24h_1uM=m(E1m, i24),
             crc_rafi_meki_24h=m(E_crc_rm, i24),
             crc_rafi_meki_biopsy_note="Corcoran 2015 (clinical_trials.json): on-treatment pERK reduction 47 +/- 24 % in CRC biopsies on dabrafenib+trametinib")
    out["i_erk_rebound"] = r
    curves["i"] = (t, E_crc, E_mel, E_crc_e)

    # (ii) paradoxical ERK activation in RAS-mutant / BRAF-wt cells
    cs = np.geomspace(0.001, 10, 13)
    def para(pl, clone_of):
        return np.array([[time_course(p["K"], clone_of(p), conc_vec(RAFi=c), hours=0)[1][0] for c in cs] for p in pl])
    P_kras = para(by("KRAS"), lambda p: p["tumour"][0])
    P_mel = para(mel, lambda p: p["tumour"][0])
    P_hea = para(patients[:200], lambda p: p["healthy"])
    # 72 h growth of KRAS-mutant tumour cells under RAFi at css relative to untreated
    g72 = []
    for p in by("KRAS"):
        o = M.simulate_wells(p["K"], p["tumour"][:1], [1.0], np.stack([conc_vec(), conc_vec(RAFi=CSS["RAFi"])]), 3.0, read_days=[3.0])
        n = o["N"][-1].sum(axis=1)
        g72.append(n[1] / n[0])
    r = {"kras_max_erk_fold": float(np.median(P_kras.max(axis=1))), "healthy_max_erk_fold": float(np.median(P_hea.max(axis=1))),
         "braf_mel_erk_at_css": float(np.median(P_mel[:, np.argmin(abs(cs - CSS["RAFi"]))])),
         "kras_viability72_rafi_css": float(np.median(g72)),
         "literature": "Poulikakos 2010, Hatzivassiliou 2010: RAF inhibitors raise pMEK/pERK in BRAF-wt cells with active RAS; fold not stated in text (unverified)"}
    r["pass"] = bool(r["kras_max_erk_fold"] > 1.2 and r["braf_mel_erk_at_css"] < 0.3)
    out["ii_paradox"] = r
    curves["ii"] = (cs, P_kras, P_hea, P_mel)

    # (iii) relief of feedback after MEK inhibition (Lito 2012); AKT up (Turke 2012); RAF+MEK deeper
    def tc(pl, conc):
        X = [time_course(p["K"], p["tumour"][0], conc) for p in pl]
        return np.array([x[1] for x in X]), np.array([x[2] for x in X]), np.array([x[3] for x in X])
    Em, Am, Rm = tc(mel, conc_vec(MEKi=CSS["MEKi"]))
    Er, _, _ = tc(mel, conc_vec(RAFi=CSS["RAFi"]))
    Erm, _, _ = tc(mel, conc_vec(RAFi=CSS["RAFi"], MEKi=CSS["MEKi"]))
    Ec, Ac, Rc = tc(crc, conc_vec(MEKi=CSS["MEKi"]))
    i12, i48 = int(np.where(t == 12)[0][0]), int(np.where(t == 48)[0][0])
    r = {"mel_meki_erk_nadir": m(Em, 0), "mel_meki_erk_24h": m(Em, i24), "mel_meki_erk_48h": m(Em, i48),
         "mel_ras_gtp_fold_12h": m(Rm, i12), "mel_ras_gtp_fold_24h": m(Rm, i24), "mel_ras_gtp_fold_48h": m(Rm, i48),
         "mel_akt_fold_24h": m(Am, i24), "crc_meki_erk_24h": m(Ec, i24), "crc_akt_fold_24h": m(Ac, i24),
         "mel_rafi_erk_24h": m(Er, i24), "mel_rafi_meki_erk_24h": m(Erm, i24),
         "literature": "Lito 2012: pERK partly rebounds by 24-48 h; Ras-GTP not raised at 12 h, induced by 24-48 h; RAF+MEK near-complete. Turke 2012: MEK inhibition raises pAKT."}
    r["pass"] = bool(r["mel_meki_erk_48h"] > 1.5 * r["mel_meki_erk_nadir"] and r["mel_ras_gtp_fold_24h"] > 1.5
                     and r["mel_akt_fold_24h"] > 1.0 and r["mel_rafi_meki_erk_24h"] < min(r["mel_rafi_erk_24h"], r["mel_meki_erk_24h"]))
    r["ras_gtp_already_up_at_12h"] = bool(r["mel_ras_gtp_fold_12h"] > 1.5)
    out["iii_feedback_relief"] = r
    # Lito 2012's own experiment: RAF inhibitor in BRAF V600E melanoma, pERK nadir then partial rebound by 24-48 h
    Er2, _, Rr2 = tc(mel, conc_vec(RAFi=CSS["RAFi"]))
    r2 = {"mel_rafi_erk_nadir": m(Er2, 0), "mel_rafi_erk_8h": m(Er2, int(np.where(t == 8)[0][0])), "mel_rafi_erk_24h": m(Er2, i24),
          "mel_rafi_erk_48h": m(Er2, i48), "mel_rafi_ras_gtp_fold_12h": m(Rr2, i12), "mel_rafi_ras_gtp_fold_24h": m(Rr2, i24),
          "literature": "Lito 2012: vemurafenib in BRAF V600E melanoma, pERK minimal at ~2 h, partial rebound by 8 h, more by 24-48 h; Ras-GTP induced by 24-48 h, not at 12 h"}
    r2["pass"] = bool(r2["mel_rafi_erk_48h"] > 1.3 * r2["mel_rafi_erk_nadir"] and r2["mel_rafi_ras_gtp_fold_24h"] > 1.5)
    r2["ras_gtp_already_up_at_12h"] = bool(r2["mel_rafi_ras_gtp_fold_12h"] > 1.5)
    out["iii_b_lito_rafi_melanoma"] = r2
    curves["iii"] = (t, Em, Rm, Am, Erm)

    # (iv) MET amplification rescues EGFR-mutant cells from EGFRi; METi restores sensitivity (Engelman 2007, Turke 2010)
    egfr = by("EGFR")
    _, lo, hi = GDSC_DRUG["EGFRi"]
    ce = np.geomspace(lo, hi, NCONC)
    ic = {"parental": [], "MET-amp": [], "MET-amp + METi": []}
    V = {k: [] for k in ic}
    for p in egfr:
        par = p["tumour"][0]
        met = dict(par, alt=1.0, name="met-amp")
        for key, clone, extra in (("parental", par, 0.0), ("MET-amp", met, 0.0), ("MET-amp + METi", met, CSS["METi"])):
            rows = [conc_vec(METi=extra)] + [conc_vec(EGFRi=c, METi=extra) for c in ce]
            o = M.simulate_wells(p["K"], [clone], [1.0], np.stack(rows), 3.0, read_days=[3.0])
            n = o["N"][-1][:, 0]
            # viability relative to untreated (no drug at all) as GDSC does
            o0 = M.simulate_wells(p["K"], [clone], [1.0], conc_vec()[None], 3.0, read_days=[3.0])["N"][-1][0, 0]
            v = n[1:] / o0
            V[key].append(v)
            ic[key].append(fit_curve(ce, v)[0])
    r = {k: round(float(np.median(v)) / np.log(10), 2) for k, v in ic.items()}
    r = {"median_log10_ic50_uM": r, "fold_shift_met_amp": round(10 ** (r["MET-amp"] - r["parental"]), 1),
         "fold_shift_with_meti": round(10 ** (r["MET-amp + METi"] - r["parental"]), 1),
         "literature": "Engelman 2007: MET amplification causes gefitinib resistance, MET inhibition restores sensitivity; Turke 2010: pre-existing MET-amplified subpopulations"}
    r["pass"] = bool(r["fold_shift_met_amp"] >= 10 and r["fold_shift_with_meti"] <= 3)
    out["iv_met_bypass"] = r
    curves["iv"] = (ce, {k: np.array(v) for k, v in V.items()})

    # (v) pre-existing KRAS-mutant minor clone under EGFR blockade: time to regrowth vs initial share
    shares = np.array([1e-4, 3e-4, 1e-3, 3e-3, 1e-2, 3e-2, 1e-1])
    plan = M.expand_regimen(C.constant(C.GUIDELINES["EGFR"][0]))
    days_pd, kras_at_pd, never = [], [], []
    for p in egfr[:80]:
        par = p["tumour"][0]
        kr = dict(par, ras=0.85, name="minor-ras")
        row, rowk = [], []
        for f in shares:
            o = M.simulate_tissue_mean(p["K"], p["healthy"], [par, kr], np.array([1 - f, f]), plan[None])
            T = o["tumour"][0]
            tot = T.sum(axis=1)
            d = progression_day(tot, tot[0])
            row.append(np.nan if d is None else d)
            rowk.append(np.nan if d is None else T[d, 1] / max(tot[d], 1e-9))
        days_pd.append(row); kras_at_pd.append(rowk)
    D = np.array(days_pd, float)
    Kf = np.array(kras_at_pd, float)
    med = [float(np.nanmedian(D[:, j])) if np.isfinite(D[:, j]).any() else None for j in range(len(shares))]
    on_tx = [float(np.mean(D[:, j] <= M.TREAT_DAYS)) for j in range(len(shares))]
    ok = [x for x in med if x is not None]
    slope = float(np.polyfit(np.log10(shares[[i for i, x in enumerate(med) if x is not None]]), ok, 1)[0]) if len(ok) > 2 else None
    r = {"shares": shares.tolist(), "median_days_to_regrowth": med, "share_regrowing_on_treatment": on_tx,
         "median_kras_share_at_regrowth": [float(np.nanmedian(Kf[:, j])) if np.isfinite(Kf[:, j]).any() else None for j in range(len(shares))],
         "days_per_10x_share": slope, "n_patients": len(D),
         "literature": "Misale 2012: KRAS-mutant cells pre-exist (~0.2 % in Lim1215), 60 % of relapses gain KRAS mutations; Diaz 2012: KRAS mutants in serum 5-6 months into panitumumab"}
    mono = all(np.diff([x for x in med if x is not None]) <= 0)
    r["pass"] = bool(mono and slope is not None and slope < 0 and np.nanmedian(Kf[:, 3]) > 0.5)
    out["v_kras_outgrowth"] = r
    curves["v"] = (shares, D)
    return out, curves


# ---------------------------------------------------------------------------
# Part c: guideline regimen in the unselected population vs trial ORRs
# ---------------------------------------------------------------------------
TRIALS = {   # pivotal trial ORR (%) for the first guideline regimen of each type; see science/validation_trials.json
    "BRAF_MEL": ("dabrafenib + trametinib, COMBI-d / COMBI-v", 67, [64, 67]),
    "BRAF_CRC": ("encorafenib + binimetinib + cetuximab, BEACON (triplet); D+T+panitumumab", 26, [21, 26.8]),
    "KRAS": ("adagrasib + cetuximab (KRYSTAL-1) / sotorasib + panitumumab (CodeBreaK 300)", 36, [26.4, 46]),
    "EGFR": ("osimertinib (FLAURA); first-generation EGFR-TKI arm", 80, [76, 80]),
}


def part_c(res, lattice_reps):
    out = {}
    for a in P.ARCHETYPES:
        sub = [r["clin"] for r in res if r["arch"] == a]
        o = {"n": len(sub), "regimen": C.GUIDELINES[a][0], "trial": TRIALS[a][0], "trial_orr_pct": TRIALS[a][1],
             "trial_orr_range_pct": TRIALS[a][2]}
        mb = np.array([s["mean_best"] for s in sub])
        o["mean_model"] = {"orr_like_30pct_burden": round(float(np.mean(mb <= -0.30)), 3),
                           "orr_like_51pct_area": round(float(np.mean(mb <= -0.51)), 3),
                           "median_best_change": round(float(np.median(mb)), 3),
                           "share_reward_ge_0.4": round(float(np.mean([s["mean_reward"] >= 0.4 for s in sub])), 3),
                           "mean_cure_probability": round(float(np.mean([s["mean_cure"] for s in sub])), 3)}
        pd_ = [s["mean_pd"] for s in sub]
        o["mean_model"]["share_progressing_on_treatment"] = round(float(np.mean([d is not None and d <= M.TREAT_DAYS for d in pd_])), 3)
        o["mean_model"]["median_days_to_regrowth"] = _median_censored(pd_)
        if lattice_reps:
            lb = np.array([np.mean(s["lat_best"]) for s in sub])
            lr = np.array([np.mean(s["lat_reward"]) for s in sub])
            o["lattice"] = {"reps_per_patient": lattice_reps,
                            "orr_like_30pct_burden": round(float(np.mean(lb <= -0.30)), 3),
                            "orr_like_51pct_area": round(float(np.mean(lb <= -0.51)), 3),
                            "median_best_change": round(float(np.median(lb)), 3),
                            "share_reward_ge_0.4": round(float(np.mean(lr >= 0.4)), 3),
                            "share_eliminated_day364": round(float(np.mean([np.mean(s["lat_elim"]) for s in sub])), 3)}
            pdl = [x for s in sub for x in s["lat_pd"]]
            o["lattice"]["share_progressing_on_treatment"] = round(float(np.mean([d is not None and d <= M.TREAT_DAYS for d in pdl])), 3)
            o["lattice"]["median_days_to_regrowth"] = _median_censored(pdl)
            o["_curves"] = np.array([s["lat_curve"] for s in sub])
            o["_best"] = lb
        out[a] = o
    return out


def _median_censored(days, horizon=M.FOLLOW_DAYS):
    x = np.array([d if d is not None else np.inf for d in days], float)
    if np.mean(np.isinf(x)) >= 0.5:
        return f">{horizon} (no regrowth in {np.mean(np.isinf(x)):.0%})"
    return float(np.median(x))


# ---------------------------------------------------------------------------
# Part d: PK
# ---------------------------------------------------------------------------
REP = {"RAFi": "dabrafenib", "MEKi": "trametinib", "EGFRi": "erlotinib", "METi": None, "KRASi": None,
       "PI3Ki": "alpelisib", "BCLi": "navitoclax (ABT-263; not approved)", "MCLi": None}
NOT_IN_TABLE = {"METi": "capmatinib (label t1/2 6.5 h: unverified, not in science/pk_table.json)",
                "KRASi": "adagrasib (t1/2 ~23-24 h: unverified, not in science/pk_table.json)",
                "MCLi": "no approved MCL-1 inhibitor; 12 h is a modelling choice"}


def part_d():
    pk = json.load(open(os.path.join(SCI, "pk_table.json")))["drugs"]
    by = {d["drug"]: d for d in pk}
    rows = []
    for d in M.DRUGS:
        row = {"class": d, "model_t_half_h": M.PK[d]["t_half_h"], "model_css_uM": M.PK[d]["css"]}
        rep = REP[d]
        if rep and rep in by:
            e = by[rep]
            th = e["terminal_half_life_h"]["value"]
            row["drug"] = rep
            row["label_t_half_h"] = th
            cav = e.get("Cavg_ss", {}).get("uM") or e.get("Cavg_day1_315mg", {}).get("uM")
            fu = e.get("plasma_protein_binding_pct", {}).get("fraction_unbound")
            row["label_cavg_total_uM"] = cav
            row["fraction_unbound"] = fu
            if cav and isinstance(fu, (int, float)):
                row["label_cavg_free_uM"] = round(cav * fu, 4)
                row["model_css_over_free_cavg"] = round(M.PK[d]["css"] / (cav * fu), 2)
            if cav:
                row["model_css_over_total_cavg"] = round(M.PK[d]["css"] / cav, 3)
        else:
            row["drug"] = NOT_IN_TABLE[d]
        rows.append(row)
    return rows


# ---------------------------------------------------------------------------
# Part e: cohort
# ---------------------------------------------------------------------------
def part_e():
    path = os.path.join(SCI, "cohort_scan_summary.json")
    if os.path.exists(path):
        return json.load(open(path))
    return {"note": "science/cohort_scan_summary.json missing; run science/validation/cohort_summary.py"}


# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------
BLUE, ORANGE, GREY, AQUA, VIOLET = "#2a78d6", "#eb6834", "#8a8984", "#1baf7a", "#6250d6"


def _style(plt):
    plt.rcParams.update({"font.size": 13, "axes.titlesize": 14, "axes.labelsize": 13, "legend.fontsize": 11,
                         "axes.spines.top": False, "axes.spines.right": False, "axes.grid": True,
                         "grid.color": "#e6e5e1", "grid.linewidth": 0.8, "figure.dpi": 130, "savefig.dpi": 130})


def fig_a(table, plt):
    fig, axes = plt.subplots(2, 4, figsize=(17, 8.5))
    for ax, d in zip(axes.flat, M.DRUGS):
        rows = table[d]
        pos, labels = [], []
        for k, r in enumerate(rows):
            x = 3 * k
            for dx, key, col in ((0, "_sim_ic", BLUE), (1, "_gdsc_ic", ORANGE)):
                vals = np.array(r[key]) / np.log(10)
                if len(vals):
                    bp = ax.boxplot(vals, positions=[x + dx], widths=0.7, patch_artist=True, showfliers=False,
                                    medianprops=dict(color="black", linewidth=2))
                    for b in bp["boxes"]:
                        b.set(facecolor=col, alpha=0.75, edgecolor=col)
            pos.append(x + 0.5)
            labels.append(r["sim_group"] + ("" if r["gdsc_drug"] == GDSC_DRUG[d][0] else f"\nvs {r['gdsc_drug']}") + f"\n(n={r['sim_ln_ic50']['n'] if r['sim_ln_ic50'] else 0}/"
                          f"{r['gdsc_ln_ic50']['n'] if r['gdsc_ln_ic50'] else 0})")
        lo, hi = GDSC_DRUG[d][1], GDSC_DRUG[d][2]
        ax.axhspan(np.log10(lo), np.log10(hi), color="#f3f2ee", zorder=0)
        ax.axhline(np.log10(M.PK[d]["css"]), color=GREY, ls="--", lw=1.5)
        ax.set_xticks(pos)
        ax.set_xticklabels(labels, fontsize=10)
        ax.set_title(f"{d} vs {GDSC_DRUG[d][0]}", fontsize=13)
        ax.set_ylabel("log10 IC50 (µM)")
    from matplotlib.patches import Patch
    from matplotlib.lines import Line2D
    fig.legend(handles=[Patch(color=BLUE, alpha=0.75, label="simulated virtual tumours"),
                        Patch(color=ORANGE, alpha=0.75, label="GDSC2 cell lines"),
                        Patch(color="#f3f2ee", label="GDSC screened range"),
                        Line2D([], [], color=GREY, ls="--", label="model css at label dose")],
               loc="lower center", ncol=4, frameon=False)
    fig.suptitle("72 h single-agent IC50: simulated tumours vs GDSC2 lines (box = IQR, n = sim/GDSC)", fontsize=15)
    fig.tight_layout(rect=[0, 0.05, 1, 0.96])
    fig.savefig(os.path.join(IMG, "val_gdsc_ic50.png"))
    plt.close(fig)


def fig_b(curves, b, plt):
    fig, axes = plt.subplots(2, 3, figsize=(17, 9.5))
    def band(ax, t, X, col, lab):
        ax.plot(t, np.median(X, axis=0), color=col, lw=2.5, label=lab)
        ax.fill_between(t, np.percentile(X, 25, axis=0), np.percentile(X, 75, axis=0), color=col, alpha=0.18, lw=0)
    t, E_crc, E_mel, E_crc_e = curves["i"]
    ax = axes[0, 0]
    band(ax, t, E_crc, ORANGE, "BRAF CRC, RAFi"); band(ax, t, E_mel, BLUE, "BRAF mel, RAFi"); band(ax, t, E_crc_e, AQUA, "BRAF CRC, RAFi+EGFRi")
    ax.axhspan(0.25, 0.5, xmin=0.45, xmax=0.55, color=ORANGE, alpha=0.25)
    ax.set(xlabel="hours", ylabel="pERK (fraction of untreated)", ylim=(0, 1.05),
           title=f"(i) ERK rebound  {'PASS' if b['i_erk_rebound']['pass'] else 'FAIL'}")
    ax.legend(loc="upper left", frameon=False)
    cs, P_kras, P_hea, P_mel = curves["ii"]
    ax = axes[0, 1]
    for X, col, lab in ((P_kras, ORANGE, "KRAS-mutant tumour"), (P_hea, GREY, "healthy (RAS/RAF wt)"), (P_mel, BLUE, "BRAF-mutant mel")):
        ax.plot(cs, np.median(X, 0), color=col, lw=2.5, label=lab)
        ax.fill_between(cs, np.percentile(X, 25, 0), np.percentile(X, 75, 0), color=col, alpha=0.18, lw=0)
    ax.axhline(1, color="black", lw=1)
    ax.set_xscale("log")
    ax.set(xlabel="RAFi (µM)", ylabel="pERK (fold of untreated), first hour",
           title=f"(ii) Paradoxical activation  {'PASS' if b['ii_paradox']['pass'] else 'FAIL'}")
    ax.legend(loc="lower left", frameon=False)
    t, Em, Rm, Am, Erm = curves["iii"]
    ax = axes[0, 2]
    band(ax, t, Em, BLUE, "pERK, MEKi"); band(ax, t, Rm, ORANGE, "RAS-GTP (receptor), MEKi"); band(ax, t, Am, VIOLET, "pAKT, MEKi")
    band(ax, t, Erm, AQUA, "pERK, RAFi+MEKi")
    ax.axhline(1, color="black", lw=1)
    ax.set(xlabel="hours", ylabel="fold of untreated", ylim=(0, None),
           title=f"(iii) Feedback relief, BRAF mel  {'PASS' if b['iii_feedback_relief']['pass'] else 'FAIL'}")
    ax.legend(loc="upper right", frameon=False, fontsize=10)
    ce, V = curves["iv"]
    ax = axes[1, 0]
    for (k, X), col in zip(V.items(), (BLUE, ORANGE, AQUA)):
        ax.plot(ce, np.median(X, 0), color=col, lw=2.5, marker="o", ms=6, label=k)
        ax.fill_between(ce, np.percentile(X, 25, 0), np.percentile(X, 75, 0), color=col, alpha=0.18, lw=0)
    ax.set_xscale("log")
    ax.set(xlabel="EGFRi (µM)", ylabel="72 h viability", ylim=(0, 1.1),
           title=f"(iv) MET bypass, EGFR lung  {'PASS' if b['iv_met_bypass']['pass'] else 'FAIL'}")
    ax.legend(loc="lower left", frameon=False)
    shares, D = curves["v"]
    ax = axes[1, 1]
    med = np.array([np.nanmedian(D[:, j]) if np.isfinite(D[:, j]).any() else np.nan for j in range(D.shape[1])])
    ax.plot(shares * 100, med, color=ORANGE, lw=2.5, marker="o", ms=8)
    ax.fill_between(shares * 100, np.nanpercentile(D, 25, 0), np.nanpercentile(D, 75, 0), color=ORANGE, alpha=0.18, lw=0)
    ax.axhline(M.TREAT_DAYS, color=GREY, ls="--", lw=1.5)
    ax.text(shares[0] * 100, M.TREAT_DAYS + 2, "treatment ends", color=GREY, fontsize=11)
    ax.set_xscale("log")
    ax.set(xlabel="initial KRAS-mutant clone share (%)", ylabel="days to regrowth (+20 % from nadir)",
           title=f"(v) KRAS clone under EGFRi  {'PASS' if b['v_kras_outgrowth']['pass'] else 'FAIL'}")
    ax = axes[1, 2]
    ax.axis("off")
    lines = ["Sources", "(i) Corcoran 2012 Cancer Discov;", "     Prahallad 2012 Nature",
             "(ii) Poulikakos 2010; Hatzivassiliou", "     2010 (Nature)", "(iii) Lito 2012 Cancer Cell;",
             "     Turke 2012 Cancer Res (pAKT)", "(iv) Engelman 2007 Science;", "     Turke 2010 Cancer Cell",
             "(v) Misale 2012; Diaz 2012 (Nature)", "", "Bands: IQR over virtual tumours", "Drugs at model css (label dose)"]
    ax.text(0, 1, "\n".join(lines), va="top", fontsize=12, family="sans-serif")
    fig.tight_layout()
    fig.savefig(os.path.join(IMG, "val_mechanisms.png"))
    plt.close(fig)


ARCH_LABEL = {"BRAF_MEL": "BRAF mel\nRAFi+MEKi", "BRAF_CRC": "BRAF CRC\nRAFi+MEKi+EGFRi", "KRAS": "KRAS G12C CRC\nKRASi+EGFRi",
              "EGFR": "EGFR lung\nEGFRi"}


def fig_c(c, plt, lattice_reps):
    src = "lattice" if lattice_reps else "mean_model"
    fig, axes = plt.subplots(1, 2, figsize=(16, 6.5), gridspec_kw={"width_ratios": [1.2, 1]})
    ax = axes[0]
    order = ["BRAF_MEL", "BRAF_CRC", "KRAS", "EGFR"]
    x = np.arange(len(order))
    sim30 = [100 * c[a][src]["orr_like_30pct_burden"] for a in order]
    sim51 = [100 * c[a][src]["orr_like_51pct_area"] for a in order]
    tr = [c[a]["trial_orr_pct"] for a in order]
    trlo = [c[a]["trial_orr_range_pct"][0] for a in order]
    trhi = [c[a]["trial_orr_range_pct"][1] for a in order]
    w = 0.27
    ax.bar(x - w, sim30, w * 0.92, color=BLUE, label="sim: burden −30 %")
    ax.bar(x, sim51, w * 0.92, color=BLUE, alpha=0.45, label="sim: area −51 % (≈ diameter −30 %)")
    ax.bar(x + w, tr, w * 0.92, color=ORANGE, label="pivotal trial ORR")
    ax.errorbar(x + w, tr, yerr=[np.array(tr) - trlo, np.array(trhi) - tr], fmt="none", ecolor="black", capsize=4)
    for xi, v in zip(x, sim30):
        ax.text(xi - w, v + 1.5, f"{v:.0f}", ha="center", fontsize=11)
    for xi, v in zip(x, tr):
        ax.text(xi + w, v + 1.5, f"{v:.0f}", ha="center", fontsize=11)
    ax.set_xticks(x)
    ax.set_xticklabels([ARCH_LABEL[a] for a in order], fontsize=11)
    ax.set(ylabel="% of patients", ylim=(0, 105), title="Response-like rate, label-dose standard regimen")
    ax.legend(frameon=False, loc="upper right", fontsize=10)
    ax = axes[1]
    if lattice_reps:
        days = np.arange(0, M.FOLLOW_DAYS + 1, 7)
        for a, col in zip(order, (BLUE, ORANGE, AQUA, VIOLET)):
            X = c[a]["_curves"]
            ax.plot(days[:X.shape[1]], 100 * np.median(X, 0), color=col, lw=2.5, label=a)
            ax.fill_between(days[:X.shape[1]], 100 * np.percentile(X, 25, 0), 100 * np.percentile(X, 75, 0), color=col, alpha=0.12, lw=0)
        ax.axvline(M.TREAT_DAYS, color=GREY, ls="--", lw=1.5)
        ax.text(M.TREAT_DAYS + 3, 5, "treatment ends", color=GREY, fontsize=11)
        ax.set(xlabel="day", ylabel="lesion tumour cells (% of day 0)", yscale="symlog", ylim=(0, 600),
               title="Tumour burden on the lattice (median, IQR)")
        ax.legend(frameon=False, loc="upper left", fontsize=11)
    fig.tight_layout()
    fig.savefig(os.path.join(IMG, "val_clinical_orr.png"))
    plt.close(fig)


def fig_e(e, plt):
    if "by_type" not in e:
        return
    order = ["BRAF_MEL", "BRAF_CRC", "KRAS", "EGFR"]
    cats = [("eligible (5 guidelines, approx.)", "#2a78d6"), ("eligible only with 3 guidelines", "#9cc0ec"),
            ("a guideline regimen works", ORANGE), ("no grid regimen cures", GREY), ("fails lattice checks", "#d8d7d2")]
    fig, ax = plt.subplots(figsize=(12, 5.5))
    left = np.zeros(len(order))
    for (k, col) in cats:
        vals = np.array([100 * e["by_type"][a]["share"].get(k, 0.0) for a in order])
        ax.barh(np.arange(len(order)), vals, left=left, color=col, edgecolor="white", linewidth=2, label=k)
        for yi, (l0, v) in enumerate(zip(left, vals)):
            if v >= 6:
                ax.text(l0 + v / 2, yi, f"{v:.0f}", ha="center", va="center", fontsize=11,
                        color="white" if col in ("#2a78d6", GREY, ORANGE) else "black")
        left += vals
    ax.set_yticks(np.arange(len(order)))
    ax.set_yticklabels([f"{a} (n={e['by_type'][a]['n']})" for a in order])
    ax.set(xlabel="% of generated patients", xlim=(0, 100), title="Who enters the benchmark cohort (development scan, 1600 patients)")
    ax.legend(frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.15), ncol=3, fontsize=10)
    ax.grid(axis="y", visible=False)
    fig.tight_layout()
    fig.savefig(os.path.join(IMG, "val_cohort.png"))
    plt.close(fig)


def strip(o):
    if isinstance(o, dict):
        return {k: strip(v) for k, v in o.items() if not k.startswith("_")}
    if isinstance(o, list):
        return [strip(v) for v in o]
    if isinstance(o, (np.floating,)):
        return float(o)
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, float):
        return round(o, 4)
    return o


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=600)
    ap.add_argument("--procs", type=int, default=16)
    ap.add_argument("--lattice-reps", type=int, default=2)
    a = ap.parse_args()
    os.makedirs(IMG, exist_ok=True)
    t0 = time.time()
    with Pool(min(a.procs, 16)) as pool:
        res = pool.map(work, [(i, a.lattice_reps) for i in range(a.n)], chunksize=4)
    t1 = time.time()
    patients = [P.make_patient(r["seed"]) for r in res]
    A = part_a(res)
    B, curves = part_b(patients)
    Cc = part_c(res, a.lattice_reps)
    D = part_d()
    E = part_e()
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    _style(plt)
    fig_a(A, plt)
    fig_b(curves, B, plt)
    fig_c(Cc, plt, a.lattice_reps)
    fig_e(E, plt)
    from collections import Counter
    out = {"command": "science/.venv/bin/python science/validation/validate_model.py",
           "population": {"salt": "validation", "n": a.n, "by_type": dict(Counter(r["arch"] for r in res)),
                          "lattice_reps": a.lattice_reps, "selection": "none (unselected generator output)"},
           "a_gdsc": A, "b_mechanisms": B, "c_clinical": Cc, "d_pk": D, "e_cohort": E,
           "runtime_s": {"population": round(t1 - t0), "total": round(time.time() - t0)}}
    json.dump(strip(out), open(os.path.join(DOCS, "validation.json"), "w"), indent=1)
    print(json.dumps(strip({"a": {d: [(r["sim_group"], r.get("delta_log10_ic50_median")) for r in A[d]] for d in A},
                            "b": {k: v["pass"] for k, v in B.items()},
                            "c": {k: (v.get("lattice", v["mean_model"])["orr_like_30pct_burden"], v["trial_orr_pct"]) for k, v in Cc.items()}}), indent=1))
    print("done in %.0f s" % (time.time() - t0))


if __name__ == "__main__":
    main()
