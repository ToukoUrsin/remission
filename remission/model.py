"""Remission: the public cell model.

This file is the whole biology of the environment. The lab's wells and the
patient's tissue are both simulated with exactly these equations; what the
agent does not know is one patient's parameter values. A copy of this file is
in the agent's container as /app/cellmodel.py.

Scope: human tumour cells and matched healthy cells, eight targeted-drug
classes, receptor -> MAPK / PI3K signalling with ERK-dependent feedback,
apoptotic priming, clonal heterogeneity, one-compartment pharmacokinetics and
an on-lattice tissue. See docs/MODEL.md for the equations and the literature
each term comes from.

Everything is plain numpy. Concentrations are free drug in micromolar, time is
in days, rates are per day.
"""
from __future__ import annotations

import numpy as np

DRUGS = ("RAFi", "MEKi", "EGFRi", "METi", "KRASi", "PI3Ki", "BCLi", "MCLi")
ND = len(DRUGS)
IDX = {d: i for i, d in enumerate(DRUGS)}
iRAF, iMEK, iEGFR, iMET, iKRAS, iPI3K, iBCL, iMCL = range(ND)

# ---------------------------------------------------------------------------
# Fixed constants of the model family (the same for every patient).
# ---------------------------------------------------------------------------
ARAS = 0.25   # receptor -> RAS-GTP gain
KE = 0.60     # half-saturation of the RAF/MEK -> ERK cascade
NE = 1.5      # steepness of the cascade
AP = 0.35     # receptor -> PI3K gain
WE = 0.50     # share of EGFR-family signal that reaches PI3K (the bypass receptor couples fully)
KA = 0.50     # half-saturation of PI3K -> AKT
GE = 0.12     # ERK level for half-maximal proliferation
G0 = 0.40     # proliferation left without AKT
GA = 0.30     # AKT level for half of the AKT-dependent proliferation
HE = 0.50     # ERK level for half-maximal survival signalling
HA = 0.20     # AKT level for half-maximal survival signalling
MH = 8.0      # steepness of the apoptosis switch
TAU_F = 0.75  # days; time constant of feedback relief (receptor re-expression)
KD0 = 0.01    # per day; death unrelated to apoptotic priming
PX = 0.8      # paradoxical activation of wild-type RAF dimers at full RAF-inhibitor occupancy

# Pharmacokinetics at dose 1.0 (= the labelled dose), once daily. One
# compartment, first-order elimination. `css` is the average steady-state
# concentration at the tumour in micromolar, in the units of the wells (what a
# well at that concentration does, the tissue at that level does too).
# Half-lives follow the drug labels of one representative per class; see
# docs/MODEL.md and science/pk_table.json.
PK = {
    "RAFi":  {"t_half_h": 8.0,   "css": 0.30},
    "MEKi":  {"t_half_h": 96.0,  "css": 0.012},
    "EGFRi": {"t_half_h": 36.0,  "css": 0.24},
    "METi":  {"t_half_h": 6.5,   "css": 0.060},
    "KRASi": {"t_half_h": 24.0,  "css": 0.50},
    "PI3Ki": {"t_half_h": 9.0,   "css": 0.50},
    "BCLi":  {"t_half_h": 17.0,  "css": 0.175},
    "MCLi":  {"t_half_h": 12.0,  "css": 0.075},
}
CSS_LABEL = np.array([PK[d]["css"] for d in DRUGS])
DOSE_MAX = 2.0          # highest allowed daily dose, as a multiple of the labelled dose
TREAT_DAYS = 84         # the regimen may dose on days 0..83
FOLLOW_DAYS = 364       # the tissue is read one year after the start
WASHOUT_DAYS = 42       # drug levels are followed this long after the last dosing day; then they are zero
SHED = 0.06             # per day; normal tissue is a renewing tissue that sheds this share of its cells
DT = 0.25               # days; integration step everywhere

# Tissue lattice
LATTICE = 64
TUMOUR_SHARE0 = 0.22    # share of sites that are tumour on day 0
EMPTY_SHARE0 = 0.04
P_INVADE = 0.50         # chance that a dividing tumour cell displaces a healthy neighbour

# Wells
SEED_CELLS = 3000       # cells seeded per well
CAP_TUMOUR = 12.0       # well capacity as a multiple of the seeded count (tumour and derived lines)
CAP_HEALTHY = 5.0       # healthy cells stop earlier (contact inhibition)

CLONE_KEYS = ("braf", "ras", "pi3k", "pten", "egfr", "alt", "fb",
              "kp", "kdmax", "bcl", "bfrac", "mfrac", "sE", "sA")



def make_clone(**kw) -> dict:
    """A cell population. Defaults are a healthy epithelial cell."""
    c = dict(braf=0.0, ras=0.0, pi3k=0.0, pten=1.0, egfr=1.0, alt=0.1, fb=3.0,
             kp=0.35, kdmax=1.2, bcl=0.9, bfrac=0.4, mfrac=0.3, sE=0.7, sA=0.7,
             kshift=np.ones(ND), tox=np.zeros(ND), healthy=False, name="clone")
    c.update(kw)
    c["kshift"] = np.asarray(c["kshift"], dtype=float)
    c["tox"] = np.asarray(c["tox"], dtype=float)
    return c


def pack(clones: list[dict]) -> dict:
    """Stack clone dicts into arrays of shape [C] (kshift: [C, ND])."""
    p = {k: np.array([float(c[k]) for c in clones]) for k in CLONE_KEYS}
    p["kshift"] = np.stack([np.asarray(c["kshift"], dtype=float) for c in clones])
    p["tox"] = np.stack([np.asarray(c.get("tox", np.zeros(ND)), dtype=float) for c in clones])
    p["healthy"] = np.array([bool(c.get("healthy", False)) for c in clones])
    return p


def occupancy(conc, K, kshift):
    """Target occupancy of each drug. conc [..., ND], K [ND], kshift [..., ND]."""
    conc = np.maximum(np.asarray(conc, dtype=float), 0.0)
    return conc / (conc + K * kshift)


def signalling(p: dict, occ, F):
    """ERK and AKT activity (0..1) for clones `p` at occupancy `occ` and feedback state `F`.

    Shapes broadcast: p[...] with occ[..., ND] and F[...].
    """
    Re = p["egfr"] * F * (1.0 - occ[..., iEGFR])          # EGFR-family receptor output
    Ra = p["alt"] * F * (1.0 - occ[..., iMET])            # bypass receptor (MET), not hit by EGFRi
    ras = p["ras"] * (1.0 - occ[..., iKRAS]) + ARAS * (Re + Ra)   # RAS-GTP: mutant KRAS + receptor-driven
    U = p["braf"] * (1.0 - occ[..., iRAF]) + ras * (1.0 + PX * occ[..., iRAF])
    M = U * (1.0 - occ[..., iMEK])
    Mn = np.power(np.maximum(M, 0.0), NE)
    E = Mn / (Mn + KE ** NE)
    P = (p["pi3k"] + AP * (WE * Re + Ra)) * (1.0 - occ[..., iPI3K]) / p["pten"]
    A = P / (P + KA)
    return E, A


def baseline(p: dict):
    """Untreated activity (feedback state 1)."""
    z = np.zeros(p["kp"].shape + (ND,))
    return signalling(p, z, np.ones(p["kp"].shape))


def feedback_target(p: dict, E, E0):
    """Where the feedback state relaxes to. 1 at baseline, up to 1 + fb when ERK is off."""
    x = E / np.maximum(E0, 1e-9)
    return (1.0 + p["fb"]) / (1.0 + p["fb"] * x * x)


def fate_rates(p: dict, E, A, occ, conc=None):
    """Proliferation and death rate (per day).

    conc [..., ND] (broadcast against the clone axis) adds off-target toxicity:
    every drug kills the patient's normal cells at a rate proportional to its
    concentration (clone["tox"], per day at the labelled average concentration;
    zero for tumour clones).
    """
    gE = E * E / (E * E + GE * GE)
    gA = G0 + (1.0 - G0) * A / (A + GA)
    kp = p["kp"] * gE * gA
    hE = E / (E + HE)
    hA = A / (A + HA)
    D = (1.0 - p["sE"] * hE) * (1.0 - p["sA"] * hA)          # death drive
    # anti-apoptotic buffer left: a share bfrac of it is BCL-2/BCL-xL, a share mfrac MCL-1,
    # the rest is out of reach of both drugs
    B = p["bcl"] * (1.0 - p["bfrac"] * occ[..., iBCL] - p["mfrac"] * occ[..., iMCL])
    Dm = np.power(D, MH)
    kd = KD0 + p["kdmax"] * Dm / (Dm + np.power(np.maximum(B, 1e-6), MH))
    if conc is not None:
        kd = kd + (p["tox"] * (np.asarray(conc, dtype=float) / CSS_LABEL)).sum(axis=-1)
    return kp, kd


def step_feedback(p, F, E, E0, dt=DT):
    return F + (feedback_target(p, E, E0) - F) * (1.0 - np.exp(-dt / TAU_F))


# ---------------------------------------------------------------------------
# What is known before any experiment (public; the same for every patient).
# A patient's own values are hidden. They lie inside these ranges.
# ---------------------------------------------------------------------------
# Typical potency of each drug (K, micromolar: the concentration that occupies half of the target).
# Between patients K varies about 2.2-fold (1 sd of ln K = 0.8). The same K holds for all of a
# patient's cells, times the clone's `kshift` for that drug.
K_TYPICAL = {"RAFi": 0.050, "MEKi": 0.0040, "EGFRi": 0.060, "METi": 0.020, "KRASi": 0.10, "PI3Ki": 0.20,
             "BCLi": 0.070, "MCLi": 0.030}
K_LN_SD = 0.8
# Off-target toxicity of each drug to the patient's healthy cells: extra death rate (per day) at the
# concentration `css`. Median 0.02, and it varies a lot between patients and drugs (1 sd of ln = 1.2).
TOX_TYPICAL, TOX_LN_SD = 0.02, 1.2

# Healthy cells. Receptor levels are in units of the healthy cell's EGFR.
HEALTHY_FIXED = dict(braf=0.0, ras=0.0, pi3k=0.0, pten=1.0, egfr=1.0, alt=0.3)
HEALTHY_RANGE = dict(fb=(1.0, 5.0), kp=(0.24, 0.50), kdmax=(0.6, 0.9), bfrac=(0.15, 0.60), mfrac=(0.15, 0.60),
                     sE=(0.45, 0.65), sA=(0.45, 0.65),
                     margin=(1.6, 2.6))      # bcl = margin x the cell's untreated death drive: not primed

# The dominant tumour clone. `margin` as above: tumour cells are primed for apoptosis.
TUMOUR_RANGE = dict(egfr=(0.005, 6.0), alt=(0.005, 1.5), fb=(0.5, 6.0), kp=(0.50, 0.90), kdmax=(0.9, 1.5),
                    bfrac=(0.10, 0.40), mfrac=(0.10, 0.40), sE=(0.70, 0.95), sA=(0.25, 0.90), margin=(1.30, 2.00))

# What each alteration in a sequencing report does to a clone. (parameter, low, high, how)
#   set: the parameter takes a value in [low, high]      mul: it is multiplied by a value in [low, high]
#   add: a value in [low, high] is added (sA is capped at 0.9)      kshift:<drug>: K of that drug is multiplied
ALTERATIONS = {
    "BRAF V600E": [("braf", 1.0, 1.0, "set")],
    "KRAS G12C": [("ras", 0.7, 1.1, "set")],
    "EGFR L858R": [("egfr", 2.5, 5.0, "set"), ("EGFRi", 0.1, 0.1, "kshift")],   # the inhibitor binds the mutant 10x better
    "PIK3CA H1047R": [("pi3k", 0.4, 0.8, "set"), ("sA", 0.2, 0.2, "add")],
    "PTEN loss": [("pten", 0.20, 0.35, "set"), ("sA", 0.2, 0.2, "add")],
    "MAP2K1 (MEK1) K57N": [("MEKi", 15.0, 40.0, "kshift")],
    "MET amplification": [("alt", 0.6, 1.5, "set")],
    "EGFR T790M": [("EGFRi", 3.0, 10.0, "kshift")],
    "MCL1 amplification": [("bcl", 1.6, 2.4, "mul"), ("bfrac", 0.05, 0.05, "set"), ("mfrac", 0.80, 0.90, "set")],
    "BCL2L1 (BCL-xL) amplification": [("bcl", 1.6, 2.4, "mul"), ("bfrac", 0.80, 0.90, "set"), ("mfrac", 0.05, 0.05, "set")],
}
# A minor clone is the dominant clone plus one more alteration, at 0.4 % to 6 % of the tumour cells.
MINOR_SHARE = (0.004, 0.06)


def death_drive0(clone: dict) -> float:
    """Untreated death drive D of a clone (see fate_rates). bcl = margin x this."""
    p = pack([clone])
    E, A = baseline(p)
    return float(((1.0 - p["sE"] * E / (E + HE)) * (1.0 - p["sA"] * A / (A + HA)))[0])


# ---------------------------------------------------------------------------
# Pharmacokinetics
# ---------------------------------------------------------------------------
def pk_constants():
    k = np.array([np.log(2.0) / (PK[d]["t_half_h"] / 24.0) for d in DRUGS])   # per day
    css = np.array([PK[d]["css"] for d in DRUGS])
    bolus = css * k * 1.0          # concentration added by a dose of 1.0 given once a day
    return k, bolus


def expand_regimen(regimen: dict) -> np.ndarray:
    """Regimen -> planned daily doses, shape [TREAT_DAYS, ND].

    regimen = {"phases": [{"days": 14, "doses": {"RAFi": 1.0, "MEKi": 0.5}}, ...]}
    Phases run back to back from day 0. Drugs not named in a phase get 0. Days
    after the last phase, and all days from TREAT_DAYS on, are untreated.
    """
    plan = np.zeros((TREAT_DAYS, ND))
    day = 0
    phases = regimen.get("phases", []) if isinstance(regimen, dict) else []
    for ph in phases:
        n = int(ph.get("days", 0))
        if n <= 0:
            continue
        doses = np.zeros(ND)
        for name, v in (ph.get("doses") or {}).items():
            if name in IDX:
                doses[IDX[name]] = min(max(float(v), 0.0), DOSE_MAX)
        plan[day:min(day + n, TREAT_DAYS)] = doses
        day += n
        if day >= TREAT_DAYS:
            break
    return plan


def validate_regimen(regimen) -> str | None:
    """None if the regimen is well formed, else what is wrong with it."""
    if not isinstance(regimen, dict) or not isinstance(regimen.get("phases"), list):
        return 'regimen must be {"phases": [{"days": int, "doses": {drug: dose}}, ...]}'
    if not 1 <= len(regimen["phases"]) <= 12:
        return "1 to 12 phases"
    total = 0
    for ph in regimen["phases"]:
        if not isinstance(ph, dict):
            return "each phase is an object"
        try:
            n = int(ph.get("days", 0))
        except (TypeError, ValueError):
            return "days must be an integer"
        if n < 1:
            return "each phase lasts at least 1 day"
        total += n
        doses = ph.get("doses") or {}
        if not isinstance(doses, dict):
            return "doses must be an object"
        for name, v in doses.items():
            if name not in IDX:
                return f"unknown drug {name!r}; drugs are {', '.join(DRUGS)}"
            try:
                v = float(v)
            except (TypeError, ValueError):
                return f"dose of {name} must be a number"
            if not (0.0 <= v <= DOSE_MAX) or v != v:
                return f"dose of {name} must be between 0 and {DOSE_MAX}"
    if total > TREAT_DAYS:
        return f"phases add up to {total} days; the limit is {TREAT_DAYS}"
    return None


def pk_profile(plan: np.ndarray, days: int = FOLLOW_DAYS, dt: float = DT) -> np.ndarray:
    """Step-average concentrations [days/dt, ND] for a daily dose plan."""
    k, bolus = pk_constants()
    spd = int(round(1.0 / dt))
    decay = np.exp(-k * dt)
    mean_factor = (1.0 - decay) / (k * dt)
    c = np.zeros(ND)
    out = np.zeros((days * spd, ND))
    for s in range(days * spd):
        d = s // spd
        if s % spd == 0 and d < plan.shape[0]:
            c = c + plan[d] * bolus
        out[s] = c * mean_factor
        c = c * decay
    return out


# ---------------------------------------------------------------------------
# Wells (in vitro). Constant concentrations, stochastic births and deaths.
# ---------------------------------------------------------------------------
def simulate_wells(K, clones: list[dict], frac0, conc, days: float, rng=None,
                   read_days=None, seed_cells: int = SEED_CELLS, cap_mult: float = CAP_TUMOUR,
                   dt: float = DT):
    """Grow W wells.

    K        [ND] drug potencies of the patient (micromolar)
    clones   list of C clone dicts
    frac0    [C] or [W, C] starting share of each clone
    conc     [W, ND] constant concentrations
    days     how long the wells run
    rng      numpy Generator for stochastic growth, or None for the mean
    read_days  times (days) at which to return counts; default: every step

    Returns dict(t [T], N [T, W, C] cell counts, E [T, W, C], A [T, W, C]).
    E and A are relative to the clone's own untreated level.
    """
    p = pack(clones)
    conc = np.atleast_2d(np.asarray(conc, dtype=float))
    W, C = conc.shape[0], len(clones)
    occ = occupancy(conc[:, None, :], np.asarray(K, dtype=float), p["kshift"][None, :, :])   # [W, C, ND]
    E0, A0 = baseline(p)
    frac0 = np.broadcast_to(np.asarray(frac0, dtype=float), (W, C))
    if rng is None:
        N = seed_cells * frac0.copy()
    else:
        N = rng.poisson(seed_cells * frac0).astype(float)
    F = np.ones((W, C))
    cap = seed_cells * cap_mult
    nstep = int(round(days / dt))
    if read_days is None:
        read_steps = np.arange(nstep + 1)
    else:
        read_steps = np.clip(np.round(np.asarray(read_days, dtype=float) / dt).astype(int), 0, nstep)
    want = {int(s) for s in read_steps}
    rec = {}
    for s in range(nstep + 1):
        E, A = signalling(p, occ, F)
        if s in want:
            rec[s] = (N.copy(), E / E0, A / A0)
        if s == nstep:
            break
        kp, kd = fate_rates(p, E, A, occ, conc[:, None, :])
        crowd = np.clip(1.0 - N.sum(axis=1, keepdims=True) / cap, 0.0, 1.0)
        if rng is None:
            N = N * np.exp((kp * crowd - kd) * dt)
        else:
            births = rng.poisson(N * (np.expm1(kp * crowd * dt)))
            deaths = rng.binomial(N.astype(np.int64), 1.0 - np.exp(-kd * dt))
            N = np.maximum(N + births - deaths, 0.0)
        F = step_feedback(p, F, E, E0, dt)
    t = np.array([s * dt for s in read_steps])
    Ns = np.stack([rec[int(s)][0] for s in read_steps])
    Es = np.stack([rec[int(s)][1] for s in read_steps])
    As = np.stack([rec[int(s)][2] for s in read_steps])
    return {"t": t, "N": Ns, "E": Es, "A": As}


# ---------------------------------------------------------------------------
# The patient. Three pieces of the same biology:
#   regimen_course        drug levels, signalling and per-clone rates over the course (deterministic)
#                         and the patient's normal tissue
#   simulate_tissue       the on-lattice stochastic lesion the regimen is scored on
#   simulate_tissue_mean  a fast well-mixed approximation of the lesion, for planning
# ---------------------------------------------------------------------------
_DIRS = [(-1, -1), (-1, 0), (-1, 1), (0, -1), (0, 1), (1, -1), (1, 0), (1, 1)]


def regimen_course(K, clones: list[dict], plans, days: int = FOLLOW_DAYS, dt: float = DT, full: bool = True):
    """Everything about a course that does not depend on cell numbers.

    clones[0] must be the patient's healthy cell type. plans [Q, TREAT_DAYS, ND].

    The patient's normal tissue (healthy cells elsewhere in the body, no tumour)
    is followed as one well-mixed compartment that starts at 1.0. It is a
    renewing tissue: it sheds SHED of its cells per day and replaces them by
    division. `harm` is its largest fractional loss.

    Returns dict(conc [Q, S, ND], kp [Q, S, C], kd [Q, S, C], E [Q, S, C], A [Q, S, C],
    normal [Q, days+1], harm [Q], kp0 [C], kd0 [C], steps S). S covers the days with drug on
    board (TREAT_DAYS + WASHOUT_DAYS); after that every clone has its untreated rates kp0, kd0
    (conc, E, A only if `full`).
    """
    plans = np.asarray(plans, dtype=float)
    if plans.ndim == 2:
        plans = plans[None]
    Q = plans.shape[0]
    p = pack(clones)
    C = len(clones)
    K = np.asarray(K, dtype=float)
    E0, A0 = baseline(p)
    kp0, kd0 = fate_rates(p, E0, A0, np.zeros((C, ND)))
    free0 = min(max((kd0[0] + SHED) / kp0[0], 1e-6), 1.0)     # free share of normal tissue at homeostasis
    kel, bolus = pk_constants()
    spd = int(round(1.0 / dt))
    active = min(days, TREAT_DAYS + WASHOUT_DAYS)             # after this: no drug left, untreated rates
    S = active * spd
    decay = np.exp(-kel * dt)
    mean_factor = (1.0 - decay) / (kel * dt)
    conc = np.zeros((Q, ND))
    F = np.ones((Q, C))
    h = np.ones(Q)
    o_kp = np.zeros((Q, S, C))
    o_kd = np.zeros((Q, S, C))
    if full:
        o_conc = np.zeros((Q, S, ND))
        o_E = np.zeros((Q, S, C))
        o_A = np.zeros((Q, S, C))
    normal = np.ones((Q, days + 1))
    kshift = p["kshift"][None, :, :]
    for day in range(active):
        if day < TREAT_DAYS:
            conc = conc + plans[:, day] * bolus
        for j in range(spd):
            s = day * spd + j
            cmean = conc * mean_factor
            occ = occupancy(cmean[:, None, :], K, kshift)
            E, A = signalling(p, occ, F)
            kp, kd = fate_rates(p, E, A, occ, cmean[:, None, :])
            o_kp[:, s] = kp
            o_kd[:, s] = kd
            if full:
                o_conc[:, s] = cmean
                o_E[:, s] = E / E0
                o_A[:, s] = A / A0
            F = step_feedback(p, F, E, E0, dt)
            conc = conc * decay
            room = np.clip(1.0 - h * (1.0 - free0), 0.0, 1.0)
            h = h * np.exp((kp[:, 0] * room - kd[:, 0] - SHED) * dt)
        normal[:, day + 1] = h
    for day in range(active, days):                           # recovery of the normal tissue, no drug
        for j in range(spd):
            room = np.clip(1.0 - h * (1.0 - free0), 0.0, 1.0)
            h = h * np.exp((kp0[0] * room - kd0[0] - SHED) * dt)
        normal[:, day + 1] = h
    harm = np.clip(1.0 - normal.min(axis=1), 0.0, 1.0)
    out = {"kp": o_kp, "kd": o_kd, "kp0": kp0, "kd0": kd0, "steps": S, "normal": normal, "harm": harm}
    if full:
        out.update(conc=o_conc, E=o_E, A=o_A)
    return out


def init_tissue(frac_tumour, rng, reps: int, L: int = LATTICE):
    """Lattices [reps, L, L]: 0 empty, 1 healthy, 2+i tumour clone i.

    The tumour starts as one lesion in the middle; its clones are mixed inside it.
    """
    frac_tumour = np.asarray(frac_tumour, dtype=float)
    frac_tumour = frac_tumour / frac_tumour.sum()
    yy, xx = np.mgrid[0:L, 0:L]
    grid = np.ones((reps, L, L), dtype=np.int8)
    r_les = np.sqrt(TUMOUR_SHARE0 * L * L / np.pi)
    for r in range(reps):
        cy, cx = L / 2 + rng.normal(0, 2, 2)
        ang = np.arctan2(yy - cy, xx - cx)       # a ragged lesion edge
        rad = r_les * (1 + 0.12 * np.sin(3 * ang + rng.uniform(0, 6.3)) + 0.08 * np.sin(5 * ang + rng.uniform(0, 6.3)))
        inside = (yy - cy) ** 2 + (xx - cx) ** 2 <= rad ** 2
        n_in = int(inside.sum())
        labels = 2 + rng.choice(len(frac_tumour), size=n_in, p=frac_tumour)
        g = grid[r]
        g[inside] = labels.astype(np.int8)
        empty = rng.random((L, L)) < EMPTY_SHARE0
        g[empty & ~inside] = 0
    return grid


def simulate_tissue(K, healthy: dict, tumour: list[dict], frac_tumour, regimen: dict, seed: int = 0,
                    reps: int = 8, L: int = LATTICE, days: int = FOLLOW_DAYS, record_every: int | None = None,
                    dt: float = DT):
    """Treat the patient with `regimen` and follow `reps` copies of the lesion to `days`.

    Returns dict:
      healthy [reps, days+1], tumour [reps, days+1, n_clones]  cell counts in the lesion each day
      normal [days+1]   the patient's normal tissue (1.0 = day 0)
      conc [days, ND]   mean drug concentration each day
      healed [reps]     healthy share of the lesion's sites on the last day
      harm              largest fractional loss of normal tissue
      reward [reps]     max(0, healed - harm)
      frames            (if record_every) list of (day, grid copy of lesion 0)
    """
    rng = np.random.default_rng(seed)
    clones = [healthy] + list(tumour)
    C = len(clones)
    course = regimen_course(K, clones, expand_regimen(regimen)[None], days, dt)
    kp, kd, S = course["kp"][0], course["kd"][0], course["steps"]
    pdie0 = np.concatenate([[0.0], 1.0 - np.exp(-course["kd0"] * dt)])
    pdiv0 = np.concatenate([[0.0], 1.0 - np.exp(-course["kp0"] * dt)])
    spd = int(round(1.0 / dt))
    grid = init_tissue(frac_tumour, rng, reps, L)
    sites = L * L
    H = np.zeros((reps, days + 1))
    T = np.zeros((reps, days + 1, C - 1))
    frames = []

    def count(day):
        H[:, day] = (grid == 1).sum(axis=(1, 2))
        for i in range(C - 1):
            T[:, day, i] = (grid == 2 + i).sum(axis=(1, 2))

    count(0)
    if record_every:
        frames.append((0, grid[0].copy()))
    for day in range(days):
        for j in range(spd):
            s = day * spd + j
            if s < S:
                pdie = np.concatenate([[0.0], 1.0 - np.exp(-kd[s] * dt)])   # by cell type; 0 = empty
                pdiv = np.concatenate([[0.0], 1.0 - np.exp(-kp[s] * dt)])
            else:
                pdie, pdiv = pdie0, pdiv0
            u = rng.random(grid.shape)
            grid[u < pdie[grid]] = 0
            u = rng.random(grid.shape)
            dividing = u < pdiv[grid]
            direction = rng.integers(0, 8, size=grid.shape)
            inv = rng.random(grid.shape) < P_INVADE
            for kdir in rng.permutation(8):
                dy, dx = _DIRS[kdir]
                src = dividing & (direction == kdir)
                if not src.any():
                    continue
                parent = np.where(src, grid, 0)
                parent_at_target = np.roll(parent, (dy, dx), axis=(1, 2))
                inv_at_target = np.roll(inv & src, (dy, dx), axis=(1, 2))
                ok = (parent_at_target > 0) & ((grid == 0) | ((grid == 1) & (parent_at_target >= 2) & inv_at_target))
                grid[ok] = parent_at_target[ok]
        count(day + 1)
        if record_every and (day + 1) % (record_every if day < TREAT_DAYS + WASHOUT_DAYS else 3 * record_every) == 0:
            frames.append((day + 1, grid[0].copy()))

    healed = H[:, -1] / sites
    harm = float(course["harm"][0])
    conc_day = np.zeros((days, ND))
    conc_day[:S // spd] = course["conc"][0].reshape(S // spd, spd, ND).mean(axis=1)
    out = {"healthy": H, "tumour": T, "normal": course["normal"][0], "conc": conc_day,
           "healed": healed, "harm": harm, "reward": np.maximum(0.0, healed - harm), "sites": sites}
    if record_every:
        out["frames"] = frames
    return out


def simulate_tissue_mean(K, healthy: dict, tumour: list[dict], frac_tumour, plans, days: int = FOLLOW_DAYS,
                         dt: float = DT, sites: int = LATTICE * LATTICE):
    """Well-mixed approximation of simulate_tissue for many dose plans at once.

    plans [Q, TREAT_DAYS, ND] daily doses (see expand_regimen). Returns dict of
    arrays over the Q plans: healed, harm, reward, cure (probability that every
    tumour clone goes extinct), kill (0..1, how far the hardest clone was pushed
    towards extinction) and the trajectories healthy [Q, days+1],
    tumour [Q, days+1, n_clones] in cells.

    Differences from the lattice: no space (tumour growth is not limited to the
    lesion edge), counts are expectations, and extinction is a probability:
    exp(-smallest expected count) per clone.
    """
    plans = np.asarray(plans, dtype=float)
    if plans.ndim == 2:
        plans = plans[None]
    Q = plans.shape[0]
    clones = [healthy] + list(tumour)
    C = len(clones)
    course = regimen_course(K, clones, plans, days, dt, full=False)
    KP, KD, S = course["kp"], course["kd"], course["steps"]
    kp_end, kd_end = course["kp0"][None, :], course["kd0"][None, :]
    spd = int(round(1.0 / dt))
    frac_tumour = np.asarray(frac_tumour, dtype=float)
    frac_tumour = frac_tumour / frac_tumour.sum()

    N = np.zeros((Q, C))
    N[:, 0] = sites * (1.0 - TUMOUR_SHARE0 - EMPTY_SHARE0)
    N[:, 1:] = sites * TUMOUR_SHARE0 * frac_tumour[None, :]
    N0 = N.copy()
    Hs = np.zeros((Q, days + 1))
    Ts = np.zeros((Q, days + 1, C - 1))
    Hs[:, 0] = N[:, 0]
    Ts[:, 0] = N[:, 1:]
    nmin = N[:, 1:].copy()
    for day in range(days):
        for j in range(spd):
            s = day * spd + j
            kp, kd = (KP[:, s], KD[:, s]) if s < S else (kp_end, kd_end)
            tot = N.sum(axis=1)
            empty = np.clip(1.0 - tot / sites, 0.0, 1.0)
            hshare = N[:, 0] / sites
            # a dividing cell needs room: healthy cells need an empty site,
            # tumour cells an empty site or (with P_INVADE) a healthy one
            room_t = empty + P_INVADE * hshare
            births_t = N[:, 1:] * kp[:, 1:] * room_t[:, None] * dt
            lost_h = P_INVADE * hshare / np.maximum(room_t, 1e-12) * births_t.sum(axis=1)
            dH = N[:, 0] * (kp[:, 0] * empty - kd[:, 0]) * dt - lost_h
            dT = births_t - N[:, 1:] * kd[:, 1:] * dt
            N[:, 0] = np.maximum(N[:, 0] + dH, 0.0)
            N[:, 1:] = np.maximum(N[:, 1:] + dT, 0.0)
            nmin = np.minimum(nmin, N[:, 1:])
        Hs[:, day + 1] = N[:, 0]
        Ts[:, day + 1] = N[:, 1:]
    cure_clone = np.exp(-nmin)                       # P(no cell of the clone survives its lowest point)
    cure = cure_clone.prod(axis=1)
    # how far each clone was pushed down, on a log scale from its start to 0.01 expected cells
    span = np.log(np.maximum(N0[:, 1:], 1e-9) / 0.01)
    kill = np.clip(np.log(np.maximum(N0[:, 1:], 1e-9) / np.maximum(nmin, 0.01)) / np.maximum(span, 1e-9), 0.0, 1.0).min(axis=1)
    harm = course["harm"]
    healed_if_relapse = Hs[:, -1] / sites
    healed_if_cure = _regrow_healthy(pack(clones), Hs, Ts, sites, days)
    healed = cure * healed_if_cure + (1.0 - cure) * healed_if_relapse
    reward = cure * np.maximum(0.0, healed_if_cure - harm) + (1.0 - cure) * np.maximum(0.0, healed_if_relapse - harm)
    return {"healed": healed, "harm": harm, "reward": reward, "cure": cure, "kill": kill,
            "healthy": Hs, "tumour": Ts, "normal": course["normal"]}


def _regrow_healthy(p, Hs, Ts, sites, days):
    """Healthy share on the last day if the tumour had gone extinct at its lowest point."""
    Q = Hs.shape[0]
    tsum = Ts.sum(axis=2)
    t_low = tsum.argmin(axis=1)
    E0, A0 = baseline(p)
    kp0, kd0 = fate_rates(p, E0, A0, np.zeros(p["kp"].shape + (ND,)))
    r = kp0[0]
    k_eff = max(1.0 - kd0[0] / r, 1e-6)
    h = np.minimum(Hs[np.arange(Q), t_low] / sites, k_eff)
    tt = np.maximum(days - t_low - 14, 0)            # a conservative 14-day lag for washout
    return k_eff / (1.0 + (k_eff / np.maximum(h, 1e-9) - 1.0) * np.exp(-(r - kd0[0]) * tt))
