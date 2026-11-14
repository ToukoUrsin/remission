"""Who is in the benchmark: patients that standard treatment fails, who can still be cured.

Hidden side. A generated patient joins the cohort when all of this holds
(checked on the fast model and again on the tissue lattice):

  1. refractory: each of the five guideline regimens for the tumour type
     scores below REFRACTORY_BELOW. Two of them are the label-dose standards of
     care. The other three are the constant regimens that cure the largest
     share of the patients the earlier ones leave uncured, picked one after
     the other by exhaustive search over the grid in a development population
     of 1600 generated patients;
  2. curable: at least one constant-dose regimen out of a grid of 3112 scores
     CURABLE_FROM or more;
  3. progressive: left untreated, the tumour takes the lesion (healthy share
     at one year at most UNTREATED_BELOW), so doing nothing earns nothing;
  4. no simple answer: every single drug and every pair of drugs at the label
     dose scores below REFRACTORY_BELOW, so a quick screen read off at label
     doses does not find the cure;
  5. a hidden clone decides: the tumour holds at least one resistant clone
     below the detection limit of sequencing, and the three regimens that
     would be best if the tumour were only what sequencing of the biopsy
     shows all score below REFRACTORY_BELOW on the real tumour;
  6. hedging does not work: take every resistant clone the tumour type can
     carry and sequencing did not show, one hypothetical tumour per clone, and
     search for the regimen that is best on average over all of them (with the
     patient's true drug potencies and true healthy cells, which no agent
     has). Done twice, for typical clones at 1 % and for the strongest clones
     at 6 %. Both regimens score below REFRACTORY_BELOW on the real tumour.

So a guess from the chart fails, a screen fails, the answer for the visible
tumour fails, covering everything that might be there fails, and a cure exists
for every patient. What is left is to find out which clone is hiding.
Functional precision oncology is used for exactly these patients: the ones
whose tumour comes back from cells nobody saw.

Rule 6 costs a few minutes of one core per patient that passes rules 1 to 5,
and about 1 generated patient in 400 passes all six.
"""
from __future__ import annotations

import itertools

import numpy as np

from . import model as M
from . import patients as P

REFRACTORY_BELOW = 0.10
CURABLE_FROM = 0.50
UNTREATED_BELOW = 0.02
VISIBLE_TOP = 3
SEQ_LOD = 0.02          # the detection limit of sequencing (lab.SEQ_LOD)
HEDGES = (("typical", 0.5, 0.01), ("strong", 1.0, 0.06))    # (name, strength within the generator's range, share)

GUIDELINES = {
    "BRAF_CRC": [{"RAFi": 1.0, "MEKi": 1.0, "EGFRi": 1.0},
                 {"RAFi": 1.0, "EGFRi": 1.0},
                 {"RAFi": 0.5, "MEKi": 0.5, "KRASi": 0.5, "MCLi": 0.5},
                 {"RAFi": 0.5, "MEKi": 0.5, "BCLi": 1.0},
                 {"RAFi": 0.5, "MEKi": 1.0, "MCLi": 2.0}],
    "BRAF_MEL": [{"RAFi": 1.0, "MEKi": 1.0},
                 {"RAFi": 1.0},
                 {"RAFi": 0.5, "MEKi": 0.5, "METi": 0.5, "BCLi": 0.5},
                 {"RAFi": 1.5, "MEKi": 0.5, "KRASi": 0.5, "MCLi": 0.5},
                 {"RAFi": 2.0, "MEKi": 1.0, "METi": 0.5}],
    "KRAS": [{"KRASi": 1.0, "EGFRi": 1.0},
             {"KRASi": 1.0, "MEKi": 1.0},
             {"MEKi": 0.5, "KRASi": 0.5, "MCLi": 1.0},
             {"MEKi": 0.5, "KRASi": 0.5, "BCLi": 2.0},
             {"KRASi": 0.5, "BCLi": 1.0, "MCLi": 2.0}],
    "EGFR": [{"EGFRi": 1.0},
             {"EGFRi": 1.0, "METi": 1.0},
             {"EGFRi": 0.25, "BCLi": 0.5},
             {"EGFRi": 0.5, "KRASi": 0.5, "MCLi": 0.5},
             {"EGFRi": 1.0, "KRASi": 0.5, "BCLi": 0.5}],
}


def constant(doses: dict, days: int = M.TREAT_DAYS) -> dict:
    return {"phases": [{"days": int(days), "doses": {k: float(v) for k, v in doses.items() if v > 0}}]}


def grid() -> np.ndarray:
    """3112 constant-dose regimens [n, ND]: every set of 1 to 4 drugs at a few dose levels."""
    regs = []
    for k in range(1, 5):
        levels = [0.25, 0.5, 1.0, 2.0] if k < 3 else ([0.5, 1.0, 2.0] if k == 3 else [0.5, 1.5])
        for sub in itertools.combinations(range(M.ND), k):
            for lv in itertools.product(levels, repeat=k):
                d = np.zeros(M.ND)
                d[list(sub)] = lv
                regs.append(d)
    return np.array(regs)


_GRID = None


def grid_rewards(patient: dict) -> np.ndarray:
    """Reward of every grid regimen (given for all 84 days) on the fast model."""
    global _GRID
    if _GRID is None:
        _GRID = grid()
    plans = np.repeat(_GRID[:, None, :], M.TREAT_DAYS, axis=1)
    return M.simulate_tissue_mean(patient["K"], patient["healthy"], patient["tumour"], patient["frac"], plans)["reward"]


def doses_of(g: np.ndarray) -> dict:
    return {M.DRUGS[i]: float(g[i]) for i in range(M.ND) if g[i] > 0}


def _lattice(patient: dict, regimen: dict, salt: str, reps: int) -> dict:
    return M.simulate_tissue(patient["K"], patient["healthy"], patient["tumour"], patient["frac"], regimen,
                             seed=P.seed_from(salt, patient["seed"]), reps=reps)


def visible(patient: dict) -> dict:
    """The patient as sequencing of the biopsy shows it: only the clones at or above the detection limit."""
    frac = np.asarray(patient["frac"], dtype=float)
    keep = [0] + [k for k in range(1, len(frac)) if frac[k] >= SEQ_LOD]
    f = frac[keep]
    return dict(patient, tumour=[patient["tumour"][k] for k in keep], frac=f / f.sum())


def hypothetical_clone(major: dict, kind: str, s: float) -> dict:
    """The dominant clone with one more resistance alteration, at strength s (0 weakest, 1 strongest generated)."""
    c = {k: (v.copy() if isinstance(v, np.ndarray) else v) for k, v in major.items() if k not in ("name", "alterations")}
    lerp = lambda lo, hi: lo + s * (hi - lo)                            # noqa: E731
    if kind == "ras":
        c["ras"] = lerp(0.7, 1.0)
    elif kind == "mek1":
        c["kshift"][M.iMEK] = lerp(15.0, 40.0)
    elif kind == "t790m":
        c["kshift"][M.iEGFR] = lerp(3.0, 10.0)
    elif kind == "pten":
        c["pten"] = lerp(0.35, 0.20)
        c["sA"] = min(c["sA"] + 0.2, 0.9)
    elif kind == "alt":
        c["alt"] = lerp(0.6, 1.5)
    elif kind == "mcl1":
        c["bcl"] = c["bcl"] * lerp(1.6, 2.4)
        c["bfrac"], c["mfrac"] = 0.05, lerp(0.80, 0.90)
    elif kind == "bclx":
        c["bcl"] = c["bcl"] * lerp(1.6, 2.4)
        c["bfrac"], c["mfrac"] = lerp(0.80, 0.90), 0.05
    else:
        raise ValueError(kind)
    return c


def hedge_worlds(patient: dict, s: float, share: float) -> list[dict]:
    """One tumour per resistant clone that could be hiding: the visible clones plus that clone at `share`."""
    vis = visible(patient)
    known = {a for c in vis["tumour"] for a in c["alterations"]}
    worlds = []
    for kind in P.MINORS[patient["archetype"]]:
        if kind in known:
            continue
        f = np.r_[vis["frac"] * (1.0 - share), share]
        worlds.append(dict(K=patient["K"], healthy=patient["healthy"], frac=f,
                           tumour=[*vis["tumour"], hypothetical_clone(patient["tumour"][0], kind, s)]))
    return worlds


def hedges(patient: dict, lattice_reps: int = 6, seed: int = 3) -> list[dict]:
    """Rule 6. For each of HEDGES: the regimen that is best on average over hedge_worlds, its reward there,
    and its reward on the real tumour (fast model and lattice)."""
    from . import oracle as RS
    out = []
    for name, s, share in HEDGES:
        worlds = hedge_worlds(patient, s, share)
        x, _, _ = RS.search(worlds, iters=14, pop=256, seed=seed)
        x, r = RS.polish(worlds, x)
        reg = RS.regimen_from_x(x)
        fast = float(M.simulate_tissue_mean(patient["K"], patient["healthy"], patient["tumour"], patient["frac"],
                                            RS.plans_from_x(x[None]))["reward"][0])
        lat = float(_lattice(patient, reg, "cohort-h", lattice_reps)["reward"].mean())
        out.append({"hedge": name, "regimen": reg, "expected": round(float(r), 3), "reward": round(fast, 3),
                    "lattice": round(lat, 3)})
    return out


def screen(patient: dict, lattice_reps: int = 6, keep_rewards: bool = False, hedge: bool = True) -> dict:
    """Cohort test for one patient. Returns dict(eligible, why, guideline, best, best_reward, ...).
    hedge=False stops after rule 5 (out["rules_1_to_5"] says whether those hold)."""
    global _GRID
    out = {"seed": patient["seed"], "archetype": patient["archetype"], "eligible": False,
           "minors": [c["name"][6:] for c in patient["tumour"][1:]],
           "hidden": [c["name"][6:] for c, f in zip(patient["tumour"][1:], patient["frac"][1:]) if f < SEQ_LOD]}
    args = (patient["K"], patient["healthy"], patient["tumour"], patient["frac"])
    plans = np.stack([M.expand_regimen(constant(g)) for g in GUIDELINES[patient["archetype"]]])
    gl = M.simulate_tissue_mean(*args, plans)["reward"]
    out["guideline"] = [round(float(v), 3) for v in gl]
    if gl.max() >= REFRACTORY_BELOW:
        out["why"] = "a guideline regimen works"
        if not keep_rewards:
            return out
    R = grid_rewards(patient)
    G = _GRID
    if keep_rewards:
        out["grid_rewards"] = R.astype(np.float16)
    out["grid_best"] = round(float(R.max()), 3)
    out["window"] = int((R >= 0.8 * max(R.max(), 1e-9)).sum())
    if "why" in out:
        return out
    if R.max() < CURABLE_FROM:
        out["why"] = "no grid regimen cures"
        return out
    un = M.simulate_tissue_mean(*args, np.zeros((1, M.TREAT_DAYS, M.ND)))
    out["untreated_healed"] = round(float(un["healed"][0]), 3)
    if un["healed"][0] > UNTREATED_BELOW:
        out["why"] = "the untreated tumour does not take the lesion"
        return out
    label = np.flatnonzero(((G > 0).sum(axis=1) <= 2) & np.all((G == 0) | (G == 1.0), axis=1))
    jl = int(label[np.argmax(R[label])])
    out["label_best"] = [doses_of(G[jl]), round(float(R[jl]), 3)]
    if R[jl] >= REFRACTORY_BELOW:
        out["why"] = "a label-dose drug or pair works"
        return out
    if not out["hidden"]:
        out["why"] = "sequencing of the biopsy shows every clone"
        return out
    Rv = grid_rewards(visible(patient))
    top = np.argsort(-Rv)[:VISIBLE_TOP]
    out["visible_best"] = [[doses_of(G[j]), round(float(Rv[j]), 3), round(float(R[j]), 3)] for j in top]
    if R[top].max() >= REFRACTORY_BELOW:
        out["why"] = "the answer for the visible clones also cures the hidden one"
        return out

    # the same on the lattice; keep the best of the top few cures
    best, best_r = None, -1.0
    for j in np.argsort(-R)[:3]:
        reg = constant(doses_of(G[j]))
        r = float(_lattice(patient, reg, "cohort", lattice_reps)["reward"].mean())
        if r > best_r:
            best, best_r = reg, r
    out["best"], out["best_reward"] = best, round(best_r, 3)
    gl_lat = [float(_lattice(patient, constant(g), "cohort-g", 8)["reward"].mean()) for g in GUIDELINES[patient["archetype"]]]
    out["guideline_lattice"] = [round(v, 3) for v in gl_lat]
    fails = [float(_lattice(patient, constant(doses_of(G[j])), "cohort-f", lattice_reps)["reward"].mean()) for j in [jl, *top]]
    none = _lattice(patient, {"phases": [{"days": 1, "doses": {}}]}, "cohort-u", lattice_reps)
    out["fail_lattice"] = [round(v, 3) for v in fails]
    out["untreated_lattice"] = round(float(none["healed"].mean()), 3)
    if best_r < CURABLE_FROM:
        out["why"] = "the grid's best regimen fails on the lattice"
    elif max(gl_lat) >= REFRACTORY_BELOW:
        out["why"] = "a guideline regimen works on the lattice"
    elif max(fails) >= REFRACTORY_BELOW:
        out["why"] = "a label-dose pair or the answer for the visible clones works on the lattice"
    elif none["healed"].mean() > 2.5 * UNTREATED_BELOW:
        out["why"] = "the untreated tumour does not take the lesion on the lattice"
    out["rules_1_to_5"] = "why" not in out
    if "why" in out or not hedge:
        return out
    out["hedges"] = hedges(patient, lattice_reps)
    if max(max(h["reward"], h["lattice"]) for h in out["hedges"]) >= REFRACTORY_BELOW:
        out["why"] = "a regimen that hedges against every clone that could be hiding works"
    else:
        out["eligible"] = True
    return out


def find(salt: str, n: int, start: int = 0, max_tries: int = 100000):
    """The first n cohort patients of a salt: yields (index, seed, screen result)."""
    found = 0
    for i in range(start, start + max_tries):
        seed = P.seed_from(salt, i)
        res = screen(P.make_patient(seed))
        if res["eligible"]:
            found += 1
            yield i, seed, res
            if found >= n:
                return
