"""Who is in the benchmark: patients that guideline regimens fail, but who can be cured.

Hidden side. A generated patient joins the cohort when

  1. refractory: each of the five guideline regimens for the tumour type
     scores below REFRACTORY_BELOW. Two of them are the label-dose standards of
     care. The other three are the constant regimens that cure the largest
     share of the patients the earlier ones leave uncured, picked one after
     the other by exhaustive search over the grid in a development population
     of 1600 generated patients;
  2. curable: at least one constant-dose regimen out of a grid of 3112 scores
     CURABLE_FROM or more, on the fast model and again on the tissue lattice.

So a guess from the chart alone fails by construction, and a cure exists for
every patient. Functional precision oncology is used for exactly these
patients: the ones standard treatment does not help.
"""
from __future__ import annotations

import itertools

import numpy as np

from . import model as M
from . import patients as P

REFRACTORY_BELOW = 0.10
CURABLE_FROM = 0.50

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


def screen(patient: dict, lattice_reps: int = 6, keep_rewards: bool = False) -> dict:
    """Cohort test for one patient. Returns dict(eligible, why, guideline, best, best_reward, ...)."""
    out = {"seed": patient["seed"], "archetype": patient["archetype"], "eligible": False,
           "minors": [c["name"][6:] for c in patient["tumour"][1:]]}
    plans = np.stack([M.expand_regimen(constant(g)) for g in GUIDELINES[patient["archetype"]]])
    gl = M.simulate_tissue_mean(patient["K"], patient["healthy"], patient["tumour"], patient["frac"], plans)["reward"]
    out["guideline"] = [round(float(v), 3) for v in gl]
    if gl.max() >= REFRACTORY_BELOW:
        out["why"] = "a guideline regimen works"
        if not keep_rewards:
            return out
    R = grid_rewards(patient)
    if keep_rewards:
        out["grid_rewards"] = R.astype(np.float16)
    out["grid_best"] = round(float(R.max()), 3)
    out["window"] = int((R >= 0.8 * max(R.max(), 1e-9)).sum())
    if "why" in out:
        return out
    if R.max() < CURABLE_FROM:
        out["why"] = "no grid regimen cures"
        return out
    # confirm on the lattice; keep the best of the top few
    best, best_r = None, -1.0
    for j in np.argsort(-R)[:3]:
        reg = constant(doses_of(_GRID[j]))
        r = float(M.simulate_tissue(patient["K"], patient["healthy"], patient["tumour"], patient["frac"], reg,
                                    seed=P.seed_from("cohort", patient["seed"]), reps=lattice_reps)["reward"].mean())
        if r > best_r:
            best, best_r = reg, r
    out["best"], out["best_reward"] = best, round(best_r, 3)
    # the guidelines must fail on the lattice too
    gl_lat = [float(M.simulate_tissue(patient["K"], patient["healthy"], patient["tumour"], patient["frac"], constant(g),
                                      seed=P.seed_from("cohort-g", patient["seed"]), reps=8)["reward"].mean())
              for g in GUIDELINES[patient["archetype"]]]
    out["guideline_lattice"] = [round(v, 3) for v in gl_lat]
    if best_r < CURABLE_FROM:
        out["why"] = "the grid's best regimen fails on the lattice"
    elif max(gl_lat) >= REFRACTORY_BELOW:
        out["why"] = "a guideline regimen works on the lattice"
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
