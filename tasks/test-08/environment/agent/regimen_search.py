"""Regimen search on the public model for a given set of parameters.

Used three ways:
  * hidden side: with the patient's true parameters it gives the ceiling
    ("the best regimen we know for this patient");
  * reference scientist: with parameters it fitted from its own lab data;
  * any agent: the same code is in /app (cellmodel.py + regimen_search.py).

The search is a cross-entropy method over two-phase regimens on the fast
well-mixed model, then the best few are re-scored on the lattice.
"""
from __future__ import annotations

import numpy as np

try:                                    # inside the package
    from . import model as M
except ImportError:                     # as /app/regimen_search.py next to /app/cellmodel.py
    import cellmodel as M

DIM = 2 * M.ND + 2
MIN_DOSE = 0.1                          # doses below this are rounded down to no drug


def plans_from_x(x: np.ndarray) -> np.ndarray:
    """x [Q, 2*ND + 2] -> daily dose plans [Q, TREAT_DAYS, ND].

    x[:, :ND] phase-1 doses, x[:, ND:2ND] phase-2 doses (multiples of the label
    dose, 0..DOSE_MAX), x[:, -2] total length / TREAT_DAYS, x[:, -1] share of
    that length spent in phase 1.
    """
    x = np.atleast_2d(np.asarray(x, dtype=float))
    d1 = np.where(x[:, :M.ND] < MIN_DOSE, 0.0, np.minimum(x[:, :M.ND], M.DOSE_MAX))
    d2 = np.where(x[:, M.ND:2 * M.ND] < MIN_DOSE, 0.0, np.minimum(x[:, M.ND:2 * M.ND], M.DOSE_MAX))
    total = np.clip(np.round(x[:, -2] * M.TREAT_DAYS), 1, M.TREAT_DAYS).astype(int)
    t1 = np.clip(np.round(x[:, -1] * total), 0, total).astype(int)
    days = np.arange(M.TREAT_DAYS)[None, :]
    in1 = days < t1[:, None]
    in2 = (days >= t1[:, None]) & (days < total[:, None])
    return in1[:, :, None] * d1[:, None, :] + in2[:, :, None] * d2[:, None, :]


def regimen_from_x(x: np.ndarray) -> dict:
    x = np.asarray(x, dtype=float)
    d1 = np.where(x[:M.ND] < MIN_DOSE, 0.0, np.minimum(x[:M.ND], M.DOSE_MAX))
    d2 = np.where(x[M.ND:2 * M.ND] < MIN_DOSE, 0.0, np.minimum(x[M.ND:2 * M.ND], M.DOSE_MAX))
    total = int(np.clip(np.round(x[-2] * M.TREAT_DAYS), 1, M.TREAT_DAYS))
    t1 = int(np.clip(np.round(x[-1] * total), 0, total))
    phases = []
    for n, d in ((t1, d1), (total - t1, d2)):
        if n > 0:
            phases.append({"days": n, "doses": {M.DRUGS[i]: round(float(d[i]), 2) for i in range(M.ND) if d[i] > 0}})
    if len(phases) == 2 and phases[0]["doses"] == phases[1]["doses"]:
        phases = [{"days": phases[0]["days"] + phases[1]["days"], "doses": phases[0]["doses"]}]
    return {"phases": phases or [{"days": 1, "doses": {}}]}


def x_from_regimen(regimen: dict) -> np.ndarray:
    """Inverse of regimen_from_x for regimens of one or two phases (others: first two phases)."""
    x = np.zeros(DIM)
    ph = [p for p in regimen.get("phases", []) if int(p.get("days", 0)) > 0][:2]
    if not ph:
        x[-2], x[-1] = 1.0 / M.TREAT_DAYS, 1.0
        return x
    for j, p in enumerate(ph):
        for name, v in (p.get("doses") or {}).items():
            x[j * M.ND + M.IDX[name]] = float(v)
    n1 = int(ph[0]["days"])
    total = min(n1 + (int(ph[1]["days"]) if len(ph) > 1 else 0), M.TREAT_DAYS)
    x[-2] = total / M.TREAT_DAYS
    x[-1] = min(n1, total) / total
    if len(ph) == 1:
        x[M.ND:2 * M.ND] = x[:M.ND]
    return x


def score_plans(param_sets: list[dict], weights, plans: np.ndarray, shaped: bool = True) -> tuple[np.ndarray, np.ndarray]:
    """(objective, reward), each [Q]: weighted mean over the parameter sets.

    The objective adds a small bonus for pushing the hardest clone towards
    extinction with little harm, so that the search has a slope to follow
    before it has found a cure.
    """
    obj = np.zeros(plans.shape[0])
    rew = np.zeros(plans.shape[0])
    for wi, ps in zip(weights, param_sets):
        out = M.simulate_tissue_mean(ps["K"], ps["healthy"], ps["tumour"], ps["frac"], plans)
        rew += wi * out["reward"]
        obj += wi * (out["reward"] + (0.05 * out["kill"] * (1.0 - out["harm"]) if shaped else 0.0))
    return obj, rew


def _one_phase(x: np.ndarray) -> np.ndarray:
    x[:, M.ND:2 * M.ND] = x[:, :M.ND]
    x[:, -1] = 1.0
    return x


def search(param_sets: list[dict], weights=None, iters: int = 16, pop: int = 320, elite: int = 32,
           seed: int = 0, extra_x: np.ndarray | None = None, single_phase: bool = False):
    """Cross-entropy search for the regimen with the best mean reward over `param_sets`.

    Each param set is dict(K, healthy, tumour, frac). Returns (best x, its mean
    reward, the final elite xs sorted best first). single_phase=True searches
    constant-dose regimens only (one phase, any length).
    """
    rng = np.random.default_rng(seed)
    w = np.ones(len(param_sets)) / len(param_sets) if weights is None else np.asarray(weights, float) / np.sum(weights)
    q = np.full(2 * M.ND, 0.30)                 # chance that a drug is in a phase
    mu = np.full(2 * M.ND, 1.0)                 # its dose when it is
    sd = np.full(2 * M.ND, 0.6)
    tmu, tsd = np.array([0.8, 0.7]), np.array([0.25, 0.3])
    best_x, best_o, best_r = None, -1.0, 0.0
    keep = np.zeros((0, DIM))
    for it in range(iters):
        inc = rng.random((pop, 2 * M.ND)) < q
        dose = np.clip(mu + sd * rng.normal(size=(pop, 2 * M.ND)), MIN_DOSE, M.DOSE_MAX)
        x = np.zeros((pop, DIM))
        x[:, :2 * M.ND] = np.where(inc, dose, 0.0)
        x[:, -2:] = np.clip(tmu + tsd * rng.normal(size=(pop, 2)), 0.0, 1.0)
        same = rng.random(pop) < 0.4            # many good regimens are one phase
        x[same, M.ND:2 * M.ND] = x[same, :M.ND]
        if it == 0 and extra_x is not None and len(extra_x):
            n = min(len(extra_x), pop)
            x[:n] = extra_x[:n]
        if single_phase:
            x = _one_phase(x)
        x = np.vstack([x, keep])
        o, r = score_plans(param_sets, w, plans_from_x(x))
        order = np.argsort(-o)
        keep = x[order[:elite]]
        if o[order[0]] > best_o:
            best_o, best_r, best_x = float(o[order[0]]), float(r[order[0]]), x[order[0]].copy()
        on = keep[:, :2 * M.ND] >= MIN_DOSE
        q = np.clip(0.5 * q + 0.5 * on.mean(axis=0), 0.04, 0.96)
        cnt = np.maximum(on.sum(axis=0), 1)
        m = (keep[:, :2 * M.ND] * on).sum(axis=0) / cnt
        v = (((keep[:, :2 * M.ND] - m) ** 2) * on).sum(axis=0) / cnt
        has = on.sum(axis=0) >= 3
        mu = np.where(has, 0.5 * mu + 0.5 * m, mu)
        sd = np.where(has, np.maximum(0.5 * sd + 0.5 * np.sqrt(v), 0.05), sd)
        tmu = 0.5 * tmu + 0.5 * keep[:, -2:].mean(axis=0)
        tsd = np.maximum(0.5 * tsd + 0.5 * keep[:, -2:].std(axis=0), 0.03)
    return best_x, best_r, keep


def polish(param_sets: list[dict], x0: np.ndarray, weights=None, rounds: int = 3, seed: int = 0,
           single_phase: bool = False):
    """Coordinate refinement of one regimen on the mean model. Returns (x, mean reward)."""
    rng = np.random.default_rng(seed)
    w = np.ones(len(param_sets)) / len(param_sets) if weights is None else np.asarray(weights, float) / np.sum(weights)
    x = np.asarray(x0, dtype=float).copy()
    o, r = score_plans(param_sets, w, plans_from_x(x[None]))
    best_o, best_r = float(o[0]), float(r[0])
    step = np.r_[np.full(2 * M.ND, 0.25), 0.08, 0.08]
    for _ in range(rounds):
        dims = np.r_[np.arange(M.ND), DIM - 2] if single_phase else np.arange(DIM)
        for j in rng.permutation(dims):
            hi = M.DOSE_MAX if j < 2 * M.ND else 1.0
            cand = np.repeat(x[None], 8, axis=0)
            cand[:7, j] = np.clip(x[j] + step[j] * np.array([-2, -1, -0.5, 0.5, 1, 2, 3]), 0.0, hi)
            cand[7, j] = 0.0
            if single_phase:
                cand = _one_phase(cand)
            o, r = score_plans(param_sets, w, plans_from_x(cand))
            k = int(np.argmax(o))
            if o[k] > best_o + 1e-6:
                best_o, best_r, x = float(o[k]), float(r[k]), cand[k].copy()
        step = step * 0.5
    return x, best_r


def lattice_reward(ps: dict, regimen: dict, seed: int = 0, reps: int = 6) -> float:
    res = M.simulate_tissue(ps["K"], ps["healthy"], ps["tumour"], ps["frac"], regimen, seed=seed, reps=reps)
    return float(res["reward"].mean())


def best_regimen(param_sets: list[dict], lattice_check: list[dict] | dict | None = None, seed: int = 0,
                 top: int = 4, reps: int = 6, restarts: int = 1, extra_x: np.ndarray | None = None,
                 single_phase: bool = False, **kw):
    """Search, polish, and (optionally) pick among the top candidates on the lattice.

    lattice_check: one param set or a list of them to re-score the top `top`
    distinct candidates on with simulate_tissue (None = trust the mean model).
    Returns dict(regimen, x, mean_reward, lattice_reward or None).
    """
    cands = []
    for rs in range(restarts):
        x, r, elite = search(param_sets, seed=seed + 1000 * rs, extra_x=extra_x, single_phase=single_phase, **kw)
        x, r = polish(param_sets, x, seed=seed + rs, single_phase=single_phase)
        cands.append((x, r))
        for e in elite:
            if len(cands) >= top * (rs + 1):
                break
            if all(np.abs(e - c[0]).max() > 0.15 for c in cands):
                cands.append(polish(param_sets, e, rounds=2, seed=seed + rs + 1, single_phase=single_phase))
    cands.sort(key=lambda c: -c[1])
    cands = cands[:top]
    x, r = cands[0]
    out = {"regimen": regimen_from_x(x), "x": x, "mean_reward": r, "lattice_reward": None}
    if lattice_check is not None:
        checks = [lattice_check] if isinstance(lattice_check, dict) else list(lattice_check)
        best = -1.0
        for cx, cr in cands:
            reg = regimen_from_x(cx)
            lr = float(np.mean([lattice_reward(ps, reg, seed=seed + 101 + i, reps=reps) for i, ps in enumerate(checks)]))
            if lr > best:
                best = lr
                out = {"regimen": reg, "x": cx, "mean_reward": cr, "lattice_reward": lr}
    return out
