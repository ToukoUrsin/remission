"""The reference scientist for Remission.

It plays through the lab's public tools only (the `call` function it is given)
and uses nothing but what every agent has in /app: the public cell model
(cellmodel.py), the regimen search (regimen_search.py), numpy and scipy.

What it does, in lab days:

  day 0    deep-sequence the biopsy; pERK/pAKT panels and long viability
           panels on tumour and on healthy cells
  day 10   fit the model to everything; wells where cells grow back although
           the fit says they should be dead hold a resistant clone the fit does
           not know: expand those survivors into lines. Search regimens on an
           ensemble of fits and culture the best ones
  day 13   sequence the lines and test drugs on them
  day 18   refit with the resistant clones named; a candidate whose culture
           grew back gives one more line; search again; culture the new
           candidates
  day 27   final fit; submit the regimen with the best expected reward among
           those whose culture behaved as the fit says
"""
from __future__ import annotations

import time

import numpy as np
from scipy.optimize import least_squares

try:                                    # in the agent container: /app/cellmodel.py, /app/regimen_search.py
    import cellmodel as M
    import regimen_search as RS
except ImportError:                     # in the repository
    from remission import model as M
    from remission import oracle as RS

CSS = M.CSS_LABEL
PATHWAY = ("RAFi", "MEKi", "EGFRi", "METi", "KRASi", "PI3Ki")
PRIMARY = {"BRAF V600E": "RAFi", "KRAS G12C": "KRASi", "EGFR L858R": "EGFRi"}
SOFT = 30.0                             # cells; counts are compared as ln(count + SOFT)
LINE_COST = 8 + 16 + 16                 # wells to expand a line, sequence it and test the drugs on it


def conc(**mult) -> dict:
    """{drug: multiple of its average concentration at the label dose} -> micromolar."""
    return {d: round(float(m * M.PK[d]["css"]), 6) for d, m in mult.items() if m > 0}


def conc_vec(well: dict) -> np.ndarray:
    v = np.zeros(M.ND)
    for d, c in well.items():
        v[M.IDX[d]] = c
    return v


# ---------------------------------------------------------------------------
# Parameters: fit vector <-> clones
# ---------------------------------------------------------------------------
class Structure:
    """Which parameters are free for a tumour with these alterations and these minor clones."""

    def __init__(self, major_labels, minor_labels, share_obs=None):
        self.major_labels = list(major_labels)
        self.minor_labels = list(minor_labels)
        self.names, self.lo, self.hi, self.mu, self.sd, self.is_log = [], [], [], [], [], []
        for d in M.DRUGS:
            mu = np.log(M.K_TYPICAL[d])
            self._add(f"K:{d}", mu - 3.2 * M.K_LN_SD, mu + 3.2 * M.K_LN_SD, mu, M.K_LN_SD, True)
        for d in M.DRUGS:
            mu = np.log(M.TOX_TYPICAL)
            self._add(f"tox:{d}", mu - 3.5 * M.TOX_LN_SD, mu + 3.5 * M.TOX_LN_SD, mu, M.TOX_LN_SD, True)
        for k, (lo, hi) in M.HEALTHY_RANGE.items():
            self._add(f"h:{k}", lo, hi)
        rng = dict(M.TUMOUR_RANGE)
        self.major_fixed = dict(braf=0.0, ras=0.0, pi3k=0.0, pten=1.0)
        self.major_kshift = np.ones(M.ND)
        for lab in self.major_labels:
            for prm, lo, hi, how in M.ALTERATIONS.get(lab, []):
                if how == "kshift":
                    self.major_kshift[M.IDX[prm]] = lo
                elif how == "set" and lo == hi:
                    self.major_fixed[prm] = lo
                elif how == "set":
                    rng[prm] = (lo, hi)
                    self.major_fixed.pop(prm, None)
        for k, (lo, hi) in rng.items():
            if k in ("egfr", "alt") and lo < 0.1:
                self._add(f"t:{k}", np.log(lo), np.log(hi), is_log=True)
            else:
                self._add(f"t:{k}", lo, hi)
        self.major_free = list(rng)
        share_obs = share_obs or {}
        self.minor_free = []
        for i, lab in enumerate(self.minor_labels):
            f = share_obs.get(lab)
            if f:
                self._add(f"m{i}:share", np.log(0.01), np.log(0.12), np.log(f), 0.25, True)
            else:
                self._add(f"m{i}:share", np.log(0.002), np.log(0.03), np.log(0.008), 0.7, True)
            free = []
            for prm, lo, hi, how in M.ALTERATIONS[lab]:
                if lo < hi:
                    self._add(f"m{i}:{prm}", np.log(lo) if how == "kshift" else lo, np.log(hi) if how == "kshift" else hi,
                              is_log=how == "kshift")
                    free.append(prm)
            self.minor_free.append(free)
        self.lo, self.hi = np.array(self.lo), np.array(self.hi)
        self.mu, self.sd = np.array(self.mu), np.array(self.sd)
        self.is_log = np.array(self.is_log)
        self.ix = {n: i for i, n in enumerate(self.names)}
        self.n = len(self.names)
        self.C = 2 + len(self.minor_labels)

    def _add(self, name, lo, hi, mu=None, sd=None, is_log=False):
        self.names.append(name)
        self.lo.append(lo)
        self.hi.append(hi)
        self.mu.append(0.5 * (lo + hi) if mu is None else mu)
        self.sd.append(np.inf if sd is None else sd)
        self.is_log.append(is_log)

    def centre(self) -> np.ndarray:
        return np.clip(self.mu, self.lo, self.hi)

    def draw(self, rng) -> np.ndarray:
        z = np.where(np.isfinite(self.sd), self.mu + np.where(np.isfinite(self.sd), self.sd, 0.0) * rng.normal(size=self.n),
                     rng.uniform(self.lo, self.hi))
        return np.clip(z, self.lo, self.hi)

    def build(self, Z: np.ndarray):
        """Z [B, n] -> K [B, ND], P {key: [B, C]}, kshift [B, C, ND], tox [B, C, ND], tumour shares [B, C-1].

        Clone 0 is the healthy cell, clone 1 the dominant tumour clone, then the minor clones.
        """
        Z = np.atleast_2d(Z)
        B, C = Z.shape[0], self.C

        def v(name):
            x = Z[:, self.ix[name]]
            return np.exp(x) if self.is_log[self.ix[name]] else x

        K = np.exp(Z[:, :M.ND])
        P = {k: np.zeros((B, C)) for k in M.CLONE_KEYS}
        ksh = np.ones((B, C, M.ND))
        tox = np.zeros((B, C, M.ND))
        tox[:, 0, :] = np.exp(Z[:, M.ND:2 * M.ND])
        for k, val in M.HEALTHY_FIXED.items():
            P[k][:, 0] = val
        for k in M.HEALTHY_RANGE:
            if k != "margin":
                P[k][:, 0] = v(f"h:{k}")
        for k, val in self.major_fixed.items():
            P[k][:, 1] = val
        for k in self.major_free:
            if k != "margin":
                P[k][:, 1] = v(f"t:{k}")
        ksh[:, 1, :] = self.major_kshift
        P["bcl"][:] = 1.0
        D0 = death_drive0({k: P[k][:, :2] for k in M.CLONE_KEYS})
        P["bcl"][:, 0] = v("h:margin") * D0[:, 0]
        P["bcl"][:, 1] = v("t:margin") * D0[:, 1]
        shares = []
        for i, lab in enumerate(self.minor_labels):
            c = 2 + i
            for k in M.CLONE_KEYS:
                P[k][:, c] = P[k][:, 1]
            ksh[:, c, :] = ksh[:, 1, :]
            for prm, lo, hi, how in M.ALTERATIONS[lab]:
                val = lo if lo == hi else v(f"m{i}:{prm}")
                if how == "set":
                    P[prm][:, c] = val
                elif how == "mul":
                    P[prm][:, c] = P[prm][:, c] * val
                elif how == "add":
                    P[prm][:, c] = np.minimum(P[prm][:, c] + val, 0.9)
                elif how == "kshift":
                    ksh[:, c, M.IDX[prm]] = val
            shares.append(v(f"m{i}:share"))
        sh = np.stack(shares, axis=1) if shares else np.zeros((B, 0))
        frac = np.concatenate([1.0 - sh.sum(axis=1, keepdims=True), sh], axis=1)
        return K, P, ksh, tox, frac

    def param_set(self, z: np.ndarray) -> dict:
        """One fit vector -> the dict(K, healthy, tumour, frac) the regimen search takes."""
        K, P, ksh, tox, frac = self.build(z[None])
        clones = [M.make_clone(name="healthy" if c == 0 else f"clone{c}", healthy=c == 0,
                               kshift=ksh[0, c], tox=tox[0, c], **{k: float(P[k][0, c]) for k in M.CLONE_KEYS})
                  for c in range(self.C)]
        return {"K": K[0], "healthy": clones[0], "tumour": clones[1:], "frac": frac[0]}


def death_drive0(P: dict) -> np.ndarray:
    E, A = M.baseline(P)
    return (1.0 - P["sE"] * E / (E + M.HE)) * (1.0 - P["sA"] * A / (A + M.HA))


def sim_wells(K, P, ksh, tox, cols, frac, cw, nstep, reads, cap):
    """Deterministic wells for B parameter sets at once.

    cols: which clones are in the wells; frac [B, len(cols)]; cw [W, ND].
    Returns N [R, B, W, C] at the read steps and absolute ERK, AKT activity [B, W, C] at the last step.
    """
    p = {k: P[k][:, None, cols] for k in M.CLONE_KEYS}
    p["kshift"] = ksh[:, None, cols, :]
    p["tox"] = tox[:, None, cols, :]
    c4 = cw[None, :, None, :]
    occ = M.occupancy(c4, K[:, None, None, :], p["kshift"])
    E0, A0 = M.baseline(p)
    N = M.SEED_CELLS * np.broadcast_to(frac[:, None, :], (K.shape[0], cw.shape[0], len(cols))).copy()
    F = np.ones_like(N)
    want = set(int(r) for r in reads)
    rec = {}
    E = A = None
    for s in range(nstep + 1):
        E, A = M.signalling(p, occ, F)
        if s in want:
            rec[s] = N.copy()
        if s == nstep:
            break
        kp, kd = M.fate_rates(p, E, A, occ, c4)
        crowd = np.clip(1.0 - N.sum(axis=2, keepdims=True) / cap, 0.0, 1.0)
        N = N * np.exp((kp * crowd - kd) * M.DT)
        F = M.step_feedback(p, F, E, E0, M.DT)
    return np.stack([rec[int(r)] for r in reads]), E, A, N


# ---------------------------------------------------------------------------
# The scientist
# ---------------------------------------------------------------------------
class Scientist:
    def __init__(self, call, log=print, seed: int = 0, budget: str = "full"):
        self.call, self.log = call, log
        self.rng = np.random.default_rng(seed)
        self.budget = budget
        self.data: list[dict] = []           # finished experiments, in the form the fit uses
        self.lines: dict[str, dict] = {}     # derived lines: parent, sequencing
        self.st: Structure | None = None
        self.z: np.ndarray | None = None
        self.ens: list[np.ndarray] = []
        self.t0 = time.time()
        self.notes: list[str] = []

    # ------------------------------------------------------------------ lab helpers
    def say(self, msg: str):
        self.notes.append(msg)
        self.log(f"[{time.time() - self.t0:6.1f}s] {msg}")

    def lab(self, tool: str, **args) -> dict:
        r = self.call(tool, **args)
        if isinstance(r, dict) and "error" in r:
            self.say(f"{tool} refused: {r['error']}")
        return r

    def start(self, kind: str, cells: str, wells=None, days=None, read_days=None, tag: str = "") -> dict | None:
        args = {"kind": kind, "cells": cells}
        if wells is not None:
            args["wells"] = wells
        if days is not None:
            args["days"] = days
        if read_days is not None:
            args["read_days"] = read_days
        r = self.lab("start_experiment", **args)
        if "error" in r:
            return None
        return {"id": r["id"], "kind": kind, "cells": cells, "wells": wells, "days": days, "tag": tag}

    def collect(self, exp: dict | None) -> dict | None:
        """Read a finished experiment and file it for the fit."""
        if exp is None:
            return None
        r = self.lab("get_results", experiment_id=exp["id"])
        if "error" in r or r.get("status") != "done":
            return None
        exp = dict(exp)
        exp["wells"] = r.get("wells", exp.get("wells"))
        if exp["kind"] == "sequencing":
            exp["alterations"] = {a["alteration"]: a["cell_fraction"] for a in r["alterations"]}
            return exp
        exp["conc"] = np.array([conc_vec(w) for w in r["wells"]])
        if exp["kind"] == "viability":
            exp["reads"] = [int(round(d / M.DT)) for d in r["read_days"]]
            exp["nstep"] = int(round(exp["days"] / M.DT))
            exp["counts"] = np.array(r["counts"], dtype=float)           # [W, R]
        else:
            exp["nstep"] = int(round(exp["days"] / M.DT))
            exp["pERK"], exp["pAKT"] = np.array(r["pERK"]), np.array(r["pAKT"])
        self.data.append(exp)
        return exp

    # ------------------------------------------------------------------ model of the data
    def line_frac(self, name: str, frac: np.ndarray) -> np.ndarray | None:
        """Clone mix of a cell source for B parameter sets [B, C-1], or None if unknown."""
        if name == "tumour":
            return frac
        line = self.lines.get(name)
        if line is None or line.get("alterations") is None:
            return None
        f = np.zeros(self.st.C - 1)
        for i, lab in enumerate(self.st.minor_labels):
            f[1 + i] = min(line["alterations"].get(lab, 0.0), 1.0)
        f[1:] *= min(1.0, 0.999 / max(f[1:].sum(), 1e-9))
        f[0] = 1.0 - f[1:].sum()
        return np.broadcast_to(f, frac.shape)

    def predict(self, Z: np.ndarray, only: dict | None = None) -> list:
        """Predictions for every filed experiment: viability -> N [R, B, W]; phospho -> (pERK, pAKT) [B, W]."""
        K, P, ksh, tox, frac = self.st.build(Z)
        tum = list(range(1, self.st.C))
        out = []
        for e in self.data if only is None else [only]:
            if e["cells"] == "healthy":
                cols, fr, cap = [0], np.ones((K.shape[0], 1)), M.CAP_HEALTHY
            else:
                cols, cap = tum, M.CAP_TUMOUR
                fr = self.line_frac(e["cells"], frac)
                if fr is None:
                    out.append(None)
                    continue
            if e["kind"] == "viability":
                N, _, _, _ = sim_wells(K, P, ksh, tox, cols, fr, e["conc"], e["nstep"], e["reads"], M.SEED_CELLS * cap)
                out.append(N.sum(axis=3))
            else:
                cw = np.vstack([e["conc"], np.zeros((1, M.ND))])
                _, E, A, N = sim_wells(K, P, ksh, tox, cols, fr, cw, e["nstep"], [e["nstep"]], M.SEED_CELLS * cap)
                w = N / np.maximum(N.sum(axis=2, keepdims=True), 1e-9)
                pe, pa = (w * E).sum(axis=2), (w * A).sum(axis=2)
                out.append((pe[:, :-1] / np.maximum(pe[:, -1:], 1e-9), pa[:, :-1] / np.maximum(pa[:, -1:], 1e-9)))
        return out

    def residuals(self, Z: np.ndarray) -> np.ndarray:
        Z = np.atleast_2d(Z)
        parts = []
        for e, pred in zip(self.data, self.predict(Z)):
            if pred is None:
                continue
            if e["kind"] == "viability":
                obs = e["counts"].T[:, None, :]                                    # [R, 1, W]
                sig = np.sqrt(0.10 ** 2 + 0.06 ** 2 + 1.0 / (obs + SOFT))
                r = (np.log(pred + SOFT) - np.log(obs + SOFT)) / sig
                parts.append(np.moveaxis(r, 1, 0).reshape(Z.shape[0], -1))
            else:
                for obs, pr in ((e["pERK"], pred[0]), (e["pAKT"], pred[1])):
                    sig = np.sqrt((0.10 * obs) ** 2 + 0.03 ** 2)
                    parts.append((pr - obs[None, :]) / sig[None, :])
        has = np.isfinite(self.st.sd)
        parts.append((Z[:, has] - self.st.mu[has]) / self.st.sd[has])
        return np.concatenate(parts, axis=1)

    def jac(self, z: np.ndarray) -> np.ndarray:
        st = self.st
        h = np.where(st.is_log, 0.02, 0.01 * (st.hi - st.lo))
        h = np.where(z + h > st.hi, -h, h)
        Z = np.vstack([z[None], z[None] + np.diag(h)])
        R = self.residuals(Z)
        return ((R[1:] - R[0]) / h[:, None]).T

    def fit(self, starts: int = 3, warm: dict | None = None, max_nfev: int = 40):
        """Least-squares fit of every free parameter to all filed data. Keeps the best of a few starts."""
        st = self.st
        x0s = [st.centre()]
        if warm:
            x = st.centre()
            for name, val in warm.items():
                if name in st.ix:
                    x[st.ix[name]] = val
            x0s.insert(0, np.clip(x, st.lo, st.hi))
        while len(x0s) < starts:
            x0s.append(0.5 * (st.centre() + st.draw(self.rng)))
        best = None
        for x0 in x0s[:starts]:
            x0 = np.clip(x0, st.lo + 1e-6, st.hi - 1e-6)
            try:
                res = least_squares(lambda x: self.residuals(x)[0], x0, jac=self.jac, bounds=(st.lo, st.hi),
                                    loss="soft_l1", f_scale=2.0, max_nfev=max_nfev, x_scale="jac")
            except Exception as exc:  # noqa: BLE001
                self.say(f"fit start failed: {exc!r}")
                continue
            if best is None or res.cost < best.cost:
                best = res
        self.z = best.x
        J = self.jac(self.z)
        H = J.T @ J + 1e-3 * np.eye(st.n)
        try:
            self.cov = np.linalg.inv(H)
            self.chol = np.linalg.cholesky(self.cov + 1e-12 * np.eye(st.n))
        except np.linalg.LinAlgError:
            self.chol = np.diag(0.05 * (st.hi - st.lo))
        n_ens = 4 if self.budget == "fast" else 6
        self.ens = [self.z] + [np.clip(self.z + 1.3 * self.chol @ self.rng.normal(size=st.n), st.lo, st.hi)
                               for _ in range(n_ens - 1)]
        r = self.residuals(self.z)[0]
        self.say(f"fit: {st.n} parameters, {len(r)} residuals, rms {np.sqrt(np.mean(r ** 2)):.2f}")
        return best

    def warm(self) -> dict:
        return {} if self.z is None else {n: self.z[i] for n, i in self.st.ix.items()}

    # ------------------------------------------------------------------ regimens
    def sets(self) -> list[dict]:
        return [self.st.param_set(z) for z in self.ens]

    def expected(self, xs: np.ndarray) -> np.ndarray:
        sets = self.sets()
        _, rew = RS.score_plans(sets, np.ones(len(sets)) / len(sets), RS.plans_from_x(xs), shaped=False)
        return rew

    def candidates(self, n: int, extra: list | None = None) -> list[np.ndarray]:
        """Constant-dose regimens with the best mean reward over the ensemble of fits, best first."""
        sets = self.sets()
        fast = self.budget == "fast"
        ex = np.array(extra) if extra else None
        x, r, elite = RS.search(sets, iters=9 if fast else 14, pop=160 if fast else 256, elite=24, seed=int(self.rng.integers(1 << 30)),
                                extra_x=ex, single_phase=True)
        x, r = RS.polish(sets, x, rounds=2 if fast else 3, single_phase=True)
        out = [x]
        for e in elite:
            if len(out) >= n:
                break
            if all(np.abs(e[:M.ND] - c[:M.ND]).max() > 0.25 for c in out):
                out.append(e)
        rew = self.expected(np.array(out))
        order = np.argsort(-rew)
        self.say("candidates: " + "; ".join(f"{self.show(out[i])} -> {rew[i]:.2f}" for i in order[:4]))
        return [out[i] for i in order]

    @staticmethod
    def show(x: np.ndarray) -> str:
        reg = RS.regimen_from_x(x)["phases"][0]
        return "+".join(f"{d[:-1]} {v:g}" for d, v in reg["doses"].items()) + f" x{reg['days']}d"

    @staticmethod
    def well_of(x: np.ndarray, scale: float = 1.0) -> dict:
        d = RS.regimen_from_x(x)["phases"][0]["doses"]
        return conc(**{k: v * scale for k, v in d.items()})

    # ------------------------------------------------------------------ designs
    def tumour_panel(self, primary: str) -> list[dict]:
        w = [{}, {}, {}]
        for d in M.DRUGS:
            w += [conc(**{d: m}) for m in (0.3, 1.0, 3.0)]
        for i, a in enumerate(M.DRUGS):
            for b in M.DRUGS[i + 1:]:
                w.append(conc(**{a: 1.0, b: 1.0}))
        others = [d for d in M.DRUGS if d != primary]
        w += [conc(**{primary: 1.0, d: 3.0}) for d in others]
        second = "MEKi" if primary != "MEKi" else "RAFi"
        w += [conc(**{primary: 1.0, second: 1.0, d: 1.0}) for d in others if d != second]
        w += [conc(**{primary: 3.0, second: 3.0}), conc(**{primary: 0.3, second: 0.3}), conc(BCLi=3.0, MCLi=3.0),
              conc(**{primary: 1.0, "BCLi": 1.0, "MCLi": 1.0}), conc(**{primary: 3.0, "BCLi": 3.0}),
              conc(**{primary: 3.0, "MCLi": 3.0}), conc(MEKi=3.0, BCLi=3.0), conc(MEKi=3.0, MCLi=3.0),
              conc(PI3Ki=3.0, BCLi=3.0), conc(PI3Ki=3.0, MCLi=3.0)]
        rng = np.random.default_rng(11)
        while len(w) < 96:
            k = int(rng.integers(3, 6))
            ds = rng.choice(M.ND, size=k, replace=False)
            w.append(conc(**{M.DRUGS[j]: float(rng.choice([0.3, 1.0, 3.0])) for j in ds}))
        return w[:96]

    def healthy_panel(self, primary: str) -> list[dict]:
        w = [{}, {}, {}]
        for d in M.DRUGS:
            w += [conc(**{d: m}) for m in (1.0, 4.0, 12.0)]
        for i, a in enumerate(M.DRUGS):
            for b in M.DRUGS[i + 1:]:
                w.append(conc(**{a: 2.0, b: 2.0}))
        w += [conc(BCLi=4.0, MCLi=4.0), conc(MEKi=2.0, BCLi=2.0, MCLi=2.0), conc(**{primary: 2.0, "MEKi": 2.0, "BCLi": 2.0}),
              conc(**{primary: 2.0, "MEKi": 2.0, "MCLi": 2.0}), conc(EGFRi=2.0, MEKi=2.0, PI3Ki=2.0),
              conc(PI3Ki=2.0, BCLi=2.0, MCLi=2.0), conc(**{d: 1.0 for d in PATHWAY}), conc(**{d: 1.0 for d in M.DRUGS})]
        rng = np.random.default_rng(12)
        while len(w) < 72:
            k = int(rng.integers(3, 6))
            ds = rng.choice(M.ND, size=k, replace=False)
            w.append(conc(**{M.DRUGS[j]: float(rng.choice([0.5, 1.5, 4.0])) for j in ds}))
        return w[:72]

    # ------------------------------------------------------------------ resistant clones
    def outgrowth(self, exp: dict) -> list[tuple[float, int]]:
        """Wells of a tumour culture that hold far more cells at the end than the fit explains."""
        pred = self.predict(self.z[None], only=exp)[0][-1, 0]            # [W] at the last read
        obs = exp["counts"][:, -1]
        prev = exp["counts"][:, -2] if exp["counts"].shape[1] > 1 else obs
        excess = np.log((obs + SOFT) / (pred + SOFT))
        found = [(float(excess[i]), int(i)) for i in range(len(obs))
                 if excess[i] > 1.0 and obs[i] >= 80 and obs[i] > 1.15 * prev[i]]
        return sorted(found, reverse=True)

    def pick_wells(self, exp: dict, found: list[tuple[float, int]], k: int = 2) -> list[int]:
        """Up to k outgrowth wells whose drug sets differ, so that they can hold different clones."""
        picks: list[int] = []
        for _, i in found:
            drugs = set(exp["wells"][i])
            if all(len(drugs ^ set(exp["wells"][j])) >= 2 for j in picks):
                picks.append(i)
            if len(picks) >= k:
                break
        return picks

    # ------------------------------------------------------------------ the session
    def run(self) -> dict:
        chart = self.lab("patient_chart")
        major = list(chart["bulk_sequencing"]["alterations"])
        primary = next((PRIMARY[a] for a in major if a in PRIMARY), "MEKi")
        self.major, self.minors, self.share_obs = major, [], {}
        self.say(f"chart: {chart['tissue']}; {', '.join(major)}; driver drug {primary}")

        # ---- day 0: everything that needs time starts now
        seq = self.start("sequencing", "tumour")
        ph = [conc(**{d: m}) for d in PATHWAY for m in (0.1, 0.5, 2.5)]
        ph += [conc(RAFi=0.5, MEKi=0.5), conc(EGFRi=2.5, METi=2.5), conc(MEKi=0.5, PI3Ki=0.5),
               conc(**{primary: 2.5, "EGFRi" if primary != "EGFRi" else "METi": 2.5}), conc(KRASi=2.5, RAFi=2.5),
               conc(EGFRi=2.5, METi=2.5, KRASi=2.5)]
        ph2 = [conc(MEKi=0.5), conc(MEKi=2.5), conc(**{primary: 2.5}), conc(EGFRi=2.5) if primary != "EGFRi" else conc(RAFi=2.5),
               conc(PI3Ki=2.5), conc(MEKi=2.5, EGFRi=2.5), conc(**{primary: 1.0, "MEKi": 1.0}) if primary != "MEKi" else conc(RAFi=1.0, MEKi=1.0),
               conc(METi=2.5)]
        hp = [conc(**{d: m}) for d in PATHWAY for m in (0.3, 3.0)]
        hp2 = [conc(MEKi=1.0), conc(EGFRi=3.0), conc(RAFi=3.0), conc(PI3Ki=3.0)]
        first = [self.start("phospho", "tumour", ph, 0.25), self.start("phospho", "tumour", ph2, 2.0),
                 self.start("phospho", "healthy", hp, 0.25), self.start("phospho", "healthy", hp2, 2.0),
                 self.start("viability", "healthy", self.healthy_panel(primary), 7, [1, 2, 3, 5, 7])]
        tpanel = self.start("viability", "tumour", self.tumour_panel(primary), 10, [1, 2, 3, 4, 6, 8, 10], tag="panel")
        self.lab("advance_time", days=10)

        # ---- day 10: first fit; survivors the fit cannot explain become lines; culture the first candidates
        seq = self.collect(seq)
        for e in first:
            self.collect(e)
        tpanel = self.collect(tpanel)
        self.share_obs = {a: f for a, f in (seq or {}).get("alterations", {}).items() if a not in major and a in M.ALTERATIONS}
        self.minors = list(self.share_obs)
        self.say(f"deep sequencing of the biopsy: minor clones {self.share_obs or 'none above 2 %'}")
        self.refit(starts=2)
        self.expand(tpanel, 2)
        pool: list[dict] = []
        verA = self.verify(self.candidates(6), days=8, reads=[2, 4, 6, 8], low=True)
        self.lab("advance_time", days=3)

        # ---- day 13: name the resistant clones (sequencing) and test drugs on them
        pending = self.study_lines()
        self.lab("advance_time", days=5)

        # ---- day 18: refit with the clones named; a candidate whose culture grew back points to one more clone
        self.collect_lines(pending)
        pool += self.collect_verify(verA)
        self.refit(starts=2)
        if verA["tumour"] is not None:
            self.expand(next(e for e in self.data if e["id"] == verA["tumour"]["id"]), 1)
        cand = self.candidates(6, extra=[c["x"] for c in pool])
        cand = [x for x in cand if all(np.abs(x[:M.ND] - c["x"][:M.ND]).max() > 0.12 for c in pool)][:5]
        verB = self.verify(cand, days=9, reads=[2, 4, 6, 8, 9], low=False)
        self.lab("advance_time", days=3)

        # ---- day 21: study a late line, if any
        pending = self.study_lines()
        self.lab("advance_time", days=6)

        # ---- day 27: final fit and choice
        self.collect_lines(pending)
        pool += self.collect_verify(verB)
        self.refit(starts=1, max_nfev=30)
        tum = {e["id"]: e for e in self.data}
        for c in pool:
            exp = tum.get(c["exp"])
            bad = {i for _, i in self.outgrowth(exp)} if exp is not None else set()
            c["explained"] = exp is not None and not (bad & set(c["wells"]))
        rew = self.expected(np.array([c["x"] for c in pool])) if pool else np.zeros(0)
        for c, r in zip(pool, rew):
            c["expected"] = float(r)
        self.say("cultured regimens: " + "; ".join(f"{self.show(c['x'])} explained={c['explained']} exp={c['expected']:.2f}" for c in pool))
        ok = [c for c in pool if c["explained"]] or pool
        best = max(ok, key=lambda c: c["expected"]) if ok else None
        # the last fit may have changed the picture (a clone named on day 26): search once more
        if best is None or best["expected"] < 0.35:
            x = self.candidates(1, extra=[c["x"] for c in pool])[0]
            r = float(self.expected(x[None])[0])
            if best is None or r > best["expected"] + 0.15:
                self.say("no cultured regimen is expected to cure under the final fit; taking the final search's best")
                best = {"x": x, "expected": r, "explained": None}
        regimen = RS.regimen_from_x(best["x"])
        self.say(f"submitting {regimen} (expected reward {best['expected']:.2f})")
        r = self.lab("submit_regimen", regimen=regimen)
        return {"regimen": regimen, "expected": best["expected"], "explained": best["explained"],
                "accepted": bool(r.get("accepted", False)), "minors": self.minors, "notes": self.notes}

    # ------------------------------------------------------------------ steps
    def refit(self, starts: int = 2, max_nfev: int = 40):
        warm = self.warm()
        self.st = Structure(self.major, self.minors, self.share_obs)
        self.fit(starts=starts, warm=warm, max_nfev=max_nfev)

    def wells_left(self) -> int:
        return int(self.lab("lab_status").get("wells_left", 0))

    def expand(self, exp: dict | None, k: int):
        """Expand the survivors of up to k wells that grew back although the fit says they should not."""
        if exp is None:
            return
        found = self.outgrowth(exp)
        picks = self.pick_wells(exp, found, k) if found else []
        self.say(f"{exp['id']}: unexplained regrowth in {len(found)} wells; expanding wells {picks}")
        for i in picks:
            if self.wells_left() < LINE_COST + 12:
                self.say("not enough wells left to study another line")
                break
            r = self.lab("derive_line", experiment_id=exp["id"], well=i)
            if r.get("ok"):
                self.lines[r["line"]] = {"well": i, "alterations": None, "sel": exp["wells"][i], "ready": r["ready_day"]}

    def study_lines(self) -> dict:
        """Sequence every new line and test the drugs on it. Returns the running experiments by line."""
        day = float(self.lab("lab_status").get("day", 0))
        out = {}
        for name, line in self.lines.items():
            if line.get("studied") or line["ready"] > day + 1e-9:
                continue
            line["studied"] = True
            sel = line["sel"]
            wells = [{}, sel] + [conc(**{d: 1.0}) for d in M.DRUGS] + [dict(sel, **conc(**{d: 1.0})) for d in M.DRUGS if d not in sel][:6]
            out[name] = (self.start("sequencing", name), self.start("viability", name, wells, 5, [1, 2, 3, 4, 5]))
        return out

    def collect_lines(self, pending: dict):
        for name, (s, _) in pending.items():
            s = self.collect(s)
            if s:
                self.lines[name]["alterations"] = s["alterations"]
                for a in s["alterations"]:
                    if a not in self.major and a in M.ALTERATIONS and a not in self.minors:
                        self.minors.append(a)
                self.say(f"line {name} (survivors of {self.lines[name]['sel']}): {s['alterations']}")
        for name, (_, v) in pending.items():
            self.collect(v)

    def verify(self, xs: list[np.ndarray], days: float, reads: list, low: bool) -> dict:
        """Culture tumour cells (two wells each) and healthy cells at the regimens' average concentrations."""
        status = self.lab("lab_status")
        left = int(status.get("wells_left", 0))
        days = min(days, float(status.get("days_left", days)))
        xs = xs[:max(0, min(len(xs), (left - 4) // 3))]
        wells, owner = [], []
        for i, x in enumerate(xs):
            wells += [self.well_of(x), self.well_of(x)]
            owner += [i, i]
        if low:
            for i, x in enumerate(xs[:3]):
                if left - len(wells) - len(xs) > LINE_COST + 30:
                    wells.append(self.well_of(x, 0.6))
                    owner.append(-1 - i)
        reads = [r for r in reads if r <= days]
        t = self.start("viability", "tumour", wells, days, reads, tag="verify") if wells and days >= 2 else None
        h = self.start("viability", "healthy", [self.well_of(x) for x in xs], min(days, 7), [r for r in (1, 2, 3, 5, 7) if r <= days],
                       tag="verify") if xs and days >= 2 else None
        return {"xs": xs, "owner": owner, "tumour": t, "healthy": h}

    def collect_verify(self, v: dict) -> list[dict]:
        t = self.collect(v["tumour"])
        self.collect(v["healthy"])
        return [{"x": x, "exp": t["id"] if t else None, "wells": [j for j, o in enumerate(v["owner"]) if o == i]}
                for i, x in enumerate(v["xs"])]


class LocalLab:
    """The same tool calls, straight into a Lab object (for tests in the repository; not available to agents)."""

    def __init__(self, lab):
        self.lab = lab

    def __call__(self, tool: str, **a):
        from remission import lab as L
        from remission import patients as P
        lab = self.lab
        try:
            if tool == "patient_chart":
                return P.chart(lab.patient)
            if tool == "lab_status":
                return lab.status()
            if tool == "start_experiment":
                if a["kind"] == "viability":
                    return lab.start_viability(a["cells"], a.get("wells"), a.get("days"), a.get("read_days"))
                if a["kind"] == "phospho":
                    return lab.start_phospho(a["cells"], a.get("wells"), a.get("days"))
                return lab.start_sequencing(a["cells"])
            if tool == "derive_line":
                return lab.derive_line(a["experiment_id"], a["well"])
            if tool == "advance_time":
                return lab.advance_time(a["days"])
            if tool == "get_results":
                return lab.results(a["experiment_id"])
            if tool == "submit_regimen":
                return lab.submit(a["regimen"])
        except L.LabError as exc:
            return {"error": str(exc)}
        return {"error": f"unknown tool {tool}"}
