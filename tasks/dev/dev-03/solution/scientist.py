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
SOFT = 10.0                             # cells; counts are compared as ln(count + SOFT)
HIDDEN_GAIN = 8.0                       # a hidden clone is accepted when it lowers the fit's cost by this much
LINE_COST = 8 + 16 + 16                 # wells to expand a line, sequence it and test the drugs on it
PROBE_DAYS = 8                          # length of the cultures that look for clones nobody has seen
SEL_SHARE0 = 0.005                      # share a clone nobody has seen is assumed to start from
SEL_MIN = 0.10                          # share among a well's survivors from which sequencing their line would show it
HYP_PRIOR = 0.12                        # weight of "a clone with this alteration is there" while no line has ruled it out


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
                sig = np.sqrt(0.10 ** 2 + 0.06 ** 2 + (obs + 1.0) / (obs + SOFT) ** 2)
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

    def fit_subset(self, z0: np.ndarray, free: np.ndarray, max_nfev: int = 12):
        """Least-squares fit of the parameters `free` only (indices), the others held at z0. Returns (z, cost)."""
        st = self.st
        h = np.where(st.is_log, 0.02, 0.01 * (st.hi - st.lo))[free]

        def full(x):
            Z = np.repeat(z0[None], np.atleast_2d(x).shape[0], axis=0)
            Z[:, free] = x
            return Z

        def res(x):
            return self.residuals(full(x[None]))[0]

        def jac(x):
            hh = np.where(x + h > st.hi[free], -h, h)
            R = self.residuals(full(np.vstack([x[None], x[None] + np.diag(hh)])))
            return ((R[1:] - R[0]) / hh[:, None]).T

        x0 = np.clip(z0[free], st.lo[free] + 1e-6, st.hi[free] - 1e-6)
        r = least_squares(res, x0, jac=jac, bounds=(st.lo[free], st.hi[free]), loss="soft_l1", f_scale=2.0, max_nfev=max_nfev)
        return full(r.x[None])[0], float(r.cost)

    def cost(self, z: np.ndarray) -> float:
        f = self.residuals(z[None])[0] / 2.0
        return float(0.5 * 4.0 * np.sum(2.0 * (np.sqrt(1.0 + f * f) - 1.0)))

    def scan_hidden(self, rounds: int = 2):
        """Is there a clone nobody has seen? Try every known resistance alteration as a small hidden clone.

        Each hypothesis adds one minor clone (share under 3 %) to the fitted model and fits only that clone's
        share and parameters. The best one is kept when it explains the data much better: survivors that
        neither die nor grow in wells that should be empty are the usual trace.
        """
        for _ in range(rounds):
            base = self.cost(self.z)
            warm = self.warm()
            trials = []
            for lab in M.ALTERATIONS:
                if lab not in self.unseen():
                    continue
                st = Structure(self.major, self.minors + [lab], self.share_obs)
                z0 = st.centre()
                for name, val in warm.items():
                    if name in st.ix:
                        z0[st.ix[name]] = val
                free = np.array([i for i, n in enumerate(st.names) if n.startswith(f"m{len(self.minors)}:")])
                keep, self.st = self.st, st
                try:
                    best = None
                    for share in (np.log(0.004), np.log(0.015)):
                        z1 = z0.copy()
                        z1[st.ix[f"m{len(self.minors)}:share"]] = share
                        z, c = self.fit_subset(z1, free)
                        if best is None or c < best[1]:
                            best = (z, c)
                finally:
                    self.st = keep
                trials.append((base - best[1], lab))
            trials.sort(reverse=True)
            self.say("hidden-clone scan (cost gain): " + ", ".join(f"{lab} {g:.1f}" for g, lab in trials[:4]))
            if not trials or trials[0][0] < HIDDEN_GAIN:
                return
            self.minors.append(trials[0][1])
            self.hidden.append(trials[0][1])
            self.say(f"a hidden clone with {trials[0][1]} explains the data better; added to the model")
            self.refit(starts=1, scan=False)

    def hypothesis(self, lab: str, share: float):
        """The fitted tumour plus one more minor clone carrying `lab` at `share`: (structure, fit vector)."""
        st = Structure(self.major, self.minors + [lab], self.share_obs)
        z = st.centre()
        for name, val in self.warm().items():
            if name in st.ix:
                z[st.ix[name]] = val
        z[st.ix[f"m{len(self.minors)}:share"]] = np.log(share)
        return st, z

    def untested(self) -> list[str]:
        """Alterations a clone nobody has seen could carry and no sequenced line of survivors has ruled out."""
        return [a for a in self.unseen() if a not in self.tested]

    def selection(self, exp: dict) -> dict:
        """What would the survivors of each well of a finished tumour culture hold, if a rare clone were there?

        For every untested alteration: the share a clone carrying it (0.5 % of the biopsy) would have among the
        cells alive at the end of each well, by the fitted model. -> {alteration: share [W]}
        """
        out = {}
        for lab in self.untested():
            st, z = self.hypothesis(lab, SEL_SHARE0)
            K, P, ksh, tox, frac = st.build(z[None])
            N = sim_wells(K, P, ksh, tox, list(range(1, st.C)), frac, exp["conc"], exp["nstep"], exp["reads"],
                          M.SEED_CELLS * M.CAP_TUMOUR)[0][-1, 0]            # [W, clones] at the last read
            out[lab] = N[:, -1] / np.maximum(N.sum(axis=1), 1e-9)
        return out

    def unseen(self) -> list[str]:
        """Alterations a clone nobody has seen could carry: every known one the tumour is not known to have."""
        return [a for a in M.ALTERATIONS if a not in self.major and a not in self.minors and a not in self.rejected
                and not (a in ("BRAF V600E", "EGFR L858R") or (a == "EGFR T790M" and "EGFR L858R" not in self.major))]

    def probe_wells(self, regimens: list[np.ndarray], per: int = 2, most: int = 18) -> list[dict]:
        """Cultures whose survivors would be rich in a clone nobody has seen.

        For every untested alteration the fitted tumour gets a 0.5 % clone carrying it, and a pool of a few
        hundred drug conditions is simulated. Kept: the conditions that kill the known clones but leave that
        clone alive, so that a line grown from the survivors and sequenced would show it.
        """
        labels = self.untested()
        if not labels or most <= 0:
            return []
        rng = np.random.default_rng(23)
        pool = []
        for _ in range(480):
            ds = rng.choice(M.ND, size=int(rng.integers(1, 5)), replace=False)
            pool.append(conc(**{M.DRUGS[j]: float(rng.choice([0.5, 1.0, 2.0, 4.0])) for j in ds}))
        for x in regimens:
            pool += [self.well_of(x, 1.0), self.well_of(x, 2.0)]
        cw = np.array([conc_vec(w) for w in pool])
        nstep = int(round(PROBE_DAYS / M.DT))
        score = {}
        for lab in labels:
            st, z = self.hypothesis(lab, SEL_SHARE0)
            K, P, ksh, tox, frac = st.build(z[None])
            N = sim_wells(K, P, ksh, tox, list(range(1, st.C)), frac, cw, nstep, [nstep], M.SEED_CELLS * M.CAP_TUMOUR)[0][-1, 0]
            total = N.sum(axis=1)
            share = N[:, -1] / np.maximum(total, 1e-9)
            score[lab] = np.where((total >= 60) & (N[:, -1] >= 25), share, 0.0)        # enough alive to grow a line from
        chosen: dict[int, list[str]] = {}
        left = set(labels)
        while left and len(chosen) < most:                                 # fewest wells that cover every alteration
            gain = np.array([sum(min(score[a][j], 0.5) for a in left if score[a][j] >= 2 * SEL_MIN) for j in range(len(pool))])
            for j in chosen:
                gain[j] = 0.0
            j = int(np.argmax(gain))
            if gain[j] <= 0:
                break
            chosen[j] = [a for a in left if score[a][j] >= 2 * SEL_MIN]
            left -= set(chosen[j])
        order = list(chosen)
        for lab in labels:                                                 # a second, different condition for each
            if len(order) >= most:
                break
            for j in np.argsort(-score[lab]):
                if score[lab][j] < 2 * SEL_MIN:
                    break
                if j not in chosen and all(set(pool[j]) != set(pool[i]) for i in chosen if lab in chosen[i]):
                    chosen[int(j)] = [lab]
                    order.append(int(j))
                    break
        self.say("selection wells for clones nobody has seen: " + ("; ".join(
            "+".join(f"{d[:-1]} {v / CSS[M.IDX[d]]:g}x" for d, v in pool[j].items()) + " (" + ", ".join(a.split()[0] for a in chosen[j]) + ")"
            for j in order[:8]) or "no condition would enrich one") + (f"; nothing for {sorted(a.split()[0] for a in left)}" if left else ""))
        return [pool[j] for j in order]

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
    def sets(self) -> tuple[list[dict], np.ndarray]:
        """The worlds a regimen has to work in, with weights: the ensemble of fits, plus (once `robust` is on)
        the fitted tumour with one more 1 % clone for every alteration that no line of survivors has ruled out."""
        sets = [self.st.param_set(z) for z in self.ens]
        hyp = []
        if self.robust:
            for lab in self.untested():
                st, z = self.hypothesis(lab, 0.01)
                hyp.append(st.param_set(z))
        wh = min(HYP_PRIOR, 0.6 / max(len(hyp), 1))
        w = np.array([(1.0 - wh * len(hyp)) / len(sets)] * len(sets) + [wh] * len(hyp))
        return sets + hyp, w

    def expected(self, xs: np.ndarray) -> np.ndarray:
        sets, w = self.sets()
        _, rew = RS.score_plans(sets, w, RS.plans_from_x(xs), shaped=False)
        return rew

    def candidates(self, n: int, extra: list | None = None) -> list[np.ndarray]:
        """Constant-dose regimens with the best mean reward over the ensemble of fits, best first."""
        sets, w = self.sets()
        fast = self.budget == "fast"
        ex = np.array(extra) if extra else None
        x, r, elite = RS.search(sets, weights=w, iters=9 if fast else 14, pop=160 if fast else 256, elite=24,
                                seed=int(self.rng.integers(1 << 30)), extra_x=ex, single_phase=True)
        x, r = RS.polish(sets, x, weights=w, rounds=2 if fast else 3, single_phase=True)
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
        excess = np.log((obs + 30.0) / (pred + 30.0))
        found = [(float(excess[i]), int(i)) for i in range(len(obs))
                 if excess[i] > 1.0 and obs[i] >= 80 and obs[i] > 1.15 * prev[i]]
        return sorted(found, reverse=True)

    def enriched(self, exp: dict, label: str) -> list[tuple[float, int]]:
        """Wells of a tumour culture whose survivors the fit says are mostly the inferred clone `label`."""
        K, P, ksh, tox, frac = self.st.build(self.z[None])
        N = sim_wells(K, P, ksh, tox, list(range(1, self.st.C)), frac, exp["conc"], exp["nstep"], exp["reads"],
                      M.SEED_CELLS * M.CAP_TUMOUR)[0][-1, 0]            # [W, clones] at the last read
        share = N[:, 1 + self.minors.index(label)] / np.maximum(N.sum(axis=1), 1e-9)
        obs = exp["counts"][:, -1]
        found = [(float(share[i]) + min(1.0, obs[i] / 300.0), int(i)) for i in range(len(obs)) if obs[i] >= 45 and share[i] >= 0.6]
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
        self.major, self.minors, self.share_obs, self.hidden, self.sought, self.rejected = major, [], {}, [], set(), set()
        self.tested, self.running, self.robust = set(), [], False
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
        self.expand(tpanel, 2, sel=2)
        pool: list[dict] = []
        first_choice = self.candidates(5)
        verA = self.verify(first_choice, days=8, reads=[2, 4, 6, 8], low=True)
        spare = self.wells_left() - LINE_COST * (1 + sum(1 for ln in self.lines.values() if not ln.get("studied"))) - 16
        pw = self.probe_wells(first_choice[:3], most=max(0, min(10, spare)))
        probe = self.start("viability", "tumour", pw, PROBE_DAYS, [4, 6, PROBE_DAYS], tag="probe") if pw else None
        self.lab("advance_time", days=3)

        # ---- day 13: sequence the lines; test the drugs on those that grew back
        pending = self.study_lines()
        self.lab("advance_time", days=2)

        # ---- day 15: the lines are named; a line that holds a new clone gets its drug test
        self.collect_lines(pending)
        self.lab("advance_time", days=3)

        # ---- day 18: refit with the clones named; from here a regimen must also hold against every clone that
        #      no line has ruled out. Survivors of the selection wells and of a candidate that grew back: more lines
        self.sweep()
        pool += self.collect_verify(verA)
        probe = self.collect(probe)
        self.robust = True
        self.refit(starts=2)
        self.expand(probe, 1, sel=2)
        if verA["tumour"] is not None and not any(not ln.get("studied") for ln in self.lines.values()):
            self.expand(next(e for e in self.data if e["id"] == verA["tumour"]["id"]), 1)
        cand = self.candidates(6, extra=[c["x"] for c in pool])
        cand = [x for x in cand if all(np.abs(x[:M.ND] - c["x"][:M.ND]).max() > 0.12 for c in pool)][:5]
        verB = self.verify(cand, days=9, reads=[2, 4, 6, 8, 9], low=False)
        self.lab("advance_time", days=3)

        # ---- day 21: late lines: sequencing and the drug test at once
        self.sweep()
        late = self.study_lines(drugs=True)
        self.lab("advance_time", days=2)

        # ---- day 23: a clone named this late changes the choice; culture what the new fit prefers
        verC = None
        if self.collect_lines(late):
            self.refit(starts=1, max_nfev=30)
            seen = pool + [{"x": x} for x in verB["xs"]]
            cand = self.candidates(4, extra=[c["x"] for c in seen])
            cand = [x for x in cand if all(np.abs(x[:M.ND] - c["x"][:M.ND]).max() > 0.12 for c in seen)][:3]
            verC = self.verify(cand, days=4, reads=[1, 2, 3, 4], low=False)
        self.lab("advance_time", days=4)

        # ---- day 27: final fit and choice
        self.sweep()
        pool += self.collect_verify(verB)
        if verC:
            pool += self.collect_verify(verC)
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
        best["x"], best["expected"] = self.step_down(best["x"], best["expected"])
        regimen = RS.regimen_from_x(best["x"])
        self.say(f"submitting {regimen} (expected reward {best['expected']:.2f})")
        r = self.lab("submit_regimen", regimen=regimen)
        return {"regimen": regimen, "expected": best["expected"], "explained": best["explained"],
                "accepted": bool(r.get("accepted", False)), "minors": self.minors, "inferred": self.hidden, "notes": self.notes}

    def step_down(self, x: np.ndarray, expected: float) -> tuple[np.ndarray, float]:
        """Induction, then a lower maintenance dose: normal tissue is lost slowly, the tumour is killed fastest
        while it is large. Kept only when every fitted world agrees that it is at least as good."""
        sets, w = self.sets()
        plans = RS.plans_from_x(x[None])
        base = np.array([RS.score_plans([ps], [1.0], plans, shaped=False)[1][0] for ps in sets])
        best = (x, expected, 0.0)
        for share in (0.1, 0.2, 0.35):
            x0 = x.copy()
            x0[M.ND:2 * M.ND], x0[-1] = x0[:M.ND], share
            x1, _ = RS.polish(sets, x0, weights=w, rounds=2, seed=int(self.rng.integers(1 << 30)))
            new = np.array([RS.score_plans([ps], [1.0], RS.plans_from_x(x1[None]), shaped=False)[1][0] for ps in sets])
            gain = float(w @ (new - base))
            if gain > best[2] + 0.005 and np.all(new >= base - 0.02):
                best = (x1, float(w @ new), gain)
        if best[2] > 0:
            # the fast model is well mixed; on the lattice a few cells of a resistant clone can outlast a short
            # induction. Keep the step-down only if lesions simulated from the fits heal as often as before.
            worlds = sets[:3] + sets[len(self.ens):][:2]
            def healed(xx):
                reg = RS.regimen_from_x(xx)
                return float(np.mean([RS.lattice_reward(ps, reg, seed=17 + i, reps=4) for i, ps in enumerate(worlds)]))
            was, now = healed(x), healed(best[0])
            if now < was - 0.01:
                self.say(f"step-down dropped: lesions simulated from the fits score {now:.2f} with it, {was:.2f} without")
                return x, expected
            self.say(f"stepping the dose down after an induction phase adds {best[2]:.2f}: {RS.regimen_from_x(best[0])}")
        return best[0], best[1]

    # ------------------------------------------------------------------ steps
    def refit(self, starts: int = 2, max_nfev: int = 40, scan: bool = True):
        warm = self.warm()
        self.st = Structure(self.major, self.minors, self.share_obs)
        self.fit(starts=starts, warm=warm, max_nfev=max_nfev)
        if scan:
            self.scan_hidden()

    def wells_left(self) -> int:
        return int(self.lab("lab_status").get("wells_left", 0))

    @staticmethod
    def wells_text(well: dict) -> str:
        return "+".join(f"{d[:-1]} {v / CSS[M.IDX[d]]:.2g}x" for d, v in well.items()) or "no drug"

    def expand(self, exp: dict | None, k: int, sel: int = 0):
        """Expand the survivors of up to k wells that grew back although the fit says they should not, and of up
        to `sel` wells whose survivors would show a clone nobody has seen."""
        if exp is None:
            return
        found = self.outgrowth(exp)
        picks = [(i, None) for i in (self.pick_wells(exp, found, k) if found else [])]
        self.say(f"{exp['id']}: unexplained regrowth in {len(found)} wells; expanding wells {[i for i, _ in picks]}")
        # a clone the fit only infers is isolated too: its line names it and measures what it responds to
        if exp["cells"] == "tumour":
            for label in self.hidden:
                if label in self.sought or len(picks) >= k:
                    continue
                wells = [i for _, i in self.enriched(exp, label) if all(i != j for j, _ in picks)]
                if wells:
                    self.sought.add(label)
                    picks.append((wells[0], label))
                    self.say(f"{exp['id']}: well {wells[0]} should hold mostly the inferred {label} clone; expanding it")
        picks = [(i, label, [], True) for i, label in picks]
        sel = max(0, min(sel, 3 - len(picks)))
        # survivors of the hardest selection: a rare resistant clone is a large share of them, so sequencing
        # their line names it, or rules it out, even where the counts show nothing unusual
        if exp["cells"] == "tumour" and sel > 0:
            share, obs = self.selection(exp), exp["counts"][:, -1]
            left, taken = set(share), {i for i, *_ in picks}
            for _ in range(sel):
                gain = np.array([sum(min(share[a][i], 0.5) for a in left if share[a][i] >= SEL_MIN) if obs[i] >= 45 and i not in taken else 0.0
                                 for i in range(len(obs))])
                i = int(np.argmax(gain)) if len(gain) else 0
                if not left or gain[i] <= 0:
                    break
                covers = [a for a in left if share[a][i] >= SEL_MIN]
                picks.append((i, None, covers, False))
                taken.add(i)
                left -= set(covers)
                self.say(f"{exp['id']}: the survivors of well {i} ({self.wells_text(exp['wells'][i])}) would show a clone with "
                         + ", ".join(a.split()[0] for a in covers) + "; expanding them")
            if left:
                self.say(f"{exp['id']}: no well would enrich {sorted(a.split()[0] for a in left)}")
        for i, label, covers, study in picks:
            if self.wells_left() < LINE_COST + 12:
                self.say("not enough wells left to study another line")
                break
            r = self.lab("derive_line", experiment_id=exp["id"], well=i)
            if r.get("ok"):
                self.lines[r["line"]] = {"well": i, "alterations": None, "sel": exp["wells"][i], "ready": r["ready_day"],
                                         "sought": label, "covers": covers, "drugs": study}
                self.tested |= set(covers)

    def study_lines(self, drugs: bool = False) -> dict:
        """Sequence every new line and test the drugs on it. Returns the running experiments by line."""
        day = float(self.lab("lab_status").get("day", 0))
        out = {}
        for name, line in self.lines.items():
            if line.get("studied") or line["ready"] > day + 1e-9:
                continue
            line["studied"] = True
            out[name] = (self.start("sequencing", name), self.test_line(name) if line.get("drugs", True) or drugs else None)
        return out

    def test_line(self, name: str) -> dict | None:
        """The drugs, alone and on top of what the line survived, on a derived line."""
        line = self.lines[name]
        sel = line["sel"]
        days = min(5.0, float(self.lab("lab_status").get("days_left", 5)))
        if line.get("tested") or days < 3 or self.wells_left() < 30:
            return None
        line["tested"] = True
        wells = [{}, sel] + [conc(**{d: 1.0}) for d in M.DRUGS] + [dict(sel, **conc(**{d: 1.0})) for d in M.DRUGS if d not in sel][:6]
        e = self.start("viability", name, wells, days, [r for r in (1, 2, 3, 4, 5) if r <= days])
        if e:
            self.running.append(e)
        return e

    def sweep(self):
        """File every drug test on a line that has finished."""
        for e in list(self.running):
            if self.collect(e):
                self.running.remove(e)

    def collect_lines(self, pending: dict) -> bool:
        """Read the lines' sequencing. Returns True when it named a clone that was not in the model."""
        news = False
        for name, (s, _) in pending.items():
            line = self.lines[name]
            if line.get("alterations") is not None:
                continue
            s = self.collect(s)
            if not s:
                continue
            line["alterations"] = s["alterations"]
            new = [a for a in s["alterations"] if a not in self.major and a in M.ALTERATIONS and a not in self.minors]
            self.minors += new
            for a in s["alterations"]:
                if a in self.hidden:
                    self.hidden.remove(a)
            self.say(f"line {name} (survivors of {self.wells_text(line['sel'])}): {s['alterations']}")
            for want in ([line["sought"]] if line.get("sought") else []) + list(line.get("covers") or []):
                if want in s["alterations"]:
                    continue
                self.rejected.add(want)
                if want in self.hidden:
                    self.say(f"the line does not carry {want}; dropping that hypothesis")
                    self.hidden.remove(want)
                    self.minors.remove(want)
            if line.get("covers"):
                self.say(f"ruled out by line {name}: " + (", ".join(a.split()[0] for a in line["covers"] if a not in s["alterations"]) or "nothing"))
            if new:
                news = True
                self.test_line(name)
        self.sweep()
        return news

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
