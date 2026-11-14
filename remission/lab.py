"""The lab: what can be done with one patient's cells in 28 days and 384 wells.

Hidden side. Every assay is an abstract call that returns numbers simulated
with the public model (model.py) and the patient's hidden parameters, plus
measurement noise. Nothing here is a laboratory protocol.
"""
from __future__ import annotations

import numpy as np

from . import model as M
from . import patients as P

WELL_BUDGET = 384            # four 96-well plates
LAB_DAYS = 28                # the regimen must be submitted by the end of day 28
MAX_WELLS = 96               # per experiment
MAX_READS = 8                # imaging time points per viability experiment
COUNT_CV = 0.10              # cell-count error (sd of ln) per read
PHOSPHO_CV = 0.10            # pERK / pAKT error (sd of ln) ...
PHOSPHO_FLOOR = 0.02         # ... plus this much additive error
SEQ_COST, SEQ_DAYS, SEQ_LOD, SEQ_CV = 16, 2.0, 0.02, 0.15
DERIVE_COST, DERIVE_DAYS, DERIVE_MIN_CELLS = 8, 3.0, 30
CONC_MAX = 100.0             # micromolar


class LabError(ValueError):
    """Raised for requests the lab refuses; the message is shown to the agent."""


def _q(x: float) -> float:
    """Clock resolution is a quarter day."""
    return round(float(x) * 4.0) / 4.0


class Lab:
    def __init__(self, patient: dict, seed: int | None = None, well_budget: int = WELL_BUDGET,
                 lab_days: float = LAB_DAYS):
        self.patient = patient
        self.seed = int(patient["seed"] if seed is None else seed)
        self.day = 0.0
        self.lab_days = float(lab_days)
        self.wells_left = int(well_budget)
        self.well_budget = int(well_budget)
        self.experiments: list[dict] = []
        self.lines: dict[str, dict] = {
            "tumour": {"kind": "biopsy", "ready_day": 0.0, "frac": np.asarray(patient["frac"], dtype=float),
                       "note": "the patient's tumour biopsy"},
            "healthy": {"kind": "healthy", "ready_day": 0.0, "frac": None,
                        "note": "matched healthy cells from the same patient"},
        }
        self.regimen: dict | None = None
        self.submitted_day: float | None = None
        self.log: list[dict] = []

    # ------------------------------------------------------------------ helpers
    def _rng(self, n: int):
        return np.random.default_rng([self.seed, 104729, n])

    def _check_open(self):
        if self.regimen is not None:
            raise LabError("the regimen has been submitted; the lab is closed")

    def _conc_matrix(self, wells) -> np.ndarray:
        if not isinstance(wells, list) or not wells:
            raise LabError("wells must be a non-empty list, one {drug: micromolar} object per well ({} = no drug)")
        if len(wells) > MAX_WELLS:
            raise LabError(f"at most {MAX_WELLS} wells per experiment")
        conc = np.zeros((len(wells), M.ND))
        for i, w in enumerate(wells):
            if w is None:
                w = {}
            if not isinstance(w, dict):
                raise LabError(f"well {i}: expected an object {{drug: micromolar}}")
            for name, v in w.items():
                if name not in M.IDX:
                    raise LabError(f"well {i}: unknown drug {name!r}; drugs are {', '.join(M.DRUGS)}")
                try:
                    v = float(v)
                except (TypeError, ValueError):
                    raise LabError(f"well {i}: concentration of {name} must be a number (micromolar)")
                if not (0.0 <= v <= CONC_MAX):
                    raise LabError(f"well {i}: concentration of {name} must be between 0 and {CONC_MAX} micromolar")
                conc[i, M.IDX[name]] = v
        return conc

    def _cells(self, name: str):
        """(clones, frac0, capacity multiple) for a cell source that is ready."""
        if name not in self.lines:
            raise LabError(f"unknown cells {name!r}; available: {', '.join(self.lines)}")
        line = self.lines[name]
        if line["ready_day"] > self.day + 1e-9:
            raise LabError(f"{name} is still expanding; ready on day {line['ready_day']:g}")
        if line["kind"] == "healthy":
            return [self.patient["healthy"]], np.ones(1), M.CAP_HEALTHY
        return list(self.patient["tumour"]), line["frac"], M.CAP_TUMOUR

    def _spend(self, n: int):
        if n > self.wells_left:
            raise LabError(f"that needs {n} wells; {self.wells_left} left")
        self.wells_left -= n

    def _duration(self, days) -> float:
        try:
            d = _q(days)
        except (TypeError, ValueError):
            raise LabError("days must be a number")
        if d < 0.25:
            raise LabError("an experiment lasts at least 0.25 days")
        if self.day + d > self.lab_days + 1e-9:
            raise LabError(f"it is day {self.day:g}; an experiment of {d:g} days would end after day {self.lab_days:g}")
        return d

    def _record(self, tool: str, **kw):
        self.log.append({"day": self.day, "tool": tool, **kw})

    # ------------------------------------------------------------------ assays
    def start_viability(self, cells: str, wells, days, read_days=None) -> dict:
        """Grow cells in wells at constant drug concentrations; count live cells by imaging."""
        self._check_open()
        clones, frac0, cap = self._cells(cells)
        conc = self._conc_matrix(wells)
        d = self._duration(days)
        if read_days is None:
            reads = [d]
        else:
            if not isinstance(read_days, list) or not read_days:
                raise LabError("read_days must be a list of days after the start, e.g. [1, 3, 7]")
            reads = sorted({_q(r) for r in read_days})
            if len(reads) > MAX_READS:
                raise LabError(f"at most {MAX_READS} read days per experiment")
            if reads[0] < 0 or reads[-1] > d + 1e-9:
                raise LabError("read_days must lie between 0 and the experiment's length")
        self._spend(conc.shape[0])
        n = len(self.experiments)
        rng = self._rng(n)
        sim = M.simulate_wells(self.patient["K"], clones, frac0, conc, d, rng=rng, read_days=reads + [d], cap_mult=cap)
        total = sim["N"][:-1].sum(axis=2)                                    # [reads, wells]
        obs = rng.poisson(total) * np.exp(rng.normal(0.0, COUNT_CV, total.shape))
        exp = {"id": f"E{n + 1}", "kind": "viability", "cells": cells, "start_day": self.day, "end_day": self.day + d,
               "wells": [{M.DRUGS[j]: float(conc[i, j]) for j in range(M.ND) if conc[i, j] > 0} for i in range(conc.shape[0])],
               "read_days": reads, "counts": np.rint(obs).astype(int).T.tolist(),   # [wells][reads]
               "_final": sim["N"][-1], "_tumour": self.lines[cells]["kind"] != "healthy"}
        self.experiments.append(exp)
        self._record("start_experiment", id=exp["id"], kind="viability", cells=cells, wells=conc.shape[0], days=d,
                     drugs=sorted({k for w in exp["wells"] for k in w}))
        return self._public(exp)

    def start_phospho(self, cells: str, wells, days) -> dict:
        """Treat cells, then measure pERK and pAKT at one time point (relative to untreated cells)."""
        self._check_open()
        clones, frac0, cap = self._cells(cells)
        conc = self._conc_matrix(wells)
        d = self._duration(days)
        self._spend(conc.shape[0])
        n = len(self.experiments)
        rng = self._rng(n)
        full = np.vstack([conc, np.zeros((1, M.ND))])                        # last well: untreated control
        sim = M.simulate_wells(self.patient["K"], clones, frac0, full, d, rng=rng, read_days=[d], cap_mult=cap)
        p = M.pack(clones)
        E0, A0 = M.baseline(p)
        N = np.maximum(sim["N"][0], 1e-9)
        w = N / N.sum(axis=1, keepdims=True)
        perk = (w * sim["E"][0] * E0).sum(axis=1)
        pakt = (w * sim["A"][0] * A0).sum(axis=1)

        def noisy(x):
            rel = x[:-1] / max(x[-1], 1e-9)
            return np.maximum(rel * np.exp(rng.normal(0.0, PHOSPHO_CV, rel.shape)) + rng.normal(0.0, PHOSPHO_FLOOR, rel.shape), 0.0)

        exp = {"id": f"E{n + 1}", "kind": "phospho", "cells": cells, "start_day": self.day, "end_day": self.day + d,
               "wells": [{M.DRUGS[j]: float(conc[i, j]) for j in range(M.ND) if conc[i, j] > 0} for i in range(conc.shape[0])],
               "pERK": np.round(noisy(perk), 3).tolist(), "pAKT": np.round(noisy(pakt), 3).tolist()}
        self.experiments.append(exp)
        self._record("start_experiment", id=exp["id"], kind="phospho", cells=cells, wells=conc.shape[0], days=d,
                     drugs=sorted({k for w_ in exp["wells"] for k in w_}))
        return self._public(exp)

    def start_sequencing(self, cells: str) -> dict:
        """Deep sequencing of a cell source: alterations carried by at least 2 % of its cells."""
        self._check_open()
        if cells not in self.lines:
            raise LabError(f"unknown cells {cells!r}; available: {', '.join(self.lines)}")
        if self.lines[cells]["ready_day"] > self.day + 1e-9:
            raise LabError(f"{cells} is still expanding; ready on day {self.lines[cells]['ready_day']:g}")
        d = self._duration(SEQ_DAYS)
        self._spend(SEQ_COST)
        n = len(self.experiments)
        rng = self._rng(n)
        found = []
        if self.lines[cells]["kind"] != "healthy":
            frac = np.asarray(self.lines[cells]["frac"], dtype=float)
            share: dict[str, float] = {}
            for f, clone in zip(frac, self.patient["tumour"]):
                for a in clone["alterations"]:
                    share[a] = share.get(a, 0.0) + float(f)
            for a, f in share.items():
                if f >= SEQ_LOD:
                    est = min(1.0, f * float(np.exp(rng.normal(0.0, SEQ_CV)))) if f < 0.9 else 1.0
                    found.append({"alteration": P.ALTERATION_LABEL[a], "cell_fraction": round(est, 3)})
        exp = {"id": f"E{n + 1}", "kind": "sequencing", "cells": cells, "start_day": self.day, "end_day": self.day + d,
               "alterations": found, "detection_limit": SEQ_LOD}
        self.experiments.append(exp)
        self._record("start_experiment", id=exp["id"], kind="sequencing", cells=cells, wells=SEQ_COST, days=d, drugs=[])
        return self._public(exp)

    def derive_line(self, experiment_id: str, well: int) -> dict:
        """Expand the cells that are alive in one well at the end of a finished viability experiment."""
        self._check_open()
        exp = self._find(experiment_id)
        if exp["kind"] != "viability":
            raise LabError("lines can only be derived from viability experiments")
        if exp["end_day"] > self.day + 1e-9:
            raise LabError(f"{exp['id']} ends on day {exp['end_day']:g}; it is day {self.day:g}")
        if not exp["_tumour"]:
            raise LabError("healthy cells are already available as 'healthy'")
        try:
            well = int(well)
        except (TypeError, ValueError):
            raise LabError("well must be an index into the experiment's wells")
        if not 0 <= well < len(exp["wells"]):
            raise LabError(f"well must be between 0 and {len(exp['wells']) - 1}")
        if self.day + DERIVE_DAYS > self.lab_days + 1e-9:
            raise LabError(f"expansion takes {DERIVE_DAYS:g} days; it would end after day {self.lab_days:g}")
        self._spend(DERIVE_COST)
        alive = np.asarray(exp["_final"][well], dtype=float)
        if alive.sum() < DERIVE_MIN_CELLS:
            self._record("derive_line", source=exp["id"], well=well, ok=False)
            return {"ok": False, "reason": f"fewer than {DERIVE_MIN_CELLS} live cells in that well; nothing grew out",
                    "wells_left": self.wells_left}
        # drug-free expansion shifts the mix by each clone's untreated growth rate
        p = M.pack(self.patient["tumour"])
        E0, A0 = M.baseline(p)
        kp0, kd0 = M.fate_rates(p, E0, A0, np.zeros((len(alive), M.ND)))
        grown = alive * np.exp((kp0 - kd0) * DERIVE_DAYS)
        name = f"L{sum(1 for v in self.lines.values() if v['kind'] == 'derived') + 1}"
        self.lines[name] = {"kind": "derived", "ready_day": self.day + DERIVE_DAYS, "frac": grown / grown.sum(),
                            "note": f"survivors of {exp['id']} well {well} ({exp['cells']}; "
                                    + (", ".join(f"{k} {v:g} uM" for k, v in exp['wells'][well].items()) or "no drug") + ")"}
        self._record("derive_line", source=exp["id"], well=well, ok=True, line=name)
        return {"ok": True, "line": name, "ready_day": self.day + DERIVE_DAYS, "wells_left": self.wells_left}

    # ------------------------------------------------------------------ clock, results, submission
    def advance_time(self, days) -> dict:
        self._check_open()
        try:
            d = _q(days)
        except (TypeError, ValueError):
            raise LabError("days must be a number")
        if d <= 0:
            raise LabError("days must be positive (the clock moves in quarter days)")
        self.day = min(self.day + d, self.lab_days)
        self._record("advance_time", to_day=self.day)
        return self.status()

    def _find(self, experiment_id: str) -> dict:
        for e in self.experiments:
            if e["id"] == str(experiment_id):
                return e
        raise LabError(f"unknown experiment {experiment_id!r}")

    def _public(self, exp: dict) -> dict:
        """What the agent may see of an experiment at the current day."""
        out = {"id": exp["id"], "kind": exp["kind"], "cells": exp["cells"], "start_day": exp["start_day"],
               "end_day": exp["end_day"], "status": "done" if exp["end_day"] <= self.day + 1e-9 else "running",
               "wells_left": self.wells_left, "day": self.day}
        if exp["kind"] == "viability":
            seen = [i for i, r in enumerate(exp["read_days"]) if exp["start_day"] + r <= self.day + 1e-9]
            out["wells"] = exp["wells"]
            out["read_days"] = [exp["read_days"][i] for i in seen]
            out["counts"] = [[row[i] for i in seen] for row in exp["counts"]]
            out["note"] = "counts[well][read]: live cells per well; about 3000 seeded"
        elif exp["kind"] == "phospho":
            out["wells"] = exp["wells"]
            if out["status"] == "done":
                out["pERK"], out["pAKT"] = exp["pERK"], exp["pAKT"]
                out["note"] = "relative to untreated cells of the same source harvested at the same time"
        elif exp["kind"] == "sequencing" and out["status"] == "done":
            out["alterations"] = exp["alterations"]
            out["note"] = f"alterations carried by at least {SEQ_LOD:.0%} of the cells; cell_fraction is approximate"
        return out

    def results(self, experiment_id: str) -> dict:
        return self._public(self._find(experiment_id))

    def status(self) -> dict:
        return {"day": self.day, "days_left": self.lab_days - self.day, "wells_left": self.wells_left,
                "cells": {k: {"ready_day": v["ready_day"], "note": v["note"]} for k, v in self.lines.items()},
                "experiments": [{"id": e["id"], "kind": e["kind"], "cells": e["cells"], "end_day": e["end_day"],
                                 "status": "done" if e["end_day"] <= self.day + 1e-9 else "running"} for e in self.experiments],
                "regimen_submitted": self.regimen is not None}

    def submit(self, regimen) -> dict:
        self._check_open()
        err = M.validate_regimen(regimen)
        if err:
            raise LabError(err)
        self.regimen = {"phases": [{"days": int(ph["days"]),
                                    "doses": {k: float(v) for k, v in (ph.get("doses") or {}).items() if float(v) > 0}}
                                   for ph in regimen["phases"]]}
        self.submitted_day = self.day
        self._record("submit_regimen", regimen=self.regimen)
        return {"accepted": True, "note": "The regimen is final. The patient is treated after the session; no outcome is reported here."}
