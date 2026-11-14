"""Remission environment server: an MCP sidecar that owns all hidden state.

Runs as its own container next to the agent's `main` container. The agent only
sees the MCP tools (streamable HTTP at /mcp). Nothing secret is reachable from
the agent side:

  REMISSION_SEED          env var, set for this container only (compose file)
  REMISSION_ADMIN_TOKEN   env var, set for this container only; the verifier
                          has a copy in tests/ (uploaded after the agent stopped)

HTTP surface
  POST /mcp              MCP streamable HTTP (stateless, JSON responses)
  GET  /health           liveness, no state
  GET  /admin/episode    final score + full episode log. Needs
                         `Authorization: Bearer <token>`; 403 otherwise.
                         The first authorised call closes the lab.
"""
from __future__ import annotations

import base64
import hmac
import json
import os
import threading
import time
from typing import Annotated, Any

import anyio
import numpy as np
from mcp.server.fastmcp import FastMCP
from pydantic import BeforeValidator
from starlette.requests import Request
from starlette.responses import JSONResponse

from . import lab as L
from . import model as M
from . import patients as P

SEED = int(os.environ.get("REMISSION_SEED", "0"))
ADMIN_TOKEN = os.environ.get("REMISSION_ADMIN_TOKEN", "")
HOST = os.environ.get("REMISSION_HOST", "0.0.0.0")
PORT = int(os.environ.get("REMISSION_PORT", "8000"))
KEEP_ALIVE_S = int(os.environ.get("REMISSION_KEEP_ALIVE_S", "7200"))
SCORE_REPS = int(os.environ.get("REMISSION_SCORE_REPS", "16"))
try:                                    # a regimen known to cure this patient; shown in the replay, never to the agent
    BEST_KNOWN = json.loads(os.environ.get("REMISSION_BEST_KNOWN", "") or "{}")
except ValueError:
    BEST_KNOWN = {}
LOG_SUMMARY_CHARS = 1500


def brief() -> dict:
    k, _ = M.pk_constants()
    return {
        "mission": (
            "You run a precision-oncology lab for one patient. You have a biopsy of the tumour and matched healthy "
            "cells from the same patient, 28 lab days and 384 wells. Work out how this patient's cells respond, then "
            "submit ONE treatment regimen. After the session the patient is treated with it and followed for one year (364 days). "
            "Cure the tumour without wrecking the patient's normal tissue. Every patient here is refractory: the standard "
            "regimens for the tumour type do not cure them, and at least one regimen of these drugs does."),
        "scoring": {
            "reward": "max(0, healed - harm), between 0 and 1",
            "healed": "share of the lesion's tissue sites that hold healthy cells on day 364 (about 0.74 on day 0; "
                      "an untreated tumour takes the tissue over, healed goes to 0)",
            "harm": "largest fractional loss of the patient's normal tissue at any time during the course "
                    "(drug toxicity to healthy cells elsewhere in the body)",
            "notes": ["The regimen is scored on the patient's tissue, which you never see; no score is shown to you.",
                      "No regimen submitted = untreated = 0.",
                      "A tumour clone that is not driven extinct grows back and takes the tissue over well before day 364."],
        },
        "regimen_format": {
            "example": {"phases": [{"days": 28, "doses": {"RAFi": 1.0, "MEKi": 0.5}}, {"days": 14, "doses": {"RAFi": 1.0}}]},
            "rules": [f"doses are multiples of the labelled daily dose, 0 to {M.DOSE_MAX:g}, given once a day",
                      f"phases run back to back from day 0 and may cover at most {M.TREAT_DAYS} days in total; 1 to 12 phases",
                      "a drug not named in a phase is not given in that phase"],
        },
        "drugs": [{"name": d, "class": CLASS[d], "half_life_h": M.PK[d]["t_half_h"],
                   "avg_concentration_uM_at_dose_1": M.PK[d]["css"]} for d in M.DRUGS],
        "drug_note": "avg_concentration_uM_at_dose_1 is the average concentration the tissue sees at the labelled daily "
                     "dose, in the same units as the wells: a well at that concentration matches dose 1.0 on average.",
        "lab": {
            "clock": f"{L.LAB_DAYS} days; time passes only when you call advance_time; quarter-day resolution",
            "wells": f"{L.WELL_BUDGET} in total; at most {L.MAX_WELLS} per experiment",
            "cells": "'tumour' (the biopsy), 'healthy' (matched healthy cells), and lines you derive ('L1', 'L2', ...)",
            "assays": {
                "viability": f"cells in wells at constant drug concentrations (micromolar), about {M.SEED_CELLS} cells seeded; "
                             f"live-cell counts by imaging on up to {L.MAX_READS} read days; 1 well each; count error about "
                             f"{L.COUNT_CV:.0%}",
                "phospho": f"pERK and pAKT at one time point, relative to untreated cells; 1 well each; error about {L.PHOSPHO_CV:.0%}",
                "sequencing": f"alterations carried by at least {L.SEQ_LOD:.0%} of the cells of a cell source; costs "
                              f"{L.SEQ_COST} wells and takes {L.SEQ_DAYS:g} days",
                "derive_line": f"expand the survivors of one well of a finished viability experiment into a new cell "
                               f"source; costs {L.DERIVE_COST} wells and takes {L.DERIVE_DAYS:g} days",
            },
        },
        "model": {
            "file": "/app/cellmodel.py",
            "what": "The exact equations the wells and the patient's tissue are simulated with: signalling, cell fate, "
                    "pharmacokinetics, wells, tissue. /app/regimen_search.py searches regimens for parameters you give it.",
            "hidden": "This patient's parameter values: the potency of each drug (K), the parameters of the healthy "
                      "cells, of the dominant tumour clone and of any minor clones, and the minor clones' shares.",
            "priors": "cellmodel.py also holds what is known about patients in general: typical potencies and their "
                      "spread (K_TYPICAL, K_LN_SD), toxicity (TOX_TYPICAL), parameter ranges (HEALTHY_RANGE, "
                      "TUMOUR_RANGE) and what each alteration changes (ALTERATIONS).",
            "known": "Tumours often contain minor clones (well under 10 % of the cells) that already carry a resistance "
                     "alteration. The chart's bulk panel does not report clones under 10 %.",
        },
        "tools": ["get_brief", "patient_chart", "lab_status", "start_experiment", "derive_line", "advance_time",
                  "get_results", "submit_regimen"],
    }


CLASS = {
    "RAFi": "BRAF V600 inhibitor",
    "MEKi": "MEK1/2 inhibitor",
    "EGFRi": "EGFR inhibitor",
    "METi": "MET inhibitor",
    "KRASi": "KRAS G12C inhibitor",
    "PI3Ki": "PI3K inhibitor",
    "BCLi": "BCL-2/BCL-xL inhibitor (BH3 mimetic)",
    "MCLi": "MCL-1 inhibitor (BH3 mimetic)",
}


class Episode:
    """One episode per sidecar lifetime: the lab plus an append-only call log."""

    def __init__(self, seed: int) -> None:
        self.patient = P.make_patient(seed)
        self.lab = L.Lab(self.patient)
        self.calls: list[dict[str, Any]] = []
        self.lock = threading.Lock()
        self.started_wall = time.time()
        self.closed_by: str | None = None
        self._final: dict | None = None

    @staticmethod
    def _summarise(result: Any) -> Any:
        text = json.dumps(result, default=str)
        if len(text) <= LOG_SUMMARY_CHARS:
            return result
        return {"_truncated": True, "chars": len(text), "head": text[:LOG_SUMMARY_CHARS]}

    def call(self, tool: str, args: dict[str, Any], fn) -> dict[str, Any]:
        with self.lock:
            day_before = self.lab.day
            t_wall = round(time.time() - self.started_wall, 3)
            if self.closed_by == "verifier":
                result, ok = {"error": "the session is over"}, False
            else:
                try:
                    result, ok = fn(), True
                except L.LabError as exc:
                    result, ok = {"error": str(exc)}, False
                except Exception as exc:  # noqa: BLE001 - a malformed request must not take the lab down
                    result, ok = {"error": f"the lab could not read that request ({type(exc).__name__})"}, False
            self.calls.append({"i": len(self.calls), "tool": tool, "args": self._summarise(args), "ok": ok,
                               "day_before": day_before, "day": self.lab.day, "t_wall": t_wall,
                               "wells_left": self.lab.wells_left, "result": self._summarise(result)})
            if self.lab.regimen is not None and self.closed_by is None:
                self.closed_by = "submit"
            return result

    def finalise(self) -> dict[str, Any]:
        """Verifier hand-off: close the lab, treat the patient, return score + log."""
        with self.lock:
            if self._final is None:
                if self.closed_by is None:
                    self.closed_by = "verifier"
                self._final = {"score": self._score(), "episode": self._episode()}
            return self._final

    def _score(self) -> dict:
        pt = self.patient
        regimen = self.lab.regimen or {"phases": [{"days": 1, "doses": {}}]}
        res = M.simulate_tissue(pt["K"], pt["healthy"], pt["tumour"], pt["frac"], regimen,
                                seed=P.seed_from("score", pt["seed"]), reps=SCORE_REPS, record_every=2)
        self._replay = res
        cured = (res["tumour"][:, -1].sum(axis=1) == 0)
        return {"reward": round(float(res["reward"].mean()), 4), "healed": round(float(res["healed"].mean()), 4),
                "harm": round(float(res["harm"]), 4), "cured": round(float(cured.mean()), 4),
                "submitted": 1.0 if self.lab.regimen is not None else 0.0,
                "wells_used": float(self.lab.well_budget - self.lab.wells_left),
                "lab_days_used": float(self.lab.day), "n_experiments": float(len(self.lab.experiments))}

    def _episode(self) -> dict:
        pt, res = self.patient, self._replay
        frames = [{"day": int(d), "grid": base64.b64encode(g.astype(np.uint8).tobytes()).decode()} for d, g in res["frames"]]
        clones = [{"name": c["name"], "alterations": [P.ALTERATION_LABEL[a] for a in c["alterations"]],
                   "share": float(f)} for c, f in zip(pt["tumour"], pt["frac"])]
        experiments = []
        for e in self.lab.experiments:
            pub = {k: v for k, v in e.items() if not k.startswith("_")}
            experiments.append(pub)
        best = None
        if BEST_KNOWN.get("regimen") and M.validate_regimen(BEST_KNOWN["regimen"]) is None:
            b = M.simulate_tissue(pt["K"], pt["healthy"], pt["tumour"], pt["frac"], BEST_KNOWN["regimen"],
                                  seed=P.seed_from("score", pt["seed"]), reps=4, record_every=2)
            best = {"regimen": BEST_KNOWN["regimen"], "reward": round(float(b["reward"].mean()), 4),
                    "healed": round(float(b["healed"].mean()), 4), "harm": round(float(b["harm"]), 4),
                    "frames": [{"day": int(d), "grid": base64.b64encode(g.astype(np.uint8).tobytes()).decode()}
                               for d, g in b["frames"]],
                    "healthy": b["healthy"].mean(axis=0).round(1).tolist(),
                    "tumour": b["tumour"].mean(axis=0).round(1).tolist(),
                    "normal": np.round(b["normal"], 4).tolist()}
        return {
            "best_known": best,
            "closed_by": self.closed_by, "n_calls": len(self.calls),
            "wall_seconds": round(time.time() - self.started_wall, 3),
            "seed": pt["seed"], "chart": P.chart(pt), "regimen": self.lab.regimen,
            "submitted_day": self.lab.submitted_day,
            "truth": {"archetype": pt["archetype"], "clones": clones, "K": dict(zip(M.DRUGS, map(float, pt["K"]))),
                      "tox": dict(zip(M.DRUGS, map(float, pt["healthy"]["tox"])))},
            "replay": {"lattice": M.LATTICE, "frames": frames,
                       "healthy": res["healthy"].mean(axis=0).round(1).tolist(),
                       "tumour": res["tumour"].mean(axis=0).round(1).tolist(),
                       "normal": np.round(res["normal"], 4).tolist(),
                       "conc": np.round(res["conc"], 5).tolist(),
                       "reward_reps": np.round(res["reward"], 4).tolist()},
            "lines": {k: {"kind": v["kind"], "ready_day": v["ready_day"], "note": v["note"],
                          "frac": None if v["frac"] is None else np.round(v["frac"], 4).tolist()}
                      for k, v in self.lab.lines.items()},
            "experiments": experiments, "lab_log": self.lab.log, "calls": list(self.calls),
        }


EPISODE = Episode(SEED)

mcp = FastMCP(
    "remission",
    instructions=("Remission precision-oncology lab. Call get_brief first. All interaction happens through these "
                  "tools; submit_regimen is final and can be called once."),
    host=HOST, port=PORT, stateless_http=True, json_response=True,
)
try:  # the agent reaches us as http://remission:8000, so the Host header is not localhost
    from mcp.server.transport_security import TransportSecuritySettings
    mcp.settings.transport_security = TransportSecuritySettings(enable_dns_rebinding_protection=False)
except Exception:  # noqa: BLE001
    pass


async def _play(tool: str, args: dict[str, Any], fn) -> dict[str, Any]:
    args = {k: v for k, v in args.items() if v is not None}
    return await anyio.to_thread.run_sync(EPISODE.call, tool, args, fn)


def _number_only(v: Any) -> Any:
    if isinstance(v, bool):
        raise ValueError("expected a number, not true/false")
    return v


Num = Annotated[float, BeforeValidator(_number_only)]
Int = Annotated[int, BeforeValidator(_number_only)]
LAB = EPISODE.lab


@mcp.tool()
async def get_brief() -> dict:
    """The mission, how the regimen is scored, the regimen format, the drugs, the lab's assays and budget, and
    where the cell model is. Free; call it first."""
    return await _play("get_brief", {}, brief)


@mcp.tool()
async def patient_chart() -> dict:
    """The patient's chart: tumour type and the bulk sequencing report of the biopsy. Free."""
    return await _play("patient_chart", {}, lambda: P.chart(EPISODE.patient))


@mcp.tool()
async def lab_status() -> dict:
    """Lab day, days and wells left, the cell sources you can use and every experiment with its status. Free."""
    return await _play("lab_status", {}, LAB.status)


@mcp.tool()
async def start_experiment(kind: str, cells: str, wells: list[dict[str, float]] | None = None, days: Num | None = None,
                           read_days: list[float] | None = None) -> dict:
    """Start an experiment on a cell source ('tumour', 'healthy' or a derived line such as 'L1').

    kind = "viability": `wells` is a list with one {drug: micromolar} object per well ({} = no drug); the cells
        grow for `days` at those constant concentrations and live cells are counted on each of `read_days`
        (days after the start, at most 8; default: the last day). Costs 1 well per well.
    kind = "phospho": same `wells`; pERK and pAKT are measured once, after `days` (0.25 = 6 hours), relative to
        untreated cells. Costs 1 well per well.
    kind = "sequencing": alterations carried by at least 2 % of the cells of `cells`. Costs 16 wells, takes 2 days.

    Results appear as the lab clock passes the read days (see advance_time, get_results)."""
    def run():
        if kind == "viability":
            return LAB.start_viability(cells, wells, days, read_days)
        if kind == "phospho":
            return LAB.start_phospho(cells, wells, days)
        if kind == "sequencing":
            return LAB.start_sequencing(cells)
        raise L.LabError("kind must be 'viability', 'phospho' or 'sequencing'")
    return await _play("start_experiment", {"kind": kind, "cells": cells, "wells": wells, "days": days,
                                            "read_days": read_days}, run)


@mcp.tool()
async def derive_line(experiment_id: str, well: Int) -> dict:
    """Expand the cells still alive in one well (index into the experiment's `wells`) at the end of a finished
    viability experiment on tumour cells into a new cell source. Costs 8 wells; the line is ready 3 days later."""
    return await _play("derive_line", {"experiment_id": experiment_id, "well": well},
                       lambda: LAB.derive_line(experiment_id, well))


@mcp.tool()
async def advance_time(days: Num) -> dict:
    """Let the lab clock run for `days` (quarter-day steps). The only way time passes. Returns lab_status."""
    return await _play("advance_time", {"days": days}, lambda: LAB.advance_time(days))


@mcp.tool()
async def get_results(experiment_id: str) -> dict:
    """Everything measured so far in one experiment. Free."""
    return await _play("get_results", {"experiment_id": experiment_id}, lambda: LAB.results(experiment_id))


@mcp.tool()
async def submit_regimen(regimen: dict[str, Any]) -> dict:
    """Submit the patient's treatment regimen, e.g.
    {"phases": [{"days": 28, "doses": {"RAFi": 1.0, "MEKi": 0.5}}, {"days": 14, "doses": {"RAFi": 1.0}}]}.
    Final: it can be called once, closes the lab and returns no outcome."""
    return await _play("submit_regimen", {"regimen": regimen}, lambda: LAB.submit(regimen))


@mcp.custom_route("/health", methods=["GET"])
async def health(_: Request) -> JSONResponse:
    return JSONResponse({"ok": True})


@mcp.custom_route("/admin/episode", methods=["GET"])
async def admin_episode(request: Request) -> JSONResponse:
    auth = request.headers.get("authorization", "")
    supplied = auth[7:] if auth.lower().startswith("bearer ") else ""
    if not ADMIN_TOKEN or not hmac.compare_digest(supplied.encode(), ADMIN_TOKEN.encode()):
        return JSONResponse({"error": "forbidden"}, status_code=403)
    return JSONResponse(await anyio.to_thread.run_sync(EPISODE.finalise))


def main() -> None:
    import uvicorn
    uvicorn.run(mcp.streamable_http_app(), host=mcp.settings.host, port=mcp.settings.port,
                log_level=mcp.settings.log_level.lower(), timeout_keep_alive=KEEP_ALIVE_S)


if __name__ == "__main__":
    main()
