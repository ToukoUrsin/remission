#!/usr/bin/env python3
"""Remission viewer server: the static replay viewer plus "Treat a patient", a public sandbox.

    python3 viewer/serve.py [--port 8377] [--host 127.0.0.1]

Serves viewer/ as static files and two endpoints:

  GET  /api/patients   the 12 development patients of tasks/suites.json ("remission-v3-dev"):
                       chart, guideline regimen, drugs and limits. No hidden truth.
  POST /api/treat      {"patient": "dev-04", "regimen": {"phases": [{"days": 84, "doses": {"RAFi": 1.0}}]}}
                       treats the simulated patient, follows the lesion to day 364 and returns the
                       score, the replay of lesion 0, the daily curves, and the truth about the clones.
  GET  /healthz        liveness

Stdlib + numpy + the repo's `remission` package. Public on purpose, so it is defensive: small bodies
only, at most 3 simulations at once, an LRU cache of results, no stack traces to the client, and
nothing outside viewer/ is ever served. The test patients ("remission-v3-test") are never loaded.
"""
from __future__ import annotations

import argparse
import base64
import gzip
import json
import logging
import mimetypes
import os
import sys
import threading
import time
import urllib.parse
import zlib
from collections import OrderedDict
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path

import numpy as np

VIEWER = Path(__file__).resolve().parent
REPO = Path(os.environ.get("REMISSION_REPO", VIEWER.parent)).resolve()
sys.path.insert(0, str(REPO))

from remission import model as M  # noqa: E402
from remission import patients as P  # noqa: E402

SUITE = "remission-v3-dev"                 # the only suite this server ever reads
MAX_BODY = 16 * 1024
MAX_SIMS = 3
SLOT_WAIT_S = 6.0
CACHE_SIZE = 160
REPS = 4
SEQ_LOD = 0.02                             # the lab's deep sequencing reports alterations in >= 2 % of cells
UI_PHASES = 3

log = logging.getLogger("remission.serve")

CLASS = {
    "RAFi": "BRAF V600 inhibitor",
    "MEKi": "MEK1/2 inhibitor",
    "EGFRi": "EGFR inhibitor",
    "METi": "MET inhibitor",
    "KRASi": "KRAS G12C inhibitor",
    "PI3Ki": "PI3K inhibitor",
    "BCLi": "BCL-2/BCL-xL inhibitor",
    "MCLi": "MCL-1 inhibitor",
}
TUMOUR_TYPE = {
    "BRAF_CRC": "colorectal tumour, BRAF-mutant",
    "BRAF_MEL": "melanoma, BRAF-mutant",
    "KRAS": "colorectal tumour, KRAS-mutant",
    "EGFR": "lung tumour, EGFR-mutant",
}


# ----------------------------------------------------------------------------- patients
def clone_label(alterations: list[str], dominant: list[str], index: int) -> str:
    """Short name of a clone by what sets it apart from the dominant clone (as in scripts/build_viewer.py)."""
    if index == 0:
        return "dominant clone"
    own = [a for a in alterations if a not in dominant] or alterations
    a = own[0]
    if "(" in a and ")" in a:
        gene = a[a.index("(") + 1:a.index(")")]
    else:
        gene = a.split()[0]
    if "amplification" in a:
        return f"{gene}-amplified clone"
    if a.endswith(" loss"):
        return f"{gene}-loss clone"
    if "(" in a:
        return f"{gene} {a.split()[-1]} clone"
    return f"{a} clone"


def constant(doses: dict) -> dict:
    return {"phases": [{"days": M.TREAT_DAYS, "doses": dict(doses)}]}


def load_patients() -> tuple[dict, list[dict]]:
    suites = json.loads((REPO / "tasks" / "suites.json").read_text())
    fixed = json.loads((REPO / "reference" / "fixed_regimens.json").read_text())
    table, public = {}, []
    for i, sp in enumerate(suites[SUITE]["patients"]):
        pid = f"dev-{i + 1:02d}"
        pt = P.make_patient(int(sp["seed"]))
        guide = fixed[pt["archetype"]]["guideline"]
        # lesion 0 on day 0, exactly as the simulation starts it, with every tumour cell drawn as one colour
        grid0 = M.init_tissue(pt["frac"], np.random.default_rng(P.seed_from("score", pt["seed"])), REPS)[0]
        masked = np.where(grid0 >= 2, 2, grid0).astype(np.uint8)
        dom = [P.ALTERATION_LABEL[a] for a in pt["tumour"][0]["alterations"]]
        clones = []
        for k, (c, f) in enumerate(zip(pt["tumour"], pt["frac"])):
            alts = [P.ALTERATION_LABEL[a] for a in c["alterations"]]
            clones.append({"label": clone_label(alts, dom, k), "alterations": alts,
                           "own": [a for a in alts if a not in dom] if k else alts,
                           "share": round(float(f), 5), "sequencing_sees": bool(f >= SEQ_LOD), "on_chart": k == 0})
        table[pid] = {"pt": pt, "best": sp["best"], "best_reward": float(sp["best_reward"]), "clones": clones,
                      "guideline": constant(guide)}
        public.append({
            "id": pid, "archetype": pt["archetype"], "tumour_type": TUMOUR_TYPE[pt["archetype"]],
            "chart": P.chart(pt),
            "guideline": [{"name": "Guideline regimen", "regimen": constant(guide)}],
            "day0": {"L": M.LATTICE, "grid": base64.b64encode(zlib.compress(masked.tobytes(), 9)).decode()},
        })
    return table, public


PATIENTS, PUBLIC = load_patients()
PATIENTS_JSON = json.dumps({
    "patients": PUBLIC,
    "drugs": [{"name": d, "class": CLASS[d]} for d in M.DRUGS],
    "limits": {"dose_min": 0.0, "dose_max": M.DOSE_MAX, "dose_step": 0.25, "treat_days": M.TREAT_DAYS,
               "follow_days": M.FOLLOW_DAYS, "max_phases": UI_PHASES, "api_max_phases": 12},
    "full_dose": constant({d: 1.0 for d in M.DRUGS}),
    "reps": REPS, "sequencing_lod": SEQ_LOD,
    "about": "Simulated patients: a computational model of tumour and healthy cell populations. No real patient data.",
}, separators=(",", ":")).encode()


# ----------------------------------------------------------------------------- treatment
class BadRequest(Exception):
    pass


def canonical(regimen: dict) -> dict:
    """The validated regimen in one spelling: drug order fixed, zero doses dropped, numbers normalised."""
    phases = []
    for ph in regimen["phases"]:
        doses = (ph.get("doses") or {})
        phases.append({"days": int(ph["days"]),
                       "doses": {d: round(float(doses[d]), 4) for d in M.DRUGS if d in doses and float(doses[d]) > 0}})
    return {"phases": phases}


def first_gone(counts: np.ndarray) -> int | None:
    """counts [reps, days]: first day from which no replicate lesion holds a cell, else None."""
    alive = counts.sum(axis=0) > 0
    if alive[-1]:
        return None
    idx = np.nonzero(alive)[0]
    return 0 if idx.size == 0 else int(idx[-1]) + 1


def simulate(pid: str, regimen: dict) -> dict:
    rec = PATIENTS[pid]
    pt = rec["pt"]
    t0 = time.time()
    res = M.simulate_tissue(pt["K"], pt["healthy"], pt["tumour"], pt["frac"], regimen,
                            seed=P.seed_from("score", pt["seed"]), reps=REPS, record_every=2)
    T = res["tumour"]                                  # [reps, days+1, clones]
    H = res["healthy"]
    days = [int(d) for d, _ in res["frames"]]
    blob = zlib.compress(b"".join(g.astype(np.uint8).tobytes() for _, g in res["frames"]), 9)
    mean_t = T.mean(axis=0)                            # [days+1, clones]
    total = mean_t.sum(axis=1)
    nadir = int(total.argmin())
    normal = res["normal"]
    cured = T[:, -1].sum(axis=1) == 0
    treat_end = min(M.TREAT_DAYS, sum(int(ph["days"]) for ph in regimen["phases"]))
    dosed = M.expand_regimen(regimen).sum(axis=1) > 0
    last_dose = int(np.nonzero(dosed)[0][-1]) + 1 if dosed.any() else 0
    clones = []
    for k, c in enumerate(rec["clones"]):
        col = mean_t[:, k]
        clones.append(dict(c, gone_day=first_gone(T[:, :, k]), end=round(float(col[-1]), 1),
                           nadir=round(float(col.min()), 1), nadir_day=int(col.argmin()),
                           cured_reps=int((T[:, -1, k] == 0).sum())))
    out = {
        "patient": pid, "regimen": regimen,
        "reward": round(float(res["reward"].mean()), 4), "healed": round(float(res["healed"].mean()), 4),
        "harm": round(float(res["harm"]), 4), "cured": round(float(cured.mean()), 4), "reps": REPS,
        "reward_reps": np.round(res["reward"], 4).tolist(),
        "treat_end": treat_end, "last_dose_day": last_dose,
        "replay": {
            "L": M.LATTICE, "days": days, "grids": base64.b64encode(blob).decode(),
            "healthy": [round(float(v), 1) for v in H.mean(axis=0)],
            "tumour": [[round(float(v), 1) for v in mean_t[:, k]] for k in range(mean_t.shape[1])],
            "normal": [round(float(v), 4) for v in normal],
            "conc": np.round(res["conc"], 5).tolist(),
        },
        "facts": {
            "tumour_gone_day": first_gone(T.sum(axis=2)), "healthy_gone_day": first_gone(H),
            "nadir_day": nadir, "nadir": round(float(total[nadir]), 1), "start": round(float(total[0]), 1),
            "end": round(float(total[-1]), 1), "harm_day": int(np.argmin(normal)),
            "min_normal": round(float(normal.min()), 4),
        },
        "truth": {"archetype": pt["archetype"], "clones": clones, "sequencing_lod": SEQ_LOD,
                  "chart_lod": 0.10},
        # the best regimen's reward on these same 4 lesion copies (what a visitor gets when they try it),
        # and its reward in the benchmark (16 copies), from tasks/suites.json
        "best": {"regimen": rec["best"], "reward": rec.get("best_here", round(float(res["reward"].mean()), 4)),
                 "benchmark_reward": rec["best_reward"]},
    }
    log.info("simulated %s in %.1f s", pid, time.time() - t0)
    return out


class Cache:
    """LRU of response bodies; pinned entries (the presets) are never evicted."""

    def __init__(self, size: int) -> None:
        self.size, self.data, self.pinned, self.lock = size, OrderedDict(), {}, threading.Lock()

    def get(self, key):
        with self.lock:
            if key in self.pinned:
                return self.pinned[key]
            if key in self.data:
                self.data.move_to_end(key)
                return self.data[key]
        return None

    def put(self, key, value, pin: bool = False) -> None:
        with self.lock:
            if pin:
                self.pinned[key] = value
                self.data.pop(key, None)
                return
            self.data[key] = value
            self.data.move_to_end(key)
            while len(self.data) > self.size:
                self.data.popitem(last=False)


CACHE = Cache(CACHE_SIZE)
SLOTS = threading.BoundedSemaphore(MAX_SIMS)
INFLIGHT_LOCK = threading.Lock()
INFLIGHT: dict[str, threading.Event] = {}


class Busy(Exception):
    pass


def ensure_best(pid: str, wait: float = SLOT_WAIT_S) -> None:
    """Score the patient's best known regimen on the sandbox's lesion copies once, before anything else."""
    rec = PATIENTS[pid]
    if "best_here" not in rec:
        body = treat(pid, canonical(rec["best"]), wait=wait, pin=True)
        rec["best_here"] = json.loads(gzip.decompress(body))["reward"]


def treat(pid: str, regimen: dict, wait: float = SLOT_WAIT_S, pin: bool = False) -> bytes:
    """gzip-compressed JSON of the treatment, from the cache when it was run before."""
    key = pid + "|" + json.dumps(regimen, sort_keys=True, separators=(",", ":"))
    while True:
        hit = CACHE.get(key)
        if hit is not None:
            return hit
        with INFLIGHT_LOCK:                           # the same request twice at once runs once
            ev = INFLIGHT.get(key)
            mine = ev is None
            if mine:
                ev = INFLIGHT[key] = threading.Event()
        if not mine:
            if not ev.wait(timeout=30):
                raise Busy()
            continue
        try:
            if not SLOTS.acquire(timeout=wait):
                raise Busy()
            try:
                body = gzip.compress(json.dumps(simulate(pid, regimen), separators=(",", ":")).encode(), 6)
            finally:
                SLOTS.release()
            CACHE.put(key, body, pin=pin)
            return body
        finally:
            with INFLIGHT_LOCK:
                INFLIGHT.pop(key, None)
            ev.set()


def parse_treat(raw: bytes) -> tuple[str, dict]:
    try:
        req = json.loads(raw.decode("utf-8"))
    except (UnicodeDecodeError, ValueError):
        raise BadRequest("the body must be JSON") from None
    if not isinstance(req, dict):
        raise BadRequest('send {"patient": "dev-01", "regimen": {"phases": [...]}}')
    pid, regimen = req.get("patient"), req.get("regimen")
    if not isinstance(pid, str) or pid not in PATIENTS:
        raise BadRequest(f"patient must be one of {', '.join(PATIENTS)}")
    err = M.validate_regimen(regimen)
    if err:
        raise BadRequest(err)
    return pid, canonical(regimen)


def prewarm() -> None:
    """Run the presets once at start, so the first visitor's guideline click answers at once."""
    t0 = time.time()
    for pid, rec in PATIENTS.items():
        try:
            ensure_best(pid, wait=600)
            for reg in (rec["guideline"], constant({d: 1.0 for d in M.DRUGS})):
                treat(pid, canonical(reg), wait=600, pin=True)
        except Exception:  # noqa: BLE001
            log.exception("prewarm %s failed", pid)
    log.info("prewarm done in %.0f s: %d presets cached", time.time() - t0, len(CACHE.pinned))


# ----------------------------------------------------------------------------- HTTP
TYPES = {".html": "text/html; charset=utf-8", ".js": "text/javascript; charset=utf-8",
         ".css": "text/css; charset=utf-8", ".woff2": "font/woff2", ".json": "application/json",
         ".png": "image/png", ".svg": "image/svg+xml", ".txt": "text/plain; charset=utf-8",
         ".ico": "image/x-icon", ".md": "text/plain; charset=utf-8"}
STATIC_DENY = {".py", ".pyc"}


class Handler(BaseHTTPRequestHandler):
    server_version = "remission"
    sys_version = ""
    protocol_version = "HTTP/1.1"
    timeout = 20

    def log_message(self, fmt, *args):  # one short line per request
        log.info("%s %s", self.address_string(), fmt % args)

    def _send(self, code: int, body: bytes, ctype: str, headers: dict | None = None, gz: bool = False) -> None:
        if gz and "gzip" not in (self.headers.get("Accept-Encoding") or ""):
            body, gz = gzip.decompress(body), False
        self.send_response(code)
        self.send_header("Content-Type", ctype)
        self.send_header("Content-Length", str(len(body)))
        self.send_header("X-Content-Type-Options", "nosniff")
        self.send_header("Referrer-Policy", "no-referrer")
        if gz:
            self.send_header("Content-Encoding", "gzip")
            self.send_header("Vary", "Accept-Encoding")
        for k, v in (headers or {}).items():
            self.send_header(k, v)
        self.end_headers()
        if self.command != "HEAD":
            self.wfile.write(body)

    def _json(self, code: int, obj, headers: dict | None = None) -> None:
        h = {"Cache-Control": "no-store"}
        h.update(headers or {})
        self._send(code, json.dumps(obj).encode(), "application/json", h)

    def _guard(self, fn) -> None:
        try:
            fn()
        except (BrokenPipeError, ConnectionResetError):
            pass
        except Exception:  # noqa: BLE001 - never a stack trace to the client
            log.exception("request failed: %s %s", self.command, self.path[:200])
            try:
                self._json(500, {"error": "the server could not handle that request"})
            except Exception:  # noqa: BLE001
                pass

    def do_HEAD(self):
        self._guard(self._get)

    def do_GET(self):
        self._guard(self._get)

    def do_POST(self):
        self._guard(self._post)

    def _get(self) -> None:
        path = urllib.parse.urlsplit(self.path).path
        if path == "/healthz":
            return self._json(200, {"ok": True})
        if path == "/api/patients":
            return self._send(200, PATIENTS_JSON, "application/json", {"Cache-Control": "no-cache"})
        if path.startswith("/api/"):
            return self._json(404, {"error": "not found"})
        self._static(path)

    def _static(self, path: str) -> None:
        rel = urllib.parse.unquote(path).lstrip("/")
        if rel == "" or rel.endswith("/"):
            rel += "index.html"
        if "\x00" in rel or "\\" in rel or any(part.startswith(".") for part in rel.split("/") if part):
            return self._not_found()
        root = VIEWER
        try:
            target = (root / rel).resolve()
            target.relative_to(root)
        except (ValueError, OSError):
            return self._not_found()
        if not target.is_file() or target.suffix in STATIC_DENY:
            return self._not_found()
        ctype = TYPES.get(target.suffix) or mimetypes.guess_type(target.name)[0] or "application/octet-stream"
        body = target.read_bytes()
        cache = "public, max-age=604800, immutable" if "/vendor/" in "/" + rel else "no-cache"
        gz = False
        if ctype.startswith(("text/", "application/json")) and len(body) > 1024 \
                and "gzip" in (self.headers.get("Accept-Encoding") or ""):
            body, gz = gzip.compress(body, 6), True
        self._send(200, body, ctype, {"Cache-Control": cache}, gz=gz)

    def _not_found(self) -> None:
        self._send(404, b"not found\n", "text/plain; charset=utf-8", {"Cache-Control": "no-store"})

    def _post(self) -> None:
        path = urllib.parse.urlsplit(self.path).path
        if path != "/api/treat":
            return self._json(404, {"error": "not found"})
        try:
            n = int(self.headers.get("Content-Length") or "-1")
        except ValueError:
            n = -1
        if n < 0:
            return self._json(411, {"error": "Content-Length is required"}, {"Connection": "close"})
        if n > MAX_BODY:
            self.close_connection = True
            return self._json(413, {"error": f"the body may be at most {MAX_BODY} bytes"}, {"Connection": "close"})
        raw = self.rfile.read(n)
        try:
            pid, regimen = parse_treat(raw)
        except BadRequest as exc:
            return self._json(400, {"error": str(exc)})
        try:
            ensure_best(pid)
            body = treat(pid, regimen)
        except Busy:
            return self._json(429, {"error": "the simulator is busy; try again in a few seconds"}, {"Retry-After": "5"})
        self._send(200, body, "application/json", {"Cache-Control": "no-store"}, gz=True)


class Server(ThreadingHTTPServer):
    daemon_threads = True
    request_queue_size = 64


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--host", default=os.environ.get("HOST", "127.0.0.1"))
    ap.add_argument("--port", type=int, default=int(os.environ.get("PORT", "8377")))
    ap.add_argument("--no-prewarm", action="store_true", help="do not run the preset regimens at start")
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
    srv = Server((args.host, args.port), Handler)
    if not args.no_prewarm:
        threading.Thread(target=prewarm, name="prewarm", daemon=True).start()
    log.info("Remission viewer on http://%s:%d/  (Treat a patient: #view=treat)", args.host, args.port)
    try:
        srv.serve_forever()
    except KeyboardInterrupt:
        pass


if __name__ == "__main__":
    main()
