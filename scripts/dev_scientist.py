"""Run the reference scientist in-process on patients and score its regimen on the lattice.

    python scripts/dev_scientist.py SEED [SEED ...] [--fast] [--quiet] [--json OUT]

Development tool: it builds the Lab object directly, so it needs the hidden side
of the repository. The benchmark itself runs the scientist through MCP.
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "reference"))

from remission import cohort as C  # noqa: E402
from remission import lab as L  # noqa: E402
from remission import model as M  # noqa: E402
from remission import patients as P  # noqa: E402
import scientist as S  # noqa: E402


def lattice(pt, regimen, reps=8):
    res = M.simulate_tissue(pt["K"], pt["healthy"], pt["tumour"], pt["frac"], regimen,
                            seed=P.seed_from("score", pt["seed"]), reps=reps)
    return float(res["reward"].mean()), float(res["healed"].mean()), float(res["harm"])


def one(seed: int, fast: bool = False, quiet: bool = False, grid: bool = True) -> dict:
    pt = P.make_patient(seed)
    lab = L.Lab(pt)
    t = time.time()
    lines = []
    sci = S.Scientist(S.LocalLab(lab), log=(lines.append if quiet else print), seed=seed, budget="fast" if fast else "full")
    try:
        out = sci.run()
    except Exception as exc:  # noqa: BLE001
        import traceback
        traceback.print_exc()
        out = {"regimen": None, "error": repr(exc), "minors": []}
    took = time.time() - t
    reg = lab.regimen
    rew, healed, harm = lattice(pt, reg) if reg else (0.0, 0.0, 0.0)
    res = {"seed": seed, "archetype": pt["archetype"], "true_minors": [c["name"][6:] for c in pt["tumour"][1:]],
           "shares": [round(float(f), 4) for f in pt["frac"][1:]], "found_minors": out.get("minors"),
           "regimen": reg, "reward": round(rew, 3), "healed": round(healed, 3), "harm": round(harm, 3),
           "expected": out.get("expected"), "sterile": out.get("sterile"), "wells_left": lab.wells_left,
           "seconds": round(took, 1), "error": out.get("error"), "log": lines}
    if grid:
        R = C.grid_rewards(pt)
        j = int(R.argmax())
        res["grid_best"] = round(float(R[j]), 3)
        res["grid_regimen"] = C.doses_of(C._GRID[j])
        if reg:
            plans = M.expand_regimen(reg)[None]
            res["mean_model"] = round(float(M.simulate_tissue_mean(pt["K"], pt["healthy"], pt["tumour"], pt["frac"], plans)["reward"][0]), 3)
    return res


def _one_quiet(seed):
    return one(seed, False, True)


def _one_quiet_fast(seed):
    return one(seed, True, True)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("seeds", nargs="*", type=int)
    ap.add_argument("--salt", help="take the patients from this salt's cohort instead (needs --pkl from a cohort scan)")
    ap.add_argument("--pkl", help="pickle of cohort.screen results (list of dicts with i, seed, eligible)")
    ap.add_argument("--first", type=int, default=24)
    ap.add_argument("--skip", type=int, default=0)
    ap.add_argument("--procs", type=int, default=1)
    ap.add_argument("--fast", action="store_true")
    ap.add_argument("--quiet", action="store_true")
    ap.add_argument("--json")
    a = ap.parse_args()
    seeds = list(a.seeds)
    if a.pkl:
        import pickle
        el = sorted((o["i"], o["seed"]) for o in pickle.load(open(a.pkl, "rb")) if o["eligible"])
        seeds = [s for _, s in el[a.skip:a.skip + a.first]]
    rows = []
    out = open(a.json, "w") if a.json else None

    def done(r):
        rows.append(r)
        if out:
            out.write(json.dumps(r) + "\n")
            out.flush()
        print(json.dumps({k: v for k, v in r.items() if k != "log"}), flush=True)

    if a.procs > 1:
        from multiprocessing import Pool
        with Pool(a.procs) as pool:
            for r in pool.imap_unordered(_one_quiet_fast if a.fast else _one_quiet, seeds):
                done(r)
    else:
        for s in seeds:
            done(one(s, a.fast, a.quiet))
    rw = np.array([r["reward"] for r in rows]); gb = np.array([r.get("grid_best", np.nan) for r in rows])
    print(f"n {len(rows)}  reward mean {rw.mean():.3f}  cured(>=0.4) {(rw >= 0.4).mean():.2f}  grid best mean {np.nanmean(gb):.3f}  "
          f"seconds median {np.median([r['seconds'] for r in rows]):.0f}")
