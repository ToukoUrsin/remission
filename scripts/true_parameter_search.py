#!/usr/bin/env python3
"""What the regimen search reaches with a patient's true parameters and no lab work: the ceiling for a scientist
who learned everything about the patient. Not a run; it reads the hidden parameters.

    scripts/true_parameter_search.py tasks/test-01 [tasks/test-02 ...] [--out results/true_parameter_search.jsonl]

One patient takes about 3 minutes on one core. The regimen is the best of 4 searches (two-phase regimens, the same
search the agent gets as /app/regimen_search.py), scored like a submission: 16 replicates of the tissue lattice.
"""
from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from remission import model as M, oracle as O, patients as P  # noqa: E402


def one(task: Path) -> dict:
    seed = int(re.search(r'REMISSION_SEED: "(\d+)"', (task / "environment" / "docker-compose.yaml").read_text()).group(1))
    pt = P.make_patient(seed)
    best = None
    for restart in range(4):
        r = O.best_regimen([pt], lattice_check=pt, seed=restart * 17, top=6, reps=6, restarts=2)
        if best is None or r["lattice_reward"] > best["lattice_reward"]:
            best = r
    res = M.simulate_tissue(pt["K"], pt["healthy"], pt["tumour"], pt["frac"], best["regimen"],
                            seed=P.seed_from("score", pt["seed"]), reps=16)
    return {"task": task.name, "seed": seed, "reward": round(float(res["reward"].mean()), 4),
            "healed": round(float(np.mean(res["healed"])), 4), "harm": round(float(np.mean(res["harm"])), 4),
            "regimen": best["regimen"]}


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("tasks", nargs="+")
    p.add_argument("--out", default=str(ROOT / "results" / "true_parameter_search.jsonl"))
    a = p.parse_args()
    rows = []
    for t in a.tasks:
        rows.append(one(Path(t)))
        print(rows[-1]["task"], rows[-1]["reward"])
    Path(a.out).write_text("".join(json.dumps(r) + "\n" for r in rows))
    print(f"mean {sum(r['reward'] for r in rows) / len(rows):.3f} -> {a.out}")


if __name__ == "__main__":
    main()
