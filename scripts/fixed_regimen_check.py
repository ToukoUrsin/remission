#!/usr/bin/env python3
"""Can one fixed regimen, submitted without any experiment, score on the benchmark's patients?

    scripts/fixed_regimen_check.py [--procs 16] [--out results/fixed_regimens.json]

For every cohort patient recorded in tasks/suites.json (36 of the test salt, 12 of the development salt) the
reward of each of the 3112 constant-dose regimens of remission/cohort.py is computed on the fast model, with
the patient's hidden parameters. Then:

  blind      per tumour type, the regimen with the best mean over the development patients (this is what
             reference/fixed_regimens.json "hindsight" holds and what the `blind` baseline submits), scored on
             the test-salt patients, which it never saw
  hindsight  the same choice made on the test-salt patients themselves: an upper bound no agent can reach,
             because it needs the hidden parameters of the very patients it is scored on

This needs the repository and the seeds; nothing here is available to an agent.
"""
from __future__ import annotations

import argparse
import json
import sys
from collections import defaultdict
from multiprocessing import Pool
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from remission import cohort as C  # noqa: E402
from remission import patients as P  # noqa: E402


def rewards(seed: int) -> np.ndarray:
    return C.grid_rewards(P.make_patient(seed))


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--procs", type=int, default=16)
    p.add_argument("--out", default=str(ROOT / "results" / "fixed_regimens.json"))
    a = p.parse_args()
    suites = json.loads((ROOT / "tasks" / "suites.json").read_text())
    dev, test = suites["remission-v3-dev"]["patients"], suites["remission-v3-test"]["patients"]
    with Pool(a.procs) as pool:
        R = dict(zip([q["seed"] for q in dev + test], pool.map(rewards, [q["seed"] for q in dev + test])))
    grid = C.grid()

    def by_type(patients):
        out = defaultdict(list)
        for q in patients:
            out[q["archetype"]].append(q)
        return out

    blind = {}
    for kind, qs in by_type(dev).items():
        mean = np.mean([R[q["seed"]] for q in qs], axis=0)
        j = int(np.argmax(mean))
        blind[kind] = {"index": j, "regimen": C.doses_of(grid[j]), "development_patients": len(qs),
                       "mean_on_development": round(float(mean[j]), 3)}
    per_patient = []
    for q in test:
        b = blind.get(q["archetype"])
        per_patient.append({"index": q["index"], "seed": q["seed"], "tumour_type": q["archetype"],
                            "blind_reward": round(float(R[q["seed"]][b["index"]]), 3) if b else 0.0,
                            "best_constant_regimen_for_this_patient": round(float(R[q["seed"]].max()), 3)})
    hindsight = {}
    for kind, qs in by_type(test).items():
        mean = np.mean([R[q["seed"]] for q in qs], axis=0)
        j = int(np.argmax(mean))
        hindsight[kind] = {"regimen": C.doses_of(grid[j]), "test_patients": len(qs), "mean_on_test": round(float(mean[j]), 3)}
    all_mean = np.mean([R[q["seed"]] for q in test], axis=0)
    j = int(np.argmax(all_mean))
    blind_rewards = [x["blind_reward"] for x in per_patient]
    out = {
        "_about": "scripts/fixed_regimen_check.py: constant-dose regimens (84 days) on the fast model, hidden parameters",
        "grid_regimens": int(len(grid)),
        "blind": {"chosen_on": "the 12 development patients, per tumour type", "regimens": blind,
                  "scored_on": "the 36 test-salt cohort patients",
                  "mean_reward": round(float(np.mean(blind_rewards)), 3),
                  "patients_above_0.10": int(sum(r > 0.10 for r in blind_rewards)),
                  "rewards_above_0.10": sorted((r for r in blind_rewards if r > 0.10), reverse=True)},
        "hindsight_on_the_test_patients": {
            "per_tumour_type": hindsight,
            "per_tumour_type_mean_over_patients": round(float(sum(h["mean_on_test"] * h["test_patients"] for h in hindsight.values())
                                                              / len(test)), 3),
            "one_regimen_for_all": {"regimen": C.doses_of(grid[j]), "mean_on_test": round(float(all_mean[j]), 3)}},
        "per_patient": per_patient,
    }
    for b in blind.values():
        del b["index"]
    Path(a.out).write_text(json.dumps(out, indent=1) + "\n")
    print(json.dumps({k: out[k] for k in ("blind", "hindsight_on_the_test_patients")}, indent=1))


if __name__ == "__main__":
    main()
