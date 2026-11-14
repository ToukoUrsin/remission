"""Summarise the existing development cohort scan into science/cohort_scan_summary.json (cheap; no re-screening).

    science/.venv/bin/python science/validation/cohort_summary.py [/tmp/cohort_dev.pkl] [/tmp/cohort_v2.json]

Input: the pickle written by the owner's scan (cohort.screen(make_patient(seed_from("dev", i)), keep_rewards=True)
for i < 1600, with the 3-guideline cohort.py of 2026-11-14 ~10:35 UTC). The 5-guideline refinement now in
remission/cohort.py is approximated from the stored grid rewards (fast model only: guidelines 3-5 are grid points);
its lattice re-check was not rerun here.
"""
import json
import os
import pickle
import sys
from collections import Counter

import numpy as np

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, ROOT)
from remission import cohort as C  # noqa: E402
from remission import model as M  # noqa: E402

pkl = sys.argv[1] if len(sys.argv) > 1 else "/tmp/cohort_dev.pkl"
v2 = sys.argv[2] if len(sys.argv) > 2 else "/tmp/cohort_v2.json"
out_rows = pickle.load(open(pkl, "rb"))
G = C.grid()


def gidx(doses):
    v = np.zeros(M.ND)
    for k, x in doses.items():
        v[M.IDX[k]] = x
    hit = np.where((np.abs(G - v) < 1e-9).all(axis=1))[0]
    return int(hit[0]) if len(hit) else None


keep = set(json.load(open(v2))["keep_i"]) if os.path.exists(v2) else None
by = {}
for a in ("BRAF_MEL", "BRAF_CRC", "KRAS", "EGFR"):
    rows = [o for o in out_rows if o["archetype"] == a]
    gl_idx = [gidx(g) for g in C.GUIDELINES[a]]
    cnt = Counter()
    for o in rows:
        why = o.get("why", "eligible")
        if why == "eligible":
            extra = [float(o["grid_rewards"][j]) for j in gl_idx[3:] if j is not None and "grid_rewards" in o]
            if extra and max(extra) >= C.REFRACTORY_BELOW:
                cnt["eligible only with 3 guidelines"] += 1
            else:
                cnt["eligible (5 guidelines, approx.)"] += 1
        elif why.startswith("a guideline regimen works"):
            cnt["a guideline regimen works"] += 1
        elif why == "no grid regimen cures":
            cnt["no grid regimen cures"] += 1
        else:
            cnt["fails lattice checks"] += 1
    n = len(rows)
    by[a] = {"n": n, "counts": dict(cnt), "share": {k: round(v / n, 4) for k, v in cnt.items()},
             "eligible_3_guidelines": sum(1 for o in rows if o.get("eligible")),
             "eligible_v2_keep_list": (sum(1 for o in rows if o.get("eligible") and o["i"] in keep) if keep is not None else None),
             "grid_best_median_all": round(float(np.median([o.get("grid_best", np.nan) for o in rows if "grid_best" in o])), 3)}
tot = Counter(o.get("why", "eligible") for o in out_rows)
summary = {
    "source": f"{pkl} (salt 'dev', i < {len(out_rows)}), {v2}",
    "n_generated": len(out_rows),
    "n_eligible_3_guidelines": tot["eligible"],
    "why_3_guidelines": dict(tot),
    "n_eligible_5_guidelines_fast_model_approx": sum(v["counts"].get("eligible (5 guidelines, approx.)", 0) for v in by.values()),
    "n_eligible_v2_keep_list": len(keep) if keep is not None else None,
    "by_type": by,
    "note": "The 3-guideline numbers are exact cohort.screen results. The 5-guideline count only re-checks guidelines 3-5 on the fast model.",
}
json.dump(summary, open(os.path.join(ROOT, "science", "cohort_scan_summary.json"), "w"), indent=1)
print(json.dumps(summary, indent=1))
