"""Summarise the cohort scans into science/cohort_scan_summary.json (cheap; screens at most a few missing patients).

    science/.venv/bin/python science/validation/cohort_summary.py [TEST_PKL [DEV_PKL]]

Inputs: lists of remission.cohort.screen() results (all six rules, current model) plus "i", the patient's index in
its salt: patients make_patient(seed_from(salt, i)) for i = 0..n-1. Defaults: /tmp/v3/merged_scrT.pkl (salt
remission-v3-test) and /tmp/v3/merged_scrD.pkl (salt remission-v3-dev). An index missing from a pickle is taken
from the row with the same seed when there is one (two indices can draw the same seed), else screened here.

The per-patient record of the cohort is tasks/suites.json; this file holds only counts.
"""
import json
import os
import pickle
import sys
from collections import Counter
from datetime import datetime, timezone

import numpy as np

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, ROOT)
from remission import cohort as C  # noqa: E402
from remission import patients as P  # noqa: E402

TYPES = ("BRAF_MEL", "BRAF_CRC", "KRAS", "EGFR")
# screen()'s reasons of exclusion, in the order screen() checks them, and the rule each one belongs to.
RULES = [
    ("rule 1", "a guideline regimen works"),
    ("rule 2", "no grid regimen cures"),
    ("rule 3", "the untreated tumour does not take the lesion"),
    ("rule 4", "a label-dose drug or pair works"),
    ("rule 5", "sequencing of the biopsy shows every clone"),
    ("rule 5", "the answer for the visible clones also cures the hidden one"),
]
LATTICE = [
    ("rule 2", "the grid's best regimen fails on the lattice"),
    ("rule 1", "a guideline regimen works on the lattice"),
    ("rule 4 or 5", "a label-dose pair or the answer for the visible clones works on the lattice"),
    ("rule 3", "the untreated tumour does not take the lesion on the lattice"),
]
RULE6 = "a regimen that hedges against every clone that could be hiding works"


def load(path, salt):
    rows = pickle.load(open(path, "rb"))
    by_i = {r["i"]: r for r in rows}
    by_seed = {r["seed"]: r for r in rows}
    n = max(by_i) + 1
    filled = []
    for i in range(n):
        seed = P.seed_from(salt, i)
        if i in by_i:
            assert by_i[i]["seed"] == seed, (salt, i)
            continue
        if seed in by_seed:
            by_i[i] = dict(by_seed[seed], i=i)
            filled.append({"i": i, "how": f"same seed as index {by_seed[seed]['i']}"})
        else:
            by_i[i] = dict(C.screen(P.make_patient(seed)), i=i)
            filled.append({"i": i, "how": "screened here"})
    return [by_i[i] for i in range(n)], filled, os.path.getmtime(path)


def stats(x):
    x = np.asarray(x, dtype=float)
    if not len(x):
        return None
    return {"median": round(float(np.median(x)), 3), "min": round(float(x.min()), 3), "max": round(float(x.max()), 3)}


def summarise(rows):
    why = Counter(r.get("why", "eligible") for r in rows)
    known = {w for _, w in RULES + LATTICE} | {RULE6, "eligible"}
    assert set(why) <= known, set(why) - known
    left = len(rows)
    funnel = [{"step": "generated", "left": left}]
    for rule, w in RULES:
        left -= why[w]
        funnel.append({"step": rule, "why": w, "excluded": why[w], "left": left})
    for rule, w in LATTICE:
        left -= why[w]
        funnel.append({"step": "lattice re-check, " + rule, "why": w, "excluded": why[w], "left": left})
    r5 = [r for r in rows if r.get("rules_1_to_5")]
    cohort = [r for r in rows if r["eligible"]]
    assert left == len(r5) and len(cohort) == len(r5) - why[RULE6]
    funnel.append({"step": "rule 6", "why": RULE6, "excluded": why[RULE6], "left": len(cohort)})

    def real(h):           # what rule 6 tests: the better of the fast model and the lattice on the real tumour
        return max(h["reward"], h["lattice"])

    rule6 = {}
    for group, rr in (("cohort", cohort), ("rejected", [r for r in r5 if not r["eligible"]])):
        g = {"n": len(rr)}
        for hedge in ("typical", "strong"):
            hs = [h for r in rr for h in r["hedges"] if h["hedge"] == hedge]
            g[hedge] = {"expected": stats([h["expected"] for h in hs]), "real": stats([real(h) for h in hs]),
                        "real_at_least_0.10": sum(real(h) >= C.REFRACTORY_BELOW for h in hs)}
        rule6[group] = g
    hidden = Counter(a for r in cohort for a in r["hidden"])
    return {
        "generated": len(rows),
        "rules_1_to_5": len(r5),
        "cohort": len(cohort),
        "one_in": round(len(rows) / len(cohort)) if cohort else None,
        "funnel": funnel,
        "by_type": {a: {"generated": sum(r["archetype"] == a for r in rows),
                        "rules_1_to_5": sum(r["archetype"] == a for r in r5),
                        "cohort": sum(r["archetype"] == a for r in cohort)} for a in TYPES},
        "hidden_clone": dict(hidden.most_common()),
        "two_hidden_clones": sum(len(r["hidden"]) >= 2 for r in cohort),
        "best_reward": stats([r["best_reward"] for r in cohort]),
        "rule6": rule6,
    }


def main():
    srcs = [("remission-v3-test", sys.argv[1] if len(sys.argv) > 1 else "/tmp/v3/merged_scrT.pkl"),
            ("remission-v3-dev", sys.argv[2] if len(sys.argv) > 2 else "/tmp/v3/merged_scrD.pkl")]
    out = {"about": "Cohort scans with all six rules of remission/cohort.py on the current model. funnel: patients "
                    "excluded at each step, in the order cohort.screen() checks them; the lattice re-check repeats "
                    "rules 1 to 5 on the tissue lattice. rule6: the two hedge regimens' reward in their own "
                    "hypothetical tumours (expected) and on the real tumour (real, the better of fast model and "
                    "lattice), for cohort patients and for those rule 6 rejects. Per patient: tasks/suites.json.",
           "thresholds": {"refractory_below": C.REFRACTORY_BELOW, "curable_from": C.CURABLE_FROM,
                          "untreated_below": C.UNTREATED_BELOW, "seq_lod": C.SEQ_LOD}}
    every, scanned = [], []
    for salt, path in srcs:
        rows, filled, mtime = load(path, salt)
        scanned.append(mtime)
        out[salt] = dict(summarise(rows), source=path, filled_in=filled)
        every += rows
    out["overall"] = summarise(every)
    out["scanned_utc"] = datetime.fromtimestamp(min(scanned), timezone.utc).strftime("%Y-%m-%d %H:%M")
    out["scanned_mtime"] = min(scanned)
    json.dump(out, open(os.path.join(ROOT, "science", "cohort_scan_summary.json"), "w"), indent=1)
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
