"""Reference solution and the scripted baselines. Runs inside the agent container, public tools only.

    REMISSION_POLICY=scientist   the reference scientist (default)
    REMISSION_POLICY=guideline   the standard-of-care regimen for the chart's tumour type, no experiment
    REMISSION_POLICY=blind       the constant regimen that is best on average for the chart's tumour type
                                 (chosen with hindsight over a development cohort), no experiment
    REMISSION_POLICY=screen      a three-day viability screen, then the combination that kills most
                                 tumour cells while sparing healthy cells, at label doses
    REMISSION_POLICY=hedge       sequencing of the biopsy and nothing else, then the regimen that the public
                                 model and regimen search rate best on average over every tumour it could
                                 be: what sequencing shows plus, in turn, each resistant clone that could
                                 hide below the detection limit, at the public typical parameter values
"""
from __future__ import annotations

import json
import os
import sys

sys.path.insert(0, "/app")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from labclient import call  # noqa: E402

import cellmodel as M  # noqa: E402

# By the chart's driver alteration: (guideline regimen, hindsight-best constant regimen on the development cohort).
FIXED = json.loads(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "fixed_regimens.json")).read())


def tumour_type(chart: dict) -> str:
    alts = chart["bulk_sequencing"]["alterations"]
    if "BRAF V600E" in alts:
        return "BRAF_MEL" if "melanoma" in chart["tissue"] else "BRAF_CRC"
    return "KRAS" if "KRAS G12C" in alts else "EGFR"


def constant(doses: dict, days: int = M.TREAT_DAYS) -> dict:
    return {"phases": [{"days": days, "doses": doses}]}


def guideline():
    kind = tumour_type(call("patient_chart"))
    print(call("submit_regimen", regimen=constant(FIXED[kind]["guideline"])))


def blind():
    kind = tumour_type(call("patient_chart"))
    print(call("submit_regimen", regimen=constant(FIXED[kind]["hindsight"])))


def screen():
    """What a quick drug screen suggests: three days, every single drug and pair at the label concentration."""
    css = {d: M.PK[d]["css"] for d in M.DRUGS}
    combos = [()] + [(d,) for d in M.DRUGS] + [(a, b) for i, a in enumerate(M.DRUGS) for b in M.DRUGS[i + 1:]]
    wells = [{d: css[d] for d in c} for c in combos]
    t = call("start_experiment", kind="viability", cells="tumour", wells=wells, days=3, read_days=[3])
    h = call("start_experiment", kind="viability", cells="healthy", wells=wells, days=3, read_days=[3])
    call("advance_time", days=3)
    tc = [row[-1] for row in call("get_results", experiment_id=t["id"])["counts"]]
    hc = [row[-1] for row in call("get_results", experiment_id=h["id"])["counts"]]
    best, best_kill = (), 0.0
    for c, tn, hn in zip(combos[1:], tc[1:], hc[1:]):
        kill = 1.0 - tn / max(tc[0], 1)
        if hn >= 0.7 * hc[0] and kill > best_kill:
            best, best_kill = c, kill
    print("screen picks", best, "tumour kill", round(best_kill, 3))
    print(call("submit_regimen", regimen=constant({d: 1.0 for d in best})))


def hedge():
    """Cover everything that might be there instead of finding out what is."""
    import regimen_search as RS
    import scientist as S
    major = list(call("patient_chart")["bulk_sequencing"]["alterations"])
    seq = call("start_experiment", kind="sequencing", cells="tumour")
    call("advance_time", days=2)
    found = {a["alteration"]: a["cell_fraction"] for a in call("get_results", experiment_id=seq["id"])["alterations"]}
    sci = S.Scientist(call, seed=0)
    sci.major, sci.rejected, sci.tested = major, set(), set()
    sci.share_obs = {a: f for a, f in found.items() if a not in major and a in M.ALTERATIONS}
    sci.minors = list(sci.share_obs)
    worlds = [st.param_set(z) for st, z in (sci.hypothesis(lab, 0.01) for lab in sci.unseen())]
    x, _, _ = RS.search(worlds, iters=14, pop=256, seed=0)
    x, r = RS.polish(worlds, x)
    reg = RS.regimen_from_x(x)
    print("hedge over", len(worlds), "possible tumours; expected reward", round(float(r), 3), json.dumps(reg))
    print(call("submit_regimen", regimen=reg))


def scientist():
    import scientist as S
    out = S.Scientist(call, seed=0).run()
    print(json.dumps({k: v for k, v in out.items() if k != "notes"}, default=str))


if __name__ == "__main__":
    policy = os.environ.get("REMISSION_POLICY", "scientist")
    {"scientist": scientist, "guideline": guideline, "blind": blind, "screen": screen, "hedge": hedge}[policy]()
