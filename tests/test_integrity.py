"""Reward integrity and lab rules. Fast, no Docker:  .venv/bin/python -m pytest tests -q"""
from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from remission import cohort as C  # noqa: E402
from remission import lab as L  # noqa: E402
from remission import model as M  # noqa: E402
from remission import patients as P  # noqa: E402

SEEDS = [P.seed_from("unit-test", i) for i in range(4)]


def score(pt, regimen, reps=2):
    return M.simulate_tissue(pt["K"], pt["healthy"], pt["tumour"], pt["frac"], regimen, seed=1, reps=reps)


@pytest.mark.parametrize("seed", SEEDS)
def test_untreated_scores_zero(seed):
    pt = P.make_patient(seed)
    res = score(pt, {"phases": [{"days": 1, "doses": {}}]})
    assert float(res["reward"].max()) == 0.0
    assert float(res["harm"]) == 0.0


@pytest.mark.parametrize("seed", SEEDS)
def test_everything_at_full_dose_is_too_toxic(seed):
    pt = P.make_patient(seed)
    res = score(pt, {"phases": [{"days": M.TREAT_DAYS, "doses": {d: M.DOSE_MAX for d in M.DRUGS}}]})
    assert float(res["reward"].mean()) < 0.05


def test_reward_is_bounded_and_harm_counts():
    pt = P.make_patient(SEEDS[0])
    res = score(pt, {"phases": [{"days": 84, "doses": {"BCLi": 2.0, "MCLi": 2.0}}]})
    assert np.all((res["reward"] >= 0) & (res["reward"] <= 1))
    assert np.allclose(res["reward"], np.maximum(0.0, res["healed"] - res["harm"]))


@pytest.mark.parametrize("bad", [
    None, {}, {"phases": []}, {"phases": [{"days": 0, "doses": {}}]},
    {"phases": [{"days": 85, "doses": {"MEKi": 1}}]},
    {"phases": [{"days": 50, "doses": {"MEKi": 1}}, {"days": 50, "doses": {"MEKi": 1}}]},
    {"phases": [{"days": 10, "doses": {"MEKi": 2.01}}]},
    {"phases": [{"days": 10, "doses": {"MEKi": -1}}]},
    {"phases": [{"days": 10, "doses": {"MEKi": float("nan")}}]},
    {"phases": [{"days": 10, "doses": {"MEKi": float("inf")}}]},
    {"phases": [{"days": 10, "doses": {"aspirin": 1}}]},
    {"phases": [{"days": 10, "doses": {"MEKi": "lots"}}]},
    {"phases": [{"days": True, "doses": {"MEKi": 1}}]},
    {"phases": [{"days": 1, "doses": {"MEKi": 1}}] * 13},
])
def test_bad_regimens_are_refused(bad):
    assert M.validate_regimen(bad) is not None
    lab = L.Lab(P.make_patient(SEEDS[0]))
    with pytest.raises(L.LabError):
        lab.submit(bad)
    assert lab.regimen is None


def test_budget_clock_and_single_submission():
    lab = L.Lab(P.make_patient(SEEDS[1]))
    with pytest.raises(L.LabError):
        lab.start_viability("tumour", [{}] * 97, 3)                    # too many wells in one experiment
    for _ in range(4):
        lab.start_viability("tumour", [{}] * 96, 3)
    assert lab.wells_left == 0
    with pytest.raises(L.LabError):
        lab.start_viability("tumour", [{}], 3)                         # budget spent
    with pytest.raises(L.LabError):
        lab.start_sequencing("tumour")
    lab.advance_time(27.5)
    with pytest.raises(L.LabError):
        lab.start_phospho("tumour", [{}], 1)                           # would end after day 28
    lab.advance_time(100)
    assert lab.day == 28
    lab.submit({"phases": [{"days": 10, "doses": {"MEKi": 1}}]})
    with pytest.raises(L.LabError):
        lab.submit({"phases": [{"days": 84, "doses": {"MEKi": 2}}]})   # final means final
    with pytest.raises(L.LabError):
        lab.advance_time(1)


def test_results_appear_only_when_the_clock_passes_them():
    lab = L.Lab(P.make_patient(SEEDS[2]))
    e = lab.start_viability("tumour", [{}, {"MEKi": 0.01}], 6, [2, 4, 6])
    assert e["counts"] == [[], []] and e["status"] == "running"
    lab.advance_time(2)
    assert lab.results(e["id"])["read_days"] == [2]
    s = lab.start_sequencing("tumour")
    assert "alterations" not in lab.results(s["id"])
    with pytest.raises(L.LabError):
        lab.derive_line(e["id"], 0)                                    # not finished yet
    lab.advance_time(4)
    assert lab.results(e["id"])["read_days"] == [2, 4, 6]
    assert "alterations" in lab.results(s["id"])


def test_nothing_hidden_in_what_the_agent_sees():
    pt = P.make_patient(SEEDS[3])
    lab = L.Lab(pt)
    seen = [P.chart(pt), lab.status()]
    e = lab.start_viability("tumour", [{}, {"BCLi": 1.0}], 5, [1, 5])
    ph = lab.start_phospho("tumour", [{"MEKi": 0.01}], 0.25)
    sq = lab.start_sequencing("tumour")
    lab.advance_time(5)
    seen += [lab.results(e["id"]), lab.results(ph["id"]), lab.results(sq["id"]), lab.derive_line(e["id"], 0), lab.status()]
    text = json.dumps(seen)
    for secret in ("_final", "_tumour", '"frac"', "seed", "kshift", "tox", str(pt["seed"]), "archetype"):
        assert secret not in text
    for value in list(pt["K"]) + list(pt["frac"]):
        assert f"{value:.6g}" not in text


def test_minor_clones_below_the_detection_limit_are_not_reported():
    for i in range(40):
        pt = P.make_patient(P.seed_from("unit-seq", i))
        lab = L.Lab(pt)
        s = lab.start_sequencing("tumour")
        lab.advance_time(2)
        found = {a["alteration"] for a in lab.results(s["id"])["alterations"]}
        for clone, f in zip(pt["tumour"][1:], pt["frac"][1:]):
            label = P.ALTERATION_LABEL[clone["alterations"][-1]]
            assert (label in found) == (f >= L.SEQ_LOD)


def test_same_seed_same_patient_same_data_same_score():
    a, b = P.make_patient(SEEDS[0]), P.make_patient(SEEDS[0])
    assert np.allclose(a["K"], b["K"]) and np.allclose(a["frac"], b["frac"])
    la, lb = L.Lab(a), L.Lab(b)
    wells = [{}, {"MEKi": 0.004}, {"BCLi": 0.07, "MCLi": 0.03}]
    assert la.start_viability("tumour", wells, 3, [1, 3]) == lb.start_viability("tumour", wells, 3, [1, 3])
    reg = {"phases": [{"days": 84, "doses": {"MEKi": 1.0, "BCLi": 0.5}}]}
    assert np.array_equal(score(a, reg)["reward"], score(b, reg)["reward"])


def test_cohort_guidelines_fail_and_a_cure_exists():
    """The cohort test does what it says, on the first cohort patient of a test salt."""
    for i in range(400):
        pt = P.make_patient(P.seed_from("unit-cohort", i))
        res = C.screen(pt, lattice_reps=2)
        if res["eligible"]:
            break
    else:
        pytest.fail("no cohort patient in 400 generated patients")
    for g in C.GUIDELINES[pt["archetype"]]:
        assert float(score(pt, C.constant(g), reps=4)["reward"].mean()) < 0.25
    assert float(score(pt, res["best"], reps=4)["reward"].mean()) >= 0.4
    assert not math.isnan(res["best_reward"])
