#!/usr/bin/env python3
"""Why did each run score what it scored? Reads Harbor job directories, rebuilds every
patient from its seed and names the decision that cost the reward.

    scripts/analyze_runs.py "LABEL=JOB_DIR" ["LABEL=JOB_DIR" ...]
        [--reference "Reference scientist"] [--out results/run_analysis.jsonl]
        [--summary results/summary.json] [--md results/summary.md]

A run's outcome is first sorted into one of six classes, from what happened in the scored
lesion (not from what the agent said):

    solved             within 0.08 of the best regimen known for the patient, or better
    no_regimen         nothing was submitted: the patient went untreated
    too_toxic          the tumour is gone, but the drugs cost so much normal tissue that less than half of the
                       best known reward is left
    missed_clone       the tumour came back from a minor resistant clone
    tumour_survived    the main clone was never cleared (too short, too little, wrong drugs)
    suboptimal         cured with avoidable damage: at least half of the best known reward, more than 0.08 below it

Runs that say nothing about the agent are left out and listed under `excluded` in the summary: trials that
never reached the verifier (a harness failure) and sessions that ended on a model API error before a regimen
was submitted.

Then the class is tied to what the agent did in the lab: did it sequence, did it grow out
survivors and name them, how long did it watch a culture, did the submitted combination
ever meet healthy cells in a well. Counterfactuals ("the same drugs at half the dose would
have scored 0.61") are run on the same lattice the score comes from, with the hidden
parameters: this script needs the repository, it is not something an agent can run.
"""
from __future__ import annotations

import argparse
import json
import sys
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from remission import lab as L  # noqa: E402
from remission import model as M  # noqa: E402
from remission import patients as P  # noqa: E402

CSS = dict(zip(M.DRUGS, M.CSS_LABEL))
SOLVED_WITHIN = 0.08
TITLES = {
    "solved": "Cured the patient",
    "no_regimen": "No regimen submitted",
    "too_toxic": "Cleared the tumour, poisoned the patient",
    "missed_clone": "Relapse from a resistant clone",
    "tumour_survived": "The main clone was never cleared",
    "suboptimal": "Cured, with avoidable damage",
}


# ----------------------------------------------------------------------------- scoring helpers
def lattice(pt: dict, regimen: dict, reps: int = 4) -> dict:
    res = M.simulate_tissue(pt["K"], pt["healthy"], pt["tumour"], pt["frac"], regimen,
                            seed=P.seed_from("analysis", pt["seed"]), reps=reps)
    return {"reward": float(res["reward"].mean()), "healed": float(res["healed"].mean()), "harm": float(res["harm"])}


def mean_model(pt: dict, regimens: list[dict]) -> dict:
    plans = np.stack([M.expand_regimen(r) for r in regimens])
    return M.simulate_tissue_mean(pt["K"], pt["healthy"], pt["tumour"], pt["frac"], plans)


def scaled(regimen: dict, s: float) -> dict:
    return {"phases": [{"days": ph["days"], "doses": {d: round(min(M.DOSE_MAX, v * s), 3) for d, v in ph["doses"].items()}}
                       for ph in regimen["phases"]]}


def extended(regimen: dict) -> dict:
    """The same phases, the last one held to the end of the treatment window."""
    phases = [dict(ph) for ph in regimen["phases"]]
    used = sum(ph["days"] for ph in phases)
    phases[-1] = {"days": phases[-1]["days"] + M.TREAT_DAYS - used, "doses": phases[-1]["doses"]}
    return {"phases": phases}


def without(regimen: dict, drug: str) -> dict:
    return {"phases": [{"days": ph["days"], "doses": {d: v for d, v in ph["doses"].items() if d != drug}} for ph in regimen["phases"]]}


def main_phase(regimen: dict) -> dict:
    """Average daily dose per drug over the treated days (what a constant culture can mimic)."""
    days = sum(ph["days"] for ph in regimen["phases"])
    out: dict[str, float] = defaultdict(float)
    for ph in regimen["phases"]:
        for d, v in ph["doses"].items():
            out[d] += v * ph["days"] / days
    return {d: v for d, v in out.items() if v >= 0.03}


def show(regimen: dict | None) -> str:
    if not regimen:
        return "nothing"
    parts = []
    for ph in regimen["phases"]:
        doses = " + ".join(f"{d} {v:g}" for d, v in ph["doses"].items() if v > 0) or "no drug"
        parts.append(f"{doses} for {ph['days']} d")
    return "; then ".join(parts)


def pct(x: float) -> str:
    return f"{100 * x:.1f} %" if x < 0.1 else f"{100 * x:.0f} %"


# ----------------------------------------------------------------------------- what the agent did
def lab_facts(ep: dict) -> dict:
    exps = ep.get("experiments") or []
    calls = ep.get("calls") or []
    read = {}                                   # experiment id -> lab day of the first read after it ended
    for c in calls:
        if c.get("tool") == "get_results" and c.get("ok"):
            eid = (c.get("args") or {}).get("experiment_id")
            e = next((e for e in exps if e["id"] == eid), None)
            if e is not None and c.get("day", 0) >= e["end_day"] - 1e-9:
                read.setdefault(eid, c["day"])
    seen = {}                                   # alteration -> (experiment, cells, day) the agent was told about
    for e in exps:
        if e["kind"] == "sequencing" and e["id"] in read:
            for a in e.get("alterations") or []:
                seen.setdefault(a["alteration"], {"experiment": e["id"], "cells": e["cells"], "day": read[e["id"]]})
    via = [e for e in exps if e["kind"] == "viability"]
    tum = [e for e in via if e["cells"] != "healthy"]
    hea = [e for e in via if e["cells"] == "healthy"]
    lines = {k: v for k, v in (ep.get("lines") or {}).items() if v.get("kind") == "derived"}
    return {
        "n_calls": len(calls), "wall_minutes": round(float(ep.get("wall_seconds") or 0) / 60, 1),
        "lab_days_used": ep.get("submitted_day") if ep.get("regimen") else max([c.get("day", 0) for c in calls] or [0]),
        "wells_used": sum(len(e.get("wells") or []) or L.SEQ_COST for e in exps),
        "n_experiments": len(exps),
        "biopsy_sequenced": any(e["kind"] == "sequencing" and e["cells"] == "tumour" and e["id"] in read for e in exps),
        "lines_derived": len(lines),
        "lines_sequenced": sum(1 for e in exps if e["kind"] == "sequencing" and e["cells"] in lines and e["id"] in read),
        "phospho_wells": sum(len(e.get("wells") or []) for e in exps if e["kind"] == "phospho"),
        "tumour_wells": sum(len(e["wells"]) for e in tum), "healthy_wells": sum(len(e["wells"]) for e in hea),
        "healthy_combination_wells": sum(1 for e in hea for w in e["wells"] if len(w) >= 2),
        "longest_tumour_culture_days": max([e["end_day"] - e["start_day"] for e in tum] or [0]),
        "seen_alterations": seen, "read": read,
    }


def cultured(ep: dict, regimen: dict, cells: str, lo: float, hi: float, min_days: float, read: dict) -> bool:
    """Did a well hold this regimen's drugs, each at lo..hi times its average level, on these cells?"""
    want = {d: v * CSS[d] for d, v in main_phase(regimen).items()}
    if not want:
        return False
    for e in ep.get("experiments") or []:
        if e["kind"] != "viability" or e["id"] not in read or e["end_day"] - e["start_day"] < min_days:
            continue
        if (e["cells"] == "healthy") != (cells == "healthy"):
            continue
        for w in e["wells"]:
            if set(want) <= set(w) and all(lo * c <= w[d] <= hi * c for d, c in want.items()):
                if cells == "healthy" or set(w) == set(want):
                    return True
    return False


# ----------------------------------------------------------------------------- the diagnosis
def diagnose(ep: dict, score: dict) -> dict:
    pt = P.make_patient(int(ep["seed"]))
    truth = ep["truth"]
    bk = ep.get("best_known") or {}
    facts = lab_facts(ep)
    if score.get("wells_used") is not None:
        facts["wells_used"] = int(score["wells_used"])          # the lab's own count, which includes derived lines
    read = facts.pop("read")
    seen = facts["seen_alterations"]
    regimen = ep.get("regimen")
    reward, healed, harm, cured = (float(score.get(k, 0.0)) for k in ("reward", "healed", "harm", "cured"))
    best = float(bk.get("reward") or 0.0)
    major = set(truth["clones"][0]["alterations"])
    minors = [{"index": j, "alteration": next((a for a in c["alterations"] if a not in major), c["name"]), "share": c["share"]}
              for j, c in enumerate(truth["clones"]) if j > 0]
    hidden = [m for m in minors if m["share"] < L.SEQ_LOD]
    out = {"reward": round(reward, 4), "healed": round(healed, 3), "harm": round(harm, 3), "cured_share": round(cured, 3),
           "best_known": round(best, 3), "regimen": regimen, "facts": facts,
           "minor_clones": [{"alteration": m["alteration"], "share": round(m["share"], 4), "seen": m["alteration"] in seen,
                             "below_sequencing_limit": m["share"] < L.SEQ_LOD} for m in minors],
           "counterfactuals": {}, "mistakes": []}
    cf, mistakes = out["counterfactuals"], out["mistakes"]

    def add(key: str, title: str, detail: str):
        mistakes.append({"key": key, "title": title, "detail": detail})

    # ---- nothing submitted
    if not regimen:
        how = {"timeout": "the 60 minutes ran out", "agent": "the session ended"}.get(ep.get("closed_by") or "", "the session ended")
        add("no_regimen", "No regimen submitted",
            f"{how} after {facts['n_calls']} tool calls, {facts['n_experiments']} experiments and lab day "
            f"{facts['lab_days_used']:g} of 28; the patient went untreated.")
        out["class"] = "no_regimen"
        out["summary"] = mistakes[0]["detail"]
        return out

    # ---- what happened in the lesion
    tum_end = np.array(ep["replay"]["tumour"][-1], dtype=float)           # mean cells per clone on day 364
    relapse = cured < 0.5
    tested_h = cultured(ep, regimen, "healthy", 0.5, 1e9, 3.0, read)
    tested_t = cultured(ep, regimen, "tumour", 0.5, 2.0, 5.0, read)
    facts["final_regimen_cultured_on_healthy_cells"] = tested_h
    facts["final_regimen_cultured_on_tumour_cells"] = tested_t
    treated = sum(ph["days"] for ph in regimen["phases"] if any(v > 0 for v in ph["doses"].values()))
    facts["days_treated"] = treated

    # ---- counterfactuals on the mean model (ranking), the chosen ones re-scored on the lattice
    scales = [0.25, 0.35, 0.5, 0.7, 1.4, 2.0]
    variants = [scaled(regimen, s) for s in scales] + [extended(regimen)] + [scaled(extended(regimen), s) for s in scales]
    names = [f"same drugs, doses x{s:g}" for s in scales] + ["same doses, held to day 84"] + \
            [f"held to day 84, doses x{s:g}" for s in scales]
    mm = mean_model(pt, variants)
    j = int(np.argmax(mm["reward"]))
    best_var = lattice(pt, variants[j])
    cf["best_rescaling"] = {"change": names[j], "regimen": variants[j], **{k: round(v, 3) for k, v in best_var.items()}}
    drugs = sorted({d for ph in regimen["phases"] for d, v in ph["doses"].items() if v > 0})
    if drugs:
        drop = mean_model(pt, [without(regimen, d) for d in drugs])
        base = mean_model(pt, [regimen])
        cf["harm_without_drug"] = {d: round(float(h), 3) for d, h in zip(drugs, drop["harm"])}
        cf["harm_mean_model"] = round(float(base["harm"][0]), 3)

    solved = reward >= best - SOLVED_WITHIN
    if solved:
        out["class"] = "solved"
        why = []
        named = [m for m in minors if m["alteration"] in seen]
        for m in named:
            s = seen[m["alteration"]]
            where = "the biopsy" if s["cells"] == "tumour" else f"line {s['cells']}, grown from surviving cells"
            why.append(f"named the {m['alteration']} clone ({pct(m['share'])} of the tumour) by sequencing {where} (result read on day {s['day']:g})")
        if tested_t:
            why.append("cultured the submitted combination on tumour cells before submitting")
        if tested_h:
            why.append("checked it on the patient's healthy cells")
        add("solved", TITLES["solved"], (f"Reward {reward:.2f} against {best:.2f} for the best regimen known. " +
                                         ("It " + "; ".join(why) + "." if why else "")).strip())
        out["summary"] = mistakes[0]["detail"]
        return out

    if relapse:
        k = int(np.argmax(tum_end))
        if k > 0:
            m = next(x for x in minors if x["index"] == k)
            out["class"] = "missed_clone"
            out["relapse_clone"] = m["alteration"]
            if m["alteration"] in seen:
                s = seen[m["alteration"]]
                detail = (f"The tumour grew back from the {m['alteration']} clone ({pct(m['share'])} of the tumour). The agent had it "
                          f"named on day {s['day']:g} (sequencing of {'the biopsy' if s['cells'] == 'tumour' else 'line ' + s['cells']}) "
                          f"and still submitted {show(regimen)}, which does not kill it.")
                add("clone_not_covered", f"Knew about the {m['alteration']} clone, did not cover it", detail)
            else:
                if m["share"] >= L.SEQ_LOD and not facts["biopsy_sequenced"]:
                    how = "The biopsy was never sequenced; one sequencing run would have shown it"
                elif facts["lines_derived"] == 0:
                    how = ("It is below the 2 % sequencing limit, so the biopsy read shows nothing. The agent never grew out "
                           "the cells that survived its drug wells, the only way to name it")
                elif facts["lines_sequenced"] == 0:
                    how = f"The agent derived {facts['lines_derived']} line(s) from survivors but never sequenced one"
                else:
                    how = (f"The agent sequenced {facts['lines_sequenced']} line(s) of survivors, none of them from a well where "
                           f"this clone had been selected")
                detail = (f"The tumour grew back from the {m['alteration']} clone ({pct(m['share'])} of the tumour on day 0). {how}. "
                          f"Submitted: {show(regimen)}.")
                add("missed_clone", f"Never found the {m['alteration']} clone", detail)
        else:
            out["class"] = "tumour_survived"
            ext = lattice(pt, extended(regimen)) if treated < M.TREAT_DAYS else None
            if ext:
                cf["held_to_day_84"] = {k2: round(v, 3) for k2, v in ext.items()}
            if ext and ext["reward"] >= reward + 0.2:
                add("stopped_too_early", f"Stopped treatment on day {treated}",
                    f"The main clone was still there when dosing stopped on day {treated} of 84 and grew back. The same doses held "
                    f"to day 84 would have scored {ext['reward']:.2f}. Longest tumour culture the agent watched: "
                    f"{facts['longest_tumour_culture_days']:g} days.")
            elif best_var["reward"] >= reward + 0.2:
                add("wrong_dose", "Right drugs, wrong dose",
                    f"The main clone survived {show(regimen)}. {names[j].capitalize()} would have scored {best_var['reward']:.2f}.")
            else:
                add("wrong_drugs", "These drugs do not clear the main clone",
                    f"The main clone survived {show(regimen)}; no rescaling of these drugs within the allowed doses cures "
                    f"(best: {best_var['reward']:.2f}). The best regimen known uses {show(bk.get('regimen'))}.")
            if not tested_t:
                mistakes[-1]["detail"] += " The submitted combination was never cultured on tumour cells."
    else:
        out["class"] = "too_toxic" if reward < 0.5 * best else "suboptimal"

    # ---- harm, for every unsolved run that hurt the patient
    if harm >= 0.25 and (not relapse or harm >= 0.4):
        worst = None
        if cf.get("harm_without_drug"):
            d, h = min(cf["harm_without_drug"].items(), key=lambda kv: kv[1])
            if cf["harm_mean_model"] - h >= 0.1:
                worst = f" Most of it comes from {d}: without it the loss would be {pct(h)}."
        check = ("The submitted combination was never put on the patient's healthy cells in a well"
                 if not tested_h else "The agent had this combination on healthy cells in a well and submitted it anyway")
        if facts["healthy_wells"] == 0:
            check = "The agent never ran a single well of the patient's healthy cells"
        better = ""
        if best_var["reward"] >= reward + 0.1:
            better = f" {names[j].capitalize()} would have scored {best_var['reward']:.2f}."
        add("avoidable_harm" if out.get("class") == "suboptimal" else "too_toxic", f"Normal tissue fell to {pct(max(0.0, 1 - harm))}",
            f"{'The tumour was gone, but the' if not relapse else 'The'} regimen ({show(regimen)}) cost {pct(harm)} of the normal "
            f"tissue at its lowest point.{worst or ''} {check}.{better}")
    if out["class"] == "suboptimal" and not mistakes:
        add("suboptimal", TITLES["suboptimal"],
            f"Cured with {show(regimen)} for a reward of {reward:.2f}; the best regimen known scores {best:.2f} "
            f"({show(bk.get('regimen'))}).")
    out["summary"] = mistakes[0]["detail"] if mistakes else ""
    return out


# ----------------------------------------------------------------------------- job directories
EXCLUDED: list[dict] = []     # runs that say nothing about the agent: the harness or the model API failed


def api_error(trial_dir: Path) -> str | None:
    """The model API's error when the agent's session ended on one (OpenCode's last event is an error)."""
    log = trial_dir / "agent" / "opencode.txt"
    if not log.is_file():
        return None
    last = None
    with open(log, errors="replace") as f:
        for line in f:
            if line.startswith("{"):
                try:
                    last = json.loads(line)
                except ValueError:
                    pass
    if last and last.get("type") == "error":
        return json.dumps(last.get("error"))[:240]
    return None


def trials(label: str, job_dir: Path):
    for result in sorted(job_dir.glob("*/result.json")):
        r = json.loads(result.read_text())
        task = (r.get("task_name") or result.parent.name).split("/")[-1]
        row = {"label": label, "job": job_dir.name, "trial": result.parent.name, "task": task}
        try:
            data = json.loads((result.parent / "verifier" / "episode.json").read_text())
        except (OSError, ValueError):
            exc = (r.get("exception_info") or {}).get("exception_type") or "no verifier result"
            EXCLUDED.append({**row, "reason": f"the trial never reached the verifier ({exc})"})
            continue
        if data["episode"].get("closed_by") != "submit":
            err = api_error(result.parent)
            if err:
                EXCLUDED.append({**row, "reason": f"the model API failed before a regimen was submitted: {err}"})
                continue
        row.update(diagnose(data["episode"], data["score"]))
        if row["facts"].get("n_calls", 0) == 0 and label not in ("Do nothing",):
            row["idle"] = True
        yield row


def bootstrap(by_task: dict[str, list[float]], n: int = 20000, seed: int = 0) -> tuple[float, float, float]:
    """Mean over tasks of the per-task mean, with a 95 % interval from resampling tasks and runs."""
    rng = np.random.default_rng(seed)
    tasks = sorted(by_task)
    if not tasks:
        return 0.0, 0.0, 0.0
    point = float(np.mean([np.mean(by_task[t]) for t in tasks]))
    stats = np.empty(n)
    arrs = [np.asarray(by_task[t], dtype=float) for t in tasks]
    for b in range(n):
        pick = rng.integers(0, len(arrs), len(arrs))
        stats[b] = np.mean([rng.choice(arrs[i], len(arrs[i])).mean() for i in pick])
    lo, hi = np.percentile(stats, [2.5, 97.5])
    return point, float(lo), float(hi)


def summarise(rows: list[dict], reference: str) -> dict:
    labels = list(dict.fromkeys(r["label"] for r in rows))
    ref = defaultdict(list)
    for r in rows:
        if r["label"] == reference:
            ref[r["task"]].append(r["reward"])
    out = {"reference": reference, "labels": {}}
    for lab in labels:
        mine = [r for r in rows if r["label"] == lab]
        by_task = defaultdict(list)
        for r in mine:
            by_task[r["task"]].append(r["reward"])
        mean, lo, hi = bootstrap(by_task)
        s = {"runs": len(mine), "tasks": len(by_task), "reward_mean": round(mean, 3), "reward_ci95": [round(lo, 3), round(hi, 3)],
             "cured_runs": sum(1 for r in mine if r.get("cured_share", 0) >= 0.5),
             "solved_runs": sum(1 for r in mine if r.get("class") == "solved"),
             "classes": dict(Counter(r.get("class", "error") for r in mine).most_common()),
             "mistakes": dict(Counter(m["key"] for r in mine for m in r["mistakes"] if m["key"] != "solved").most_common()),
             "wall_minutes_median": round(float(np.median([r["facts"].get("wall_minutes", 0) for r in mine])), 1) if mine else 0,
             "tool_calls_median": float(np.median([r["facts"].get("n_calls", 0) for r in mine])) if mine else 0}
        if ref and lab != reference:
            gap = {t: [float(np.mean(ref[t])) - v for v in vs] for t, vs in by_task.items() if t in ref}
            g, glo, ghi = bootstrap(gap, seed=1)
            s["gap_to_reference"] = round(g, 3)
            s["gap_ci95"] = [round(glo, 3), round(ghi, 3)]
            s["runs_at_reference_level"] = sum(1 for t, vs in by_task.items() if t in ref for v in vs if v >= np.mean(ref[t]) - 0.05)
        out["labels"][lab] = s
    return out


def markdown(summary: dict) -> str:
    lines = ["| run | tasks x runs | reward (95 % CI) | tumour cleared | gap to reference (95 % CI) | tags (a run can carry several) |",
             "|---|---|---|---|---|---|"]
    for lab, s in summary["labels"].items():
        gap = "-" if "gap_to_reference" not in s else f"{s['gap_to_reference']:.2f} ({s['gap_ci95'][0]:.2f} to {s['gap_ci95'][1]:.2f})"
        wrong = ", ".join(f"{k} {v}" for k, v in s["mistakes"].items()) or "-"
        per = s["runs"] // max(1, s["tasks"])
        lines.append(f"| {lab} | {s['tasks']} x {per} | {s['reward_mean']:.2f} ({s['reward_ci95'][0]:.2f} to {s['reward_ci95'][1]:.2f}) "
                     f"| {s['cured_runs']}/{s['runs']} | {gap} | {wrong} |")
    return "\n".join(lines) + "\n"


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("columns", nargs="+", help="LABEL=JOB_DIR")
    ap.add_argument("--reference", default="Reference scientist")
    ap.add_argument("--out", default=str(ROOT / "results" / "run_analysis.jsonl"))
    ap.add_argument("--summary", default=str(ROOT / "results" / "summary.json"))
    ap.add_argument("--md", default=str(ROOT / "results" / "summary.md"))
    a = ap.parse_args()
    rows = []
    for col in a.columns:
        label, _, path = col.partition("=")
        job = Path(path)
        if not job.is_dir():
            print(f"skipping {label}: no such job directory {job}", file=sys.stderr)
            continue
        rows += list(trials(label, job))
    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    with open(a.out, "w") as f:
        for r in rows:
            f.write(json.dumps(r) + "\n")
    summary = summarise(rows, a.reference)
    summary["excluded"] = EXCLUDED
    for e in EXCLUDED:
        print(f"excluded {e['label']} {e['task']} ({e['trial']}): {e['reason']}", file=sys.stderr)
    for path, text in ((a.summary, json.dumps(summary, indent=1) + "\n"), (a.md, markdown(summary))):
        Path(path).parent.mkdir(parents=True, exist_ok=True)
        Path(path).write_text(text)
    print(markdown(summary))
    for r in rows:
        if r.get("class") not in ("solved",) and r["label"] != "Do nothing":
            print(f"{r['label']:22s} {r['task']:8s} {r['reward']:.2f}  {r.get('class', '?'):16s} {r['summary'][:230]}")
