#!/usr/bin/env python3
"""Collect finished Remission runs into compact data for the replay viewer.

    python3 scripts/build_viewer.py "Reference scientist=/path/to/job" "Claude Fable 5.2=/path/to/job" ... \
        [--out viewer/data] [--mistakes results/run_analysis.jsonl]

Each JOB_DIR is a Harbor job directory: <job_dir>/<trial>/verifier/episode.json and
<job_dir>/<trial>/result.json. Output, loaded by viewer/index.html through <script> tags
(so the page also works from file://):

    <out>/index.js              patients -> runs (label, reward, healed, harm, cured, file)
    <out>/<patient>.<run>.js    one run: replay frames, curves, lab, diagnosis   (<= 250 KB)
    <out>/<patient>.best.js     the best-known regimen for that patient, when the episode has one

Only the standard library is used.
"""
from __future__ import annotations

import argparse
import base64
import json
import re
import sys
import zlib
from pathlib import Path

DRUGS = ("RAFi", "MEKi", "EGFRi", "METi", "KRASi", "PI3Ki", "BCLi", "MCLi")
TREAT_DAYS = 84
CLASS_TAG = {"solved": "healed", "suboptimal": "healed, extra harm", "too_toxic": "normal tissue lost"}   # by analysis class
LAB_DAYS = 28
WELL_BUDGET = 384
MAX_FRAMES = 72
NBSP = " "


# ----------------------------------------------------------------------------- small helpers
def slug(text: str) -> str:
    return re.sub(r"[^a-z0-9]+", "-", text.lower()).strip("-") or "run"


def sig(x: float, n: int = 3) -> float:
    return float(f"{float(x):.{n}g}")


def pct(x: float, digits: int = 0) -> str:
    v = round(100.0 * float(x), 4)
    if digits == 0 and 0.05 <= v < 10:
        digits = 1
    return f"{v:.{digits}f}{NBSP}%"


def num(x: float) -> str:
    x = float(x)
    if x >= 10 or x == int(x):
        return f"{x:,.0f}".replace(",", NBSP)
    return f"{x:.1f}"


def plural(n, word: str) -> str:
    n = int(n) if float(n) == int(n) else n
    return f"{n} {word}" if n == 1 else f"{n} {word}s"


def clone_label(clone: dict, dominant: dict | None, index: int) -> str:
    """Short name of a clone by what sets it apart from the dominant clone."""
    if index == 0 or dominant is None:
        return "dominant clone"
    own = [a for a in clone["alterations"] if a not in dominant["alterations"]] or clone["alterations"]
    a = own[0]
    m = re.search(r"\(([^)]+)\)", a)
    gene = m.group(1) if m else a.split()[0]
    if "amplification" in a:
        return f"{gene}-amplified clone"
    if a.endswith(" loss"):
        return f"{gene}-loss clone"
    if m:                                   # e.g. "MAP2K1 (MEK1) K57N"
        return f"{gene} {a.split()[-1]} clone"
    return f"{a} clone"


def own_alterations(clone: dict, dominant: dict) -> list[str]:
    return [a for a in clone["alterations"] if a not in dominant["alterations"]]


def drugs_text(doses: dict, unit: str = "") -> str:
    return " + ".join(f"{d} {v:g}{unit}" for d, v in doses.items()) or "no drug"


# ----------------------------------------------------------------------------- replay packing
def unpack_frames(frames: list[dict]) -> tuple[list[int], list[bytes]]:
    return [int(f["day"]) for f in frames], [base64.b64decode(f["grid"]) for f in frames]


def pick_frames(days: list[int], grids: list[bytes], keep: int = MAX_FRAMES, must: tuple[int, ...] = ()) -> list[int]:
    """Indices of at most `keep` frames, spent where the tissue changes most."""
    n = len(grids)
    if n <= keep:
        return list(range(n))
    change = [0.0]
    for i in range(1, n):
        diff = sum(1 for a, b in zip(grids[i - 1], grids[i]) if a != b)
        change.append(diff + 12.0 * (days[i] - days[i - 1]))      # a floor so quiet stretches still get frames
    cum, total = [], 0.0
    for c in change:
        total += c
        cum.append(total)
    chosen = {0, n - 1}
    for d in must:                                                # e.g. the last treatment day
        chosen.add(min(range(n), key=lambda i: abs(days[i] - d)))
    targets = [total * k / (keep - 1) for k in range(keep)]
    j = 0
    for t in targets:
        while j < n - 1 and cum[j] < t:
            j += 1
        chosen.add(j)
    out = sorted(chosen)
    while len(out) > keep:                                        # drop the frame whose removal leaves the smallest gap
        k = min(range(1, len(out) - 1), key=lambda i: cum[out[i + 1]] - cum[out[i - 1]])
        out.pop(k)
    return out


def pack_replay(rep: dict, treat_end: int) -> dict:
    days, grids = unpack_frames(rep["frames"])
    idx = pick_frames(days, grids, must=(treat_end,))
    blob = zlib.compress(b"".join(grids[i] for i in idx), 9)
    tumour = rep["tumour"]
    n_clones = len(tumour[0]) if tumour else 0
    out = {
        "L": int(rep.get("lattice", 64)),
        "days": [days[i] for i in idx],
        "grids": base64.b64encode(blob).decode(),
        "healthy": [int(round(v)) for v in rep["healthy"]],
        "tumour": [[round(float(row[c]), 1) if row[c] < 100 else int(round(row[c])) for row in tumour] for c in range(n_clones)],
        "normal": [round(float(v), 3) for v in rep["normal"]],
    }
    if "conc" in rep:
        conc = {}
        for j, d in enumerate(DRUGS):
            col = [float(r[j]) for r in rep["conc"]]
            last = max((i for i, v in enumerate(col) if v > 0), default=-1)
            if last >= 0:
                conc[d] = [sig(v) for v in col[:last + 1]]
        out["conc"] = conc
    if "reward_reps" in rep:
        out["reward_reps"] = [round(float(v), 3) for v in rep["reward_reps"]]
    return out


# ----------------------------------------------------------------------------- lab
def pack_experiments(experiments: list[dict]) -> list[dict]:
    out = []
    for e in experiments:
        x = {"id": e["id"], "kind": e["kind"], "cells": e["cells"], "start": e["start_day"], "end": e["end_day"]}
        if "wells" in e:
            x["wells"] = [{d: sig(v, 4) for d, v in w.items()} for w in e["wells"]]
        if e["kind"] == "viability":
            x["read_days"] = e["read_days"]
            x["counts"] = e["counts"]
        elif e["kind"] == "phospho":
            x["pERK"] = [round(v, 2) for v in e["pERK"]]
            x["pAKT"] = [round(v, 2) for v in e["pAKT"]]
        elif e["kind"] == "sequencing":
            x["alterations"] = e.get("alterations", [])
            x["lod"] = e.get("detection_limit", 0.02)
            x["n_wells"] = 16
        out.append(x)
    return out


def call_summary(c: dict) -> str:
    a = c.get("args") or {}
    if not isinstance(a, dict) or a.get("_truncated"):
        head = a.get("head", "") if isinstance(a, dict) else ""
        kind = re.search(r'"kind":\s*"(\w+)"', head)
        cells = re.search(r'"cells":\s*"(\w+)"', head)
        days = re.search(r'"days":\s*([\d.]+)', head)
        bits = [m.group(1) for m in (kind, cells) if m]
        if days:
            bits.append(f"{float(days.group(1)):g} d")
        return " · ".join(bits)
    t = c["tool"]
    if t == "start_experiment":
        wells = a.get("wells") or []
        drugs = sorted({d for w in wells if isinstance(w, dict) for d in w})
        bits = [str(a.get("kind", "")), str(a.get("cells", ""))]
        if wells:
            bits.append(plural(len(wells), "well"))
        if a.get("days") is not None:
            bits.append(f"{float(a['days']):g} d")
        if drugs:
            bits.append(", ".join(drugs))
        return " · ".join(b for b in bits if b)
    if t == "derive_line":
        return f"{a.get('experiment_id', '?')} well {a.get('well', '?')}"
    if t == "advance_time":
        return f"+{a.get('days', '?')} d"
    if t == "get_results":
        return str(a.get("experiment_id", ""))
    if t == "submit_regimen":
        ph = (a.get("regimen") or {}).get("phases") or []
        return "; ".join(f"{p.get('days')} d: {drugs_text(p.get('doses') or {})}" for p in ph if isinstance(p, dict))[:160]
    return ""


def pack_calls(calls: list[dict], experiments: list[dict]) -> list[list]:
    by_start: dict[tuple, list[dict]] = {}
    out = []
    n_exp = 0
    for c in calls:
        text = call_summary(c)
        if c["tool"] == "start_experiment" and c.get("ok"):
            e = experiments[n_exp] if n_exp < len(experiments) else None
            n_exp += 1
            if e is not None:
                wells = len(e.get("wells", [])) or 16
                drugs = sorted({d for w in e.get("wells", []) for d in w})
                text = " · ".join([e["id"], e["kind"], e["cells"], plural(wells, "well"),
                                   f"{e['end_day'] - e['start_day']:g} d"] + ([", ".join(drugs)] if drugs else []))
        out.append([c["tool"], c.get("day_before", c.get("day", 0)), round(float(c.get("t_wall", 0)), 1),
                    1 if c.get("ok") else 0, text])
    del by_start
    return out


# ----------------------------------------------------------------------------- diagnosis
def diagnose(score: dict, ep: dict, best: dict | None) -> dict:
    """Facts computed from the episode, and from them a verdict and up to three callouts."""
    truth = ep["truth"]
    clones = truth["clones"]
    dom = clones[0]
    labels = [clone_label(c, dom, i) for i, c in enumerate(clones)]
    rep = ep["replay"]
    T = rep["tumour"]
    H = rep["healthy"]
    normal = rep["normal"]
    n_days = len(T) - 1
    regimen = ep.get("regimen")
    submitted = regimen is not None
    treat_end = min(TREAT_DAYS, sum(int(p["days"]) for p in regimen["phases"])) if submitted else 0
    last_dose = 0
    if submitted:
        d = 0
        for p in regimen["phases"]:
            d += int(p["days"])
            if p.get("doses"):
                last_dose = min(d, TREAT_DAYS)
    reps = len(rep.get("reward_reps") or []) or 16
    cured = float(score.get("cured", 0.0))
    relapsed = round((1.0 - cured) * reps)

    per_clone = []
    for i, c in enumerate(clones):
        col = [float(r[i]) for r in T]
        nadir_day = min(range(len(col)), key=lambda k: col[k])
        gone = next((k for k in range(len(col)) if all(v == 0 for v in col[k:])), None)
        per_clone.append({"label": labels[i], "share": c["share"], "start": col[0], "at_stop": col[min(last_dose, n_days)],
                          "nadir": col[nadir_day], "nadir_day": nadir_day, "end": col[-1], "gone_day": gone})
    total = [sum(r) for r in T]
    gone_day = next((k for k in range(len(total)) if all(v == 0 for v in total[k:])), None)
    takeover = next((k for k in range(len(total)) if total[k] > H[k]), None)
    healthy_gone = next((k for k in range(len(H)) if all(v < 1 for v in H[k:])), None)
    regrew = [i for i, pc in enumerate(per_clone) if pc["end"] > 0.5]
    main = max(regrew, key=lambda i: per_clone[i]["end"]) if regrew else None
    harm = float(score.get("harm", 0.0))
    harm_day = min(range(len(normal)), key=lambda k: normal[k])
    best_harm = float(best["harm"]) if best else None
    best_reward = float(best["reward"]) if best else None

    # what happened in the lab
    exps = ep.get("experiments") or []
    lines = ep.get("lines") or {}
    seq = [e for e in exps if e["kind"] == "sequencing"]
    seq_biopsy = [e for e in seq if e["cells"] == "tumour"]
    seq_lines = [e for e in seq if lines.get(e["cells"], {}).get("kind") == "derived"]
    derived = {k: v for k, v in lines.items() if v.get("kind") == "derived"}
    found: dict[int, str] = {}                # minor clone index -> experiment that reported its alteration
    for i, c in enumerate(clones[1:], start=1):
        own = own_alterations(c, dom)
        for e in seq:
            seen = {a["alteration"] for a in e.get("alterations", [])}
            if own and all(a in seen for a in own):
                found[i] = e["id"]
                break
    via = [e for e in exps if e["kind"] == "viability"]
    via_tumour = [e for e in via if lines.get(e["cells"], {}).get("kind") != "healthy"]
    via_healthy = [e for e in via if lines.get(e["cells"], {}).get("kind") == "healthy"]
    longest = max((e["end_day"] - e["start_day"] for e in via_tumour), default=0.0)
    long_ids = [e["id"] for e in via_tumour if e["end_day"] - e["start_day"] > 5]
    combos = [frozenset(p.get("doses") or {}) for p in (regimen["phases"] if submitted else []) if p.get("doses")]
    final_combo = max(combos, key=len) if combos else frozenset()

    def tested(es):
        return [e["id"] for e in es if any(frozenset(w) == final_combo for w in e.get("wells", []))] if final_combo else []

    tested_tumour, tested_healthy = tested(via_tumour), tested(via_healthy)
    wells_used = int(score.get("wells_used", 0))
    days_used = float(score.get("lab_days_used", 0))
    wells_left, days_left = WELL_BUDGET - wells_used, LAB_DAYS - days_used
    n_exp = len(exps)

    facts = {
        "submitted": submitted, "treat_end": treat_end, "last_dose_day": last_dose, "clones": per_clone,
        "tumour_gone_day": gone_day, "takeover_day": takeover, "healthy_gone_day": healthy_gone,
        "regrew": [labels[i] for i in regrew], "harm_day": harm_day,
        "sequenced_biopsy": [e["id"] for e in seq_biopsy], "sequenced_lines": [e["id"] for e in seq_lines],
        "derived_lines": list(derived), "found_clones": {labels[i]: e for i, e in found.items()},
        "longest_tumour_viability_days": longest, "viability_over_5_days": long_ids,
        "final_combo": sorted(final_combo), "final_combo_tested_on_tumour": tested_tumour,
        "final_combo_tested_on_healthy": tested_healthy, "wells_left": wells_left, "lab_days_left": days_left,
    }

    toxic = harm > 0.3 and (best_harm is None or harm > best_harm + 0.1)
    reward = float(score.get("reward", 0.0))
    success = submitted and cured >= 0.99 and not toxic
    callouts: list[dict] = []

    def add(key, title, detail, tone="bad", day=None, exp=None):
        c = {"key": key, "title": title, "detail": detail, "tone": tone}
        if day is not None:
            c["day"] = int(day)
        if exp:
            c["exp"] = list(exp)
        callouts.append(c)

    # ---- verdict: one line under the lattice
    if not submitted:
        verdict = "Nothing submitted: the untreated tumour took the tissue over"
    elif main is not None and cured < 0.5:
        verdict = f"Tumour regrew from the {labels[main]}"
    elif main is not None:
        verdict = f"Tumour gone in {reps - relapsed} of {reps} lesions; regrew from the {labels[main]} in {relapsed}"
    elif toxic:
        verdict = f"Tumour gone, but normal tissue fell to {pct(1 - harm)}"
    else:
        verdict = f"Tumour gone by day {gone_day}" if gone_day is not None else "Tumour gone"
    if main is not None and toxic and submitted:
        verdict += f"; normal tissue fell to {pct(1 - harm)}"

    lab_none = n_exp == 0
    combo_text = " + ".join(sorted(final_combo, key=DRUGS.index))
    # ---- failure callouts, each tied to a decision; the first three are shown
    if not submitted:
        used = "No experiment was run." if lab_none else f"{plural(n_exp, 'experiment')} and {wells_used} wells were used, but the lab closed with no regimen."
        add("not_submitted", "No regimen submitted",
            f"{used} No regimen means an untreated patient: tumour cells outnumber healthy cells from day {takeover}"
            + (f" and the last healthy cells in the lesion are gone by day {healthy_gone}." if healthy_gone else "."),
            day=takeover)
    if submitted and main is not None:
        pc = per_clone[main]
        share = f"{pct(pc['share'])} of the tumour on day 0" if main else "The bulk of the tumour"
        if pc["at_stop"] >= pc["start"]:
            how = (f"It was never pushed back: {num(pc['start'])} cells on day 0, {num(pc['at_stop'])} on day {last_dose}, "
                   f"the last day of treatment.")
            title = f"The regimen did not kill the {labels[main]}"
        elif pc["nadir"] < 1.0 and relapsed < reps:
            how = (f"Pushed below one cell per lesion (day {pc['nadir_day']}) but it survived in {relapsed} of {reps} lesions "
                   f"and regrew after the last dose on day {last_dose}.")
            title = f"The {labels[main]} was held, not killed"
        else:
            how = (f"Down to {num(pc['nadir'])} cells on day {pc['nadir_day']}, never extinct; it regrew after the last dose "
                   f"on day {last_dose}.")
            title = f"The {labels[main]} was held, not killed"
        tail = f" Tumour cells outnumber healthy cells from day {takeover}." if takeover is not None else ""
        add("regrowth", title, f"{share}. {how}{tail}", day=takeover if takeover is not None else pc["nadir_day"])
    if submitted and toxic:
        why = ""
        if not lab_none and final_combo and not tested_healthy:
            why = f" {combo_text} was never run together on the patient's healthy cells."
        add("toxic", f"Too toxic: normal tissue fell to {pct(1 - harm)}",
            f"Harm {harm:.2f} on day {harm_day}"
            + (f"; a regimen that cures this patient costs {best_harm:.2f}." if best_harm is not None else ".")
            + why + (f" Reward = healed {score['healed']:.2f} − harm {harm:.2f}." if main is None else ""),
            day=harm_day, exp=[e["id"] for e in via_healthy])
    if submitted and main is not None and main != 0:
        share_v = per_clone[main]["share"]
        untested = "" if (lab_none or tested_tumour or not final_combo) else f" The submitted combination ({combo_text}) was never run together on tumour cells."
        if main in found:
            add("found_not_killed", f"Found the {labels[main]}, then did not test the answer",
                f"Sequencing {found[main]} reported it." + (untested or " The submitted regimen still left it alive."),
                exp=[found[main]])
        elif lab_none:
            add("missed_clone", f"Never looked for the {labels[main]}",
                f"0 of {WELL_BUDGET} wells and 0 of {LAB_DAYS} lab days used: the regimen was chosen from the chart alone, "
                f"and the chart's bulk panel does not report clones under 10{NBSP}%.")
        elif share_v >= 0.02 and not seq_biopsy:
            add("missed_clone", f"Never found the {labels[main]}",
                f"The biopsy was never deep-sequenced (16 wells, 2 days), which reports alterations in at least 2{NBSP}% of "
                f"cells; this clone is {pct(share_v)}." + ("" if derived else " No line was derived from surviving wells either."))
        elif not derived:
            add("missed_clone", f"Never found the {labels[main]}",
                (f"At {pct(share_v)} it is under the 2{NBSP}% limit of biopsy sequencing. " if share_v < 0.02 else "")
                + "It only shows as wells that regrow under drug; no line was derived from such a well and sequenced.",
                exp=[e["id"] for e in seq_biopsy])
        else:
            add("missed_clone", f"Never found the {labels[main]}",
                f"{plural(len(derived), 'line')} derived ({', '.join(derived)}), but none was shown by sequencing to carry "
                f"{' + '.join(own_alterations(clones[main], dom))}.", exp=[e["id"] for e in seq_lines])
    if submitted and not success and not lab_none:
        if final_combo and not tested_tumour and not any(c["key"] == "found_not_killed" for c in callouts):
            add("untested", "Submitted a combination it never tested",
                f"{combo_text} was not run together in any well on the patient's tumour cells.")
        elif longest <= 5 and main is not None:
            add("short_assays", "No experiment long enough to see regrowth",
                f"The longest viability run on tumour cells was {longest:g} days; a resistant clone of a few per cent needs "
                f"more than 5 days to show as regrowth.")
    if not success and not lab_none and (wells_left >= WELL_BUDGET / 2 or days_left >= 7):
        add("unused", f"Left {wells_left} of {WELL_BUDGET} wells and {days_left:g} of {LAB_DAYS} lab days unused",
            "The budget allowed more experiments before committing to a regimen.")
    if submitted and not success and lab_none and not any(c["key"] == "missed_clone" for c in callouts):
        add("blind", "Chosen without a single experiment", f"0 of {WELL_BUDGET} wells and 0 of {LAB_DAYS} lab days used.")

    # ---- what made it work
    if success:
        for i, eid in found.items():
            e = next(x for x in exps if x["id"] == eid)
            own = " + ".join(own_alterations(clones[i], dom))
            if lines.get(e["cells"], {}).get("kind") == "derived":
                add("found_clone", f"Derived a line from surviving wells and sequenced it: found the {labels[i]}",
                    f"{e['cells']} = {lines[e['cells']]['note'].split(' (')[0]}; sequencing {eid} reported {own}. "
                    f"The clone is {pct(clones[i]['share'])} of the tumour, under the limit of biopsy sequencing.",
                    tone="good", exp=[eid])
            else:
                frac = next((a["cell_fraction"] for a in e["alterations"] if a["alteration"] in own_alterations(clones[i], dom)), None)
                add("found_clone", f"Sequenced the biopsy: found the {labels[i]}",
                    f"{eid} reported {own} in {pct(frac) if frac is not None else 'a few per cent'} of cells; the chart's bulk "
                    f"panel does not show it.", tone="good", exp=[eid])
        if tested_tumour:
            eid = tested_tumour[-1]
            e = next(x for x in exps if x["id"] == eid)
            add("tested", "Tested the final combination before submitting",
                f"{' + '.join(sorted(final_combo, key=DRUGS.index))} on tumour cells for {e['end_day'] - e['start_day']:g} days ({', '.join(tested_tumour)})"
                + (f" and on healthy cells ({', '.join(tested_healthy)})." if tested_healthy else "."),
                tone="good", exp=tested_tumour + tested_healthy)
        elif lab_none:
            add("blind", "Chosen without a single experiment", f"0 of {WELL_BUDGET} wells used; the regimen happened to fit this patient.",
                tone="note")
        if gone_day is not None:
            add("cleared", f"Every tumour cell gone by day {gone_day}",
                f"All {reps} lesions are tumour-free on day {n_days}; healthy cells refilled {pct(score['healed'])} of the "
                f"lesion. Normal tissue dipped to {pct(1 - harm)} (harm {harm:.2f}).", tone="good", day=gone_day)

    summary = verdict + "."
    if not submitted:
        tag = "nothing submitted"
    elif main is not None and cured < 0.5:
        tag = "regrew: " + labels[main].replace(" clone", "") + (", toxic" if toxic else "")
    elif main is not None:
        tag = f"cured {reps - relapsed} of {reps}" + (", toxic" if toxic else "")
    elif toxic:
        tag = "tumour gone, too toxic"
    else:
        tag = "tumour gone"
    return {"verdict": verdict, "tag": tag, "success": bool(success), "callouts": callouts[:3], "summary": summary,
            "facts": facts, "source": "derived"}


# ----------------------------------------------------------------------------- runs
def run_kind(label: str, agent_cfg: dict) -> str:
    name = (agent_cfg.get("name") or "").lower()
    if agent_cfg.get("model_name"):
        return "agent"
    if re.search(r"reference|oracle scientist", label, re.I):
        return "reference"
    if name == "oracle" and not (agent_cfg.get("env") or {}):
        return "reference"
    return "baseline"


def load_mistakes(path: Path | None) -> list[dict]:
    if not path:
        return []
    if not path.exists():
        print(f"warning: {path} not found; diagnosis is derived from the episodes", file=sys.stderr)
        return []
    out = []
    for line in path.read_text().splitlines():
        line = line.strip()
        if line:
            try:
                out.append(json.loads(line))
            except ValueError:
                print(f"warning: skipped a bad line in {path}", file=sys.stderr)
    return out


def find_mistakes(rows: list[dict], job: str, trial: str, label: str, task: str) -> dict | None:
    def base(p):
        return Path(str(p or "")).name
    for r in rows:
        if base(r.get("job")) == job and base(r.get("trial")) == trial:
            return r
    hits = [r for r in rows if r.get("label") == label and str(r.get("task", "")).split("/")[-1] == task.split("/")[-1]
            and (not r.get("trial") or base(r.get("trial")) == trial)]
    return hits[0] if len(hits) == 1 else None


def write_js(path: Path, call: str, obj: dict) -> int:
    text = f"REMISSION_DATA.{call}(" + json.dumps(obj, separators=(",", ":"), ensure_ascii=False) + ");\n"
    path.write_text(text, encoding="utf-8")
    return len(text.encode("utf-8"))


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("runs", nargs="+", metavar="LABEL=JOB_DIR")
    ap.add_argument("--out", default=str(Path(__file__).resolve().parent.parent / "viewer" / "data"))
    ap.add_argument("--mistakes", default=None, help="JSONL with one object per run (job, trial, label, task, mistakes, summary)")
    ap.add_argument("--open", default=None, metavar="PATIENT:RUN_ID:DAY",
                    help="what the viewer opens with, e.g. test-07:gpt-6-1-sol:24 (default: the largest agent failure)")
    args = ap.parse_args()

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    mistakes = load_mistakes(Path(args.mistakes) if args.mistakes else None)

    patients: dict[str, dict] = {}
    labels: list[dict] = []
    written: list[str] = []
    for order, spec in enumerate(args.runs):
        if "=" not in spec:
            ap.error(f"expected LABEL=JOB_DIR, got {spec!r}")
        label, job_dir = spec.split("=", 1)
        label, job = label.strip(), Path(job_dir.strip())
        if not job.is_dir():
            print(f"warning: {job} is not a directory; skipped {label!r}", file=sys.stderr)
            continue
        kind_of_label = None
        seen: dict[str, int] = {}
        for trial in sorted(p for p in job.iterdir() if p.is_dir()):
            res_path, ep_path = trial / "result.json", trial / "verifier" / "episode.json"
            result = {}
            if res_path.exists():
                try:
                    result = json.loads(res_path.read_text())
                except ValueError:
                    result = {}
            task = result.get("task_name") or trial.name.split("__")[0]
            pid = task.split("/")[-1]
            agent_cfg = (result.get("config") or {}).get("agent") or {}
            kind = run_kind(label, agent_cfg)
            kind_of_label = kind_of_label or kind
            seen[pid] = seen.get(pid, 0) + 1
            rid = slug(label) + (f"-{seen[pid]}" if seen[pid] > 1 else "")
            pat = patients.setdefault(pid, {"id": pid, "task": task, "runs": []})
            entry = {"id": rid, "label": label + (f" · run {seen[pid]}" if seen[pid] > 1 else ""), "kind": kind,
                     "order": order, "job": job.name, "trial": trial.name}
            if agent_cfg.get("model_name"):
                entry["model"] = agent_cfg["model_name"]
            if not ep_path.exists():
                exc = (result.get("exception_info") or {}).get("exception_message")
                if not res_path.exists():
                    entry.update({"missing": True, "running": True, "reason": "still running"})
                    pat["runs"].append(entry)
                    print(f"  {label:<24} {pid:<8} still running (no result.json)")
                    continue
                entry.update({"missing": True, "reason": (exc or "the run produced no episode")[:200]})
                pat["runs"].append(entry)
                print(f"  {label:<24} {pid:<8} no episode: {entry['reason'][:70]}")
                continue
            data = json.loads(ep_path.read_text())
            score, ep = data["score"], data["episode"]
            truth = ep["truth"]
            dom = truth["clones"][0]
            best = ep.get("best_known")
            regimen = ep.get("regimen")
            treat_end = min(TREAT_DAYS, sum(int(p["days"]) for p in regimen["phases"])) if regimen else 0

            if "chart" not in pat:
                pat.update({
                    "chart": ep["chart"], "archetype": truth["archetype"],
                    "clones": [{"name": c["name"], "label": clone_label(c, dom, i), "alterations": c["alterations"],
                                "own": own_alterations(c, dom) if i else c["alterations"], "share": round(c["share"], 4)}
                               for i, c in enumerate(truth["clones"])],
                })
            if best and "best" not in pat:
                b_end = min(TREAT_DAYS, sum(int(p["days"]) for p in best["regimen"]["phases"]))
                bscore = {"reward": best["reward"], "healed": best["healed"], "harm": best["harm"]}
                tot = [sum(r) for r in best["tumour"]]
                gone = next((k for k in range(len(tot)) if all(v == 0 for v in tot[k:])), None)

                def gone_at(col):
                    return next((k for k in range(len(col)) if all(v == 0 for v in col[k:])), None)
                bfacts = {"tumour_gone_day": gone, "healthy_gone_day": gone_at([float(v) for v in best["healthy"]]),
                          "clones": [{"gone_day": gone_at([float(r[i]) for r in best["tumour"]])}
                                     for i in range(len(best["tumour"][0]))]}
                bobj = {"id": "best", "patient": pid, "label": "Best-known regimen", "kind": "best", "score": bscore,
                        "regimen": best["regimen"], "treat_end": b_end, "replay": pack_replay(best, b_end),
                        "diag": {"verdict": f"Tumour gone by day {gone}" if gone is not None else "Tumour not cleared",
                                 "success": gone is not None, "source": "derived", "facts": bfacts,
                                 "callouts": [{"key": "best", "tone": "note", "title": "Found with the hidden parameters",
                                               "detail": "A regimen known to cure this patient, searched with the true "
                                                         "parameter values the lab never shows. No experiments: it is the "
                                                         "yardstick, not a run."}]}}
                fname = f"{pid}.best.js"
                size = write_js(out / fname, "add", bobj)
                written.append(fname)
                pat["best"] = {"id": "best", "label": "Best-known regimen", "kind": "best", "file": fname, **bscore,
                               "verdict": bobj["diag"]["verdict"], "tag": "tumour gone" if gone is not None else "",
                               "success": gone is not None}
                print(f"  {'Best-known regimen':<24} {pid:<8} reward {best['reward']:.2f}  {size / 1024:6.0f} KB")

            diag = diagnose(score, ep, best)
            row = find_mistakes(mistakes, job.name, trial.name, label, task)
            if row and row.get("mistakes"):
                tone = "good" if diag["success"] else "bad"
                diag["callouts"] = [{"key": m.get("key", ""), "title": m.get("title", ""), "detail": m.get("detail", ""),
                                     "tone": m.get("tone", tone), **({"day": m["day"]} if "day" in m else {}),
                                     **({"exp": m["exp"]} if "exp" in m else {})} for m in row["mistakes"][:3]]
                diag["source"] = "analysis"
            if row and row.get("summary"):
                diag["summary"] = row["summary"]
            if row and row.get("class") in CLASS_TAG and not diag["tag"].startswith(("cured", "regrew")):
                diag["tag"] = CLASS_TAG[row["class"]]          # the same words as the README's tables

            lines = {k: {"kind": v["kind"], "ready": v["ready_day"], "note": v["note"], "frac": v.get("frac")}
                     for k, v in (ep.get("lines") or {}).items()}
            derives = [{"day": l["day"], "source": l.get("source"), "well": l.get("well"), "ok": bool(l.get("ok")),
                        "line": l.get("line")} for l in ep.get("lab_log") or [] if l.get("tool") == "derive_line"]
            obj = {
                "id": rid, "patient": pid, "label": entry["label"], "kind": kind, "job": job.name, "trial": trial.name,
                "model": entry.get("model"),
                "score": {k: (round(float(v), 4) if isinstance(v, (int, float)) else v) for k, v in score.items()},
                "closed_by": ep.get("closed_by"), "submitted_day": ep.get("submitted_day"),
                "n_calls": ep.get("n_calls"), "wall_seconds": ep.get("wall_seconds"),
                "regimen": regimen, "treat_end": treat_end, "replay": pack_replay(ep["replay"], treat_end),
                "experiments": pack_experiments(ep.get("experiments") or []), "lines": lines, "derives": derives,
                "calls": pack_calls(ep.get("calls") or [], ep.get("experiments") or []), "diag": diag,
            }
            fname = f"{pid}.{rid}.js"
            size = write_js(out / fname, "add", obj)
            written.append(fname)
            entry.update({"file": fname, "reward": score["reward"], "healed": score["healed"], "harm": score["harm"],
                          "cured": score["cured"], "submitted": bool(score.get("submitted")),
                          "wells_used": int(score.get("wells_used", 0)), "lab_days_used": score.get("lab_days_used", 0),
                          "n_experiments": int(score.get("n_experiments", 0)), "verdict": diag["verdict"],
                          "tag": diag["tag"], "success": diag["success"]})
            pat["runs"].append(entry)
            flag = "  OVER 250 KB" if size > 250 * 1024 else ""
            print(f"  {label:<24} {pid:<8} reward {score['reward']:.2f}  {size / 1024:6.0f} KB{flag}")
        if kind_of_label:
            labels.append({"label": label, "kind": kind_of_label, "order": order})

    patients = {k: v for k, v in sorted(patients.items()) if "chart" in v}
    if not patients:
        print("no finished runs found (no <trial>/verifier/episode.json)", file=sys.stderr)
        return 1
    index = {"patients": list(patients.values()), "labels": labels, "files": sorted(set(written)),
             "lab": {"days": LAB_DAYS, "wells": WELL_BUDGET, "treat_days": TREAT_DAYS}, "drugs": list(DRUGS)}
    if args.open:
        pid, rid, day = (args.open.split(":") + ["", ""])[:3]
        if any(p["id"] == pid and any(r.get("id") == rid for r in p["runs"]) for p in patients.values()):
            index["open"] = {"patient": pid, "left": rid, **({"day": int(day)} if day.isdigit() else {})}
        else:
            print(f"--open {args.open}: no such run, ignored", file=sys.stderr)
    old = out / "index.js"
    stale = []
    if old.exists():                                   # drop data files of an earlier build that this one did not write
        m = re.search(r'"files":\[(.*?)\]', old.read_text())
        if m:
            stale = [f for f in re.findall(r'"([\w.\-]+\.js)"', m.group(1)) if f not in written]
    for f in stale:
        p = out / f
        if p.is_file() and p.parent == out:
            p.unlink()
    write_js(old, "index", index)
    n_runs = sum(1 for p in patients.values() for r in p["runs"] if not r.get("missing"))
    print(f"{n_runs} runs, {len(patients)} patients -> {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
