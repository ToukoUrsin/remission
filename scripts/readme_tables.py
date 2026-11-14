#!/usr/bin/env python3
"""Write the result tables of the README from the analysed runs, so the README cannot drift from the data.

    scripts/readme_tables.py [--analysis results/run_analysis.jsonl] [--summary results/summary.json]
                             [--readme README.md] [--print]

Replaces what is between `<!-- results:start -->` / `<!-- results:end -->` (one row per policy) and between
`<!-- patients:start -->` / `<!-- patients:end -->` (one row per patient, every run's reward) and between
`<!-- failures:start -->` / `<!-- failures:end -->` (how each agent's lost runs tie to its lab work) in the README.
"""
from __future__ import annotations

import argparse
import json
import re
import statistics
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SHORTCUTS = ("Hedge", "Blind guess", "Guideline regimen", "Three-day screen", "Stop early", "Do nothing")
WHAT = {  # what each scripted policy does
    "Reference scientist": "scripted scientist, public tools only",
    "Hedge": "sequences the biopsy, then covers every clone that could hide",
    "Blind guess": "best fixed regimen for the tumour type, no experiment",
    "Guideline regimen": "standard of care for the tumour type, no experiment",
    "Three-day screen": "one short drug screen, best-looking pair",
    "Stop early": "reads the brief and quits",
    "Do nothing": "no tool call",
}
CLASS = {  # how a run ended, in the README's words
    "solved": "healed",
    "suboptimal": "healed with avoidable harm",
    "too_toxic": "tumour gone, normal tissue destroyed",
    "missed_clone": "relapse from a hidden clone",
    "tumour_survived": "main tumour never cleared",
    "no_regimen": "no regimen",
}


def experience_arm(label: str) -> bool:
    """The runs with the model's own notes: reported in the README's experience section, not in the main tables."""
    return " with its " in label
ALTERATION = {"BCL2L1 (BCL-xL) amplification": "BCL-xL amplified", "MCL1 amplification": "MCL-1 amplified",
              "MET amplification": "MET amplified"}
TUMOUR = {"BRAF_MEL": "melanoma, BRAF V600E", "EGFR": "lung, EGFR L858R", "KRAS": "colorectal, KRAS G12C",
          "BRAF_CRC": "colorectal, BRAF V600E"}


def results(summary: dict) -> str:
    out = ["| | Runs | Reward (95 % CI) | Below the reference by (95 % CI) | Runs that healed the patient | How the other runs ended |",
           "|---|---|---|---|---|---|"]
    for lab, s in summary["labels"].items():
        if experience_arm(lab):
            continue
        lo, hi = s["reward_ci95"]
        name = f"**{lab}**" if lab not in SHORTCUTS else lab
        if lab in WHAT:
            name += f"<br><sub>{WHAT[lab]}</sub>"
        gap = "–"
        if "gap_to_reference" in s:
            glo, ghi = s["gap_ci95"]
            gap = f"{s['gap_to_reference']:.2f} ({glo:.2f} to {ghi:.2f})"
        healed = s["classes"].get("solved", 0) + s["classes"].get("suboptimal", 0)
        rest = ", ".join(f"{n} {CLASS[k]}" for k, n in sorted(s["classes"].items(), key=lambda kv: -kv[1])
                         if k not in ("solved",) and n) or "–"
        out.append(f"| {name} | {s['runs']} | **{s['reward_mean']:.2f}** ({lo:.2f} to {hi:.2f}) | {gap} | "
                   f"{healed} of {s['runs']} | {rest} |")
    return "\n".join(out)


def patients(rows: list[dict], summary: dict) -> str:
    suites = json.loads((ROOT / "tasks" / "suites.json").read_text())
    kinds = {}
    for i, q in enumerate(suites["remission-v3-dev"]["patients"], 1):
        kinds[f"dev-{i:02d}"] = q["archetype"]
    taken = defaultdict(int)                       # the test suite: the first 10 with at most 4 of one tumour type
    for q in suites["remission-v3-test"]["patients"]:
        if taken[q["archetype"]] < 4 and sum(taken.values()) < 10:
            taken[q["archetype"]] += 1
            kinds[f"test-{sum(taken.values()):02d}"] = q["archetype"]
    labels = [lab for lab in summary["labels"] if lab not in SHORTCUTS and not experience_arm(lab)]
    by = defaultdict(lambda: defaultdict(list))
    clones, best_shortcut = {}, defaultdict(float)
    for r in rows:
        if r["label"] in SHORTCUTS:
            best_shortcut[r["task"]] = max(best_shortcut[r["task"]], r["reward"])
            continue
        if experience_arm(r["label"]):
            continue
        by[r["task"]][r["label"]].append(r["reward"])
        clones[r["task"]] = ", ".join(f"{ALTERATION.get(c['alteration'], c['alteration'])} {100 * c['share'] + 1e-6:.1f} %"
                                      for c in r.get("minor_clones", []) if c.get("below_sequencing_limit", True))
    head = ["Patient", "Tumour", "Clone below the sequencing limit"] + labels + ["Best shortcut"]
    out = ["| " + " | ".join(head) + " |", "|" + "---|" * len(head)]
    for task in sorted(by):
        cells = [task, TUMOUR.get(kinds.get(task, ""), kinds.get(task, "")), clones.get(task, "")]
        for lab in labels:
            v = by[task].get(lab, [])
            if lab == summary["reference"] and len(set(round(x, 3) for x in v)) == 1:
                cells.append(f"{v[0]:.2f}")
            else:
                cells.append(" ".join(f"{x:.2f}" for x in v) or "–")
        cells.append(f"{best_shortcut[task]:.2f}")
        out.append("| " + " | ".join(cells) + " |")
    return "\n".join(out)


def failure_counts(rows: list[dict], summary: dict, experience: bool = False) -> dict:
    """Per agent: how the lost runs tie to what was done in the lab. experience=True adds the runs with notes."""
    out = {}
    for lab in summary["labels"]:
        if lab in SHORTCUTS or lab == summary["reference"] or (experience_arm(lab) and not experience):
            continue
        mine = [r for r in rows if r["label"] == lab]
        toxic = [r for r in mine if r["class"] == "too_toxic"]
        sub = [r for r in mine if r["class"] == "suboptimal"]
        med = lambda k: statistics.median(r["facts"].get(k) or 0 for r in mine) if mine else 0  # noqa: E731
        missed = [r for r in mine if r["class"] == "missed_clone"]
        out[lab] = {
            "runs": len(mine),
            "healed": sum(r["class"] in ("solved", "suboptimal") for r in mine),
            "suboptimal": len(sub),
            "suboptimal_untested": sum(not r["facts"].get("final_regimen_cultured_on_healthy_cells") for r in sub),
            "harm_mean": round(sum(r["harm"] for r in mine) / len(mine), 2) if mine else 0,
            "lab_days_median": med("lab_days_used"), "wall_minutes_median": med("wall_minutes"),
            "too_toxic": len(toxic),
            "too_toxic_untested": sum(not r["facts"].get("final_regimen_cultured_on_healthy_cells") for r in toxic),
            "missed_clone": len(missed),
            "missed_clone_never_named": sum(not all(c.get("seen") for c in r.get("minor_clones", [])
                                                    if c.get("below_sequencing_limit", True)) for r in missed),
            "tumour_survived": sum(r["class"] == "tumour_survived" for r in mine),
            "no_regimen": sum(r["class"] == "no_regimen" for r in mine),
            "no_line_sequenced": sum(not r["facts"].get("lines_sequenced") for r in mine),
            "wells_used_median": med("wells_used"),
        }
    return out


def failures(rows: list[dict], summary: dict) -> str:
    fc = failure_counts(rows, summary)
    a = ["| | Runs | Healed | of those: with avoidable harm | Tumour gone, normal tissue destroyed "
         "| Relapse from a hidden clone | of those: the clone was never named in the lab | Main tumour not cleared |",
         "|---|---|---|---|---|---|---|---|"]
    for lab, c in fc.items():
        a.append(f"| {lab} | {c['runs']} | {c['healed']} | {c['suboptimal']} | {c['too_toxic']} | {c['missed_clone']} | "
                 f"{c['missed_clone_never_named']} | {c['tumour_survived']} |")
    b = ["| | Mean harm | Runs that lost reward to harm: final combination never on healthy cells in a well "
         "| Median lab day at submission (of 28) | Median wells used (of 384) | Median minutes (of 60) |",
         "|---|---|---|---|---|---|"]
    ref = [r for r in rows if r["label"] == summary["reference"]]
    if ref:
        med = lambda k: statistics.median(r["facts"].get(k) or 0 for r in ref)  # noqa: E731
        b.append(f"| {summary['reference']} | {sum(r['harm'] for r in ref) / len(ref):.2f} | – | {med('lab_days_used'):g} | "
                 f"{med('wells_used'):g} | {med('wall_minutes'):.0f} |")
    for lab, c in fc.items():
        lost = c["too_toxic"] + c["suboptimal"]
        b.append(f"| {lab} | {c['harm_mean']:.2f} | {c['too_toxic_untested'] + c['suboptimal_untested']} of {lost} | "
                 f"{c['lab_days_median']:g} | {c['wells_used_median']:g} | {c['wall_minutes_median']:.0f} |")
    return "\n".join(a) + "\n\nHarm and how much of the budget was used:\n\n" + "\n".join(b)


def experience(rows: list[dict], summary: dict) -> str:
    """The same small model on the test patients: first time, and with the notes it wrote from its own runs."""
    fc = failure_counts(rows, summary, experience=True)
    arms = [lab for lab in fc if experience_arm(lab)]
    base = [lab for lab in fc if any(a.startswith(lab + " with its ") for a in arms)]
    out = ["| | Runs | Reward (95 % CI) | Healed | Tumour gone, normal tissue destroyed | Relapse from a hidden clone | Mean harm "
           "| Median lab day at submission | Median wells used |",
           "|---|---|---|---|---|---|---|---|---|"]
    for lab in base + arms:
        s, c = summary["labels"][lab], fc[lab]
        lo, hi = s["reward_ci95"]
        name = lab if lab in base else lab.split(" with its ")[1].join(["with its ", ""])
        out.append(f"| {name} | {c['runs']} | **{s['reward_mean']:.2f}** ({lo:.2f} to {hi:.2f}) | {c['healed']} | {c['too_toxic']} | "
                   f"{c['missed_clone']} | {c['harm_mean']:.2f} | {c['lab_days_median']:g} | {c['wells_used_median']:g} |")
    return "\n".join(out)


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--analysis", default=str(ROOT / "results" / "run_analysis.jsonl"))
    p.add_argument("--summary", default=str(ROOT / "results" / "summary.json"))
    p.add_argument("--readme", default=str(ROOT / "README.md"))
    p.add_argument("--print", action="store_true", help="print the tables, leave the README alone")
    a = p.parse_args()
    rows = [json.loads(x) for x in Path(a.analysis).read_text().splitlines() if x.strip()]
    summary = json.loads(Path(a.summary).read_text())
    blocks = {"results": results(summary), "patients": patients(rows, summary), "failures": failures(rows, summary),
              "experience": experience(rows, summary)}
    if a.print:
        print("\n\n".join(blocks.values()))
        return
    text = Path(a.readme).read_text()
    for key, table in blocks.items():
        text, n = re.subn(rf"(<!-- {key}:start -->\n).*?(<!-- {key}:end -->)", lambda m: m.group(1) + table + "\n" + m.group(2),
                          text, flags=re.S)
        if n != 1 and key != "experience":
            raise SystemExit(f"README has no <!-- {key}:start --> ... <!-- {key}:end --> block")
    Path(a.readme).write_text(text)
    print("README tables written")


if __name__ == "__main__":
    main()
