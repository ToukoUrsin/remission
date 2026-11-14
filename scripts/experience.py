#!/usr/bin/env python3
"""Learning from experience: a model plays the development patients, sees how its regimens
scored, writes down what it learned, and plays the test patients again with those notes.

No weights change. The experience is the model's own scored episodes on tasks/dev (never a
test patient); the notes are written by the same model; they reach the model as text at the
end of the task instruction.

    # 1. play the development patients, then explain every run (needs the repository)
    scripts/run_eval.sh --agent-only -k 1 -m anthropic/claude-haiku-5-5 -t exp-dev tasks/dev/dev-*
    scripts/analyze_runs.py "Haiku 5.5=jobs/exp-dev-agent" --out experience/dev_runs.jsonl

    # 2. the model reads its own post-mortems and writes the notes
    scripts/experience.py notes experience/dev_runs.jsonl --model claude-haiku-5-5 --out experience/notes.md

    # 3. copies of the test tasks with the notes appended to the instruction; play them
    scripts/experience.py tasks experience/notes.md tasks/test-* --out jobs/exp-tasks
    scripts/run_eval.sh --agent-only -k 2 -m anthropic/claude-haiku-5-5 -t exp-test jobs/exp-tasks/test-*

Step 2 calls the Anthropic Messages API (ANTHROPIC_API_KEY, optional ANTHROPIC_BASE_URL).
"""
from __future__ import annotations

import argparse
import json
import os
import shutil
import sys
import urllib.request
from pathlib import Path

HEADING = "## Notes from earlier patients"
ASK = """You are an AI agent that has just treated {n} patients in a simulated precision-oncology lab. This was the task:

<task>
{task}
</task>

After each patient the regimen you submitted was scored, and a post-mortem was written with what only the
simulator knows: which resistant clones were really there, what the best regimen found for that patient scores,
and what the same drugs would have scored at other doses. Here are your {n} episodes, one JSON object per line:

<episodes>
{episodes}
</episodes>

You will now treat new patients of the same kind, with the same tools, budget and reward. You will not see these
episodes again; only the notes you write now. Write the notes for yourself: what to do, in which order, what to
check before submitting, and which of your own mistakes to avoid, with the numbers that matter (doses, days,
wells). Use only what the episodes show. At most 450 words, plain text with short bullet points, no preamble."""


def episode(row: dict) -> dict:
    """What the model gets to see of one of its own runs."""
    f = row.get("facts", {})
    cf = row.get("counterfactuals", {})
    return {
        "patient": row["task"], "reward": row["reward"], "healed": row["healed"], "harm": row["harm"],
        "regimen_submitted": row.get("regimen"),
        "lab": {k: f.get(k) for k in ("lab_days_used", "wells_used", "n_experiments", "biopsy_sequenced", "lines_derived",
                                      "lines_sequenced", "tumour_wells", "healthy_wells", "healthy_combination_wells",
                                      "longest_tumour_culture_days", "final_regimen_cultured_on_tumour_cells",
                                      "final_regimen_cultured_on_healthy_cells")},
        "clones_really_there": [{k: c.get(k) for k in ("alteration", "share", "seen")} for c in row.get("minor_clones", [])],
        "best_regimen_known_scores": row.get("best_known"),
        "same_drugs_other_doses": {k: cf["best_rescaling"].get(k) for k in ("change", "reward", "healed", "harm")}
        if cf.get("best_rescaling") else None,
        "harm_if_one_drug_is_left_out": cf.get("harm_without_drug"),
        "post_mortem": row.get("summary"),
    }


def ask(model: str, prompt: str, max_tokens: int = 16000) -> str:   # room for the model's thinking and the notes
    base = os.environ.get("ANTHROPIC_BASE_URL", "https://api.anthropic.com").rstrip("/")
    base = base[:-3] if base.endswith("/v1") else base
    req = urllib.request.Request(
        base + "/v1/messages", method="POST",
        data=json.dumps({"model": model, "max_tokens": max_tokens, "messages": [{"role": "user", "content": prompt}]}).encode(),
        headers={"content-type": "application/json", "anthropic-version": "2023-06-01",
                 "x-api-key": os.environ["ANTHROPIC_API_KEY"]})
    with urllib.request.urlopen(req, timeout=600) as r:
        out = json.load(r)
    return "".join(b.get("text", "") for b in out["content"] if b.get("type") == "text").strip()


def notes(a) -> None:
    rows = [json.loads(line) for line in Path(a.analysis).read_text().splitlines() if line.strip()]
    if a.label:
        rows = [r for r in rows if r["label"] == a.label]
    if not rows:
        sys.exit("no runs in " + a.analysis)
    eps = [episode(r) for r in rows]
    task = (Path(__file__).resolve().parent.parent / "tasks" / "_template" / "instruction.md").read_text()
    prompt = ASK.format(n=len(eps), task=task.strip(), episodes="\n".join(json.dumps(e) for e in eps))
    text = ask(a.model, prompt)
    if not text:
        sys.exit("the model returned no text")
    out = Path(a.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(text + "\n")
    out.with_suffix(".prompt.txt").write_text(prompt + "\n")
    print(f"{out}: {len(text.split())} words from {len(eps)} episodes (mean reward {sum(e['reward'] for e in eps) / len(eps):.3f})")


def tasks(a) -> None:
    text = Path(a.notes).read_text().strip()
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    for t in a.task_dirs:
        src = Path(t)
        dst = out / src.name
        if dst.exists():
            shutil.rmtree(dst)
        shutil.copytree(src, dst)
        ins = dst / "instruction.md"
        ins.write_text(ins.read_text().rstrip() + f"\n\n{HEADING}\n\nYou wrote these notes after treating earlier patients of the "
                       "same kind in this lab and seeing how your regimens scored.\n\n" + text + "\n")
        print(dst)


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)
    n = sub.add_parser("notes", help="the model writes its notes from its analysed development runs")
    n.add_argument("analysis", help="jsonl from scripts/analyze_runs.py --out")
    n.add_argument("--label", help="only the runs with this label")
    n.add_argument("--model", default="claude-haiku-5-5")
    n.add_argument("--out", default="experience/notes.md")
    n.set_defaults(fn=notes)
    t = sub.add_parser("tasks", help="copies of task directories with the notes appended to the instruction")
    t.add_argument("notes")
    t.add_argument("task_dirs", nargs="+")
    t.add_argument("--out", required=True)
    t.set_defaults(fn=tasks)
    a = p.parse_args()
    a.fn(a)


if __name__ == "__main__":
    main()
