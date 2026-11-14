#!/usr/bin/env python3
"""Print a per-task results table from Harbor job directories.

    scripts/results_table.py LABEL=JOB_DIR [LABEL=JOB_DIR ...] [--key reward] [--json out.json] [--agent LABEL]

Each JOB_DIR is one `harbor run` output (<jobs_dir>/<job_name>); it holds one
directory per trial with a result.json. A trial without a verifier result
(agent setup crashed, verifier crashed) is shown as `err` and counted as 0 in
the means, which is also how Harbor aggregates.

In the --agent column a run that never called a tool is shown as `idle` and
named below the table. That is a setup problem (the model could not be
reached), not a result, and it would otherwise look like a fair 0.
"""

from __future__ import annotations

import argparse
import json
from collections import defaultdict
from pathlib import Path


IDLE = "idle"   # stands in for the reward of an agent run that never called a tool


def tool_calls(trial: Path) -> int | None:
    """How many tool calls reached the lab in this trial (None when the verifier left no episode)."""
    try:
        return len(json.loads((trial / "verifier" / "episode.json").read_text())["episode"]["calls"])
    except (OSError, ValueError, KeyError, TypeError):
        return None


def load(job_dir: Path, key: str, agent: bool = False) -> tuple[dict[str, list[float | str | None]], list[str]]:
    out: dict[str, list[float | str | None]] = defaultdict(list)
    idle = []
    for result in sorted(job_dir.glob("*/result.json")):
        r = json.loads(result.read_text())
        task = (r.get("task_name") or result.parent.name).split("/")[-1]
        rewards = (r.get("verifier_result") or {}).get("rewards")
        if agent and tool_calls(result.parent) == 0:
            idle.append(str(result.parent))
            out[task].append(IDLE)
        else:
            out[task].append(None if rewards is None else float(rewards.get(key, 0.0)))
    return out, idle


def num(v: float | str | None) -> float:
    return v if isinstance(v, float) else 0.0


def fmt(values: list[float | str | None]) -> str:
    if not values:
        return "-"

    def one(v):
        return "err" if v is None else v if isinstance(v, str) else f"{v:.2f}"

    cell = f"{sum(num(v) for v in values) / len(values):.3f}"
    if len(values) > 1:
        cell += " [" + " ".join(one(v) for v in values) + "]"
    elif not isinstance(values[0], float):
        cell = one(values[0])
    return cell


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("jobs", nargs="+", help="LABEL=JOB_DIR")
    p.add_argument("--key", default="reward", help="reward key to tabulate")
    p.add_argument("--json", type=Path, help="also write the table as JSON")
    p.add_argument("--agent", help="label of the column with real agent runs (runs without a tool call are flagged)")
    args = p.parse_args()

    columns: list[tuple[str, dict[str, list[float | str | None]]]] = []
    idle: list[str] = []
    for spec in args.jobs:
        label, _, path = spec.rpartition("=")  # labels may contain "="
        col, col_idle = load(Path(path), args.key, agent=label == args.agent)
        columns.append((label, col))
        idle += col_idle

    tasks = sorted({t for _, col in columns for t in col})
    rows = [[t] + [fmt(col.get(t, [])) for _, col in columns] for t in tasks]
    means = ["MEAN"]
    for _, col in columns:
        vals = [num(v) for vs in col.values() for v in vs]
        means.append(f"{sum(vals) / len(vals):.3f}" if vals else "-")
    header = ["task"] + [label for label, _ in columns]
    table = [header, *rows, means]
    widths = [max(len(r[i]) for r in table) for i in range(len(header))]
    for i, r in enumerate(table):
        print("  ".join(c.ljust(w) for c, w in zip(r, widths)))
        if i == 0 or i == len(table) - 2:
            print("  ".join("-" * w for w in widths))

    if idle:
        print(f"\n{len(idle)} agent run(s) never called a tool (shown as {IDLE}, counted as 0). The agent did not reach the model:")
        print("check the API key and, for a sealed run, the model's host (--model-host HOST, or --open). Agent logs:")
        for trial in idle:
            print(f"  {trial}/agent")

    if args.json:
        args.json.write_text(
            json.dumps({"key": args.key, "columns": {label: col for label, col in columns}, "idle": idle}, indent=1)
        )


if __name__ == "__main__":
    main()
