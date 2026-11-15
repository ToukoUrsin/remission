#!/usr/bin/env python3
"""Did the 60-minute session limit shape the agents' runs?

For every reported agent run this reads the agent's own messages and shell commands from the Harbor trial
folder and records two things:

  mentions_limit   a message names the time limit (the sentence is kept as `quote`)
  date_calls       how many shell commands ran `date`, the only way to see elapsed time in these runs
                   (`lab_status` shows lab days and wells, not minutes)

    scripts/clock_mentions.py JOBS_DIR [--out results/clock_mentions.jsonl]

JOBS_DIR holds the Harbor job folders (<job>/<trial>/agent/trajectory.json). The trials are the ones named in
results/run_analysis.jsonl. Without JOBS_DIR it only prints the summary of the saved file.
"""
import argparse
import collections
import glob
import json
import re
import statistics
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
AGENTS = ("Claude Fable 5.2", "Claude Opus 5.5", "GPT-6.1 Sol", "Claude Haiku 5.5")
LIMIT = re.compile(r"(60[- ]minute|sixty[- ]minute|wall[- ]clock|time limit|time budget|out of time|running out of time|"
                   r"minutes? (?:left|remaining)|time (?:left|remaining)|time constraint|time pressure|"
                   r"session (?:limit|time)|short on time|running short)", re.I)
DATE = re.compile(r"(?:^|[;&|\n(]\s*)date(?:\s|$)", re.M)


def walk(o, texts, cmds):
    if isinstance(o, dict):
        if (o.get("role") or o.get("source")) in ("assistant", "agent"):
            for k in ("content", "text", "message"):
                v = o.get(k)
                if isinstance(v, str) and v.strip():
                    texts.append(v)
                elif isinstance(v, list):
                    texts.extend(q["text"] for q in v if isinstance(q, dict) and isinstance(q.get("text"), str))
        for k in ("command", "cmd"):
            if isinstance(o.get(k), str):
                cmds.append(o[k])
        for v in o.values():
            walk(v, texts, cmds)
    elif isinstance(o, list):
        for v in o:
            walk(v, texts, cmds)


def sentence(text, m):
    """the sentence around a match, one line"""
    a = max(text.rfind(c, 0, m.start()) for c in ".\n!?") + 1
    ends = [i for i in (text.find(c, m.end()) for c in ".\n!?") if i >= 0]
    b = min(ends) + 1 if ends else len(text)
    return " ".join(text[a:b].split())


def scan(jobs):
    rows = [json.loads(line) for line in open(REPO / "results" / "run_analysis.jsonl")]
    ref = collections.defaultdict(list)
    for r in rows:
        if r["label"] == "Reference scientist":
            ref[r["task"]].append(r["reward"])
    out = []
    for r in rows:
        if r["label"] not in AGENTS:
            continue
        hits = glob.glob(f"{jobs}/*/{r['trial']}/agent/trajectory.json")
        if not hits:
            raise SystemExit(f"no trajectory for {r['trial']} under {jobs}")
        texts, cmds = [], []
        walk(json.load(open(hits[0])), texts, cmds)
        quotes = []
        for t in texts:
            for m in LIMIT.finditer(t):
                s = sentence(t, m)
                if s not in quotes:
                    quotes.append(s)
        out.append({"agent": r["label"], "task": r["task"], "trial": r["trial"], "reward": r["reward"],
                    "reference_reward": round(statistics.mean(ref[r["task"]]), 4),
                    "minutes": r["facts"]["wall_minutes"], "lab_day_at_submission": r["facts"]["lab_days_used"],
                    "wells_used": r["facts"]["wells_used"], "mentions_limit": bool(quotes),
                    "date_calls": sum(1 for c in cmds if DATE.search(c)), "quotes": quotes})
    return out


def summary(rows):
    med = statistics.median
    print(f"{'':18} {'runs':>4}  {'reward':>6}  {'below ref':>9}  {'minutes (median, range)':>24}  {'lab day':>7}  {'wells':>6}")
    for a in AGENTS:
        v = [r for r in rows if r["agent"] == a]
        groups = [("all runs", v), ("name the limit", [r for r in v if r["mentions_limit"]]),
                  ("do not", [r for r in v if not r["mentions_limit"]]),
                  ("ran `date`", [r for r in v if r["date_calls"]]),
                  ("never ran `date`", [r for r in v if not r["date_calls"]])]
        print(a)
        for name, s in groups:
            if not s:
                print(f"  {name:16} {0:>4}")
                continue
            gap = statistics.mean(r["reference_reward"] - r["reward"] for r in s)
            mins = [r["minutes"] for r in s]
            print(f"  {name:16} {len(s):>4}  {statistics.mean(r['reward'] for r in s):>6.3f}  {gap:>9.3f}  "
                  f"{med(mins):>8.1f} ({min(mins):.1f} to {max(mins):.1f})".ljust(72)
                  + f"{med(r['lab_day_at_submission'] for r in s):>7g}  {med(r['wells_used'] for r in s):>6g}")
        both = sum(1 for r in v if r["mentions_limit"] and r["date_calls"])
        print(f"  name the limit and ran `date`: {both}")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("jobs", nargs="?")
    ap.add_argument("--out", default=str(REPO / "results" / "clock_mentions.jsonl"))
    a = ap.parse_args()
    if a.jobs:
        rows = scan(a.jobs)
        with open(a.out, "w") as f:
            for r in rows:
                f.write(json.dumps(r, ensure_ascii=False) + "\n")
    else:
        rows = [json.loads(line) for line in open(a.out)]
    summary(rows)


if __name__ == "__main__":
    main()
