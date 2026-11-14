#!/usr/bin/env python3
"""Stamp self-contained Harbor task directories out of tasks/_template, one patient each.

    scripts/new_task.py --salt remission-v3-test --prefix test         -> tasks/test-01 .. test-10
    scripts/new_task.py --fresh 3 --out /tmp/fresh                     -> 3 patients nobody has seen
    scripts/new_task.py --seed 123456 --name demo                      -> one given patient seed

A salt names a suite: its patients are the first `count` generated patients that
pass the cohort test (remission/cohort.py: seven rules). --fresh draws a random
salt, so the patients are new.

The cohort test is slow: about 1 generated patient in 400 passes, and the last
rule costs minutes of one core per candidate (plan on roughly one core-hour per
cohort patient passing rules 1 to 6, and about 4 in 10 of those pass rule 7).
The cohort patients of the named salts are therefore recorded in
tasks/suites.json (index in the salt's sequence, seed, the cohort test's
verdict); --salt reads them from there. --rescreen repeats the test instead.

What is different per task:
  environment/docker-compose.yaml   REMISSION_SEED, REMISSION_ADMIN_TOKEN and REMISSION_BEST_KNOWN
                                    (sidecar-only environment variables)
  tests/admin_token                 the same token, for the verifier
  task.toml                         [task].name and [metadata].seed_label

The seed is not in task.toml, instruction.md or anything else the agent's container can see.
"""
from __future__ import annotations

import argparse
import json
import re
import secrets
import shutil
import sys
from multiprocessing import Pool
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))
TEMPLATE = REPO / "tasks" / "_template"
PKG = REPO / "remission"
SIDECAR_FILES = ("__init__.py", "model.py", "patients.py", "lab.py", "server.py")
NAME_RE = re.compile(r"^[a-z0-9][a-z0-9._-]*$")


def sub_once(pattern: str, repl: str, text: str, what: str) -> str:
    new, n = re.subn(pattern, lambda m: m.group(1) + repl, text, count=1, flags=re.MULTILINE)
    if n != 1:
        sys.exit(f"template changed: could not find {what}")
    return new


def stamp(seed: int, name: str, best: dict | None, out_root: Path, force: bool, org: str = "remission",
          wells: int = 384, lab_days: int = 28) -> Path:
    if not NAME_RE.match(name):
        sys.exit(f"invalid task name {name!r}: use lowercase letters, digits, '.', '_' and '-'")
    dest = out_root / name
    if dest.exists():
        if not force:
            sys.exit(f"{dest} exists (use --force to overwrite)")
        shutil.rmtree(dest)
    shutil.copytree(TEMPLATE, dest, ignore=shutil.ignore_patterns("__pycache__", "*.pyc"))
    token = secrets.token_hex(24)

    compose_path = dest / "environment" / "docker-compose.yaml"
    compose = compose_path.read_text()
    compose = sub_once(r'^(\s*REMISSION_SEED:\s*).*$', f'"{seed}"', compose, "REMISSION_SEED")
    compose = sub_once(r'^(\s*REMISSION_ADMIN_TOKEN:\s*).*$', f'"{token}"', compose, "REMISSION_ADMIN_TOKEN")
    compose = sub_once(r'^(\s*REMISSION_BEST_KNOWN:\s*).*$', json.dumps(json.dumps(best or {})), compose, "REMISSION_BEST_KNOWN")
    compose = sub_once(r'^(\s*REMISSION_WELLS:\s*).*$', f'"{wells}"', compose, "REMISSION_WELLS")
    compose = sub_once(r'^(\s*REMISSION_LAB_DAYS:\s*).*$', f'"{lab_days}"', compose, "REMISSION_LAB_DAYS")
    compose_path.write_text(compose)
    ins = dest / "instruction.md"
    text = ins.read_text()
    if "**28 lab days**" not in text or "**384 wells**" not in text or "before day 28" not in text:
        sys.exit("template changed: could not find the lab budget in instruction.md")
    ins.write_text(text.replace("**28 lab days**", f"**{lab_days} lab days**").replace("**384 wells**", f"**{wells} wells**")
                   .replace("before day 28", f"before day {lab_days}"))

    sidecar = dest / "environment" / "sidecar" / "remission"
    sidecar.mkdir(parents=True, exist_ok=True)
    for f in SIDECAR_FILES:                                       # the hidden side: generator, lab, server
        shutil.copy(PKG / f, sidecar / f)
    agent = dest / "environment" / "agent"                        # the public side
    shutil.copy(PKG / "model.py", agent / "cellmodel.py")
    shutil.copy(PKG / "oracle.py", agent / "regimen_search.py")
    for src in sorted((REPO / "reference").glob("*.py")) + sorted((REPO / "reference").glob("*.json")):
        shutil.copy(src, dest / "solution" / src.name)            # the reference: public tools only

    (dest / "tests" / "admin_token").write_text(token + "\n")
    toml_path = dest / "task.toml"
    toml = toml_path.read_text()
    toml = sub_once(r'^(name = )".*"$', f'"{org}/{name}"', toml, "[task].name")
    toml = sub_once(r'^(seed_label = )".*"$', f'"{name}"', toml, "[metadata].seed_label")
    toml_path.write_text(toml)
    for script in (dest / "tests" / "test.sh", dest / "solution" / "solve.sh"):
        script.chmod(0o755)
    return dest


def _screen(args):
    from remission import cohort as C
    from remission import patients as P
    salt, i = args
    seed = P.seed_from(salt, i)
    res = C.screen(P.make_patient(seed))
    return i, seed, res


SUITES = REPO / "tasks" / "suites.json"


def capped(rows: list[tuple[int, int, dict]], count: int, per_type: int) -> list[tuple[int, dict]] | None:
    """The first `count` of rows [(index, seed, result)] in index order, at most `per_type` of one tumour type."""
    out, seen = [], {}
    for _, seed, res in sorted(rows, key=lambda r: r[0]):
        kind = res["archetype"]
        if seen.get(kind, 0) < per_type:
            seen[kind] = seen.get(kind, 0) + 1
            out.append((seed, res))
            if len(out) == count:
                return out
    return None


def recorded(salt: str, count: int, per_type: int) -> list[tuple[int, dict]] | None:
    """The suite of a salt as tasks/suites.json records its cohort patients, or None."""
    try:
        suite = json.loads(SUITES.read_text())[salt]
    except (OSError, KeyError, ValueError):
        return None
    return capped([(p["index"], p["seed"], p) for p in suite["patients"] if p.get("eligible", True)], count, per_type)


def cohort(salt: str, count: int, per_type: int, procs: int) -> list[tuple[int, dict]]:
    """The suite of a salt, by running the cohort test: [(seed, screen result)]."""
    found: list[tuple[int, int, dict]] = []
    start, chunk = 0, max(procs, 8) * 8
    with Pool(procs) as pool:
        while True:
            res = pool.map(_screen, [(salt, i) for i in range(start, start + chunk)])
            found += [r for r in res if r[2]["eligible"]]
            start += chunk
            print(f"  screened {start} generated patients, {len(found)} in the cohort", file=sys.stderr)
            got = capped(found, count, per_type)
            if got:
                return got


def main() -> None:
    import os
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    g = p.add_mutually_exclusive_group(required=True)
    g.add_argument("--salt", help="suite name; its patients are the first --count cohort patients of this salt")
    g.add_argument("--fresh", type=int, metavar="N", help="N new patients from a random salt")
    g.add_argument("--seed", type=int, help="one given patient seed (no cohort test)")
    p.add_argument("--count", type=int, default=10)
    p.add_argument("--per-type", type=int, default=1000, help="at most this many patients of one tumour type in a suite (default: no cap)")
    p.add_argument("--name", help="task directory name (with --seed)")
    p.add_argument("--prefix", default="patient", help="name prefix: <prefix>-01, <prefix>-02, ...")
    p.add_argument("--out", type=Path, default=REPO / "tasks", help="directory to create the tasks in")
    p.add_argument("--procs", type=int, default=min(os.cpu_count() or 4, 32))
    p.add_argument("--force", action="store_true", help="overwrite existing task directories")
    p.add_argument("--wells", type=int, default=384, help="lab budget in wells")
    p.add_argument("--lab-days", type=int, default=28, help="lab days before the regimen is due")
    p.add_argument("--rescreen", action="store_true", help="with --salt: repeat the cohort test instead of reading tasks/suites.json")
    a = p.parse_args()

    a.out.mkdir(parents=True, exist_ok=True)
    if a.seed is not None:
        print(stamp(a.seed, a.name or f"{a.prefix}-{a.seed}", None, a.out, a.force, wells=a.wells, lab_days=a.lab_days))
        return
    salt, count = (a.salt, a.count) if a.salt else ("fresh-" + secrets.token_hex(8), a.fresh)
    found = (None if a.fresh or a.rescreen else recorded(salt, count, a.per_type)) or cohort(salt, count, a.per_type, a.procs)
    for n, (seed, res) in enumerate(found, start=1):
        tol = res.get("tolerable") or {"regimen": res["best"], "reward": res["best_reward"]}
        best = {"regimen": tol["regimen"], "reward": tol["reward"]}
        print(stamp(seed, f"{a.prefix}-{n:02d}", best, a.out, a.force, wells=a.wells, lab_days=a.lab_days))


if __name__ == "__main__":
    main()
