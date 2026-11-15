#!/usr/bin/env bash
# Run the evaluation over Harbor task directories and print a results table:
#
#   reference    the scripted scientist (Harbor's oracle agent, public tools only)
#   do-nothing   Harbor's nop agent
#   stop-early   reads the brief and quits
#   guideline    submits the standard-of-care regimen for the tumour type, no experiment
#   blind        submits the constant regimen that is best on average for the tumour type, no experiment
#   screen       a three-day drug screen, then the best-looking pair at label doses
#   hedge        sequencing only, then the regimen the public model rates best on average over every
#                clone that could be hiding (regimen_search.py on the public typical values)
#   agent        k attempts of a real agent
#
#   scripts/run_eval.sh [options] [TASK_DIR ...]
#
#   -k N          agent attempts per task                (default 3)
#   -a AGENT      Harbor agent                           (default opencode)
#   -m MODEL      model                                  (default anthropic/claude-fable-5-2)
#   -n N          concurrent trials                      (default 4)
#   -o DIR        jobs directory                         (default $JOBS_DIR or ./jobs)
#   -t TAG        prefix for the job names               (default eval-<timestamp>)
#   --no-agent    reference and baselines only (no model, no API key needed)
#   --agent-only  agent runs only (baselines already done under the same TAG)
#   --fresh N     run on N new patients instead of the ones in tasks/. They are drawn
#                 from a random salt on this machine, so nobody has seen them
#                 (about a minute per patient; uses .venv or python3 with numpy,
#                 else a python:3.12 container)
#   --open        let the agent keep its internet. Without it the agent plays
#                 sealed, as the tasks are written: its network is cut down to the
#                 lab, the model's host and OpenCode's model list (Harbor
#                 allowlist, local Docker). The host is taken from the model name
#                 or *_BASE_URL; name others with --model-host HOST
#
# Without TASK_DIRs: every tasks/*/ except _template.
# The agent phase uses the credentials in the environment (ANTHROPIC_API_KEY,
# OPENAI_API_KEY, optional ANTHROPIC_BASE_URL / OPENAI_BASE_URL, ...).
#
# Examples
#   scripts/run_eval.sh --no-agent                                  # reference + baselines on all tasks
#   ANTHROPIC_API_KEY=... scripts/run_eval.sh tasks/test-04          # one task: baselines + 3 agent runs
#   OPENAI_API_KEY=... scripts/run_eval.sh -m openai/gpt-6.1-sol tasks/test-*
#   ANTHROPIC_API_KEY=... scripts/run_eval.sh --fresh 3              # patients nobody has seen
set -euo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SCRIPTS="$REPO/scripts"
K=3; AGENT=opencode; MODEL=anthropic/claude-fable-5-2; NCONC=4
JOBS_DIR="${JOBS_DIR:-$REPO/jobs}"
TAG="eval-$(date +%Y%m%d-%H%M%S)"; RUN_AGENT=1; RUN_BASE=1; FRESH=0; SEALED=1
TASKS=(); MODEL_HOSTS=()
HARBOR="${HARBOR:-harbor}"

while [ $# -gt 0 ]; do
  case "$1" in
    -k) K="$2"; shift 2 ;;
    -a) AGENT="$2"; shift 2 ;;
    -m) MODEL="$2"; shift 2 ;;
    -n) NCONC="$2"; shift 2 ;;
    -o) JOBS_DIR="$2"; shift 2 ;;
    -t) TAG="$2"; shift 2 ;;
    --no-agent) RUN_AGENT=0; shift ;;
    --agent-only) RUN_BASE=0; shift ;;
    --fresh) FRESH="$2"; shift 2 ;;
    --open) SEALED=0; shift ;;
    --model-host) MODEL_HOSTS+=("$2"); shift 2 ;;
    -h|--help) sed -n '2,43p' "$0"; exit 0 ;;
    -*) echo "unknown option $1" >&2; exit 2 ;;
    *) TASKS+=("$1"); shift ;;
  esac
done

command -v "$HARBOR" >/dev/null || { echo "harbor not found: uv tool install harbor==0.25.0" >&2; exit 1; }
mkdir -p "$JOBS_DIR"
JOBS_DIR="$(cd "$JOBS_DIR" && pwd)"

if [ "$FRESH" -gt 0 ]; then
  [ ${#TASKS[@]} -eq 0 ] || { echo "--fresh makes its own tasks; drop the TASK_DIRs" >&2; exit 2; }
  PY="${PYTHON:-python3}"; [ -z "${PYTHON:-}" ] && [ -x "$REPO/.venv/bin/python" ] && PY="$REPO/.venv/bin/python"
  FRESH_DIR="$JOBS_DIR/$TAG-fresh-tasks"; mkdir -p "$FRESH_DIR"
  echo "== drawing $FRESH new patient(s) into $FRESH_DIR" >&2
  if "$PY" -c 'import numpy' 2>/dev/null; then
    while IFS= read -r d; do TASKS+=("$d"); done \
      < <("$PY" "$SCRIPTS/new_task.py" --fresh "$FRESH" --prefix fresh --out "$FRESH_DIR")
  else
    echo "   (no numpy here, using a python:3.12 container)" >&2
    while IFS= read -r d; do TASKS+=("$FRESH_DIR/$(basename "$d")"); done \
      < <(docker run --rm -u "$(id -u):$(id -g)" -e HOME=/tmp -v "$REPO":/repo:ro -v "$FRESH_DIR":/out -w /tmp \
            python:3.12-slim sh -c 'pip install -q --disable-pip-version-check --target /tmp/deps numpy==2.5.4 >&2 \
              && PYTHONPATH=/tmp/deps PYTHONDONTWRITEBYTECODE=1 python /repo/scripts/new_task.py --fresh "$0" --prefix fresh --out /out' "$FRESH")
  fi
  [ ${#TASKS[@]} -eq "$FRESH" ] || { echo "could not make the fresh tasks" >&2; exit 1; }
fi

if [ ${#TASKS[@]} -eq 0 ]; then
  for d in "$REPO"/tasks/*/; do
    [ -f "$d/task.toml" ] && [ "$(basename "$d")" != "_template" ] && TASKS+=("${d%/}")
  done
  [ ${#TASKS[@]} -gt 0 ] || { echo "no tasks in $REPO/tasks; make some with scripts/new_task.py or pass --fresh N" >&2; exit 1; }
fi

WORK="$(mktemp -d)"
cleanup() { rm -rf "${WORK:?}"; }
trap cleanup EXIT

# A job config that only lists the tasks; everything else comes from flags.
python3 - "$WORK/tasks.json" "${TASKS[@]}" <<'PY'
import json, sys
from pathlib import Path
out, *tasks = sys.argv[1:]
paths = [str(Path(t).resolve()) for t in tasks]
missing = [p for p in paths if not (Path(p) / "task.toml").is_file()]
if missing:
    sys.exit(f"not a task directory: {missing}")
Path(out).write_text(json.dumps({"tasks": [{"path": p} for p in paths]}))
PY

# The tasks ship sealed: task.toml lets the agent phase reach the lab's MCP server
# and nothing else. The agent plays on them like that, plus the model's host and
# OpenCode's model list. Two sets of copies are made here and nothing else about
# a task changes:
#   sealed   the seal written in (so a task made without it is sealed too)
#   open     the seal taken out: for the agent with --open, and for the reference
#            and the shortcut scripts, which are ours and call only the lab (they
#            ran like this for the reported numbers)
python3 - "$WORK/tasks.json" "$WORK" <<'PY'
import json, re, shutil, sys
from pathlib import Path
cfg, work = Path(sys.argv[1]), Path(sys.argv[2])
SEAL = 'network_mode = "allowlist"\nallowed_hosts = ["remission"]\n'
AGENT = re.compile(r"^\[agent\][^\n]*\n(.*?)(?=^\[|\Z)", re.M | re.S)     # the body of the [agent] table
SEAL_LINE = re.compile(r"^(?:network_mode|allowed_hosts)[ \t]*=.*\n", re.M)
lists = {"sealed": [], "open": []}
for t in json.loads(cfg.read_text())["tasks"]:
    src = Path(t["path"])
    toml = (src / "task.toml").read_text()
    m = AGENT.search(toml)
    if not m:
        sys.exit(f"{src}: no [agent] section in task.toml")
    body = SEAL_LINE.sub("", m.group(1))
    for kind, new_body in (("sealed", SEAL + body), ("open", body)):
        dst = work / kind / src.name
        shutil.copytree(src, dst)
        (dst / "task.toml").write_text(toml[:m.start(1)] + new_body + toml[m.end(1):])
        lists[kind].append({"path": str(dst)})
for kind, tasks in lists.items():
    (work / f"tasks-{kind}.json").write_text(json.dumps({"tasks": tasks}))
PY
SEAL_ARGS=(); AGENT_CFG="$WORK/tasks-open.json"
if [ "$SEALED" = 1 ] && [ "$RUN_AGENT" = 1 ]; then
  if [ ${#MODEL_HOSTS[@]} -eq 0 ]; then
    case "$MODEL" in
      anthropic/*) url="${ANTHROPIC_BASE_URL:-https://api.anthropic.com}" ;;
      openai/*) url="${OPENAI_BASE_URL:-https://api.openai.com}" ;;
      google/*) url="https://generativelanguage.googleapis.com" ;;
      openrouter/*) url="https://openrouter.ai" ;;
      *) echo "the agent plays sealed: name the model's host with --model-host HOST, or pass --open" >&2; exit 2 ;;
    esac
    MODEL_HOSTS=("$(python3 -c 'import sys, urllib.parse as u; print(u.urlparse(sys.argv[1]).hostname)' "$url")")
  fi
  AGENT_CFG="$WORK/tasks-sealed.json"
  for h in "${MODEL_HOSTS[@]}" models.dev; do SEAL_ARGS+=(--allow-agent-host "$h"); done
  echo "== sealed: the agent reaches the lab, the model (${MODEL_HOSTS[*]}) and models.dev, nothing else. --open lifts this" >&2
fi

cd "$JOBS_DIR"   # Harbor drops incidental files into the cwd; keep them out of the repo

run_phase() {   # run_phase <label> <attempts> <harbor args...>
  local label="$1" attempts="$2" cfg="$WORK/tasks-open.json"; shift 2
  [ "$label" = agent ] && cfg="$AGENT_CFG"
  echo "== $label: ${#TASKS[@]} task(s) x $attempts attempt(s)" >&2
  PYTHONPATH="$SCRIPTS${PYTHONPATH:+:$PYTHONPATH}" "$@" \
    -c "$cfg" -o "$JOBS_DIR" --job-name "$TAG-$label" \
    -k "$attempts" -n "$NCONC" -y -q >"$WORK/$label.log" 2>&1 \
    || { echo "harbor failed for $label:" >&2; tail -30 "$WORK/$label.log" >&2; }
}

if [ "$RUN_BASE" = 1 ]; then
  run_phase oracle 1 "$HARBOR" run -a oracle
  run_phase nop 1 "$HARBOR" run -a nop
  run_phase stop-early 1 "$HARBOR" run -a stop_early_agent:StopEarlyAgent
  # Shortcuts that avoid the task. They run the same way as the reference
  # (solution/solve.py inside the agent container, public tools only).
  run_phase guideline 1 "$HARBOR" run -a oracle --ae REMISSION_POLICY=guideline
  run_phase blind 1 "$HARBOR" run -a oracle --ae REMISSION_POLICY=blind
  run_phase screen 1 "$HARBOR" run -a oracle --ae REMISSION_POLICY=screen
  run_phase hedge 1 "$HARBOR" run -a oracle --ae REMISSION_POLICY=hedge
fi

COLUMNS=("reference=$JOBS_DIR/$TAG-oracle" "do-nothing=$JOBS_DIR/$TAG-nop" "stop-early=$JOBS_DIR/$TAG-stop-early"
         "guideline=$JOBS_DIR/$TAG-guideline" "blind=$JOBS_DIR/$TAG-blind" "screen=$JOBS_DIR/$TAG-screen"
         "hedge=$JOBS_DIR/$TAG-hedge")
AGENT_COL=()

if [ "$RUN_AGENT" = 1 ]; then
  run_phase agent "$K" "$HARBOR" run -a "$AGENT" -m "$MODEL" ${SEAL_ARGS[@]+"${SEAL_ARGS[@]}"}
  COLUMNS+=("$AGENT(k=$K)=$JOBS_DIR/$TAG-agent")
  AGENT_COL=(--agent "$AGENT(k=$K)")
fi

echo
echo "reward by task  (agent: $AGENT / $MODEL; jobs: $JOBS_DIR/$TAG-*)"
python3 "$SCRIPTS/results_table.py" "${COLUMNS[@]}" --json "$JOBS_DIR/$TAG-table.json" ${AGENT_COL[@]+"${AGENT_COL[@]}"}
