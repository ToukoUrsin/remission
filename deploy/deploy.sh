#!/usr/bin/env bash
# Deploy the Remission viewer (replays + "Treat a patient") to Fly.io.
# Stages a minimal build context (viewer/, remission/, two JSON files), then `fly deploy`.
# Rebuild viewer/data first if the runs changed (scripts/build_viewer.py); this script ships whatever is there.
# Token: FLY_API_TOKEN in the environment, or through the vault: hsec exec --only FLY_API_TOKEN -- deploy/deploy.sh
set -euo pipefail
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo="$(dirname "$here")"
app="$(sed -n 's/^app = "\(.*\)"/\1/p' "$here/fly.toml")"
if [[ -z "${FLY_API_TOKEN:-}" ]] && command -v hsec >/dev/null 2>&1; then
  exec hsec exec --only FLY_API_TOKEN -- "$0" "$@"
fi
ctx="$(mktemp -d)"
trap 'rm -rf -- "${ctx:?}"' EXIT
mkdir -p "$ctx/reference" "$ctx/tasks"
tar -C "$repo" --exclude=__pycache__ --exclude='*.pyc' -cf - viewer remission | tar -C "$ctx" -xf -
cp "$repo/reference/fixed_regimens.json" "$ctx/reference/"
cp "$repo/tasks/suites.json" "$ctx/tasks/"
cp "$here/Dockerfile" "$here/fly.toml" "$ctx/"
apps="$(fly apps list 2>/dev/null || true)"   # into a variable: grep -q on a pipe can end fly early and look like "no app"
grep -q "^ *$app " <<<"$apps" || fly apps create "$app" --org "${FLY_ORG:-helios-one}"
cd "$ctx"
fly deploy --config fly.toml --app "$app" --remote-only --ha=false --yes "$@"
fly scale count 1 --app "$app" --yes >/dev/null 2>&1 || true
echo "https://$app.fly.dev/#view=treat"
