#!/bin/bash
# Verifier. Harbor uploads tests/ to /tests and runs this in `main` after the
# agent has finished. It asks the sidecar for the score and the episode log
# (admin token: /tests/admin_token) and writes
#   /logs/verifier/reward.txt    scalar reward
#   /logs/verifier/reward.json   reward + named sub-scores (Harbor prefers this)
#   /logs/verifier/episode.json  full episode log for the replay viewer
# A reward is written on every path; anything going wrong means 0.
set -u
OUT=/logs/verifier
mkdir -p "$OUT"

# Shared-verifier hardening: stop whatever the agent user left running and take
# the output directory away from it before anything is written.
if [ "$(id -u)" = "0" ] && id agent >/dev/null 2>&1; then
  pkill -KILL -u agent 2>/dev/null || true
  chown root:root "$OUT" 2>/dev/null || true
  chmod 755 "$OUT" 2>/dev/null || true
fi

rm -f "$OUT/reward.txt" "$OUT/reward.json"
echo 0 > "$OUT/reward.txt"

if ! python3 -I /tests/verify.py; then
  echo "verify.py failed; falling back to reward 0" >&2
  rm -f "$OUT/reward.json"
  echo 0 > "$OUT/reward.txt"
fi

[ -s "$OUT/reward.txt" ] || echo 0 > "$OUT/reward.txt"
echo "reward: $(cat "$OUT/reward.txt")"
exit 0
