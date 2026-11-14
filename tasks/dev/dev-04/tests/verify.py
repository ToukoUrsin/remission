"""Fetch score + episode log from the sidecar and write Harbor's reward files.

Standard library only. Exits non-zero on any failure; test.sh then writes 0.
"""

import json
import os
import sys
import time
import urllib.error
import urllib.request
from pathlib import Path

OUT = Path("/logs/verifier")
URL = os.environ.get("REMISSION_ADMIN_URL", "http://remission:8000/admin/episode")
TOKEN = (Path(__file__).parent / "admin_token").read_text().strip()


def fetch() -> dict:
    last = None
    for attempt in range(10):
        try:
            req = urllib.request.Request(URL, headers={"Authorization": f"Bearer {TOKEN}"})
            with urllib.request.urlopen(req, timeout=600) as resp:
                return json.loads(resp.read().decode())
        except urllib.error.HTTPError as exc:  # wrong token etc.: retrying won't help
            raise SystemExit(f"admin route answered HTTP {exc.code}")
        except Exception as exc:  # noqa: BLE001 - sidecar restarting / DNS hiccup
            last = exc
            time.sleep(1 + attempt)
    raise SystemExit(f"sidecar unreachable: {last!r}")


def main() -> int:
    data = fetch()
    score = data["score"]
    reward = float(score["reward"])
    if not (0.0 <= reward <= 1.0):
        raise SystemExit(f"reward out of range: {reward}")

    episode = data["episode"]
    (OUT / "episode.json").write_text(json.dumps(data))
    # Harbor wants a flat {name: number} object; `reward` is the headline score.
    flat = {k: float(v) for k, v in score.items() if isinstance(v, (int, float))}
    flat["n_tool_calls"] = float(episode["n_calls"])
    (OUT / "reward.json").write_text(json.dumps(flat))
    (OUT / "reward.txt").write_text(f"{reward}\n")
    print(json.dumps({"score": score, "closed_by": episode["closed_by"], "n_calls": episode["n_calls"]}))
    return 0


if __name__ == "__main__":
    sys.exit(main())
