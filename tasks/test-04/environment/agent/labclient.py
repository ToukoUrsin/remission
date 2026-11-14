"""labclient: call the Remission lab tools from your own Python scripts.

    from labclient import call
    brief = call("get_brief")
    e = call("start_experiment", kind="viability", cells="tumour",
             wells=[{}, {"RAFi": 0.1}, {"RAFi": 0.1, "MEKi": 0.01}], days=3, read_days=[1, 2, 3])
    call("advance_time", days=3)
    print(call("get_results", experiment_id=e["id"]))

These are exactly the MCP tools of the `remission` server (same budgets, same
rules). A rejected action returns {"error": "..."} instead of raising.
Standard library only.
"""
import itertools
import json
import os
import urllib.request

URL = os.environ.get("REMISSION_MCP_URL", "http://remission:8000/mcp")
_ids = itertools.count(1)


def call(tool, **arguments):
    body = json.dumps({"jsonrpc": "2.0", "id": next(_ids), "method": "tools/call",
                       "params": {"name": tool, "arguments": arguments}}).encode()
    req = urllib.request.Request(URL, data=body, headers={
        "Content-Type": "application/json", "Accept": "application/json, text/event-stream"})
    with urllib.request.urlopen(req, timeout=600) as resp:
        raw = resp.read().decode()
    if raw.lstrip().startswith("event:") or raw.lstrip().startswith("data:"):   # SSE framing
        raw = "".join(line[5:].strip() for line in raw.splitlines() if line.startswith("data:"))
    msg = json.loads(raw)
    if "error" in msg:
        raise RuntimeError(f"{tool}: {msg['error']}")
    result = msg["result"]
    text = "".join(c.get("text", "") for c in result.get("content", []) if c.get("type") == "text")
    if result.get("isError"):
        return {"error": text}
    try:
        return json.loads(text)
    except ValueError:
        return {"text": text}


if __name__ == "__main__":
    print(json.dumps(call("lab_status"), indent=1))
