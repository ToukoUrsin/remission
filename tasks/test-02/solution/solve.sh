#!/bin/bash
# Reference run (Harbor's oracle agent executes this inside `main`).
# It plays through the MCP tools only: same URL, same tools, same client
# library and the same public model files the agent has in /app. No seed, no
# token, no server internals.
#   REMISSION_POLICY=scientist (default) | guideline | blind | screen
set -euo pipefail
cd "$(dirname "$0")"
PYTHONPATH=/app python3 solve.py
