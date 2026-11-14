#!/usr/bin/env bash
# Make a host-local model endpoint (default 127.0.0.1:8317) reachable from
# Harbor's task containers.
#
# Containers cannot reach the host's 127.0.0.1. This forwards
#   <docker0 gateway ip>:PORT  ->  127.0.0.1:PORT
# The listener is bound to the docker bridge address only, so it is reachable
# from local containers but not from the public interface.
#
#   scripts/model_proxy_forward.sh start    # idempotent, prints the base URL
#   scripts/model_proxy_forward.sh stop
#   scripts/model_proxy_forward.sh url
set -euo pipefail

PORT="${MODEL_POOL_PORT:-8317}"
BRIDGE_IP="${MODEL_POOL_BRIDGE_IP:-$(docker network inspect bridge -f '{{(index .IPAM.Config 0).Gateway}}')}"
PIDFILE="${TMPDIR:-/tmp}/remission-model-proxy-forward.pid"

running() { [ -f "$PIDFILE" ] && kill -0 "$(cat "$PIDFILE")" 2>/dev/null; }
listening() { (exec 3<>"/dev/tcp/${BRIDGE_IP}/${PORT}") 2>/dev/null; }

case "${1:-start}" in
  start)
    if ! running && ! listening; then
      command -v socat >/dev/null || { echo "socat is required (apt-get install socat)" >&2; exit 1; }
      nohup socat "TCP-LISTEN:${PORT},bind=${BRIDGE_IP},fork,reuseaddr" "TCP:127.0.0.1:${PORT}" >/dev/null 2>&1 &
      echo $! > "$PIDFILE"
      sleep 0.5
      running || { echo "forwarder failed to start (is ${BRIDGE_IP}:${PORT} taken?)" >&2; exit 1; }
    fi
    echo "http://${BRIDGE_IP}:${PORT}"
    ;;
  stop)
    if running; then kill "$(cat "$PIDFILE")"; fi
    rm -f "$PIDFILE"
    ;;
  url) echo "http://${BRIDGE_IP}:${PORT}" ;;
  status) if running || listening; then echo "running: http://${BRIDGE_IP}:${PORT}"; else echo "stopped"; exit 1; fi ;;
  *) echo "usage: $0 start|stop|url|status" >&2; exit 2 ;;
esac
