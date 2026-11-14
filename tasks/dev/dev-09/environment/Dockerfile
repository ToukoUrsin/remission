# The agent's container (`main`). It holds no server code, no seed and no admin
# token: only a Python MCP client, the public cell model and the regimen search,
# so that anything the reference solution does, the agent can do too.
FROM python:3.12-slim

ENV PIP_NO_CACHE_DIR=1 PIP_DISABLE_PIP_VERSION_CHECK=1
RUN apt-get update \
    && apt-get install -y --no-install-recommends bash ca-certificates curl git procps \
    && rm -rf /var/lib/apt/lists/*
# numpy/scipy so the agent can fit and optimise; mcp for a full MCP client if it wants one.
RUN pip install mcp==1.30.0 numpy==2.5.4 scipy==1.18.1

# Node + npm for the agent CLIs (opencode, claude-code). Without `node` and
# `npm` on PATH Harbor apt-installs them at the start of every trial (~2 min).
COPY --from=node:22-slim /usr/local/bin/node /usr/local/bin/node
COPY --from=node:22-slim /usr/local/lib/node_modules /usr/local/lib/node_modules
RUN ln -s ../lib/node_modules/npm/bin/npm-cli.js /usr/local/bin/npm \
    && ln -s ../lib/node_modules/npm/bin/npx-cli.js /usr/local/bin/npx \
    && node --version && npm --version

RUN useradd --create-home --shell /bin/bash agent \
    && mkdir -p /app && chown agent:agent /app
# labclient.py   stdlib client for the lab's tools, so tool calls can be scripted
# cellmodel.py   the public cell model (the equations wells and tissue are simulated with)
# regimen_search.py  regimen search on that model for parameters you give it
COPY --chown=agent:agent agent/ /app/
WORKDIR /app
