# Sidecar image: the environment server and the hidden patient. Never used for
# the agent's container. Seed and admin token are NOT baked in; compose injects
# them as env vars for this service only.
# Build context: this Dockerfile, requirements.txt and the `remission` package
# (scripts/new_task.py assembles it as <task>/environment/sidecar).
FROM python:3.12-slim

ENV PYTHONUNBUFFERED=1 PIP_NO_CACHE_DIR=1 PIP_DISABLE_PIP_VERSION_CHECK=1
WORKDIR /srv

COPY requirements.txt .
RUN pip install -r requirements.txt

COPY remission ./remission

EXPOSE 8000
CMD ["python", "-m", "remission.server"]
