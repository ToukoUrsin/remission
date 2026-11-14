# Remission viewer + "Treat a patient" sandbox. Built by deploy/deploy.sh from a staged context that holds only
# viewer/, remission/, reference/fixed_regimens.json and tasks/suites.json.
FROM python:3.12-slim
ENV PYTHONDONTWRITEBYTECODE=1 PYTHONUNBUFFERED=1 PIP_NO_CACHE_DIR=1 PIP_DISABLE_PIP_VERSION_CHECK=1
RUN pip install "numpy==2.5.4" && useradd --create-home --uid 10001 app
WORKDIR /app
COPY remission/ remission/
COPY reference/fixed_regimens.json reference/fixed_regimens.json
COPY tasks/suites.json tasks/suites.json
COPY viewer/ viewer/
USER app
EXPOSE 8080
CMD ["python3", "viewer/serve.py", "--host", "0.0.0.0", "--port", "8080"]
