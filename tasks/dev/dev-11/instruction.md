# Remission: find the regimen that cures this patient

You run a precision-oncology lab for one patient. The standard regimens for this tumour type do not cure
this patient; at least one regimen of the eight available drugs does. Find it.

You have a biopsy of the tumour and matched healthy cells from the same patient, **28 lab days** and
**384 wells**. Run experiments on the cells, work out how they respond, then submit **one** treatment regimen.
After the session the patient is treated with it and followed for a year.

**Reward = max(0, healed − harm)**, between 0 and 1:

- `healed`: share of the lesion's tissue that is healthy one year after the start. A tumour clone that is not
  driven extinct grows back and takes the tissue over, so a regimen that only holds the tumour scores 0.
- `harm`: the largest loss of the patient's normal tissue during the course (drug toxicity).

No regimen submitted scores 0. You get no score back; the regimen is scored on the patient's tissue, which you
never see.

## What you have

- The MCP server `remission` is the lab. Call `get_brief` first: it holds the rules, the drugs, the assays, the
  budget and the regimen format. The other tools: `patient_chart`, `lab_status`, `start_experiment`,
  `derive_line`, `advance_time`, `get_results`, `submit_regimen`.
- `/app/cellmodel.py`: the exact equations the wells and the patient's tissue are simulated with, and what is
  known about patients in general (typical drug potencies and their spread, parameter ranges, what each
  alteration changes). This patient's parameter values are hidden; the experiments are how you learn them.
- `/app/regimen_search.py`: searches regimens for parameter values you give it.
- `/app/labclient.py`: the same lab tools from Python (`from labclient import call`), for scripts.
- Python 3.12 with numpy and scipy, 2 CPUs. Work in `/app`.

## Rules

- Lab time passes only when you call `advance_time`. Experiments take lab days; plan them so the results arrive
  before day 28.
- `submit_regimen` is final and can be called once. Submit within **60 minutes** of wall-clock time; when the
  session ends without a submitted regimen the patient goes untreated (reward 0).
- Everything happens through the lab tools. There is nothing to install and nothing to look up.
