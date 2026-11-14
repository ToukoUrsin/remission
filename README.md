# Remission

**A cell lab where an AI agent has to cure one patient.**

The agent gets a tumour biopsy and healthy cells from one virtual patient, 28 lab days and 384 wells. It runs
experiments on the cells, works out how they respond, and submits one treatment. The patient is then treated and
followed for a year on tissue the agent never sees. Reward = healed − harm.

The standard regimens fail on every patient here. A resistant clone too small for sequencing to see survives the
treatment and takes the tissue back. The agent has to find that clone in the lab, and then find doses that kill
the tumour without destroying the patient's normal tissue.

![The same patient, day 24: the agent's regimen on the left, the reference scientist's on the right](docs/img/hero.png)

**Live:** https://remission-lab.fly.dev plays every run below day by day next to the reference, and lets you treat a patient yourself.

Built at Junction 2026 for AfterQuery's "Break the Agent". A Harbor task, played by OpenCode over MCP.

## Results

Ten test patients, sealed network, 60 minutes per run. Every agent run is one OpenCode session in Harbor: 3 runs per
patient for the reference, Claude Fable 5.2, Claude Opus 5.5 and GPT-6.1 Sol, 2 for Claude Haiku 5.5.

<!-- results:start -->
| | Runs | Reward (95 % CI) | Below the reference by (95 % CI) | Runs that healed the patient | How the other runs ended |
|---|---|---|---|---|---|
| **Reference scientist**<br><sub>scripted scientist, public tools only</sub> | 30 | **0.64** (0.59 to 0.69) | – | 30 of 30 | – |
| **Claude Fable 5.2** | 30 | **0.59** (0.52 to 0.66) | 0.06 (0.01 to 0.11) | 30 of 30 | 9 healed with avoidable harm |
| **Claude Opus 5.5** | 30 | **0.56** (0.48 to 0.63) | 0.08 (0.03 to 0.16) | 28 of 30 | 4 healed with avoidable harm, 2 tumour gone, normal tissue destroyed |
| **GPT-6.1 Sol** | 30 | **0.15** (0.05 to 0.27) | 0.49 (0.39 to 0.59) | 8 of 30 | 16 tumour gone, normal tissue destroyed, 5 healed with avoidable harm, 4 relapse from a hidden clone, 2 main tumour never cleared |
| **Claude Haiku 5.5** | 20 | **0.28** (0.16 to 0.42) | 0.36 (0.24 to 0.47) | 8 of 20 | 10 tumour gone, normal tissue destroyed, 5 healed with avoidable harm, 1 main tumour never cleared, 1 relapse from a hidden clone |
| Hedge<br><sub>sequences the biopsy, then covers every clone that could hide</sub> | 10 | **0.05** (0.00 to 0.12) | 0.59 (0.52 to 0.66) | 1 of 10 | 5 relapse from a hidden clone, 3 tumour gone, normal tissue destroyed, 1 healed with avoidable harm, 1 main tumour never cleared |
| Blind guess<br><sub>best fixed regimen for the tumour type, no experiment</sub> | 10 | **0.04** (0.00 to 0.12) | 0.61 (0.53 to 0.68) | 1 of 10 | 6 relapse from a hidden clone, 3 tumour gone, normal tissue destroyed, 1 healed with avoidable harm |
| Guideline regimen<br><sub>standard of care for the tumour type, no experiment</sub> | 10 | **0.00** (0.00 to 0.01) | 0.64 (0.59 to 0.69) | 0 of 10 | 10 relapse from a hidden clone |
| Three-day screen<br><sub>one short drug screen, best-looking pair</sub> | 10 | **0.00** (0.00 to 0.01) | 0.64 (0.59 to 0.69) | 0 of 10 | 8 relapse from a hidden clone, 1 tumour gone, normal tissue destroyed, 1 main tumour never cleared |
| Stop early<br><sub>reads the brief and quits</sub> | 10 | **0.00** (0.00 to 0.00) | 0.64 (0.59 to 0.69) | 0 of 10 | 10 no regimen |
| Do nothing<br><sub>no tool call</sub> | 10 | **0.00** (0.00 to 0.00) | 0.64 (0.59 to 0.69) | 0 of 10 | 10 no regimen |
<!-- results:end -->

![Mean reward per policy with 95 % interval and per-patient means](docs/img/results_ladder.png)

Healed means the tumour is gone in at least half of the 16 replicate lesions and at least half of the best known
reward for that patient is left. Gaps are computed before rounding.

How to read it:

- **The reference heals every patient in every run** and scores 0.64. With a patient's true parameters and no lab
  work the same regimen search reaches 0.68 (`results/true_parameter_search.jsonl`), so little is left above it.
- **GPT-6.1 Sol and Claude Haiku 5.5 are far below it**, by 0.49 and 0.36. Sol scores 0.00 in 17 of its 30 runs.
- **Claude Opus 5.5 and Claude Fable 5.2 are close**, 0.08 and 0.06 below. Both intervals exclude 0 and both models
  are below the reference on 8 of the 10 patients, but not on every patient or in every run: 8 of Fable's 30 runs and
  6 of Opus's beat the reference on the same patient. What these two lose is normal tissue: mean harm 0.37 and 0.39
  against the reference's 0.29.

Every run, per patient (each number is one run's reward; the reference gave the same reward in all three runs):

<!-- patients:start -->
| Patient | Tumour | Clone below the sequencing limit | Reference scientist | Claude Fable 5.2 | Claude Opus 5.5 | GPT-6.1 Sol | Claude Haiku 5.5 | Best shortcut |
|---|---|---|---|---|---|---|---|---|
| test-01 | melanoma, BRAF V600E | BCL-xL amplified 1.1 % | 0.65 | 0.67 0.57 0.71 | 0.60 0.62 0.63 | 0.66 0.06 0.41 | 0.26 0.58 | 0.35 |
| test-02 | lung, EGFR L858R | MET amplified 1.0 % | 0.63 | 0.68 0.69 0.73 | 0.71 0.69 0.56 | 0.54 0.00 0.12 | 0.34 0.21 | 0.02 |
| test-03 | melanoma, BRAF V600E | MCL-1 amplified 0.9 % | 0.79 | 0.62 0.74 0.69 | 0.72 0.63 0.75 | 0.46 0.58 0.00 | 0.72 0.46 | 0.12 |
| test-04 | melanoma, BRAF V600E | BCL-xL amplified 1.3 % | 0.64 | 0.48 0.60 0.56 | 0.62 0.60 0.54 | 0.37 0.00 0.32 | 0.19 0.14 | 0.00 |
| test-05 | melanoma, BRAF V600E | BCL-xL amplified 0.9 % | 0.73 | 0.80 0.78 0.77 | 0.59 0.70 0.77 | 0.00 0.00 0.00 | 0.00 0.26 | 0.00 |
| test-06 | lung, EGFR L858R | BCL-xL amplified 1.2 % | 0.71 | 0.67 0.58 0.64 | 0.62 0.61 0.00 | 0.71 0.00 0.04 | 0.56 0.60 | 0.38 |
| test-07 | lung, EGFR L858R | MET amplified 1.0 % | 0.65 | 0.39 0.57 0.49 | 0.65 0.55 0.56 | 0.00 0.00 0.00 | 0.15 0.65 | 0.04 |
| test-08 | lung, EGFR L858R | PTEN loss 1.3 %, EGFR T790M 0.9 % | 0.54 | 0.48 0.48 0.34 | 0.56 0.38 0.45 | 0.00 0.20 0.00 | 0.07 0.00 | 0.00 |
| test-09 | colorectal, KRAS G12C | BCL-xL amplified 0.5 % | 0.51 | 0.49 0.50 0.50 | 0.53 0.54 0.47 | 0.00 0.00 0.00 | 0.01 0.42 | 0.00 |
| test-10 | colorectal, BRAF V600E | BCL-xL amplified 0.5 % | 0.60 | 0.54 0.40 0.42 | 0.25 0.44 0.53 | 0.08 0.00 0.00 | 0.00 0.08 | 0.01 |
<!-- patients:end -->

All runs are in `results/run_analysis.jsonl` (one line per run: regimen, lab work, what went wrong, what the same
drugs would have scored at other doses) and can be replayed in the viewer.

## Run it

Needs Docker and [Harbor](https://github.com/harbor-framework/harbor) 0.25.0 (`uv tool install harbor==0.25.0`).

```bash
# reference + do-nothing + stop-early + four shortcuts, 10 patients. No model, no API key. About 20 minutes.
scripts/run_eval.sh --no-agent

# everything above plus 3 agent runs per patient (default model: anthropic/claude-fable-5-2)
ANTHROPIC_API_KEY=... scripts/run_eval.sh
```

The second command is the whole evaluation: it prints the reward per patient for the reference, the baselines and
the agent. With 30 agent runs of up to an hour it takes about 2 hours at the default 4 runs in parallel (`-n 8` halves it).

Smaller pieces:

```bash
ANTHROPIC_API_KEY=... scripts/run_eval.sh tasks/test-07        # one patient: reference, baselines, 3 agent runs. 20 to 40 minutes.
OPENAI_API_KEY=... scripts/run_eval.sh -m openai/gpt-6.1-sol     # another model
ANTHROPIC_API_KEY=... scripts/run_eval.sh --fresh 3              # 3 patients nobody has seen (see below)
scripts/run_eval.sh --help
```

After a run:

```bash
scripts/analyze_runs.py "Reference scientist=jobs/<tag>-oracle" "My agent=jobs/<tag>-agent" --md -   # why each run scored what it did
python3 scripts/build_viewer.py "Reference scientist=jobs/<tag>-oracle" "My agent=jobs/<tag>-agent"
python3 viewer/serve.py                                                                              # http://localhost:8377
python3 -m pytest tests -q                                                                           # 32 reward-integrity tests, no Docker
```

The reference is deterministic: three Harbor runs per patient gave the same reward to four decimals
(`results/summary.json`). On a second machine, from a clean copy, `scripts/run_eval.sh --no-agent` printed the same
reference and baseline rewards in 18 min 35 s. Agent runs vary; the intervals above come from 3 runs on each of 10
patients, resampling patients and runs.

## What the agent has to do

One episode is one patient. The lab is an MCP server with eight tools: `get_brief`, `patient_chart`, `lab_status`,
`start_experiment`, `derive_line`, `advance_time`, `get_results`, `submit_regimen`.

| | |
|---|---|
| Drugs | 8 classes: RAF, MEK, EGFR, MET, KRAS G12C, PI3K, BCL-2/xL and MCL-1 inhibitors |
| Budget | 384 wells, 28 lab days, 60 minutes of wall clock |
| Viability assay | cells at constant drug concentrations, counted on up to 8 days, 1 well each, 10 % count error |
| Phospho assay | pERK and pAKT at one time point, 1 well each, 10 % error |
| Sequencing | alterations carried by at least 2 % of the cells; 16 wells, 2 days |
| Derive a line | grow the survivors of one finished well into a new cell source; 8 wells, 3 days |
| Regimen | up to 12 phases of daily doses, 0 to 2× the label dose, at most 84 days |
| Score | the regimen on a 64 × 64 tissue lattice, read on day 364, mean of 16 replicates |

The agent also gets, as files in its container:

- `/app/cellmodel.py`: the exact equations the wells and the tissue are simulated with, and what is known about
  patients in general (typical potencies and their spread, parameter ranges, what each alteration does).
- `/app/regimen_search.py`: searches regimens for parameter values the agent supplies.
- `/app/labclient.py`: the lab tools from Python.

What is hidden is this patient's parameter values: how potent each drug is in this patient, how the healthy cells
and the tumour respond, and which minor clones are there at what share. The experiments are the only way to learn
them. Nothing else is withheld: the brief tells the agent that the standard regimens fail, that a clone below the
detection limit decides the outcome, and that a regimen covering every possible clone without finding out which
one is there does not cure either (`brief()` in `remission/server.py`).

The reference scientist (`reference/scientist.py`, 49 tool calls, about 4 minutes) uses only these tools and files:

1. Day 0: sequence the biopsy; long viability panels on tumour and healthy cells; pERK/pAKT panels.
2. Day 10: fit the model. Wells where cells grow back although the fit says they should be dead hold a clone the
   fit does not know. Expand those survivors into lines.
3. Day 13: sequence the lines, test drugs on them.
4. Day 18: refit with the clones named, search regimens over an ensemble of fits, culture the best candidates on
   tumour and healthy cells.
5. Day 27: submit the regimen with the best expected reward among those whose culture behaved as the fit says.

## Where agents lose the reward

Every run is first sorted by what happened in the scored tissue, then tied to what the agent did in the lab
(`scripts/analyze_runs.py`). The text for each run is in `results/run_analysis.jsonl` and in the viewer.

<!-- failures:start -->
| | Runs | Healed | of those: with avoidable harm | Tumour gone, normal tissue destroyed | Relapse from a hidden clone | of those: the clone was never named in the lab | Main tumour not cleared |
|---|---|---|---|---|---|---|---|
| Claude Fable 5.2 | 30 | 30 | 9 | 0 | 0 | 0 | 0 |
| Claude Opus 5.5 | 30 | 28 | 4 | 2 | 0 | 0 | 0 |
| GPT-6.1 Sol | 30 | 8 | 5 | 16 | 4 | 1 | 2 |
| Claude Haiku 5.5 | 20 | 8 | 5 | 10 | 1 | 1 | 1 |

Harm and how much of the budget was used:

| | Mean harm | Runs that lost reward to harm: final combination never on healthy cells in a well | Median lab day at submission (of 28) | Median wells used (of 384) | Median minutes (of 60) |
|---|---|---|---|---|---|
| Reference scientist | 0.29 | – | 27 | 322.5 | 4 |
| Claude Fable 5.2 | 0.37 | 2 of 9 | 25 | 280 | 12 |
| Claude Opus 5.5 | 0.39 | 4 of 6 | 19 | 229.5 | 10 |
| GPT-6.1 Sol | 0.69 | 14 of 21 | 19 | 329.5 | 11 |
| Claude Haiku 5.5 | 0.63 | 4 of 15 | 22.5 | 225 | 5 |
<!-- failures:end -->

![How each agent's runs ended](docs/img/results_outcomes.png)

Three decisions account for the lost reward.

**1. The dose has to be tested on the patient's healthy cells.** The reward subtracts the worst loss of normal
tissue during the course. Every drug also kills healthy cells, by an amount that differs per patient and per drug,
and the only way to learn it is to put the drugs on the healthy cells in a well. GPT-6.1 Sol destroyed the normal
tissue in 16 of 30 runs; in 13 of the 16 the submitted combination had never been on the patient's healthy cells.
Claude Haiku 5.5 did it in 10 of 20 runs, and in 7 of the 10 the combination had been on healthy cells and was
submitted anyway. Claude Opus 5.5: 2 of 30 runs. Claude Fable 5.2: none, but 9 of its 30 runs carry avoidable harm.

**2. Survivors have to be grown out, sequenced, and then killed.** The clone that decides the outcome is under the
2 % sequencing limit in the biopsy. It shows up as wells that regrow under a drug that should have emptied them.
`derive_line` expands those survivors and sequencing the line names the clone. The shortcuts cannot do this: the
guideline regimen relapses on all 10 patients, the three-day screen on 8. The agents mostly can: Fable sequenced a
derived line in 30 of 30 runs, Opus in 29, Sol in 27, Haiku in 10 of 20. 4 Sol runs and 1 Haiku run relapsed from
the hidden clone.

**3. The course has to drive every clone extinct.** A tumour that is held and not eradicated regrows after the last
dose. Short courses at high dose leave survivors and still cost more normal tissue than a low dose held to day 84.
The reference, Fable and Opus treat for 84 days at the median. Sol's median course is 28 days and Haiku's 24.5,
at doses high enough to cost the tissue.

In every one of these runs the agent had what it needed. The brief says that a clone below the detection limit
decides the outcome. Healthy cells, long cultures, `derive_line` and sequencing fit in the budget: the reference
uses about 320 of the 384 wells and 27 of the 28 days. `cellmodel.py` holds the exact equations, including how harm
is computed. Nothing is hidden behind wording, and no run was lost to the clock: the longest agent run took 36 of
the 60 minutes, and the median is 5 to 12 minutes depending on the model. Opus and Sol submit on lab day 19 of 28.

## Why the reward is hard to game

Doing less than the task scores about nothing. Every shortcut below is a script in `solution/solve.py` that plays
through the same public tools, and `scripts/run_eval.sh` reruns all of them:

| Shortcut | What it does | Mean reward | Best patient |
|---|---|---|---|
| Do nothing | no tool call | 0.00 | 0.00 |
| Stop early | reads the brief and quits | 0.00 | 0.00 |
| Guideline regimen | standard of care for the tumour type, no experiment | 0.00 | 0.04 |
| Three-day screen | one short screen of every drug and pair, then the pair that kills most | 0.00 | 0.04 |
| Blind guess | the fixed regimen that is best for the tumour type on the 12 development patients, no experiment | 0.04 | 0.38 |
| Hedge | sequences the biopsy, then the regimen the public model rates best on average over every clone that could be hiding | 0.05 | 0.35 |

Blind guess and hedge are not exactly zero: each gets lucky on one patient (0.38 and 0.35) and scores 0.00 to 0.12
on the other nine.

- **The reward is the task.** `reward = max(0, healed − harm)`, read from the simulated tissue on day 364. No
  regimen is 0. A tumour that is held but not driven extinct regrows and takes the tissue, so `healed` is 0. Killing
  the tumour with everything at full dose destroys the normal tissue, so `harm` is about 1. There is no partial
  credit for tool calls, for using the budget, or for anything the agent says.
- **No fixed answer.** `scripts/fixed_regimen_check.py` scores 3112 constant regimens on all 36 held-out cohort
  patients with their hidden parameters. The single best regimen for all of them averages 0.09. Chosen per tumour
  type, still with hindsight: 0.14 for melanoma (22 patients) and 0.14 for EGFR lung (12 patients). Chosen on the
  development patients, as an agent with memory of earlier patients could: 0.04 (`results/fixed_regimens.json`).
- **Nothing to look up.** By default the agent's container can reach the lab, the model's API host and OpenCode's
  model list, nothing else (`--open` lifts this). The patient's seed, parameters and the scoring code are in the lab
  sidecar. The agent's container holds the public model, the search and the lab client.
- **One submission, no score back.** `submit_regimen` works once. Results appear only when the lab clock passes
  them. Wells, days and the regimen format are enforced by the lab, not by the prompt.
- **The verifier does not trust the agent's container.** It runs as root after the agent's processes are killed,
  takes the output directory, and fetches the score from the sidecar with a token that is uploaded only after
  the agent has finished. Any failure writes 0.
- **Tested.** `tests/test_integrity.py` (32 tests) checks the points above in-process, among them: untreated
  scores 0, everything at full dose scores under 0.05, malformed regimens are refused, nothing hidden leaks into
  any tool result, clones below the detection limit are not reported.
- **Attacked.** A separate GPT-6.1 Sol session read the repository and then tried, inside real Harbor trials, to earn
  reward without treating the patient or to read something hidden: the environment, `/proc`, mounts and the Docker
  socket, the lab's admin routes, the reward file, malformed regimens, the clock and the well budget, racing
  submissions, a second MCP session, call timing. No attack scored above 0 and none exposed a secret
  ([docs/REDTEAM.md](docs/REDTEAM.md)).

One thing a judge should know: the seeds of the ten test patients are in this repository
(`tasks/test-*/environment/docker-compose.yaml`), because Harbor builds the sidecar from the task directory. An agent
inside the task cannot read them. Someone with the repository can compute a good regimen offline. For a run where
that is impossible, `scripts/run_eval.sh --fresh N` draws N new cohort patients from a random salt on your machine
(about one core-hour of screening per patient).

## Is the cell model believable?

The model is a signalling ODE (receptor → RAS → RAF → MEK → ERK and PI3K → AKT, one feedback state), a cell-fate
rule with an apoptotic threshold that BH3 mimetics lower, one-compartment pharmacokinetics, and a tissue lattice.
[docs/MODEL.md](docs/MODEL.md) checks it against public data on 600 unselected virtual patients:

| Check | Against | Result |
|---|---|---|
| 72 h IC50 per drug class and tumour type | GDSC2 cell lines (14 comparisons) | 7 within half a log10; the 7 misses are listed |
| Half-lives | FDA and EMA labels | match |
| Response to standard regimens | COMBI-d, BEACON, CodeBreaK 300, FLAURA and others | BRAF melanoma 74 % vs 64–67 %, BRAF colorectal 29 % vs 21–27 %, KRAS 21 % vs 26–46 %, EGFR lung 75 % vs 76–80 % |
| Resistance mechanisms | Engelman 2007, Misale 2012, Lito 2012 and others | MET bypass, minor KRAS clone outgrowth and RAF-inhibitor feedback pass; ERK rebound in colorectal, paradoxical activation and MEK feedback relief fail |
| Toxicity at label dose | grade ≥ 3 adverse events in the same trials | single drugs believable; combinations 19–47 points too toxic |

Three parameters were moved towards the published data during validation, before the benchmark patients were
drawn. Nothing was changed after that and nothing was tuned on agent results. The misses matter for how to read
the benchmark and are spelled out in [docs/MODEL.md](docs/MODEL.md), section 6. The main one: full-dose
combinations cost more normal tissue in the model than in the clinic, so an agent gains more from lowering
combination doses than a clinician would. The agent's lab runs on
the same model as the scored tissue, so the healthy cells in its wells show the harm it will be scored on.

This is a simulator for testing agents. It is not a clinical tool and its regimens are not treatment advice.

## How the patients were chosen

About 1 generated patient in 400 enters the benchmark. Six rules (`remission/cohort.py`), checked on the tissue
lattice:

1. Refractory: each of the five guideline regimens for the tumour type scores below 0.10.
2. Curable: at least one of 3112 constant-dose regimens scores 0.50 or more.
3. Progressive: untreated, the tumour takes the lesion.
4. No simple answer: every single drug and every pair at label dose scores below 0.10.
5. A hidden clone decides: a resistant clone is below the detection limit of sequencing (2 % of cells), and the
   best regimens for the visible tumour score below 0.10.
6. Hedging fails: the regimen that is best on average over every clone that could be hiding scores below 0.10.

20,000 generated patients → 48 in the cohort: 36 under the test salt, 12 under the development salt
([docs/MODEL.md](docs/MODEL.md), section 5, has the funnel). The test suite is the first 10 test-salt patients
with at most 4 of one tumour type: 4 melanoma, 4 EGFR lung, 1 KRAS, 1 BRAF colorectal. The development patients
(`tasks/dev/`) are for building agents and for the experience runs below; no reported number comes from them
except where stated.

The cohort is hard on purpose and is not a random sample of patients. The selection uses the simulator itself, so
the lab's brief tells the agent what the selection was.

## Learning from experience

The challenge's tiebreaker asks whether a smaller model improves after training on experience from the environment.
We ran the loop that needs no weights (`scripts/experience.py`). Claude Haiku 5.5 treats the 12 development patients
once. It then reads the analysis of each of its own runs (the text of `scripts/analyze_runs.py`: what its lab work
found, which clone regrew, what other doses would have scored) and writes notes for itself
([experience/notes.md](experience/notes.md), 462 words). The notes are appended to the task text. No test patient is
in the experience.

The honest result: it gets better on patients it has treated before, and not measurably on new ones.

- **Patients it had treated:** the same 12 development patients go from 0.09 to 0.34 with the notes. The notes name
  those patients' clones and doses, so part of this is memory.
- **New patients:** on the ten test patients the reward goes from 0.28 to 0.30, a difference of +0.01 (95 % interval
  −0.14 to +0.16), better on 5 patients and worse on 5.
- **A second round** (it treats the development patients again with its notes, then rewrites them:
  [experience/notes_round2.md](experience/notes_round2.md)) scores 0.21 on the test patients, −0.08 (−0.27 to +0.10).

<!-- experience:start -->
| | Runs | Reward (95 % CI) | Healed | Tumour gone, normal tissue destroyed | Relapse from a hidden clone | Mean harm | Median lab day at submission | Median wells used |
|---|---|---|---|---|---|---|---|---|
| Claude Haiku 5.5 | 20 | **0.28** (0.16 to 0.42) | 8 | 10 | 1 | 0.63 | 22.5 | 225 |
| with its notes | 20 | **0.30** (0.17 to 0.43) | 10 | 6 | 4 | 0.51 | 25 | 233 |
| with its second notes | 20 | **0.21** (0.09 to 0.33) | 8 | 6 | 6 | 0.52 | 23 | 227.5 |
<!-- experience:end -->

The notes do change how it works. It destroys the normal tissue in 6 of 20 runs, down from 10, it sequences a
derived line in 18 runs, up from 10, and its median course goes from 24.5 to 84 days. The reward does not follow,
because the lower doses it now gives let the hidden clone back: relapses go from 1 to 4, and to 6 in the second
round. Fine-tuning on the trajectories was not tried.

## Limits

- **It is a simulator.** The model is calibrated against public data and misses some of it
  ([docs/MODEL.md](docs/MODEL.md), section 6). Full-dose combinations are too toxic compared with the clinic. Lesions
  are about 900 cells and time is compressed.
- **The cohort is selected with the simulator itself**, 1 generated patient in 400. It says nothing about how
  common such patients are.
- **The reference is not perfect.** On all 36 held-out cohort patients, run in-process, it averages 0.63 and
  fails one (`results/reference_all_test_salt_patients.jsonl`). On the ten test patients it heals all ten.
- **Ten patients, three runs each.** The intervals are wide and the tumour mix is 4 melanoma, 4 EGFR lung,
  1 KRAS, 1 BRAF colorectal.
- **One harness.** Agent numbers are OpenCode with its default prompt and tools in Harbor 0.25.0. No prompt was
  tuned for any model.
- **The two strongest agents are close to the reference.** Fable is 0.06 below it and Opus 0.08, with intervals of
  0.01 to 0.11 and 0.03 to 0.16, and single runs beat it. The gap that is large and holds in every comparison is
  GPT-6.1 Sol's and Haiku's.
- **On two patients the reference leaves something.** On test-01 and test-02 its regimen clears 14 of the 16
  replicate lesions (rewards 0.65 and 0.63).
- **A stricter variant did not widen the gap, so it is not the benchmark.** Branch `v4` adds a hard tolerability
  limit (harm above 0.30 scores 0) and a session clock. On five pilot patients Fable scored 0.70 against the
  reference's 0.67, using 26 to 31 of the 60 minutes and nearly all wells; Opus 0.66 on three patients; Sol 0.00 on
  all five. Two things changed at once there, so it is a pilot and no reported number comes from it.
- **Runs left out** are listed with the reason in `results/excluded_runs.jsonl`: 9 that say nothing about the agent
  (4 where the model API failed before a regimen was submitted, 5 reference trials that never reached the
  verifier), 5 we cancelled mid-run when the runs moved to another machine, and 5 valid reference runs beyond the
  first three per patient. The reported runs are the first valid ones per patient in the order they were started.
  Nothing was left out for its score.

## Repository

| Path | What |
|---|---|
| `remission/` | the cell model (`model.py`), patient generator (`patients.py`), lab (`lab.py`), MCP server (`server.py`), regimen search (`oracle.py`), cohort rules (`cohort.py`) |
| `tasks/test-01..10`, `tasks/dev/` | Harbor tasks: one patient each (agent container, lab sidecar, verifier, reference solution) |
| `reference/scientist.py` | the reference scientist |
| `scripts/run_eval.sh` | the evaluation: reference, baselines, agent runs, results table |
| `scripts/analyze_runs.py` | explains every run: what the lab work found, which clone regrew, what other doses would have scored |
| `scripts/new_task.py` | makes task directories; `--fresh` draws new cohort patients |
| `scripts/fixed_regimen_check.py`, `scripts/experience.py`, `scripts/make_figures.py` | the fixed-regimen check, the experience loop, the figures |
| `viewer/` | the replay viewer and "treat a patient" (`viewer/serve.py`) |
| `results/` | every reported number: `summary.json`, `run_analysis.jsonl`, `fixed_regimens.json` |
| `science/`, `docs/MODEL.md` | public data, sources and the validation of the model |
| `docs/REDTEAM.md` | what an attacker tried against the reward, and what happened |
| `deck/` | the pitch deck (`Remission-deck.pdf`), built from `results/` by `deck/build.py` |
| `tests/` | reward-integrity tests |

MIT licence. All code was written during Junction 2026.
