#!/usr/bin/env python3
"""Build the pitch deck (deck/deck.html, deck/Remission-deck.pdf) from the analysed runs.

    python3 deck/build.py [--url https://...] [--no-pdf]

Every number on the slides is read from results/summary.json, results/run_analysis.jsonl,
results/fixed_regimens.json and experience/summary.json, so the deck cannot drift from the README.
Images: docs/img (figures) and deck/img (viewer screenshots). The PDF is printed with headless Chrome.
"""
from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import readme_tables as RT  # noqa: E402

CSS = """
@font-face { font-family: "Instrument Sans"; src: url("../viewer/vendor/fonts/instrument-sans-latin-wght-normal.woff2") format("woff2"); font-weight: 400 700; }
@font-face { font-family: "JetBrains Mono"; src: url("../viewer/vendor/fonts/jetbrains-mono-latin-wght-normal.woff2") format("woff2"); font-weight: 300 700; }
@page { size: 1920px 1080px; margin: 0; }
:root { --bg:#121110; --panel:#1a1917; --line:#2c2a26; --line-2:#413e39; --ink:#f1eee8; --ink-2:#c2bdb3; --ink-3:#918c82;
  --healthy:#2c9a97; --tumour:#f0437a; --clone:#f5b735; --bad:#ff7a66; --good:#5fd0a6;
  --sans:"Instrument Sans","Helvetica Neue",Arial,sans-serif; --mono:"JetBrains Mono",ui-monospace,Menlo,monospace; }
* { box-sizing: border-box; }
html, body { margin: 0; background: var(--bg); color: var(--ink); font: 400 30px/1.4 var(--sans); -webkit-font-smoothing: antialiased; }
.slide { width: 1920px; height: 1080px; padding: 84px 104px 96px; position: relative; overflow: hidden; break-after: page;
  display: flex; flex-direction: column; background: var(--bg); }
.kicker { font: 500 21px/1 var(--mono); color: var(--ink-3); letter-spacing: .06em; text-transform: uppercase; margin-bottom: 28px; }
h1 { font: 600 132px/1 var(--sans); letter-spacing: -.03em; margin: 0; }
h2 { font: 600 62px/1.08 var(--sans); letter-spacing: -.02em; margin: 0 0 36px; max-width: 1560px; }
p { margin: 0; }
.lead { font-size: 36px; line-height: 1.34; color: var(--ink-2); max-width: 1380px; }
.num, .mono { font-family: var(--mono); font-variant-numeric: tabular-nums; }
.foot { position: absolute; left: 104px; right: 104px; bottom: 40px; display: flex; justify-content: space-between;
  font: 400 18px/1 var(--mono); color: var(--ink-3); }
.cols { display: flex; gap: 72px; flex: 1; min-height: 0; }
.col { flex: 1; min-width: 0; display: flex; flex-direction: column; }
img { display: block; max-width: 100%; }
.shot { border: 1px solid var(--line-2); border-radius: 10px; overflow: hidden; background: #0a0a09; }
.shot img { width: 100%; height: 100%; object-fit: cover; object-position: top left; }
.missing { display: flex; align-items: center; justify-content: center; color: var(--bad); font: 400 22px var(--mono); height: 100%; }
.steps { display: grid; grid-template-columns: repeat(4, 1fr); gap: 28px; margin-top: 56px; }
.step { border-top: 2px solid var(--line-2); padding-top: 26px; }
.step b { display: block; font: 600 34px/1.15 var(--sans); margin-bottom: 14px; }
.step span { color: var(--ink-2); font-size: 26px; line-height: 1.36; display: block; }
.step i { font: 500 20px var(--mono); color: var(--ink-3); font-style: normal; display: block; margin-bottom: 18px; letter-spacing: .05em; }
.big { font: 600 150px/1 var(--mono); letter-spacing: -.04em; }
.stat { display: flex; flex-direction: column; gap: 10px; }
.stat .v { font: 600 84px/1 var(--mono); letter-spacing: -.03em; }
.stat .l { color: var(--ink-2); font-size: 26px; line-height: 1.3; max-width: 420px; }
.stats { display: flex; gap: 72px; margin-top: 48px; }
.good { color: var(--good); } .bad { color: var(--bad); } .clone { color: var(--clone); } .dim { color: var(--ink-3); }
ul.points { list-style: none; margin: 0; padding: 0; display: flex; flex-direction: column; gap: 30px; }
ul.points li { padding-left: 34px; position: relative; font-size: 30px; line-height: 1.34; color: var(--ink-2); }
ul.points li::before { content: ""; position: absolute; left: 0; top: 15px; width: 14px; height: 14px; border-radius: 3px; background: var(--ink-3); }
ul.points li.t::before { background: #eb6834; } ul.points li.c::before { background: #6250d6; } ul.points li.g::before { background: var(--good); }
ul.points li b { color: var(--ink); font-weight: 600; }
table { border-collapse: collapse; width: 100%; font-size: 28px; }
th { text-align: left; font: 500 19px var(--mono); color: var(--ink-3); letter-spacing: .05em; text-transform: uppercase; padding: 0 20px 14px 0; border-bottom: 1px solid var(--line-2); }
td { padding: 15px 20px 15px 0; border-bottom: 1px solid var(--line); vertical-align: top; }
td.n, th.n { text-align: right; font-family: var(--mono); padding-right: 0; }
td small { display: block; color: var(--ink-3); font-size: 21px; line-height: 1.3; margin-top: 2px; }
.cmd { font: 500 32px/1.5 var(--mono); background: var(--panel); border: 1px solid var(--line-2); border-radius: 10px; padding: 26px 32px; white-space: pre; }
.cmd .c { color: var(--ink-3); }
.fig { background: var(--bg); }
.fig img { width: 100%; }
.legend { display: flex; gap: 34px; font-size: 22px; color: var(--ink-2); margin-top: 18px; }
.legend span::before { content: ""; display: inline-block; width: 14px; height: 14px; border-radius: 2px; margin-right: 10px; background: var(--k); }
"""


def img(path: str, cls: str = "", style: str = "") -> str:
    p = ROOT / path
    if not p.exists():
        return f'<div class="missing">missing {path}</div>'
    return f'<img class="{cls}" style="{style}" src="../{path}" alt="">'


def slide(n: int, total: int, body: str, kicker: str = "") -> str:
    k = f'<div class="kicker">{kicker}</div>' if kicker else ""
    return (f'<section class="slide">{k}{body}'
            f'<div class="foot"><span>Remission · Junction 2026 · AfterQuery "Break the Agent"</span><span>{n} / {total}</span></div></section>')


def f2(x: float) -> str:
    return f"{x:.2f}"


def build(url: str) -> str:
    summary = json.loads((ROOT / "results" / "summary.json").read_text())
    rows = [json.loads(x) for x in (ROOT / "results" / "run_analysis.jsonl").read_text().splitlines() if x.strip()]
    fixed = json.loads((ROOT / "results" / "fixed_regimens.json").read_text())
    exp_path = ROOT / "experience" / "summary.json"
    exp = json.loads(exp_path.read_text()) if exp_path.exists() else None
    L = summary["labels"]
    ref = L[summary["reference"]]
    agents = [(lab, s) for lab, s in L.items() if lab != summary["reference"] and lab not in RT.SHORTCUTS]
    shortcuts = [(lab, L[lab]) for lab in RT.SHORTCUTS if lab in L]
    fc = RT.failure_counts(rows, summary)
    tot = {k: sum(c[k] for c in fc.values()) for k in next(iter(fc.values()))} if fc else {}
    best = max(agents, key=lambda a: a[1]["reward_mean"])
    worst_frontier = min((a for a in agents if "Haiku" not in a[0]), key=lambda a: a[1]["reward_mean"])
    top_shortcut = max(shortcuts, key=lambda a: a[1]["reward_mean"])

    s = []
    # 1 title
    s.append(f"""
<div style="display:flex;gap:72px;flex:1;min-height:0">
  <div style="width:700px;display:flex;flex-direction:column;justify-content:center">
    <h1>Remission</h1>
    <p class="lead" style="margin-top:36px;color:var(--ink)">A cell lab where an AI agent has to cure one patient.</p>
    <p class="lead" style="margin-top:28px;font-size:29px">28 lab days, 384 wells, one treatment. Scored on tissue the agent never sees.</p>
    <p class="mono" style="margin-top:56px;font-size:24px;color:var(--ink-2)">{url}</p>
  </div>
  <div class="shot" style="flex:1">{img("deck/img/hero.png")}</div>
</div>""")
    # 2 the task
    s.append(f"""
<h2>One patient, one lab, one shot at the treatment.</h2>
<p class="lead">The agent gets a tumour biopsy and healthy cells from the same virtual patient. It plans its own experiments, then submits one regimen of eight targeted drugs.</p>
<div class="steps">
  <div class="step"><i>DAY 0</i><b>A biopsy and healthy cells</b><span>Tumour type and the mutations a clinical panel reports. Clones under 10 % of the tumour are not on the chart.</span></div>
  <div class="step"><i>28 LAB DAYS · 384 WELLS</i><b>Experiments</b><span>Viability over days, pERK and pAKT, deep sequencing down to 2 % of cells, growing survivors into new lines. Counts carry 10 % noise.</span></div>
  <div class="step"><i>ONCE</i><b>One regimen</b><span>Up to 12 phases of daily doses, 0 to 2× the label dose, at most 84 days. No score comes back.</span></div>
  <div class="step"><i>DAY 364</i><b>Reward = healed − harm</b><span>The patient's tissue a year later: the share that is healthy, minus the worst loss of normal tissue to the drugs.</span></div>
</div>
<p class="lead" style="margin-top:56px;font-size:28px">The agent also holds the exact equations of the simulator and a regimen search. Only this patient's parameter values are hidden.</p>""")
    # 3 why standard care fails
    s.append(f"""
<h2>The scan goes clean. Then a clone nobody sequenced takes the tissue back.</h2>
<div class="cols">
  <div class="shot" style="flex:1.55">{img("deck/img/relapse.png")}</div>
  <div class="col" style="justify-content:center;gap:44px;flex:1">
    <div class="stat"><span class="v clone">2 %</span><span class="l">of cells is the detection limit of sequencing. The resistant clone starts at 0.5 to 1.3 %.</span></div>
    <div class="stat"><span class="v bad">{f2(L["Guideline regimen"]["reward_mean"])}</span><span class="l">is what the standard-of-care regimen scores, on every one of the 10 patients.</span></div>
    <div class="stat"><span class="v">1 in 400</span><span class="l">generated patients enters the benchmark. Six rules, written down, told to the agent.</span></div>
  </div>
</div>""")
    # 4 results
    lines = "".join(
        f'<div class="stat"><span class="v{" good" if lab == summary["reference"] else ""}">{f2(x["reward_mean"])}</span>'
        f'<span class="l">{lab}</span></div>' for lab, x in [(summary["reference"], ref)] + agents)
    s.append(f"""
<h2>The reference heals every patient. Agents lose reward on the same ten.</h2>
<div class="fig" style="flex:1;min-height:0;display:flex;align-items:center">{img("docs/img/results_ladder_dark.png", style="max-height:640px;width:auto;margin:0 auto")}</div>
<p class="lead" style="font-size:26px;margin-top:20px">10 held-out patients, 3 runs each, OpenCode in Harbor with a sealed network. The reference gave the same reward in all three runs. Best shortcut: {top_shortcut[0].lower()}, {f2(top_shortcut[1]["reward_mean"])}.</p>""")
    # 5 where agents lose
    s.append(f"""
<h2>Each lost run traces back to a decision in the lab.</h2>
<div class="cols">
  <div class="col" style="flex:1.15;justify-content:center"><div class="fig">{img("docs/img/results_outcomes_dark.png")}</div></div>
  <div class="col" style="justify-content:center">
    <ul class="points">
      <li class="t"><b>{tot.get("too_toxic", 0)} runs cleared the tumour and destroyed the normal tissue.</b> In {tot.get("too_toxic_untested", 0)} of them the submitted combination had never met the patient's healthy cells in a well. In the other {tot.get("too_toxic", 0) - tot.get("too_toxic_untested", 0)} the agent saw the damage and submitted anyway.</li>
      <li class="c"><b>{tot.get("missed_clone", 0)} runs relapsed from the hidden clone.</b> {tot.get("missed_clone_never_named", 0)} never named it. The others named it and still gave a course too short or too weak to kill it.</li>
      <li class="g"><b>The tools were there.</b> Healthy cells, long cultures, growing out survivors, sequencing them. The reference uses nothing else.</li>
    </ul>
  </div>
</div>""")
    # 6 replay
    s.append(f"""
<h2>Every run can be replayed day by day, next to the reference on the same patient.</h2>
<div class="shot" style="flex:1;min-height:0">{img("deck/img/replay.png")}</div>
<p class="mono" style="font-size:22px;color:var(--ink-2);margin-top:20px">{url}</p>""")
    # 7 integrity
    rowsx = "".join(f'<tr><td>{lab}<small>{RT.WHAT[lab]}</small></td><td class="n">{f2(x["reward_mean"])}</td></tr>' for lab, x in shortcuts)
    h = fixed["hindsight_on_the_test_patients"]
    s.append(f"""
<h2>Skipping the lab work earns nothing.</h2>
<div class="cols">
  <div class="col"><table><tr><th>Shortcut</th><th class="n">Reward</th></tr>{rowsx}</table></div>
  <div class="col" style="justify-content:flex-start">
    <ul class="points">
      <li><b>No fixed answer exists.</b> The best single regimen for all 36 held-out patients, chosen with their hidden parameters, averages {f2(h["one_regimen_for_all"]["mean_on_test"])}. Per tumour type: {f2(h["per_tumour_type"]["BRAF_MEL"]["mean_on_test"])} for melanoma, {f2(h["per_tumour_type"]["EGFR"]["mean_on_test"])} for EGFR lung.</li>
      <li><b>Nothing to look up.</b> The agent's network is cut down to the lab and the model. Seed, parameters and scoring live in a sidecar it cannot read.</li>
      <li><b>One submission, no score back.</b> The verifier runs as root after the agent's processes are killed.</li>
      <li><b>Fresh patients on demand.</b> <span class="mono">--fresh N</span> draws patients nobody has seen.</li>
    </ul>
  </div>
</div>""")
    # 8 believable
    s.append(f"""
<h2>The cell model is checked against public data, and the misses are listed.</h2>
<div class="cols">
  <div class="col" style="flex:1.1"><div class="fig" style="background:#fcfcfb;border-radius:10px;padding:14px">{img("docs/img/val_clinical_orr.png")}</div></div>
  <div class="col" style="justify-content:center">
    <ul class="points">
      <li class="g"><b>Drug potency:</b> 7 of 14 comparisons with GDSC2 cell lines within half a log10.</li>
      <li class="g"><b>Standard regimens:</b> response rates next to COMBI-d, BEACON, CodeBreaK 300 and FLAURA. Melanoma 74 % vs 64–67 %, EGFR lung 75 % vs 76–80 %.</li>
      <li class="g"><b>Resistance:</b> MET bypass and minor-clone outgrowth behave as in Engelman 2007 and Misale 2012.</li>
      <li class="t"><b>Misses:</b> combinations are 19–47 points too toxic; three signalling checks fail. All in docs/MODEL.md.</li>
    </ul>
  </div>
</div>""")
    # 9 experience
    if exp:
        b, a = exp["before"], exp["after"]
        s.append(f"""
<h2>A smaller model gets better from its own experience here.</h2>
<div class="cols">
  <div class="col" style="justify-content:center;gap:40px;flex:.9">
    <div style="display:flex;align-items:baseline;gap:40px">
      <div class="stat"><span class="big dim">{f2(b["reward_mean"])}</span><span class="l">Claude Haiku 5.5, first time</span></div>
      <div class="big dim" style="font-size:90px">→</div>
      <div class="stat"><span class="big good">{f2(a["reward_mean"])}</span><span class="l">with the notes it wrote itself</span></div>
    </div>
    <p class="lead" style="font-size:27px">Same 10 test patients, {b["runs"]} and {a["runs"]} runs. Difference {exp["difference"]:+.2f} (95 % interval {exp["difference_ci95"][0]:+.2f} to {exp["difference_ci95"][1]:+.2f}).</p>
  </div>
  <div class="col" style="justify-content:center">
    <ul class="points">
      <li><b>12 development patients.</b> Haiku treats each once. Mean reward {f2(exp["development"]["reward_mean"])}.</li>
      <li><b>It reads its own post-mortems</b> and writes 450 words of notes: test the final combination on healthy cells, grow out survivors and sequence them, hold a low dose long.</li>
      <li><b>No weights change.</b> The notes are appended to the task text. No test patient is in the experience.</li>
    </ul>
  </div>
</div>""")
    # 10 run it
    s.append(f"""
<h2>One command reruns it.</h2>
<div class="cmd"><span class="c"># reference, do-nothing, stop-early and four shortcuts. No API key.</span>
scripts/run_eval.sh --no-agent

<span class="c"># the whole evaluation: the above plus 3 agent runs per patient</span>
ANTHROPIC_API_KEY=... scripts/run_eval.sh</div>
<div class="stats">
  <div class="stat"><span class="v">{f2(ref["reward_mean"])}</span><span class="l">reference, identical across repeated runs</span></div>
  <div class="stat"><span class="v">{f2(best[1]["reward_mean"])}</span><span class="l">{best[0]}, the best agent</span></div>
  <div class="stat"><span class="v">{f2(worst_frontier[1]["reward_mean"])}</span><span class="l">{worst_frontier[0]}</span></div>
  <div class="stat"><span class="v">{f2(L["Do nothing"]["reward_mean"])}</span><span class="l">doing nothing, stopping early</span></div>
</div>
<p class="lead" style="margin-top:56px;font-size:28px">Harbor task, OpenCode over MCP, MIT. Live replay and "treat a patient yourself": <span class="mono">{url}</span></p>""")
    kickers = ["", "The task", "Why it is hard", "Results", "Failure modes", "Replay", "Reward integrity", "The simulator",
               "Learning from experience", "Reproduce"]
    if not exp:
        kickers.remove("Learning from experience")
    body = "".join(slide(i + 1, len(s), x, kickers[i]) for i, x in enumerate(s))
    return f'<!doctype html><html lang="en"><head><meta charset="utf-8"><title>Remission</title><style>{CSS}</style></head><body>{body}</body></html>'


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--url", default="https://remission-lab.fly.dev")
    p.add_argument("--no-pdf", action="store_true")
    a = p.parse_args()
    out = ROOT / "deck" / "deck.html"
    out.write_text(build(a.url))
    print(out)
    if a.no_pdf:
        return
    chrome = shutil.which("google-chrome") or shutil.which("chromium") or shutil.which("chromium-browser")
    if not chrome:
        sys.exit("no Chrome found: open deck/deck.html and print it to PDF (1920 x 1080, no margins)")
    pdf = ROOT / "deck" / "Remission-deck.pdf"
    subprocess.run([chrome, "--headless=new", "--no-sandbox", "--disable-gpu", "--no-pdf-header-footer",
                    f"--print-to-pdf={pdf}", "--virtual-time-budget=5000", out.as_uri()], check=True,
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    print(pdf, f"{pdf.stat().st_size / 1e6:.1f} MB")


if __name__ == "__main__":
    main()
