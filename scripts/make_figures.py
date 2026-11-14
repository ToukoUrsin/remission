#!/usr/bin/env python3
"""Draw the two result figures of the README from the analysed runs.

    scripts/make_figures.py [--analysis results/run_analysis.jsonl] [--summary results/summary.json] [--out docs/img]
                            [--dark]   (dark surface, files end in _dark.png: for the deck)

    docs/img/results_ladder.png     mean reward per policy with its 95 % interval and the per-patient means
    docs/img/results_outcomes.png   how each agent's runs ended (the six classes of scripts/analyze_runs.py)

Needs matplotlib (science/.venv has it). Colours: one hue per group (reference, agents, shortcuts), never per
rank; text stays in ink.
"""
from __future__ import annotations

import argparse
import json
from collections import defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib import font_manager  # noqa: E402

SURFACE, INK, INK2, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#e4e3df"
REFERENCE, AGENT, SHORTCUT = "#1baf7a", "#2a78d6", "#8a8984"
SHORTCUTS = ("Hedge", "Blind guess", "Guideline regimen", "Three-day screen", "Stop early", "Do nothing")
CLASSES = [  # key, label, colour
    ("solved", "healed, near the best regimen known", "#1baf7a"),
    ("suboptimal", "healed, with avoidable harm", "#2a78d6"),
    ("too_toxic", "tumour gone, normal tissue destroyed", "#eb6834"),
    ("missed_clone", "relapse from a hidden clone", "#6250d6"),
    ("tumour_survived", "main tumour never cleared", "#e34948"),
    ("no_regimen", "no regimen submitted", "#8a8984"),
]


def font() -> str:
    for f in sorted((Path.home() / ".local/share/fonts/inter").glob("Inter-*.otf")):
        font_manager.fontManager.addfont(str(f))
    names = {f.name for f in font_manager.fontManager.ttflist}
    return "Inter" if "Inter" in names else "DejaVu Sans"


def style(ax) -> None:
    ax.set_facecolor(SURFACE)
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)
    ax.spines["bottom"].set_color(GRID)
    ax.tick_params(axis="both", length=0, colors=INK2, labelsize=10)
    ax.xaxis.grid(True, color=GRID, linewidth=0.8)
    ax.set_axisbelow(True)


def ladder(rows: list[dict], summary: dict, out: Path) -> None:
    labels = [lab for lab in summary["labels"] if " with its " not in lab]      # the experience arms have their own section
    per_task = defaultdict(lambda: defaultdict(list))
    for r in rows:
        per_task[r["label"]][r["task"]].append(r["reward"])
    fig, ax = plt.subplots(figsize=(9.2, 0.46 * len(labels) + 1.5), dpi=200)
    fig.patch.set_facecolor(SURFACE)
    style(ax)
    for i, lab in enumerate(labels):
        s = summary["labels"][lab]
        colour = REFERENCE if lab == summary["reference"] else SHORTCUT if lab in SHORTCUTS else AGENT
        y = len(labels) - 1 - i
        m, (lo, hi) = s["reward_mean"], s["reward_ci95"]
        ax.barh(y, m, height=0.46, color=colour, alpha=0.28, linewidth=0, zorder=2)
        pts = sorted(sum(v) / len(v) for v in per_task[lab].values())
        offs = [((k % 5) - 2) * 0.055 for k in range(len(pts))]
        ax.scatter(pts, [y + o for o in offs], s=13, color=colour, alpha=0.55, linewidths=0, zorder=3)
        ax.plot([lo, hi], [y, y], color=INK, linewidth=1.6, solid_capstyle="round", zorder=4)
        ax.scatter([m], [y], s=62, color=colour, edgecolors=SURFACE, linewidths=1.6, zorder=5)
        ax.text(1.0, y, f"{m:.2f}", va="center", ha="right", color=INK, fontsize=11.5, fontweight="bold")
        ax.text(1.012, y, f"{lo:.2f} to {hi:.2f}", va="center", ha="left", color=INK2, fontsize=10)
    ax.set_yticks(range(len(labels)))
    ax.set_yticklabels(list(reversed(labels)), color=INK, fontsize=11)
    ax.set_xlim(0, 1.0)
    ax.set_xticks([0, 0.2, 0.4, 0.6, 0.8])
    ax.set_ylim(-0.7, len(labels) - 0.3)
    ax.text(1.0, len(labels) - 0.42, "mean", va="bottom", ha="right", color=INK2, fontsize=10)
    ax.text(1.012, len(labels) - 0.42, "95 % interval", va="bottom", ha="left", color=INK2, fontsize=10)
    ax.set_xlabel("reward = healed − harm   (dot: mean, line: 95 % interval, small dots: the patients)", color=INK2,
                  fontsize=10, labelpad=8)
    fig.tight_layout()
    fig.savefig(out, facecolor=SURFACE, bbox_inches="tight", pad_inches=0.15)
    plt.close(fig)


def outcomes(summary: dict, out: Path) -> None:
    labels = [lab for lab in summary["labels"] if lab not in SHORTCUTS and " with its " not in lab]
    fig, ax = plt.subplots(figsize=(9.2, 0.5 * len(labels) + 1.9), dpi=200)
    fig.patch.set_facecolor(SURFACE)
    style(ax)
    ax.xaxis.grid(False)
    for i, lab in enumerate(labels):
        s = summary["labels"][lab]
        y, left, n = len(labels) - 1 - i, 0.0, s["runs"]
        for key, _, colour in CLASSES:
            k = s["classes"].get(key, 0)
            if not k:
                continue
            w = k / n
            ax.barh(y, w - 0.004, left=left + 0.002, height=0.5, color=colour, linewidth=0, zorder=2)
            if w >= 0.06:
                ax.text(left + w / 2, y, str(k), va="center", ha="center", color="white", fontsize=10, fontweight="bold")
            left += w
        ax.text(1.012, y, f"{n} runs", va="center", ha="left", color=INK2, fontsize=10)
    ax.set_yticks(range(len(labels)))
    ax.set_yticklabels(list(reversed(labels)), color=INK, fontsize=11)
    ax.set_xlim(0, 1.0)
    ax.set_xticks([0, 0.25, 0.5, 0.75, 1.0])
    ax.set_xticklabels(["0", "25 %", "50 %", "75 %", "100 %"])
    ax.set_xlabel("share of runs", color=INK2, fontsize=10, labelpad=6)
    handles = [plt.Rectangle((0, 0), 1, 1, color=c) for _, _, c in CLASSES]
    leg = ax.legend(handles, [t for _, t, _ in CLASSES], loc="upper center", bbox_to_anchor=(0.5, -0.32), ncol=2,
                    frameon=False, fontsize=10, handlelength=1.0, handleheight=1.0, columnspacing=2.4)
    for t in leg.get_texts():
        t.set_color(INK)
    fig.tight_layout()
    fig.savefig(out, facecolor=SURFACE, bbox_inches="tight", pad_inches=0.15)
    plt.close(fig)


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--analysis", default="results/run_analysis.jsonl")
    p.add_argument("--summary", default="results/summary.json")
    p.add_argument("--out", default="docs/img")
    p.add_argument("--dark", action="store_true")
    a = p.parse_args()
    if a.dark:
        global SURFACE, INK, INK2, GRID
        SURFACE, INK, INK2, GRID = "#121110", "#f1eee8", "#c2bdb3", "#2c2a26"
    plt.rcParams["font.family"] = font()
    rows = [json.loads(x) for x in Path(a.analysis).read_text().splitlines() if x.strip()]
    summary = json.loads(Path(a.summary).read_text())
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    end = "_dark.png" if a.dark else ".png"
    ladder(rows, summary, out / ("results_ladder" + end))
    outcomes(summary, out / ("results_outcomes" + end))
    print(out / ("results_ladder" + end), out / ("results_outcomes" + end))


if __name__ == "__main__":
    main()
