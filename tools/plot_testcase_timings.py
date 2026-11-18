#!/usr/bin/env python3
"""
Bar plot of the execution times of one test case across (platform, precision) builds,
from the CSV written by tools/run_testcase.sl:

    platform,precision,disp_errinf,disp_errl2,att_errinf,att_errl2,step_seconds,wall_seconds,status

Bars are grouped by platform, one bar per precision. Builds that are missing or
failed are left out.

Examples:
    plot_testcase_timings.py gausshill-att-3-1234567.csv
    plot_testcase_timings.py refine-el-2-*.csv --metric step -o refine-el-2-step.png
"""
import argparse
import csv
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

METRICS = {
    "wall": ("wall_seconds", "wall-clock time (s)"),
    "step": ("step_seconds", "time stepping phase (s)"),
}


def read_times(path, column):
    """Return {(platform, precision): seconds} for the runs that finished."""
    times = {}
    with open(path, newline="") as f:
        for row in csv.DictReader(f):
            if row["status"] == "ok" and row[column]:
                times[(row["platform"], row["precision"])] = float(row[column])
    return times


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("csvfile", help="CSV written by tools/run_testcase.sl")
    parser.add_argument("-o", "--output", help="image file (default: <csvfile>.png)")
    parser.add_argument("--metric", choices=METRICS, default="wall",
                        help="which time to plot (default: wall)")
    parser.add_argument("--title", help="plot title (default: the CSV file name)")
    args = parser.parse_args()

    column, ylabel = METRICS[args.metric]
    times = read_times(args.csvfile, column)
    if not times:
        sys.exit(f"no finished runs with {column} in {args.csvfile}")

    # keep the order in which the runs appear in the CSV
    platforms = list(dict.fromkeys(p for p, _ in times))
    precs = list(dict.fromkeys(q for _, q in times))
    width = 0.8 / len(precs)

    # centre each platform's bars on its tick, however many precisions it has
    xpos = {}
    for i, p in enumerate(platforms):
        present = [q for q in precs if (p, q) in times]
        for j, q in enumerate(present):
            xpos[(p, q)] = i + (j - (len(present) - 1) / 2) * width

    fig, ax = plt.subplots(figsize=(1.8 * len(platforms) + 2, 4.5))
    for prec in precs:
        keys = [(p, prec) for p in platforms if (p, prec) in times]
        bars = ax.bar([xpos[k] for k in keys], [times[k] for k in keys], width, label=prec)
        ax.bar_label(bars, fmt="%.1f", padding=2, fontsize=9)

    ax.set_xticks(range(len(platforms)), platforms)
    ax.set_ylabel(ylabel)
    ax.set_title(args.title or os.path.basename(args.csvfile), fontsize=10)
    ax.legend(title="precision")
    ax.margins(y=0.12)
    ax.grid(axis="y", alpha=0.3)
    ax.set_axisbelow(True)
    fig.tight_layout()

    output = args.output or os.path.splitext(args.csvfile)[0] + ".png"
    fig.savefig(output, dpi=150)
    print(f"bar plot written to {output}")


if __name__ == "__main__":
    main()
