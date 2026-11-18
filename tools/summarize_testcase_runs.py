#!/usr/bin/env python3
"""
Collect the per-build results of several tools/run_testcase.sl jobs from their Slurm
logs (testcase-<jobid>.out) into one CSV and a grouped bar plot of the wall-clock
times: one panel per precision, one bar per (build, job).

The logs are read rather than the <name>-<jobid>.csv files so that a job that is
still running contributes the runs it has finished.

Runs that used Intel MPI's mpirun with I_MPI_PIN_DOMAIN=omp but without
I_MPI_PIN_CELL=core (every job before tools/run_testcase.sl set it) had both OpenMP
threads of a rank on one physical core. Those bars are hatched: they are only
comparable with other Intel MPI runs of the same job.

Example:
    summarize_testcase_runs.py testcase-9305990.out testcase-9307163.out \\
        -o gausshill-att-3-all.png --csv gausshill-att-3-all.csv
"""
import argparse
import csv
import os
import re
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HEAD_RE = re.compile(r"^test case: (?P<case>\S+), node: (?P<node>[^,.]+)")
LAUNCH_RE = re.compile(r"^(?P<cmd>mpirun|srun)\b(?P<args>.*?)\s(?P<sw4>\S+/build-(?P<build>[^/]+)/bin/sw4)\s")
RESULT_RE = re.compile(
    r"^(?P<build>\S+): disp errInf=(?P<errinf>\S+) errL2=(?P<errl2>\S+), .*"
    r"time stepping=(?P<step>\S+?)s wall=(?P<wall>\S+?)s \((?P<status>\w+)\)")


def parse_log(path):
    """Yield one dict per finished run in a run_testcase.sl log."""
    job = re.search(r"(\d+)", os.path.basename(path)).group(1)
    node, launch = "?", {}
    with open(path) as f:
        text = f.read()
    # run_testcase.sl logs "pinning: ... I_MPI_PIN_CELL=core" since it started setting it
    pin_cell = "I_MPI_PIN_CELL=core" in text
    for line in text.splitlines():
        m = HEAD_RE.match(line)
        if m:
            node = m["node"]
            continue
        m = LAUNCH_RE.match(line)
        if m:
            launch[m["build"]] = m["cmd"] + m["args"]
            continue
        m = RESULT_RE.match(line)
        if m and m["wall"] != "-":   # "-": the run left no timing behind
            build = m["build"]
            cmd = launch.get(build, "")
            intel_mpi = cmd.startswith("mpirun") and "--map-by" not in cmd
            yield {
                "job": job, "node": node, "build": build,
                "platform": build.rsplit("-", 1)[0], "precision": build.rsplit("-", 1)[1],
                "launch": cmd, "half_cores": intel_mpi and not pin_cell,
                "errinf": m["errinf"], "errl2": m["errl2"],
                "step_seconds": m["step"], "wall_seconds": m["wall"], "status": m["status"],
            }


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("logs", nargs="+", help="testcase-<jobid>.out files")
    parser.add_argument("-o", "--output", default="testcase-runs.png", help="image file")
    parser.add_argument("--csv", help="also write the collected runs to this CSV")
    parser.add_argument("--title", help="figure title")
    args = parser.parse_args()

    runs = [r for path in args.logs for r in parse_log(path)]
    if args.csv:
        with open(args.csv, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(runs[0]))
            w.writeheader()
            w.writerows(runs)
        print(f"{len(runs)} runs written to {args.csv}")
    ok = [r for r in runs if r["status"] == "ok"]
    if not ok:
        sys.exit("no finished runs")

    jobs = list(dict.fromkeys(r["job"] for r in ok))
    nodes = {r["job"]: r["node"] for r in ok}
    colors = dict(zip(jobs, plt.rcParams["axes.prop_cycle"].by_key()["color"]))
    precs = [p for p in ("default", "single") if any(r["precision"] == p for r in ok)]

    fig, axes = plt.subplots(len(precs), 1, figsize=(14, 4.2 * len(precs)), squeeze=False)
    for ax, prec in zip(axes[:, 0], precs):
        sel = [r for r in ok if r["precision"] == prec]
        platforms = list(dict.fromkeys(r["platform"] for r in sel))
        width = 0.8 / len(jobs)
        for i, p in enumerate(platforms):
            here = [r for r in sel if r["platform"] == p]
            for j, r in enumerate(here):
                x = i + (j - (len(here) - 1) / 2) * width
                wall = float(r["wall_seconds"])
                ax.bar(x, wall, width, color=colors[r["job"]],
                       hatch="///" if r["half_cores"] else None,
                       edgecolor="white" if r["half_cores"] else None, linewidth=0)
                ax.text(x, wall, f"{wall:.0f}", ha="center", va="bottom", fontsize=7)
        ax.set_xticks(range(len(platforms)), platforms, rotation=20, ha="right")
        ax.set_ylabel(f"wall-clock time (s), {prec}")
        ax.grid(axis="y", alpha=0.3)
        ax.set_axisbelow(True)
        ax.margins(y=0.1)

    handles = [plt.Rectangle((0, 0), 1, 1, color=colors[j]) for j in jobs]
    labels = [f"job {j} ({nodes[j]})" for j in jobs]
    handles.append(plt.Rectangle((0, 0), 1, 1, facecolor="grey", hatch="///", edgecolor="white"))
    labels.append("Intel MPI, 1 core per rank (not comparable)")
    axes[0, 0].legend(handles, labels, fontsize=8, ncol=2, loc="upper right")
    fig.suptitle(args.title or "SW4 test case timings, 4 MPI ranks x 2 OpenMP threads", fontsize=11)
    fig.tight_layout()
    fig.savefig(args.output, dpi=150)
    print(f"bar plot written to {args.output}")


if __name__ == "__main__":
    main()
