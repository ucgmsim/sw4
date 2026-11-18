#!/usr/bin/env python3
"""Compare the sampled (tau_exec -ebs) time of two or more TAU profile directories.

    tools/compare_tau_profiles.py gausshill-att-3-tau/9352716 gausshill-att-3-tau/9352857 \
        --labels generic hpc3-genoa

Runs "pprof -s -p" in each directory and takes the flat [SAMPLE] entries plus the MPI
calls of the total summary, so the numbers are thread-seconds summed over all ranks and
threads (sample entries are not subtracted from their enclosing timers, so the TAU
thread/region timers are left out to avoid counting time twice). Functions are grouped
into categories; the first run is the baseline for the speed-up column.
"""
import argparse
import re
import subprocess
import sys

CATEGORIES = [
    ("curvilinear4sg stencil", r"curvilinear4sg"),
    ("sin/cos (libm, svml)", r"__libm_|__svml_|\bsinf\b|\bcosf\b|sincos"),
    ("twilight forcing bodies", r"forcing|twilight|exactSol|exactAcc"),
    ("rhs4th3fort", r"rhs4th3"),
    ("OpenMP runtime wait", r"__kmp|kmp_"),
    ("memory variables", r"[mM]emVar|memvar"),
    ("MPI", r"^MPI_"),
]


def read_profile(path, pprof):
    out = subprocess.run([pprof, "-s", "-p"], cwd=path, capture_output=True,
                         text=True, check=True).stdout
    total = out.split("FUNCTION SUMMARY (total)")[1].split("FUNCTION SUMMARY (mean)")[0]
    funcs = {}
    for line in total.splitlines():
        # %Time Exclusive Inclusive #Call #Subrs Inclusive/call Name  (-p: plain msec)
        m = re.match(r"\s*[\d.]+\s+([\d.,]+)\s+[\d.,]+\s+\d+\s+\d+\s+\d+\s+(.*\S)", line)
        if not m:
            continue
        name = m.group(2)
        if name.startswith("[SAMPLE] "):
            name = name[len("[SAMPLE] "):]
        elif not name.startswith("MPI_"):
            continue
        # anonymous-namespace hashes differ between builds
        name = re.sub(r"_INTERNAL[0-9a-f]+::", "", name)
        funcs[name] = funcs.get(name, 0.0) + float(m.group(1).replace(",", "")) / 1000.0
    return funcs


def category(name):
    for cat, pat in CATEGORIES:
        if re.search(pat, name):
            return cat
    return "other"


def short(name, width=70):
    # drop the trailing argument list only, keeping template arguments like <(char)61>
    name = name.strip()
    if name.endswith(")"):
        depth = 0
        for i in range(len(name) - 1, -1, -1):
            depth += {")": 1, "(": -1}.get(name[i], 0)
            if depth == 0:
                name = name[:i]
                break
    name = re.sub(r"^void |_INTERNAL[0-9a-f]+::", "", name)
    return name if len(name) <= width else name[:width - 3] + "..."


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("dirs", nargs="+", help="TAU profile directories")
    ap.add_argument("--labels", nargs="+", help="column labels (default: directory names)")
    ap.add_argument("--top", type=int, default=15, help="number of functions to list")
    ap.add_argument("--pprof", default="pprof", help="pprof executable")
    args = ap.parse_args()
    labels = args.labels or args.dirs
    if len(labels) != len(args.dirs):
        sys.exit("need one label per directory")

    runs = [read_profile(d, args.pprof) for d in args.dirs]

    def table(title, rows):
        head = "| " + title + " | " + " | ".join(f"{l} (s)" for l in labels)
        head += "".join(f" | speed-up {l}" for l in labels[1:]) + " |"
        print(head)
        print("|" + "---|" * (1 + 2 * len(labels) - 1))
        for key, vals in rows:
            cells = " | ".join(f"{v:.1f}" for v in vals)
            ups = "".join(f" | {vals[0] / v:.2f}x" if v > 0 else " | -" for v in vals[1:])
            print(f"| {key} | {cells}{ups} |")
        print()

    cats = {}
    for i, run in enumerate(runs):
        for name, t in run.items():
            cats.setdefault(category(name), [0.0] * len(runs))[i] += t
    totals = [sum(r.values()) for r in runs]
    rows = sorted(cats.items(), key=lambda kv: -kv[1][0]) + [("**total sampled**", totals)]
    table("category", rows)

    names = sorted(set().union(*runs), key=lambda n: -max(r.get(n, 0.0) for r in runs))
    table("function", [(short(n), [r.get(n, 0.0) for r in runs]) for n in names[:args.top]])


if __name__ == "__main__":
    main()
