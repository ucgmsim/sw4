#!/usr/bin/env python3
"""
Extract per-test wall-clock timings from saved ctest output files and print them
as an aligned table (default on the terminal), a Markdown table, or CSV.

The input files are the captured stdout of `ctest`, named
res-<platform>-<precision>.txt (e.g. res-intel-single.txt). Each ctest result line
looks like

      5/24 Test  #5: Run_attenuation/tw-topo-att-1 ......   Passed    3.30 sec

Examples:
    # selected tests, all res-*.txt files in the current directory
    ctest_timings.py Run_attenuation/tw-topo-att-1 Run_meshrefine/refine-el-1

    # every Run_ test as a Markdown table, ready to paste into a PR or issue
    ctest_timings.py --format markdown

    # every Run_ test as CSV, one column per platform-precision
    ctest_timings.py --format csv --wide -o timings.csv

    # explicit files; the Run_ prefix may be omitted
    ctest_timings.py --files res-aocc-*.txt -- meshrefine/refine-el-1
"""
import argparse
import csv
import glob
import os
import re
import sys

LINE_RE = re.compile(
    r"^\s*\d+/\d+\s+Test\s+#\d+:\s+(?P<name>\S+)\s+\.*\s*"
    r"(?P<status>.*?)\s+(?P<secs>[0-9.]+)\s+sec\s*$"
)
FILE_RE = re.compile(r"^res-(?P<platform>.+)-(?P<prec>[^-]+)\.txt$")


def parse_file(path):
    """Return {test name: (status, seconds)} for every result line in a ctest log."""
    results = {}
    with open(path, errors="replace") as f:
        for line in f:
            m = LINE_RE.match(line)
            if m:
                status = m.group("status").lstrip("*").strip()
                results[m.group("name")] = (status, float(m.group("secs")))
    return results


def label_for(path):
    """(platform, precision) from res-<platform>-<precision>.txt, else (stem, '')."""
    base = os.path.basename(path)
    m = FILE_RE.match(base)
    if m:
        return m.group("platform"), m.group("prec")
    return os.path.splitext(base)[0], ""


def cell(res, t):
    """Human-readable table entry: seconds if passed, otherwise the status."""
    if t not in res:
        return "-"
    status, secs = res[t]
    if status == "Passed":
        return "%.2f" % secs
    return "%s (%.2f)" % (status.upper(), secs) if secs else status.upper()


def print_table(runs, tests, out, markdown):
    """One row per test, one column per platform-precision, times in seconds."""
    header = ["test"] + ["-".join(x for x in (pl, pr) if x) for pl, pr, _ in runs]
    rows = [[t[4:] if t.startswith("Run_") else t] + [cell(res, t) for _, _, res in runs]
            for t in tests]
    widths = [max(len(r[i]) for r in [header] + rows) for i in range(len(header))]

    def fmt(r):
        # test names left-aligned, timings right-aligned
        parts = [r[0].ljust(widths[0])] + [c.rjust(w) for c, w in zip(r[1:], widths[1:])]
        return ("| " + " | ".join(parts) + " |") if markdown else "  ".join(parts)

    print(fmt(header), file=out)
    if markdown:
        print("|" + "|".join([":" + "-" * (widths[0] + 1)] +
                             ["-" * (w + 1) + ":" for w in widths[1:]]) + "|", file=out)
    else:
        print("  ".join("-" * w for w in widths), file=out)
    for r in rows:
        print(fmt(r), file=out)
    if not markdown:
        print("\n(wall-clock seconds per Run_ test; '-' = not in that build's output)", file=out)


def main():
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("tests", nargs="*",
                   help="test names (Run_ prefix optional); default: all Run_ tests found")
    p.add_argument("--dir", default=".",
                   help="directory searched for res-*.txt (default: current directory)")
    p.add_argument("--files", nargs="+", help="explicit ctest output files instead of --dir")
    p.add_argument("--format", choices=["table", "markdown", "csv"],
                   help="output format (default: table on stdout, csv with -o)")
    p.add_argument("--wide", action="store_true",
                   help="csv only: one row per test, one column per platform-precision")
    p.add_argument("-o", "--output", help="file to write (default: stdout)")
    args = p.parse_args()

    files = args.files or sorted(glob.glob(os.path.join(args.dir, "res-*.txt")))
    if not files:
        sys.exit("no ctest output files found")

    runs = []  # [(platform, prec, {name: (status, secs)})]
    for path in files:
        platform, prec = label_for(path)
        runs.append((platform, prec, parse_file(path)))

    if args.tests:
        tests = [t if t.startswith(("Run_", "Check_")) else "Run_" + t for t in args.tests]
    else:
        seen = {}
        for _, _, res in runs:
            for name in res:
                if name.startswith("Run_"):
                    seen.setdefault(name, None)
        tests = list(seen)

    fmt = args.format or ("csv" if args.output else "table")
    out = open(args.output, "w", newline="") if args.output else sys.stdout
    w = csv.writer(out)
    if fmt in ("table", "markdown"):
        print_table(runs, tests, out, fmt == "markdown")
    elif args.wide:
        cols = ["-".join(x for x in (pl, pr) if x) for pl, pr, _ in runs]
        w.writerow(["test"] + cols)
        for t in tests:
            # a timing from a failed run is not comparable, so leave it blank
            w.writerow([t] + ["%.2f" % res[t][1] if t in res and res[t][0] == "Passed" else ""
                              for _, _, res in runs])
    else:
        w.writerow(["platform", "precision", "test", "status", "seconds"])
        for platform, prec, res in runs:
            for t in tests:
                status, secs = res.get(t, ("missing", None))
                w.writerow([platform, prec, t, status, "" if secs is None else "%.2f" % secs])
    if args.output:
        out.close()
        print("wrote %s" % args.output, file=sys.stderr)


if __name__ == "__main__":
    main()
