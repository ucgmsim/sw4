#!/usr/bin/env python3
"""Compare two directories of SAC station files and report the difference in ULP.

A hash tells you only that two runs differ, which on a single-precision branch
is nearly useless: changing vectorisation changes FMA contraction and summation
order, so bit-identity is not the right expectation. What matters is the
magnitude. This reports the worst |A-B| relative to the trace maximum, in units
of float epsilon, so "6.7 ULP over 200 timesteps" can be told apart from a real
solver difference.

The 632-byte SAC header is skipped: it carries an NZSEC reference-second
stamped from wall-clock at write time (verified -- offset 297).
"""
import glob, os, struct, sys

EPS = {4: 1.1920929e-07, 8: 2.220446049250313e-16}

def read(path):
    d = open(path, "rb").read()
    n = (len(d) - 632) // 4
    return struct.unpack("<%df" % n, d[632:632 + 4 * n])

def main(a_dir, b_dir, word=4):
    files = sorted(glob.glob(os.path.join(a_dir, "**", "st*.[xyz]"), recursive=True))
    if not files:
        print("nofiles"); return 2
    npts = ndiff = 0
    worst = 0.0
    worst_file = None
    for fa in files:
        fb = os.path.join(b_dir, os.path.relpath(fa, a_dir))
        if not os.path.exists(fb):
            print("missing:" + os.path.basename(fb)); return 2
        a, b = read(fa), read(fb)
        if len(a) != len(b):
            print("length-mismatch:" + os.path.basename(fa)); return 2
        amax = max((abs(x) for x in a), default=0.0) or 1.0
        for x, y in zip(a, b):
            npts += 1
            if x != y:
                ndiff += 1
                r = abs(x - y) / amax
                if r > worst:
                    worst, worst_file = r, os.path.basename(fa)
    if ndiff == 0:
        print(f"identical:{npts}")
    else:
        ulp = worst / EPS[word]
        print(f"differ:{ndiff}/{npts}:{ulp:.1f}ulp:{worst:.3e}:{worst_file}")
    return 0

if __name__ == "__main__":
    sys.exit(main(sys.argv[1], sys.argv[2],
                  int(sys.argv[3]) if len(sys.argv) > 3 else 4))
