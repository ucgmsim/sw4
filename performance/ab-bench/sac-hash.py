#!/usr/bin/env python3
"""Hash the waveform data of every SAC station file under a directory.

Production benchmark cases have no analytic solution, so the recorded
waveforms are the correctness oracle. The 632-byte SAC header is skipped
deliberately: it carries an NZSEC reference-second stamped from wall-clock at
write time, so byte-comparing whole files reports a spurious one-byte
difference on every run (verified -- offset 297).
"""
import glob, hashlib, os, sys

if len(sys.argv) < 2:
    sys.exit(0)
h, n = hashlib.sha1(), 0
for f in sorted(glob.glob(os.path.join(sys.argv[1], "**", "st*.[xyz]"), recursive=True)):
    try:
        with open(f, "rb") as fh:
            fh.seek(632)
            h.update(fh.read())
            n += 1
    except OSError:
        pass
print(f"sac{n}:{h.hexdigest()[:16]}" if n else "")
