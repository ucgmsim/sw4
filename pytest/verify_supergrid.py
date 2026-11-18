"""Checks for the supergrid absorbing-layer guards.

Both test harnesses compare only named result files, and every assertion here is
about what sw4 printed, so this module reads the captured stdout instead. It is
wired into test_sw4.py's dispatch on the 'supergrid.log' result name.

The cases come in pairs on purpose: sg-source-1 aborts and sg-source-2 is the
same geometry with the developer override, sg-margin-1 aborts at W+2h and
sg-margin-2 passes at W+6h. A single case could be satisfied by a predicate that
is always true or always false; the pairs cannot.
"""

import math
import os

# Vs of the single block in every supergrid test input, m/s.
TEST_VS = 3464.0
# Supergrid width in those inputs: gp=10 at h=100.
TEST_WIDTH = 1000.0
# (2772/1024)/(2*pi): the maximum slope of SW4's own stretching polynomial,
# Psi0'(xi) = 2772 xi^5 (1-xi)^5, divided by 2*pi. Written out here rather than
# imported so an edit to Psi0 in SuperGrid.C makes this test fail loudly.
C_ADIAB = (2772.0 / 1024.0) / (2.0 * math.pi)


def _read(pytest_dir, case_dir):
    path = os.path.join(pytest_dir, "supergrid", case_dir + ".out")
    with open(path) as f:
        return f.read()


def _require(cond, msg):
    if not cond:
        print("ERROR: verify_supergrid:", msg)
    return cond


def _no_error_word(text):
    # pytest/test_sw4.py fails any hdf5 test whose stdout contains the substring
    # 'Error'. The supergrid block lands in that stdout, so pin the rule here
    # too: no message this feature adds may use that word.
    bad = [l for l in text.splitlines()
           if "Error" in l and "Errorcode" not in l]
    return _require(not bad,
                    "stdout contains 'Error': %r" % (bad[:3],))


def _count(text, needle):
    return sum(1 for l in text.splitlines() if needle in l)


def verify(pytest_dir, case_dir):
    text = _read(pytest_dir, case_dir)
    ok = True

    if case_dir == "sg-source-1":
        # Source 200 m inside the x=0 layer: abort, naming the face and the
        # penetration, and offering the override.
        ok &= _require("FATAL:" in text, "no FATAL line")
        ok &= _require("INSIDE the layer at face x=0" in text,
                       "did not name the x=0 face")
        ok &= _require("allowsourceinsupergrid" in text,
                       "did not mention the override")
        # A source in the sponge must never reach the solver.
        ok &= _require("Supergrid absorbing layer" not in text,
                       "ran past the abort into printPreamble")

    elif case_dir == "sg-source-2":
        # Identical geometry, override set: warn, print the whole diagnostic,
        # and run.
        ok &= _require("WARNING: 1 of 1 source(s)" in text,
                       "no source warning")
        ok &= _require("FATAL:" not in text, "aborted despite the override")
        ok &= _require("INSIDE the layer at face x=0" in text,
                       "override suppressed the diagnostic")
        ok &= _require("Supergrid absorbing layer" in text,
                       "did not reach printPreamble")
        ok &= _no_error_word(text)

    elif case_dir == "sg-margin-1":
        # W + 2h from the boundary: outside the layer, inside the stencil
        # margin. Pins the lower half of margin_pts = 5.
        ok &= _require("FATAL:" in text, "no FATAL line")
        ok &= _require("within the stencil margin of face x=0" in text,
                       "did not report a margin violation")
        ok &= _require("5 grid points required" in text,
                       "margin is not 5 grid points at 4th order")

    elif case_dir == "sg-margin-2":
        # W + 6h: must pass. Pins the upper half of margin_pts = 5.
        ok &= _require("FATAL:" not in text and "source(s) are inside" not in text,
                       "flagged a source that clears the layer by 6 points")
        # Receiver flag: one of the two stations is in the layer.
        ok &= _require("WARNING: receiver inlayer is" in text,
                       "did not flag the station inside the layer")
        ok &= _require("WARNING: receiver interior" not in text,
                       "flagged the interior station")
        ok &= _require("1 of 2 receivers are inside the supergrid" in text,
                       "wrong receiver summary")
        # The metric block, exactly once - printPreamble is re-entered per event
        # and per inversion iteration and pytest-sw4mopt reads inversion stdout
        # at a fixed offset from the end.
        n = _count(text, "Supergrid absorbing layer")
        ok &= _require(n == 1, "metric block printed %d times, expected 1" % n)
        # The printed T_max must equal W/(C*Vs) analytically, which pins the
        # formula against edits to Psi0.
        expect = TEST_WIDTH / (C_ADIAB * TEST_VS)
        got = None
        for line in text.splitlines():
            if "in the fastest material" in line:
                got = float(line.split()[0])
        ok &= _require(got is not None, "no fastest-material T_max line")
        if got is not None:
            ok &= _require(abs(got - expect) / expect < 1e-3,
                           "T_max = %g, expected W/(%.7f*Vs) = %g"
                           % (got, C_ADIAB, expect))
        # fc1 = 0.5 Hz asks for 2 s against 0.67 s absorbable, so the shortfall
        # warning must fire and must come from fc1, not fc2.
        ok &= _require("prefilter bandpass low corner fc1" in text,
                       "band edge not taken from fc1")
        ok &= _require("cannot absorb the source band" in text,
                       "no shortfall warning")
        ok &= _no_error_word(text)

    elif case_dir == "sg-rupture-1":
        # A rupture-derived source in the layer. If the guard only looked at
        # 'source' commands this case would run to completion.
        ok &= _require("FATAL:" in text, "no FATAL line")
        ok &= _require("of 1 source(s) are inside" in text,
                       "did not count the SRF subfault as a source")
        ok &= _require("INSIDE the layer at face y=0" in text,
                       "did not name the y=0 face")

    elif case_dir == "sg-topface-1":
        # lz=2 puts a sponge on the z=0 face. The old index-based predicate
        # never looked at it.
        ok &= _require("WARNING: receiver intop is" in text,
                       "did not flag the station in the z=0 layer")
        ok &= _require("z=0(top)" in text,
                       "did not attribute the flag to the top face")
        ok &= _require("WARNING: receiver interior" not in text,
                       "flagged the interior station")
        ok &= _require("1 of 2 receivers are inside the supergrid" in text,
                       "wrong receiver summary")
        # Adequate-sponge negative: a 0.31 s pulse against 0.67 s absorbable.
        # The shortfall warning must stay quiet, or it fires on every run and
        # stops meaning anything.
        ok &= _require("cannot absorb the source band" not in text,
                       "warned about a band the layer can absorb")
        ok &= _no_error_word(text)

    else:
        return _require(False, "unknown case " + case_dir)

    return bool(ok)
