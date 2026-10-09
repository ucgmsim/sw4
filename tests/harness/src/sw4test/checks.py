"""Check kinds. Each takes a Ctx and returns True on pass, printing why."""

import json
import math
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable

import numpy as np

from . import outputs
from .cases import Case, Check
from .config import Config

CHECKS: dict[str, Callable[["Ctx"], bool]] = {}
BLESS: dict[str, Callable[["Ctx"], dict]] = {}


def check(kind: str, bless: Callable[["Ctx"], dict] | None = None):
    def deco(fn):
        CHECKS[kind] = fn
        if bless:
            BLESS[kind] = bless
        return fn

    return deco


@dataclass
class Ctx:
    cfg: Config
    case: Case
    check: Check

    @property
    def opts(self) -> dict[str, Any]:
        return self.check.options

    @property
    def run_names(self) -> list[str]:
        return self.check.runs or list(self.case.runs)

    def run_dir(self, run: str, out_dir: str | None = None) -> Path:
        return self.cfg.run_dir(self.case.id, run, out_dir)

    def params(self, run: str) -> dict[str, Any]:
        p = self.run_dir(run) / "params.json"
        return json.loads(p.read_text()) if p.exists() else dict(self.case.runs[run].params)

    def tol(self, key: str, default: float | dict) -> float:
        """Option that may be a scalar or {double = .., float = ..}."""
        v = self.opts.get(key, default)
        if isinstance(v, dict):
            return float(v[self.cfg.precision])
        return float(v)

    # Goldens -------------------------------------------------------------
    @property
    def golden_path(self) -> Path:
        return self.cfg.golden_dir / self.case.id / f"{self.check.name}.json"

    def golden(self) -> dict | None:
        p = self.golden_path
        if not p.exists():
            print(f"FAIL: no golden file {p}. Generate one with:\n"
                  f"  uv run --project tests/harness sw4test bless --config <build>/tests/sw4test.json {self.case.id}")
            return None
        return json.loads(p.read_text())


def run_check(cfg: Config, case: Case, check: Check) -> bool:
    fn = CHECKS.get(check.kind)
    if fn is None:
        raise ValueError(f"unknown check kind {check.kind!r}")
    if check.options.get("description"):
        print(check.options["description"])
    ok = fn(Ctx(cfg, case, check))
    xfail = check.options.get("xfail")
    if isinstance(xfail, dict):
        xfail = xfail.get(cfg.precision)
    if xfail:
        if ok:
            print(f"XPASS: check now passes; remove xfail ({xfail})")
            return False
        print(f"XFAIL (known defect, still present): {xfail}")
        return True
    return ok


def relerr(a: float, b: float) -> float:
    d = max(abs(a), abs(b))
    return abs(a - b) / d if d > 0 else 0.0


def compare_values(label: str, got: list[float], ref: list[float], rtol: float, atol: float = 0.0) -> bool:
    ok = True
    if len(got) != len(ref):
        print(f"FAIL {label}: {len(got)} values, golden has {len(ref)}")
        return False
    for i, (g, r) in enumerate(zip(got, ref)):
        e = relerr(g, r)
        good = e <= rtol or abs(g - r) <= atol
        flag = "ok " if good else "BAD"
        print(f"  {flag} {label}[{i}]: got {g:.8e} golden {r:.8e} rel {e:.2e} (rtol {rtol:.1e})")
        ok &= good
    return ok


# --------------------------------------------------------------------------
# Error norms written by SW4's exact-solution modes
# --------------------------------------------------------------------------
NORMS = ("errInf", "errL2", "solInf")


def _twilight_values(ctx: Ctx, run: str) -> dict[str, list[float]]:
    return outputs.read_twilight_err(ctx.run_dir(run) / "TwilightErr.txt")


def _bless_twilight(ctx: Ctx) -> dict:
    return {"runs": {r: _twilight_values(ctx, r) for r in ctx.run_names}}


@check("twilight_golden", bless=_bless_twilight)
def twilight_golden(ctx: Ctx) -> bool:
    """TwilightErr.txt error norms match the golden values for this precision."""
    gold = ctx.golden()
    if gold is None:
        return False
    rtol = ctx.tol("rtol", {"double": 1e-6, "float": 1e-3})
    ok = True
    for r in ctx.run_names:
        got = _twilight_values(ctx, r)
        ref = gold["runs"][r]
        for key in ref:
            if key not in got:
                print(f"FAIL {r}: no {key} block in output")
                ok = False
                continue
            ok &= compare_values(f"{r}.{key}", got[key], ref[key], rtol)
    return ok


@check("twilight_agree")
def twilight_agree(ctx: Ctx) -> bool:
    """TwilightErr.txt norms are the same in every run (e.g. across MPI rank counts)."""
    rtol = ctx.tol("rtol", {"double": 1e-10, "float": 1e-4})
    runs = ctx.run_names
    ref = _twilight_values(ctx, runs[0])
    ok = True
    for r in runs[1:]:
        got = _twilight_values(ctx, r)
        for key in ref:
            ok &= compare_values(f"{r} vs {runs[0]} {key}", got[key], ref[key], rtol)
    return ok


def _series_summary(ctx: Ctx, run: str) -> dict[str, list[float]]:
    a = outputs.read_err_series(ctx.run_dir(run) / ctx.opts["file"])
    # Final errInf/errL2/solInf, and the largest errInf/errL2 over the run.
    return {"final": list(map(float, a[-1, 1:])), "max": [float(a[:, 1].max()), float(a[:, 2].max())]}


def _bless_series(ctx: Ctx) -> dict:
    return {"runs": {r: _series_summary(ctx, r) for r in ctx.run_names}}


@check("err_golden", bless=_bless_series)
def err_golden(ctx: Ctx) -> bool:
    """LambErr/PointSourceErr/RayleighErr series match golden final and peak norms."""
    gold = ctx.golden()
    if gold is None:
        return False
    rtol = ctx.tol("rtol", {"double": 1e-6, "float": 1e-3})
    ok = True
    for r in ctx.run_names:
        got = _series_summary(ctx, r)
        for key in ("final", "max"):
            ok &= compare_values(f"{r}.{key}", got[key], gold["runs"][r][key], rtol)
    return ok


@check("err_ratio")
def err_ratio(ctx: Ctx) -> bool:
    """errInf/errL2 of each run (`stat` final or max) are at most max_ratio x the reference run's."""
    ref_run = ctx.opts["reference"]
    max_ratio = float(ctx.opts.get("max_ratio", 2.0))
    stat = ctx.opts.get("stat", "final")
    ref = _series_summary(ctx, ref_run)[stat]
    ok = True
    for r in ctx.run_names:
        e = _series_summary(ctx, r)[stat]
        for i, norm in enumerate(("errInf", "errL2")):
            ratio = e[i] / ref[i]
            good = ratio <= max_ratio
            print(f"  {'ok ' if good else 'BAD'} {r}: {stat} {norm} {e[i]:.4e} = {ratio:.2f} x {ref_run} (limit {max_ratio})")
            ok &= good
    return ok


def _error_for_rate(ctx: Ctx, run: str, norm: str) -> tuple[float, float]:
    """Return (error, solution magnitude) used in the convergence check."""
    src = ctx.opts.get("source", "twilight")
    idx = NORMS.index(norm)
    if src == "twilight":
        v = _twilight_values(ctx, run)[ctx.opts.get("field", "disp")]
        return v[idx], v[2]
    a = outputs.read_err_series(ctx.run_dir(run) / src)
    stat = ctx.opts.get("stat", "final")
    if stat == "final":
        return float(a[-1, 1 + idx]), float(a[-1, 3])
    if stat == "max":
        return float(a[:, 1 + idx].max()), float(a[:, 3].max())
    raise ValueError(f"unknown stat {stat}")


@check("rate")
def rate(ctx: Ctx) -> bool:
    """Observed order of accuracy between successive resolutions.

    Each run must define param `h`. The order on the finest pair must lie in
    [expect - band, expect + band]; coarser pairs are printed for information.
    `band_<norm>` overrides the band for one norm: at PR resolutions the max
    norm is often still pre-asymptotic (e.g. 3.5 with attenuation) while L2 is
    already at 4. Pairs whose fine-grid error is below `floor` x solution magnitude are
    skipped: there the error is round-off, not truncation (matters in float).
    """
    expect = float(ctx.opts["expect"])
    # float_sw4 round-off perturbs the observed order at PR resolutions: wider default band
    band = ctx.tol("band", {"double": 0.4, "float": 0.8})
    if ctx.cfg.precision == "float":
        band = max(band, 0.8)
    floor = ctx.tol("floor", {"double": 1e-11, "float": 2e-5})
    norms = ctx.opts.get("norms", ["errInf", "errL2"])
    runs = sorted(ctx.run_names, key=lambda r: -float(ctx.params(r)["h"]))
    if len(runs) < 2:
        print("FAIL: a rate check needs at least two runs")
        return False
    ok = True
    for norm in norms:
        rows = []
        for r in runs:
            e, s = _error_for_rate(ctx, r, norm)
            rows.append((r, float(ctx.params(r)["h"]), e, s))
        print(f"{norm}:")
        last = None
        for (r1, h1, e1, s1), (r2, h2, e2, s2) in zip(rows, rows[1:]):
            if e2 <= floor * max(s2, 1e-300):
                print(f"  {r1}->{r2}: fine error {e2:.3e} is below round-off floor {floor:.0e}*{s2:.3e}; skipped")
                continue
            p = math.log(e1 / e2) / math.log(h1 / h2)
            print(f"  {r1} (h={h1:g}, e={e1:.4e}) -> {r2} (h={h2:g}, e={e2:.4e}): p = {p:.3f}")
            last = p
        if last is None:
            # Every fine-grid error is already at the round-off floor: as accurate
            # as this precision allows; the golden check guards the values.
            print(f"  ok: all fine-grid {norm} errors are at the {ctx.cfg.precision} round-off floor")
        else:
            b = ctx.tol(f"band_{norm}", band) if f"band_{norm}" in ctx.opts else band
            if ctx.cfg.precision == "float":
                b = max(b, band)
            if abs(last - expect) > b:
                print(f"FAIL {norm}: finest-pair order {last:.3f} outside {expect} +/- {b}")
                ok = False
            else:
                print(f"  ok: finest-pair order {last:.3f} within {expect} +/- {b}")
    return ok


# --------------------------------------------------------------------------
# Energy (testenergy) and stdout
# --------------------------------------------------------------------------
@check("energy")
def energy(ctx: Ctx) -> bool:
    """Discrete energy from `testenergy`.

    mode = "conserve": |E_n - E_skip| / E_skip <= tol for all n >= skip
    mode = "decay":    E_{n+1} <= E_n (1 + tol) for all n >= skip
    mode = "grow":     E_end / E_skip >= factor (proves the detector sees an instability)
    """
    mode = ctx.opts.get("mode", "decay")
    skip = int(ctx.opts.get("skip", 2))
    ok = True
    for r in ctx.run_names:
        e = outputs.read_energy(ctx.run_dir(r) / ctx.opts.get("file", "energy.log"))
        print(f"{r}: {len(e)} samples, E[{skip}]={e[skip]:.10e}, E[end]={e[-1]:.10e}")
        if not np.all(np.isfinite(e)):
            print(f"FAIL {r}: non-finite energy at sample {int(np.argmin(np.isfinite(e)))}")
            ok = False
            continue
        e0 = e[skip]
        if mode == "conserve":
            tol = ctx.tol("tol", {"double": 1e-7, "float": 1e-3})
            drift = float(np.max(np.abs(e[skip:] - e0)) / e0)
            print(f"  max relative drift {drift:.3e} (tol {tol:.1e})")
            ok &= drift <= tol
        elif mode == "decay":
            tol = ctx.tol("tol", {"double": 1e-8, "float": 2e-3})
            rel = np.diff(e[skip:]) / e[skip:-1]
            worst = int(np.argmax(rel))
            print(f"  largest step-to-step relative increase {rel[worst]:.3e} at step {worst + skip} (tol {tol:.1e});"
                  f" total change {(e[-1] - e0) / e0:.3e}")
            ok &= bool(rel[worst] <= tol)
        elif mode == "grow":
            factor = float(ctx.opts.get("factor", 10.0))
            g = e[-1] / e0
            print(f"  growth E_end/E_skip = {g:.3e} (need >= {factor:g})")
            ok &= g >= factor
        else:
            raise ValueError(f"unknown energy mode {mode}")
    return ok


def _output_text(ctx: Ctx, run: str) -> str:
    d = ctx.run_dir(run)
    return (d / "sw4.out").read_text(errors="replace") + (d / "sw4.err").read_text(errors="replace")


@check("stdout_absent")
def stdout_absent(ctx: Ctx) -> bool:
    """None of `patterns` appears in any run's output."""
    pats = ctx.opts.get("patterns", [ctx.opts.get("pattern")])
    ok = True
    for r in ctx.run_names:
        text = _output_text(ctx, r)
        for p in pats:
            hits = [ln for ln in text.splitlines() if re.search(p, ln)]
            if hits:
                print(f"FAIL {r}: {len(hits)} lines match {p!r}, first: {hits[0].strip()}")
                ok = False
    if ok:
        print(f"ok: no output matches {pats}")
    return ok


@check("stdout_present")
def stdout_present(ctx: Ctx) -> bool:
    pats = ctx.opts.get("patterns", [ctx.opts.get("pattern")])
    ok = True
    for r in ctx.run_names:
        text = _output_text(ctx, r)
        for p in pats:
            if not re.search(p, text, re.MULTILINE):
                print(f"FAIL {r}: no output matches {p!r}")
                ok = False
    return ok


# Further check kinds register themselves on import.
from . import checks_waveforms  # noqa: E402,F401
