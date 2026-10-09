"""Waveform, HDF5 output and cross-run check kinds."""

import glob
import json
import math
import re
from pathlib import Path

import h5py
import numpy as np

from . import outputs
from .checks import Ctx, check

COMPONENTS = ("NS", "EW", "UP")


# --------------------------------------------------------------------------
# rechdf5 output, read the way ~/src/workflow's lf_to_xarray reads it
# --------------------------------------------------------------------------
def read_rechdf5(path: Path) -> dict:
    out = {"stations": {}}
    with h5py.File(path, "r") as f:
        out["attrs"] = {k: (v.decode() if isinstance(v, bytes) else v) for k, v in f.attrs.items()}
        for k in ("DELTA", "DOWNSAMPLE", "ORIGINTIME", "SGWIDTH", "SGWIDTHGP"):
            if k in f:
                out[k] = float(np.asarray(f[k][()]).ravel()[0])
        for name, g in f.items():
            if not isinstance(g, h5py.Group):
                continue
            st = {k: np.asarray(g[k][()]) for k in g}
            # variables=velocity names the components Vns, Vew, Vup
            for v, c in (("Vns", "NS"), ("Vew", "EW"), ("Vup", "UP")):
                if v in st and c not in st:
                    st[c] = st[v]
            out["stations"][name] = st
    return out


def _station_inputs(ctx: Ctx, run: str) -> dict[str, tuple[float, float]]:
    p = ctx.params(run)
    return {k: tuple(v) for k, v in p.get("stations", {}).items()}


@check("rechdf5_contract")
def rechdf5_contract(ctx: Ctx) -> bool:
    """The properties of out.h5 that ~/src/workflow relies on (lf_to_xarray, im_calc)."""
    ok = True

    def need(cond: bool, msg: str):
        nonlocal ok
        print(f"  {'ok ' if cond else 'BAD'} {msg}")
        ok &= bool(cond)

    for run in ctx.run_names:
        p = ctx.params(run)
        d = read_rechdf5(ctx.run_dir(run) / ctx.opts.get("file", "out.h5"))
        print(f"{run}:")
        need("DELTA" in d and d["DELTA"] > 0, f"DELTA dataset present and positive ({d.get('DELTA')})")
        need(d.get("ORIGINTIME", 0.0) == 0.0, "ORIGINTIME is 0: recordings start at t=0")
        sts = d["stations"]
        inputs = _station_inputs(ctx, run)
        need(set(sts) == set(inputs), f"one group per input station ({sorted(sts)})")
        npts = {n: int(s["NPTS"].ravel()[0]) for n, s in sts.items()}
        need(len(set(npts.values())) == 1, f"NPTS equal across stations {set(npts.values())}")
        n = next(iter(npts.values()))
        expect_n = round(float(p["t"]) / d["DELTA"]) + 1
        need(abs(n - expect_n) <= 1, f"NPTS {n} ~ t/DELTA + 1 = {expect_n}")
        for name, s in sts.items():
            for c in COMPONENTS:
                arr = s.get(c)
                need(arr is not None and arr.shape == (n,) and np.all(np.isfinite(arr)),
                     f"{name}/{c}: dataset of {n} finite samples")
            need(int(s.get("ISNSEW", [1]).ravel()[0]) == 1, f"{name}: ISNSEW=1 (geographic NS/EW/UP)")
            if name in inputs:
                lat, lon = inputs[name]
                got = s["STLA,STLO,STDP"].ravel()
                dlon = (got[1] - lon + 180.0) % 360.0 - 180.0
                etol = ctx.tol("echo_tol_deg", {"double": 1e-9, "float": 1e-4})
                need(abs(got[0] - lat) < etol and abs(dlon) < etol and got[2] == 0.0,
                     f"{name}: STLA,STLO,STDP echo the station file ({got[0]:.9f}, {got[1]:.9f}, {got[2]:g})")
        if "sponge" in p:
            need("SGWIDTH" in d and abs(d["SGWIDTH"] - float(p["sponge"])) < 1e-6,
                 f"SGWIDTH = {d.get('SGWIDTH')} (sponge {p['sponge']})")
            need("SGWIDTHGP" in d and abs(d["SGWIDTHGP"] - float(p["sponge"]) / float(p["hcoarse"])) < 1e-6,
                 f"SGWIDTHGP = {d.get('SGWIDTHGP')} (sponge/hcoarse {float(p['sponge']) / float(p['hcoarse'])})")
            for name, s in sts.items():
                sg = float(s["SGDEPTH"].ravel()[0]) if "SGDEPTH" in s else math.nan
                sggp = float(s["SGDEPTHGP"].ravel()[0]) if "SGDEPTHGP" in s else math.nan
                inside = name in ctx.opts.get("in_sponge", ["SPNG"])
                need((sg > 0) == inside and (sggp > 0) == inside and np.isfinite(sg) and np.isfinite(sggp),
                     f"{name}: SGDEPTH={sg:g} SGDEPTHGP={sggp:g} ({'inside' if inside else 'outside'} the sponge)")
    return ok


@check("receiver_placement")
def receiver_placement(ctx: Ctx) -> bool:
    """Each receiver records at the grid node nearest to its requested position:
    horizontal distance <= h_surface * sqrt(2) / 2."""
    ok = True
    for run in ctx.run_names:
        p = ctx.params(run)
        h = float(p[ctx.opts.get("h_param", "hfine")])
        limit = h * math.sqrt(2) / 2 + 1e-6
        d = read_rechdf5(ctx.run_dir(run) / "out.h5")
        for name, s in d["stations"].items():
            want = s["STX,STY,STZ"].ravel()
            got = s["ACTUALSTX,STY,STZ"].ravel()
            dist = math.hypot(got[0] - want[0], got[1] - want[1])
            good = dist <= limit
            print(f"  {'ok ' if good else 'BAD'} {run}/{name}: requested ({want[0]:.1f}, {want[1]:.1f}) recorded at "
                  f"({got[0]:.1f}, {got[1]:.1f}); {dist:.1f} m away (nearest node is <= {limit:.1f} m)")
            ok &= good
    return ok


@check("projection")
def projection(ctx: Ctx) -> bool:
    """SW4's Cartesian position of each lat/lon station matches an independent
    pyproj calculation in the projection named by `frame`:
      frame = "sw4":  the PROJ string SW4 builds (+scale=...),
      frame = "nztm": what the workflow means (tmerc with +k = scale)."""
    from . import mininz

    frame = ctx.opts.get("frame", "sw4")
    tol = ctx.tol("tol_m", 0.01)
    ok = True
    for run in ctx.run_names:
        p = ctx.params(run)
        proj = (mininz.sw4_proj_string(173.0, 0.0, 0.9996) if frame == "sw4"
                else mininz.nztm_like_string(173.0, 0.0, 0.9996))
        grid = mininz.Grid(float(p["lon0"]), float(p["lat0"]), float(p["az"]), proj)
        d = read_rechdf5(ctx.run_dir(run) / "out.h5")
        for name, s in d["stations"].items():
            lat, lon, _ = s["STLA,STLO,STDP"].ravel()
            x, y = grid.to_xy(lon, lat)
            sx, sy, _ = s["STX,STY,STZ"].ravel()
            err = math.hypot(float(x) - sx, float(y) - sy)
            good = err <= tol
            print(f"  {'ok ' if good else 'BAD'} {run}/{name}: SW4 ({sx:.3f}, {sy:.3f}) pyproj[{frame}] "
                  f"({float(x):.3f}, {float(y):.3f}) differ by {err:.3f} m (tol {tol} m)")
            ok &= good
    return ok


def _traces(ctx: Ctx, run: str, out_dir: str | None = None) -> dict[str, dict[str, np.ndarray]]:
    d = read_rechdf5(ctx.run_dir(run, out_dir) / ctx.opts.get("file", "out.h5"))
    return {n: {c: s[c].astype(np.float64) for c in COMPONENTS if c in s} for n, s in d["stations"].items()}


def _bless_traces(ctx: Ctx) -> dict:
    return {"runs": {r: {n: {c: v.tolist() for c, v in comps.items()} for n, comps in _traces(ctx, r).items()}
                     for r in ctx.run_names}}


def trace_misfit(a: np.ndarray, b: np.ndarray) -> float:
    """Relative L2 misfit ||a - b|| / max(||a||, ||b||)."""
    n = max(np.linalg.norm(a), np.linalg.norm(b))
    return float(np.linalg.norm(a - b) / n) if n > 0 else 0.0


@check("traces_golden", bless=_bless_traces)
def traces_golden(ctx: Ctx) -> bool:
    """Every station trace matches the golden for this precision (relative L2)."""
    gold = ctx.golden()
    if gold is None:
        return False
    rtol = ctx.tol("rtol", {"double": 1e-6, "float": 2e-3})
    ok = True
    for r in ctx.run_names:
        got = _traces(ctx, r)
        ref = gold["runs"][r]
        if set(got) != set(ref):
            print(f"FAIL {r}: stations {sorted(got)} != golden {sorted(ref)}")
            ok = False
            continue
        worst = (0.0, "")
        for n in ref:
            for c, v in ref[n].items():
                a = np.asarray(v)
                b = got[n][c]
                if a.shape != b.shape:
                    print(f"FAIL {r}/{n}/{c}: {b.shape} samples, golden {a.shape}")
                    ok = False
                    continue
                m = trace_misfit(b, a)
                worst = max(worst, (m, f"{n}/{c}"))
                if m > rtol:
                    print(f"  BAD {r}/{n}/{c}: misfit {m:.2e} > {rtol:.1e}")
                    ok = False
        print(f"{r}: worst trace misfit {worst[0]:.2e} ({worst[1]}), tolerance {rtol:.1e}")
    return ok


# --------------------------------------------------------------------------
# imagehdf5
# --------------------------------------------------------------------------
def read_image(path: Path) -> dict:
    with h5py.File(path, "r") as f:
        d = {k: np.asarray(f[k][()]) for k in f}
    ni, nj = int(d["ni"][0]), int(d["nj"][0])
    d["data"] = d["patches"].reshape(nj, ni).T  # Fortran order: i fastest
    return d


def _images(run_dir: Path, file: str, mode_suffix: str) -> list[Path]:
    return sorted(Path(p) for p in glob.glob(str(run_dir / f"{file}.cycle=*.z=0.{mode_suffix}.sw4img.h5")))


# file prefix (from mininz.in) -> SW4 mode suffix in the file name
WORKFLOW_IMAGES = {
    "topo": "topo", "grid": "gridx", "surf_vp": "p", "surf_vs": "s", "surf_rho": "rho",
    "surf_mag": "mag", "surf_velmag": "magdudt", "surf_uz": "uz", "surf_hmax": "hmax", "surf_vmax": "vmax",
}


@check("images_contract")
def images_contract(ctx: Ctx) -> bool:
    """The 10 imagehdf5 outputs the workflow requests exist with the right times,
    sizes and spacing; hmax/vmax agree with the station traces; the topo image
    is the sfile topography; the surface material images match the model."""
    from . import mininz

    ok = True

    def need(cond: bool, msg: str):
        nonlocal ok
        print(f"  {'ok ' if cond else 'BAD'} {msg}")
        ok &= bool(cond)

    for run in ctx.run_names:
        p = ctx.params(run)
        rd = ctx.run_dir(run)
        t, dt_img, h = float(p["t"]), float(p["image_dt"]), float(p["hfine"])
        print(f"{run}:")
        imgs = {}
        for prefix, suffix in WORKFLOW_IMAGES.items():
            files = _images(rd, prefix, suffix)
            need(bool(files), f"{prefix}: {len(files)} file(s)")
            if not files:
                continue
            imgs[prefix] = [read_image(f) for f in files]
        if not ok:
            continue
        for prefix in ("topo", "grid", "surf_vp", "surf_vs", "surf_rho"):
            need(len(imgs[prefix]) == 1 and float(imgs[prefix][0]["time"][0]) == 0.0, f"{prefix}: one image at t=0")
        n_expect = int(round(t / dt_img)) + 1
        for prefix in ("surf_mag", "surf_velmag", "surf_uz"):
            times = [float(i["time"][0]) for i in imgs[prefix]]
            need(len(times) == n_expect and abs(times[-1] - t) < 1e-4 * max(1.0, t),
                 f"{prefix}: {len(times)} images at t={times[0]:g}..{times[-1]:g} every {dt_img:g} s (expect {n_expect})")
        for prefix in ("surf_hmax", "surf_vmax"):
            need(len(imgs[prefix]) == 1 and abs(float(imgs[prefix][0]["time"][0]) - t) < 1e-4 * max(1.0, t),
                 f"{prefix}: one image at the final time t={t:g}")
        for prefix, lst in imgs.items():
            need(all(abs(float(i["grid_size"][0]) - h) < 1e-9 for i in lst), f"{prefix}: grid_size = {h:g}")
        # Topography image vs the sfile surface, and the computational surface
        side = float(p["lx"])
        hill = (side / 2 - 1000.0, side / 2 + 500.0)
        topo = imgs["topo"][0]["data"]
        X, Y = np.meshgrid(np.arange(topo.shape[0]) * h, np.arange(topo.shape[1]) * h, indexing="ij")
        ana = mininz.topography(X, Y, *hill) if p.get("elevation", True) else np.zeros_like(X)
        if ctx.opts.get("check_topo", True):
            need(np.abs(topo - ana).max() < 0.05, f"topo image = sfile topography (max diff {np.abs(topo - ana).max():.3g} m)")
        # Surface material at depth 0. Density is the model's; with attenuation
        # SW4 stores (and images) unrelaxed moduli, so velocities come out above
        # the model's phasefreq velocities by the dispersion correction (~6% at Qs=25).
        rho0, vp0, vs0, _, _ = mininz.material(0.0)
        lo, hi = ctx.opts.get("velocity_ratio_range", [1.0, 1.08])
        img = imgs["surf_rho"][0]["data"]
        rel = float(np.abs(img / rho0 - 1).min())
        need(rel < 1e-4, f"surf_rho: {img.min():.1f}..{img.max():.1f}, model {rho0:.1f} at the surface (min rel diff {rel:.1e})")
        for prefix, want in (("surf_vp", vp0), ("surf_vs", vs0)):
            img = imgs[prefix][0]["data"]
            r = float(img.min() / want)
            need(lo <= r <= hi,
                 f"{prefix}: minimum {img.min():.1f} = {r:.4f} x the model's {want:.1f} (unrelaxed; expect {lo}..{hi})")
        # hmax / vmax at station nodes vs the station traces (sampled every step)
        d = read_rechdf5(rd / "out.h5")
        hmax = imgs["surf_hmax"][0]["data"]
        vmax = imgs["surf_vmax"][0]["data"]
        for name, s in d["stations"].items():
            ax, ay, az = s["ACTUALSTX,STY,STZ"].ravel()
            i, j = int(round(ax / h)), int(round(ay / h))
            th = float(np.max(np.hypot(s["NS"], s["EW"])))
            tv = float(np.max(np.abs(s["UP"])))
            rh = abs(hmax[i, j] - th) / max(th, 1e-30)
            rv = abs(vmax[i, j] - tv) / max(tv, 1e-30)
            need(rh < 2e-3 and rv < 2e-3,
                 f"{name}: hmax image {hmax[i, j]:.4e} vs trace {th:.4e}; vmax {vmax[i, j]:.4e} vs {tv:.4e}")
            if p.get("elevation", True):
                surf = -az
                need(abs(surf - ana[i, j]) < ctx.opts.get("surface_tol_m", 30.0),
                     f"{name}: computational surface {surf:.1f} m vs sfile {ana[i, j]:.1f} m (smoothing allowance)")
    return ok


# --------------------------------------------------------------------------
# Comparisons between runs (metamorphic tests)
# --------------------------------------------------------------------------
def _load_series(ctx: Ctx, run: str, spec: dict) -> dict[str, np.ndarray]:
    """Traces from either a rechdf5 file or `rec usgsformat=1` text files."""
    rd = ctx.run_dir(run)
    if spec.get("format", "rechdf5") == "rechdf5":
        return {f"{n}/{c}": v for n, comps in _traces(ctx, run).items() for c, v in comps.items()
                if not n.startswith(spec.get("ignore_prefix", "\0"))}
    out = {}
    for f in sorted(rd.glob(spec.get("glob", "*.txt"))):
        if f.name.endswith("Err.txt"):
            continue
        t, data, names = outputs.read_usgs(f)
        for k in range(data.shape[1]):
            comp = names[k] if k < len(names) else str(k)
            stem = f.stem[len(spec.get("strip", "")):] if f.stem.startswith(spec.get("strip", "")) else f.stem
            out[f"{stem}/{_canon(comp)}"] = data[:, k]
    return out


def _canon(name: str) -> str:
    """Map USGS column titles to NS/EW/UP (or X/Y/Z)."""
    n = name.lower()
    for key, short in (("north", "NS"), ("east", "EW"), ("up", "UP"), ("x ", "X"), ("y ", "Y"), ("z ", "Z")):
        if key in n + " ":
            return short
    return name


@check("traces_agree")
def traces_agree(ctx: Ctx) -> bool:
    """Traces of `runs[0]` and `runs[1]` agree (relative L2 per trace <= rtol).

    Options: `a_format`/`b_format` ("rechdf5" | "usgs"), `transform` applied to
    run b before comparing: "none" | "derivative" | "integral", `scale_b`,
    `components` to restrict to, `skip` initial samples, `bitwise`.
    """
    ra, rb = ctx.run_names[:2]
    a = _load_series(ctx, ra, {"format": ctx.opts.get("a_format", "rechdf5"), "ignore_prefix": ctx.opts.get("a_ignore", "\0")})
    b = _load_series(ctx, rb, {"format": ctx.opts.get("b_format", "rechdf5"), "strip": ctx.opts.get("b_strip", "")})
    comps = ctx.opts.get("components")
    if comps:
        a = {k: v for k, v in a.items() if k.split("/")[-1] in comps}
    keys = sorted(set(a) & set(b))
    if not keys:
        print(f"FAIL: no common traces between {ra} ({sorted(a)[:4]}..) and {rb} ({sorted(b)[:4]}..)")
        return False
    if set(a) != set(b):
        print(f"  note: comparing {len(keys)} common traces of {len(a)} / {len(b)}")
    rtol = ctx.tol("rtol", {"double": 1e-6, "float": 1e-3})
    bw = ctx.opts.get("bitwise", False)
    bitwise = bool(bw.get(ctx.cfg.precision, False) if isinstance(bw, dict) else bw)
    skip = int(ctx.opts.get("skip", 0))
    # normalize = "station": divide by the station's largest component norm, so a
    # weak component near a nodal plane is judged against the station's motion.
    station_norm = {}
    if ctx.opts.get("normalize") == "station":
        for k in keys:
            st = k.split("/")[0]
            station_norm[st] = max(station_norm.get(st, 0.0), float(np.linalg.norm(b[k][skip:])))
    ok = True
    worst = (0.0, "")
    for k in keys:
        x, y = a[k], b[k]
        if x.shape != y.shape:
            print(f"  BAD {k}: {x.shape} vs {y.shape} samples")
            ok = False
            continue
        if bitwise:
            same = np.array_equal(x, y)
            if not same:
                print(f"  BAD {k}: not bitwise identical (max diff {np.abs(x - y).max():.3e})")
            ok &= same
            continue
        if station_norm:
            m = float(np.linalg.norm(x[skip:] - y[skip:]) / station_norm[k.split("/")[0]])
        else:
            m = trace_misfit(x[skip:], y[skip:])
        worst = max(worst, (m, k))
        if m > rtol:
            print(f"  BAD {k}: misfit {m:.3e} > {rtol:.1e}")
            ok = False
    if bitwise:
        print(f"{len(keys)} traces {'bitwise identical' if ok else 'differ'}")
    else:
        print(f"{len(keys)} traces; worst misfit {worst[0]:.3e} ({worst[1]}), tolerance {rtol:.1e}")
    return ok


@check("peer_traces_agree")
def peer_traces_agree(ctx: Ctx) -> bool:
    """Float vs double: the same run in the other precision's build (cfg.peer_out_dir)."""
    if not ctx.cfg.peer_out_dir:
        print("FAIL: no peer build configured (SW4TEST_PEER_OUT_DIR)")
        return False
    rtol = float(ctx.opts.get("rtol", 1e-4))
    ok = True
    for r in ctx.run_names:
        a = _traces(ctx, r)
        b = _traces(ctx, r, ctx.cfg.peer_out_dir)
        worst = (0.0, "")
        for n in a:
            for c in a[n]:
                m = trace_misfit(a[n][c], b[n][c])
                worst = max(worst, (m, f"{n}/{c}"))
                if m > rtol:
                    print(f"  BAD {r}/{n}/{c}: {ctx.cfg.precision} vs {ctx.cfg.peer_precision} misfit {m:.3e} > {rtol:.1e}")
                    ok = False
        print(f"{r}: worst {ctx.cfg.precision}-vs-{ctx.cfg.peer_precision} misfit {worst[0]:.3e} ({worst[1]})")
    return ok


@check("stdout_count")
def stdout_count(ctx: Ctx) -> bool:
    """`pattern` occurs exactly `count` times in each run's output."""
    ok = True
    for r in ctx.run_names:
        text = (ctx.run_dir(r) / "sw4.out").read_text(errors="replace")
        n = len(re.findall(ctx.opts["pattern"], text))
        good = n == int(ctx.opts["count"])
        print(f"  {'ok ' if good else 'BAD'} {r}: {n} matches of {ctx.opts['pattern']!r} (expect {ctx.opts['count']})")
        ok &= good
    return ok


@check("station_relations")
def station_relations(ctx: Ctx) -> bool:
    """Symmetry: trace A/comp == sign x trace B/comp for each [A, compA, B, compB, sign]."""
    rtol = ctx.tol("rtol", {"double": 1e-8, "float": 1e-4})
    ok = True
    for r in ctx.run_names:
        tr = _traces(ctx, r)
        for a, ca, b, cb, sign in ctx.opts["relations"]:
            x, y = tr[a][ca], float(sign) * tr[b][cb]
            m = trace_misfit(x, y)
            good = m <= rtol
            print(f"  {'ok ' if good else 'BAD'} {r}: {a}/{ca} vs {sign:+g} x {b}/{cb}: misfit {m:.2e} (tol {rtol:.0e})")
            ok &= good
    return ok


@check("rechdf5_only_requested")
def rechdf5_only_requested(ctx: Ctx) -> bool:
    """out.h5 holds exactly the stations of its infile, each with data (NPTS > 0)."""
    ok = True
    for r in ctx.run_names:
        with h5py.File(ctx.run_dir(r) / "stations.h5") as f:
            want = {k for k in f if isinstance(f[k], h5py.Group)}
        d = read_rechdf5(ctx.run_dir(r) / "out.h5")
        got = set(d["stations"])
        extra = sorted(got - want)
        empty = sorted(n for n, s in d["stations"].items() if int(s["NPTS"].ravel()[0]) == 0)
        print(f"  {r}: stations {sorted(got)}; not in infile {extra}; empty (NPTS=0) {empty}")
        ok &= not extra and not empty
    return ok


@check("component_orientation")
def component_orientation(ctx: Ctx) -> bool:
    """The NS/EW components point to true north/east.

    runs = [nsew_run, xyz_run]: the same problem recorded as NS/EW/UP and as raw
    grid X/Y/Z. For each station, fit [NS, EW] = M [X, Y] and compare the
    directions of M's rows with true north/east in the grid frame (pyproj, in
    the projection SW4 uses)."""
    from . import mininz

    tol = float(ctx.opts.get("tol_deg", 0.05))
    rn, rx = ctx.run_names[:2]
    p = ctx.params(rn)
    grid = mininz.Grid(float(p["lon0"]), float(p["lat0"]), float(p["az"]), mininz.sw4_proj_string(173.0, 0.0, 0.9996))
    a = read_rechdf5(ctx.run_dir(rn) / "out.h5")["stations"]
    with h5py.File(ctx.run_dir(rx) / "out.h5") as f:
        b = {n: {c: np.asarray(g[c][()], float) for c in ("X", "Y")} for n, g in f.items() if isinstance(g, h5py.Group)}
    ok = True
    for name, s in a.items():
        A = np.vstack([b[name]["X"], b[name]["Y"]]).T
        M = np.vstack([np.linalg.lstsq(A, s[c].astype(float), rcond=None)[0] for c in ("NS", "EW")])
        lat, lon, _ = s["ACTUALSTLA,STLO,STDP"].ravel() if "ACTUALSTLA,STLO,STDP" in s else s["STLA,STLO,STDP"].ravel()
        x0, y0 = grid.to_xy(lon, lat)
        xn, yn = grid.to_xy(lon, lat + 1e-4)
        xe, ye = grid.to_xy(lon + 1e-4, lat)
        north = np.array([xn - x0, yn - y0]); north /= np.linalg.norm(north)
        east = np.array([xe - x0, ye - y0]); east /= np.linalg.norm(east)
        ang = lambda u, v: math.degrees(math.atan2(u[0] * v[1] - u[1] * v[0], u @ v))
        ens = ang(north, M[0] / np.linalg.norm(M[0]))
        eew = ang(east, M[1] / np.linalg.norm(M[1]))
        good = abs(ens) <= tol and abs(eew) <= tol
        print(f"  {'ok ' if good else 'BAD'} {name} (lon {lon:.3f}): NS axis off true north by {ens:+.3f} deg, "
              f"EW axis off true east by {eew:+.3f} deg (tol {tol})")
        ok &= good
    return ok


@check("restart_agree")
def restart_agree(ctx: Ctx) -> bool:
    """runs = [straight, restarted]: the restarted run's samples equal the
    straight run's samples at the same times (aligned on the trace ends)."""
    rs, rr = ctx.run_names[:2]
    a, b = _traces(ctx, rs), _traces(ctx, rr)
    bw = ctx.opts.get("bitwise", True)
    bitwise = bool(bw.get(ctx.cfg.precision, True) if isinstance(bw, dict) else bw)
    rtol = ctx.tol("rtol", {"double": 1e-12, "float": 1e-6})
    ok = True
    for n in a:
        for c in a[n]:
            x, y = a[n][c], b[n][c]
            k = min(len(x), len(y))
            x, y = x[-k:], y[-k:]
            good = np.array_equal(x, y) if bitwise else trace_misfit(x, y) <= rtol
            if not good:
                print(f"  BAD {n}/{c}: last {k} samples differ (max {np.abs(x - y).max():.3e})")
            ok &= good
    print(f"compared the last {k} samples of {len(a)} stations ({'bitwise' if bitwise else f'rtol {rtol}'})")
    return ok


@check("prefilter_response")
def prefilter_response(ctx: Ctx) -> bool:
    """runs = [unfiltered, filtered]: the filtered seismograms equal the unfiltered
    ones passed through a zero-phase Butterworth lowpass (scipy butter + sosfiltfilt
    for passes=2), so SW4's filter has the intended corner, order and no shift.
    Misfit is normalised per station (largest component norm)."""
    from scipy.signal import butter, sosfilt, sosfiltfilt

    ru, rf = ctx.run_names[:2]
    fc, order, passes = float(ctx.opts["fc"]), int(ctx.opts.get("order", 4)), int(ctx.opts.get("passes", 2))
    rtol = float(ctx.opts.get("rtol", 0.02))
    a, b = _traces(ctx, ru), _traces(ctx, rf)
    dt = read_rechdf5(ctx.run_dir(ru) / "out.h5")["DELTA"]
    sos = butter(order, fc, btype="low", fs=1.0 / dt, output="sos")
    ok = True
    for n in a:
        norm = max(np.linalg.norm(b[n][c]) for c in b[n])
        for c in ("NS", "EW", "UP"):
            x = a[n][c]
            # zeros before the source starts; hold the final (static) value after the record ends
            xx = np.concatenate([np.zeros(len(x)), x, np.full(len(x), x[-1])])
            y = (sosfiltfilt(sos, xx, padtype=None) if passes == 2 else sosfilt(sos, xx))[len(x):2 * len(x)]
            # The backward pass needs data past the record end, which the reference
            # does not have: compare only up to `compare_until` seconds.
            k = int(float(ctx.opts.get("compare_until", 1e30)) / dt)
            m = float(np.linalg.norm(b[n][c][:k] - y[:k]) / norm)
            good = m <= rtol
            print(f"  {'ok ' if good else 'BAD'} {n}/{c}: SW4 prefilter vs scipy Butterworth (order {order}, "
                  f"{passes} pass{'es' if passes > 1 else ''}, fc {fc} Hz): misfit {m:.2e} (tol {rtol})")
            ok &= good
    return ok


@check("shifted_agree")
def shifted_agree(ctx: Ctx) -> bool:
    """runs = [a, b]: trace b(t) == trace a(t + shift) over b's record (linear
    interpolation of a), relative L2 <= rtol. Used to compare one source at two
    start times."""
    ra, rb = ctx.run_names[:2]
    shift = float(ctx.opts["shift"])
    a, b = _traces(ctx, ra), _traces(ctx, rb)
    da = read_rechdf5(ctx.run_dir(ra) / "out.h5")["DELTA"]
    db = read_rechdf5(ctx.run_dir(rb) / "out.h5")["DELTA"]
    rtol = float(ctx.opts.get("rtol", 1e-2))
    ok = True
    worst = (0.0, "")
    for n in b:
        for c in b[n]:
            tb = np.arange(len(b[n][c])) * db
            ta = np.arange(len(a[n][c])) * da
            sel = tb + shift <= ta[-1]
            x = np.interp(tb[sel] + shift, ta, a[n][c])
            m = trace_misfit(b[n][c][sel], x)
            worst = max(worst, (m, f"{n}/{c}"))
            ok &= m <= rtol
    print(f"worst misfit {worst[0]:.3e} ({worst[1]}) between {rb}(t) and {ra}(t + {shift:g} s); tolerance {rtol:g}")
    return ok


@check("traces_bounded")
def traces_bounded(ctx: Ctx) -> bool:
    """Long-time stability: every trace is finite and its largest value over the
    last `tail` fraction of the record is <= `ratio` x its peak over the first half."""
    tail, ratio = float(ctx.opts.get("tail", 0.25)), float(ctx.opts.get("ratio", 1.0))
    ok = True
    for r in ctx.run_names:
        res = json.loads((ctx.run_dir(r) / "result.json").read_text())
        if res.get("returncode", 0) != 0:
            print(f"  BAD {r}: SW4 aborted (exit {res['returncode']}): the solution did not stay bounded")
            ok = False
            continue
        worst = (0.0, "")
        for n, comps in _traces(ctx, r).items():
            for c, v in comps.items():
                if not np.all(np.isfinite(v)):
                    print(f"  BAD {r}/{n}/{c}: non-finite samples from t index {int(np.argmin(np.isfinite(v)))}")
                    ok = False
                    continue
                k = len(v)
                early = np.abs(v[: k // 2]).max()
                late = np.abs(v[int(k * (1 - tail)):]).max()
                q = late / early if early > 0 else 0.0
                worst = max(worst, (q, f"{n}/{c}"))
                ok &= q <= ratio
        print(f"{r}: late/early amplitude ratio worst {worst[0]:.3e} ({worst[1]}), limit {ratio}")
    return ok


@check("loh_reference")
def loh_reference(ctx: Ctx) -> bool:
    """Station 10 of LOH.1/LOH.3 against the PROSE reference (tools/LOH.*_prose*).

    The SW4 run records grid x/y/z displacement at (6 km, 8 km) from the source
    with a Gaussian moment rate, which equals ground velocity for the reference's
    step moment. Radial = (6x + 8y)/10, transverse = (-8x + 6y)/10 (the
    receiver azimuth), vertical = z (positive down). Misfit is relative L2 over
    [0, compare_until] after resampling the reference to SW4's time step."""
    from . import loh

    ref = loh.exact(Path(ctx.cfg.source_dir) / "tools" / ctx.opts["reference"], float(ctx.opts["sigma"]))
    rtol = float(ctx.opts.get("rtol", 0.05))
    until = float(ctx.opts.get("compare_until", 9.0))
    ok = True
    for r in ctx.run_names:
        t, d, names = outputs.read_usgs(ctx.run_dir(r) / ctx.opts.get("file", "sta10.txt"))
        x, y, z = d[:, 0], d[:, 1], d[:, 2]
        sw4 = {"radial": (6 * x + 8 * y) / 10, "transverse": (-8 * x + 6 * y) / 10, "vertical": z}
        sel = t <= until
        for comp in ("radial", "transverse", "vertical"):
            rr = np.interp(t[sel], ref["t"], ref[comp])
            m = trace_misfit(sw4[comp][sel], rr)
            peak = np.abs(sw4[comp][sel]).max() / max(np.abs(rr).max(), 1e-30)
            good = m <= rtol
            print(f"  {'ok ' if good else 'BAD'} {r} {comp}: misfit vs PROSE {m:.3e} (tol {rtol}); peak ratio {peak:.3f}")
            ok &= good
    return ok
