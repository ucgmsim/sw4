"""Inputs for the physics-equivalence (metamorphic) tests in tests/cases/equivalence.

One point source in a homogeneous or two-layer half-space, written several ways
(SRF-HDF5, ASCII SRF, `source` command) with stations given as lat/lon or as
grid x,y. Positions are defined in a projected frame around a fixed centre so
the same physical problem can be posed on grids with different azimuths.
"""

import math
from pathlib import Path

import h5py
import numpy as np

from . import mininz


def _frame(params: dict) -> tuple[mininz.Grid, float, float]:
    """Grid whose square domain (side L) is centred on (clon, clat) for any azimuth."""
    L = float(params["L"])
    az = float(params.get("az", 0.0))
    proj = mininz.sw4_proj_string(173.0, 0.0, 0.9996)
    base = mininz.Grid(float(params.get("clon", 172.6)), float(params.get("clat", -43.5)), 0.0, proj)
    a = math.radians(az)
    # origin = centre - L/2 (xhat + yhat), xhat = (sin a, cos a), yhat = (cos a, -sin a) in (E, N)
    de = -L / 2 * (math.sin(a) + math.cos(a))
    dn = -L / 2 * (math.cos(a) - math.sin(a))
    # base frame has az=0: x = north, y = east
    lon0, lat0 = base.to_lonlat(dn, de)
    grid = mininz.Grid(float(lon0), float(lat0), az, proj)
    return grid, float(lon0), float(lat0)


def _xy_from_en(params: dict, de: float, dn: float) -> tuple[float, float]:
    """Frame coordinates of the point (de, dn) metres east/north of the centre."""
    L = float(params["L"])
    a = math.radians(float(params.get("az", 0.0)))
    # Round away trig noise (cos 270 = -1.8e-16): SW4 snaps receivers with floor(),
    # so 5999.9999999 would land a whole cell away from 6000.
    x = L / 2 + de * math.sin(a) + dn * math.cos(a)
    y = L / 2 + de * math.cos(a) - dn * math.sin(a)
    return round(x, 6), round(y, 6)


def write_ascii_srf(path: Path, pts: list[mininz.SrfPoint]) -> None:
    lines = ["2.0", "PLANE 1"]
    p0 = pts[0]
    lines.append(f"{p0.lon:.8f} {p0.lat:.8f} 1 {len(pts)} 1.0 1.0")
    lines.append(f"{p0.stk:g} {p0.dip:g} {p0.dep_km:g} 0.0 0.5")
    lines.append(f"POINTS {len(pts)}")
    f32 = lambda v: float(np.float32(v))  # the SRF-HDF5 POINTS fields are float32
    for p in pts:
        lines.append(f"{f32(p.lon):.8f} {f32(p.lat):.8f} {f32(p.dep_km):.6f} {p.stk:g} {p.dip:g} {p.area_cm2:.6e} {p.tinit:g} {p.dt:g} "
                     f"{p.vs_cm_s:.6e} {p.den_g_cm3:g}")
        lines.append(f"{p.rake:g} {p.slip_cm:.6e} {len(p.slip_rate_cm_s)} 0.0 0 0.0 0")
        sr = list(p.slip_rate_cm_s)
        for k in range(0, len(sr), 6):
            lines.append(" ".join(f"{v:.8e}" for v in sr[k:k + 6]))
    path.write_text("\n".join(lines) + "\n")


def write_layered_sfile(path: Path, grid: mininz.Grid, L: float, Z: float, h: float, layers: list[dict],
                        az_offset: float = 0.0) -> None:
    """Flat sfile with one patch per constant-property layer.

    layers: [{"bottom": z, "vp":, "vs":, "rho":, "qp":, "qs":}, ...] top to bottom.
    """
    ni = int(round(L / h)) + 4
    with h5py.File(path, "w") as f:
        f.attrs["Origin longitude, latitude, azimuth"] = np.array([grid.lon0, grid.lat0, grid.az + az_offset],
                                                                  dtype=np.float64)
        f.attrs["Attenuation"] = np.int32(1)
        f.attrs["ngrids"] = np.int32(len(layers))
        f.attrs["Min, max depth"] = np.array([0.0, layers[-1]["bottom"]], dtype=np.float64)
        f.attrs["Coarsest horizontal grid spacing"] = np.float64(h)
        zi = f.create_group("Z_interfaces")
        zi.create_dataset("z_values_0", data=np.zeros((ni, ni), np.float32))
        mm = f.create_group("Material_model")
        top = 0.0
        for p, lay in enumerate(layers):
            zi.create_dataset(f"z_values_{p + 1}", data=np.full((ni, ni), lay["bottom"], np.float32))
            nk = max(2, int(round((lay["bottom"] - top) / h)) + 1)
            g = mm.create_group(f"grid_{p}")
            g.attrs["Horizontal grid size"] = np.float64(h)
            g.attrs["Number of components"] = np.int32(5)
            for name, key in (("Rho", "rho"), ("Cp", "vp"), ("Cs", "vs"), ("Qp", "qp"), ("Qs", "qs")):
                g.create_dataset(name, data=np.full((ni, ni, nk), lay[key], np.float32))
            top = lay["bottom"]


def _grid_moment_tensor(grid: mininz.Grid, lon: float, lat: float, stk: float, dip: float, rake: float) -> dict:
    """Unit moment tensor of a true-north strike/dip/rake in SW4's grid frame.

    The strike is converted to the grid with pyproj alone: true north at
    (lon, lat) has the grid direction (cos b, -sin b), b = atan2(-y_n, x_n),
    so the strike from the grid x-axis is stk - b (b = az + the projection's
    meridian convergence). Components follow Aki & Richards with x, y, z =
    grid x, grid y, down, as SW4's `source mxx=...` expects. Returned as
    gmxx..gmyz (strings; SW4 multiplies by m0), plus stk_grid and conv (b - az) in degrees.
    """
    x0, y0 = grid.to_xy(lon, lat)
    xn, yn = grid.to_xy(lon, lat + 1e-4)
    b = math.degrees(math.atan2(-(float(yn) - float(y0)), float(xn) - float(x0)))
    s_grid = stk - b
    S, D, R = (math.radians(v) for v in (s_grid, dip, rake))
    mxx = -(math.sin(D) * math.cos(R) * math.sin(2 * S) + math.sin(2 * D) * math.sin(R) * math.sin(S) ** 2)
    myy = math.sin(D) * math.cos(R) * math.sin(2 * S) - math.sin(2 * D) * math.sin(R) * math.cos(S) ** 2
    mxy = math.sin(D) * math.cos(R) * math.cos(2 * S) + 0.5 * math.sin(2 * D) * math.sin(R) * math.sin(2 * S)
    mxz = -(math.cos(D) * math.cos(R) * math.cos(S) + math.cos(2 * D) * math.sin(R) * math.sin(S))
    myz = -(math.cos(D) * math.cos(R) * math.sin(S) - math.cos(2 * D) * math.sin(R) * math.cos(S))
    conv = (b - grid.az + 180.0) % 360.0 - 180.0
    # 10 significant digits as strings: SW4 reads input lines into a 256-char
    # buffer and hangs on longer lines, and full float reprs overflow it.
    g = {"gmxx": mxx, "gmyy": myy, "gmzz": -(mxx + myy), "gmxy": mxy, "gmxz": mxz, "gmyz": myz}
    out = {k: f"{v:.10g}" for k, v in g.items()}
    out.update({"stk_grid": s_grid, "conv": conv})
    return out


def build(run_dir: Path, params: dict) -> dict:
    grid, lon0, lat0 = _frame(params)
    out: dict = {"lon0": lon0, "lat0": lat0}
    h = float(params["h"])
    mu = float(params["rho"]) * float(params["vs"]) ** 2
    # Point source at (src_de, src_dn) from the centre.
    sx, sy = _xy_from_en(params, float(params.get("src_de", 0.0)), float(params.get("src_dn", 0.0)))
    slon, slat = grid.to_lonlat(sx, sy)
    dep = float(params.get("src_dep_km", 2.0))
    sigma = float(params.get("sigma", 0.3))
    t0 = float(params.get("t0", 1.5))
    dt = float(params.get("srf_dt", 0.01))
    nt = int(round(2 * t0 / dt)) + 1
    slip, area = float(params.get("slip_cm", 50.0)), float(params.get("area_cm2", 1.0e10))
    pt = mininz.SrfPoint(lon=float(slon), lat=float(slat), dep_km=dep, stk=float(params.get("stk", 30.0)),
                         dip=float(params.get("dip", 70.0)), rake=float(params.get("rake", 45.0)), area_cm2=area,
                         slip_cm=slip, tinit=float(params.get("srf_tinit", 0.0)), dt=dt,
                         slip_rate_cm_s=mininz.gaussian_slip_rate(slip, dt, t0, sigma, nt))
    mininz.write_srf_hdf5(run_dir / "src.srf.h5", [pt])
    write_ascii_srf(run_dir / "src.srf", [pt])
    out.update({
        "src_lat": float(slat), "src_lon": float(slon), "src_dep": dep * 1000.0, "src_x": sx, "src_y": sy,
        "m0": mu * area * 1e-4 * slip * 1e-2, "src_t0": t0, "src_freq": 1.0 / sigma,
        "stk": pt.stk, "dip": pt.dip, "rake": pt.rake,
    })
    out.update(_grid_moment_tensor(grid, float(slon), float(slat), pt.stk, pt.dip, pt.rake))
    # Stations: (name, de, dn) offsets from the centre, metres.
    stations = params.get("stations", [["N", 0, 3000], ["E", 3000, 0], ["S", 0, -3000], ["W", -3000, 0],
                                       ["NE", 2000, 2000]])
    mode = params.get("station_mode", "latlon")
    rec_lines = []
    with h5py.File(run_dir / "stations.h5", "w") as f:
        for name, de, dn in stations:
            x, y = _xy_from_en(params, float(de), float(dn))
            g = f.create_group(name)
            g.create_dataset("ISNSEW", data=np.array([int(params.get("isnsew", 1))], np.int32))
            if mode == "xy":
                g.create_dataset("STX,STY,STZ", data=np.array([x, y, 0.0], np.float64))
                rec_lines.append(f"rec x={x:.6f} y={y:.6f} z=0 file=txt{name} usgsformat=1 sacformat=0 nsew=1")
            else:
                lon, lat = grid.to_lonlat(x, y)
                g.create_dataset("STLA,STLO,STDP", data=np.array([lat, lon, 0.0], np.float64))
                rec_lines.append(f"rec lat={float(lat):.10f} lon={float(lon):.10f} depth=0 file=txt{name} "
                                 f"usgsformat=1 sacformat=0 nsew=1")
    out["rec_lines"] = "\n".join(rec_lines)
    if params.get("layers"):
        write_layered_sfile(run_dir / "layers.sfile", grid, float(params.get("sfile_L", params["L"])),
                            float(params["Z"]), h, params["layers"], float(params.get("sfile_az_offset", 0.0)))
    if params.get("srf_drop"):
        with h5py.File(run_dir / "src.srf.h5", "a") as f:
            del f[params["srf_drop"]]
    return out
