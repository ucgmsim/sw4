"""A miniature of the SW4 problem ~/src/workflow generates for New Zealand events.

Same structure as `create-sw4-input` output, shrunk to run in seconds:

* NZTM-style `grid proj=tmerc ...` with an origin in Canterbury and azimuth 329.29
* an HDF5 sfile with topography, two material patches and attenuation (Q from
  the workflow's rule Qs = 0.05 Vs, Qp = 2 Qs)
* curvilinear top grid (mapping order 3) over two Cartesian refinement levels
  in the ratio 1:2:4, each band at least 12 cells thick
* a multi-point SRF v2.0 converted to SW4's SRF-HDF5 layout (SR1 only)
* an HDF5 station file with STLA,STLO,STDP (STDP=0, i.e. on the surface)

Everything here is deterministic.
"""

from dataclasses import dataclass, field
from pathlib import Path

import h5py
import numpy as np
from pyproj import Transformer

# SW4 projects from this geographic CRS (GeographicProjection.C).
SW4_GEOGRAPHIC_CRS = "+proj=latlong +datum=NAD83"


def sw4_proj_string(lon_p: float, lat_p: float, scale: float | None, ellps: str = "GRS80") -> str:
    """The PROJ string SW4 builds from `grid proj=tmerc ellps= lon_p= lat_p= scale=`."""
    s = f"+proj=tmerc +ellps={ellps} +lon_0={lon_p:g} +lat_0={lat_p:g}"
    if scale is not None:
        s += f" +scale={scale:g}"
    return s + " +units=m"


def nztm_like_string(lon_p: float, lat_p: float, scale: float, ellps: str = "GRS80") -> str:
    """What the workflow means by those parameters: NZTM without false origin."""
    return f"+proj=tmerc +ellps={ellps} +lon_0={lon_p:g} +lat_0={lat_p:g} +k={scale:g} +units=m"


@dataclass
class Grid:
    """SW4 grid frame: origin (lon, lat), x along azimuth az, y along az + 90."""

    lon0: float
    lat0: float
    az: float
    proj: str

    def __post_init__(self):
        self._fwd = Transformer.from_crs(SW4_GEOGRAPHIC_CRS, self.proj, always_xy=True)
        self._inv = Transformer.from_crs(self.proj, SW4_GEOGRAPHIC_CRS, always_xy=True)
        self._e0, self._n0 = self._fwd.transform(self.lon0, self.lat0)
        a = np.deg2rad(self.az)
        self._sin, self._cos = np.sin(a), np.cos(a)

    def to_xy(self, lon, lat):
        e, n = self._fwd.transform(lon, lat)
        e, n = np.asarray(e) - self._e0, np.asarray(n) - self._n0
        return e * self._sin + n * self._cos, e * self._cos - n * self._sin

    def to_lonlat(self, x, y):
        x, y = np.asarray(x, float), np.asarray(y, float)
        e = x * self._sin + y * self._cos
        n = x * self._cos - y * self._sin
        return self._inv.transform(e + self._e0, n + self._n0)


# --------------------------------------------------------------------------
# Material: a 1D profile hung from the topographic surface
# --------------------------------------------------------------------------
def vs_profile(depth_m):
    """Vs (m/s) as a function of depth below the surface: 500 m/s at the
    surface (the workflow's Vs floor) rising linearly to 3500 m/s at 6 km."""
    d = np.clip(np.asarray(depth_m, float), 0.0, None)
    return np.minimum(500.0 + 0.5 * d, 3500.0)


def material(depth_m):
    vs = vs_profile(depth_m)
    # Vp/Vs from 2.5 near the surface down to 1.75 at depth (inside NZCVM's 1.73..4)
    ratio = 1.75 + 0.75 * np.exp(-np.asarray(depth_m, float).clip(0) / 1500.0)
    vp = vs * ratio
    rho = 1700.0 + 0.2 * vs  # kg/m^3
    qs = 0.05 * vs  # workflow Q rule
    qp = 2.0 * qs
    return rho, vp, vs, qp, qs


def topography(x, y, hill_x, hill_y, amp=400.0, width=3000.0):
    """Elevation (m, positive up): a Gaussian hill plus a gentle tilt."""
    return amp * np.exp(-((x - hill_x) ** 2 + (y - hill_y) ** 2) / (2 * width**2)) + 0.004 * x


@dataclass
class Model:
    grid: Grid
    lx: float  # domain extent along x (m)
    ly: float
    lz: float
    hh: float  # sfile horizontal spacing (m)
    patch_bottoms: list[float]  # depth of the bottom of each sfile patch (m)
    nk: list[int]
    hill: tuple[float, float]
    elevation: bool = True
    pad: int = 3  # extra sfile cells beyond the domain on each high side
    attenuation: bool = True

    def surface(self, x, y):
        if not self.elevation:
            return np.zeros(np.broadcast(x, y).shape)
        return topography(x, y, *self.hill)


def write_sfile(path: Path, m: Model) -> dict:
    """Write an SW4 sfile (MaterialSfile.C layout): all patches share `hh`."""
    ni = int(np.ceil(m.lx / m.hh)) + 1 + m.pad
    nj = int(np.ceil(m.ly / m.hh)) + 1 + m.pad
    x = np.arange(ni) * m.hh
    y = np.arange(nj) * m.hh
    X, Y = np.meshgrid(x, y, indexing="ij")
    elev = m.surface(X, Y)
    tops = [-elev] + [np.full_like(elev, b) for b in m.patch_bottoms[:-1]]
    bottoms = [np.full_like(elev, b) for b in m.patch_bottoms]
    with h5py.File(path, "w") as f:
        f.attrs["Origin longitude, latitude, azimuth"] = np.array([m.grid.lon0, m.grid.lat0, m.grid.az], dtype=np.float64)
        f.attrs["Attenuation"] = np.int32(1 if m.attenuation else 0)
        f.attrs["ngrids"] = np.int32(len(m.patch_bottoms))
        f.attrs["Min, max depth"] = np.array([float(-elev.max()), float(m.patch_bottoms[-1])], dtype=np.float64)
        f.attrs["Coarsest horizontal grid spacing"] = np.float64(m.hh)
        zi = f.create_group("Z_interfaces")
        zi.create_dataset("z_values_0", data=tops[0].astype(np.float32))
        for p, b in enumerate(m.patch_bottoms):
            zi.create_dataset(f"z_values_{p + 1}", data=np.full((ni, nj), b, dtype=np.float32))
        mm = f.create_group("Material_model")
        for p, nk in enumerate(m.nk):
            g = mm.create_group(f"grid_{p}")
            g.attrs["Horizontal grid size"] = np.float64(m.hh)
            g.attrs["Number of components"] = np.int32(5 if m.attenuation else 3)
            t = np.linspace(0.0, 1.0, nk)
            z = tops[p][..., None] + (bottoms[p] - tops[p])[..., None] * t  # (ni, nj, nk) z coordinate
            depth = z + elev[..., None]  # below the local surface
            rho, vp, vs, qp, qs = material(depth)
            for name, arr in (("Rho", rho), ("Cp", vp), ("Cs", vs)) + ((("Qp", qp), ("Qs", qs)) if m.attenuation else ()):
                g.create_dataset(name, data=arr.astype(np.float32))
    return {"sfile_ni": ni, "sfile_nj": nj, "max_elevation": float(elev.max())}


# --------------------------------------------------------------------------
# Sources: SRF v2.0 in SW4's SRF-HDF5 layout (source_modelling.srf.write_sw4_hdf5)
# --------------------------------------------------------------------------
PLANE_DTYPE = np.dtype([("ELON", "<f4"), ("ELAT", "<f4"), ("NSTK", "<i4"), ("NDIP", "<i4"), ("LEN", "<f4"),
                        ("WID", "<f4"), ("STK", "<f4"), ("DIP", "<f4"), ("DTOP", "<f4"), ("SHYP", "<f4"),
                        ("DHYP", "<f4")])
POINT_DTYPE = np.dtype([("LON", "<f4"), ("LAT", "<f4"), ("DEP", "<f4"), ("STK", "<f4"), ("DIP", "<f4"),
                        ("AREA", "<f4"), ("TINIT", "<f4"), ("DT", "<f4"), ("VS", "<f4"), ("DEN", "<f4"),
                        ("RAKE", "<f4"), ("SLIP1", "<f4"), ("NT1", "<i4"), ("SLIP2", "<f4"), ("NT2", "<i4"),
                        ("SLIP3", "<f4"), ("NT3", "<i4")])


@dataclass
class SrfPoint:
    lon: float
    lat: float
    dep_km: float
    stk: float
    dip: float
    rake: float
    area_cm2: float
    slip_cm: float
    tinit: float
    dt: float
    slip_rate_cm_s: np.ndarray  # SR1 samples
    vs_cm_s: float = 3.0e5
    den_g_cm3: float = 2.7


def gaussian_slip_rate(slip_cm: float, dt: float, t0: float, sigma: float, n: int) -> np.ndarray:
    """Slip rate whose integral is slip_cm: a sampled Gaussian centred on t0."""
    t = np.arange(n) * dt
    g = np.exp(-0.5 * ((t - t0) / sigma) ** 2)
    return (slip_cm * g / (g.sum() * dt)).astype(np.float32)


def write_srf_hdf5(path: Path, points: list[SrfPoint], plane: dict | None = None) -> None:
    pts = np.zeros(len(points), dtype=POINT_DTYPE)
    sr = []
    for i, p in enumerate(points):
        pts[i] = (p.lon, p.lat, p.dep_km, p.stk, p.dip, p.area_cm2, p.tinit, p.dt, p.vs_cm_s, p.den_g_cm3,
                  p.rake, p.slip_cm, len(p.slip_rate_cm_s), 0.0, 0, 0.0, 0)
        sr.append(np.asarray(p.slip_rate_cm_s, np.float32))
    pl = np.zeros(1, dtype=PLANE_DTYPE)
    p0 = points[0]
    plane = plane or {}
    pl[0] = (plane.get("elon", p0.lon), plane.get("elat", p0.lat), plane.get("nstk", 1), plane.get("ndip", len(points)),
             plane.get("len", 1.0), plane.get("wid", 1.0), p0.stk, p0.dip, plane.get("dtop", p0.dep_km),
             0.0, plane.get("dhyp", 0.5))
    with h5py.File(path, "w") as f:
        f.attrs["VERSION"] = np.float32(2.0)
        f.attrs["PLANE"] = pl
        f.create_dataset("POINTS", data=pts)
        f.create_dataset("SR1", data=np.concatenate(sr))


def write_stations(path: Path, stations: dict[str, tuple[float, float]]) -> None:
    """Station file as workflow/scripts/generate_station_coordinates.py writes it."""
    with h5py.File(path, "w") as f:
        for name, (lat, lon) in stations.items():
            g = f.create_group(name)
            g.create_dataset("STLA,STLO,STDP", data=np.array([lat, lon, 0.0], dtype=np.float64))


# --------------------------------------------------------------------------
# The mini-NZ problem
# --------------------------------------------------------------------------
@dataclass
class MiniNZ:
    """Geometry in SW4 frame metres. `h` is the finest spacing (top grid)."""

    lon0: float = 172.30
    lat0: float = -43.70
    az: float = 329.29
    h: float = 200.0
    lx: float = 20000.0
    ly: float = 20000.0
    sponge: float = 3200.0
    elevation: bool = True
    lon_p: float = 173.0
    lat_p: float = 0.0
    scale: float = 0.9996
    stations: dict[str, tuple[float, float]] = field(default_factory=dict)  # name -> (x, y)

    @property
    def grid(self) -> Grid:
        return Grid(self.lon0, self.lat0, self.az, sw4_proj_string(self.lon_p, self.lat_p, self.scale))

    def bands(self) -> dict:
        """Depths in the workflow layout: curvilinear top, then h, 2h, 4h Cartesian bands of >= 12 cells."""
        h = self.h
        topo_zmax = 12 * h
        r1 = topo_zmax + 12 * h
        r2 = r1 + 12 * 2 * h
        lz = r2 + 8 * 4 * h
        return {"topo_zmax": topo_zmax, "r1": r1, "r2": r2, "lz": lz}


STATION_OFFSETS = {
    # name: (dx, dy) from the domain centre in the SW4 frame, metres
    "HILL": (-1000.0, 500.0),  # on the hilltop
    "NEAR": (1500.0, -1500.0),  # near the source
    "WEST": (-2500.0, 2500.0),
    "EAST": (2800.0, -2800.0),
    "FAR1": (-2800.0, -2800.0),
    "SPNG": (0.0, 0.0),  # x is set to the middle of the low-x sponge: SGDEPTH > 0
}


def build(run_dir: Path, params: dict) -> dict:
    """Write sfile, SRF and station file for a mini-NZ variant; return template params.

    The sponge is 12 coarse cells: with two refinement interfaces, 6 cells
    diverges, 10 is marginal (the refinement interface iteration stops
    converging) and 12 is clean; see tests/cases/stability/narrow-sponge and sources sit at least 5
    coarse cells inside it (the fork's margin_pts check), so the domain grows
    with the coarse spacing: lx = 2 * (sponge + 5 hc) + 6 km.
    """
    levels = int(params.get("levels", 3))  # 3 grids (1 Hz shape), 2 (0.5 Hz), 1 (0.25 Hz)
    hfine = float(params.get("hfine", 200.0))
    hc = hfine * 2 ** (levels - 1)
    sponge = float(params.get("sponge", 12 * hc))
    side = float(params.get("lx", 2 * (sponge + 5 * hc) + 6000.0))
    nz = MiniNZ(
        lon0=float(params.get("lon0", 172.30)), lat0=float(params.get("lat0", -43.70)),
        az=float(params.get("az", 329.29)), h=hfine, lx=side, ly=side, sponge=sponge,
        elevation=bool(params.get("elevation", True)), lon_p=float(params.get("lon_p", 173.0)),
    )
    b = nz.bands()
    # The bottom grid must hold the sponge (> sponge/hc points) plus slack.
    lz = max(b["lz"], b["r2"] + sponge + 4 * hc) if levels >= 2 else max(b["lz"], b["r1"] + sponge + 4 * hc)
    b["lz"] = lz
    grid = nz.grid
    cx0, cy0 = side / 2, side / 2
    hill = (cx0 - 1000.0, cy0 + 500.0)
    model = Model(grid=grid, lx=side, ly=side, lz=lz, hh=hfine,
                  patch_bottoms=[b["r1"], lz + 4 * hc], nk=[25, 33], hill=hill,
                  elevation=nz.elevation)
    info = write_sfile(run_dir / "mininz.sfile", model)

    # A 2 x 3 point fault, strike 45 (geographic), dip 60, at 3-4.2 km depth,
    # 1.5 km from the domain centre; tinit increases away from the first point.
    pts = []
    for i in range(2):
        for j in range(3):
            x, y = cx0 + 500.0 + 600.0 * i, cy0 - 1000.0 + 400.0 * j
            lon, lat = grid.to_lonlat(x, y)
            dep = 3.0 + 0.6 * j
            slip = 60.0 + 20.0 * i + 10.0 * j
            tinit = 0.25 * (i + j)
            pts.append(SrfPoint(lon=float(lon), lat=float(lat), dep_km=dep, stk=45.0, dip=60.0, rake=90.0 - 15 * i,
                                area_cm2=4.0e9, slip_cm=slip, tinit=tinit, dt=0.05,
                                slip_rate_cm_s=gaussian_slip_rate(slip, 0.05, 0.8, 0.25, 40)))
    write_srf_hdf5(run_dir / "mininz.srf.h5", pts)

    stations = {}
    xy = {}
    for name, (dx, dy) in STATION_OFFSETS.items():
        x, y = (cx0 + dx, cy0 + dy) if name != "SPNG" else (0.5 * sponge, cy0 + dy)
        lon, lat = grid.to_lonlat(x, y)
        lon = ((float(lon) + 180.0) % 360.0) - 180.0
        stations[name] = (float(lat), lon)
        xy[name] = (x, y)
    write_stations(run_dir / "stations.h5", stations)

    refinements = []
    hcoarse = hc
    if levels >= 3:
        refinements.append(f"refinement zmax={b['r1']:.1f}")
    if levels >= 2:
        refinements.append(f"refinement zmax={b['r2']:.1f}")
    out = {
        "lon0": nz.lon0, "lat0": nz.lat0, "az": nz.az, "lx": side, "ly": side, "lz": lz,
        "hcoarse": hcoarse, "topo_zmax": b["topo_zmax"],
        "refinements": "\n".join(refinements),
        "sponge": nz.sponge,
        "topography": (f"topography input=sfile zmax={b['topo_zmax']:.1f} order=3 file=mininz.sfile"
                       if nz.elevation else ""),
        "station_xy": xy, "stations": stations,
        "srf_points": [{"lon": p.lon, "lat": p.lat, "dep_km": p.dep_km} for p in pts],
        **info,
    }
    return out
