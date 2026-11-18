"""PEER LOH.1 / LOH.3 reference solutions (Day et al. 2001), from the PROSE
frequency-wavenumber files bundled in tools/ (Apsel & Luco 1983).

Port of tools/ReadUHS.m, tools/loh1exact.m and tools/loh3exact.m: the files
hold station-10 velocities for PROSE's own source time function; two
convolutions replace it with SW4's `type=Gaussian` moment rate of spread
`sigma`, centred at t0 = 6 sigma.
"""

from pathlib import Path

import numpy as np


def read_prose(path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """-> t, vertical, radial, transverse (scaled as in ReadUHS.m)."""
    a = np.loadtxt(path)
    return a[:, 0], a[:, 1] * -1e5, a[:, 2] * 1e5, a[:, 3] * 1e5


def exact(path: Path, sigma: float) -> dict[str, np.ndarray]:
    """Reference radial/transverse/vertical velocities for a Gaussian of spread
    sigma. Vertical is positive down (SW4's z), as loh3exact.m notes."""
    t, ve1, ra1, tr1 = read_prose(path)
    nt = len(t)
    dt = t[1] - t[0]
    T = 0.1
    f1 = (1 / T**2) * t * np.exp(-t / T)
    ts = 6 * sigma
    tau = t - ts
    factor = 1 - (2 * T / sigma**2) * tau - ((T / sigma) ** 2) * (1 - (tau / sigma) ** 2)
    f2 = (1 / np.sqrt(2 * np.pi) / sigma) * factor * np.exp(-0.5 * (tau / sigma) ** 2)
    out = {}
    for name, x in (("radial", ra1), ("transverse", tr1), ("vertical", ve1)):
        y = dt * np.convolve(x, f1)[:nt]
        y = dt * np.convolve(y, f2)[:nt]
        out[name] = y
    out["t"] = t
    return out
