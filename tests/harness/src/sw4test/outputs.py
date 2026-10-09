"""Readers for the files SW4 writes."""

import re
from pathlib import Path

import numpy as np


def read_twilight_err(path: Path) -> dict[str, list[float]]:
    """TwilightErr.txt -> {"disp": [errInf, errL2, solInf], "atten": [...]}."""
    lines = path.read_text().splitlines()
    out: dict[str, list[float]] = {}
    for i, line in enumerate(lines):
        key = "disp" if line.startswith("Displacement") else "atten" if line.startswith("Attenuation") else None
        if key and i + 1 < len(lines):
            out[key] = [float(x) for x in lines[i + 1].split()]
    if "disp" not in out:
        raise ValueError(f"{path}: no 'Displacement variables' block")
    return out


def read_err_series(path: Path) -> np.ndarray:
    """LambErr/PointSourceErr/RayleighErr.txt -> (n, 4) array of t, errInf, errL2, solInf."""
    rows = []
    for line in path.read_text().splitlines():
        parts = line.split()
        if len(parts) == 4:
            try:
                rows.append([float(p) for p in parts])
            except ValueError:
                continue
    if not rows:
        raise ValueError(f"{path}: no numeric rows")
    return np.array(rows)


def read_energy(path: Path) -> np.ndarray:
    return np.array([float(x) for x in path.read_text().split()])


def read_usgs(path: Path) -> tuple[np.ndarray, np.ndarray, list[str]]:
    """SW4 `rec usgsformat=1` text file -> (t, data[n, 3], component names)."""
    header: list[str] = []
    rows = []
    for line in path.read_text().splitlines():
        if line.startswith("#"):
            header.append(line)
            continue
        parts = line.split()
        if parts:
            rows.append([float(p) for p in parts])
    arr = np.array(rows)
    names = []
    for h in header:
        m = re.match(r"#\s*Column\s+(\d+):\s*(.*)", h)
        if m and int(m.group(1)) > 1:
            names.append(m.group(2).strip())
    return arr[:, 0], arr[:, 1:], names
