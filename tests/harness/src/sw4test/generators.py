"""Input asset generators (`generate = "<name>"` in a case or run).

A generator writes files into the run directory and returns extra template
parameters. Generators must be deterministic so goldens stay valid.
"""

from pathlib import Path
from typing import Any, Callable

GENERATORS: dict[str, Callable[[Path, dict[str, Any], Path], dict[str, Any]]] = {}


def generator(name: str):
    def deco(fn):
        GENERATORS[name] = fn
        return fn

    return deco


def generate(name: str, run_dir: Path, params: dict[str, Any], case_dir: Path) -> dict[str, Any]:
    if name not in GENERATORS:
        raise ValueError(f"unknown generator {name!r}; known: {sorted(GENERATORS)}")
    return GENERATORS[name](run_dir, params, case_dir)


@generator("mininz")
def _mininz(run_dir: Path, params: dict[str, Any], case_dir: Path) -> dict[str, Any]:
    from . import mininz

    return mininz.build(run_dir, params)


@generator("eq")
def _eq(run_dir: Path, params: dict[str, Any], case_dir: Path) -> dict[str, Any]:
    from . import eqgen

    return eqgen.build(run_dir, params)
