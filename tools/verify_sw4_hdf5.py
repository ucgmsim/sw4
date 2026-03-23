"""Verify that write_sw4_hdf5 produces output identical to SW4's srf2hdf5.py.

Uses rtest.srf from this SW4 repo as a Version 1.0 SRF test case.

Usage
-----
    python tools/verify_sw4_hdf5.py

Requires the source_modelling package to be installed (pip install source_modelling).
"""

import sys
import tempfile
from pathlib import Path

from source_modelling import srf

SW4_REPO = Path(__file__).parent.parent


def run_srf2hdf5(srf_path: Path, output_path: Path) -> None:
    """Run SW4's srf2hdf5.py as a subprocess."""
    import subprocess

    script = SW4_REPO / "tools" / "srf2hdf5.py"
    result = subprocess.run(
        [sys.executable, str(script), str(srf_path), str(output_path)],
        capture_output=True,
        text=True,
    )
    if result.returncode != 0:
        print(result.stdout)
        print(result.stderr)
        raise RuntimeError(f"srf2hdf5.py failed with return code {result.returncode}")


def compare(ref_path: Path, ours_path: Path) -> bool:
    import hashlib

    ref_md5 = hashlib.md5(ref_path.read_bytes()).hexdigest()
    ours_md5 = hashlib.md5(ours_path.read_bytes()).hexdigest()
    print(f"srf2hdf5.py:     {ref_md5}")
    print(f"write_sw4_hdf5:  {ours_md5}")
    return ref_md5 == ours_md5


def main():
    srf_path = SW4_REPO / "examples" / "rupture" / "rtest.srf"

    print(f"SRF file: {srf_path}")

    with tempfile.TemporaryDirectory() as tmp_dir:
        ref_path = Path(tmp_dir) / "ref.h5"
        ours_path = Path(tmp_dir) / "ours.h5"

        print("\n--- Running SW4's srf2hdf5.py ---")
        run_srf2hdf5(srf_path, ref_path)

        print("\n--- Running write_sw4_hdf5 ---")
        srf_file = srf.read_srf(srf_path)
        srf_file.write_sw4_hdf5(ours_path)
        print(f"Written: {len(srf_file.header)} plane(s), {len(srf_file.points)} points")

        print("\n--- Comparing outputs ---")
        ok = compare(ref_path, ours_path)

    if ok:
        print("\nChecksums match.")
    else:
        print("\nChecksums differ.")
        sys.exit(1)


if __name__ == "__main__":
    main()
