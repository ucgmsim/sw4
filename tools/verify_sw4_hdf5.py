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

import h5py
import numpy as np

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
    all_ok = True

    def fail(msg):
        nonlocal all_ok
        all_ok = False
        print(f"FAIL: {msg}")

    with h5py.File(ref_path, "r") as ref, h5py.File(ours_path, "r") as ours:
        # VERSION
        if ref.attrs["VERSION"] != ours.attrs["VERSION"]:
            fail(f"VERSION: ref={ref.attrs['VERSION']}, ours={ours.attrs['VERSION']}")
        else:
            print(f"OK   VERSION = {ours.attrs['VERSION']}")

        # PLANE
        ref_plane = ref.attrs["PLANE"]
        ours_plane = ours.attrs["PLANE"]
        if ref_plane.shape != ours_plane.shape:
            fail(f"PLANE shape: ref={ref_plane.shape}, ours={ours_plane.shape}")
        else:
            for field in ref_plane.dtype.names:
                ref_vals = np.array([ref_plane[i][field] for i in range(len(ref_plane))])
                ours_vals = np.array([ours_plane[i][field] for i in range(len(ours_plane))])
                if np.array_equal(ref_vals, ours_vals):
                    print(f"OK   PLANE.{field}")
                else:
                    fail(f"PLANE.{field}: max diff = {np.max(np.abs(ref_vals - ours_vals))}")

        # POINTS
        ref_pts = ref["POINTS"]
        ours_pts = ours["POINTS"]
        if ref_pts.shape != ours_pts.shape:
            fail(f"POINTS shape: ref={ref_pts.shape}, ours={ours_pts.shape}")
        else:
            for field in ref_pts.dtype.names:
                ref_vals = ref_pts[field][:]
                ours_vals = ours_pts[field][:]
                if np.array_equal(ref_vals, ours_vals):
                    print(f"OK   POINTS.{field}")
                else:
                    mismatches = np.sum(ref_vals != ours_vals)
                    max_diff = np.max(np.abs(ref_vals.astype(np.float64) - ours_vals.astype(np.float64)))
                    fail(f"POINTS.{field}: {mismatches}/{len(ref_vals)} mismatches, max diff = {max_diff}")

        # SR1
        if "SR1" in ref and "SR1" in ours:
            ref_sr1 = ref["SR1"][:]
            ours_sr1 = ours["SR1"][:]
            if ref_sr1.shape != ours_sr1.shape:
                fail(f"SR1 shape: ref={ref_sr1.shape}, ours={ours_sr1.shape}")
            elif np.array_equal(ref_sr1, ours_sr1):
                print(f"OK   SR1 ({len(ref_sr1)} values)")
            else:
                mismatches = np.sum(ref_sr1 != ours_sr1)
                max_diff = np.max(np.abs(ref_sr1.astype(np.float64) - ours_sr1.astype(np.float64)))
                fail(f"SR1: {mismatches}/{len(ref_sr1)} mismatches, max diff = {max_diff}")
        elif "SR1" not in ref and "SR1" not in ours:
            print("OK   SR1 (both absent)")
        else:
            fail(f"SR1 presence: ref={'yes' if 'SR1' in ref else 'no'}, ours={'yes' if 'SR1' in ours else 'no'}")

    return all_ok


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
        print("\nAll fields match exactly.")
    else:
        print("\nSome fields differ.")
        sys.exit(1)


if __name__ == "__main__":
    main()
