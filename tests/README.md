# SW4 test suite

Every test is a ctest test. A case is a TOML file under `cases/<group>/`. The
harness in `harness/` (a uv project, the `sw4test` package) turns each case's
runs and checks into ctest tests:

- `<group>/<case>:run:<run>` runs SW4 once. Each run is a ctest fixture.
- `<group>/<case>:check:<check>` asserts something about one or more runs and
  requires their fixtures, so ctest runs them first.

## Running

```sh
cmake -S . -B build -G Ninja -DUSE_HDF5=ON -DUSE_PROJ=ON   # needs `uv` on PATH
cmake --build build
ctest --test-dir build -L '^pr$' -j "$(nproc)"     # pull-request tier
ctest --test-dir build -j "$(nproc)"                # everything (nightly tier)
ctest --test-dir build -R workflow/mininz           # one case, with its runs
```

Labels are regular expressions, so anchor them: `-L pr` also matches `proj`.
Each run declares its MPI ranks as ctest `PROCESSORS`, so `-j N` never runs
more than N ranks at once. The harness passes `--oversubscribe --bind-to none`
to Open MPI so that concurrent tests do not all pin to the same cores, and sets
`OMP_NUM_THREADS` from the run (default 1).

Outputs go to `build/tests/out/<group>/<case>/<run>/`: the rendered input
`sw4.in`, `sw4.out`/`sw4.err`, `params.json` and SW4's own output files.
To rerun a check without rerunning SW4:

```sh
uv run --project tests/harness sw4test check --config build/tests/sw4test.json workflow/mininz images-contract
```

## Writing a case

```toml
description = "What is tested and why; measured values for tolerances."
labels = ["mms", "curvi"]          # extra ctest labels
input = "topo.in"                  # SW4 input template, {{name}} placeholders
generate = "mininz"                # optional: writes input files (sw4test/generators.py)
requires = ["hdf5", "proj"]        # build features the case needs
precisions = ["double", "float"]   # default: both
ranks = 4                          # default 4; per-run override
timeout = 300

[params]                           # template values shared by all runs
t = 0.8

[[run]]
name = "h1"
tier = "pr"                        # "pr" (default) or "nightly"
params = { nx = 31, h = 0.0333 }

[[check]]
name = "rate"
kind = "rate"                      # see sw4test/checks*.py
runs = ["h1", "h2"]                # default: all runs
expect = 4.0
```

A run can expect SW4 to fail (`expect = "failure"`, with `stdout_regex`), wait
for another run (`after = [...]`) and copy files from it (`copy_from`).

**Known defects** are marked `xfail = "reason"` on a run or check. The test
passes while the defect is present (XFAIL) and fails once it is fixed (XPASS),
so the marker gets removed and the test starts guarding the fix.

## Goldens

Checks such as `twilight_golden`, `err_golden` and `traces_golden` compare with
`golden/<precision>/<case>/<check>.json`. Each golden records the commit and
build it came from. Regenerate a golden only for an intended change, and review
the diff like code:

```sh
ctest --test-dir build -R 'mms/flat:run'
uv run --project tests/harness sw4test bless --config build/tests/sw4test.json 'mms/flat'
```

Tolerances are per precision: `rtol = { double = 1e-6, float = 1e-3 }`.

## Groups

| Group | What it checks |
|---|---|
| `mms` | Order of accuracy with manufactured solutions (twilight): Cartesian, supergrid, topography, refinement, the workflow's curvilinear and two-refinement stack, attenuation, high Vp/Vs, 2nd-order time, rank invariance |
| `analytic` | Lamb's problem, a whole-space point source, a source near a refinement interface |
| `energy` | Discrete energy of random data: conserved or decaying on the workflow grid stack, with attenuation and supergrid; the attenuation CFL regression (6f62621b) |
| `workflow` | Mini-NZ, a scaled-down copy of the input ~/src/workflow generates: the rechdf5 and imagehdf5 contract, projection, golden traces, float vs double |
| `equivalence` | One physical problem posed in different ways: SRF vs `source`, SRF HDF5 vs ASCII, slip-rate output is velocity, sfile vs `block`, flat curvilinear vs Cartesian, rechdf5 vs text, `width` vs `gp`, mirror symmetry, grid rotation, NS/EW orientation |
| `supergrid` | Reflection from the absorbing layer |
| `prefilter` | The prefilter is a zero-phase Butterworth; precursor truncation |
| `parallel` | Bitwise identical results across MPI ranks and threads; checkpoint and restart |
| `stability` | Long-time stability at the default CFL; narrow sponges with refinement |
| `negative` | Inputs SW4 must reject or warn about |
| `benchmark` | LOH.1 and LOH.3 against the PROSE references in `tools/` (nightly) |

## CI

`.github/workflows/tests.yml` runs the `pr` tier in double and float on every
pull request, then compares the float and double mini-NZ results
(`sw4test crosscheck`). The nightly run executes every tier and
~/src/workflow's end-to-end test with this SW4.
