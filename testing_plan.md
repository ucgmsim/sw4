# SW4 test expansion plan

**Goal:** check SW4 completely on every pull request, with emphasis on the features `~/src/workflow` actually uses. Every test runs under `ctest`, and a GitHub Actions workflow gates PRs on them.

Status: draft, 2026-10-08. Branch `single_precision`.

---

## 0. Implementation status (2026-10-08, branch `test-suite`)

Delivered (see `tests/README.md`):

- **Harness:** `tests/harness`, a uv project (`sw4test`). Cases are TOML under `tests/cases`, registered with ctest as run fixtures plus checks, labelled `pr` or `nightly`. It supports per-precision goldens, xfail markers for known defects, and a `crosscheck` command for float vs double.
- **Size:** 45 cases. In double that is **239 ctest tests: 164 in the PR tier**, which takes about 10 min on 6 local cores. The old 24-case ctest block is replaced. Phase 0 fixes to `check_results.py` and `test_sw4.py` are included.
- **CI:** `.github/workflows/tests.yml` runs the PR tier in double and float, the float-vs-double crosscheck, and nightly every tier plus the workflow e2e test with the real SW4 (needs the `EMOD3D_TOKEN` secret).
- **Workflow e2e with real SW4 (XR-01, run locally):** it fails. The e2e test's shrunk geometry (`nz_min=4`, 1 km grid) makes `create-sw4-input` emit `refinement zmax=2000` inside `topography zmax=6167`, giving a 5-point grid that real SW4 rejects. The fake SW4 stub accepts it. This needs fixing in ~/src/workflow.
- **testenergy:** the `rhobase/mubase/qsbase/basement` options from `testenergy-master.patch` are merged with the defaults unchanged. EN-05 uses them, and it fails on the build before 6f62621b.

Plan items mapped to what was built:

| Plan | Status |
|---|---|
| MMS-01..09, 11, 12 | Done (`mms/`). MMS-10 dropped: SW4 twilight forcing supports nmech=1 only. Topography `order=3` (the workflow's setting) converges at about 3rd order: `mms/topo-order3`. |
| AN-01..04 | Done (`analytic/`); AN-03 is nightly (needs h=0.05). AN-05 (Rayleigh) not done. |
| EN-01..07 | Done (`energy/`). EN-03 is Cartesian only, because the curvilinear energy diagnostic omits supergrid stretching. EN-06 (float long run) is covered by the nightly `energy-long` checks. |
| WF-01..10 | Done (`workflow/mininz`); WF-10 is the nightly `1hz-full` run. |
| EQ-01..07, 09, 10 (+ EQ-11 orientation) | Done (`equivalence/`). EQ-08 (reciprocity) not done. |
| SG-01..04 | SG-01/02 done (`supergrid/reflection`). SG-04 covered by `negative/inputs` (source in sponge) and `stability/narrow-sponge`. |
| PF-01/02 (+ PF-03 precursor) | Done (`prefilter/response`). |
| AT-01/02 (Q fit, plane-wave decay) | Not done. They need C++ unit tests or a 1D setup. AT-03 (LOH.3) is done. |
| PI-01..03, PX-01/02 | Done (`parallel/`, `workflow/mininz:float-vs-double`, per-precision goldens). |
| BM-01 (+ AT-03) | Done (`benchmark/loh1`, `loh3`, nightly, against the PROSE references). BM-02..05 not done. |
| Negative tests | Done (`negative/inputs`, 10 runs). |
| IO-01..03 | Not done (HDF5 image writer race). |
| C++ unit tests | Not done. |

Known defects encoded as xfail (each turns red when fixed, so the marker gets removed):
NS/EW orientation and skew; the az=180/270 rotation checks; receiver floor snapping; `grid scale=` ignored by PROJ; the prefilter precursor truncated at t=0; failonnan missing Inf; default CFL with a 12-point sponge; SRF missing SR1/POINTS or outside the grid not stopping the run; an sfile smaller than the grid silently clamped; text `rec` leaking empty groups into rechdf5; the uninitialised topo image without topography.

---

## 1. Where we are today

### 1.1 Existing tests

| Item | State | Problems |
|---|---|---|
| ctest (`CMakeLists.txt:286-428`) | 24 cases (twilight, attenuation, meshrefine, lamb, pointsource) × Run + Check | Check uses `DEPENDS ${Test_Name}` rather than `Run_${Test_Name}` (l.425), so ordering is not enforced. No labels, fixtures or timeouts. Float tolerance is 5e-1 (l.348-358), which is close to no check at all. |
| `pytest/check_results.py` | compare mode | For Lamb and PointSource it compares only columns 0–1, and column 0 is *time*, so errL2 and solInf are never checked. `max(t,b)` without `abs`. Undefined `printf` at l.84 turns an exception into a NameError. |
| `pytest/test_sw4.py` | 26 entries: twilight, attenuation, meshrefine, curvimeshrefine, energy, lamb, pointsource, hdf5, geodynbc, supergrid | `main_test` returns True even when tests fail (l.410-412), so it exits 0. MPI launcher is chosen by hostname. Not precision-aware: energy check is 1e-10 against float energy sums. Not in ctest at all. |
| curvimeshrefine, energy, hdf5, geodynbc, supergrid | Only in pytest | Never run under ctest or in PR CI. |
| `pytest/scec`, `pytest/loh3` | Inputs only | Never run. |
| C++ unit tests | None | — |
| CI `compilers.yml` | gcc/clang/intel/aocc × double/float, HDF5+FFTW+PROJ | **Builds only; runs no tests.** |
| CI `linux.yml` | Makefile build + `test_sw4.py` on master/developer | Cannot fail because of the exit-code bug. Does not trigger on `single_precision` or on PRs. |
| Convergence-rate checks | None | The reference ladders do show 4th order (flat-twi errInf 6.85e-4 → 3.99e-5 → 2.01e-6), but nothing asserts it. |

**What this means in practice:** a PR can currently break the numerics, the HDF5 I/O, or every feature the workflow uses, and CI stays green. The tests that do exist mostly check that an error norm is unchanged against a stored reference. They do not check that the answer is correct.

### 1.2 How the workflow drives SW4 (the contract we must protect)

Sources: `workflow/scripts/sw4_template.py`, `workflow/default_parameters/v26_7_*Hz/defaults.yaml`, `workflow/scripts/lf_to_xarray.py`.

| Feature | What the workflow emits or depends on |
|---|---|
| `fileio` | `verbose=2 printcycle=10` |
| `grid` | `proj=tmerc ellps=GRS80 lon_p=173 lat_p=0 scale=0.9996`, plus `x y z h az lat lon`. az is about 329°; x is north. Domains can cross the 180° antimeridian. |
| `time` | `t=<duration>` (about 100 s) |
| `refinement zmax=` | 1Hz: 100 m to 5 km, 200 m to 25 km, 400 m below. 0.5Hz: 2 levels. 0.25Hz: none. |
| `topography` | `input=sfile order=3 zmax=3·e_max`. Each layer has at least `nz_min=12` cells. |
| `sfile` | NZCVM HDF5 sfile with attenuation (Qs=0.05·Vs, Qp=2·Qs) and multiple Z interfaces |
| `attenuation` | `nmech=3 phasefreq=0.5 maxfreq=10` |
| `supergrid` | `width=12000` (metres, not `gp`). The sponge sits inside the bottom refinement. Sources must be ≥ 5 coarse cells from the sponge. |
| `prefilter` | `type=lowpass order=4 passes=2 fc2≈1.116/0.558/0.279` |
| `developer` | `cfl=0.9 reporttiming=1 failonnan=1` |
| `rupturehdf5` | SRF v2.0 → HDF5 (`PLANE` attribute, `POINTS`, `SR1`; SLIP2/3 = 0). This is the **only** source path; point sources are 1-point SRFs. |
| `rechdf5` | HDF5 station file (STLA/STLO, STDP=0, so topodepth puts the receiver at the surface). Defaults: displacement mode, NSEW, writeEvery=1000. |
| `imagehdf5` (×10) | z=0, float. topo/grid/p/s/rho at cycle 0. mag/velmag/uz every 0.5 s. hmax/vmax at `time=t_end`. |

What downstream code (`lf_to_xarray`, `bb_sim`, `im_calc`) relies on in the output:

- `DELTA` is correct, and `NPTS` is the same for every station.
- `NS`/`EW`/`UP` are h5py *datasets* (`USE_DSET_ATTR`, `src/sachdf5.C:53`), with UP positive up.
- `STLA`/`STLO` are echoed back.
- `SGWIDTH`/`SGWIDTHGP`/`SGDEPTH`/`SGDEPTHGP` (fork-only).
- The time axis starts at t=0.
- **"Displacement" output from an SRF slip-rate source is really ground velocity in m/s** (`lf_to_xarray.py:~61`).

Features the workflow does **not** use: `source`, `rupture` (ASCII), `block`/`rfile`/`pfile`/`efile`/`ifile`, `rec`/`sac` text output, `image`/`volimage`, `checkpoint`, `geodynbc`, anisotropy, sw4mopt. These keep their current coverage but get no new investment beyond the basic smoke tests in §3.9.

### 1.3 Known open defects the suite must cover (from project memory)

- **Curvilinear interface non-convergence in float** (`CurvilinearInterface2::impose_ic`). A production float run blew up at about 7470 steps (`convergence_plan.md`).
- **Attenuation + high CFL instability.** Confirmed with an energy test at CFL 1.3; CFL ≤ 1.1 is stable. The workflow uses 0.9, but the margin is not tested.
- **HDF5 image writer race.** Large z-plane images abort with lock errno 11.
- Two float cases with suspiciously large errors: `attenuation/tw-att-2` and `meshrefine/refine-el-1`, both about 3.6e-1 relative.

---

## 2. Testing principles (from the literature)

1. **Verification, not just regression.**
   - Use the Method of Manufactured Solutions (SW4 calls it "twilight") with *observed order of accuracy* assertions: Salari & Knupp 2000 (SAND2000-1444); Roache 2002; Oberkampf & Roy 2010.
   - A stored-norm comparison only says "unchanged". p_obs ≈ 4 says "correct".
   - Formal orders to assert: 4 for smooth interior, curvilinear, refinement and attenuation (Sjögreen & Petersson 2012; Zhang, Wang & Petersson 2021); **2 for singular point sources** (Petersson & Sjögreen 2010).
2. **Energy as an oracle.**
   - The SBP schemes are provably energy-conserving (elastic, no damping) or energy-dissipating (supergrid or attenuation): Nilsson et al. 2007; Petersson & Sjögreen 2012, 2014.
   - Random initial data excites grid-scale modes and finds slow instabilities that smooth sources hide. The attenuation/CFL instability in §1.3 was found exactly this way.
3. **Metamorphic and invariance tests where no exact solution exists.** This addresses the "oracle problem" in Kanewala & Bieman 2014. Examples:
   - MPI decomposition invariance
   - restart equivalence
   - two input paths that must give the same physics (sfile vs block, 1-point SRF vs `source`)
   - reciprocity (Eisner & Clayton 2001)
   - symmetry
4. **Community benchmarks for validation of the assembled code.** LOH.1 and LOH.3 (Day et al. 2001/2003) against frequency-wavenumber references, scored with Kristeková et al. 2006/2009 time-frequency envelope and phase misfits rather than raw RMS.
5. **Each precision has its own tolerances and goldens.** Assert convergence only while the error is about 100× above the round-off floor. Keep an explicit float-vs-double cross-check (single precision is the production configuration here).
6. **Small and fast in PR, broad nightly.** Coarse grids still exercise every code path (Dubey et al. FLASH process; xSDK policies). Two-level convergence runs per PR; three or four levels and the long energy runs nightly.

---

## 3. Test catalogue

Every test is a ctest test. Each has `LABELS` from {`pr`, `nightly`}, a feature label (`twilight`, `curvi`, `mr`, `atten`, `sfile`, `rupture`, `rechdf5`, `image`, `proj`, `energy`, `mpi`, `restart`, `float-xcheck`, `workflow`, `negative`, `unit`), and `TIMEOUT`, `PROCESSORS` and `RESOURCE_LOCK`/`FIXTURES_*` where needed.

Columns:

- **WF** = relevance to the workflow: ★★★ direct contract, ★★ underlying numerics it relies on, ★ general.
- **Tier**: when the test runs.

### 3.1 Order-of-accuracy verification (MMS / twilight)

Each case runs at h and h/2 in PR and at h, h/2, h/4 nightly, with CFL fixed at the workflow value of 0.9. Assertions:

- (a) error norms within tolerance of the golden values, for each precision;
- (b) `p_obs` between the two finest levels in [3.6, 4.4] (L2 and L∞; [1.7, 2.3] where 2nd order is expected);
- (c) memory-variable errors also converge when attenuation is on;
- (d) stdout contains no "no convergence in curvilinear interface".

| ID | Configuration | WF | Tier |
|---|---|---|---|
| MMS-01 | Cartesian, free surface top, Dirichlet sides | ★★ | pr |
| MMS-02 | MMS-01 + `sgstretching` supergrid on 4 sides + bottom | ★★★ | pr |
| MMS-03 | Gaussian topography, curvilinear, `order=3` (workflow's order) | ★★★ | pr |
| MMS-04 | One refinement interface, Cartesian | ★★★ | pr |
| MMS-05 | **Two** refinement interfaces (1:2:4 ratio as in 1Hz) | ★★★ | pr |
| MMS-06 | Curvilinear + 2 refinements (the full 1Hz stack shape) | ★★★ | pr |
| MMS-07 | MMS-06 + attenuation nmech=3 | ★★★ | pr |
| MMS-08 | MMS-07 + supergrid stretching | ★★★ | pr |
| MMS-09 | Heterogeneous material with high Vp/Vs (≈4, the NZCVM max) | ★★ | pr |
| MMS-10 | Attenuation nmech=1,3,5 sweep | ★★ | nightly |
| MMS-11 | 2nd-order time stepping (`time_order=2`): asserts p≈2 in time | ★ | nightly |
| MMS-12 | Rank-count sweep of MMS-06 at 1/3/4 ranks (decomposition across interfaces) | ★★ | pr |

The existing twilight/attenuation/meshrefine/curvimeshrefine cases map onto MMS-01, 03, 04, 07 and so on. They become rate tests instead of being duplicated.

### 3.2 Analytic solutions

| ID | Test | Assertion | WF | Tier |
|---|---|---|---|---|
| AN-01 | Lamb's problem (`testlamb`), 2 levels | errInf/errL2 golden values; p_obs ≥ 3.5 (surface) | ★★ | pr |
| AN-02 | Whole-space point source (`testpointsource`), smooth time function | p_obs ≈ 2 for the singular source (Petersson & Sjögreen 2010) | ★★ | pr |
| AN-03 | AN-02 with the source placed 0, ½ and 1 cells from a refinement interface | Error independent of position within 2× | ★★★ (fault planes cross refinement boundaries) | pr |
| AN-04 | AN-02 at 3 levels | p_obs ≈ 2 | ★★ | nightly |
| AN-05 | Rayleigh wave (`testrayleigh`) | Golden error norms; p_obs | ★★ | nightly |

### 3.3 Energy and stability

These use `testenergy` (random data) and must log energy in **double**, even in a float build. This needs a small code change: accumulate `compute_energy` in double. Assertions:

- monotone non-increase (or conservation to drift ε) with a precision-aware ε;
- no NaN;
- for the "unstable" cases, a growth rate > 0 detected. These are expected-failure tests that prove the detector works.

| ID | Configuration | Steps (PR / nightly) | WF | Tier |
|---|---|---|---|---|
| EN-01 | Cartesian, no damping: conservation | 2k / 20k | ★★ | pr |
| EN-02 | Curvilinear + 2 refinements, CFL 0.9 | 2k / 20k | ★★★ | pr |
| EN-03 | EN-02 + supergrid: monotone decay | 2k / 20k | ★★★ | pr |
| EN-04 | EN-02 + attenuation nmech=3 (needs the `testenergy` attenuation patch from the CFL experiment) | 2k / 20k | ★★★ | pr |
| EN-05 | CFL margin sweep with attenuation: 0.9, 1.0, 1.1 must be stable; 1.3 is recorded as known-unstable (`WILL_FAIL` until fixed, then flipped) | — / 20k | ★★★ | nightly |
| EN-06 | Float build: EN-02/03 run to 20k steps with no curvilinear-interface non-convergence warning | — / 20k | ★★★ | nightly |
| EN-07 | Existing `energy/*` and `curvimeshrefine/energy-1`, ported to ctest | as now | ★★ | pr |

### 3.4 Workflow-shaped integration test: "mini-NZ"

This is the centrepiece. It is a tiny problem generated to match `sw4_template.py`'s output **exactly in structure**, with shrunk numbers.

**Inputs**, built by `tests/sw4test/make_mininz.py` (h5py + numpy, deterministic seed) at configure time, or committed (< 2 MB):

- **sfile:** 3 grids, attenuation on, Gaussian-hill topography with a max elevation of about 500 m, layered Vs 500–3500 m/s, Qs = 0.05·Vs, Qp = 2·Qs.
- **grid:** origin at NZ lat/lon, `az=329.29`, NZTM `proj=tmerc` parameters, 3-level refinement with h_coarse=400 m scaled down (e.g. 200/400/800 m) on a domain of about 30×30×15 km.
- **SRF-HDF5:** 2×3-point SRF with v2.0 VS/DEN, mixed rakes, SR1 only, written with the same layout as `source_modelling.srf.write_sw4_hdf5` (or produced by importing it, if available in CI).
- **stations.h5:** 6 stations, including one on the hilltop and one inside the sponge (to exercise SGDEPTH > 0).
- **Input file:** the full workflow command set (`fileio`, `grid`, `time`, `rupturehdf5`, `refinement`, `supergrid width=`, `attenuation`, `developer cfl=0.9 failonnan=1`, `prefilter`, `topography input=sfile order=3`, all 10 `imagehdf5`, `sfile`, `rechdf5`).
- **Runtime:** about 20–60 s on 4 ranks.

Tests (all share one run through a `FIXTURES_SETUP mininz_run`, except where noted):

| ID | Assertion | WF | Tier |
|---|---|---|---|
| WF-01 | Run exits 0; no NaN; no "no convergence" warning | ★★★ | pr |
| WF-02 | **rechdf5 output contract:** `DELTA` equals the SW4 dt × downsample; NPTS is equal across stations and equals ceil(t/dt); `NS`/`EW`/`UP` are datasets; STLA/STLO echo the input to 1e-9°; `SGWIDTH(GP)` present at the root; SGDEPTH/SGDEPTHGP datasets, > 0 only for the sponge station | ★★★ | pr |
| WF-03 | **Golden waveform regression** per precision. Relative L2 per trace below 1e-10 (double) or 1e-5 (float), *and* Kristeková envelope/phase misfit below 1e-3, so failures are diagnosable | ★★★ | pr |
| WF-04 | **imagehdf5 contract:** 10 files with expected names/modes and cycle/time stamps (cycle 0 for topo/grid/p/s/rho, every 0.5 s for mag/velmag/uz, t_end for hmax/vmax); dimensions match the top grid; `topo` equals the sfile elevation interpolated (≤ 1 cm); p/s/rho equal the sfile surface values at the nodes | ★★★ | pr |
| WF-05 | **Image/trace consistency:** at stations that coincide with grid nodes, `hmax` ≥ max horizontal amplitude of the rechdf5 trace and matches it within 1%; `vmax` likewise | ★★★ | pr |
| WF-06 | **Projection:** the x,y SW4 reports for each lat/lon station (parsed from stdout with `verbose=2`) match `pyproj` NZTM-without-false-origin + azimuth rotation within h/100 | ★★★ | pr |
| WF-07 | Antimeridian variant: a domain whose corners straddle 180° (as realisation 2012p578973 does). Same assertions as WF-02/06 | ★★★ | pr |
| WF-08 | 0.5Hz-shape (2 levels) and 0.25Hz-shape (no refinement) variants: WF-01/02/03 | ★★★ | pr |
| WF-09 | No-topography variant (NZCVM `decay` set, sfile min depth 0, no `topography` command, so SW4 runs Cartesian): WF-01/02/03 | ★★★ | pr |
| WF-10 | Larger mini-NZ at the full 1Hz ratio (100/200/400 m), 40 s, 8 ranks (oversubscribed) | ★★★ | nightly |

### 3.5 Physics-equivalence (metamorphic) tests for the workflow's input paths

Each test is two or more runs (fixtures) followed by a compare test.

| ID | Equivalence | Why it matters | Tier |
|---|---|---|---|
| EQ-01 | **1-point `rupturehdf5` SRF ≡ `source` moment tensor** with the matching strike/dip/rake/M0 and time function. Traces agree within discretisation round-off | Workflow point sources are 1-point SRFs; this checks the SRF → moment conversion, AREA·μ·slip, and the rake/strike conventions | pr |
| EQ-02 | **Slip-rate semantics:** with SRF slip *rate* input, rechdf5 "displacement" equals the velocity output (`variables=velocity`) of the same model driven by the integrated slip, i.e. the downstream assumption "output is velocity in m/s". Check against User Guide §11.2.2 when building the test | The whole downstream unit chain (×100, differentiate) rests on this | pr |
| EQ-03 | `rupturehdf5` ≡ ASCII `rupture` for the same SRF | Isolates HDF5 reader bugs | pr |
| EQ-04 | **sfile ≡ `block`** for a layered model aligned with the grid (flat, no topography) | Isolates sfile reader/interpolation bugs | pr |
| EQ-05 | Curvilinear grid with **flat** topography ≡ Cartesian grid (within round-off) | Checks that the curvilinear machinery is consistent | pr |
| EQ-06 | `rechdf5` ≡ `rec usgsformat=1` text output for the same stations | Isolates the HDF5 writer | pr |
| EQ-07 | `supergrid width=W` ≡ `supergrid gp=W/h_coarse` | The workflow uses `width`; SW4 measures `gp` on the coarsest grid | pr |
| EQ-08 | Reciprocity: point force at A recorded at B ≡ the reverse, in the mini-NZ heterogeneous, attenuating, curvilinear model (Eisner & Clayton 2001) | Solution-free check of operator self-adjointness and of source/receiver interpolation near interfaces | nightly |
| EQ-09 | Mirror symmetry: a symmetric model and source give mirrored traces | Catches index and orientation errors | pr |
| EQ-10 | **Azimuth rotation:** the same physical problem at az=0 vs az=329.29 (homogeneous half-space, point source), with NS/EW traces equal after rotation | The workflow relies on SW4's NSEW output with no rotation step of its own | pr |

### 3.6 Absorbing boundary (supergrid) quality

| ID | Test | WF | Tier |
|---|---|---|---|
| SG-01 | Reflection measurement: run domain A (workflow-like `width`), then domain B (3× larger, so the reflection arrives after t_end). Assert ‖A−B‖/‖B‖ over the station window is below the threshold implied by the logged absorbed period (`sw4.py:226`) | ★★★ | pr |
| SG-02 | Reflection vs width ladder: the error decreases monotonically with width | ★★ | nightly |
| SG-03 | Sponge inside the bottom refinement: the configuration the workflow enforces runs and passes SG-01 | ★★★ | pr |
| SG-04 | Existing `pytest/supergrid` guard cases (source in sponge, margin, rupture in sponge) as `WILL_FAIL`/`PASS_REGULAR_EXPRESSION` tests. Includes the 5-cell `margin_pts` that `workflow/sw4.py:20` mirrors | ★★★ | pr |

### 3.7 Prefilter and attenuation response

| ID | Test | WF | Tier |
|---|---|---|---|
| PF-01 | Spectral ratio of a prefilter run to an unfiltered run at a station equals \|H_butter(f)\|² for order 4 with 2 passes (zero-phase): within 1 dB below 0.8·fc2, half-power at fc2 within 2% | ★★★ | pr |
| PF-02 | Zero-phase check: cross-correlation peak lag = 0 | ★★★ | pr |
| AT-01 | Q(ω) fit unit test: the mechanisms fitted for `nmech=3 phasefreq=0.5 maxfreq=10` give Q within 5% of target over [maxfreq/100, maxfreq] (Emmerich & Korn 1987; Liu et al. 1976) | ★★★ | pr |
| AT-02 | 1D plane-wave decay: amplitude decay vs exp(−πfx/(Qc)) and the dispersion at phasefreq | ★★ | nightly |
| AT-03 | LOH.3 (`pytest/loh3/LOH.3-h100.in`, `phasefreq=2.5`) at station 10 vs the bundled PROSE reference `tools/LOH.3_prose_corrected`, post-processed as in `tools/loh3exact.m`: Kristeková EM/PM < 5% | ★★★ | nightly |

### 3.8 Parallel, restart and precision invariance

| ID | Test | WF | Tier |
|---|---|---|---|
| PI-01 | mini-NZ at 1, 2, 3 and 4 ranks: traces **bitwise** identical (or ≤ 1 ulp if reductions are involved; establish which) | ★★★ (HPC runs at many rank counts) | pr |
| PI-02 | OpenMP 1 vs 2 vs 4 threads: bitwise | ★★ | pr |
| PI-03 | Checkpoint at N/2 + restart ≡ straight run, bitwise, including memory variables and supergrid state. Ports `loh1-h100-mr-restart-hdf5` | ★ | pr |
| PX-01 | **Float vs double cross-check** on mini-NZ: per-trace relative L2 < 1e-4 and Kristeková EM/PM < 1e-3 | ★★★ (production is float) | pr (needs both builds; see §5) |
| PX-02 | Float twilight errors within 10× of double; replaces the blanket 5e-1. Triages the 3.6e-1 cases | ★★ | pr |

### 3.9 Benchmarks and legacy coverage

| ID | Test | Tier |
|---|---|---|
| BM-01 | LOH.1 (`examples/scec/LOH.1-h50.in` run at h=100) at station 10 vs the bundled PROSE reference `tools/LOH.1_prose3`, post-processed as in `tools/loh1exact.m`: Kristeková misfit < 5%. Also asserts the vertical sign convention (the two MATLAB scripts disagree in their comments) | nightly |
| BM-05 | *Optional:* NZ-like layered 1D model (workflow Q rules Qs=0.05·Vs, Qp=2·Qs, `phasefreq=0.5`), several stations, point moment tensor. Compared against wavenumber-integration references generated once offline (CPS `hspec96` or Axitra) and committed with a provenance header | nightly |
| BM-02 | Existing hdf5 cases (`loh1-h100-mr-hdf5`, `-sfile`, `-restart`) ported to ctest with `verify_hdf5.py` | pr |
| BM-03 | Existing geodynbc cases | nightly |
| BM-04 | Smoke tests (runs, finite output, golden L2) for features the workflow does not use: `rfile`, `pfile`, `ifile`, `image`, `volimage`, `ssioutput`, `sfileoutput`, anisotropy | nightly |

### 3.10 Negative and input-validation tests

These use `WILL_FAIL` or `PASS_REGULAR_EXPRESSION` and run in a few seconds, all at tier pr.

- `cfl > 1.5` is rejected.
- A source inside the sponge aborts.
- The sfile does not cover the grid → clear error.
- Domain bottom below the sfile zmax → error.
- A refinement level with fewer than the minimum gridpoints → error (protects the `nz_min=12` assumption).
- A NaN injected with `failonnan=1` gives a non-zero exit (the workflow relies on this to fail fast on HPC).
- Malformed SRF-HDF5 (missing `SR1` or `PLANE`) → clear error, not a segfault.
- A station outside the domain → warned/skipped, and NPTS stays consistent.

### 3.11 HDF5 I/O robustness

| ID | Test | Tier |
|---|---|---|
| IO-01 | z-plane image at more ranks than writers (lock errno-11 race). `WILL_FAIL` until fixed, then a regression test | pr |
| IO-02 | `fileio nwriters` workaround path exercised | pr |
| IO-03 | rechdf5 with `writeEvery` smaller than steps and with the default 1000 > steps: the final file is complete in both cases | pr |

### 3.12 C++ unit tests (new, doctest, single-header)

These are fast and run without MPI where possible, labelled `unit`, tier pr:

- Butterworth prefilter coefficients and frequency response.
- Attenuation mechanism fit (the AT-01 core, without a solve).
- Projection helpers (grid lat/lon ↔ x/y, azimuth).
- sfile header parsing and interface interpolation.
- SRF-HDF5 point parsing (AREA, rake, unit conversions).
- `Sarray` indexing and halo-exchange packing.
- Supergrid damping profile, `gp` ↔ `width` conversion.

### 3.13 Cross-repo contract test (workflow e2e with a real SW4)

| ID | Test | Tier |
|---|---|---|
| XR-01 | Check out `ucgmsim/workflow` (pinned ref, bumped by PR) and run its e2e (`WORKFLOW_E2E=1 SW4=<built sw4> SW4_LAUNCHER="mpirun -n 4"`) instead of `fake_sw4.py`. Then run `lf-to-xarray --format sw4` on the output. This checks the real generator → real SW4 → real reader chain end to end | nightly (and on PRs labelled `workflow`) |

This is the only test that catches a change in what `create-sw4-input` emits that SW4 no longer accepts. Because `SW4_COMMAND_SCHEMA` allows arbitrary commands, the risk is real.

### Count

| Section | Test specs | PR | Nightly-only | Ports of existing |
|---|---|---|---|---|
| 3.1 MMS | 12 | 10 | 2 | ~4 |
| 3.2 Analytic | 5 | 3 | 2 | 2 |
| 3.3 Energy | 7 | 5 | 2 | 1 (5 inputs) |
| 3.4 mini-NZ | 10 | 9 | 1 | 0 |
| 3.5 Equivalence | 10 | 9 | 1 | 0 |
| 3.6 Supergrid | 4 | 3 | 1 | 1 |
| 3.7 Prefilter/atten | 5 | 3 | 2 | 0 |
| 3.8 Invariance/precision | 5 | 5 | 0 | 1 |
| 3.9 Benchmarks | 5 | 1 | 4 | 2 |
| 3.10 Negative | 8 | 8 | 0 | 0 |
| 3.11 HDF5 I/O | 3 | 3 | 0 | 0 |
| 3.12 Unit | 7 | 7 | 0 | 0 |
| 3.13 Cross-repo | 1 | 0 | 1 | 0 |
| **Total** | **82** | **66** | **16** | **~12** |

- **About 70 of the 82 are new.** The other ~12 replace existing cases with stronger assertions.
- **Each spec expands into several ctest entries:** one run per resolution level or comparison variant, plus a check. The PR tier is therefore about **160 ctest entries per precision** (≈ 320 per PR across double and float), against 48 today. Nightly adds about 60 more.
- Every feature in §1.2 has at least one test that would fail if it broke.

---

## 4. Harness design

### 4.1 Layout

```
tests/
  CMakeLists.txt            # add_subdirectory from root; registers everything below
  sw4test/                  # python package, stdlib + numpy + h5py (+ scipy for PF/Kristeková)
    run.py                  # ctest entry point: `python -m sw4test.run <case.toml> --stage run|check`
    checks.py               # norm compare, rate fit, energy monotone, trace compare, misfits
    rechdf5.py, images.py   # readers for SW4 HDF5 outputs (contract checks)
    kristekova.py           # TF envelope/phase misfit (port of the ObsPy tf_misfit algorithm)
    make_mininz.py          # deterministic sfile/SRF/stations generator
  cases/
    mms/mms-06.toml + mms-06.in.j2 ...
    workflow/mininz.toml ...
  golden/
    double/<case>/...       # small: norms, decimated traces (.npz); < 20 MB total
    float/<case>/...
  unit/                     # doctest C++ tests
```

### 4.2 Case spec (TOML, one per case)

```toml
name = "mms-06"
labels = ["pr", "twilight", "curvi", "mr"]
ranks = 4
threads = 1
timeout = 120
input = "mms-06.in.j2"            # templated: {{h}}, {{nx}}, ...
levels = [{h=0.04}, {h=0.02}]     # nightly adds {h=0.01}
[[check]]
kind = "golden_norms"             # TwilightErr.txt last block, all columns
rtol = { double = 1e-8, float = 1e-3 }
[[check]]
kind = "rate"
expect = 4.0
band = 0.4
floor = { float = 1e-5 }          # skip the rate check once errors near the round-off floor
[[check]]
kind = "stdout_absent"
pattern = "no convergence"
```

`tests/CMakeLists.txt` globs `cases/**/*.toml` (via a tiny CMake-side parser or `execute_process(python -m sw4test.list)` at configure time). For each case it emits:

- `run:<case>:<level>` tests with `FIXTURES_SETUP <case>-<level>` and `PROCESSORS ranks*threads`;
- one `check:<case>` test with `FIXTURES_REQUIRED` on all its levels. This fixes the current `DEPENDS` bug structurally;
- labels: `pr` cases at levels ≤ PR depth; `nightly` adds the extra levels and long cases. Selection uses `ctest -L pr` or `-L nightly`. `TESTING_LEVEL` is kept for compatibility and maps onto the labels.

### 4.3 Goldens and tolerances

- Goldens are generated **only** by `python -m sw4test.bless <case> --precision double|float`. It writes a header with commit, compiler and flags, and the PR diff must show the change. Blessing is reviewed like code.
- Starting tolerances:
  - double: 1e-8 relative on norms, 1e-10 on traces;
  - float: 1e-4 on norms, 1e-5 on traces;
  - convergence and energy assertions are precision-independent physics checks.
- These replace the blanket 5e-1.
- **Why bitwise is out:** goldens cannot be bitwise because compiler and FMA differences show up across the gcc/clang/intel/aocc matrix. Bitwise comparison is used only *within a single build* (PI-01..03). `SW4_STRICT_FP=ON` is available if cross-compiler bitwise checks ever become desirable.

### 4.4 Retire or fix the old harness

- Fix `check_results.py`: compare all columns by name, use `abs` in the denominator, fix `printf`. Then fold it into `checks.py` and keep a shim.
- Make `test_sw4.py` exit non-zero on failure. Add a deprecation note pointing to ctest; `linux.yml` is replaced by the new workflow.
- Delete the old CMake test block once its cases are ported (same inputs, now with rate checks).
- Add the generated `pytest/<dir>/` output directories to `.gitignore`.

---

## 5. GitHub Actions: `.github/workflows/tests.yml`

```
on: pull_request, push to master/single_precision, schedule (nightly 13:00 UTC), workflow_dispatch
concurrency: cancel in-progress per ref

jobs:
  build (matrix precision: [double, float]):
    ubuntu-24.04 (4 vCPU / 16 GB)
    apt: gfortran, libopenmpi-dev, libhdf5-openmpi-dev, libproj-dev, libfftw3-mpi-dev, liblapack-dev, ninja, ccache
    ccache via actions/cache keyed on compiler + precision + hashFiles(src/**)
    reuse .github/ci/build.sh gcc <precision>  (+ -DSW4_TESTS=ON)
    upload build/bin/sw4 + build/tests as an artifact

  test (needs build; matrix precision; tier = pr on PRs, nightly on schedule):
    pip install numpy h5py scipy pyproj  (cached)
    env: OMPI_MCA_rmaps_base_oversubscribe=1, PRTE_MCA_rmaps_default_mapping_policy=:oversubscribe,
         OMP_NUM_THREADS=1, OMPI_ALLOW_RUN_AS_ROOT* if containerised
    ctest --test-dir build -L ${tier} -j4 --resource-spec-file .github/ci/resources.json
          --timeout 600 --output-on-failure --output-junit junit.xml
    on failure: upload tests/out/** (stdout, .h5, diff plots) as an artifact
    publish JUnit summary to the PR checks (e.g. mikepenz/action-junit-report)

  cross-precision (needs both test jobs; pr tier):
    download the float+double outputs of mini-NZ → run PX-01

  workflow-e2e (nightly, or PRs labelled `workflow`):
    checkout ucgmsim/workflow@<pinned>, uv sync, run e2e with the real sw4 (XR-01)
```

- **Budget:** the PR tier should finish in **≤ 25 min wall per precision** on a 4-vCPU runner. The level-0 runs today take about 6–7 min on 4×2 cores locally, so the expanded catalogue needs careful grid sizing; `tools/ctest_timings.py` reports per-test times so slow cases can be moved to nightly.
- The CI test environment is fixed to gcc. The existing `compilers.yml` still builds every compiler. Optionally, a nightly job runs the `pr` tier on clang and intel, with per-compiler golden tolerances (not separate goldens).
- Make `tests / double` and `tests / float` **required status checks** on master and single_precision.

---

## 6. Phased delivery

| Phase | Content | Exit criterion |
|---|---|---|
| **0. Make failures visible** (≈ 1 day) | Fix the ctest `DEPENDS` bug, the `check_results.py` bugs and the `test_sw4.py` exit code. Add `tests.yml` running the existing ctest suite at double and float. Add the gitignore for outputs. | The CI is red when an existing check fails. Confirm by deliberately perturbing one golden. |
| **1. Harness + rates** | `sw4test` package, TOML cases, fixtures and labels. Port all existing pytest cases (curvimeshrefine, energy, hdf5, supergrid) into ctest. MMS-01..09, AN-01..03, rate checks. Per-precision goldens and real float tolerances; triage the two 3.6e-1 cases. | ≈ 60 PR tests green on both precisions. |
| **2. Workflow contract** | `make_mininz.py`; WF-01..09, EQ-01..07 + 09/10, SG-01/03/04, PF-01/02, PI-01..03, PX-01/02, negative tests. | Every row in §1.2 is covered by a failing-if-broken test. |
| **3. Stability + I/O** | Double energy accumulation. EN-01..07 (incl. the attenuation `testenergy` patch). IO-01..03. | Known defects from §1.3 encoded as `WILL_FAIL` tests that flip when fixed. |
| **4. Validation + units** | doctest unit tests, AT-01..03, BM-01 LOH.1/LOH.3 against the bundled PROSE references (port `tools/loh{1,3}exact.m` + `ReadUHS.m` to `sw4test/loh.py`), optional BM-05, EQ-08, XR-01 nightly. | Nightly green; LOH misfits recorded. |

Each phase is a separate PR, and each PR in phases 1–4 runs under the CI from phase 0.

---

## 7. Open questions / decisions

1. **Bitwise across rank counts:** does SW4 already achieve this for traces? Run PI-01 once to find out. If not, decide between a tolerance and fixing the reductions.
2. **LOH references: resolved.** The repo already ships the original PEER frequency-wavenumber references (Apsel & Luco PROSE): `tools/LOH.1_prose3` and `tools/LOH.3_prose_corrected`. Each is 2048 samples at dt=0.008 s for station 10, with columns t, vertical, radial, transverse. `tools/loh{1,3}exact.m` convolves them to SW4's Gaussian source. No new solver is needed for BM-01 or AT-03.
   - *Open:* whether to add BM-05, a layered-model reference beyond the LOH geometry. If yes, use CPS (strongest causal-dispersion/Q handling, explicit reference frequency) or Axitra (GitHub, Python wrapper, forces and moment tensors), run offline, and commit the results. pyfk or fk only if their Q and reference-frequency conventions check out. Running another 3D code (SPECFEM3D, SeisSol) is not worth it except for a one-off topography benchmark.
3. **EQ-02 semantics:** confirm against User Guide §11.2.2 exactly what rechdf5 "displacement" means for `rupturehdf5`. If the downstream assumption is wrong, this test is the one that will reveal it.
4. **Is `ucgmsim/workflow` reachable from SW4's Actions** (public, or needs a token) for XR-01?
5. **The `testenergy`-with-attenuation patch** (`testenergy-master.patch`) needs upstreaming before EN-04/05 can run in CI.
6. **Precision on the HPC target:** the plan treats float as production. If the workflow switches to double for some runs, the PX-01 tolerance becomes a release criterion.

---

## References

- Salari & Knupp (2000), *Code Verification by the Method of Manufactured Solutions*, SAND2000-1444.
- Roache (2002), *J. Fluids Eng.* 124; Oberkampf & Roy (2010), *Verification and Validation in Scientific Computing*, CUP.
- Nilsson, Petersson, Sjögreen & Kreiss (2007), *SIAM J. Numer. Anal.* 45(5).
- Sjögreen & Petersson (2012), *J. Sci. Comput.* 52(1).
- Appelö & Petersson (2009), *CiCP* 5(1).
- Petersson & Sjögreen (2010), *CiCP* 8(5) — grid refinement and singular sources.
- Petersson & Sjögreen (2012), *CiCP* 12(1) — attenuation.
- Petersson & Sjögreen (2014), *CiCP* 16(4) — supergrid.
- Zhang, Wang & Petersson (2021), *SIAM J. Sci. Comput.* 43(2) — curvilinear + refinement.
- Day et al. (2001, 2003), PEER Lifelines 1A01/1A02 — LOH.1/2/3.
- Bielak et al. (2010), *GJI* 180 — ShakeOut verification.
- Chaljub et al. (2010) *BSSA* 100(4); Chaljub et al. (2015) *GJI* 201(1) — E2VP.
- Kristeková et al. (2006) *BSSA* 96(5); Kristeková, Kristek & Moczo (2009) *GJI* 178(2) — misfit criteria.
- Anderson (2004), 13th WCEE Paper 243 — goodness of fit.
- Emmerich & Korn (1987) *Geophysics* 52(9); Liu, Anderson & Kanamori (1976) *GJRAS* 47; Day & Bradley (2001) *BSSA* 91(3).
- Eisner & Clayton (2001), *BSSA* 91(3) — reciprocity.
- Apsel & Luco (1983), *BSSA* 73(4) — PROSE layered half-space Green's functions (the LOH reference solutions in `tools/`).
- Zhu & Rivera (2002) *GJI* — fk; Herrmann (2013) *SRL* — CPS.
- Kanewala & Bieman (2014), *IST* 56(10); Hook & Kelly (2009), SECSE; Dubey et al. (2013), FLASH; xSDK Community Package Policies.
