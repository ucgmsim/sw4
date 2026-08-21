# ab-bench — A/B performance harness for SW4

Builds two commits from clean git worktrees, runs a fixed case set through both,
and reports the per-phase timing breakdown side by side.

```bash
# defaults: base = 23a3410, head = HEAD, double precision, size M, 5 reps
./performance/ab-bench/ab-bench.sh --target mn5-gpp
```

Load your compiler and MPI modules first — the script deliberately does not
touch them, so that `platform.txt` records what you actually built with.

## Why this exists

Vectorisation improvements in this code have been easy to verify statically and
hard to predict dynamically. Over the 2026-08-21 optimisation series, every
prediction of *codegen* held (zmm counts, packed/scalar ratios, vectoriser
diagnostics) and several predictions of *time* were wrong by an order of
magnitude:

| prediction | basis | outcome |
|---|---|---|
| removing 35% of the RHS kernel's compulsory traffic | traffic counting | 1.4% |
| `__restrict__` on the BC kernels | aliasing analysis | 0% |
| reciprocal `strx` array, ~2% | latency/throughput ratio | 0.3%, reverted |
| halo-only zeroing, 15-25% on a saturated socket | roofline | 0% locally, unresolved |

The pattern is consistent: instruction-level wins do not convert to time in a
way that can be reasoned about from a machine with a different vector width,
core count and bandwidth-per-core. Hence: measure, on the target, with this.

## Reading the output

`summary.txt` gives min-of-reps per phase and a speedup column. The phases come
from SW4's own `developer reporttiming=1` breakdown:

| column | what dominates it |
|---|---|
| `divstress` | `rhs4th3fort*_ci`, `curvilinear4sg*_ci` — the RHS stencils |
| `forcing` | point-source / analytic forcing. Twilight cases inflate this |
| `bc` | boundary conditions, mostly data movement |
| `sg` | supergrid damping |
| `comm` | halo exchange. Only meaningful with `--ranks` > 1 across nodes |
| `mr` | mesh-refinement interfaces, incl. the *windowed* 3D RHS |
| `updates` | predictor/corrector/dpdmt |

Watch `divstress` for stencil work and `comm` for the halo exchange. `total` is
diluted by `forcing`, which is unrepresentative in twilight cases.

Correctness is checked on **both** Linf and L2. Linf is a max and misses ~1e-8
relative changes that the L2 sum catches — this actually happened during the
series. If the harness reports `NON-DETERMINISTIC across reps`, that is
expected for anything touching an OpenMP reduction (SW4's energy diagnostic is
not reproducible above one thread) and means timings are still valid but the
correctness check is not.

Note `--strict-fp` builds both sides with `SW4_STRICT_FP=ON`, which disables FMA
contraction. Two changes in the series (the peel-param unroll and the
curvilinear fission) are bit-identical only under that flag, because extra
unrolling exposes extra FMA-contraction opportunities. Use it to separate "did
the numbers change because of my change" from "did they change because the
compiler re-associated a multiply-add".

## Asymmetry you should know about

The default base commit `23a3410` predates `cmake/SW4Optimization.cmake`. On
that side the harness falls back to the old `SW4_ARCH_FLAGS` interface and
passes an equivalent `-march` by hand, so the comparison is not accidentally
"tuned flags vs no flags". It is still not perfectly symmetric — the base side
lacks `-fno-math-errno`, `-fipa-pta` and the alignment flags, which are part of
what is being tested. `platform.txt` records the exact flags for both sides.

To isolate a single change instead, pass explicit commits:

```bash
# just the ivdep -> omp simd fix
./ab-bench.sh --base 23a3410 --head 5374160
# just the curvilinear closure unroll
./ab-bench.sh --base 86fd071 --head 5c1140d
# just the Cartesian closure unroll
./ab-bench.sh --base 5c0c2e5 --head bcafa0a --cases cart-mr
```

## If a case fails

A run that produces no timing table is recorded in `failures.txt`, marked with
`!` instead of `.` in the progress line, excluded from the summary, and warned
about at the end — the run continues. The usual cause is a geometry SW4 rejects;
it requires **at least 12 z-points (excluding ghosts) per grid**, so a
refinement interface placed too close to the curvilinear bottom aborts with
`Precondition violated: The number of grid points ... must be >= 12`. If you add
a case or a size, check it at the *coarsest* size you intend to use, since `h`
grows as `--size` shrinks.

## Submitting on a cluster without giving anyone shell access

`slurm/hpc3-ab.sl` and `slurm/hpc3-comm.sl` are self-contained: copy the repo
across, `sbatch`, collect the output file. Nothing needs an interactive session,
which matters where login requires a portal/MFA round trip.

```bash
# Copy INCLUDING .git -- the harness uses `git worktree` to check out both
# commits, so a source-only copy cannot work.
rsync -az --exclude 'build*' --exclude 'ab-*' sw4/ hpc3:~/sw4-ab/

# State the layout explicitly. --exclusive is NOT in the script headers,
# because whether you want it depends on what you are measuring:
#
#   Shared node -- schedules in minutes rather than a day, since a whole free
#   node is rare. Good enough for the vectorisation and thread-scaling
#   questions. Confirm core counts first: sinfo -p <part> -o '%n %c %m'
sbatch -p genoa --ntasks=2 --cpus-per-task=32 --mem-per-cpu=2G \
       --hint=nomultithread performance/ab-bench/slurm/hpc3-ab.sl
sbatch -p milan --ntasks=2 --cpus-per-task=32 --mem-per-cpu=2G \
       --hint=nomultithread performance/ab-bench/slurm/hpc3-ab.sl

#   Whole node -- needed for the bandwidth-saturation question, and only that.
sbatch -p genoa --exclusive --mem=0 --ntasks=2 --cpus-per-task=84 \
       --hint=nomultithread performance/ab-bench/slurm/hpc3-ab.sl

#   Multi-node halo exchange. Defaults to isolating 5374160 -> 668bd37, where
#   nothing but parallelStuff.C differs.
sbatch -p genoa --nodes=4 --ntasks-per-node=2 --cpus-per-task=32 \
       --mem-per-cpu=2G --hint=nomultithread \
       performance/ab-bench/slurm/hpc3-comm.sl
```

The summary is echoed into the job's `.out` file, so the whole result arrives
without needing to fetch anything. Knobs are environment variables:

```bash
SIZE=XL REPS=7 sbatch -p genoa performance/ab-bench/slurm/hpc3-ab.sl
PRECISION=both sbatch -p genoa performance/ab-bench/slurm/hpc3-ab.sl
STRICT=1 sbatch -p genoa performance/ab-bench/slurm/hpc3-ab.sl   # bit-reproducible
BASE=23a3410 HEAD_REF=HEAD sbatch -p genoa performance/ab-bench/slurm/hpc3-comm.sl
```

Module names are the one thing not knowable from outside the site. The scripts
probe a list of plausible sets (`foss CMake`, `GCC OpenMPI OpenBLAS CMake`, …)
and verify by looking for `cmake`/`mpicxx`/`gfortran`/`git`/`python3` on PATH
rather than trusting `module load` exit codes. If none work they print the
relevant slice of `module avail` and tell you to re-submit with:

```bash
MODULES='GCC/13.2.0 OpenMPI/4.1.6 CMake/3.27' sbatch -p genoa ...
```

Uncomment the `--account` line in the header if your site requires one.

## Platform recipes

### Target clusters

```bash
# MareNostrum 5 GPP — 2x Xeon Platinum 8480+, 112 cores
sbatch --qos=gp_debug -N1 -n1 -c112 -t 1:00:00 --wrap \
  "module load intel impi cmake; \
   ./performance/ab-bench/ab-bench.sh --target mn5-gpp --size L --threads 112"

# REANNZ HPC3 — heterogeneous, so pin the partition
sbatch -p genoa -N1 -n1 -c168 -t 1:00:00 --wrap \
  "./performance/ab-bench/ab-bench.sh --target hpc3-genoa --size L --threads 168"
sbatch -p milan -N1 -n1 -c126 -t 1:00:00 --wrap \
  "./performance/ab-bench/ab-bench.sh --target hpc3-milan --size L --threads 126"

# TACC Frontera / Stampede3 / Vista
idev -N 1 -n 1 -t 1:00:00
./performance/ab-bench/ab-bench.sh --target frontera --size L --threads 56
./performance/ab-bench/ab-bench.sh --target stampede3-spr --size L --threads 112
./performance/ab-bench/ab-bench.sh --target vista --size L --threads 144
```

For the halo-exchange work, `comm` only becomes meaningful across nodes:

```bash
sbatch -N 4 -n 4 --wrap "./performance/ab-bench/ab-bench.sh --ranks 4 --size L"
```

### Cheap cloud, when you do not have an allocation yet

Spot instances give the exact microarchitectures for roughly the price of a
coffee. Use a full-socket size if you want the bandwidth-saturation and
thread-scaling answers; a 2-4 vCPU instance will not show either.

| target | AWS | GCP |
|---|---|---|
| MareNostrum 5, Stampede3 SPR | `c7i.24xlarge` | `c3-highcpu` |
| Cascade, HPC3 genoa | `c7a.24xlarge` | `c3d-highcpu` |
| HPC3 milan | `c6a.24xlarge` | `n2d` |
| Frontera CLX | `c6i` / `c5n` | `c2` |
| Vista Grace | `c7g` / `c8g.16xlarge` | `c4a` (Axion) |

```bash
sudo apt-get install -y build-essential gfortran cmake libopenmpi-dev \
     openmpi-bin libopenblas-dev liblapack-dev
git clone <repo> sw4 && cd sw4
./performance/ab-bench/ab-bench.sh --target cascade --size L
```

### Free options

- **Chameleon Cloud / CloudLab / Jetstream2** — NSF-funded, free for academics,
  **bare metal** (no virtualisation noise, which matters for bandwidth
  measurements). Chameleon has Sapphire Rapids and Milan nodes.
- **GitHub Actions** — free on public repos, and runners are commonly Ice Lake
  (AVX-512) or EPYC Milan. Only 4 vCPU, so no thread-scaling answer, but enough
  to prove the vectorisation pays on a real AVX-512 part and to make it a
  permanent regression check. Use `--size S --reps 3`.
- **Oracle Cloud always-free** gives 4 Ampere A1 cores, but that is Neoverse
  **N1** — no SVE, so not a Grace proxy.
- **Google Colab** is the weakest option: 2 vCPUs, no control over which Xeon
  you land on, no MPI scaling. It can confirm AVX-512 code runs and little else.

`SW4CK` (github.com/LLNL/SW4CK) is a mini-app of five SW4 stencil kernels
covering roughly half of solve time. On a small instance it isolates the kernels
without needing a full simulation, which may be a better fit than this harness.

## Reading a shared-node result

A shared-node run gives each core several times the memory bandwidth it would
get on a full node. That matters because the two kernels sit on opposite sides
of the roofline:

* **Curvilinear** (`curvi`, `curvi-mr`) is compute-bound at every occupancy
  (AI 32.3 flop/B against a machine balance of 14-25). Its speedups transfer
  directly to a full node.
* **Cartesian** (`cart`, `cart-mr`) is memory-bound on x86 by 1.15-1.40x at full
  occupancy and comfortably compute-bound at one-third occupancy. Its
  shared-node speedups are therefore an **upper bound** -- expect them to
  compress on a full node.

Neighbouring jobs also contend for bandwidth. The A/B interleaving and
min-of-reps mitigate that but cannot remove it, so treat a shared-node number as
a range rather than a point. The job output prints which mode it ran in.

## What to expect, and what would falsify it

The changes fall into four classes that scale very differently, which is why a
single number is not meaningful:

1. **Bandwidth-capped, no gain.** The Cartesian *interior* loop was already
   vectorised at 512-bit. AI 12.34 flop/B against an effective machine balance
   of 16.5 on Sapphire Rapids: memory-bound by 1.34x, so better code cannot
   help it.
2. **Scales with vector width.** The closure loops went scalar -> vectorised.
   4 float lanes on the development box, 16 on SPR/Genoa. Measured 1.87-2.52x
   locally on thin-nk Cartesian geometries; should be larger on target.
3. **Scales with thread count.** `collapse(2)` on six-iteration loops was worth
   2-4% at 4 threads. At 56-112 cores it is addressing 57-89% idle threads.
4. **Width- and thread-independent.** The curvilinear fission (1.175x) works via
   dependency chains and addressing and should transfer as-is.

Class 2 and 3 gains are partly absorbed at full socket occupancy, where
bandwidth per core drops to ~5.5 GB/s on SPR. That is exactly the regime a
4-core development box cannot reproduce, and exactly what this harness is for.
