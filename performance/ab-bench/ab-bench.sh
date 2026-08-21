#!/usr/bin/env bash
#
# A/B performance harness for SW4.
#
# Builds two commits from clean worktrees, runs a fixed case set through both,
# and reports the per-phase timing breakdown side by side. Designed to be copied
# onto an unfamiliar machine and run without setup: it records what it built
# with, interleaves the A and B runs so thermal drift cannot bias one side, and
# checks that the two binaries still agree numerically before believing any
# timing.
#
# Usage:  ./ab-bench.sh [--base COMMIT] [--head COMMIT] [options]
# Help:   ./ab-bench.sh --help
#
set -euo pipefail

# ---------------------------------------------------------------------------
# defaults
# ---------------------------------------------------------------------------
BASE="23a3410"          # last commit before the 2026-08-21 optimisation series
HEAD_REF="HEAD"
PRECISIONS="double"
TARGET="auto"
STRICT_FP="OFF"
SIZE="M"
REPS=5
RANKS=1
THREADS=""
CASES="cart,cart-mr,curvi,curvi-mr"
STEPS=200
NX_OVERRIDE=""
T_OVERRIDE=""
OUTDIR=""
LAUNCHER=""
DO_BUILD=1
DO_RUN=1
JOBS=""
EXTRA_CMAKE=""

usage() {
  sed -n '3,13p' "$0" | sed 's/^# \{0,1\}//'
  cat <<EOF

Options
  --base COMMIT      baseline commit           (default: $BASE)
  --head COMMIT      commit under test         (default: $HEAD_REF)
  --precision LIST   single|double|both        (default: $PRECISIONS)
  --target NAME      SW4_TARGET preset, e.g. mn5-gpp, cascade, hpc3-genoa,
                     hpc3-milan, frontera, stampede3, stampede3-spr, vista,
                     native, generic, auto                 (default: $TARGET)
  --strict-fp        build both sides with SW4_STRICT_FP=ON (bit-reproducible)
  --size S|M|L|XL    points PER RANK: 200k/800k/3.2M/12.8M (default: $SIZE)
  --steps N          timesteps per run        (default: $STEPS)
  --nx N             fixed grid instead of per-rank scaling
  --time T           simulated time, with --nx
  --reps N           timed repetitions, interleaved         (default: $REPS)
  --ranks N          MPI ranks                 (default: $RANKS)
  --threads N        OMP_NUM_THREADS           (default: cores/ranks)
  --cases LIST       comma list from: cart,cart-mr,curvi,curvi-mr,prod,prod-curvi,prod-mr,aniso
                     (default: $CASES)
  --launcher CMD     override the MPI launcher, e.g. "srun" or
                     "mpirun --allow-run-as-root --bind-to core"
  --jobs N           parallel compile jobs     (default: nproc)
  --cmake "ARGS"     extra args passed to both cmake configures
  --outdir DIR       results directory         (default: ab-results-<stamp>)
  --build-only       build, do not run
  --skip-build       reuse an existing --outdir's builds
  -h, --help         this text

Notes
  * Run from anywhere inside the SW4 git repository.
  * Load your compiler/MPI modules first; this script does not touch them.
  * The base commit predates cmake/SW4Optimization.cmake, so on that side the
    harness falls back to the old SW4_ARCH_FLAGS interface and derives
    equivalent -march flags itself. That asymmetry is the point: it is what the
    two commits actually build as, and it is recorded in platform.txt.
  * Correctness is checked on BOTH Linf and L2. Linf alone is a max and misses
    ~1e-8 relative changes that show up in the L2 sum.
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --base) BASE="$2"; shift 2;;
    --head) HEAD_REF="$2"; shift 2;;
    --precision) PRECISIONS="$2"; shift 2;;
    --target) TARGET="$2"; shift 2;;
    --strict-fp) STRICT_FP="ON"; shift;;
    --size) SIZE="$2"; shift 2;;
    --reps) REPS="$2"; shift 2;;
    --ranks) RANKS="$2"; shift 2;;
    --threads) THREADS="$2"; shift 2;;
    --cases) CASES="$2"; shift 2;;
    --steps) STEPS="$2"; shift 2;;
    --nx) NX_OVERRIDE="$2"; shift 2;;
    --time) T_OVERRIDE="$2"; shift 2;;
    --launcher) LAUNCHER="$2"; shift 2;;
    --jobs) JOBS="$2"; shift 2;;
    --cmake) EXTRA_CMAKE="$2"; shift 2;;
    --outdir) OUTDIR="$2"; shift 2;;
    --build-only) DO_RUN=0; shift;;
    --skip-build) DO_BUILD=0; shift;;
    -h|--help) usage; exit 0;;
    *) echo "unknown option: $1" >&2; usage >&2; exit 2;;
  esac
done
[[ "$PRECISIONS" == "both" ]] && PRECISIONS="single double"

REPO="$(git rev-parse --show-toplevel)"
cd "$REPO"
NPROC="$( (command -v nproc >/dev/null && nproc) || sysctl -n hw.ncpu 2>/dev/null || echo 4)"
: "${JOBS:=$NPROC}"
: "${THREADS:=$(( NPROC / RANKS > 0 ? NPROC / RANKS : 1 ))}"
[[ -z "$OUTDIR" ]] && OUTDIR="$REPO/ab-results-$(git rev-parse --short "$HEAD_REF")-$$"
mkdir -p "$OUTDIR"
OUTDIR="$(cd "$OUTDIR" && pwd)"

BASE_SHA="$(git rev-parse --short "$BASE")"
HEAD_SHA="$(git rev-parse --short "$HEAD_REF")"

# ---------------------------------------------------------------------------
# launcher
# ---------------------------------------------------------------------------
if [[ -z "$LAUNCHER" ]]; then
  if [[ -n "${SLURM_JOB_ID:-}" ]] && command -v srun >/dev/null; then
    LAUNCHER="srun -n $RANKS"
  elif command -v mpirun >/dev/null; then
    LAUNCHER="mpirun -np $RANKS"
    # OpenMPI refuses to run as root without this; harmless otherwise.
    [[ "$(id -u)" == "0" ]] && LAUNCHER="$LAUNCHER --allow-run-as-root"
    mpirun --version 2>&1 | grep -qi 'open mpi' && LAUNCHER="$LAUNCHER --oversubscribe"
  else
    echo "no mpirun or srun found; pass --launcher" >&2; exit 1
  fi
else
  case "$LAUNCHER" in *" -n "*|*" -np "*) :;; *) LAUNCHER="$LAUNCHER -n $RANKS";; esac
fi

# ---------------------------------------------------------------------------
# case definitions.  Sizes chosen so M finishes in ~1 min/side on a laptop core
# and XL saturates a server socket.  Every case turns on reporttiming, and the
# twilight forcing gives an analytic error norm to check correctness against.
# ---------------------------------------------------------------------------
# Sizes express points PER RANK, not a fixed global grid. A fixed grid divided
# over a growing rank count collapses the per-rank work: an earlier 16/32/64-rank
# sweep at a fixed nx=161 produced total runtimes of 2.2s, 0.73s and 0.48s, where
# startup and jitter dominate and nothing is measurable. The untouched `forcing`
# phase varied +/-8% across those rounds, which is the noise floor they were
# being read against.
case "$SIZE" in
  S)  PPR=200000    ;;
  M)  PPR=800000    ;;
  L)  PPR=3200000   ;;
  XL) PPR=12800000  ;;
  *)  echo "bad --size: $SIZE" >&2; exit 2;;
esac
if [[ -n "$NX_OVERRIDE" ]]; then
  NX="$NX_OVERRIDE"; T="${T_OVERRIDE:-0.15}"
else
  # dt/h ~= 0.349 for the twilight material (measured: nx=101 t=0.15 -> 43
  # steps), so holding the step count fixed means t scales as 1/nx.
  NX=$(python3 -c "print(round((${PPR}*${RANKS})**(1/3))+1)")
  T=$(python3 -c "print(f'{${STEPS}*0.349/(${NX}-1):.6g}')")
fi
EST_PT=$(( PPR * STEPS ))
echo "  derived: nx=$NX t=$T  (~${PPR} pts/rank x ${STEPS} steps"
echo "           ~$(python3 -c "print(f'{$EST_PT*8.75e-9:.1f}')")s/run/side estimated at 4 threads/rank)"

emit_case() {  # $1=name  -> writes $OUTDIR/cases/$1.in
  local n="$1" f="$OUTDIR/cases/$1.in"
  mkdir -p "$OUTDIR/cases"
  {
    case "$n" in
      cart)
        # Cartesian, twilight. NOTE: `supergrid gp=20` only TUNES the taper --
        # it does not enable supergrid. m_use_supergrid flips only when a
        # boundary condition is bSuperGrid (setupRun.C:2066), and twilight
        # cases leave the BCs at their all-Dirichlet default. So this case
        # exercises rhs4th3fort_ci, NOT rhs4th3fortsgstr_ci, and never calls
        # addSuperGridDamping. Kept as-is because twilight's analytic solution
        # is the only exact correctness oracle we have; use the prod* cases
        # below for anything performance-related.
        echo "grid nx=$NX x=1.0 y=1.0 z=1.0"
        echo "time t=$T"
        echo "supergrid gp=20"
        echo "twilight omega=6.28 phase=0.8 momega=6.28 errorlog=1"
        ;;
      cart-mr)
        # Cartesian + mesh refinement.  Every refinement interface sets
        # bRefInterface, which turns on a one-sided closure -- so the thin
        # Cartesian levels here are where the closure vectorisation pays.
        echo "grid nx=$NX x=1.0 y=1.0 z=1.0"
        echo "refinement zmax=0.35"
        echo "time t=$T"
        echo "supergrid gp=20"
        echo "twilight omega=6.28 phase=0.8 momega=6.28 errorlog=1"
        ;;
      curvi)
        # Topography: curvilinear4sg_ci.
        echo "grid nx=$NX x=1.0 y=1.0 z=1.0"
        echo "time t=$T"
        echo "topography input=gaussian zmax=0.25 order=4 gaussianAmp=0.05"
        echo "twilight omega=6.28 phase=0.8 momega=6.28 errorlog=1"
        ;;
      curvi-mr)
        # Topography + refinement: curvilinear4sgwind and the interface solve.
        echo "grid nx=$NX x=1.0 y=1.0 z=1.0"
        echo "time t=$T"
        echo "topography input=gaussian zmax=0.25 order=4 gaussianAmp=0.05"
        # SW4 requires >= 12 z-points (excluding ghosts) per grid. The Cartesian
        # level between the curvilinear bottom (0.25) and this interface must
        # therefore be at least 12*h deep, and h grows as --size shrinks, so
        # this is set for the coarsest size (S, h=0.025) and is comfortable at
        # the finer ones. At 0.15 the level gets 9 points and SW4 aborts.
        echo "refinement zmax=0.65"
        echo "twilight omega=6.28 phase=0.8 momega=6.28 errorlog=1"
        ;;
      prod|prod-curvi|prod-mr)
        # Production-shaped: point source, real material, and REAL supergrid
        # absorbing boundaries (a bare `boundary_conditions` line defaults to
        # bSuperGrid on five sides with a free surface on top -- see
        # processBoundaryConditions). This is the configuration the twilight
        # cases silently are not: with it, addSuperGridDamping runs and
        # evalRHS takes the rhs4th3fortsgstr_ci branch.
        #
        # Measured phase split on this case vs a twilight one:
        #   production: Div-stress 72%  SG 15.5%  Updates 8%   Forcing 2.7%
        #   twilight:   Div-stress 50%  SG  0.0%  Updates 6%   Forcing 40%
        # so every phase share quoted from a twilight case is against the
        # wrong denominator, and two phases production has are absent.
        #
        # dt/h = 2.5e-4 for vp=4000 over a 30 km cube (measured: nx=61,
        # t=1.0 -> 8 steps), hence t = STEPS * 7.5/(nx-1).
        echo "grid nx=$NX x=30000 y=30000 z=30000"
        echo "time t=$(python3 -c "print(f'{${STEPS}*7.5/(${NX}-1):.6g}')")"
        echo "boundary_conditions"
        echo "supergrid gp=20"
        [ "$n" = prod-curvi ] && echo "topography input=gaussian zmax=6000 order=4 gaussianAmp=1500 gaussianXc=15000 gaussianYc=15000 gaussianLx=6000 gaussianLy=6000"
        [ "$n" = prod-mr ] && echo "refinement zmax=12000"
        echo "block vp=4000 vs=2000 rho=2600"
        echo "source x=15000 y=15000 z=8000 mxy=1e18 t0=0.36 freq=16.6667 type=Gaussian"
        for r in 1 2 3; do
          echo "rec x=$((15000+r*2500)) y=$((15000+r*2000)) z=0 file=st0$r writeEvery=1000000"
        done
        ;;
      aniso)
        # nc=21 anisotropic path.  No analytic solution, so correctness is
        # checked by comparing station output rather than an error norm.
        echo "grid h=$(python3 -c "print(30000/$NX)") x=30000 y=30000 z=17000"
        echo "time t=$(python3 -c "print(round($T*10,3))")"
        echo "anisotropy"
        echo "ablock c11=9.72e10 c22=9.72e10 c33=9.72e10 c12=3.2404e10 c13=3.2404e10 c23=3.2404e10 c44=3.2398e10 c55=3.2398e10 c66=3.2398e10 rho=2700"
        echo "supergrid gp=10"
        echo "source x=15000 y=15000 z=2000 mxy=1e18 t0=0.36 freq=16.6667 type=Gaussian"
        echo "rec x=17600 y=17800 z=0 file=sta01 writeEvery=100000"
        ;;
      *) echo "unknown case: $n" >&2; return 1;;
    esac
    echo "fileio path=OUTPATH printcycle=100000"
    echo "developer reporttiming=1"
  } > "$f"
}

# ---------------------------------------------------------------------------
# build
# ---------------------------------------------------------------------------
build_side() {  # $1=label  $2=commit  $3=single|double
  local label="$1" sha="$2" prec="$3"
  local wt="$OUTDIR/wt-$label" bd="$OUTDIR/build-$label-$prec"
  if [[ ! -d "$wt" ]]; then
    git worktree add --detach -f "$wt" "$sha" >/dev/null
  fi
  local args=(-DCMAKE_BUILD_TYPE=Release -DBUILD_SW4MOPT=OFF)
  [[ "$prec" == "single" ]] && args+=(-DUSE_DOUBLE=OFF) || args+=(-DUSE_DOUBLE=ON)
  if [[ -f "$wt/cmake/SW4Optimization.cmake" ]]; then
    args+=(-DSW4_TARGET="$TARGET" -DSW4_STRICT_FP="$STRICT_FP")
  else
    # Pre-SW4Optimization checkout: no SW4_TARGET. Give it the closest
    # equivalent by hand so the comparison is not merely "new flags vs none".
    local a; a="$(arch_flags_for "$TARGET")"
    args+=(-DSW4_ARCH_FLAGS="$a" -DSW4_EXTRA_RELEASE_FLAGS="-funroll-loops")
    echo "  note: $sha predates SW4Optimization.cmake; using SW4_ARCH_FLAGS='$a'"
  fi
  [[ -n "$EXTRA_CMAKE" ]] && args+=($EXTRA_CMAKE)
  echo "  configuring $label/$prec ..."
  cmake -S "$wt" -B "$bd" "${args[@]}" > "$OUTDIR/logs/cmake-$label-$prec.log" 2>&1 \
    || { echo "CONFIGURE FAILED, see $OUTDIR/logs/cmake-$label-$prec.log" >&2; return 1; }
  echo "  building    $label/$prec (-j$JOBS) ..."
  local t0 t1; t0=$(date +%s)
  cmake --build "$bd" -j"$JOBS" > "$OUTDIR/logs/build-$label-$prec.log" 2>&1 \
    || { echo "BUILD FAILED, see $OUTDIR/logs/build-$label-$prec.log" >&2; return 1; }
  t1=$(date +%s)
  echo "    -> $((t1-t0))s   $(du -sh "$bd" | cut -f1)"
  grep -m1 'CXX_FLAGS' "$bd/CMakeFiles/sw4_common.dir/flags.make" 2>/dev/null \
    | sed "s|^|flags[$label/$prec] |" >> "$OUTDIR/platform.txt" || true
}

arch_flags_for() {  # best-effort -march for the pre-module checkout
  case "$1" in
    mn5-gpp|stampede3-spr) echo "-march=sapphirerapids -mtune=sapphirerapids -mprefer-vector-width=512";;
    cascade|hpc3-genoa)    echo "-march=znver4 -mtune=znver4";;
    hpc3-milan)            echo "-march=znver3 -mtune=znver3";;
    hpc3-portable)         echo "-march=x86-64-v3 -mtune=znver4";;
    frontera)              echo "-march=cascadelake -mtune=cascadelake";;
    stampede3)             echo "-march=skylake-avx512 -mtune=sapphirerapids";;
    vista)                 echo "-mcpu=neoverse-v2";;
    native)                echo "-march=native -mtune=native";;
    *)                     echo "";;
  esac
}

# ---------------------------------------------------------------------------
# run
# ---------------------------------------------------------------------------
run_one() {  # $1=build dir  $2=case file  $3=outpath -> stdout: the raw sw4 output
  local bd="$1" cf="$2" op="$3"
  rm -rf "$op"; mkdir -p "$op"
  sed "s|OUTPATH|$op|" "$cf" > "$op/run.in"
  ( cd "$op" && OMP_NUM_THREADS="$THREADS" OMP_PROC_BIND=close OMP_PLACES=cores \
      $LAUNCHER "$bd/bin/sw4" "$op/run.in" 2>&1 ) || true
}

# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------
mkdir -p "$OUTDIR/logs"
{
  echo "=== ab-bench $(date -u +%Y-%m-%dT%H:%M:%SZ) ==="
  echo "repo        $REPO"
  echo "base        $BASE_SHA  $(git log -1 --format=%s "$BASE")"
  echo "head        $HEAD_SHA  $(git log -1 --format=%s "$HEAD_REF")"
  echo "precisions  $PRECISIONS"
  echo "target      $TARGET     strict_fp=$STRICT_FP"
  echo "size        $SIZE (nx=$NX t=$T)   reps=$REPS"
  echo "ranks       $RANKS   threads=$THREADS   launcher='$LAUNCHER'"
  echo "cases       $CASES"
  echo
  echo "--- cpu ---"
  (lscpu 2>/dev/null | grep -iE 'model name|^cpu\(s\)|core\(s\) per socket|socket|thread\(s\) per core|^flags' | cut -c1-300) || true
  echo "--- memory ---"; (free -h 2>/dev/null | head -2) || true
  echo "--- compilers ---"
  for c in cc gcc g++ icpx gfortran mpicxx; do
    command -v $c >/dev/null && echo "$c: $($c --version 2>&1 | head -1)"
  done
  command -v cmake >/dev/null && echo "cmake: $(cmake --version | head -1)"
  echo "--- modules ---"; (module list 2>&1 | head -20) || echo "(none)"
  echo
} > "$OUTDIR/platform.txt"
sed -n '1,12p' "$OUTDIR/platform.txt"

if [[ "$DO_BUILD" == "1" ]]; then
  echo
  echo "== building =="
  for p in $PRECISIONS; do
    build_side base "$BASE_SHA" "$p"
    build_side head "$HEAD_SHA" "$p"
  done
fi
[[ "$DO_RUN" == "0" ]] && { echo; echo "built only; results dir: $OUTDIR"; exit 0; }

IFS=',' read -ra CASE_ARR <<< "$CASES"
for c in "${CASE_ARR[@]}"; do emit_case "$c"; done

FAILED_CASES=""
rm -f "$OUTDIR/failures.txt"
CSV="$OUTDIR/results.csv"
echo "case,precision,side,rep,total,divstress,forcing,bc,sg,comm,mr,img_tseries,updates,essi,wall,linf,l2" > "$CSV"

echo
echo "== running (A and B interleaved, $REPS reps) =="
for p in $PRECISIONS; do
  for c in "${CASE_ARR[@]}"; do
    printf "  %-10s %-6s " "$c" "$p"
    for r in $(seq 1 "$REPS"); do
      for side in base head; do
        bd="$OUTDIR/build-$side-$p"
        out="$(run_one "$bd" "$OUTDIR/cases/$c.in" "$OUTDIR/run/$c-$p-$side")"
        echo "$out" > "$OUTDIR/logs/run-$c-$p-$side-r$r.log"
        # phase row is the line after the column header
        # Every extraction below ends in `|| true`. Under pipefail a grep that
        # finds nothing returns 1, and set -e would then abort the whole run
        # before the summary is written -- which is exactly what a single bad
        # case did during development. A failing case must be recorded and
        # skipped, not fatal.
        phases="$(echo "$out" | grep -A1 'Total *Div-stress' | tail -1 \
                  | awk '{printf "%s,%s,%s,%s,%s,%s,%s,%s,%s,%s",$1,$2,$3,$4,$5,$6,$7,$8,$9,$10}' || true)"
        # SW4 switches format past a minute: "... phase 2 minutes 2.25e+01 seconds".
        # Taking $(NF-1) grabs only the seconds remainder and silently drops the
        # minutes, which inverted the A/B comparison on runs over 60s. Sum the
        # hours/minutes/seconds fields by name instead.
        wall="$(echo "$out" | grep 'time stepping phase' | tail -1 \
                | awk '{h=0;m=0;s=0;
                        for(i=1;i<=NF;i++){
                          if($i ~ /^hours?$/)        h=$(i-1);
                          else if($i ~ /^minutes?$/) m=$(i-1);
                          else if($i ~ /^seconds?$/) s=$(i-1);
                        }
                        printf "%.6f", h*3600+m*60+s}' || true)"
        linf="$(echo "$out" | grep -o 'Linf = *[0-9.eE+-]*' | tail -1 | awk '{print $NF}' || true)"
        l2="$(  echo "$out" | grep -o 'L2 = *[0-9.eE+-]*'   | tail -1 | awk '{print $NF}' || true)"
        # Production cases have no analytic solution, so hash the recorded
        # waveforms instead (see sac-hash.py for why the header is skipped).
        if [ -z "$linf" ]; then
          linf="$(python3 "$REPO/performance/ab-bench/sac-hash.py" \
                  "$OUTDIR/run/$c-$p-$side" 2>/dev/null || true)"
        fi
        if [[ -z "$phases" ]]; then
          FAILED_CASES="${FAILED_CASES}${FAILED_CASES:+ }$c/$p/$side"
          reason="$(echo "$out" | grep -m1 -iE 'precondition violated|error|abort' | cut -c1-140 || true)"
          echo "$c/$p/$side rep$r: no timing table. ${reason:-<no diagnostic>}" >> "$OUTDIR/failures.txt"
          printf "!"
        else
          printf "."
        fi
        echo "$c,$p,$side,$r,${phases:-,,,,,,,,,},${wall:-},${linf:-},${l2:-}" >> "$CSV"
      done
    done
    echo " done"
  done
done

python3 - "$CSV" "$OUTDIR/summary.txt" "$BASE_SHA" "$HEAD_SHA" <<'PYEOF'
import csv, sys, collections
csvp, outp, base_sha, head_sha = sys.argv[1:5]
COLS = ["total","divstress","forcing","bc","sg","comm","mr","img_tseries","updates","essi","wall"]
# forcing and bc are not touched by any change under test, so their measured
# A/B ratio IS the noise floor for that configuration. Reporting it turns
# "is 1.05x real?" from a judgement call into arithmetic.
UNTOUCHED = ["forcing", "bc"]
MIN_TRUSTWORTHY_SECONDS = 5.0

rows = list(csv.DictReader(open(csvp)))
def f(x):
    try: return float(x)
    except: return None

best = collections.defaultdict(dict)
allreps = collections.defaultdict(lambda: collections.defaultdict(list))
norm = collections.defaultdict(list)
for r in rows:
    k = (r["case"], r["precision"], r["side"])
    for c in COLS:
        v = f(r[c])
        if v is not None and v > 0:
            best[k][c] = min(best[k].get(c, 1e18), v)
            allreps[k][c].append(v)
    if r["linf"] or r["l2"]:
        norm[k].append((r["linf"], r["l2"]))

out = [f"A = {base_sha} (base)      B = {head_sha} (head)",
       "min-of-reps per phase, seconds.  speedup = A/B, >1 means B is faster.", ""]
hdr = f"{'case':<10} {'prec':<7} {'phase':<12} {'A':>10} {'B':>10} {'speedup':>8}"

for (case, prec) in sorted({(r['case'], r['precision']) for r in rows}):
    a = best.get((case,prec,'base'), {}); b = best.get((case,prec,'head'), {})
    if not a or not b: continue

    # --- noise floor from the phases nothing under test touches --------------
    # Every phase gets its OWN error bar, because noise scales with how big
    # the phase is: on a short run a 0.19s phase genuinely varies ~100% between
    # reps while a 2.6s phase varies a few percent. A single global "noise
    # floor" taken from the untouched phases was therefore either far too tight
    # (when they were large) or absurd (+/-99.7% when they were small), and in
    # both cases it mislabelled real results. The spread below is measured, not
    # assumed: max over both sides of (max-min)/min across reps.
    def spread_of(k, c):
        v = allreps[k].get(c, [])
        return (max(v)-min(v))/min(v) if len(v) > 1 and min(v) > 0 else None

    out.append(hdr); out.append("-"*len(hdr))
    for c in COLS:
        if c in a and c in b and a[c] > 1e-4:
            sp = a[c]/b[c]
            sa = spread_of((case,prec,'base'), c)
            sb = spread_of((case,prec,'head'), c)
            band = max([x for x in (sa, sb) if x is not None], default=None)
            if band is None:
                out.append(f"{case:<10} {prec:<7} {c:<12} {a[c]:>10.4g} {b[c]:>10.4g} {sp:>8.3f}")
            else:
                tag = "  <- within its own +/-%.0f%% spread" % (band*100) \
                      if abs(sp-1.0) <= band else "   (spread +/-%.0f%%)" % (band*100)
                out.append(f"{case:<10} {prec:<7} {c:<12} {a[c]:>10.4g} {b[c]:>10.4g} {sp:>8.3f}{tag}")

    # --- is the run even long enough to mean anything? -----------------------
    shortest = min([x for x in (a.get("total"), b.get("total")) if x] or [0])
    if shortest and shortest < MIN_TRUSTWORTHY_SECONDS:
        out.append(f"{'':<10} {'':<7} WARNING: shortest run {shortest:.2f}s is under "
                   f"{MIN_TRUSTWORTHY_SECONDS}s -- startup and jitter dominate. "
                   f"Re-run with a larger --size.")

    # --- correctness, tolerant of a non-deterministic threaded reduction -----
    na, nb = norm.get((case,prec,'base'), []), norm.get((case,prec,'head'), [])
    if na and nb:
        sa, sb = set(na), set(nb)
        def spread(vals):
            xs = [f(v[1]) for v in vals if f(v[1]) is not None]
            if len(xs) < 2 or not max(xs): return 0.0
            return (max(xs)-min(xs))/abs(max(xs))
        ra, rb = spread(na), spread(nb)
        if sa == sb and len(sa) == 1:
            out.append(f"{'':<10} {'':<7} correctness: BIT-EXACT (Linf and L2)")
        elif ra or rb:
            # Each side varies run to run: the error-norm reduction is threaded,
            # so in single precision the summation order is not reproducible.
            # Compare the difference between sides against that spread instead
            # of demanding equality it cannot deliver.
            la = [f(v[1]) for v in na if f(v[1]) is not None]
            lb = [f(v[1]) for v in nb if f(v[1]) is not None]
            gap = abs(sum(la)/len(la) - sum(lb)/len(lb)) / max(abs(sum(la)/len(la)), 1e-300)
            band = max(ra, rb)
            verdict = ("CONSISTENT within run-to-run spread" if gap <= band
                       else "DIFFERS BEYOND SPREAD -- investigate")
            out.append(f"{'':<10} {'':<7} correctness: {verdict} "
                       f"(A/B L2 gap {gap:.1e}, own spread {band:.1e})")
            out.append(f"{'':<10} {'':<7}   note: the L2/Linf reduction is threaded, so it is "
                       f"not bit-reproducible in single precision; re-run with "
                       f"--threads 1 for an exact check.")
        else:
            (la,l2a), (lb,l2b) = next(iter(sa)), next(iter(sb))
            try:
                d = max(abs(float(la)-float(lb))/abs(float(la)),
                        abs(float(l2a)-float(l2b))/abs(float(l2a)))
                out.append(f"{'':<10} {'':<7} correctness: DIFFERS, max rel {d:.1e} "
                           f"(A Linf={la} L2={l2a} / B Linf={lb} L2={l2b})")
            except Exception:
                out.append(f"{'':<10} {'':<7} correctness: DIFFERS A=({la},{l2a}) B=({lb},{l2b})")

    # --- wall vs total consistency (catches the minutes-format class of bug) --
    for side, d in (("A", a), ("B", b)):
        if "total" in d and "wall" in d and d["total"] > 1e-6:
            skew = abs(d["wall"] - d["total"]) / d["total"]
            if skew > 0.05:
                out.append(f"{'':<10} {'':<7} WARNING: {side} wall={d['wall']:.4g}s vs "
                           f"total={d['total']:.4g}s, {skew*100:.0f}% skew -- parse error, "
                           f"do not trust these figures")
    out.append("")

txt = "\n".join(out)
open(outp,"w").write(txt + "\n")
print(txt)
PYEOF

if [[ -n "$FAILED_CASES" ]]; then
  echo
  echo "WARNING: these runs produced no timing table and are absent from the"
  echo "summary: $FAILED_CASES"
  echo "See $OUTDIR/failures.txt and $OUTDIR/logs/."
fi

echo
echo "results:   $OUTDIR/summary.txt"
echo "raw csv:   $OUTDIR/results.csv"
echo "platform:  $OUTDIR/platform.txt"
echo
echo "To clean up the worktrees this created:"
echo "  git -C $REPO worktree remove --force $OUTDIR/wt-base"
echo "  git -C $REPO worktree remove --force $OUTDIR/wt-head"
