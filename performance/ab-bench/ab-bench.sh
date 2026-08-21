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
  --size S|M|L|XL    case size                 (default: $SIZE)
  --reps N           timed repetitions, interleaved         (default: $REPS)
  --ranks N          MPI ranks                 (default: $RANKS)
  --threads N        OMP_NUM_THREADS           (default: cores/ranks)
  --cases LIST       comma list from: cart,cart-mr,curvi,curvi-mr,aniso
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
case "$SIZE" in
  S)  NX=41;  T=0.10;;
  M)  NX=81;  T=0.15;;
  L)  NX=161; T=0.15;;
  XL) NX=241; T=0.20;;
  *)  echo "bad --size: $SIZE" >&2; exit 2;;
esac

emit_case() {  # $1=name  -> writes $OUTDIR/cases/$1.in
  local n="$1" f="$OUTDIR/cases/$1.in"
  mkdir -p "$OUTDIR/cases"
  {
    case "$n" in
      cart)
        # Cartesian + supergrid: rhs4th3fortsgstr_ci interior.
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
        echo "refinement zmax=0.15"
        echo "twilight omega=6.28 phase=0.8 momega=6.28 errorlog=1"
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
        phases="$(echo "$out" | grep -A1 'Total *Div-stress' | tail -1 \
                  | awk '{printf "%s,%s,%s,%s,%s,%s,%s,%s,%s,%s",$1,$2,$3,$4,$5,$6,$7,$8,$9,$10}')"
        wall="$(echo "$out" | awk '/time stepping phase/{print $(NF-1)}' | tail -1)"
        linf="$(echo "$out" | grep -o 'Linf = *[0-9.eE+-]*' | tail -1 | awk '{print $NF}')"
        l2="$(  echo "$out" | grep -o 'L2 = *[0-9.eE+-]*'   | tail -1 | awk '{print $NF}')"
        echo "$c,$p,$side,$r,${phases:-,,,,,,,,,},${wall:-},${linf:-},${l2:-}" >> "$CSV"
      done
      printf "."
    done
    echo " done"
  done
done

python3 - "$CSV" "$OUTDIR/summary.txt" "$BASE_SHA" "$HEAD_SHA" <<'PYEOF'
import csv, sys, collections
csvp, outp, base_sha, head_sha = sys.argv[1:5]
COLS = ["total","divstress","forcing","bc","sg","comm","mr","img_tseries","updates","essi","wall"]
rows = list(csv.DictReader(open(csvp)))
def f(x):
    try: return float(x)
    except: return None
best = collections.defaultdict(dict)     # (case,prec,side) -> col -> min
norm = collections.defaultdict(set)      # (case,prec,side) -> {(linf,l2)}
for r in rows:
    k = (r["case"], r["precision"], r["side"])
    for c in COLS:
        v = f(r[c])
        if v is not None and v > 0:
            best[k][c] = min(best[k].get(c, 1e18), v)
    if r["linf"] or r["l2"]:
        norm[k].add((r["linf"], r["l2"]))

out = []
out.append(f"A = {base_sha} (base)      B = {head_sha} (head)")
out.append("min-of-reps per phase, seconds.  speedup = A/B, >1 means B is faster.")
out.append("")
hdr = f"{'case':<10} {'prec':<7} {'phase':<12} {'A':>10} {'B':>10} {'speedup':>8}"
for (case, prec) in sorted({(r['case'], r['precision']) for r in rows}):
    a, b = best.get((case,prec,'base'),{}), best.get((case,prec,'head'),{})
    if not a or not b: continue
    out.append(hdr); out.append("-"*len(hdr))
    for c in COLS:
        if c in a and c in b and a[c] > 1e-4:
            out.append(f"{case:<10} {prec:<7} {c:<12} {a[c]:>10.4g} {b[c]:>10.4g} {a[c]/b[c]:>8.3f}")
    # correctness
    na, nb = norm.get((case,prec,'base'), set()), norm.get((case,prec,'head'), set())
    if na and nb:
        if len(na) > 1 or len(nb) > 1:
            verdict = "NON-DETERMINISTIC across reps (thread reduction order?)"
        elif na == nb:
            verdict = "BIT-EXACT (Linf and L2)"
        else:
            (la,l2a), (lb,l2b) = next(iter(na)), next(iter(nb))
            try:
                d = max(abs(float(la)-float(lb))/abs(float(la)),
                        abs(float(l2a)-float(l2b))/abs(float(l2a)))
                verdict = f"DIFFERS, max rel {d:.1e}  (A Linf={la} L2={l2a} / B Linf={lb} L2={l2b})"
            except Exception:
                verdict = f"DIFFERS  A=({la},{l2a})  B=({lb},{l2b})"
        out.append(f"{'':<10} {'':<7} correctness: {verdict}")
    out.append("")
txt = "\n".join(out)
open(outp,"w").write(txt + "\n")
print(txt)
PYEOF

echo
echo "results:   $OUTDIR/summary.txt"
echo "raw csv:   $OUTDIR/results.csv"
echo "platform:  $OUTDIR/platform.txt"
echo
echo "To clean up the worktrees this created:"
echo "  git -C $REPO worktree remove --force $OUTDIR/wt-base"
echo "  git -C $REPO worktree remove --force $OUTDIR/wt-head"
