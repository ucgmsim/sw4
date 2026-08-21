#!/bin/bash
#SBATCH --job-name=sw4-ab
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --time=03:00:00
#SBATCH --output=sw4-ab-%x-%j.out
##SBATCH --account=CHANGE_ME          # uncomment and set if your site requires it
#
# Single-node A/B benchmark of the 2026-08-21 optimisation series on REANNZ HPC3.
#
#   sbatch -p genoa  performance/ab-bench/slurm/hpc3-ab.sl
#   sbatch -p milan  performance/ab-bench/slurm/hpc3-ab.sl
#
# The partition comes from -p and the script derives everything else from it,
# because the two are NOT interchangeable: genoa is Zen 4 with AVX-512, milan is
# Zen 3 with AVX2 only, and a -march=znver4 binary raises SIGILL on milan. Node
# geometry is read from Slurm rather than hardcoded (genoa 168 cores, milan 126).
#
# BEFORE SUBMITTING, make sure the .git directory came across. The harness uses
# `git worktree` to check out both commits, so a source-only copy cannot work:
#   rsync -az --exclude 'build*' --exclude 'ab-results-*' sw4/ hpc3:~/sw4-ab/
#
# Everything this script discovers is echoed into the job output, so if it fails
# the log should say why without needing a second attempt.

set -uo pipefail       # deliberately NOT -e: module systems return nonzero for
                       # benign reasons and we want to reach the diagnostics.

# ---------------------------------------------------------------- knobs -----
SIZE="${SIZE:-L}"                 # S | M | L | XL
REPS="${REPS:-5}"
PRECISION="${PRECISION:-double}"  # single | double | both
BASE="${BASE:-23a3410}"           # last commit before the series
HEAD_REF="${HEAD_REF:-HEAD}"
RANKS="${RANKS:-2}"               # one per socket; these nodes are dual-socket
STRICT="${STRICT:-}"              # set to 1 for a bit-reproducible comparison

echo "=================================================================="
echo " SW4 A/B benchmark"
echo " job         ${SLURM_JOB_ID:-<none>} on ${SLURMD_NODENAME:-$(hostname)}"
echo " partition   ${SLURM_JOB_PARTITION:-<unknown>}"
echo " started     $(date -u +%Y-%m-%dT%H:%M:%SZ)"
echo "=================================================================="

# ------------------------------------------------------- target selection ---
case "${SLURM_JOB_PARTITION:-}" in
  *genoa*) TARGET=hpc3-genoa ;;
  *milan*) TARGET=hpc3-milan ;;
  *)       TARGET=hpc3-portable
           echo "NOTE: partition '${SLURM_JOB_PARTITION:-unset}' not recognised."
           echo "      Falling back to the x86-64-v3 build, which runs on both"
           echo "      partitions but leaves AVX-512 unused on genoa. Pass"
           echo "      -p genoa or -p milan to get a tuned build." ;;
esac

CORES="${SLURM_CPUS_ON_NODE:-$( (command -v nproc >/dev/null && nproc) || echo 8)}"
THREADS=$(( CORES / RANKS )); [ "$THREADS" -lt 1 ] && THREADS=1
echo " target      $TARGET"
echo " geometry    ${CORES} cores -> ${RANKS} ranks x ${THREADS} threads"
echo " case size   $SIZE   reps=$REPS   precision=$PRECISION"
echo

# ------------------------------------------------------------- modules ------
# Site module names are not knowable from here, so try the plausible sets in
# order and verify by looking for the tools rather than trusting exit codes.
have_toolchain() {
  command -v cmake    >/dev/null && \
  command -v mpicxx   >/dev/null && \
  command -v gfortran >/dev/null && \
  command -v git      >/dev/null && \
  command -v python3  >/dev/null
}

CANDIDATES=(
  "foss CMake"
  "foss/2023a CMake"
  "foss/2022a CMake"
  "GCC OpenMPI OpenBLAS CMake"
  "gcc openmpi openblas cmake"
  "GCCcore OpenMPI OpenBLAS CMake"
  "cmake gcc openmpi"
)

LOADED=""
if command -v module >/dev/null 2>&1; then
  for spec in "${CANDIDATES[@]}"; do
    module purge >/dev/null 2>&1
    # shellcheck disable=SC2086
    module load $spec >/dev/null 2>&1
    if have_toolchain; then LOADED="$spec"; break; fi
  done
  if [ -z "$LOADED" ]; then
    module purge >/dev/null 2>&1
    have_toolchain && LOADED="(none needed - tools already on PATH)"
  fi
else
  have_toolchain && LOADED="(no module command; tools already on PATH)"
fi

if [ -z "$LOADED" ]; then
  echo "FAILED: could not assemble a toolchain (need cmake, mpicxx, gfortran,"
  echo "git, python3). Tried: ${CANDIDATES[*]}"
  echo
  echo "--- module avail (grep'd) ---------------------------------------"
  module avail 2>&1 | grep -iE 'cmake|gcc|openmpi|foss|openblas|intel|oneapi' | head -60
  echo "----------------------------------------------------------------"
  echo "Re-submit with the right names, e.g.:"
  echo "  MODULES='GCC/13.2.0 OpenMPI/4.1.6 CMake/3.27' sbatch -p genoa hpc3-ab.sl"
  exit 1
fi
echo " modules     $LOADED"
[ -n "${MODULES:-}" ] && { module purge >/dev/null 2>&1; module load $MODULES; echo " modules     (overridden) $MODULES"; }

echo
echo "--- toolchain ----------------------------------------------------"
for t in cmake mpicxx gfortran git python3; do
  printf "  %-9s %s\n" "$t" "$(command -v $t) :: $($t --version 2>&1 | head -1)"
done
echo "--- cpu ----------------------------------------------------------"
lscpu 2>/dev/null | grep -iE 'model name|^cpu\(s\)|socket|core\(s\) per socket|thread\(s\) per core' | sed 's/^/  /'
echo "  avx512: $(grep -o 'avx512[a-z]*' /proc/cpuinfo 2>/dev/null | sort -u | tr '\n' ' ')"
echo "------------------------------------------------------------------"
echo

# --------------------------------------------------------------- checks -----
cd "${SLURM_SUBMIT_DIR:-$PWD}" || exit 1
REPO="$(git rev-parse --show-toplevel 2>/dev/null)"
if [ -z "$REPO" ]; then
  echo "FAILED: not inside a git repository, or .git was not copied across."
  echo "The harness uses 'git worktree' to check out both commits, so a"
  echo "source-only copy cannot work. Re-sync including .git:"
  echo "  rsync -az --exclude 'build*' --exclude 'ab-results-*' sw4/ <host>:~/sw4-ab/"
  exit 1
fi
cd "$REPO" || exit 1
for c in "$BASE" "$HEAD_REF"; do
  git rev-parse --verify --quiet "$c" >/dev/null || {
    echo "FAILED: commit '$c' not present in this clone."
    echo "rsync may have excluded git objects, or the branch was not fetched."
    exit 1; }
done
echo " repo        $REPO"
echo " base        $(git rev-parse --short "$BASE")  $(git log -1 --format=%s "$BASE")"
echo " head        $(git rev-parse --short "$HEAD_REF")  $(git log -1 --format=%s "$HEAD_REF")"
echo

# ------------------------------------------------------------------ run -----
OUT="$REPO/ab-${SLURM_JOB_PARTITION:-local}-${SLURM_JOB_ID:-$$}"
ARGS=(--base "$BASE" --head "$HEAD_REF" --target "$TARGET"
      --size "$SIZE" --reps "$REPS" --precision "$PRECISION"
      --ranks "$RANKS" --threads "$THREADS"
      --jobs "$CORES" --outdir "$OUT")
[ -n "$STRICT" ] && ARGS+=(--strict-fp)

# Prefer srun inside a Slurm allocation; the harness detects this itself, but
# being explicit avoids surprises if the site wraps mpirun.
# Be explicit about the layout. A bare `srun -n N` inherits cpus-per-task=1
# from the allocation on many sites, which pins each rank's OpenMP threads onto
# a single core and makes the whole measurement meaningless.
if command -v srun >/dev/null; then
  ARGS+=(--launcher "srun --ntasks=$RANKS --cpus-per-task=$THREADS --cpu-bind=cores")
fi
export OMP_PROC_BIND=close OMP_PLACES=cores

echo "+ ./performance/ab-bench/ab-bench.sh ${ARGS[*]}"
echo
./performance/ab-bench/ab-bench.sh "${ARGS[@]}"
RC=$?

echo
echo "=================================================================="
echo " finished $(date -u +%Y-%m-%dT%H:%M:%SZ)   exit=$RC"
echo " results  $OUT/summary.txt"
echo "=================================================================="
echo
[ -f "$OUT/summary.txt" ] && { echo "--- summary ---"; cat "$OUT/summary.txt"; }
[ -f "$OUT/failures.txt" ] && { echo "--- failures ---"; cat "$OUT/failures.txt"; }

# Leave the worktrees in place so a follow-up job can reuse the builds with
# --skip-build; remove them when you are done with $OUT:
echo
echo "To reclaim the worktrees:"
echo "  git -C $REPO worktree remove --force $OUT/wt-base"
echo "  git -C $REPO worktree remove --force $OUT/wt-head"
exit $RC
