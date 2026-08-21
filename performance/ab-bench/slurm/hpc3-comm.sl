#!/bin/bash
#SBATCH --job-name=sw4-comm
#SBATCH --nodes=4
# Deliberately NOT --exclusive, and no --mem here. Both are passed at submit
# time, because the right choice depends on what you are measuring and on how
# long you are willing to queue:
#
#   Shared node, fast to schedule (recommended for a first look):
#     sbatch -p genoa --ntasks=2 --cpus-per-task=32 --mem-per-cpu=2G \\
#            --hint=nomultithread <this script>
#
#   Whole node, slow to schedule, needed for the bandwidth question:
#     sbatch -p genoa --exclusive --mem=0 --ntasks=2 --cpus-per-task=84 \\
#            --hint=nomultithread <this script>
#
# A shared-node run gives each core several times more memory bandwidth than a
# full node does, so for the memory-bound Cartesian kernel its speedups are an
# UPPER bound on what a production full-node run will show. The curvilinear
# kernel is compute-bound at every occupancy, so its numbers transfer directly.
# The job output records which mode it ran in so the caveat travels with the
# result.
#
# --hint=nomultithread keeps threads on physical cores; these nodes present 2
# hardware threads per core and letting OpenMP land on siblings halves the
# effective vector throughput.
#
# The script refuses to run an allocation below 8 CPUs -- see the guard below.
#SBATCH --time=02:00:00
#SBATCH --output=sw4-comm-%x-%j.out
##SBATCH --account=CHANGE_ME
#
# Multi-node A/B of the halo exchange on REANNZ HPC3.
#
#   sbatch -p genoa performance/ab-bench/slurm/hpc3-comm.sl
#
# This exists to answer one question the single-node job cannot. Commit 668bd37
# replaced four blocking MPI_Sendrecv calls with two phases of
# Irecv/Isend/Waitall, halving the serialisation points per halo exchange. On a
# single node that change is unmeasurable -- shared-memory MPI has no network
# latency to overlap, and the run-to-run spread in the Comm. column exceeded the
# difference. It was committed on the strength of being bit-exact and
# structurally better, with the performance claim explicitly unverified.
#
# So by default this compares ONLY that commit against its parent:
#     5374160  ivdep -> omp simd            (parent)
#     668bd37  non-blocking halo exchange   (the change under test)
# Nothing else differs between those two trees, so anything that shows up in the
# Comm. column is attributable.
#
# To instead get the whole series at multi-node scale:
#     BASE=23a3410 HEAD_REF=HEAD sbatch -p genoa hpc3-comm.sl
#
# Watch the Comm. column, not total. SW4 decomposes in i and j only, so halo
# traffic grows with rank count while the per-rank subdomain shrinks -- the
# effect is a larger fraction of runtime at high rank counts, which is exactly
# the regime this measures and a workstation cannot.

set -uo pipefail

BASE="${BASE:-5374160}"
HEAD_REF="${HEAD_REF:-668bd37}"
SIZE="${SIZE:-L}"
REPS="${REPS:-5}"
PRECISION="${PRECISION:-double}"
# Ranks per node: one per socket keeps NUMA sane. The rank COUNT is what matters
# here, so the sweep below varies nodes rather than packing more ranks per node.
RANKS_PER_NODE="${RANKS_PER_NODE:-2}"

NODES="${SLURM_JOB_NUM_NODES:-1}"
ALLOC="${SLURM_CPUS_ON_NODE:-0}"
PHYS="$( (command -v nproc >/dev/null && nproc) || echo 8)"
CORES="$ALLOC"; [ "$CORES" -lt 1 ] && CORES="$PHYS"
# Prefer what Slurm was actually asked for. Deriving THREADS from
# SLURM_CPUS_ON_NODE / ranks is wrong whenever --cpus-per-task is given: with
# SMT enabled (these nodes report 336 and 256 CPUs for 168 and 128 physical
# cores) that division can come out double the intended thread count and
# oversubscribe every rank.
if [ -n "${SLURM_CPUS_PER_TASK:-}" ]; then
  THREADS="$SLURM_CPUS_PER_TASK"
else
  THREADS=$(( CORES / RANKS_PER_NODE ))
fi
[ -n "${THREADS_OVERRIDE:-}" ] && THREADS="$THREADS_OVERRIDE"
[ "$THREADS" -lt 1 ] && THREADS=1

if [ -n "${SLURM_JOB_ID:-}" ] && [ "$ALLOC" -gt 0 ] && [ "$ALLOC" -lt 8 ]; then
  echo
  echo "REFUSING TO RUN: only $ALLOC CPU(s) allocated per node (physical: $PHYS)."
  echo "Re-submit stating the layout explicitly:"
  echo "  sbatch -p \${SLURM_JOB_PARTITION:-genoa} --nodes=$NODES \\"
  echo "         --ntasks-per-node=$RANKS_PER_NODE \\"
  echo "         --cpus-per-task=\$(( $PHYS / $RANKS_PER_NODE )) --mem=0 \\"
  echo "         performance/ab-bench/slurm/hpc3-comm.sl"
  echo "Set FORCE=1 to override."
  [ -z "${FORCE:-}" ] && exit 1
fi

echo "=================================================================="
echo " SW4 halo-exchange A/B"
echo " job        ${SLURM_JOB_ID:-<none>}  partition ${SLURM_JOB_PARTITION:-?}"
echo " nodes      $NODES x ${RANKS_PER_NODE} ranks x ${THREADS} threads"
echo " base       $BASE"
echo " head       $HEAD_REF"
echo "=================================================================="

case "${SLURM_JOB_PARTITION:-}" in
  *genoa*) TARGET=hpc3-genoa ;;
  *milan*) TARGET=hpc3-milan ;;
  *)       TARGET=hpc3-portable ;;
esac

have_toolchain() {
  command -v cmake >/dev/null && command -v mpicxx >/dev/null && \
  command -v gfortran >/dev/null && command -v git >/dev/null && \
  command -v python3 >/dev/null
}
if [ -n "${MODULES:-}" ]; then
  module purge >/dev/null 2>&1; module load $MODULES >/dev/null 2>&1
elif command -v module >/dev/null 2>&1; then
  for spec in "foss CMake" "foss/2023a CMake" "GCC OpenMPI OpenBLAS CMake" \
              "gcc openmpi openblas cmake" "cmake gcc openmpi"; do
    module purge >/dev/null 2>&1
    module load $spec >/dev/null 2>&1
    have_toolchain && { echo " modules    $spec"; break; }
  done
fi
have_toolchain || {
  echo "FAILED: no usable toolchain. Try MODULES='...' and see:"
  module avail 2>&1 | grep -iE 'cmake|gcc|openmpi|foss|openblas' | head -40
  exit 1; }

cd "${SLURM_SUBMIT_DIR:-$PWD}" || exit 1
REPO="$(git rev-parse --show-toplevel 2>/dev/null)" || true
[ -z "$REPO" ] && { echo "FAILED: .git missing -- the harness needs git worktree."; exit 1; }
cd "$REPO" || exit 1
for c in "$BASE" "$HEAD_REF"; do
  git rev-parse --verify --quiet "$c" >/dev/null || { echo "FAILED: no such commit $c"; exit 1; }
done

export OMP_PROC_BIND=close OMP_PLACES=cores

# One output directory for all rank counts, so the two builds happen once and
# every later round reuses them via --skip-build. Copying CMake build trees
# between directories would work for running (the binary is self-contained) but
# is wasteful and easy to get wrong; keeping one directory avoids the question.
OUT="$REPO/ab-comm-${SLURM_JOB_ID:-$$}"
mkdir -p "$OUT"
FIRST=1

for n in 1 2 4; do
  [ "$n" -gt "$NODES" ] && continue
  RANKS=$(( n * RANKS_PER_NODE ))
  echo
  echo "================ $n node(s), $RANKS ranks ================"
  ARGS=(--base "$BASE" --head "$HEAD_REF" --target "$TARGET"
        --size "$SIZE" --reps "$REPS" --precision "$PRECISION"
        --ranks "$RANKS" --threads "$THREADS" --jobs "$CORES"
        --cases cart,curvi --outdir "$OUT"
        --launcher "srun --nodes=$n --ntasks=$RANKS --cpus-per-task=$THREADS --cpu-bind=cores")
  [ "$FIRST" -eq 0 ] && ARGS+=(--skip-build)
  ./performance/ab-bench/ab-bench.sh "${ARGS[@]}"
  FIRST=0
  # results.csv and summary.txt are rewritten each round, so keep a copy
  for f in summary.txt results.csv failures.txt; do
    [ -f "$OUT/$f" ] && cp "$OUT/$f" "$OUT/${f%.*}-n${n}.${f##*.}"
  done
  echo "--- $n node(s) ---"
  [ -f "$OUT/summary.txt" ] && grep -E 'case|comm|total|correctness|-----' "$OUT/summary.txt"
done

echo
echo "Per-rank-count results kept as $OUT/summary-n{1,2,4}.txt"

echo
echo "=================================================================="
echo " done $(date -u +%Y-%m-%dT%H:%M:%SZ)"
echo " Compare the comm rows across the 1/2/4-node sections. If the"
echo " speedup there is ~1.0 at every rank count, the non-blocking"
echo " exchange in 668bd37 is not paying and the comm/compute overlap"
echo " restructure is not worth attempting either."
echo "=================================================================="
