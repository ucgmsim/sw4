#!/bin/bash
#SBATCH --job-name=sw4-comm
#SBATCH --nodes=4
#SBATCH --exclusive
#SBATCH --mem=0
# --mem=0 means "all memory on the node". Without it Slurm applies the site
# default (512M on NeSI), which is below what a single compile of
# rhs4th3fortc.C needs at --param=max-completely-peeled-insns=4000 (~513MB
# peak RSS) -- the build OOMs before any measurement happens.
#
# --exclusive alone does NOT get you the node's cores in the allocation on
# every Slurm configuration: it grants exclusive ACCESS while the allocation
# still reflects what was requested, so squeue shows CPUS=1 and
# SLURM_CPUS_ON_NODE reports 1. Pass the core count explicitly at submit time,
# because it differs per partition and cannot be baked in here:
#
#   sbatch -p genoa --ntasks=2 --cpus-per-task=84 ...   # 168-core nodes
#   sbatch -p milan --ntasks=2 --cpus-per-task=63 ...   # 126-core nodes
#
# The script cross-checks and refuses to run a degenerate allocation.
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
THREADS=$(( CORES / RANKS_PER_NODE )); [ "$THREADS" -lt 1 ] && THREADS=1

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
