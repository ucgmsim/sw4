#!/bin/bash
#SBATCH --job-name=sw4-ab
#SBATCH --nodes=1
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
CASES="${CASES:-}"                # override the case list, e.g. prod,prod-mr
THREADS_OVERRIDE="${THREADS_OVERRIDE:-}"

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

ALLOC="${SLURM_CPUS_ON_NODE:-0}"
# nproc reports the cgroup's CPU budget inside a Slurm job, not the node's core
# count, so on a 64-of-168 allocation it returned 4 and the occupancy line read
# "64 of 4 cores -- whole node". Ask Slurm for the node's real geometry and keep
# nproc only as a fallback outside a job.
PHYS=""
if command -v scontrol >/dev/null 2>&1 && [ -n "${SLURMD_NODENAME:-}" ]; then
  PHYS=$(scontrol show node "$SLURMD_NODENAME" 2>/dev/null | tr ' ' '\n' \
         | awk -F= '/^CoresPerSocket=/{c=$2} /^Sockets=/{s=$2} END{if(c&&s) print c*s}')
fi
[ -z "$PHYS" ] && PHYS="$( (command -v nproc >/dev/null && nproc) || echo 8)"
CORES="$ALLOC"; [ "$CORES" -lt 1 ] && CORES="$PHYS"
# Prefer what Slurm was actually asked for. Deriving THREADS from
# SLURM_CPUS_ON_NODE / ranks is wrong whenever --cpus-per-task is given: with
# SMT enabled (these nodes report 336 and 256 CPUs for 168 and 128 physical
# cores) that division can come out double the intended thread count and
# oversubscribe every rank.
if [ -n "${SLURM_CPUS_PER_TASK:-}" ]; then
  THREADS="$SLURM_CPUS_PER_TASK"
else
  THREADS=$(( CORES / RANKS ))
fi
[ -n "${THREADS_OVERRIDE:-}" ] && THREADS="$THREADS_OVERRIDE"
[ "$THREADS" -lt 1 ] && THREADS=1

# A benchmark that silently runs on one core is worse than one that refuses to
# start: the numbers look plausible and mean nothing. Bail out loudly instead.
if [ -n "${SLURM_JOB_ID:-}" ] && [ "$ALLOC" -gt 0 ] && [ "$ALLOC" -lt 8 ]; then
  echo
  echo "REFUSING TO RUN: Slurm allocated only $ALLOC CPU(s) on a node with"
  echo "$PHYS. --exclusive grants exclusive access but does not necessarily put"
  echo "the node's cores into the allocation, so timings here would be"
  echo "meaningless. Re-submit with the core count stated explicitly:"
  echo
  echo "  sbatch -p \${SLURM_JOB_PARTITION:-genoa} --ntasks=$RANKS \\"
  echo "         --cpus-per-task=\$(( $PHYS / $RANKS )) --mem=0 \\"
  echo "         performance/ab-bench/slurm/hpc3-ab.sl"
  echo
  echo "Set FORCE=1 to override and measure on $ALLOC CPU(s) anyway."
  [ -z "${FORCE:-}" ] && exit 1
fi
if [ "$ALLOC" -gt 0 ] && [ "$PHYS" -gt 0 ] && [ "$ALLOC" -lt "$PHYS" ]; then
  echo " occupancy   $ALLOC of $PHYS cores -- SHARED NODE"
  echo "             Each core has ~$(( PHYS / ALLOC ))x the memory bandwidth it"
  echo "             would get on a full node, so speedups for the"
  echo "             memory-bound Cartesian cases are an UPPER bound on what a"
  echo "             full-node production run will show. The curvilinear cases"
  echo "             are compute-bound at any occupancy and transfer directly."
  echo "             Neighbouring jobs also contend for bandwidth, which the"
  echo "             A/B interleaving and min-of-reps mitigate but cannot remove."
elif [ "$ALLOC" -gt 0 ] && [ "$PHYS" -gt 0 ]; then
  echo " occupancy   $ALLOC of $PHYS cores -- whole node (bandwidth-saturated;"
  echo "             these figures should transfer to production)"
else
  echo " occupancy   $ALLOC cpus allocated; node core count undetermined, so the"
  echo "             shared-vs-saturated caveat cannot be stated -- check by hand"
  echo "             with: scontrol show node \$SLURMD_NODENAME | grep -o 'CoresPerSocket=[0-9]*'"
fi
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
# An explicit MODULES= override is honoured FIRST and short-circuits the probe.
# It has to be here rather than after the loop: when the probe fails the script
# exits, so an override read afterwards would never run -- which is exactly the
# situation the override exists for.
if [ -n "${MODULES:-}" ] && command -v module >/dev/null 2>&1; then
  module purge >/dev/null 2>&1
  # shellcheck disable=SC2086
  module load $MODULES >/dev/null 2>&1
  if have_toolchain; then
    LOADED="$MODULES (from \$MODULES)"
  else
    echo "WARNING: MODULES='$MODULES' loaded but the toolchain is still"
    echo "         incomplete; falling through to the probe."
  fi
fi

if [ -z "$LOADED" ] && command -v module >/dev/null 2>&1; then
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
[ -n "$CASES" ] && ARGS+=(--cases "$CASES")

# Walltime sanity. The per-rank sizing means nx grows as cbrt(ranks), so
# SIZE=L at 16 ranks is a 372^3 grid -- four of six jobs in the previous round
# died two reps into curvi-mr at exactly that size, having built in 14-89s and
# completed the other three cases. Warn rather than fail: the estimate is crude
# and a big allocation may well be fine.
NCASES=$(echo "${CASES:-cart,cart-mr,curvi,curvi-mr}" | tr ',' '\n' | grep -c .)
EST=$(python3 -c "
ppr={'S':2e5,'M':8e5,'L':3.2e6,'XL':1.28e7}['$SIZE']
# 8.75 ns/pt/step per thread, x2 sides, x reps, x cases; mr cases cost ~4x
print(int(ppr*$STEPS*8.75e-9*2*$REPS*$NCASES*2.5))" 2>/dev/null || echo 0)
WALLSEC=$(python3 -c "
p='${SLURM_JOB_END_TIME:-0}'
import os,subprocess
try:
    o=subprocess.run(['squeue','-h','-j','${SLURM_JOB_ID:-0}','-o','%L'],capture_output=True,text=True).stdout.strip()
    h,m,sec=0,0,0
    parts=o.replace('-',':').split(':')
    parts=[int(x) for x in parts if x.isdigit()]
    print(sum(v*f for v,f in zip(reversed(parts),[1,60,3600,86400])))
except Exception: print(0)" 2>/dev/null || echo 0)
if [ "${EST:-0}" -gt 0 ] && [ "${WALLSEC:-0}" -gt 0 ] && [ "$EST" -gt "$WALLSEC" ]; then
  echo
  echo " WARNING: rough estimate ${EST}s of run time against ${WALLSEC}s of walltime"
  echo "          left. SIZE=$SIZE at $RANKS ranks gives nx~$(python3 -c "print(round((${EST}/1)**0))" 2>/dev/null; true)"
  echo "          Consider SIZE=M, fewer --reps, or a shorter CASES list. The"
  echo "          harness writes results.csv incrementally, so a job killed on"
  echo "          time still yields every case that finished."
  echo
fi

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
