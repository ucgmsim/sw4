#!/bin/bash
#SBATCH --job-name=sw4-prof
#SBATCH --nodes=1
#SBATCH --output=sw4-prof-%x-%j.out
#
# Sampling profile of SW4 on HPC3. Submit exactly like the benchmark jobs:
#
#   sbatch -A nesi00213 -p genoa -t 1:00:00 --ntasks=16 --cpus-per-task=4 \
#          --mem=96G --hint=nomultithread performance/ab-bench/slurm/hpc3-profile.sl
#
# Probes for perf / eu-stack / gdb / vendor profilers, picks the best the site
# actually permits, and records which one it used next to the numbers. It never
# hard-fails on a missing profiler: SW4's own per-phase timer needs no
# privileges and is always collected.
#
# PROBE_ONLY=1 reports capabilities and exits without running SW4 -- worth doing
# first, since it costs seconds and tells us whether perf is usable at all.
set -uo pipefail

CASE="${CASE:-curvi}"
SIZE="${SIZE:-L}"
STEPS="${STEPS:-100}"
HZ="${HZ:-199}"
PRECISION="${PRECISION:-single}"

have(){ command -v "$1" >/dev/null 2>&1; }
if [ -n "${MODULES:-}" ]; then module purge >/dev/null 2>&1; module load $MODULES >/dev/null 2>&1
elif command -v module >/dev/null 2>&1; then
  for m in "foss CMake" "foss/2026 CMake" "GCC OpenMPI OpenBLAS CMake"; do
    module purge >/dev/null 2>&1; module load $m >/dev/null 2>&1
    have cmake && have mpicxx && { echo "modules: $m"; break; }
  done
fi
# perf often lives in a module of its own on Cray/EB stacks
have perf || module load perf 2>/dev/null || module load linux-tools 2>/dev/null || true

cd "${SLURM_SUBMIT_DIR:-$PWD}" || exit 1
REPO="$(git rev-parse --show-toplevel 2>/dev/null)"; [ -z "$REPO" ] && { echo "not a git repo"; exit 1; }
cd "$REPO" || exit 1
PROF=performance/ab-bench/sample-profile.sh

echo "===================== capability probe ====================="
$PROF --probe
[ -n "${PROBE_ONLY:-}" ] && exit 0

case "${SLURM_JOB_PARTITION:-}" in
  *genoa*) TGT=hpc3-genoa;; *milan*) TGT=hpc3-milan;; *) TGT=hpc3-portable;;
esac
RANKS="${SLURM_NTASKS:-1}"; THREADS="${SLURM_CPUS_PER_TASK:-4}"
BD="$REPO/prof-build-$PRECISION"
echo; echo "===================== build ($TGT, $PRECISION) ====================="
# -g so the sampler can resolve symbols; RelWithDebInfo keeps the Release flags.
cmake -S "$REPO" -B "$BD" -DCMAKE_BUILD_TYPE=RelWithDebInfo -DBUILD_SW4MOPT=OFF \
      -DSW4_TARGET="$TGT" \
      -DUSE_DOUBLE=$([ "$PRECISION" = single ] && echo OFF || echo ON) \
      >/dev/null 2>&1 || { echo "configure failed"; exit 1; }
cmake --build "$BD" -j"${SLURM_CPUS_ON_NODE:-16}" >/dev/null 2>&1 || { echo "build failed"; exit 1; }

echo; echo "===================== case ====================="
OUT="$REPO/prof-${SLURM_JOB_ID:-$$}"; mkdir -p "$OUT"
case "$SIZE" in S) PPR=200000;; M) PPR=800000;; L) PPR=3200000;; XL) PPR=12800000;; esac
NX=$(python3 -c "print(round(($PPR*$RANKS)**(1/3))+1)")
T=$(python3 -c "print(f'{$STEPS*0.349/($NX-1):.6g}')")
{
  echo "grid nx=$NX x=1.0 y=1.0 z=1.0"
  echo "time t=$T"
  [ "$CASE" = curvi ] && echo "topography input=gaussian zmax=0.25 order=4 gaussianAmp=0.05"
  echo "supergrid gp=20"
  echo "twilight omega=6.28 phase=0.8 momega=6.28 errorlog=1"
  echo "fileio path=$OUT/run"
  echo "developer reporttiming=1"
} > "$OUT/case.in"
echo "  $CASE  nx=$NX t=$T  $RANKS ranks x $THREADS threads"

echo; echo "===================== profile ====================="
$PROF --bin "$BD/bin/sw4" --input "$OUT/case.in" --ranks "$RANKS" \
      --threads "$THREADS" --hz "$HZ" --outdir "$OUT"
echo; echo "results: $OUT"
