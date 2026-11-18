#!/bin/bash -e
#SBATCH --job-name=sw4-tau
#SBATCH --account=nesi99999
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=4              # MPI ranks
#SBATCH --cpus-per-task=2       # OpenMP threads per rank, one physical core each
#SBATCH --hint=nomultithread    # no SMT: a CPU is a physical core
#SBATCH --mem-per-cpu=2G
#SBATCH --partition=genoa
#SBATCH --output=tau-%x-%j.out
#
# Profile one pytest/reference test case with TAU's tau_exec, using the sw4 of
# build-intel-default (intel/2022a, Intel MPI). No re-build is needed: tau_exec
# preloads TAU, which intercepts MPI and OpenMP (OMPT) and, with event-based
# sampling, attributes the samples to sw4's functions (TAU has BFD + libunwind).
#
# Submit from the top of the sw4 source tree (or set SW4_ROOT):
#
#     sbatch tools/profile_testcase_tau.sl                      # curvimeshrefine/gausshill-att-3
#     sbatch tools/profile_testcase_tau.sl meshrefine/refine-el-2
#
# Knobs (environment):
#     SW4_BUILD=intel-single     build-<suffix> to profile (must use the intel toolchain); a list,
#                                e.g. "intel-single intel-genoa-single", profiles each in turn on
#                                the same node, for a like-for-like comparison
#     TAU_EBS=0                  no sampling: MPI and OpenMP regions only (less overhead)
#     TAU_EBS_PERIOD=<n>         sampling period (tau_exec default 1000)
#     TAU_CALLPATH=1             also record callpaths (TAU_CALLPATH_DEPTH, default 2)
#
# The profile.<rank>.<ctx>.<thread> files land in <name>-tau/<jobid>/<build>/, next to
# sw4.out; a pprof summary is printed at the end of each run. To compare builds:
#     tools/compare_tau_profiles.py <name>-tau/<jobid>/* --labels ...
# For the GUI, on a node with X:
#     paraprof <name>-tau/<jobid>/<build>
# or pack it first and copy the .ppk anywhere:
#     cd <name>-tau/<jobid>/<build> && paraprof --pack <name>.ppk

testcase=${1:-curvimeshrefine/gausshill-att-3}
name=$(basename "$testcase")

tau_root=/nesi/project/nesi99999/pletzera/tau-intel
tau_bin=$tau_root/x86_64/bin
tau_makefile=$tau_root/x86_64/lib/Makefile.tau-icpc-ompt-mpi-pdt-openmp

root=${SW4_ROOT:-$SLURM_SUBMIT_DIR}
builds=${SW4_BUILD:-intel-single}
infile="$root/pytest/reference/${testcase}.in"
reffile="$root/pytest/reference/${testcase}/TwilightErr.txt"
jobid=${SLURM_JOB_ID:-local}
ntasks=${SLURM_NTASKS:-4}

[ -f "$infile" ] || { echo "input file $infile not found; submit from the sw4 source tree or set SW4_ROOT"; exit 1; }
[ -f "$tau_makefile" ] || { echo "$tau_makefile not found"; exit 1; }
for build in $builds; do
    [ -x "$root/build-$build/bin/sw4" ] || { echo "$root/build-$build/bin/sw4 not found"; exit 1; }
done

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-2}
export OMP_PROC_BIND=close
export OMP_PLACES=cores
# Intel MPI: give each rank OMP_NUM_THREADS cores (see tools/run_testcase.sl)
export I_MPI_PIN_DOMAIN=omp
export I_MPI_PIN_CELL=core

export TAU_MAKEFILE=$tau_makefile
export TAU_PROFILE=1
if [ "${TAU_CALLPATH:-0}" = 1 ]; then
    export TAU_CALLPATH=1
    export TAU_CALLPATH_DEPTH=${TAU_CALLPATH_DEPTH:-2}
fi

tau_opts="-T icpc,ompt,mpi,pdt,openmp -ompt"
if [ "${TAU_EBS:-1}" = 1 ]; then
    tau_opts="$tau_opts -ebs"
    [ -n "$TAU_EBS_PERIOD" ] && tau_opts="$tau_opts -ebs_period=$TAU_EBS_PERIOD"
fi
outpath=$(grep -o 'path=[^ ]*' "$infile" | head -1 | cut -d= -f2)

for build in $builds; do
# subshell: one build's module environment doesn't leak into the next
(
    builddir="$root/build-$build"
    sw4="$builddir/bin/sw4"
    rundir="$root/${name}-tau/${jobid}/${build}"
    export PROFILEDIR="$rundir"

    # intel/2022a: the MPI and compiler runtime TAU was built against
    source "$builddir/sourceme.sh" > /dev/null 2>&1
    export PATH="$tau_bin:$PATH"

    rm -rf "$rundir"
    mkdir -p "$rundir"
    cd "$rundir"
    echo
    echo "=========== $build ==========="
    echo "test case: $testcase, build: $build, node: $(hostname), ranks: $ntasks, threads/rank: $OMP_NUM_THREADS"
    echo "pinning: OMP_PROC_BIND=$OMP_PROC_BIND OMP_PLACES=$OMP_PLACES I_MPI_PIN_DOMAIN=$I_MPI_PIN_DOMAIN I_MPI_PIN_CELL=$I_MPI_PIN_CELL"
    echo "mpirun -np $ntasks tau_exec $tau_opts $sw4 $infile"

    start=$(date +%s.%N)
    rc=0
    mpirun -np "$ntasks" tau_exec $tau_opts "$sw4" "$infile" > sw4.out 2>&1 || rc=$?
    end=$(date +%s.%N)

    grep -E "Execution time, (solver|time stepping) phase" sw4.out | tail -1 || true
    awk -v a="$start" -v b="$end" 'BEGIN { printf "wall-clock time: %.2f s\n", b - a }'
    if [ $rc -ne 0 ]; then
        echo "sw4 failed (exit $rc); tail of $rundir/sw4.out:"
        tail -30 sw4.out
        exit $rc
    fi

    if [ -f "${outpath:-.}/TwilightErr.txt" ]; then
        echo "errors (errInf, errL2, solInf), this run:"
        cat "${outpath:-.}/TwilightErr.txt"
        echo "reference:"
        cat "$reffile"
    fi

    nprof=$(ls profile.* 2>/dev/null | wc -l)
    echo
    echo "$nprof TAU profile files in $rundir"
    [ "$nprof" -gt 0 ] || { echo "no profiles written; check sw4.out"; exit 1; }
    # total and mean over all ranks/threads, sorted by exclusive time
    pprof -s -m | head -60 > pprof-summary.txt
    echo "pprof summary (also in $rundir/pprof-summary.txt):"
    cat pprof-summary.txt
    pprof -m > pprof-full.txt
) || echo "${build}: profiling failed"
done
