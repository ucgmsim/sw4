#!/bin/bash -e
#SBATCH --job-name=sw4-refine-el-2
#SBATCH --account=nesi99999
#SBATCH --time=01:00:00         # six runs of meshrefine/refine-el-2
#SBATCH --ntasks=4              # MPI ranks per run
#SBATCH --cpus-per-task=2       # OpenMP threads per rank
#SBATCH --mem-per-cpu=2G
#SBATCH --output=refine-el-2-%j.out
#SBATCH --partition=genoa
#SBATCH --nodes=1
#SBATCH --hint=nomultithread
#
# Run the meshrefine/refine-el-2 test case with every (platform, precision) build
# on the same node and report the final twilight errors and the execution times.
#
# Submit from the top of the sw4 source tree, where the build-<platform>-<precision>
# directories live (or set SW4_ROOT):
#
#     sbatch tools/run_refine_el_2.sl
#
# Each run gets its own directory refine-el-2-runs/<platform>-<precision>/ and the
# summary table is written to refine-el-2-<jobid>.md.

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
export OMP_PROC_BIND=close
export OMP_PLACES=cores
# Intel MPI (gimkl, intel): give each rank OMP_NUM_THREADS cores
export I_MPI_PIN_DOMAIN=omp

root=${SW4_ROOT:-$SLURM_SUBMIT_DIR}
infile="$root/pytest/reference/meshrefine/refine-el-2.in"
reffile="$root/pytest/reference/meshrefine/refine-el-2/TwilightErr.txt"
rundir="$root/refine-el-2-runs"
jobid=${SLURM_JOB_ID:-local}
tablefile="$root/refine-el-2-${jobid}.md"
ntasks=${SLURM_NTASKS:-4}

platforms="gimkl intel aocc"
precs="default single"

if [ ! -f "$infile" ]; then
    echo "input file $infile not found; submit from the sw4 source tree or set SW4_ROOT"
    exit 1
fi

# errInf and errL2 are the first two numbers on the last line of TwilightErr.txt
read_errors() {
    [ -f "$1" ] && tail -1 "$1" | awk '{ print $1, $2 }'
}

# "Execution time, solver|time stepping phase [h hours] [m minutes] s seconds" -> seconds
solver_seconds() {
    grep -hE "Execution time, (solver|time stepping) phase" "$1" 2>/dev/null | tail -1 | awk '{
        t = 0
        for (i = 1; i < NF; i++) {
            if ($(i+1) ~ /^hour/)   t += 3600 * $i
            if ($(i+1) ~ /^minute/) t += 60 * $i
            if ($(i+1) ~ /^second/) t += $i
        }
        printf "%.2f", t
    }'
}

mkdir -p "$rundir"
echo "node: $(hostname), ranks: $ntasks, threads/rank: $OMP_NUM_THREADS"

rows=()
for platform in $platforms; do
    for prec in $precs; do
        combo="${platform}-${prec}"
        builddir="$root/build-${combo}"
        sw4="$builddir/bin/sw4"
        if [ ! -x "$sw4" ]; then
            echo "skipping ${combo}: $sw4 not found"
            rows+=("| $platform | $prec | - | - | - | - | missing |")
            continue
        fi
        workdir="$rundir/$combo"
        rm -rf "$workdir"
        mkdir -p "$workdir"
        logfile="$workdir/sw4.out"

        echo "=========== running ${combo} ==========="
        status=ok
        # subshell: the module environment doesn't leak into the next build
        (
            cd "$workdir"
            source "$builddir/sourceme.sh" > /dev/null 2>&1
            if mpirun --version 2>&1 | grep -q "Open MPI"; then
                launch="mpirun -np $ntasks --map-by slot:PE=$OMP_NUM_THREADS --bind-to core"
            else
                launch="mpirun -np $ntasks"
            fi
            echo "$launch $sw4 $infile"
            start=$(date +%s.%N)
            rc=0
            $launch "$sw4" "$infile" > "$logfile" 2>&1 || rc=$?
            end=$(date +%s.%N)
            awk -v a="$start" -v b="$end" 'BEGIN { printf "%.2f\n", b - a }' > wall.txt
            exit $rc
        ) || status=failed

        wall=$(cat "$workdir/wall.txt" 2>/dev/null || echo -)
        solver=$(solver_seconds "$logfile")
        read -r errinf errl2 <<< "$(read_errors "$workdir/refine-el-2/TwilightErr.txt")"
        rows+=("| $platform | $prec | ${errinf:--} | ${errl2:--} | ${solver:--} | $wall | $status |")
        echo "${combo}: errInf=${errinf:--} errL2=${errl2:--} solver=${solver:--}s wall=${wall}s ($status)"
    done
done

read -r refinf refl2 <<< "$(read_errors "$reffile")"
{
    echo "meshrefine/refine-el-2 on $(hostname), job $jobid," \
         "$ntasks MPI ranks x $OMP_NUM_THREADS OpenMP threads"
    echo
    echo "| platform | precision | errInf | errL2 | time stepping (s) | wall (s) | status |"
    echo "|----------|-----------|--------|-------|-------------------|----------|--------|"
    echo "| reference | double | $refinf | $refl2 | - | - | - |"
    printf '%s\n' "${rows[@]}"
} > "$tablefile"

echo
cat "$tablefile"
echo
echo "table written to $tablefile"
