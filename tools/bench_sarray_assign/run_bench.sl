#!/bin/bash -e
#SBATCH --job-name=bench-assign
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=4              # independent copies, to share memory bandwidth as the 4 sw4 ranks do
#SBATCH --cpus-per-task=2       # OpenMP threads per copy
#SBATCH --hint=nomultithread
#SBATCH --mem-per-cpu=4G
#SBATCH --output=bench-assign-%j.out
#
# Build bench_sarray_assign.C with classic icpc (intel/2022a), once without arch flags and
# once with AVX-512 (-march=icelake-server -qopt-zmm-usage=high), then run both on the same
# node: first one copy, then 4 copies at once (like lamb-3's 4 ranks x 2 threads).
# Add --account/--partition for your site, and adjust the module line:
#
#     sbatch [--account=... --partition=...] tools/bench_sarray_assign/run_bench.sl [n ...]

here=${SLURM_SUBMIT_DIR:-$PWD}/tools/bench_sarray_assign
sizes="$*"
ml purge > /dev/null 2>&1
ml intel/2022a > /dev/null 2>&1

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-2}
export OMP_PROC_BIND=close
export OMP_PLACES=cores

cd "$here"
icpc -O3 -DNDEBUG -qopenmp -o bench-generic bench_sarray_assign.C
icpc -O3 -DNDEBUG -march=icelake-server -qoverride-limits -qopt-zmm-usage=high -qopenmp \
     -o bench-avx512 bench_sarray_assign.C

echo "node: $(hostname)"
for flavour in generic avx512; do
    echo
    echo "=========== $flavour, 1 copy ==========="
    srun -n 1 --cpu-bind=cores ./bench-$flavour $sizes
    echo
    echo "=========== $flavour, 4 copies (rank 0 shown) ==========="
    srun -n 4 --cpu-bind=cores bash -c "./bench-$flavour $sizes > out-$flavour.\$SLURM_PROCID"
    cat out-$flavour.0
    grep -h MISMATCH out-$flavour.* || true
    rm -f out-$flavour.*
done
