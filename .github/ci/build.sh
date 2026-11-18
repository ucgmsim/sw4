#!/usr/bin/env bash
# Configure and build SW4 with one compiler toolchain through CMake.
#
#   .github/ci/build.sh gcc|clang|intel|aocc [double|float] [build-dir]
#
# Run install-compiler.sh with the same compiler first.
#
# The OpenMPI wrappers stay the project compilers (CMakeLists.txt prefers them);
# OMPI_CC/OMPI_CXX/OMPI_FC point them at the toolchain under test, so CMake
# identifies e.g. mpicxx as IntelLLVM and applies that compiler's flags.
set -eu -o pipefail

compiler=${1:?usage: $0 gcc|clang|intel|aocc [double|float] [build-dir]}
precision=${2:-double}
builddir=${3:-build-$compiler-$precision}

case "$compiler" in
  gcc)
    export OMPI_CC=gcc OMPI_CXX=g++ OMPI_FC=gfortran ;;
  clang)
    export OMPI_CC=clang OMPI_CXX=clang++ OMPI_FC=gfortran ;;
  intel)
    # setvars.sh reads unset variables and returns non-zero on benign warnings.
    set +eu; source /opt/intel/oneapi/setvars.sh >/dev/null; set -eu
    export OMPI_CC=icx OMPI_CXX=icpx OMPI_FC=ifx ;;
  aocc)
    aocc_root=$(ls -d /opt/AMD/aocc-compiler-* | sort -V | tail -1)
    set +eu; source "$aocc_root/setenv_AOCC.sh" >/dev/null; set -eu
    export OMPI_CC=clang OMPI_CXX=clang++ OMPI_FC=flang ;;
  *)
    echo "unknown compiler: $compiler" >&2; exit 2 ;;
esac

case "$precision" in
  double) use_double=ON ;;
  float)  use_double=OFF ;;
  *) echo "unknown precision: $precision" >&2; exit 2 ;;
esac

echo "== $compiler: $(command -v "$OMPI_CXX")"
"$OMPI_CXX" --version | head -1
"$OMPI_FC" --version | head -1

cmake -S . -B "$builddir" -G Ninja \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_Fortran_COMPILER="$OMPI_FC" \
  -DSW4_TARGET=generic \
  -DUSE_DOUBLE="$use_double" \
  -DUSE_HDF5=ON -DUSE_FFTW3=ON -DUSE_PROJ=ON \
  -DHDF5_PREFER_PARALLEL=ON

# -k 0: report every file that fails, not just the first.
cmake --build "$builddir" --parallel "${CMAKE_BUILD_PARALLEL_LEVEL:-$(nproc)}" -- -k 0
