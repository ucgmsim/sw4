#!/usr/bin/env bash
# Install the libraries SW4 links against, plus one compiler toolchain.
#
#   .github/ci/install-compiler.sh gcc|clang|intel|aocc
#
# MPI, HDF5, FFTW and PROJ come from Ubuntu and are built with GCC; they expose
# C interfaces only, so every compiler below can link against them. The
# compiler-specific part is which C/C++/Fortran front ends the OpenMPI wrappers
# drive -- see ci-env.sh.
set -eu -o pipefail

compiler=${1:?usage: $0 gcc|clang|intel|aocc}
SUDO=$([ "$(id -u)" = 0 ] && echo "" || echo sudo)
export DEBIAN_FRONTEND=noninteractive

$SUDO apt-get update -q
$SUDO apt-get install -y -q --no-install-recommends \
    build-essential ca-certificates cmake curl gnupg gfortran ninja-build \
    libopenmpi-dev openmpi-bin libhdf5-openmpi-dev libfftw3-dev \
    libfftw3-mpi-dev libproj-dev liblapack-dev libblas-dev zlib1g-dev

case "$compiler" in
  gcc)
    ;;
  clang)
    # Fortran stays on gfortran: Ubuntu's flang is not yet a drop-in.
    $SUDO apt-get install -y -q --no-install-recommends clang libomp-dev clang-tidy
    ;;
  intel)
    curl -fsSL https://apt.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB \
      | gpg --dearmor | $SUDO tee /usr/share/keyrings/oneapi-archive-keyring.gpg >/dev/null
    echo "deb [signed-by=/usr/share/keyrings/oneapi-archive-keyring.gpg] https://apt.repos.intel.com/oneapi all main" \
      | $SUDO tee /etc/apt/sources.list.d/oneAPI.list >/dev/null
    $SUDO apt-get update -q
    $SUDO apt-get install -y -q --no-install-recommends \
      intel-oneapi-compiler-dpcpp-cpp intel-oneapi-compiler-fortran
    ;;
  aocc)
    ver=${AOCC_VERSION:-5.1.0}
    deb=aocc-compiler-${ver}_1_amd64.deb
    series=$(echo "$ver" | cut -d. -f1,2 | tr . -)   # 5.1.0 -> 5-1
    curl -fsSL -o "/tmp/$deb" \
      "https://download.amd.com/developer/eula/aocc/aocc-${series}/${deb}"
    $SUDO apt-get install -y -q "/tmp/$deb"
    rm -f "/tmp/$deb"
    ;;
  *)
    echo "unknown compiler: $compiler" >&2; exit 2 ;;
esac
