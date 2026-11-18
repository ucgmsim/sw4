#!/usr/bin/env bash
# Run clang-tidy (checks in .clang-tidy) on the lines changed since a base
# commit -- not on whole files, so existing code is left alone until touched.
#
#   .github/ci/static-analysis.sh [base-ref] [build-dir]
#
# base-ref defaults to the merge base with origin/single_precision. Run
# install-compiler.sh clang first; this only configures (for
# compile_commands.json), it does not build.
set -eu -o pipefail

base=${1:-$(git merge-base origin/single_precision HEAD)}
builddir=${2:-build-tidy}

tidy_diff=$(command -v clang-tidy-diff.py clang-tidy-diff \
              /usr/bin/clang-tidy-diff-*.py 2>/dev/null | head -1)
[ -n "$tidy_diff" ] || { echo "clang-tidy-diff not found" >&2; exit 2; }

# Plain clang, not the MPI wrappers: clang-tidy reads the compile commands
# but cannot ask mpicxx for its include paths, so FindMPI must spell them out.
export CC=clang CXX=clang++
cmake -S . -B "$builddir" -G Ninja -DCMAKE_EXPORT_COMPILE_COMMANDS=ON \
  -DCMAKE_BUILD_TYPE=Release -DCMAKE_Fortran_COMPILER=gfortran \
  -DSW4_TARGET=generic -DUSE_HDF5=ON -DUSE_FFTW3=ON -DUSE_PROJ=ON \
  -DHDF5_PREFER_PARALLEL=ON >/dev/null

# A changed header has no compile command of its own; clang-tidy borrows one
# from a similarly named file, which must not be a gfortran command.
python3 - "$builddir/compile_commands.json" <<'PY'
import json, sys
db = json.load(open(sys.argv[1]))
json.dump([e for e in db if not e["file"].lower().endswith((".f", ".f90"))],
          open(sys.argv[1], "w"), indent=1)
PY

echo "== clang-tidy on lines changed since $(git rev-parse --short "$base")"
git diff -U0 --no-color "$base" HEAD -- '*.C' '*.h' '*.cpp' '*.hpp' \
  | python3 -W ignore "$tidy_diff" -p1 -path "$builddir" -quiet \
      -iregex '.*\.(C|h|cpp|hpp)' -j "$(nproc)"
