#!/usr/bin/env bash
#
# Reproducibly build SW4 from this fork using a private, self-contained Spack.
#
# Intended usage on a cluster:
#
#   git clone git@github.com:ucgmsim/sw4.git && cd sw4
#   sbatch/qsub a job whose payload is:  ./install.sh
#
# ...and come back to spack-env/.spack-env/view/bin/sw4.
#
# What it does:
#   0. Clones Spack at a pinned tag into --prefix. ~/.spack is ignored entirely,
#      so the result does not depend on the invoking user's Spack config.
#   1. Discovers the host toolchain: compilers via `spack compiler find`, MPI via
#      the vendor's own query tools, optionally after loading environment modules
#      with --modules. Writes what it found to spack-env/host/host.yaml.
#   2. Concretizes and installs spack-env/, which builds SW4 and every dependency
#      except MPI with the host compiler.
#
# Nothing host-specific is hardcoded; it is all discovered or overridable by flag.
# Two things come from the host by design -- the compiler and MPI -- plus a list
# of build-only tools (see BUILD_TOOLS).
#
set -euo pipefail

readonly REPO_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"

# --- pinned versions ---------------------------------------------------------
# Spack itself and its package recipes are pinned so that a rebuild months from
# now resolves the same dependency versions. The recipe pin lives in
# spack-env/spack.yaml (repos: builtin: commit:); keep the two consistent.
readonly SPACK_GIT="${SPACK_GIT:-https://github.com/spack/spack.git}"
readonly SPACK_REF="${SPACK_REF:-v1.2.2}"
# Commit of https://github.com/spack/spack-packages that supplies the builtin
# recipes. Must match `repos: builtin: commit:` in spack-env/spack.yaml.
readonly SPACK_PACKAGES_COMMIT="${SPACK_PACKAGES_COMMIT:-119680aeee8ea802c6111b7167583bddef97e82f}"

# Build-only tools taken from the host rather than built by Spack. None of these
# contribute code to the sw4 binary -- they are programs Spack runs during builds
# -- so using the host's copies changes nothing about what gets linked, and it
# avoids a real failure mode: Spack's pinned m4@1.4.20 and tar@1.35 bundle a
# gnulib that does not compile against glibc >= 2.42 (it collides with glibc's
# _Generic bsearch), and the pinned recipes offer nothing newer. Host m4 1.4.21+
# has the fix. Deliberately excluded is anything that can end up linked:
# xz/bzip2/zlib/zstd/ncurses/gettext/libxml2/openssl/sqlite/curl.
readonly BUILD_TOOLS=(m4 tar gmake autoconf automake libtool diffutils findutils
                      pkgconf perl texinfo bison flex cmake ninja)

# A commit built via `sw4@git.<sha>=<version>` needs a declared version to sort
# against. Must be a version declared in the sw4 recipe.
readonly SW4_REF_VERSION="2026.07.30"

# --- defaults ----------------------------------------------------------------
PREFIX="${SW4_SPACK_PREFIX:-$HOME/sw4-spack}"
JOBS="$(getconf _NPROCESSORS_ONLN 2>/dev/null || echo 4)"
MODULES=""
MPI_SPEC=""
MPI_PREFIX=""
MPI_PKG=""
COMPILER_SPEC=""
COMPILER_PKG=""
INSTALL_ROOT=""
BUILD_STAGE=""
SW4_SPEC_VERSION=""
ALLOW_UNPUSHED=0
CONCRETIZE_ONLY=0
NO_HOST_BUILD_TOOLS=0
SW4_VARIANTS="${SW4_VARIANTS:-}"
SW4_VERSION_CONSTRAINT=""

usage() {
    cat <<EOF
Reproducibly build SW4 from this fork using a private, self-contained Spack.

The host compiler and MPI are auto-detected and used as externals; every other
dependency (PROJ, HDF5, FFTW, OpenBLAS, ...) is built by Spack.

Usage: ./install.sh [options]

  --prefix DIR         Where to put Spack and its install tree
                       (default: \$HOME/sw4-spack, or \$SW4_SPACK_PREFIX)
  --jobs N             Parallel build jobs (default: all cores, here $JOBS)
  --modules "A B .."   Environment modules to load before host discovery,
                       e.g. --modules "gcc/13.2.0 openmpi/4.1.6"
  --compiler SPEC      Compiler to build everything with, e.g. --compiler gcc@13.2.0.
                       Default: the newest GCC that \`spack compiler find\` locates.
  --mpi-spec SPEC      Skip MPI auto-detection, e.g. --mpi-spec "cray-mpich@8.1.29"
  --mpi-prefix DIR     Prefix for --mpi-spec (default: derived from mpicc)
  --variants "..."     Extra variants for the sw4 spec,
                       e.g. --variants "precision=single +native ~fftw"
  --sw4-version V      Build a declared package version (e.g. $SW4_REF_VERSION) instead
                       of the commit currently checked out here
  --allow-unpushed     Proceed even if HEAD is not on any remote branch. Spack
                       fetches from the remote, so the build will still fail
                       unless the commit is reachable there.
  --no-host-build-tools
                       Build m4/tar/perl/... with Spack instead of using the
                       host's. Note that Spack's pinned m4 and tar do not
                       compile against glibc >= 2.42.
  --install-root DIR   Spack install tree location (default: PREFIX/spack/opt/spack)
  --build-stage DIR    Scratch space for builds (default: Spack's own choice).
                       Point this at node-local or fast scratch on a cluster.
  --concretize-only    Stop after writing spack.lock; do not build
  -h, --help           This message
EOF
}

log()  { printf '\033[1;34m==>\033[0m %s\n' "$*"; }
warn() { printf '\033[1;33m==> WARNING:\033[0m %s\n' "$*" >&2; }
die()  { printf '\033[1;31m==> ERROR:\033[0m %s\n' "$*" >&2; exit 1; }

while [[ $# -gt 0 ]]; do
    case "$1" in
        --prefix)         PREFIX="$2"; shift 2 ;;
        --jobs|-j)        JOBS="$2"; shift 2 ;;
        --modules)        MODULES="$2"; shift 2 ;;
        --compiler)       COMPILER_SPEC="$2"; shift 2 ;;
        --mpi-spec)       MPI_SPEC="$2"; shift 2 ;;
        --mpi-prefix)     MPI_PREFIX="$2"; shift 2 ;;
        --variants)       SW4_VARIANTS="$2"; shift 2 ;;
        --sw4-version)    SW4_SPEC_VERSION="$2"; shift 2 ;;
        --allow-unpushed) ALLOW_UNPUSHED=1; shift ;;
        --no-host-build-tools) NO_HOST_BUILD_TOOLS=1; shift ;;
        --install-root)   INSTALL_ROOT="$2"; shift 2 ;;
        --build-stage)    BUILD_STAGE="$2"; shift 2 ;;
        --concretize-only) CONCRETIZE_ONLY=1; shift ;;
        -h|--help)        usage; exit 0 ;;
        *)                die "unknown option: $1 (try --help)" ;;
    esac
done

readonly SPACK_ROOT_DIR="$PREFIX/spack"
readonly ENV_DIR="$REPO_DIR/spack-env"
readonly HOST_CFG_DIR="$ENV_DIR/host"

[[ -f "$ENV_DIR/spack.yaml" ]] || die "$ENV_DIR/spack.yaml not found -- run this from a checkout of the sw4 fork"

# Keep the user's ~/.spack out of it: a reproducible build must not depend on
# whatever externals and preferences happen to live there.
export SPACK_DISABLE_LOCAL_CONFIG=1
export SPACK_USER_CACHE_PATH="$PREFIX/user-cache"

# -----------------------------------------------------------------------------
# 0. Bootstrap Spack
# -----------------------------------------------------------------------------
bootstrap_spack() {
    if [[ -d "$SPACK_ROOT_DIR/.git" ]]; then
        local have
        have="$(git -C "$SPACK_ROOT_DIR" describe --tags --exact-match 2>/dev/null \
                || git -C "$SPACK_ROOT_DIR" rev-parse --short HEAD)"
        log "Reusing existing Spack at $SPACK_ROOT_DIR ($have)"
    else
        log "Cloning Spack $SPACK_REF into $SPACK_ROOT_DIR"
        mkdir -p "$PREFIX"
        git clone --depth 1 --branch "$SPACK_REF" "$SPACK_GIT" "$SPACK_ROOT_DIR"
    fi

    # setup-env.sh trips over `set -u`.
    set +u
    # shellcheck disable=SC1091
    source "$SPACK_ROOT_DIR/share/spack/setup-env.sh"
    set -u

    command -v spack >/dev/null || die "sourcing setup-env.sh did not put spack on PATH"
    log "Using $(spack --version)"
}

# Spack clones the builtin recipe repo on first use. `repos: builtin: commit:` in
# spack.yaml records the intent, but Spack does not reliably honour it on that
# initial clone -- it lands on the tip of releases/v2026.06 instead. That is not
# cosmetic: between 119680ae and the current branch tip, fftw moves 3.3.10 ->
# 3.3.11 and the host gcc@16.1.1 external stops being usable at all. So enforce
# the pin here, and fail loudly rather than silently building something else.
pin_package_recipes() {
    local repo_dir root have
    repo_dir="$(spack location --repo builtin 2>/dev/null || true)"
    [[ -n "$repo_dir" && -d "$repo_dir" ]] || die "could not locate the builtin package repo"

    root="$(git -C "$repo_dir" rev-parse --show-toplevel 2>/dev/null || true)"
    [[ -n "$root" ]] || die "builtin package repo at $repo_dir is not a git checkout"

    have="$(git -C "$root" rev-parse HEAD)"
    if [[ "$have" != "$SPACK_PACKAGES_COMMIT" ]]; then
        log "Pinning builtin recipes to ${SPACK_PACKAGES_COMMIT:0:8} (was ${have:0:8})"
        git -C "$root" fetch --depth 1 origin "$SPACK_PACKAGES_COMMIT" \
            || die "could not fetch recipe commit $SPACK_PACKAGES_COMMIT"
        git -C "$root" checkout --detach FETCH_HEAD >/dev/null 2>&1 \
            || die "could not check out recipe commit $SPACK_PACKAGES_COMMIT"
        have="$(git -C "$root" rev-parse HEAD)"
    fi
    [[ "$have" == "$SPACK_PACKAGES_COMMIT" ]] \
        || die "builtin recipes are at $have, expected $SPACK_PACKAGES_COMMIT"
    log "Builtin recipes pinned at ${have:0:8}"
}

# -----------------------------------------------------------------------------
# 1. Host discovery
# -----------------------------------------------------------------------------
load_modules() {
    [[ -n "$MODULES" ]] || return 0

    # Batch shells often have no module function; source the usual initscripts.
    if ! type module >/dev/null 2>&1; then
        for init in /usr/share/lmod/lmod/init/bash \
                    /etc/profile.d/modules.sh \
                    /etc/profile.d/lmod.sh \
                    /etc/profile.d/z00_lmod.sh; do
            if [[ -r "$init" ]]; then
                set +u; # shellcheck disable=SC1090
                source "$init"; set -u
                break
            fi
        done
    fi
    type module >/dev/null 2>&1 || die "--modules given but no module command is available"

    log "Loading modules: $MODULES"
    set +u
    # shellcheck disable=SC2086
    module load $MODULES
    set -u
}

# Registers every compiler Spack can find into the *site* scope of our private
# Spack. Site scope sits below the environment, so spack-env/spack.yaml and
# host/host.yaml always win over it.
detect_host_compilers() {
    log "Detecting host compilers"
    spack compiler find --scope site >/dev/null 2>&1 || true
    narrow_compiler_languages
    local found
    found="$(spack compiler list 2>/dev/null | grep -Ev '^(--|==>|$)' | tr '\n' ' ' || true)"
    [[ -n "${found// /}" ]] || die "no host compiler found; load one with --modules"
    log "Host compilers: $found"
}

# `spack compiler find` records every language it detects. On a host with gdc
# installed that means languages:='c,c++,d,fortran', and Spack then refuses to use
# the external at all ("Cannot satisfy 'gcc@16.1.1'"). Narrow it to what SW4 needs.
narrow_compiler_languages() {
    local script="$ENV_DIR/narrow-compiler-languages.py" n
    [[ -f "$script" ]] || die "missing helper: $script"
    n="$(spack python "$script" \
            "$SPACK_ROOT_DIR/etc/spack/packages.yaml" \
            "$SPACK_ROOT_DIR/etc/spack/site/packages.yaml" 2>/dev/null | tail -1)"
    [[ "${n:-0}" == "0" ]] || log "Narrowed $n compiler external(s) to c/c++/fortran"
}

detect_host_build_tools() {
    if (( NO_HOST_BUILD_TOOLS )); then
        log "Skipping host build-tool detection (--no-host-build-tools)"
        return
    fi
    log "Detecting host build tools (${#BUILD_TOOLS[@]} candidates)"
    spack external find --scope site "${BUILD_TOOLS[@]}" 2>&1 | sed 's/^/    /' || true
}

pick_compiler() {
    if [[ -z "$COMPILER_SPEC" ]]; then
        # Newest GCC that was detected. sort -V puts the highest version last.
        local newest
        newest="$(spack compiler list 2>/dev/null | tr ' ' '\n' \
                  | grep -E '^gcc@[0-9]' | sed 's/^gcc@//' | sort -V | tail -1 || true)"
        [[ -n "$newest" ]] || die "no host GCC found; name a compiler with --compiler, e.g. --compiler llvm@18"
        COMPILER_SPEC="gcc@$newest"
    fi
    COMPILER_PKG="${COMPILER_SPEC%%@*}"

    spack compiler list 2>/dev/null | tr ' ' '\n' | grep -qx -- "$COMPILER_SPEC" \
        || warn "$COMPILER_SPEC is not in \`spack compiler list\`; concretization will fail if it is not usable"

    log "Building everything with $COMPILER_SPEC"
}

# Sets MPI_SPEC / MPI_PREFIX / MPI_PKG. Uses each vendor's own query tool rather
# than `spack external find`, whose openmpi detection is broken in Spack 1.2
# (AttributeError: 'Spec' object has no attribute 'cc').
detect_mpi() {
    if [[ -n "$MPI_SPEC" ]]; then
        log "MPI given on the command line: $MPI_SPEC"
    elif [[ -n "${CRAY_MPICH_DIR:-}" ]]; then
        MPI_SPEC="cray-mpich@${CRAY_MPICH_VERSION:-8.1.0}"
        MPI_PREFIX="${MPI_PREFIX:-$CRAY_MPICH_DIR}"
    elif command -v ompi_info >/dev/null 2>&1; then
        local ver
        ver="$(ompi_info --parsable 2>/dev/null | sed -n 's/^ompi:version:full://p' | head -1)"
        MPI_SPEC="openmpi@${ver:-5.0.0}"
        MPI_PREFIX="${MPI_PREFIX:-$(mpicc --showme:prefix 2>/dev/null || true)}"
    elif command -v mpichversion >/dev/null 2>&1; then
        local ver
        ver="$(mpichversion --version 2>/dev/null | sed -n 's/.*Version:[[:space:]]*\([0-9.]*\).*/\1/p' | head -1)"
        MPI_SPEC="mpich@${ver:-4.0.0}"
        MPI_PREFIX="${MPI_PREFIX:-$(mpichversion --prefix 2>/dev/null || true)}"
    elif command -v mpiname >/dev/null 2>&1; then
        MPI_SPEC="mvapich2@$(mpiname -v 2>/dev/null | head -1)"
    elif mpiexec --version 2>&1 | grep -qi 'Intel(R) MPI'; then
        local ver
        ver="$(mpiexec --version 2>&1 | sed -n 's/.*Version \([0-9.]*\).*/\1/p' | head -1)"
        MPI_SPEC="intel-oneapi-mpi@${ver:-2021.0.0}"
    else
        die "could not detect an MPI installation. Load one with --modules, or name it with --mpi-spec/--mpi-prefix"
    fi

    # Fall back to deriving the prefix from the compiler wrapper's location.
    if [[ -z "$MPI_PREFIX" ]]; then
        local wrapper
        wrapper="$(command -v mpicc || command -v mpicxx || true)"
        [[ -n "$wrapper" ]] || die "found $MPI_SPEC but no mpicc/mpicxx on PATH; pass --mpi-prefix"
        MPI_PREFIX="$(cd -- "$(dirname -- "$(dirname -- "$(readlink -f "$wrapper")")")" && pwd)"
    fi
    [[ -d "$MPI_PREFIX" ]] || die "MPI prefix does not exist: $MPI_PREFIX"

    MPI_PKG="$(awk '{print $1}' <<<"$MPI_SPEC")"; MPI_PKG="${MPI_PKG%%@*}"

    # sw4 links Fortran objects and the recipe asks MPI for a Fortran wrapper, so
    # record whether the host MPI has one. Only openmpi and mpich declare a
    # `fortran` variant -- naming it on e.g. cray-mpich is a concretization error.
    if [[ "$MPI_PKG" == openmpi || "$MPI_PKG" == mpich ]]; then
        if command -v mpifort >/dev/null 2>&1 || command -v mpif90 >/dev/null 2>&1; then
            MPI_SPEC="$MPI_SPEC +fortran"
        else
            MPI_SPEC="$MPI_SPEC ~fortran"
            warn "host MPI has no Fortran wrapper; recording ~fortran"
        fi
    fi

    log "Host MPI: $MPI_SPEC at $MPI_PREFIX"
}

# Which sw4 source to build: by default the commit checked out right here, so
# "clone at revision X, run install.sh" builds exactly X.
resolve_sw4_version() {
    if [[ -n "$SW4_SPEC_VERSION" ]]; then
        log "Building declared package version sw4@$SW4_SPEC_VERSION"
        SW4_VERSION_CONSTRAINT="@$SW4_SPEC_VERSION"
        return
    fi
    local sha
    sha="$(git -C "$REPO_DIR" rev-parse HEAD 2>/dev/null || true)"
    [[ -n "$sha" ]] || die "not a git checkout; pass --sw4-version to build a declared version"

    if ! git -C "$REPO_DIR" branch -r --contains "$sha" 2>/dev/null | grep -q .; then
        local msg="HEAD ($sha) is not reachable from any remote branch. Spack fetches from the fork's remote, so the build will fail. Push the commit, or pass --sw4-version."
        (( ALLOW_UNPUSHED )) && warn "$msg" || die "$msg"
    fi
    if ! git -C "$REPO_DIR" diff --quiet HEAD -- . 2>/dev/null; then
        warn "working tree has uncommitted changes; they will NOT be in the build (Spack builds commit $sha from the remote)"
    fi
    log "Building sw4 at commit $sha"
    SW4_VERSION_CONSTRAINT="@git.$sha=$SW4_REF_VERSION"
}

write_host_config() {
    mkdir -p "$HOST_CFG_DIR"
    log "Writing $HOST_CFG_DIR/host.yaml"
    {
        echo "# Generated by install.sh on $(hostname) at $(date -Is)."
        echo "# Do not edit -- regenerated on every run, and gitignored."
        echo "toolchains:"
        echo "  # Referenced as %host_toolchain by spack-env/spack.yaml, which"
        echo "  # requires it for every node in the DAG."
        echo "  host_toolchain:"
        echo "  - spec: \"%c=$COMPILER_SPEC\""
        echo "    when: \"%c\""
        echo "  - spec: \"%cxx=$COMPILER_SPEC\""
        echo "    when: \"%cxx\""
        echo "  - spec: \"%fortran=$COMPILER_SPEC\""
        echo "    when: \"%fortran\""
        echo "packages:"
        echo "  all:"
        echo "    providers:"
        echo "      mpi: [$MPI_PKG]"
        echo "  # The compiler comes from the host; do not let Spack decide to build"
        echo "  # one instead. The externals themselves were registered into the site"
        echo "  # scope by \`spack compiler find\`."
        echo "  $COMPILER_PKG:"
        echo "    buildable: false"
        echo "  # The one host runtime dependency. Spack points OMPI_CC/OMPI_CXX/OMPI_FC"
        echo "  # (or the equivalent) at its own wrappers, so these host wrappers still"
        echo "  # drive the compiler Spack chose."
        echo "  $MPI_PKG:"
        echo "    buildable: false"
        echo "    externals:"
        echo "    - spec: \"$MPI_SPEC\""
        echo "      prefix: $MPI_PREFIX"
        echo "  sw4:"
        echo "    require:"
        echo "    - \"$SW4_VERSION_CONSTRAINT\""
        if [[ -n "$SW4_VARIANTS" ]]; then
            echo "    - \"$SW4_VARIANTS\""
        fi
        if [[ -n "$INSTALL_ROOT" || -n "$BUILD_STAGE" ]]; then
            echo "config:"
            [[ -n "$INSTALL_ROOT" ]] && { echo "  install_tree:"; echo "    root: $INSTALL_ROOT"; }
            [[ -n "$BUILD_STAGE" ]]  && { echo "  build_stage: [$BUILD_STAGE]"; }
        fi
    } >"$HOST_CFG_DIR/host.yaml"
}

# -----------------------------------------------------------------------------
# 2. Build
# -----------------------------------------------------------------------------
install_sw4() {
    log "Concretizing $ENV_DIR"
    spack -e "$ENV_DIR" concretize --force 2>&1 | tee "$HOST_CFG_DIR/concretize.log"

    if (( CONCRETIZE_ONLY )); then
        log "--concretize-only: stopping. Lockfile: $ENV_DIR/spack.lock"
        return
    fi

    log "Installing (this is the long part)"
    spack -e "$ENV_DIR" install --fail-fast -j "$JOBS"
}

summary() {
    local view="$ENV_DIR/.spack-env/view"
    log "Done."
    cat <<EOF

  sw4 binary   : $view/bin/sw4
  environment  : spack env activate -d $ENV_DIR
  spack prefix : $SPACK_ROOT_DIR
  lockfile     : $ENV_DIR/spack.lock   <-- commit this to pin the exact DAG

Reproduce this build elsewhere by committing spack.lock and running:
  ./install.sh            # regenerates host.yaml, then installs from the lock
EOF
}

main() {
    bootstrap_spack
    pin_package_recipes
    load_modules
    detect_host_compilers
    detect_host_build_tools
    pick_compiler
    detect_mpi
    resolve_sw4_version
    write_host_config
    install_sw4
    (( CONCRETIZE_ONLY )) || summary
}

main "$@"
