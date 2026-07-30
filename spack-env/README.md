# Reproducible SW4 builds with Spack

`../install.sh` is the entry point. It bootstraps a private Spack, discovers the
host compiler and MPI, and builds SW4 from **this fork** plus every other
dependency from source.

```sh
./install.sh                                        # everything
./install.sh --concretize-only                      # just resolve, write spack.lock
./install.sh --variants "precision=single"
./install.sh --modules "gcc/13.2.0 openmpi/4.1.6"   # cluster with modules
./install.sh --compiler gcc@13.2.0 --build-stage /local/scratch
```

Verified on this host: `sw4` builds, links PROJ/HDF5/FFTW/OpenBLAS/libgfortran
from Spack and `libmpi.so.40` from `/usr`, and reproduces
`pytest/reference/twilight/flat-twi-1` bit-identically under `mpirun -n 2`.

## Layout

| Path                           | Purpose                                                                                                                  |
|--------------------------------|--------------------------------------------------------------------------------------------------------------------------|
| `spack.yaml`                   | The environment. Portable — no host paths.                                                                               |
| `repo/spack_repo/sw4_ucgmsim/` | Custom Spack repo; shadows `builtin.sw4` with a recipe pointing at this fork at pinned commits, built via CMake.         |
| `narrow-compiler-languages.py` | Post-processes detected compiler externals (see below).                                                                  |
| `host/host.yaml`               | **Generated** by `install.sh`: the `host_toolchain` definition, the MPI external, and the pinned sw4 commit. Gitignored. |
| `spack.lock`                   | The concrete DAG.                                                                                                        |


## What is pinned

- **Spack** — `SPACK_REF` in `install.sh` (a tag, not a branch).
- **Package recipes** — `SPACK_PACKAGES_COMMIT`, mirrored by
  `repos: builtin: commit:` in `spack.yaml`. Keep the two in sync.
- **SW4 sources** — `install.sh` builds the
  commit checked out in the working copy, and refuses to proceed if it is not on
  a remote, because Spack fetches from the remote and would otherwise fail late.
- **The whole DAG** — `spack.lock`.

### The recipe pin needs enforcing, not just declaring

`repos: builtin: commit:` records the intent but Spack does not honour it on the
*initial* clone — it lands on the tip of `releases/v2026.06`. That is not
cosmetic: between the pinned commit and the current branch tip, `fftw` moves
3.3.10 → 3.3.11 and `compiler-wrapper` 1.0 → 1.1.0. So `install.sh` checks the
clone out itself and fails loudly if it cannot.

## Host components

Three things come from the host, all deliberate:

1. **The compiler.** Auto-detected (newest GCC by default; `--compiler` to
   override). `spack.yaml` then requires it for *every* node via
   `packages: all: require: "%host_toolchain"`.
2. **MPI.** Detected via the vendor's own query tools — `ompi_info`,
   `mpichversion`, `mpiname`, `$CRAY_MPICH_DIR` — because `spack external find`
   crashes on openmpi in Spack 1.2 (`'Spec' object has no attribute 'cc'`).
   Spack points `OMPI_CC`/`OMPI_CXX`/`OMPI_FC` at its own wrappers, so the host
   `mpicxx` still drives the compiler Spack chose. On a module-based cluster the
   providing module is recorded on the external — see below.
3. **Build-only tools** — `m4`, `tar`, `perl`, `autoconf`, … (`BUILD_TOOLS` in
   `install.sh`). None contribute code to the binary. This is not just for
   speed: Spack's pinned `m4@1.4.20` and `tar@1.35` bundle a gnulib that does
   not compile against glibc ≥ 2.42 (it collides with glibc's `_Generic`
   `bsearch`), and the pinned recipes offer nothing newer — host m4 1.4.21+ has
   the fix. `--no-host-build-tools` disables this and will fail on such a host.

Anything that can be *linked* — xz, bzip2, zlib, zstd, ncurses, gettext,
libxml2, openssl, sqlite, curl — is excluded from that list and built by Spack.

### Module-provided MPI must be recorded as a module, not just a prefix

Spack's `clean_environment()` unsets `LD_LIBRARY_PATH`, `LIBRARY_PATH` and
`CPATH` before every build. Spack does pass `-L<mpi prefix>/lib` for the
external, so a system MPI under `/usr` works — but on an EasyBuild/Lmod cluster
that is not enough. `libmpi.so` carries `DT_NEEDED` entries for PMIx, hwloc,
libevent and UCX, which live in *sibling* module prefixes, so `ld` cannot
resolve them and the build fails at, for example, fftw's configure:

```
checking for mpicc... .../OpenMPI/5.0.8-GCC-14.3.0/bin/mpicc
checking for MPI_Init in -lmpi... no
configure: error: could not find mpi library for --enable-mpi
```

The fix is Spack's own `load_external_modules()`, which runs *after* the scrub.
`install.sh` therefore records the providing module on the external:

```yaml
  openmpi:
    externals:
    - spec: "openmpi@5.0.8 +fortran"
      prefix: /shared/.../OpenMPI/5.0.8-GCC-14.3.0
      modules: [OpenMPI/5.0.8-GCC-14.3.0]
```

The module name is matched out of `$LOADEDMODULES` (so it works whether you use
`--modules` or load them yourself in the job script), handling EasyBuild-style
`OpenMPI/5.0.8-GCC-14.3.0`, flat names like `openmpi-4.1.5-gcc11`, `mpich/…`,
`impi/…` and `cray-mpich/…`. Override with `--mpi-module`. If MPI is found at a
non-system prefix and no module can be identified, `install.sh` warns, because
that is exactly the combination that fails.

The compiler external needs no such treatment — Spack invokes it by absolute
path, and `gcc-runtime` copies `libstdc++`/`libgfortran` into the Spack store.

### One job at a time per checkout

`host/host.yaml`, `spack.lock` and `.spack-env/` live in the checkout, so two
concurrent jobs pointed at the same clone will fight over them even with
different `--prefix`. Use one job per checkout, or separate clones.

### Why `narrow-compiler-languages.py` exists

`spack compiler find` records every language it detects. On a host with `gdc`
installed it writes `languages:='c,c++,d,fortran'`, and Spack then cannot use
that external at all — concretizing against `%c=gcc@16.1.1` fails with a bare
`Cannot satisfy 'gcc@16.1.1'`. The identical entry with
`languages:='c,c++,fortran'` concretizes in ~20s. The helper narrows detected
externals to the three languages SW4 needs.

### A note on the compiler being external

`packages: all: require: "%host_toolchain"` works *because* the compiler is
external: an external has no dependencies, so there is no cycle. Pointing the
same line at a Spack-**built** GCC is unsatisfiable — the requirement also
applies to GCC's own `gmp`/`mpfr`/`mpc`, giving `gcc → gmp → gcc`. Spack 1.2
additionally refuses any compiler that is neither external nor already
installed (`Only external, or concrete, compilers are allowed for the cxx
language`), so a Spack-built compiler would need a separate install phase before
this environment could even be concretized.

## Fork-specific recipe notes

`repo/spack_repo/sw4_ucgmsim/packages/sw4/package.py` uses the fork's
`CMakeLists.txt`, unlike `builtin.sw4` which is a `MakefilePackage`. Variants:

| Variant                    | Default   | Notes                                                   |
|----------------------------|-----------|---------------------------------------------------------|
| `precision=double\|single` | `double`  | `-DUSE_DOUBLE`                                          |
| `proj` / `hdf5` / `fftw`   | on        | `+fftw` is needed for randomized material               |
| `zfp`                      | off       | requires `+hdf5`                                        |
| `native`                   | off       | `-march=native -mtune=native`; not relocatable          |
| `pytests`                  | off       | installs `pytest/` and the Python the ctest checks need |
| `build_type`               | `Release` | from `CMakePackage`; use `build_type=Debug`             |

`openmp` is not a variant because `CMakeLists.txt` makes OpenMP unconditionally
`REQUIRED`.
