# Copyright Spack Project Developers. See COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)
"""SW4 built from the UC/QuakeCoRE (ucgmsim) fork rather than geodynamics/sw4.

This package shadows ``builtin.sw4`` and differs from it in two substantive ways.

*Source.*  ``git`` points at https://github.com/ucgmsim/sw4 and every version is
pinned to an explicit commit, so a given spec always fetches identical sources.
No branch versions are declared on purpose -- a branch would silently change
what "the same spec" builds.

*Build system.*  ``builtin.sw4`` is a ``MakefilePackage``; this uses the fork's
modernised ``CMakeLists.txt``.  That fixes several problems that are awkward to
work around from a Makefile build:

* ``build_type`` and ``precision`` are ordinary CMake cache entries
  (``CMAKE_BUILD_TYPE``, ``USE_DOUBLE``) rather than environment variables that
  select a build directory by string concatenation.
* The Makefile assigns ``FFTWHOME = $(SW4ROOT)`` with ``=``, which silently
  overrides an environment variable of the same name -- so ``builtin.sw4``'s
  ``+fftw`` actually points FFTW at PROJ's prefix.
* The Makefile picks up ``configs/make.<site>`` based on ``hostname``, so an
  unlucky hostname changes the build.

Two variants from the builtin recipe are gone because CMake handles them
differently: ``openmp`` (CMakeLists makes OpenMP unconditionally ``REQUIRED``)
and ``debug`` (use ``build_type=Debug``, which CMakePackage provides).

*Compiler switches.*  ``cmake/SW4Optimization.cmake`` picks architecture and
optimisation flags per deployment target and compiler.  ``sw4_target``,
``strict_fp`` and ``lto`` pass its options through.  ``cxxflags`` (e.g.
``cxxflags="-mprefer-vector-width=256"``) becomes ``SW4_EXTRA_RELEASE_FLAGS``,
so it comes after the target's flags and can override them, and the module's
fast-math guard sees it.

Spack's compiler wrapper puts the spec's own ``target`` flags (``-march`` for
the build host's microarchitecture, unless ``target=`` says otherwise) *before*
CMake's.  So with the default ``sw4_target=generic``, which adds no
architecture flags, SW4 is built for Spack's ``target``.  A named
``sw4_target`` overrides that ISA, because the compiler takes the last
``-march``.
"""

from spack_repo.builtin.build_systems.cmake import CMakePackage

from spack.package import *


class Sw4(CMakePackage):
    """SW4 implements substantial capabilities for 3-D seismic modeling.

    This is the ucgmsim fork, pinned to explicit commits.
    """

    homepage = "https://github.com/ucgmsim/sw4"
    git = "https://github.com/ucgmsim/sw4.git"

    maintainers("ucgmsim")

    license("GPL-2.0-or-later")

    # Pinned commits only. To build something else, add a line here rather than
    # pointing a version at a branch.
    #
    # 2026.07.30 == origin/single_precision @ 4d8c3f0f, the first commit whose
    # CMakeLists.txt is usable from Spack.
    version(
        "2026.07.30",
        commit="4d8c3f0f8c45bab6c51f993282716b15a12bb67e",
        preferred=True,
    )

    variant("proj", default=True, description="Build with PROJ")
    variant("hdf5", default=True, description="Build with HDF5")
    variant("zfp", default=False, when="+hdf5", description="Build with ZFP compression")
    variant("fftw", default=True, description="Build with FFTW (needed for randomized material)")
    variant(
        "precision",
        default="double",
        values=("double", "single"),
        multi=False,
        description="Floating point precision of float_sw4",
    )
    variant(
        "sw4_target",
        default="generic",
        values=(
            "generic",
            "auto",
            "native",
            "mn5-gpp",
            "cascade",
            "hpc3-genoa",
            "hpc3-milan",
            "hpc3-portable",
            "frontera",
            "stampede3",
            "stampede3-spr",
            "vista",
        ),
        multi=False,
        description=(
            "SW4_TARGET of cmake/SW4Optimization.cmake. generic adds no -march, "
            "leaving the ISA to Spack's target; auto guesses from site environment "
            "variables, which may pick a portable baseline below Spack's target"
        ),
    )
    variant(
        "strict_fp",
        default=False,
        description="SW4_STRICT_FP: no FMA contraction, for bit-reproducible results across machines",
    )
    variant("lto", default=False, description="SW4_LTO: link-time optimisation")
    variant(
        "pytests",
        default=False,
        description="Install the pytest regression suite and the Python needed to run it",
    )

    depends_on("c", type="build")
    depends_on("cxx", type="build")
    depends_on("fortran", type="build")
    depends_on("cmake@3.16:", type="build")

    depends_on("mpi")
    depends_on("blas")
    depends_on("lapack")
    depends_on("proj@9:", when="+proj")
    depends_on("hdf5@1.14: +mpi", when="+hdf5")
    depends_on("fftw@3: +mpi", when="+fftw")
    depends_on("zfp", when="+zfp")
    depends_on("h5z-zfp@1.1.0:", when="+zfp")

    # The Check_Result_* ctest targets shell out to pytest/check_results.py.
    depends_on("python", type=("build", "run"), when="+pytests")
    depends_on("py-h5py", type=("build", "run"), when="+pytests+hdf5")

    def flag_handler(self, name, flags):
        # cxxflags go to CMake as SW4_EXTRA_RELEASE_FLAGS (see cmake_args)
        # rather than through the wrapper. The wrapper puts them before CMake's
        # target flags, so e.g. -mprefer-vector-width=256 would lose to the
        # target's 512.
        if name == "cxxflags":
            return (None, None, None)
        return (flags, None, None)

    def cmake_args(self):
        spec = self.spec

        args = [
            self.define_from_variant("USE_PROJ", "proj"),
            self.define_from_variant("USE_HDF5", "hdf5"),
            self.define_from_variant("USE_ZFP", "zfp"),
            self.define_from_variant("USE_FFTW3", "fftw"),
            self.define("USE_SZ", False),
            self.define("USE_DOUBLE", spec.satisfies("precision=double")),
            # Empty, so SW4Optimization.cmake chooses flags for SW4_TARGET.
            self.define("SW4_ARCH_FLAGS", ""),
            self.define_from_variant("SW4_TARGET", "sw4_target"),
            self.define_from_variant("SW4_STRICT_FP", "strict_fp"),
            self.define_from_variant("SW4_LTO", "lto"),
            # See flag_handler.
            self.define(
                "SW4_EXTRA_RELEASE_FLAGS", " ".join(spec.compiler_flags["cxxflags"])
            ),
        ]

        # Point FindMPI at the wrappers Spack chose rather than whatever is first
        # on PATH.
        args += [
            self.define("MPI_C_COMPILER", spec["mpi"].mpicc),
            self.define("MPI_CXX_COMPILER", spec["mpi"].mpicxx),
            self.define("MPI_Fortran_COMPILER", spec["mpi"].mpifc),
        ]

        # CMake's FindBLAS/FindLAPACK guess among many vendors; name the one we
        # actually depend on. spack.yaml pins the provider to OpenBLAS.
        if spec["blas"].name == "openblas":
            args.append(self.define("BLA_VENDOR", "OpenBLAS"))
        else:
            args += [
                self.define("BLAS_LIBRARIES", spec["blas"].libs.joined(";")),
                self.define("LAPACK_LIBRARIES", spec["lapack"].libs.joined(";")),
            ]

        # <Package>_ROOT is honoured by find_path/find_library/find_package under
        # CMP0074, which cmake_minimum_required(3.16) enables.
        if spec.satisfies("+proj"):
            args.append(self.define("PROJ_ROOT", spec["proj"].prefix))
        if spec.satisfies("+hdf5"):
            args.append(self.define("HDF5_ROOT", spec["hdf5"].prefix))
            args.append(self.define("HDF5_PREFER_PARALLEL", True))
        if spec.satisfies("+fftw"):
            args.append(self.define("FFTW3_ROOT", spec["fftw"].prefix))
        if spec.satisfies("+zfp"):
            args.append(self.define("ZFP_ROOT", spec["zfp"].prefix))
            args.append(self.define("H5Z_ZFP_ROOT", spec["h5z-zfp"].prefix))

        # Only build the ctest suite's Python half when we have Python for it.
        if not spec.satisfies("+pytests"):
            args.append(self.define("BUILD_TESTING", False))

        return args

    @run_after("install")
    def install_pytest_suite(self):
        """ctest drives pytest/check_results.py against pytest/reference/."""
        if self.spec.satisfies("+pytests"):
            install_tree("pytest", join_path(self.prefix.share, "sw4", "pytest"))
