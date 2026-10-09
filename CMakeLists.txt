cmake_minimum_required(VERSION 3.16)

# Prefer the MPI compiler wrappers as the project's compilers, unless the user
# already picked one explicitly (via -DCMAKE_<LANG>_COMPILER=... or the CC/CXX/FC
# environment variables). This must happen before project(), since a language's
# compiler can only be chosen before it is first detected.
#
# Why: find_package(MPI) below still works with a plain (non-wrapper) compiler,
# but then has to derive the right MPI include/library flags itself, usually by
# scraping the output of `mpicc -show` and friends. That heuristic is fragile
# across MPI implementations/versions -- e.g. on mahuika with OpenMPI 5.0.10 +
# GCC 15.2.0, it fails outright ("Could NOT find MPI_C (missing: MPI_C_WORKS)")
# even though mpicc/mpicxx/mpifort are on PATH and work fine invoked directly.
# Using the wrapper as CMAKE_<LANG>_COMPILER sidesteps that derivation entirely:
# the wrapper already knows its own flags, and FindMPI recognizes an MPI-aware
# compiler immediately.
if(NOT CMAKE_C_COMPILER AND NOT DEFINED ENV{CC})
  find_program(SW4_MPICC_EXECUTABLE NAMES mpicc)
  if(SW4_MPICC_EXECUTABLE)
    set(CMAKE_C_COMPILER "${SW4_MPICC_EXECUTABLE}")
  endif()
endif()
if(NOT CMAKE_CXX_COMPILER AND NOT DEFINED ENV{CXX})
  find_program(SW4_MPICXX_EXECUTABLE NAMES mpicxx mpic++)
  if(SW4_MPICXX_EXECUTABLE)
    set(CMAKE_CXX_COMPILER "${SW4_MPICXX_EXECUTABLE}")
  endif()
endif()
# There does not appear be any fortran code calling MPI so turning off
#if(NOT CMAKE_Fortran_COMPILER AND NOT DEFINED ENV{FC})
#  find_program(SW4_MPIFORT_EXECUTABLE NAMES mpifort mpif90)
#  if(SW4_MPIFORT_EXECUTABLE)
#    set(CMAKE_Fortran_COMPILER "${SW4_MPIFORT_EXECUTABLE}")
#  endif()
#endif()

project(SW4 C CXX Fortran)

# The sources use C++17 features (e.g. `if constexpr` in rhs4th3fortc.C). Without
# an explicit standard each compiler uses its own default, and those older than
# C++17 (e.g. Apple clang) emit a -Wc++17-extensions warning for every use.
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

list(PREPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/configs/")

include(CTest)

macro(today RESULT)
  if(WIN32)
    execute_process(COMMAND "date" "/T" OUTPUT_VARIABLE ${RESULT})
    string(REGEX REPLACE "(..)/(..)/..(..).*" "\\3\\2\\1" ${RESULT} ${${RESULT}})
  elseif(UNIX)
    execute_process(COMMAND "date" "+%d/%m/%Y" OUTPUT_VARIABLE ${RESULT})
    string(REGEX REPLACE "(..)/(..)/..(..).*" "\\3\\2\\1" ${RESULT} ${${RESULT}})
  else()
    message(SEND_ERROR "date not implemented")
    set(${RESULT} 000000)
  endif()
endmacro()

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING
      "Build type (Debug, Release, RelWithDebInfo, MinSizeRel)" FORCE)
endif()

set(CMAKE_Fortran_MODULE_DIRECTORY "${PROJECT_BINARY_DIR}/include")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/bin")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/lib")

site_name(HOSTNAME)
set(USERNAME $ENV{USER})
today(CURDATETIME)

find_package(MPI REQUIRED COMPONENTS C CXX)

if(NOT PROJ_OSGEO4W_HOME)
  set(PROJ_OSGEO4W_HOME ${PROJ_HOME})
endif()

option(USE_PROJ "Use PROJ." OFF)
if(USE_PROJ)
  find_package(PROJ REQUIRED)
  if(PROJ_VERSION VERSION_LESS "6.0.0")
    message(SEND_ERROR "PROJ library must be version 6 or higher")
  endif()
endif()

option(USE_HDF5 "Use HDF5." OFF)
if(USE_HDF5)
  find_package(ZLIB REQUIRED)
  find_package(HDF5 REQUIRED)
endif()

option(USE_ZFP "Use ZFP compression." OFF)
if(USE_ZFP)
  set(H5Z_ZFP_USE_STATIC_LIBS ON)
  find_package(ZFP REQUIRED CONFIG)
  find_package(H5Z_ZFP 1.1.0 REQUIRED CONFIG)
endif()

option(USE_SZ "Use SZ compression." OFF)
if(USE_SZ)
  find_package(SZ REQUIRED)
endif()

option(USE_FFTW3 "Use FFTW3." OFF)
if(USE_FFTW3)
  find_package(FFTW3 REQUIRED)
  if(NOT FFTW3_MPI_LIBRARIES)
    message(SEND_ERROR "FFTW3 MPI library (fftw3_mpi) is required when USE_FFTW3=ON")
  endif()
endif()

find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)
find_package(OpenMP REQUIRED)
message(STATUS "OpenMP_CXX_FLAGS = ${OpenMP_CXX_FLAGS}")

include(FortranCInterface)
FortranCInterface_HEADER(FC.h MACRO_NAMESPACE "FC_")

option(USE_DOUBLE "Use double precision" ON)

set(SW4_ARCH_FLAGS "" CACHE STRING
    "Architecture-specific compiler flags. Leave empty to let SW4_TARGET pick them.")
set(SW4_EXTRA_RELEASE_FLAGS "" CACHE STRING
    "Extra Release-mode flags, appended after the ones SW4_TARGET selects.")

# Selects numerically-inert Release flags for the deployment target and the
# compiler in use. Defines SW4_RESOLVED_ARCH_FLAGS / _CXX_OPT / _FORT_OPT.
list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")
include(SW4Optimization)

file(GLOB SW4_QUADPACK_SRC_FILES CONFIGURE_DEPENDS src/quadpack/*.f)

set(SW4_SRC_C_FILES src/EW.C src/Sarray.C src/version.C src/parseInputFile.C
    src/ForcingTwilight.C src/curvilinearGrid.C src/parallelStuff.C src/Source.C
    src/MaterialProperty.C src/MaterialData.C src/material.C src/setupRun.C
    src/solve.C src/Parallel_IO.C src/Image.C src/GridPointSource.C src/MaterialBlock.C
    src/TimeSeries.C src/sacsubc.C src/SuperGrid.C src/TestRayleighWave.C src/MaterialPfile.C
    src/Filter.C src/Polynomial.C src/SecondOrderSection.C src/time_functions.C src/Qspline.C
    src/MaterialIfile.C src/GeographicProjection.C src/Image3D.C src/ESSI3D.C src/ESSI3DHDF5.C
    src/MaterialVolimagefile.C src/MaterialRfile.C src/MaterialSfile.C src/AnisotropicMaterialBlock.C
    src/sacutils.C src/DataPatches.C src/addmemvarforcing2.C src/consintp.C src/oddIoddJinterp.C
    src/evenIoddJinterp.C src/oddIevenJinterp.C src/evenIevenJinterp.C src/CheckPoint.C src/geodyn.C
    src/AllDims.C src/Patch.C src/RandomizedMaterial.C src/MaterialInvtest.C src/sw4-prof.C
    src/sachdf5.C src/readhdf5.C src/TestTwilight.C src/TestPointSource.C src/curvilinear4sgwind.C
    src/TestEcons.C src/GridGenerator.C src/GridGeneratorGeneral.C src/GridGeneratorGaussianHill.C
    src/CurvilinearInterface2.C src/SfileOutput.C src/pseudohess.C src/fastmarching.C src/solveTT.C
    src/rhs4th3point.C src/MaterialGMG.C src/addsgdc.C src/bcfortc.C src/bcfortanisgc.C 
    src/bcfreesurfcurvanic.C src/boundaryOpc.C src/energy4c.C src/checkanisomtrlc.C src/computedtanisoc.C 
    src/curvilinear4sgc.C src/gradientsc.C src/randomfield3dc.C src/innerloop-ani-sgstr-vcc.C 
    src/ilanisocurvc.C src/rhs4curvilinearc.C src/rhs4curvilinearsgc.C src/rhs4th3fortc.C src/solerr3c.C 
    src/testsrcc.C src/rhs4th3windc.C src/tw_aniso_forcec.C src/tw_aniso_force_ttc.C src/velsumc.C 
    src/twilightfortc.C src/twilightsgfortc.C src/tw_ani_stiffc.C src/anisomtrltocurvilinearc.C 
    src/scalar_prodc.C src/updatememvarc.C src/addsg4windc.C src/bndryOpNoGhostc.C src/rhs4th3windc2.C
)

set(SW4_SRC_F_FILES src/rayleighfort.f src/lamb_exact_numquad.f)

# Sources specific to the material optimization driver (sw4mopt)
set(SW4_SRC_MOPT_FILES src/moptmain.C src/solve-backward-allpars.C src/lbfgs.C src/nlcg.C
    src/ProjectMtrl.C src/MaterialParameterization.C src/Mopt.C src/MaterialParCartesian.C
    src/InterpolateMaterial.C src/MaterialParCartesianVels.C src/MaterialParCartesianVp.C
    src/MParGridFile.C src/MaterialParCartesianVsVp.C src/MaterialParAllpts.C
    src/MaterialParCart.C src/solve-dudp.C src/MaterialParCurv.C
)

# Interface target carrying the compile/link settings shared by every SW4 target
add_library(sw4_config INTERFACE)

target_compile_definitions(sw4_config INTERFACE
    EW_MADEBY="${USERNAME}"
    EW_OPT_LEVEL="${CMAKE_BUILD_TYPE}"
    EW_COMPILER="${CMAKE_CXX_COMPILER}"
    EW_LIBDIR="NA"
    EW_INCDIR="NA"
    EW_HOSTNAME="${HOSTNAME}"
    EW_WHEN="${CURDATETIME}"
    CMAKE_FC_FUNC=1
)

if(USE_PROJ AND PROJ_FOUND)
  target_compile_definitions(sw4_config INTERFACE ENABLE_PROJ)
endif()

if(USE_HDF5 AND HDF5_FOUND)
  target_compile_definitions(sw4_config INTERFACE USE_HDF5)
endif()

if(USE_ZFP)
  target_compile_definitions(sw4_config INTERFACE USE_ZFP)
endif()

if(USE_SZ AND SZ_FOUND)
  target_compile_definitions(sw4_config INTERFACE USE_SZ)
endif()

if(USE_FFTW3 AND FFTW3_FOUND)
  target_compile_definitions(sw4_config INTERFACE USE_FFTW)
endif()

target_include_directories(sw4_config INTERFACE
    "${PROJECT_BINARY_DIR}"
    "${PROJECT_SOURCE_DIR}/src"
)

if(USE_DOUBLE)
  target_include_directories(sw4_config INTERFACE "${PROJECT_SOURCE_DIR}/src/double")
  message(STATUS "Precision: double (USE_DOUBLE=ON)")
else()
  target_include_directories(sw4_config INTERFACE "${PROJECT_SOURCE_DIR}/src/float")
  message(STATUS "Precision: single/float (USE_DOUBLE=OFF)")
endif()

if(USE_SZ AND SZ_FOUND)
  target_include_directories(sw4_config INTERFACE ${SZ_INCLUDE_DIR})
endif()

string(JOIN " " SW4_RELEASE_CXX_FLAGS
    ${SW4_RESOLVED_ARCH_FLAGS} ${SW4_RESOLVED_CXX_OPT} ${SW4_EXTRA_RELEASE_FLAGS})
string(JOIN " " SW4_RELEASE_FORT_FLAGS
    ${SW4_RESOLVED_ARCH_FLAGS} ${SW4_RESOLVED_FORT_OPT})

target_compile_options(sw4_config INTERFACE
    "$<$<AND:$<COMPILE_LANGUAGE:CXX>,$<CONFIG:Debug>>:SHELL:-O0>"
    "$<$<AND:$<COMPILE_LANGUAGE:CXX>,$<CONFIG:Release>>:SHELL:${SW4_RELEASE_CXX_FLAGS}>"
    "$<$<AND:$<COMPILE_LANGUAGE:CXX>,$<CONFIG:RelWithDebInfo>>:SHELL:-O3 ${SW4_RELEASE_CXX_FLAGS}>"
    "$<$<AND:$<COMPILE_LANGUAGE:Fortran>,$<CONFIG:Release>>:SHELL:${SW4_RELEASE_FORT_FLAGS}>"
    "$<$<AND:$<COMPILE_LANGUAGE:Fortran>,$<CONFIG:RelWithDebInfo>>:SHELL:${SW4_RELEASE_FORT_FLAGS}>"
)

target_link_libraries(sw4_config INTERFACE
    MPI::MPI_C
    MPI::MPI_CXX
    BLAS::BLAS
    LAPACK::LAPACK
    OpenMP::OpenMP_CXX
)

if(USE_PROJ AND PROJ_FOUND)
  target_link_libraries(sw4_config INTERFACE PROJ::proj)
endif()

if(USE_HDF5 AND HDF5_FOUND)
  target_link_libraries(sw4_config INTERFACE HDF5::HDF5)
endif()

if(USE_ZFP)
  target_link_libraries(sw4_config INTERFACE h5z_zfp::h5z_zfp)
endif()

if(USE_FFTW3 AND FFTW3_FOUND)
  target_link_libraries(sw4_config INTERFACE FFTW3::FFTW3 FFTW3::FFTW3_MPI)
endif()

# Objects shared by sw4 and sw4mopt; compiled once, linked into both drivers
add_library(sw4_common OBJECT ${SW4_SRC_C_FILES} ${SW4_SRC_F_FILES} ${SW4_QUADPACK_SRC_FILES})
target_link_libraries(sw4_common PUBLIC sw4_config)

add_executable(sw4 src/main.C)
target_link_libraries(sw4 PRIVATE sw4_common)

add_custom_command(TARGET sw4 POST_BUILD COMMAND cat ${CMAKE_SOURCE_DIR}/wave.txt)

# The sw4mopt sources have not been ported to single precision: they still use
# raw double* where the solver now expects float_sw4*, so they only compile when
# USE_DOUBLE=ON. Default the target on accordingly.
option(BUILD_SW4MOPT "Build the sw4mopt material optimization driver (requires USE_DOUBLE)." ${USE_DOUBLE})
if(BUILD_SW4MOPT AND NOT USE_DOUBLE)
  message(SEND_ERROR "BUILD_SW4MOPT requires USE_DOUBLE=ON: the sw4mopt sources are not single-precision clean.")
endif()

if(BUILD_SW4MOPT)
  add_executable(sw4mopt ${SW4_SRC_MOPT_FILES})
  target_link_libraries(sw4mopt PRIVATE sw4_common)
endif()

include(GNUInstallDirs)
install(TARGETS sw4 RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
if(BUILD_SW4MOPT)
  install(TARGETS sw4mopt RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
endif()

# Tests: cases live in tests/cases/*.toml and are registered by the harness in
# tests/harness (see tests/README.md). Select with labels: `ctest -L pr` runs
# the pull-request tier, plain `ctest` runs everything (the nightly tier).
if(BUILD_TESTING)
  add_subdirectory(tests)
endif()
