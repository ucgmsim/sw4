cmake_minimum_required(VERSION 3.16)
project(SW4 C CXX Fortran)

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

find_package(MPI REQUIRED)

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
    "Architecture-specific compiler flags (e.g. -march=native -mtune=native)")
set(SW4_EXTRA_RELEASE_FLAGS "-funroll-loops" CACHE STRING
    "Additional Release-mode optimization flags")

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

target_compile_options(sw4_config INTERFACE
    "$<$<AND:$<COMPILE_LANGUAGE:CXX>,$<CONFIG:Debug>>:SHELL:-O0>"
    "$<$<AND:$<COMPILE_LANGUAGE:CXX>,$<CONFIG:Release>>:SHELL:${SW4_ARCH_FLAGS} ${SW4_EXTRA_RELEASE_FLAGS}>"
    "$<$<AND:$<COMPILE_LANGUAGE:CXX>,$<CONFIG:RelWithDebInfo>>:SHELL:-O3 ${SW4_ARCH_FLAGS}>"
    "$<$<AND:$<COMPILE_LANGUAGE:Fortran>,$<CONFIG:Release>>:SHELL:${SW4_ARCH_FLAGS}>"
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

set(TESTING_LEVEL 0 CACHE STRING "Level of testing thoroughness: 0=short tests (<1 minute each), 2=all tests (may take an hour)")
set(MPI_NUM_TEST_PROCS 4 CACHE STRING "Number of MPI processes to use in tests")

set(PYTEST_DIR ${CMAKE_SOURCE_DIR}/pytest)
set(REF_DIR ${CMAKE_SOURCE_DIR}/pytest/reference)

set(TEST_DIRS       twilight               twilight         twilight
                    twilight               twilight         twilight
                    attenuation            attenuation      attenuation
                    attenuation            attenuation      attenuation
                    meshrefine             meshrefine       meshrefine
                    meshrefine             meshrefine       meshrefine
                    lamb                   lamb             lamb
                    pointsource            pointsource      pointsource)

set(RESULT_SUBDIRS  flat-twi-1             flat-twi-2       flat-twi-3
                    gauss-twi-1            gauss-twi-2      gauss-twi-3
                    tw-att-1               tw-att-2         tw-att-3
                    tw-topo-att-1          tw-topo-att-2    tw-topo-att-3
                    refine-el-1            refine-att-1     refine-att-2nd-1
                    refine-el-2            refine-att-2     refine-att-2nd-2
                    lamb-1                 lamb-2           lamb-3
                    pointsource-sg-1       pointsource-sg-2 pointsource-sg-3)

set(TEST_IN_FILES   flat-twi-1.in          flat-twi-2.in       flat-twi-3.in
                    gauss-twi-1.in         gauss-twi-2.in      gauss-twi-3.in
                    tw-att-1.in            tw-att-2.in         tw-att-3.in
                    tw-topo-att-1.in       tw-topo-att-2.in    tw-topo-att-3.in
                    refine-el-1.in         refine-att-1.in     refine-att-2nd-1.in
                    refine-el-2.in         refine-att-2.in     refine-att-2nd-2.in
                    lamb-1.in              lamb-2.in           lamb-3.in
                    pointsource-sg-1.in    pointsource-sg-2.in pointsource-sg-3.in)

set(TEST_BASE_FILES TwilightErr            TwilightErr             TwilightErr    
                    TwilightErr            TwilightErr             TwilightErr    
                    TwilightErr            TwilightErr             TwilightErr    
                    TwilightErr            TwilightErr             TwilightErr    
                    TwilightErr            TwilightErr             TwilightErr    
                    TwilightErr            TwilightErr             TwilightErr
                    LambErr                LambErr                 LambErr
                    PointSourceErr         PointSourceErr          PointSourceErr)

set(TEST_CHECKS     compare             compare                 compare 
                    compare             compare                 compare 
                    compare             compare                 compare 
                    compare             compare                 compare 
                    compare             compare                 compare 
                    compare             compare                 compare
                    compare             compare                 compare
                    compare             compare                 compare)

set(TEST_LEVELS     0                       0                    1
                    0                       0                    1
                    0                       0                    1
                    0                       1                    2
                    0                       0                    0
                    1                       1                    1
                    0                       1                    2
                    0                       1                    2)

set(TEST_ERRINF_TOL 1e-5                    1e-5                 1e-5
                    1e-5                    1e-5                 1e-5
                    1e-5                    1e-5                 1e-5
                    1e-5                    1e-5                 1e-5
                    1e-5                    1e-5                 1e-5
                    1e-5                    1e-5                 1e-5
                    1e-5                    1e-5                 1e-5
                    1e-5                    1e-5                 1e-5)

set(TEST_ERRL2_TOL  1e-5                    1e-5                 1e-5
                    1e-5                    1e-5                 1e-5
                    1e-5                    1e-5                 1e-5
                    1e-5                    1e-5                 1e-5
                    1e-5                    1e-5                 1e-5
                    1e-5                    1e-5                 1e-5
                    1e-5                    1e-5                 1e-5
                    1e-5                    1e-5                 1e-5)

set(TEST_SOLINF_TOL 1e-2                    1e-2                 1e-2
                    1e-2                    1e-2                 1e-2
                    1e-2                    1e-2                 1e-2
                    1e-2                    1e-2                 1e-2
                    1e-2                    1e-2                 1e-2
                    1e-2                    1e-2                 1e-2
                    0                       0                    0
                    0                       0                    0)

list(LENGTH TEST_DIRS N)
math(EXPR NUM_TESTS "${N}-1")

foreach(TEST_IND RANGE ${NUM_TESTS})
  list(GET TEST_LEVELS ${TEST_IND} TEST_LEVEL)
  if(NOT ${TESTING_LEVEL} LESS ${TEST_LEVEL})
    list(GET TEST_DIRS ${TEST_IND} TEST_DIR)
    list(GET RESULT_SUBDIRS ${TEST_IND} RESULT_SUBDIR)
    set(Test_Name ${TEST_DIR}/${RESULT_SUBDIR})
    list(GET TEST_IN_FILES ${TEST_IND} TEST_IN_FILE)
    list(GET TEST_BASE_FILES ${TEST_IND} TEST_BASE_FILE)
    list(GET TEST_CHECKS ${TEST_IND} TEST_CHECK_TYPE)
    list(GET TEST_ERRINF_TOL ${TEST_IND} ERRINF_TOL)
    list(GET TEST_ERRL2_TOL ${TEST_IND} ERRL2_TOL)
    list(GET TEST_SOLINF_TOL ${TEST_IND} SOLINF_TOL)
    add_test(
        NAME Run_${Test_Name}
        WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
        COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${MPI_NUM_TEST_PROCS} ${MPIEXEC_PREFLAGS} ${MPIEXEC_POSTFLAGS}
            ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/sw4 ${REF_DIR}/${TEST_DIR}/${TEST_IN_FILE}
    )

    set(TEST_REF_FILE ${REF_DIR}/${TEST_DIR}/${RESULT_SUBDIR}/${TEST_BASE_FILE}.txt)
    set(TEST_OUT_FILE ${PROJECT_BINARY_DIR}/${RESULT_SUBDIR}/${TEST_BASE_FILE}.txt)
    add_test(
        NAME Check_Result_${Test_Name}
        WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
        COMMAND ${PYTEST_DIR}/check_results.py ${TEST_CHECK_TYPE} ${TEST_REF_FILE} ${TEST_OUT_FILE} ${ERRINF_TOL} ${ERRL2_TOL} ${SOLINF_TOL}
    )
    set_tests_properties(
        Check_Result_${Test_Name}
        PROPERTIES DEPENDS ${Test_Name}
    )
  endif()
endforeach()
