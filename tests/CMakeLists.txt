# SW4 test suite. Every case is a TOML file under cases/; the harness turns each
# run into a `<case>:run:<name>` test (a ctest fixture) and each assertion into a
# `<case>:check:<name>` test that requires those fixtures.

find_program(UV_EXECUTABLE uv)
if(NOT UV_EXECUTABLE)
  message(WARNING "uv not found: SW4 tests are not registered. Install it from https://docs.astral.sh/uv/")
  return()
endif()

set(SW4TEST_PEER_OUT_DIR "" CACHE PATH
    "tests/out directory of a build with the other precision; enables float-vs-double cross checks")
set(SW4TEST_MPI_PREFLAGS "" CACHE STRING
    "Extra mpiexec flags for tests (;-separated), added after the auto-detected ones")

if(USE_DOUBLE)
  set(_prec double)
  set(_peer_prec float)
else()
  set(_prec float)
  set(_peer_prec double)
endif()
set(_features "")
if(USE_HDF5)
  list(APPEND _features hdf5)
endif()
if(USE_PROJ)
  list(APPEND _features proj)
endif()
if(USE_FFTW3)
  list(APPEND _features fftw)
endif()
string(REPLACE ";" "," _features "${_features}")

set(_harness ${UV_EXECUTABLE} run --quiet --frozen --project ${CMAKE_CURRENT_SOURCE_DIR}/harness sw4test)
set(_config ${CMAKE_CURRENT_BINARY_DIR}/sw4test.json)

# Re-run this configure step whenever a case or the harness changes.
file(GLOB_RECURSE _case_files CONFIGURE_DEPENDS
     ${CMAKE_CURRENT_SOURCE_DIR}/cases/*.toml
     ${CMAKE_CURRENT_SOURCE_DIR}/harness/src/*.py
     ${CMAKE_CURRENT_SOURCE_DIR}/harness/uv.lock)
set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS ${_case_files})

execute_process(
  COMMAND ${_harness} configure
          --out ${_config}
          --sw4 ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/sw4
          --precision ${_prec}
          --features "${_features}"
          --source-dir ${PROJECT_SOURCE_DIR}
          --out-dir ${CMAKE_CURRENT_BINARY_DIR}/out
          --mpiexec ${MPIEXEC_EXECUTABLE}
          "--numproc-flag=${MPIEXEC_NUMPROC_FLAG}"
          "--mpi-preflags=${MPIEXEC_PREFLAGS};${SW4TEST_MPI_PREFLAGS}"
          --uv ${UV_EXECUTABLE}
          "--peer-out-dir=${SW4TEST_PEER_OUT_DIR}"
          --peer-precision ${_peer_prec}
  RESULT_VARIABLE _rc)
if(NOT _rc EQUAL 0)
  message(FATAL_ERROR "sw4test configure failed (${_rc})")
endif()

execute_process(
  COMMAND ${_harness} cmake
          --config ${_config} --out ${CMAKE_CURRENT_BINARY_DIR}/sw4tests.cmake
  RESULT_VARIABLE _rc)
if(NOT _rc EQUAL 0)
  message(FATAL_ERROR "sw4test test registration failed (${_rc})")
endif()
include(${CMAKE_CURRENT_BINARY_DIR}/sw4tests.cmake)

