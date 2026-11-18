# FindFFTW3.cmake - Find FFTW3 library (including MPI support)
#
# Variables defined:
#   FFTW3_FOUND        - True if FFTW3 found
#   FFTW3_INCLUDE_DIRS - Include directories for FFTW3
#   FFTW3_LIBRARIES    - Libraries for FFTW3 (core only)
#   FFTW3_MPI_LIBRARIES - Libraries for FFTW3 MPI
#
# Imported targets:
#   FFTW3::FFTW3       - Core FFTW3 library
#   FFTW3::FFTW3_MPI   - FFTW3 MPI library (links FFTW3::FFTW3 transitively)

find_path(FFTW3_INCLUDE_DIR
  NAMES fftw3.h fftw3-mpi.h
  DOC "FFTW3 include directory"
)

find_library(FFTW3_LIBRARY
  NAMES fftw3
  DOC "FFTW3 library"
)

find_library(FFTW3_MPI_LIBRARY
  NAMES fftw3_mpi
  DOC "FFTW3 MPI library"
)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(FFTW3
  REQUIRED_VARS FFTW3_LIBRARY FFTW3_INCLUDE_DIR
)

if(FFTW3_FOUND)
  set(FFTW3_INCLUDE_DIRS ${FFTW3_INCLUDE_DIR})
  set(FFTW3_LIBRARIES ${FFTW3_LIBRARY})

  if(NOT TARGET FFTW3::FFTW3)
    add_library(FFTW3::FFTW3 UNKNOWN IMPORTED)
    set_target_properties(FFTW3::FFTW3 PROPERTIES
      IMPORTED_LOCATION "${FFTW3_LIBRARY}"
      INTERFACE_INCLUDE_DIRECTORIES "${FFTW3_INCLUDE_DIR}"
    )
  endif()

  if(FFTW3_MPI_LIBRARY AND NOT TARGET FFTW3::FFTW3_MPI)
    add_library(FFTW3::FFTW3_MPI UNKNOWN IMPORTED)
    set_target_properties(FFTW3::FFTW3_MPI PROPERTIES
      IMPORTED_LOCATION "${FFTW3_MPI_LIBRARY}"
      INTERFACE_INCLUDE_DIRECTORIES "${FFTW3_INCLUDE_DIR}"
      INTERFACE_LINK_LIBRARIES "FFTW3::FFTW3"
    )
    set(FFTW3_MPI_LIBRARIES ${FFTW3_MPI_LIBRARY})
    message(STATUS "Found FFTW3 MPI: ${FFTW3_MPI_LIBRARY}")
  endif()

  mark_as_advanced(FFTW3_INCLUDE_DIR FFTW3_LIBRARY FFTW3_MPI_LIBRARY)
endif()
