#
set(BLACS_FOUND FALSE)


# use variable
IF (BLACS_LIBRARIES)
  set(BLACS_FOUND TRUE)
  message(STATUS "set BLACS libraries from cmake variable: ${BLACS_LIBRARIES}")
  return()
endif()


# use ENV
set(HAVE_ENV_BLACS_LIBRARIES $ENV{BLACS_LIBRARIES})
if (HAVE_ENV_BLACS_LIBRARIES)
  set(BLACS_FOUND TRUE)
  set(BLACS_LIBRARIES $ENV{BLACS_LIBRARIES})
  message(STATUS "set B:ACS libraries from ENV: ${BLACS_LIBRARIES}")
  return()
endif()


# find libraries
message(STATUS "Finding BLACS")

set(BLACS_LIBRARIES_PATH
  "${BLACS_ROOT}"
  "${BLACS_ROOT}/lib"
  "$ENV{BLACS_ROOT}"
  "$ENV{BLACS_ROOT}/lib"
  "${CMAKE_SOURCE_DIR}"
  "${CMAKE_SOURCE_DIR}/lib"
  "${CMAKE_SOURCE_DIR}/blacs"
  "${CMAKE_SOURCE_DIR}/blacs/lib"
  /usr/lib
  /usr/local/lib
  )

find_library(
  BLACS_LIBRARIES
  NAMES "blacs" "blacs-mpi" "blacs-mpich" "blacs-mpich2" "blacs-openmpi" "blacs-lam" "blacs-pvm"
  PATHS "${BLACS_LIBRARIES_PATH}"
  )

find_library(
  BLACS_CINIT_LIBRARIES
  NAMES "blacsCinit" "blacsCinit-mpi" "blacsCinit-mpich" "blacsCinit-mpich2" "blacsCinit-openmpi" "blacsCinit-lam" "blacsCinit-pvm"
  PATHS "${BLACS_LIBRARIES_PATH}"
  )

find_library(
  BLACS_F77INIT_LIBRARIES
  NAMES "blacsF77init" "blacsF77init-mpi" "blacsF77init-mpich" "blacsF77init-mpich2" "blacsF77init-openmpi" "blacsF77init-lam" "blacsF77init-pvm"
  PATHS "${BLACS_LIBRARIES_PATH}"
  )

if (BLACS_LIBRARIES)
  set(BLACS_FOUND TRUE)
endif()

if (BLACS_FOUND)
  if (NOT BLACS_FIND_QUIETLY)
    message(STATUS "BLACS libraries: ${BLACS_LIBRARIES}")
  endif()
else()
  if (BLACS_FIND_REQUIRED)
    message(FATAL_ERROR "BLACS libraries not found")
  endif()
endif()

mark_as_advanced(
  BLACS_FOUND
  BLACS_LIBRARIES
  BLACS_CINIT_LIBRARIES
  BLACS_F77INIT_LIBRARIES
  )
