# 
set(SCALAPACK_FOUND FALSE)

IF (SCALAPACK_LIBRARIES)
  set(SCALAPACK_FOUND  TRUE)
  message(STATUS "set ScaLAPACK libraries: ${SCALAPACK_LIBRARIES}")
  return()
endif()

message(STATUS "Finding ScaLAPACK")

set(SCALAPACK_LIBRARIES_PATH
  "${SCALAPACK_ROOT}"
  "${SCALAPACK_ROOT}/lib"
  "$ENV{SCALAPACK_ROOT}"
  "$ENV{SCALAPACK_ROOT}/lib"
  "${CMAKE_SOURCE_DIR}"
  "${CMAKE_SOURCE_DIR}/lib"
  "${CMAKE_SOURCE_DIR}/scalapack"
  "${CMAKE_SOURCE_DIR}/scalapack/lib"
  /usr/lib
  /usr/local/lib
  )
  
  
find_library(
  SCALAPACK_LIBRARIES
  NAMES "scalapack" "scalapck-mpi" "scalapack-mpich" "scalapack-mpich2" "scalapack-openmpi" "scalapack-lam" "scalapack-pvm"
  PATHS ${SCALAPACK_LIBRARIES_PATH}
  )

if (SCALAPACK_LIBRARIES)
  set(SCALAPACK_FOUND TRUE)
endif()

if (SCALAPACK_FOUND)
  if (NOT SCALAPACK_FIND_QUIETLY)
    message(STATUS "ScaLAPACK libraries: ${SCALAPACK_LIBRARIES}")
  endif()
else()
  if (SCALAPACK_FIND_REQUIRED)
    message(FATAL_ERROR "ScaLAPACK libraries not found")
  endif()
endif()

mark_as_advanced(
  SCALAPACK_FOUND
  SCALAPACK_LIBRARIES
  )
