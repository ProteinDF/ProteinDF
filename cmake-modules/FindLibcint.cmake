# Try to find Libcint
#
# usage:
# find_package(Libcint)
#
# Once done this will define
#  LIBCINT_FOUND - system has libcint lib with correct version
#  LIBCINT_INCLUDE_DIR - the libcint include directory
#
# This module reads hints about search locations from
# the following enviroment variables:
#
# LIBCINT_ROOT

SET(LIBCINT_INCLUDE_SEARCH_PATHS
  ..
  /usr/include
  /usr/local/include
  /opt/libcint/include
  ${LIBCINT_ROOT}
  ${LIBCINT_ROOT}/include
  $ENV{LIBCINT_ROOT}
  $ENV{LIBCINT_ROOT}/include
  )

SET(LIBCINT_LIBRARIES_PATHS
  "${LIBCINT_ROOT}"
  "${LIBCINT_ROOT}/lib"
  "$ENV{LIBCINT_ROOT}"
  "$ENV{LIBCINT_ROOT}/lib"
  "${CMAKE_SOURCE_DIR}"
  "${CMAKE_SOURCE_DIR}/lib"
  "${CMAKE_SOURCE_DIR}/libcint"
  "${CMAKE_SOURCE_DIR}/libcint/lib"
  /usr/lib
  /usr/local/lib
  )


message(STATUS "Finding Libcint")

find_path(
  LIBCINT_INCLUDE_DIR
  NAMES cint.h
  PATHS ${LIBCINT_INCLUDE_SEARCH_PATHS}
  )
set(LIBCINT_INCLUDE_DIR ${LIBCINT_INCLUDE_DIR})

find_library(
  LIBCINT_LIBRARIES
  NAMES libcint cint
  PATHS ${LIBCINT_LIBRARIES_PATHS}
  )
message("lib: ${LIBCINT_LIBRARIES}")
set(LIBCINT_FOUND OFF)
if (LIBCINT_LIBRARIES)
  set(LIBCINT_FOUND ON)
endif(LIBCINT_LIBRARIES)

mark_as_advanced(
    LIBCINT_INCLUDE_DIR
    LIBCINT_LIBRARIES
    )
