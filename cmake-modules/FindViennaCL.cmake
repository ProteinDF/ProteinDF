# Try to find ViennaCL
#
# usage:
# find_package(ViennaCL)
#
# Once done this will define
#  VIENNACL_FOUND - system has ViennaCL lib with correct version
#  VIENNACL_INCLUDE_DIR - the ViennaCL include directory
#
# This module reads hints about search locations from
# the following enviroment variables:
#
# VIENNACL_ROOT

set(HAVE_ENV_VIENNACL_ROOT $ENV{VIENNACL_ROOT})
if (HAVE_ENV_VIENNACL_ROOT)
    set(VIENNACL_ROOT $ENV{VIENNACL_ROOT})
    message(STATUS "VIENNACL_ROOT: ${VIENNACL_ROOT}")
endif ()

SET(VIENNACL_INCLUDE_SEARCH_PATHS
  ..
  /usr/include
  /usr/local/include
  /opt/ViennaCL/include
  ${VIENNACL_ROOT}
  ${VIENNACL_ROOT}/include
  # $ENV{VIENNACL_ROOT}
  # $ENV{VIENNACL_ROOT}/include
  ENV VIENNACL_ROOT
  )

find_path(
  VIENNACL_INCLUDE_DIR
  NAMES viennacl/vector.hpp
  PATHS ${VIENNACL_INCLUDE_SEARCH_PATHS}
  )

set(VIENNACL_FOUND OFF)
if (VIENNACL_INCLUDE_DIR)
  message(STATUS "VIENNACL_INCLUDE_DIR: ${VIENNACL_INCLUDE_DIR}")
  set(VIENNACL_FOUND ON)
else (VIENNACL_INCLUDE_DIR)
  message(STATUS "not found VIENNACL_INCLUDE_DIR")
endif(VIENNACL_INCLUDE_DIR)

set(VIENNACL_INCLUDE_DIR ${VIENNACL_INCLUDE_DIR})
mark_as_advanced(VIENNACL_INCLUDE_DIR)
