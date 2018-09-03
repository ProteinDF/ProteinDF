# Try to find ViennaCL
#
# usage:
# find_package(ViennaCL)
#
# Once done this will define
#  EIGEN3_FOUND - system has eigen lib with correct version
#  EIGEN3_INCLUDE_DIR - the eigen include directory
#
# This module reads hints about search locations from
# the following enviroment variables:
#
# VIENNACL_ROOT

SET(VIENNACL_INCLUDE_SEARCH_PATHS
  ..
  /usr/include
  /usr/local/include
  /opt/ViennaCL/include
  ${VIENNACL_ROOT}
  ${VIENNACL_ROOT}/include
  $ENV{VIENNACL_ROOT}
  $ENV{VIENNACL_ROOT}/include
  )

find_path(
  VIENNACL_INCLUDE_DIR
  NAMES viennacl/vector.hpp
  PATHS ${VIENNACL_INCLUDE_SEARCH_PATHS}
  )

set(VIENNACL_FOUND OFF)
if (VIENNACL_INCLUDE_DIR)
  set(VIENNACL_FOUND ON)
endif(VIENNACL_INCLUDE_DIR)

set(VIENNACL_INCLUDE_DIR ${VIENNACL_INCLUDE_DIR})
mark_as_advanced(VIENNACL_INCLUDE_DIR)
