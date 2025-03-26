cmake_minimum_required(VERSION 3.5)
project(ProteinDF)
enable_language(C)
enable_language(CXX)
enable_testing()

# for CMake debug
# set(CMAKE_VERBOSE_MAKEFILE TRUE)

# Project Version
set(PROJECT_VERSION_MAJOR "2025")
set(PROJECT_VERSION_MINOR "3")
set(PROJECT_VERSION_REVISION "0")
set(PROJECT_VERSION "${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}.${PROJECT_VERSION_REVISION}")

set(PROJECT_VERSION_BUILD "")
execute_process(
  COMMAND git rev-parse --abbrev-ref HEAD
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  OUTPUT_VARIABLE GIT_BRANCH
  OUTPUT_STRIP_TRAILING_WHITESPACE
  )
execute_process(
  COMMAND git rev-parse --short HEAD
  OUTPUT_VARIABLE BUILD_GIT_REV
  OUTPUT_STRIP_TRAILING_WHITESPACE
  )
execute_process(
  COMMAND date +%Y%m%d
  OUTPUT_VARIABLE BUILD_DATE
  OUTPUT_STRIP_TRAILING_WHITESPACE
  )
message("GIT_BRANCH: ${GIT_BRANCH}")
message("BUILD_GIT_REV: ${BUILD_GIT_REV}")
message("BUILD_DATE: ${BUILD_DATE}")
set(PROJECT_VERSION_BUILD "${GIT_BRANCH}:${BUILD_GIT_REV}:${BUILD_DATE}")
message("PROJECT_VERSION_BUILD: ${PROJECT_VERSION_BUILD}")

# -----------------------------------------------------------------------------
# Option
# -----------------------------------------------------------------------------
OPTION(USE_EIGEN "use EIGEN" on)
OPTION(USE_VIENNACL "use ViennaCL" off)
OPTION(USE_HDF5 "use HDF5" off)
OPTION(USE_LIBCINT "use LIBCINT" off)

# -----------------------------------------------------------------------------
# Check
# -----------------------------------------------------------------------------
include(CheckTypeSize)
include(CheckIncludeFile)

# -----------------------------------------------------------------------------
# Modules
# -----------------------------------------------------------------------------
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake-modules")

# C++ compiler
if(NOT "${CMAKE_CXX_STANDARD}")
  set(CMAKE_CXX_STANDARD 14)
endif()

# check compiler flags
set(CMAKE_CXX_FLAGS_RELEASE "-DNDEBUG")
set(CMAKE_CXX_FLAGS_DEBUG "-Wall -g")
#message(STATUS "CMAKE_CXX_FLAGS: ${CMAKE_CXX_FLAGS}")
#message(STATUS "CMAKE_CXX_FLAGS_RELEASE: ${CMAKE_CXX_FLAGS_RELEASE}")
#message(STATUS "CMAKE_CXX_FLAGS_DEBUG: ${CMAKE_CXX_FLAGS_DEBUG}")
#message(STATUS "CMAKE_CXX_FLAGS_RELWITHDEBINFO: ${CMAKE_CXX_FLAGS_RELWITHDEBINFO}")
#message(STATUS "CMAKE_CXX_FLAGS_MINSIZEREL: ${CMAKE_CXX_FLAGS_MINSIZEREL}")


# check system
check_include_file("sys/resource.h" HAVE_SYS_RESOURCE_H)
check_include_file("sys/time.h" HAVE_SYS_TIME_H)
check_include_file("sys/types.h" HAVE_SYS_TYPES_H)
check_include_file("sys/mman.h" HAVE_SYS_MMAN_H)
check_include_file("sys/stat.h" HAVE_SYS_STAT_H)
check_include_file("unistd.h" HAVE_UNISTD_H)
check_include_file("stdint.h" HAVE_STDINT_H)


# -----------------------------------------------------------------------------
# Boost
# -----------------------------------------------------------------------------
find_package(Boost)
if (Boost_FOUND)
  message("Boost information:")
  message("  Boost_VERSION: ${Boost_VERSION}")
  message("  Boost_INCLUDE_DIRS: ${Boost_INCLUDE_DIRS}")
  message("  Boost_LIBRARIES: ${Boost_LIBRARIES}")

  include_directories(${Boost_INCLUDE_DIRS})
else (Boost_FOUND)
  message("Not Found Boost")
endif (Boost_FOUND)


# -----------------------------------------------------------------------------
# OpenMP
# -----------------------------------------------------------------------------
set(HAVE_ENV_OPENMP_HOME $ENV{OpenMP_HOME})
if (HAVE_ENV_OPENMP_HOME)
    set(OpenMP_HOME $ENV{OpenMP_HOME})
    message(STATUS "OpenMP Root : ${OpenMP_HOME}")
endif ()

set(HAVE_OpenMP_C_LIB_NAMES $ENV{OpenMP_C_LIB_NAMES})
if (HAVE_OpenMP_C_LIB_NAMES)
    message(STATUS "set OpenMP_C_LIB_NAMES: $ENV{OpenMP_C_LIB_NAMES}")
    set(OpenMP_C_LIB_NAMES $ENV{OpenMP_C_LIB_NAMES})
endif (HAVE_OpenMP_C_LIB_NAMES)
set(HAVE_OpenMP_CXX_LIB_NAMES $ENV{OpenMP_CXX_LIB_NAMES})
if (HAVE_OpenMP_CXX_LIB_NAMES)
    message(STATUS "set OpenMP_CXX_LIB_NAMES: $ENV{OpenMP_CXX_LIB_NAMES}")
    set(OpenMP_CXX_LIB_NAMES $ENV{OpenMP_CXX_LIB_NAMES})
endif (HAVE_OpenMP_CXX_LIB_NAMES)
set(HAVE_OpenMP_omp_LIBRARY $ENV{OpenMP_omp_LIBRARY})
if (HAVE_OpenMP_omp_LIBRARY)
    message(STATUS "set OpenMP_omp_LIBRARY: $ENV{OpenMP_omp_LIBRARY}")
    set(OpenMP_omp_LIBRARY $ENV{OpenMP_omp_LIBRARY})
endif (HAVE_OpenMP_omp_LIBRARY)

set(HAVE_ENV_OpenMP_C_FLAGS $ENV{OpenMP_C_FLAGS})
if (HAVE_ENV_OpenMP_C_FLAGS)
    message(STATUS "set OpenMP_C_FLAGS: $ENV{OpenMP_C_FLAGS}")
    set(OpenMP_C_FLAGS $ENV{OpenMP_C_FLAGS} CACHE STRING "" FORCE)
endif ()
set(HAVE_ENV_OpenMP_CXX_FLAGS $ENV{OpenMP_CXX_FLAGS})
if (HAVE_ENV_OpenMP_CXX_FLAGS)
    message(STATUS "set OpenMP_CXX_FLAGS: $ENV{OpenMP_CXX_FLAGS}")
    set(OpenMP_CXX_FLAGS $ENV{OpenMP_CXX_FLAGS} CACHE STRING "" FORCE)
endif ()

find_package(OpenMP)
if (OPENMP_FOUND)
  message("OpenMP information:")
  message("  OpenMP_VERSION: ${OpenMP_VERSION}")
  message("  OpenMP_CXX_FLAGS: ${OpenMP_CXX_FLAGS}")
  message("  OpenMP_CXX_LIB_NAMES: ${OpenMP_CXX_LIB_NAMES}")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  add_definitions(-DVIENNACL_WITH_OPENMP)
else (OPENMP_FOUND)
message("OpenMP not found")
endif (OPENMP_FOUND)


# -----------------------------------------------------------------------------
# MPI
# -----------------------------------------------------------------------------
find_package(MPI)
if (MPI_FOUND)
  if (DEFINED ENV{MPIEXEC_FLAGS})
    set(MPIEXEC_FLAGS $ENV{MPIEXEC_FLAGS})
  endif (DEFINED ENV{MPIEXEC_FLAGS})

  message("MPI information:")
  message("  MPI_VERSION: ${MPI_VERSION}")
  message("  MPIEXEC_EXECUTABLE: ${MPIEXEC_EXECUTABLE}")
  message("  MPIEXEC: ${MPIEXEC}")
  message("  MPIEXEC_FLAGS: ${MPIEXEC_FLAGS}")
  message("  MPIEXEC_PREFLAGS: ${MPIEXEC_PREFLAGS}")
  message("  MPIEXEC_NUMPROC_FLAG: ${MPIEXEC_NUMPROC_FLAG}")
  message("  MPIEXEC_POSTFLAGS: ${MPIEXEC_POSTFLAGS}")
  message("  MPI_CXX_COMPILER: ${MPI_CXX_COMPILER}")
  message("  MPI_CXX_LINK_FLAGS: ${MPI_CXX_LINK_FLAGS}")
  message("  MPI_CXX_LIBRARIES: ${MPI_CXX_LIBRARIES}")
  message("  MPI_CXX_INCLUDE_PATH: ${MPI_CXX_INCLUDE_PATH}")
  message("  MPI_CXX_FOUND: ${MPI_CXX_FOUND}")

  set(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER})
  message("overwrite: CMAKE_CXX_COMPILER = ${CMAKE_CXX_COMPILER}")
else(MPI_FOUND)
  message("Not Found MPI compiler")
endif(MPI_FOUND)


# check linear algebra package
if (USE_EIGEN)
  find_package(Eigen3)
  if (EIGEN3_FOUND)
    message("Eigen3 information:")
    message("  EIGEN3_INCLUDE_DIR: ${EIGEN3_INCLUDE_DIR}")
    message("  EIGEN3_VERSION: ${EIGEN3_VERSION}")

    set(HAVE_EIGEN 1) # for config.h
    include_directories(${EIGEN3_INCLUDE_DIR})
  endif (EIGEN3_FOUND)
else (USE_EIGEN)
  message("USE_EIGEN: ${USE_EIGEN}")
endif (USE_EIGEN)

find_package(BLAS REQUIRED)
if (BLAS_FOUND)
  message("BLAS information:")
  message("  BLA_VENDOR: ${BLA_VENDOR}")
  message("  BLA_STATIC: ${BLA_STATIC}")
  message("  BLAS_LINKER_FLAGS: ${BLAS_LINKER_FLAGS}")
  message("  BLAS_LIBRARIES: ${BLAS_LIBRARIES}")

  set(HAVE_BLAS 1) # for config.h
endif(BLAS_FOUND)

find_package(LAPACK REQUIRED)
if (LAPACK_FOUND)
  message("LAPACK information:")
  message("  LAPACK_LINKER_FLAGS: ${LAPACK_LINKER_FLAGS}")
  message("  LAPACK_LIBRARIES: ${LAPACK_LIBRARIES}")

  set(HAVE_LAPACK 1) # for config.h
endif(LAPACK_FOUND)


# -----------------------------------------------------------------------------
# OpenCL
# -----------------------------------------------------------------------------
find_package(OpenCL)
if (OpenCL_FOUND)
  message("OpenCL information:")
  message("  OpenCL_VERSION_STRING: ${OpenCL_VERSION_STRING}")
  message("  OpenCL_INCLUDE_DIRS: ${OpenCL_INCLUDE_DIRS}")
  message("  OpenCL_LIBRARIES: ${OpenCL_LIBRARIES}")

  include_directories(${OpenCL_INCLUDE_DIRS})
  set(HAVE_OPENCL 1) # for config.h
endif (OpenCL_FOUND)

# -----------------------------------------------------------------------------
# ViennaCL
# -----------------------------------------------------------------------------
if (USE_VIENNACL)

  # for boost
  find_package(Boost REQUIRED)
  if (Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS})
    message("Boost_INCLUDE_DIRS: ${Boost_INCLUDE_DIRS}")
    message("Boost_LIBRARIES: ${Boost_LIBRARIES}")
  endif (Boost_FOUND)

  find_package(ViennaCL)
  if (VIENNACL_FOUND)
    message("ViennaCL information:")
    message("  VIENNACL_INCLUDE_DIR: ${VIENNACL_INCLUDE_DIR}")

    if (OpenCL_FOUND)
      add_definitions(-DVIENNACL_WITH_OPENCL)
    endif (OpenCL_FOUND)

    set(HAVE_VIENNACL 1) # for config.h
    include_directories(${VIENNACL_INCLUDE_DIR})
  else (VIENNACL_FOUND)
    message("Not Found ViennaCL")
    message("  please check: VIENNACL_ROOT")
  endif (VIENNACL_FOUND)
else (USE_VIENNACL)
  message("USE_VIENNACL: ${USE_VIENNACL}")
endif (USE_VIENNACL)


# -----------------------------------------------------------------------------
# BLACS
# -----------------------------------------------------------------------------
find_package(BLACS)
if (BLACS_FOUND)
  message("BLACS information")
  message("  BLACS_LIBRARIES: ${BLACS_LIBRARIES}")
  message("  BLACS_CINIT_LIBRARIES: ${BLACS_CINIT_LIBRARIES}")
  message("  BLACS_F77INIT_LIBRARIES: ${BLACS_F77INIT_LIBRARIES}")
else(BLACS_FOUND)
  message("Not Found BLACS")
  set(BLACS_LIBRARIES "")
  set(BLACS_CINIT_LIBRARIES "")
  set(BLACS_F77INIT_LIBRARIES "")
endif (BLACS_FOUND)


# -----------------------------------------------------------------------------
# ScaLAPACK
# -----------------------------------------------------------------------------
find_package(ScaLAPACK)
if (SCALAPACK_FOUND)
  message("ScaLAPACK information")
  message("  SCALAPACK_LIBRARIES: ${SCALAPACK_LIBRARIES}")

  set(HAVE_SCALAPACK 1) # for config.h
else(SCALAPACK_FOUND)
  message("Not Found ScaLAPACK")
endif(SCALAPACK_FOUND)


# -----------------------------------------------------------------------------
# HDF5
# -----------------------------------------------------------------------------
if (USE_HDF5)
  find_package(HDF5 COMPONENTS C CXX)
  if (HDF5_FOUND)
    include_directories(${HDF5_INCLUDE_DIRS})
    message("HDF5 information:")
    message("  HDF5_VERSION: ${HDF5_VERSION}")
    message("  HDF5_IS_PARALLEL: ${HDF5_IS_PARALLEL}")
    message("  HDF5_INCLUDE_DIRS: ${HDF5_INCLUDE_DIRS}")
    message("  HDF5_CXX_INCLUDE_DIRS: ${HDF5_CXX_INCLUDE_DIRS}")
    message("  HDF5_LIBRARIES: ${HDF5_LIBRARIES}")
    message("  HDF5_HL_LIBRARIES: ${HDF5_HL_LIBRARIES}")
    message("  HDF5_CXX_LIBRARIES: ${HDF5_CXX_LIBRARIES}")

    set(HAVE_HDF5 1) # for config.h
  else (HDF5_FOUND)
    message("Not Found HDF5")
  endif (HDF5_FOUND)
else (USE_HDF5)
  message("USE_HDF5: ${USE_HDF5}")
endif (USE_HDF5)


# -----------------------------------------------------------------------------
# libcint
# -----------------------------------------------------------------------------
if (USE_LIBCINT)
  if (BLAS_FOUND)
    find_package(libcint)
    if (LIBCINT_FOUND)
      message("LIBCINT information:")
      message("  LIBCINT_INCLUDE_DIR: ${LIBCINT_INCLUDE_DIR}")
      message("  LIBCINT_LIBRARIES: ${LIBCINT_LIBRARIES}")
      set(HAVE_LIBCINT 1) # for config.h
    else()
      message("Not Found LIBCINT")
    endif()
  endif ()
endif ()
# check threads
set(CMAKE_THREAD_PREFER_PTHREAD ON)
set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)


# -----------------------------------------------------------------------------
# Google Test
# -----------------------------------------------------------------------------
#   ref.) env param: `GTEST_ROOT`
find_package(GTest)
if (GTEST_FOUND)
  include_directories(${GTEST_INCLUDE_DIRS})
  message("Google Test information:")
  message("  GTEST_INCLUDE_DIRS: ${GTEST_INCLUDE_DIRS}")
  message("  GTEST_LIBRARIES: ${GTEST_LIBRARIES}")
  message("  GTEST_MAIN_LIBRARIES: ${GTEST_MAIN_LIBRARIES}")
endif (GTEST_FOUND)


# -----------------------------------------------------------------------------
# Doxygen
# -----------------------------------------------------------------------------
find_package(Doxygen)
if (DOXYGEN_FOUND)
  set(DOXYGEN_CONF_FILE "doxygen.conf")
  set(DOXYGEN_CONF_PATH ${CMAKE_CURRENT_BINARY_DIR}/${DOXYGEN_CONF_FILE})
  set(DOXYGEN_SOURCE_DIR ${CMAKE_SOURCE_DIR})
  set(DOXYGEN_TARGET "doc")

  configure_file(${DOXYGEN_CONF_FILE}.in ${DOXYGEN_CONF_PATH})

  add_custom_target(${DOXYGEN_TARGET}
                    ${DOXYGEN_EXECUTABLE} ${DOXYGEN_CONF_PATH}
                    DEPENDS  ${DOXYGEN_CONF_PATH})
else (DOXYGEN_FOUND)
  message("Doxygen not found")
endif (DOXYGEN_FOUND)


# make config.h
configure_file(${CMAKE_SOURCE_DIR}/config.h.in ${CMAKE_BINARY_DIR}/config.h @ONLY)
include_directories(${CMAKE_BINARY_DIR})
add_definitions(-DHAVE_CONFIG_H)


# do subdirectories
add_subdirectory(src)
add_subdirectory(data)
add_subdirectory(scripts)

# package
#   use: `make package_source` to make tarball
set(CPACK_PACKAGE_VERSION_MAJOR "${PROJECT_VERSION_MAJOR}")
set(CPACK_PACKAGE_VERSION_MINOR "${PROJECT_VERSION_MINOR}")
set(CPACK_PACKAGE_VERSION_PATCH "${PROJECT_VERSION_REVISION}")
set(CPACK_PACKAGE_VENDOR "ProteinDF development team")
set(CPACK_SOURCE_GENERATOR TGZ)
set(CPACK_SOURCE_IGNORE_FILES
  "~$"
  "\.gitignore$"
  "\.git/*"
  "\.vscode/*"
  "^${PROJECT_SOURCE_DIR}/build"
  "^${PROJECT_SOURCE_DIR}/temp*"
  "^${PROJECT_SOURCE_DIR}/tmp*"
  "^${PROJECT_SOURCE_DIR}/out\.*"
)
include(CPack)
