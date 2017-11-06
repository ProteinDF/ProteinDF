cmake_minimum_required(VERSION 3.5)
project(ProteinDF)
enable_language(CXX)

include(CheckTypeSize)
include(CheckIncludeFile)

#set(CMAKE_VERBOSE_MAKEFILE 1)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake-modules")

# set(CMAKE_C_FLAGS "-Wall")
# set(CMAKE_CXX_FLAGS "-Wall")
# message(STATUS "${CMAKE_CXX_FLAGS_RELWITHDEBINFO}")

# check system
check_include_file("sys/resource.h" HAVE_SYS_RESOURCE_H)
check_include_file("sys/time.h" HAVE_SYS_TIME_H)
check_include_file("unistd.h" HAVE_UNISTD_H)
check_include_file("stdint.h" HAVE_STDINT_H)


# check OpenMP
find_package(OpenMP)
if (OPENMP_FOUND)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif(OPENMP_FOUND)


# check linear algebra package
find_package(BLAS REQUIRED)
if (BLAS_FOUND)
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${BLAS_LINKER_FLAGS}")
  set(CMAKE_MODULE_LINKER_FLAGS "${CMAKE_MODULE_LINKER_FLAGS} ${BLAS_LINKER_FLAGS}")
  set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} ${BLAS_LINKER_FLAGS}")

  add_definitions(-DHAVE_BLAS)
endif(BLAS_FOUND)

find_package(LAPACK REQUIRED)
if (LAPACK_FOUND)
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${LAPACK_LINKER_FLAGS}")
  set(CMAKE_MODULE_LINKER_FLAGS "${CMAKE_MODULE_LINKER_FLAGS} ${LAPACK_LINKER_FLAGS}")
  set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} ${LAPACK_LINKER_FLAGS}")

  add_definitions(-DHAVE_LAPACK)
endif(LAPACK_FOUND)


# check MPI
find_package(MPI)
if (MPI_FOUND)
  message("MPI information:")
  message("  MPI_INCLUDE_PATH: ${MPI_INCLUDE_PATH}")
  message("  MPI_LIBRARY: ${MPI_LIBRARY}")
  message("  MPI_LIBRARIES: ${MPI_LIBRARIES}")
  message("  MPIEXEC: ${MPIEXEC}")
  message("  MPI_C_COMPILER: ${MPI_C_COMPILER}")
  message("  MPI_CXX_COMPILER: ${MPI_CXX_COMPILER}")

  set(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER})
endif(MPI_FOUND)


# check ScaLAPACK
find_package(ScaLAPACK)
if (SCALAPACK_FOUND)
  message("ScaLAPACK information")
  message("  SCALAPACK_LIBRARIES: ${SCALAPACK_LIBRARIES}")

  add_definitions(-DHAVE_SCALAPACK)
endif(SCALAPACK_FOUND)


# check threads
set(CMAKE_THREAD_PREFER_PTHREAD ON)
set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)


# make config.h
configure_file(${CMAKE_SOURCE_DIR}/config.h.in ${CMAKE_BINARY_DIR}/config.h @ONLY)
include_directories(${CMAKE_BINARY_DIR})
add_definitions(-DHAVE_CONFIG_H)


# do subdirectories
add_subdirectory(src)
add_subdirectory(data)
add_subdirectory(scripts)

# package
set(CPACK_PACKAGE_VERSION_MAJOR 2017)
set(CPACK_PACKAGE_VERSION_MINOR 10)
set(CPACK_PACKAGE_VERSION_PATCH 0)
set(CPACK_PACKAGE_VENDOR "ProteinDF development team")
set(CPACK_SOURCE_GENERATOR TGZ)
set(CPACK_SOURCE_IGNORE_FILES
  "~$"
  "\.gitignore$"
  "\.git/*"
  "^${PROJECT_SOURCE_DIR}/build"
  "^${PROJECT_SOURCE_DIR}/temp*"
  "^${PROJECT_SOURCE_DIR}/tmp*"
  "^${PROJECT_SOURCE_DIR}/out\.*"
)
include(CPack)


