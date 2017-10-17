cmake_minimum_required(VERSION 3.5)
project(ProteinDF)
enable_language(CXX)

# Project Version
set(PROJECT_VERSION_MAJOR "2017")
set(PROJECT_VERSION_MINOR "10")
set(PROJECT_VERSION_REVISION "0")
set(PROJECT_VERSION "${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}.${PROJECT_VERSION_REVISION}")

set(PROJECT_VERSION_BUILD "")
execute_process(COMMAND git rev-parse --short HEAD OUTPUT_VARIABLE BUILD_GIT_REV OUTPUT_STRIP_TRAILING_WHITESPACE)
execute_process(COMMAND date +%Y%m%d OUTPUT_VARIABLE BUILD_DATE OUTPUT_STRIP_TRAILING_WHITESPACE)
if (DEFINED BUILD_GIT_REV)
  set(PROJECT_VERSION_BUILD "${BUILD_GIT_REV}")
else()
  set(PROJECT_VERSION_BUILD "${BUILD_DATE}")
endif()
#message("BUILD_GIT_REV: ${BUILD_GIT_REV}")
#message("BUILD_DATE: ${BUILD_DATE}")
#message("PROJECT_VERSION_BUILD: ${PROJECT_VERSION_BUILD}")

#
include(CheckTypeSize)
include(CheckIncludeFile)


#set(CMAKE_VERBOSE_MAKEFILE 1)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake-modules")


# check compiler flags
message(STATUS "CMAKE_CXX_FLAGS: ${CMAKE_CXX_FLAGS}")
#message(STATUS "CMAKE_CXX_FLAGS_RELEASE: ${CMAKE_CXX_FLAGS_RELEASE}")
#message(STATUS "CMAKE_CXX_FLAGS_DEBUG: ${CMAKE_CXX_FLAGS_DEBUG}")
#message(STATUS "CMAKE_CXX_FLAGS_RELWITHDEBINFO: ${CMAKE_CXX_FLAGS_RELWITHDEBINFO}")
#message(STATUS "CMAKE_CXX_FLAGS_MINSIZEREL: ${CMAKE_CXX_FLAGS_MINSIZEREL}")


# check system
check_include_file("sys/resource.h" HAVE_SYS_RESOURCE_H)
check_include_file("sys/time.h" HAVE_SYS_TIME_H)
check_include_file("sys/types.h" HAVE_SYS_TYPES_H)
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
  #set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${BLAS_LINKER_FLAGS}")
  #set(CMAKE_MODULE_LINKER_FLAGS "${CMAKE_MODULE_LINKER_FLAGS} ${BLAS_LINKER_FLAGS}")
  #set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} ${BLAS_LINKER_FLAGS}")

  message("BLAS_LINKER_FLAGS ${BLAS_LINKER_FLAGS}")
  message("BLAS_LIBRARIES ${BLAS_LIBRARIES}")
  set(HAVE_BLAS 1) # for config.h
endif(BLAS_FOUND)

find_package(LAPACK REQUIRED)
if (LAPACK_FOUND)
  #set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${LAPACK_LINKER_FLAGS}")
  #set(CMAKE_MODULE_LINKER_FLAGS "${CMAKE_MODULE_LINKER_FLAGS} ${LAPACK_LINKER_FLAGS}")
  #set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} ${LAPACK_LINKER_FLAGS}")

  message("LAPACK_LINKER_FLAGS ${LAPACK_LINKER_FLAGS}")
  message("LAPACK_LIBRARIES ${LAPACK_LIBRARIES}")
  set(HAVE_LAPACK 1) # for config.h
endif(LAPACK_FOUND)


# check MPI
find_package(MPI)
message("MPI information:")
message("  MPI_INCLUDE_PATH: ${MPI_INCLUDE_PATH}")
message("  MPI_LIBRARY: ${MPI_LIBRARY}")
message("  MPI_LIBRARIES: ${MPI_LIBRARIES}")
message("  MPIEXEC: ${MPIEXEC}")
message("  MPI_CXX_COMPILER: ${MPI_CXX_COMPILER}")
message("  MPI_CXX_LIBRARIES: ${MPI_CXX_LIBRARIES}")
message("  MPI_CXX_INCLUDE_PATH: ${MPI_CXX_INCLUDE_PATH}")
message("  MPI_CXX_FOUND: ${MPI_CXX_FOUND}")
if (MPI_FOUND)
  message("overwrite: CMAKE_CXX_COMPILER = ${MPI_CXX_COMPILER}")
  set(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER})
else(MPI_FOUND)
  message("Not Found MPI compiler")
endif(MPI_FOUND)

  # check ScaLAPACK
  find_package(ScaLAPACK)
  if (SCALAPACK_FOUND)
    # message("ScaLAPACK information:")
    # message("  SCALAPACK_LIBRARIES: ${SCALAPACK_LIBRARIES}")
    set(HAVE_SCALAPACK 1) # for config.h
  endif(SCALAPACK_FOUND)
endif(MPI_FOUND)

# check ScaLAPACK
find_package(ScaLAPACK)
if (SCALAPACK_FOUND)
  message("ScaLAPACK information")
  message("  SCALAPACK_LIBRARIES: ${SCALAPACK_LIBRARIES}")

  # add_definitions(-DHAVE_SCALAPACK)
  set(HAVE_SCALAPACK 1) # for config.h
else(SCALAPACK_FOUND)
  message("Not Found ScaLAPACK")
endif(SCALAPACK_FOUND)

# check HDF5
find_package(HDF5 COMPONENTS C CXX)
if (HDF5_FOUND)
  include_directories(${HDF5_INCLUDE_DIRS})
  set(HAVE_HDF5 1)
  message("HDF5 information:")
  message("  HDF5_LIBRARIES: ${HDF5_LIBRARIES}")
endif(HDF5_FOUND)

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
set(CPACK_PACKAGE_VERSION_MAJOR "${PROJECT_VERSION_MAJOR}")
set(CPACK_PACKAGE_VERSION_MINOR "${PROJECT_VERSION_MINOR}")
set(CPACK_PACKAGE_VERSION_PATCH "${PROJECT_VERSION_REVISION}")
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
