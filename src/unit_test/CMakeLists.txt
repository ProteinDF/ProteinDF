set(PDF_XTEST_SOURCES
  matrix_common.h
  vector_common.h

  DfFunctional_B3LYPTest.cpp
  DfFunctional_Becke88Test.cpp
  DfFunctional_LYPTest.cpp
  DfFunctional_SlaterTest.cpp
  DfFunctional_VWNTest.cpp
  TlFileTest.cpp
  TlMemManagerTest.cpp

  tl_dense_general_matrix_mmap_test.cc
  tl_dense_general_matrix_mmap_rowmajor_test.cc
  tl_dense_general_matrix_file_test.cc
  tl_dense_symmetric_matrix_file_test.cc

  tl_dense_general_matrix_arrays_coloriented_test.cc
  tl_dense_general_matrix_arrays_roworiented_test.cc
  tl_dense_general_matrix_arrays_mmap_roworiented_test.cc

  tl_dense_general_matrix_lapack_test.cc
  tl_dense_symmetric_matrix_lapack_test.cc

  #tl_sparse_matrix_test.cc
  #tl_sparse_symmetric_matrix_test.cc
  #tl_sparse_vector_test.cc

  TlPositionTest.cpp
  TlSerializeDataTest.cpp
  TlStlUtilsTest.cpp
  TlStringTokenizerTest.cpp
  TlUtilsTest.cpp

  dense_general_matrix_test_template.h
  test_general_matrix.cc

  dense_symmetric_matrix_test_template.h
  test_symmetric_matrix.cc

  dense_vector_test_template.h
  test_dense_vector.cc

  main.cpp
  )

set(PDF_XTEST_MPI_SOURCES
  TlCommunicateTest.cpp
  tl_dense_general_matrix_scalapack_test.cc
  tl_dense_symmetric_matrix_scalapack_test.cc
  tl_dense_vector_scalapack_test.cc
  main_MPI.cpp
  )


# objects using HDF5
if (HDF5_FOUND)
  set(PDF_XTEST_SOURCES
    ${PDF_XTEST_SOURCES}
    TlHdf5UtilsTest.cpp
    )
endif(HDF5_FOUND)


# objects using Lapack
if (LAPACK_FOUND)
  set(PDF_XTEST_MPI_SOURCES
    ${PDF_XTEST_MPI_SOURCES}
    tl_communicate_test_lapack.cc
  )
endif (LAPACK_FOUND)


# objects using Eigen3
if (EIGEN3_FOUND)
  set(PDF_XTEST_SOURCES
    ${PDF_XTEST_SOURCES}
    tl_dense_general_matrix_eigen_test.cc
    tl_dense_symmetric_matrix_eigen_test.cc
    #tl_dense_vector_eigen_test.cc

    tl_sparse_general_matrix_eigen_test.cc
    tl_sparse_symmetric_matrix_eigen_test.cc
    #tl_sparse_vector_eigen_test.cc
    )

  set(PDF_XTEST_MPI_SOURCES
    ${PDF_XTEST_MPI_SOURCES}
    tl_communicate_test_eigen.cc
  )
endif(EIGEN3_FOUND)


# objects using ViennaCL
if (EIGEN3_FOUND AND VIENNACL_FOUND)
  set(PDF_XTEST_SOURCES
    ${PDF_XTEST_SOURCES}
    tl_dense_general_matrix_viennacl_test.cc
    tl_dense_symmetric_matrix_viennacl_test.cc
    tl_dense_vector_viennacl_test.cc

    tl_sparse_general_matrix_viennacl_test.cc
    tl_sparse_symmetric_matrix_viennacl_test.cc
    #tl_sparse_vector_eigen_viennacl.cc
    )
endif(EIGEN3_FOUND AND VIENNACL_FOUND)


add_executable(pdf-xtest ${PDF_XTEST_SOURCES})
target_link_libraries(pdf-xtest libpdf libpdftl ${GTEST_LIBRARIES} ${GTEST_MAIN_LIBRARIES} Threads::Threads)

if (NOT DEFINED PDF_XTEST_FLAG)
  if (NOT DEFINED ENV{PDF_XTEST_FLAG})
      set(PDF_XTEST_FLAG "-d 0")
  else (NOT DEFINED ENV{PDF_XTEST_FLAG})
      set(PDF_XTEST_FLAG $ENV{PDF_XTEST_FLAG})
  endif (NOT DEFINED ENV{PDF_XTEST_FLAG})
endif (NOT DEFINED PDF_XTEST_FLAG)
separate_arguments(PDF_XTEST_FLAG)
message("PDF_XTEST_FLAG: ${PDF_XTEST_FLAG}")

add_test(
  NAME xtest
  COMMAND $<TARGET_FILE:pdf-xtest> ${PDF_XTEST_FLAG})
set_property(TEST xtest PROPERTY LABELS xtest)
set_property(TEST xtest PROPERTY RUN_SERIAL)
set_property(TEST xtest PROPERTY TIMEOUT 600)

if (MPI_FOUND)
  add_executable(pdf-xtest.MPI ${PDF_XTEST_MPI_SOURCES})
  target_link_libraries(pdf-xtest.MPI libppdf libppdftl libpdf libpdftl ${GTEST_LIBRARIES} ${GTEST_MAIN_LIBRARIES} Threads::Threads)

  # cf.) MPIEXEC_MAX_NUMPROCS
  if (NOT DEFINED MPI_PROCS)
    set(MPI_PROCS 4)
  endif (NOT DEFINED MPI_PROCS)

  add_test(
    NAME xtest.mpi
    COMMAND ${MPIEXEC}
      ${MPIEXEC_FLAGS}
      ${MPIEXEC_NUMPROC_FLAG}
      ${MPI_PROCS}
      ${MPIEXEC_PREFLAGS}
      $<TARGET_FILE:pdf-xtest.MPI>
      ${MPIEXEC_POSTFLAGS}
    )
    set_property(TEST xtest.mpi PROPERTY LABELS xtest.mpi)
    set_property(TEST xtest.mpi PROPERTY RUN_SERIAL)
    set_property(TEST xtest.mpi PROPERTY TIMEOUT 600)
endif(MPI_FOUND)
