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

  tl_dense_general_matrix_colvector_test.cc
  tl_dense_general_matrix_rowvector_test.cc

  tl_dense_general_matrix_lapack_test.cc
  tl_dense_symmetric_matrix_lapack_test.cc

  tl_sparse_matrix_test.cc
  tl_sparse_symmetric_matrix_test.cc
  tl_sparse_vector_test.cc

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
  main_MPI.cpp
  )


# objects using HDF5
if (HDF5_FOUND)
  set(PDF_XTEST_SOURCES
    ${PDF_XTEST_SOURCES}
    TlHdf5UtilsTest.cpp
    )
endif(HDF5_FOUND)


# objects using Eigen3
if (EIGEN3_FOUND)
  set(PDF_XTEST_SOURCES
    ${PDF_XTEST_SOURCES}
    tl_dense_general_matrix_eigen_test.cc
    tl_dense_symmetric_matrix_eigen_test.cc
    )
endif(EIGEN3_FOUND)


# objects using ViennaCL
if (VIENNACL_FOUND)
  set(PDF_XTEST_SOURCES
    ${PDF_XTEST_SOURCES}
    tl_dense_general_matrix_viennacl_test.cc
    tl_dense_symmetric_matrix_viennacl_test.cc
    )
endif(VIENNACL_FOUND)


add_executable(pdf-xtest ${PDF_XTEST_SOURCES})
target_link_libraries(pdf-xtest libpdf libpdftl ${OpenCL_LIBRARIES} ${GTEST_LIBRARIES} ${GTEST_MAIN_LIBRARIES} Threads::Threads)

add_test(
  NAME xtest
  COMMAND $<TARGET_FILE:pdf-xtest>)
set_property(TEST xtest PROPERTY LABELS xtest)
set_property(TEST xtest PROPERTY RUN_SERIAL)
set_property(TEST xtest PROPERTY TIMEOUT 120)

if (MPI_FOUND)
  add_executable(pdf-xtest.MPI ${PDF_XTEST_MPI_SOURCES})
  target_link_libraries(pdf-xtest.MPI libppdf libppdftl libpdf libpdftl ${OpenCL_LIBRARIES} ${GTEST_LIBRARIES} ${GTEST_MAIN_LIBRARIES} Threads::Threads)

  add_test(
    NAME xtest.mpi
    COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 $<TARGET_FILE:pdf-xtest.MPI>
    )
    set_property(TEST xtest.mpi PROPERTY LABELS xtest.mpi)
    set_property(TEST xtest.mpi PROPERTY RUN_SERIAL)
    set_property(TEST xtest.mpi PROPERTY TIMEOUT 120)
endif(MPI_FOUND)
