set(PDF_XTEST_SOURCES
  common.h
  common.cpp
  DfFunctional_B3LYPTest.cpp
  DfFunctional_Becke88Test.cpp
  DfFunctional_LYPTest.cpp
  DfFunctional_SlaterTest.cpp
  DfFunctional_VWNTest.cpp
  TlColVectorMatrixTest.cpp
  TlFileMatrixTest.cpp
  TlFileSymmetricMatrixTest.cpp
  TlFileTest.cpp
  TlMatrixTest.cpp
  TlMemManagerTest.cpp
  TlMmapMatrix_CSFD_Test.cpp
  TlMmapMatrix_RSFD_Test.cpp
  TlPositionTest.cpp
  TlRowVectorMatrixTest.cpp
  TlSparseMatrixTest.cpp
  TlSparseSymmetricMatrixTest.cpp
  TlSparseVectorTest.cpp
  TlSerializeDataTest.cpp
  TlStlUtilsTest.cpp
  TlStringTokenizerTest.cpp
  TlSymmetricMatrixTest.cpp
  TlUtilsTest.cpp
  TlVectorTest.cpp
  main.cpp
  )

set(PDF_XTEST_MPI_SOURCES
  TlCommunicateTest.cpp
  TlDistributeMatrixTest.cpp
  TlDistributeSymmetricMatrixTest.cpp
  main_MPI.cpp
  )


if (HDF5_FOUND)
  set(PDF_XTEST_SOURCES
    ${PDF_XTEST_SOURCES}
    TlHdf5UtilsTest.cpp
    )
endif(HDF5_FOUND)


add_executable(pdf-xtest ${PDF_XTEST_SOURCES})
target_link_libraries(pdf-xtest libpdf libpdftl ${GTEST_LIBRARIES} ${GTEST_MAIN_LIBRARIES} Threads::Threads)

add_test(
  NAME xtest
  COMMAND $<TARGET_FILE:pdf-xtest>)

if (MPI_FOUND)
  add_executable(pdf-xtest.MPI ${PDF_XTEST_MPI_SOURCES})
  target_link_libraries(pdf-xtest.MPI libppdf libppdftl libpdf libpdftl ${GTEST_LIBRARIES} ${GTEST_MAIN_LIBRARIES} Threads::Threads)

  add_test(
    NAME xtest.mpi
    COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 2 $<TARGET_FILE:pdf-xtest.MPI>
    )
endif(MPI_FOUND)
