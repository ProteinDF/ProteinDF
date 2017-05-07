cmake_minimum_required(VERSION 3.5)

set(PDF_XTEST_SOURCES
  DfFunctional_B3LYPTest.cpp
  DfFunctional_Becke88Test.cpp
  DfFunctional_LYPTest.cpp
  DfFunctional_SlaterTest.cpp
  DfFunctional_VWNTest.cpp
  TlColVectorMatrixTest.cpp
  TlMatrixTest.cpp
  TlMemManagerTest.cpp
  #  TlMsgPackTest.cpp
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
  gtest/gtest-all.cc
  gtest/gtest.h
  )


add_executable(pdf-xtest ${PDF_XTEST_SOURCES})
target_link_libraries(pdf-xtest libpdf libpdftl Threads::Threads)


