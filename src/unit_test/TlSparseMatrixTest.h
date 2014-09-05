#ifndef TLSPARSEMATRIXTEST_H
#define TLSPARSEMATRIXTEST_H

#include <cppunit/extensions/HelperMacros.h>

#include "TlSparseMatrix.h"

class TlSparseMatrixTest : public CppUnit::TestFixture{
  CPPUNIT_TEST_SUITE(TlSparseMatrixTest);
  CPPUNIT_TEST(testConstructer);
  CPPUNIT_TEST(testCopyConstructer);
  CPPUNIT_TEST(testMerge);
  CPPUNIT_TEST_SUITE_END();

public:
  void testConstructer();
  void testCopyConstructer();
  void testMerge();

public:
  TlSparseMatrixTest(){
  }

  void setUp(){
  }
  
  void tearDown(){
  }

private:
  static const double threshold;
};

//CPPUNIT_TEST_SUITE_REGISTRATION(TlSparseMatrixTest);

#endif // TLSPARSEMATRIXTEST_H

