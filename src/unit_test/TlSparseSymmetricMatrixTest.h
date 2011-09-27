#ifndef TLSPARSESYMMETRICMATRIXTEST_H
#define TLSPARSESYMMETRICMATRIXTEST_H

#include <cppunit/extensions/HelperMacros.h>

#include "TlSparseSymmetricMatrix.h"

class TlSparseSymmetricMatrixTest : public CppUnit::TestFixture{
  CPPUNIT_TEST_SUITE(TlSparseSymmetricMatrixTest);
  CPPUNIT_TEST(testConstructer);
  CPPUNIT_TEST(testMerge);
  CPPUNIT_TEST_SUITE_END();

public:
  void testConstructer();
  void testMerge();

public:
  TlSparseSymmetricMatrixTest(){
  }

  void setUp(){
  }
  
  void tearDown(){
  }

private:
  static const double threshold;
};

#endif // TLSPARSESYMMETRICMATRIXTEST_H

