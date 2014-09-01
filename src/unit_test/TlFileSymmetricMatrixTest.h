#ifndef TLFILESYMMETRICMATRIXTEST_H
#define TLFILESYMMETRICMATRIXTEST_H

#include <cppunit/extensions/HelperMacros.h>

#include "TlFileSymmetricMatrix.h"

class TlFileSymmetricMatrixTest : public CppUnit::TestFixture{
  CPPUNIT_TEST_SUITE(TlFileSymmetricMatrixTest);
  CPPUNIT_TEST(testConstructer);
  CPPUNIT_TEST(testGet);
  CPPUNIT_TEST(testSet);
  CPPUNIT_TEST(testAdd);
  CPPUNIT_TEST_SUITE_END();

public:
  void testConstructer();
  void testGet();
  void testSet();
  void testAdd();

public:
  TlFileSymmetricMatrixTest(){
  }

  void setUp(){
  }
  
  void tearDown(){
  }

private:
  static const double threshold;
};

#endif // TLFILESYMMETRICMATRIXTEST_H

