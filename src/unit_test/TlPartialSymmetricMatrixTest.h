#ifndef TLPARTIALSYMMETRICMATRIXTEST_H
#define TLPARTIALSYMMETRICMATRIXTEST_H

#include <cppunit/extensions/HelperMacros.h>

#include "TlPartialSymmetricMatrix.h"

class TlPartialSymmetricMatrixTest : public CppUnit::TestFixture {
    CPPUNIT_TEST_SUITE(TlPartialSymmetricMatrixTest);
    CPPUNIT_TEST(testConstructer);
    CPPUNIT_TEST(testSetGet);
    CPPUNIT_TEST(testAdd);
    CPPUNIT_TEST(testSetGet2);
    CPPUNIT_TEST(testAdd2);
    CPPUNIT_TEST(testSetGet3);
    CPPUNIT_TEST(testGetRowVector);
    CPPUNIT_TEST_SUITE_END();
    
public:
    void testConstructer();
    void testSetGet();
    void testAdd();
    void testSetGet2();
    void testAdd2();
    void testSetGet3();
    void testGetRowVector();

public:
  TlPartialSymmetricMatrixTest(){
  }

  void setUp(){
  }
  
  void tearDown(){
  }

private:
  static const double threshold;
};

#endif // TLPARTIALSYMMETRICMATRIXTEST_H

