#ifndef DFFUNCTIONAL_B3LYPTEST_H
#define DFFUNCTIONAL_B3LYPTEST_H

#include <cppunit/extensions/HelperMacros.h>

#include "DfFunctional_B3LYP.h"

class DfFunctional_B3LYPTest : public CppUnit::TestFixture{
  CPPUNIT_TEST_SUITE(DfFunctional_B3LYPTest);
  CPPUNIT_TEST(testConstructer);
  CPPUNIT_TEST(testPointwise1);
  CPPUNIT_TEST(testPointwise1_RKS);
  CPPUNIT_TEST(testPointwise2);
  CPPUNIT_TEST(testPointwise2_RKS);
  CPPUNIT_TEST(testPointwise3);
  CPPUNIT_TEST(testPointwise3_RKS);
  CPPUNIT_TEST(testPointwise4);
  CPPUNIT_TEST(testPointwise4_RKS);
  CPPUNIT_TEST(testPointwise12);
  CPPUNIT_TEST(testPointwise13);
  CPPUNIT_TEST_SUITE_END();

public:
  void testConstructer();
  void testPointwise1();
  void testPointwise1_RKS();
  void testPointwise2();
  void testPointwise2_RKS();
  void testPointwise3();
  void testPointwise3_RKS();
  void testPointwise4();
  void testPointwise4_RKS();
  void testPointwise12();
  void testPointwise13();

public:
  DfFunctional_B3LYPTest(){
  }

  void setUp(){
  }
  
  void tearDown(){
  }

private:
  const static double EPS;

};

#endif // DFFUNCTIONAL_B3LYPTEST_H
