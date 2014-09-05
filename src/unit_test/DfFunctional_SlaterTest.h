#ifndef DFFUNCTIONAL_SLATERTEST_H
#define DFFUNCTIONAL_SLATERTEST_H

#include <cppunit/extensions/HelperMacros.h>

#include "DfFunctional_Slater.h"

class DfFunctional_SlaterTest : public CppUnit::TestFixture{
  CPPUNIT_TEST_SUITE(DfFunctional_SlaterTest);
  CPPUNIT_TEST(testConstructer);
  CPPUNIT_TEST(testPointwise1);
  CPPUNIT_TEST(testPointwise1_RKS);
  CPPUNIT_TEST(testPointwise2);
  CPPUNIT_TEST(testPointwise2_RKS);
  CPPUNIT_TEST(testPointwise3);
  CPPUNIT_TEST(testPointwise4);
  CPPUNIT_TEST(testPointwise12);
  CPPUNIT_TEST(testPointwise13);
  CPPUNIT_TEST(testPointwise14);
  CPPUNIT_TEST(testPointwise15);
  CPPUNIT_TEST(testPointwise16);
  CPPUNIT_TEST(testPointwise17);
  CPPUNIT_TEST(testPointwise18);
  CPPUNIT_TEST(testPointwise19);
  CPPUNIT_TEST(testPointwise20);
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
  void testPointwise12();
  void testPointwise13();
  void testPointwise14();
  void testPointwise15();
  void testPointwise16();
  void testPointwise17();
  void testPointwise18();
  void testPointwise19();
  void testPointwise20();

public:
  DfFunctional_SlaterTest(){
  }

  void setUp(){
  }
  
  void tearDown(){
  }

private:
  const static double EPS;

};

#endif // DFFUNCTIONAL_SLATERTEST_H

