#ifndef DFFUNCTIONAL_PW91XTEST_H
#define DFFUNCTIONAL_PW91XTEST_H

#include <cppunit/extensions/HelperMacros.h>

#include "DfFunctional_PW91X.h"

class DfFunctional_PW91XTest : public CppUnit::TestFixture {
  CPPUNIT_TEST_SUITE(DfFunctional_PW91XTest);
  CPPUNIT_TEST(testConstructer);
  CPPUNIT_TEST(testPointwise1);
  CPPUNIT_TEST(testPointwise1_RKS);
  //   CPPUNIT_TEST(testPointwise2);
  //   CPPUNIT_TEST(testPointwise2_RKS);
  //   CPPUNIT_TEST(testPointwise3);
  //   CPPUNIT_TEST(testPointwise3_RKS);
  //   CPPUNIT_TEST(testPointwise4);
  //   CPPUNIT_TEST(testPointwise11);
  //   CPPUNIT_TEST(testPointwise12);
  CPPUNIT_TEST_SUITE_END();

 public:
  void testConstructer();
  void testPointwise1();
  void testPointwise1_RKS();
  //   void testPointwise2();
  //   void testPointwise2_RKS();
  //   void testPointwise3();
  //   void testPointwise3_RKS();
  //   void testPointwise4();
  //   void testPointwise11();
  //   void testPointwise12();

 public:
  DfFunctional_PW91XTest() {}

  void setUp() {}

  void tearDown() {}

 private:
  const static double EPS;
};

#endif  // DFFUNCTIONAL_PW91XTEST_H
