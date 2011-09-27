#ifndef TLPOSITIONTEST_H
#define TLPOSITIONTEST_H

#include <cppunit/extensions/HelperMacros.h>

#include "TlPosition.h"

class TlPositionTest : public CppUnit::TestFixture{
  CPPUNIT_TEST_SUITE(TlPositionTest);
  CPPUNIT_TEST(testConstructer);
  CPPUNIT_TEST(testCopyConstructer);
  CPPUNIT_TEST(testOperatorEqual);
  CPPUNIT_TEST_SUITE_END();

public:
  void testConstructer();
  void testCopyConstructer();
  void testOperatorEqual();
  void testSquareDistanceFrom();

public:
  TlPositionTest(){
  }

  void setUp(){
  }
  
  void tearDown(){
  }

private:
  static const double threshold;
};

#endif // TLPARAMETERTEST_H

