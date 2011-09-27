#ifndef TLUTILSTEST_H
#define TLUTILSTEST_H

#include <cppunit/extensions/HelperMacros.h>

#include "TlUtils.h"

class TlUtilsTest : public CppUnit::TestFixture{
  CPPUNIT_TEST_SUITE(TlUtilsTest);
  CPPUNIT_TEST(testPad);
  CPPUNIT_TEST(testTrim);
  CPPUNIT_TEST(testTrim_ws);
  CPPUNIT_TEST_SUITE_END();

public:
  void testPad();
  void testTrim();
  void testTrim_ws();

public:
  TlUtilsTest(){
  }

  void setUp(){
  }
  
  void tearDown(){
  }
};

#endif // TLUTILESTEST_H

