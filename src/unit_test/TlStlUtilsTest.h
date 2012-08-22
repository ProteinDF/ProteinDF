#ifndef TLSTLUTILSTEST_H
#define TLSTLUTILSTEST_H

#include <cppunit/extensions/HelperMacros.h>

#include "TlStlUtils.h"

class TlStlUtilsTest : public CppUnit::TestFixture {
  CPPUNIT_TEST_SUITE(TlStlUtilsTest);
  CPPUNIT_TEST(testCache);
  CPPUNIT_TEST_SUITE_END();

public:
  void testCache();

public:
  TlStlUtilsTest() {
  }

  void setUp() {
  }
  
  void tearDown() {
  }
};

#endif // TLSTLUTILESTEST_H

