#ifndef TLPARAMETERTEST_H
#define TLPARAMETERTEST_H

#include <cppunit/extensions/HelperMacros.h>

#include "TlParameter.h"

class TlParameterTest : public CppUnit::TestFixture{
  CPPUNIT_TEST_SUITE(TlParameterTest);
  CPPUNIT_TEST(testConstructer);
  CPPUNIT_TEST(testCopyConstructer);
  CPPUNIT_TEST(testOperatorEqual);
  CPPUNIT_TEST_SUITE_END();

public:
  void testConstructer();
  void testCopyConstructer();
  void testOperatorEqual();

public:
  TlParameterTest(){
  }

  void setUp(){
  }
  
  void tearDown(){
  }
};

#endif // TLPARAMETERTEST_H

