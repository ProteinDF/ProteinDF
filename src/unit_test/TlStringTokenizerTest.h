#ifndef TLSTRINGTOKENIZERTEST_H
#define TLSTRINGTOKENIZERTEST_H

#include <cppunit/extensions/HelperMacros.h>

#include "TlStringTokenizer.h"

class TlStringTokenizerTest : public CppUnit::TestFixture{
  CPPUNIT_TEST_SUITE(TlStringTokenizerTest);
  CPPUNIT_TEST(testConstructer);
  CPPUNIT_TEST(testCountTokens);
  CPPUNIT_TEST(testGetTokens);
  CPPUNIT_TEST_SUITE_END();

public:
  void testConstructer();
  void testCountTokens();
  void testGetTokens();

public:
  TlStringTokenizerTest(){
  }

  void setUp(){
  }
  
  void tearDown(){
  }
};

#endif // TLSTRINGTOKENIZERTEST_H

