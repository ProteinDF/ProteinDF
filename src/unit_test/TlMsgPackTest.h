#ifndef TLMSGPACKTEST_H
#define TLMSGPACKTEST_H

#include <cppunit/extensions/HelperMacros.h>

#include "TlMsgPack.h"

class TlMsgPackTest : public CppUnit::TestFixture {
  CPPUNIT_TEST_SUITE(TlMsgPackTest);
  CPPUNIT_TEST(testLoad);
  CPPUNIT_TEST(testDumpAndPack);
  CPPUNIT_TEST_SUITE_END();

 public:
  void testLoad();
  void testDumpAndPack();

 public:
  TlMsgPackTest() {}

  void setUp() {}

  void tearDown() {}
};

#endif  // TLMSGPACKTEST_H
