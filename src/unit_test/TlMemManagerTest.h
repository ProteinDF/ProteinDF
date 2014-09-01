#ifndef TLMEMMANAGERTEST_H
#define TLMEMMANAGERTEST_H

#include <cppunit/extensions/HelperMacros.h>

#include "TlMemManager.h"

class TlMemManagerTest : public CppUnit::TestFixture {
  CPPUNIT_TEST_SUITE(TlMemManagerTest);
  CPPUNIT_TEST(testAllocate);
  CPPUNIT_TEST_SUITE_END();

public:
  void testAllocate();

public:
  TlMemManagerTest(){
  }

  void setUp(){
  }
  
  void tearDown(){
  }
};

#endif // TLMEMMANAGERTEST_H

