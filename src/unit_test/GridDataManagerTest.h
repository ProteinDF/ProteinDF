#ifndef GRIDDATAMANAGERTEST_H
#define GRIDDATAMANAGERTEST_H

#include <cppunit/extensions/HelperMacros.h>

#include "GridDataManager.h"

class GridDataManagerTest : public CppUnit::TestFixture {
  CPPUNIT_TEST_SUITE(GridDataManagerTest);
  CPPUNIT_TEST(testConstructer);
  CPPUNIT_TEST(testSetDensity);
  CPPUNIT_TEST(testGetDensity);
  CPPUNIT_TEST(testUpdate);
  CPPUNIT_TEST_SUITE_END();

 public:
  void testConstructer();
  void testSetDensity();
  void testGetDensity();
  void testUpdate();

 public:
  GridDataManagerTest() {}

  void setUp() {}

  void tearDown() {}
};

#endif  // GRIDDATAMANAGERTEST_H
