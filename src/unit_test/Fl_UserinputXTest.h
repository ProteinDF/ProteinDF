#ifndef FL_USERINPUTXTEST_H
#define FL_USERINPUTXTEST_H

#include <cppunit/extensions/HelperMacros.h>

#include "Fl_UserinputX.h"

class Fl_UserinputXTest : public CppUnit::TestFixture{
  CPPUNIT_TEST_SUITE(Fl_UserinputXTest);
  CPPUNIT_TEST(testConstructer);
  CPPUNIT_TEST(testGetFlGlobalinputX);
  CPPUNIT_TEST_SUITE_END();

public:
  void testConstructer();
  void testGetFlGlobalinputX();

public:
  Fl_UserinputXTest(){
  }

  void setUp(){
  }
  
  void tearDown(){
  }
};

//CPPUNIT_TEST_SUITE_REGISTRATION(TlMatrixTest);

#endif // TLMATRIXTEST_H

