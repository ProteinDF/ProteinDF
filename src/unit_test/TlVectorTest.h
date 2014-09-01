#ifndef TLVECTORTEST_H
#define TLVECTORTEST_H

#include <cppunit/extensions/HelperMacros.h>

#include "TlVector.h"
#include "TlMatrix.h"
#include "TlSymmetricMatrix.h"

class TlVectorTest : public CppUnit::TestFixture{
  CPPUNIT_TEST_SUITE(TlVectorTest);
  CPPUNIT_TEST(testConstructer);
  CPPUNIT_TEST(testCopyConstructer);
  CPPUNIT_TEST(testOperatorCopy);
  CPPUNIT_TEST(testResize);
//   CPPUNIT_TEST(testConvertFromMatrix);
//   CPPUNIT_TEST(testConvertFromMatrix_Symmetric);
  CPPUNIT_TEST(testOperatorPlus);
  CPPUNIT_TEST(testOperatorPlusEqual);
  CPPUNIT_TEST(testOperatorMinus);
  CPPUNIT_TEST(testOperatorMinusEqual);
  CPPUNIT_TEST(testOperatorMultipleDouble);
  CPPUNIT_TEST(testOperatorMultipleEqualDouble);
  CPPUNIT_TEST(testOperatorMultipleVector);
  CPPUNIT_TEST(testGetMaxAbsoluteElement);
  CPPUNIT_TEST(testAdd);
  CPPUNIT_TEST(testSave);
  CPPUNIT_TEST(testLoad);
  CPPUNIT_TEST_SUITE_END();

public:
  void testConstructer();
  void testCopyConstructer();
  void testOperatorCopy();
  void testResize();
//   void testConvertFromMatrix();
//   void testConvertFromMatrix_Symmetric();
  void testOperatorPlus();
  void testOperatorPlusEqual();
  void testOperatorMinus();
  void testOperatorMinusEqual();
  void testOperatorMultipleDouble();
  void testOperatorMultipleEqualDouble();
  void testOperatorMultipleVector();
  void testGetMaxAbsoluteElement();
  void testAdd();
  void testSave();
  void testLoad();

public:
  TlVectorTest(){
  }

  void setUp(){
  }
  
  void tearDown(){
  }

private:
  TlVector getVectorA();
  TlVector getVectorB();

private:
  static const double threshold;
};

//CPPUNIT_TEST_SUITE_REGISTRATION(TlVectorTest);

#endif // TLVECTORTEST_H

