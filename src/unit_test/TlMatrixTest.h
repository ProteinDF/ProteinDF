#ifndef TLMATRIXTEST_H
#define TLMATRIXTEST_H

#include <cppunit/extensions/HelperMacros.h>

#include "TlMatrix.h"
#include "TlVector.h"

class TlMatrixTest : public CppUnit::TestFixture {
    CPPUNIT_TEST_SUITE(TlMatrixTest);
    CPPUNIT_TEST(testConstructer);
    //CPPUNIT_TEST(testConvertFromVector);
    CPPUNIT_TEST(testCopyConstructer);
    CPPUNIT_TEST(testOperatorEqual);
    CPPUNIT_TEST(testOperatorPlus);
    CPPUNIT_TEST(testOperatorPlusEqual);
    CPPUNIT_TEST(testSave);
    CPPUNIT_TEST(testLoad);
    CPPUNIT_TEST(testInverse);
    CPPUNIT_TEST(testOperatorMul_AB);
    CPPUNIT_TEST(testOperatorMul_AX);
    CPPUNIT_TEST(testOperatorMul_XA);
    CPPUNIT_TEST(testDot);
    CPPUNIT_TEST(testSum);
    CPPUNIT_TEST_SUITE_END();
    
public:
    void testConstructer();
    void testConvertFromVector();
    void testCopyConstructer();
    void testOperatorEqual();
    void testOperatorPlus();
    void testOperatorPlusEqual();
    void testSave();
    void testLoad();
    void testInverse();
    void testOperatorMul_AB();
    void testOperatorMul_AX();
    void testOperatorMul_XA();
    void testDot();
    void testSum();
    
public:
    TlMatrixTest(){
    }
    
    void setUp(){
    }
    
    void tearDown(){
    }
    
private:
    TlMatrix getMatrixA();
    TlMatrix getMatrixB();
    TlMatrix getMatrixC();
    
private:
    static const double threshold;
};

//CPPUNIT_TEST_SUITE_REGISTRATION(TlMatrixTest);

#endif // TLMATRIXTEST_H

