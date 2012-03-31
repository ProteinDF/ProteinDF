#ifndef TLMATRIX_SYMMTERICTEST_H
#define TLMATRIX_SYMMETRICTEST_H

#include <cppunit/extensions/HelperMacros.h>

#include "TlSymmetricMatrix.h"
#include "TlVector.h"

class TlSymmetricMatrixTest : public CppUnit::TestFixture{
    CPPUNIT_TEST_SUITE(TlSymmetricMatrixTest);
    CPPUNIT_TEST(testConstructer);
    CPPUNIT_TEST(testOperaterRoundBracket);
    CPPUNIT_TEST(testCopyConstructer);
    CPPUNIT_TEST(testConvertFromTlVector1);
    CPPUNIT_TEST(testConvertFromTlVector2);
    CPPUNIT_TEST(testOperatorEqual);
    CPPUNIT_TEST(testOperatorPlus);
    CPPUNIT_TEST(testOperatorPlusEqual);
    CPPUNIT_TEST(testOperatorAsterisk);
    CPPUNIT_TEST(testSave);
    CPPUNIT_TEST(testLoad);
    CPPUNIT_TEST(testInverse);
    CPPUNIT_TEST(testMulti1);
    CPPUNIT_TEST(testMulti2);
    CPPUNIT_TEST(testDot);
    CPPUNIT_TEST(testSum);
    CPPUNIT_TEST(testCholeskyDecomposition);
    CPPUNIT_TEST_SUITE_END();

public:
    void testConstructer();
    void testOperaterRoundBracket();
    void testCopyConstructer();
    void testConvertFromTlVector1();
    void testConvertFromTlVector2();
    void testOperatorEqual();
    void testOperatorPlus();
    void testOperatorPlusEqual();
    void testOperatorAsterisk();
    void testSave();
    void testLoad();
    void testInverse();
    void testMulti1();
    void testMulti2();
    void testDot();
    void testSum();
    void testCholeskyDecomposition();

public:
    TlSymmetricMatrixTest(){
    }

    void setUp(){
    }
  
    void tearDown(){
    }

private:
    TlSymmetricMatrix getMatrixA();
    TlSymmetricMatrix getMatrixB();
    TlSymmetricMatrix getMatrixC();

private:
    static const double threshold;
};

//CPPUNIT_TEST_SUITE_REGISTRATION(TlSymmetricMatrixTest);

#endif // TLMATRIX_SYMMETRICTEST_H

