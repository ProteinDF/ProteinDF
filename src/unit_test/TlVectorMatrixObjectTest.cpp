#include <limits>
#include <iostream>
#include "TlVectorMatrixObjectTest.h"
#include "TlMatrix.h"

// CPPUNIT_ASSERT( condition );
// conditionが偽(false,0)であったとき、失敗します。
//
// CPPUNIT_ASSERT_MESSAGE( message, condition );
// conditionが偽であったとき、失敗します。このときmessageを出力します。
//
// CPPUNIT_FAIL( message );
// 必ず失敗します。messageを出力します。
//
// CPPUNIT_ASSERT_EQUAL( expected, actual );
// 得られた結果actualが期待する値expectedでなかったとき、すなわちexpected != actualのときに失敗します。
//
// CPPUNIT_ASSERT_EQUAL_MESSAGE( message, expected, actual );
// 得られた結果actualが期待する値expectedでなかったとき、すなわちexpected != actualのときに失敗します。このときmessageを出力します。
//
// CPPUNIT_ASSERT_DOUBLES_EQUAL( expected, actual, delta );
// 得られた結果actualと期待する値expectedとの差がdeltaより大きいとき、失敗します。

const double TlVectorMatrixObjectTest::threshold = std::numeric_limits<double>::epsilon();


void TlVectorMatrixObjectTest::testConstructer()
{
    TlVectorMatrixObject A(100, 600);
    
    CPPUNIT_ASSERT_EQUAL(100, A.getNumOfVectors());
    CPPUNIT_ASSERT_EQUAL(600, A.getSizeOfVector());
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, A.get(0, 0), threshold);
}


void TlVectorMatrixObjectTest::testConstructer2()
{
    TlVectorMatrixObject A0(100, 600, 10, 0);
    TlVectorMatrixObject A1(100, 600, 10, 1);
    
    CPPUNIT_ASSERT_EQUAL(100, A0.getNumOfVectors());
    CPPUNIT_ASSERT_EQUAL(600, A0.getSizeOfVector());
    CPPUNIT_ASSERT_EQUAL( 10, A0.getNumOfSubunits());
    CPPUNIT_ASSERT_EQUAL(  0, A0.getSubunitID());

    CPPUNIT_ASSERT_EQUAL(100, A1.getNumOfVectors());
    CPPUNIT_ASSERT_EQUAL(600, A1.getSizeOfVector());
    CPPUNIT_ASSERT_EQUAL( 10, A1.getNumOfSubunits());
    CPPUNIT_ASSERT_EQUAL(  1, A1.getSubunitID());
}


void TlVectorMatrixObjectTest::testResize()
{
    TlVectorMatrixObject A(100, 100);
    A.resize(200, 100);
    CPPUNIT_ASSERT_EQUAL(200, A.getNumOfVectors());
    CPPUNIT_ASSERT_EQUAL(100, A.getSizeOfVector());

    TlVectorMatrixObject B(100, 100);
    B.resize(100, 200);
    CPPUNIT_ASSERT_EQUAL(100, B.getNumOfVectors());
    CPPUNIT_ASSERT_EQUAL(200, B.getSizeOfVector());

    TlVectorMatrixObject C(50, 100);
    C.resize(50, 100);
    CPPUNIT_ASSERT_EQUAL( 50, C.getNumOfVectors());
    CPPUNIT_ASSERT_EQUAL(100, C.getSizeOfVector());

    TlVectorMatrixObject D(100, 50);
    D.resize(100, 50);
    CPPUNIT_ASSERT_EQUAL(100, D.getNumOfVectors());
    CPPUNIT_ASSERT_EQUAL( 50, D.getSizeOfVector());
}


void TlVectorMatrixObjectTest::testContents()
{
    const int numOfVector = 80;
    const int sizeOfVector = 100;
    TlVectorMatrixObject vecA(numOfVector, sizeOfVector);
    TlMatrix matA(numOfVector, sizeOfVector);

    // setup
    int count = 0;
    for (int i = 0; i < numOfVector; ++i) {
        for (int j = 0; j < sizeOfVector; ++j) {
            double v = double(count);
            matA.set(i, j, v);
            vecA.set(i, j, v);
            
            ++count;
        }
    }

    // test
    for (int i = 0; i < numOfVector; ++i) {
        for (int j = 0; j < sizeOfVector; ++j) {
            CPPUNIT_ASSERT_DOUBLES_EQUAL(matA.get(i, j),
                                         vecA.get(i, j),
                                         TlVectorMatrixObjectTest::threshold);
        }
    }
}


