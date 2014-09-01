#include <limits>
#include <iostream>
#include "TlRowVectorMatrixTest.h"
#include "TlMatrix.h"

const double TlRowVectorMatrixTest::threshold = std::numeric_limits<double>::epsilon();


void TlRowVectorMatrixTest::testConstructer()
{
    TlRowVectorMatrix A(100, 600);
    
    CPPUNIT_ASSERT_EQUAL(100, A.getNumOfRows());
    CPPUNIT_ASSERT_EQUAL(600, A.getNumOfCols());
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, A.get(0, 0), threshold);
}


void TlRowVectorMatrixTest::testConstructer2()
{
    TlRowVectorMatrix A0(100, 600, 10, 0);
    TlRowVectorMatrix A1(100, 600, 10, 1);
    
    CPPUNIT_ASSERT_EQUAL(100, A0.getNumOfRows());
    CPPUNIT_ASSERT_EQUAL(600, A0.getNumOfCols());
    CPPUNIT_ASSERT_EQUAL( 10, A0.getNumOfSubunits());
    CPPUNIT_ASSERT_EQUAL(  0, A0.getSubunitID());

    CPPUNIT_ASSERT_EQUAL(100, A1.getNumOfRows());
    CPPUNIT_ASSERT_EQUAL(600, A1.getNumOfCols());
    CPPUNIT_ASSERT_EQUAL( 10, A1.getNumOfSubunits());
    CPPUNIT_ASSERT_EQUAL(  1, A1.getSubunitID());
}


void TlRowVectorMatrixTest::testResize()
{
    TlRowVectorMatrix A(100, 100);
    A.resize(200, 100);
    CPPUNIT_ASSERT_EQUAL(200, A.getNumOfRows());
    CPPUNIT_ASSERT_EQUAL(100, A.getNumOfCols());

    TlRowVectorMatrix B(100, 100);
    B.resize(100, 200);
    CPPUNIT_ASSERT_EQUAL(100, B.getNumOfRows());
    CPPUNIT_ASSERT_EQUAL(200, B.getNumOfCols());

    TlRowVectorMatrix C(50, 100);
    C.resize(50, 100);
    CPPUNIT_ASSERT_EQUAL( 50, C.getNumOfRows());
    CPPUNIT_ASSERT_EQUAL(100, C.getNumOfCols());

    TlRowVectorMatrix D(100, 50);
    D.resize(100, 50);
    CPPUNIT_ASSERT_EQUAL(100, D.getNumOfRows());
    CPPUNIT_ASSERT_EQUAL( 50, D.getNumOfCols());
}


void TlRowVectorMatrixTest::testContents()
{
    const int maxRow = 100;
    const int maxCol = 80;
    TlRowVectorMatrix vecA(maxRow, maxCol);
    TlMatrix matA(maxRow, maxCol);

    // setup
    int count = 0;
    for (int r = 0; r < maxRow; ++r) {
        for (int c = 0; c < maxCol; ++c) {
            double v = double(count);
            matA.set(r, c, v);
            vecA.set(r, c, v);
            
            ++count;
        }
    }

    // test
    for (int r = 0; r < maxRow; ++r) {
        for (int c = 0; c < maxCol; ++c) {
            CPPUNIT_ASSERT_DOUBLES_EQUAL(matA.get(r, c),
                                         vecA.get(r, c),
                                         TlRowVectorMatrixTest::threshold);
        }
    }
}


void TlRowVectorMatrixTest::testSaveLoad()
{
    const int maxRow = 100;
    const int maxCol = 80;
    TlRowVectorMatrix vecA(maxRow, maxCol);
    TlMatrix matA(maxRow, maxCol);

    // setup
    int count = 0;
    for (int r = 0; r < maxRow; ++r) {
        for (int c = 0; c < maxCol; ++c) {
            double v = double(count);
            matA.set(r, c, v);
            vecA.set(r, c, v);
            
            ++count;
        }
    }

    // save
    vecA.save("/tmp/vecA.mat");

    // load
    TlRowVectorMatrix vecB;
    vecB.load("/tmp/vecA.mat");

    // test
    for (int r = 0; r < maxRow; ++r) {
        for (int c = 0; c < maxCol; ++c) {
            CPPUNIT_ASSERT_DOUBLES_EQUAL(matA.get(r, c),
                                         vecB.get(r, c),
                                         TlRowVectorMatrixTest::threshold);
        }
    }
}


void TlRowVectorMatrixTest::testToTlMatrix()
{
    const int maxRow = 100;
    const int maxCol = 80;
    TlRowVectorMatrix vecA(maxRow, maxCol);
    TlMatrix matA(maxRow, maxCol);

    // setup
    int count = 0;
    for (int r = 0; r < maxRow; ++r) {
        for (int c = 0; c < maxCol; ++c) {
            double v = double(count);
            matA.set(r, c, v);
            vecA.set(r, c, v);
            
            ++count;
        }
    }

    TlMatrix matB = vecA.getTlMatrixObject();

    // test
    for (int r = 0; r < maxRow; ++r) {
        for (int c = 0; c < maxCol; ++c) {
            CPPUNIT_ASSERT_DOUBLES_EQUAL(matA.get(r, c),
                                         matB.get(r, c),
                                         TlRowVectorMatrixTest::threshold);
        }
    }
}
