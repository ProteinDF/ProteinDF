#include <limits>
#include <iostream>
#include "TlColVectorMatrixTest.h"
#include "TlMatrix.h"

const double TlColVectorMatrixTest::threshold = std::numeric_limits<double>::epsilon();


void TlColVectorMatrixTest::testConstructer()
{
    TlColVectorMatrix A(100, 600);
    
    CPPUNIT_ASSERT_EQUAL(100, A.getSizeOfVector());
    CPPUNIT_ASSERT_EQUAL(600, A.getNumOfVectors());
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, A.get(0, 0), threshold);
}


void TlColVectorMatrixTest::testConstructer2()
{
    TlColVectorMatrix A0(100, 600, 10, 0);
    TlColVectorMatrix A1(100, 600, 10, 1);
    
    CPPUNIT_ASSERT_EQUAL(100, A0.getSizeOfVector());
    CPPUNIT_ASSERT_EQUAL(600, A0.getNumOfVectors());
    CPPUNIT_ASSERT_EQUAL( 10, A0.getNumOfSubunits());
    CPPUNIT_ASSERT_EQUAL(  0, A0.getSubunitID());

    CPPUNIT_ASSERT_EQUAL(100, A1.getSizeOfVector());
    CPPUNIT_ASSERT_EQUAL(600, A1.getNumOfVectors());
    CPPUNIT_ASSERT_EQUAL( 10, A1.getNumOfSubunits());
    CPPUNIT_ASSERT_EQUAL(  1, A1.getSubunitID());
}


void TlColVectorMatrixTest::testResize()
{
    TlColVectorMatrix A(100, 100);
    A.resize(200, 100);
    CPPUNIT_ASSERT_EQUAL(200, A.getSizeOfVector());
    CPPUNIT_ASSERT_EQUAL(100, A.getNumOfVectors());

    TlColVectorMatrix B(100, 100);
    B.resize(100, 200);
    CPPUNIT_ASSERT_EQUAL(100, B.getSizeOfVector());
    CPPUNIT_ASSERT_EQUAL(200, B.getNumOfVectors());

    TlColVectorMatrix C(50, 100);
    C.resize(50, 100);
    CPPUNIT_ASSERT_EQUAL( 50, C.getSizeOfVector());
    CPPUNIT_ASSERT_EQUAL(100, C.getNumOfVectors());

    TlColVectorMatrix D(100, 50);
    D.resize(100, 50);
    CPPUNIT_ASSERT_EQUAL(100, D.getSizeOfVector());
    CPPUNIT_ASSERT_EQUAL( 50, D.getNumOfVectors());
}


void TlColVectorMatrixTest::testContents()
{
    const int maxRow = 100;
    const int maxCol = 80;
    TlColVectorMatrix vecA(maxRow, maxCol);
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
                                         TlColVectorMatrixTest::threshold);
        }
    }
}


void TlColVectorMatrixTest::testSaveLoad()
{
    const int maxRow = 100;
    const int maxCol = 80;
    TlColVectorMatrix vecA(maxRow, maxCol);
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
    TlColVectorMatrix vecB;
    vecB.load("/tmp/vecA.mat", 0);

    // test
    for (int r = 0; r < maxRow; ++r) {
        for (int c = 0; c < maxCol; ++c) {
            CPPUNIT_ASSERT_DOUBLES_EQUAL(matA.get(r, c),
                                         vecB.get(r, c),
                                         TlColVectorMatrixTest::threshold);
        }
    }
}


void TlColVectorMatrixTest::testToTlMatrix()
{
    const int maxRow = 100;
    const int maxCol = 80;
    TlColVectorMatrix vecA(maxRow, maxCol);
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
                                         TlColVectorMatrixTest::threshold);
        }
    }
}

