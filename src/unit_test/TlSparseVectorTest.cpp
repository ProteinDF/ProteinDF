#include <limits>
#include "TlSparseVectorTest.h"

const double TlSparseVectorTest::threshold = std::numeric_limits<double>::epsilon();

void TlSparseVectorTest::testConstructer()
{
    TlSparseVector a(100);

    for (int i = 0; i < 100; ++i) {
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, a[i], threshold);
    }
}


void TlSparseVectorTest::testAdd_double() {
    TlSparseVector a(100);
    
    a[10] = 10.0;
    a[10] += 20.0;
    a[51] += 3.0;

    CPPUNIT_ASSERT_DOUBLES_EQUAL(30.0, a[10], threshold);
    CPPUNIT_ASSERT_DOUBLES_EQUAL( 3.0, a[51], threshold);
    CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, a[80], threshold);
}

void TlSparseVectorTest::testAdd_TlSparseVector() {
    TlSparseVector a(100);
    TlSparseVector b(100);

    a[ 2] = 2.0;
    a[ 5] = 5.0;
    a[31] = 31.0;
    b[17] = 17.0;
    b[31] = 20.0;

    a += b;

    CPPUNIT_ASSERT_DOUBLES_EQUAL( 2.0, a[ 2], threshold);
    CPPUNIT_ASSERT_DOUBLES_EQUAL( 5.0, a[ 5], threshold);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(17.0, a[17], threshold);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(51.0, a[31], threshold);
    CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, a[99], threshold);
}

