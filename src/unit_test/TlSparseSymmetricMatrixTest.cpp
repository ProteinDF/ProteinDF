#include <limits>
#include "TlSparseSymmetricMatrixTest.h"

const double TlSparseSymmetricMatrixTest::threshold = std::numeric_limits<double>::epsilon();

void TlSparseSymmetricMatrixTest::testConstructer(){
  TlSparseSymmetricMatrix a(3);

  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, a(0, 0), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, a(0, 1), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, a(0, 2), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, a(1, 0), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, a(1, 1), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, a(1, 2), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, a(2, 0), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, a(2, 1), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, a(2, 2), threshold);
}

void TlSparseSymmetricMatrixTest::testMerge(){
  TlSparseSymmetricMatrix a(5);
  TlSparseSymmetricMatrix b(5);

  a(0, 0) = 1.0;
  a(1, 0) = 2.0;
  b(2, 0) = 3.0;
  b(3, 3) = 4.0;

  a.merge(b);

  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, a(0, 0), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(2.0, a(0, 1), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(3.0, a(0, 2), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, a(0, 3), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, a(0, 4), threshold);

  CPPUNIT_ASSERT_DOUBLES_EQUAL(2.0, a(1, 0), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, a(1, 1), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, a(1, 2), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, a(1, 3), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, a(1, 4), threshold);

  CPPUNIT_ASSERT_DOUBLES_EQUAL(3.0, a(2, 0), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, a(2, 1), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, a(2, 2), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, a(2, 3), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, a(2, 4), threshold);

  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, a(3, 0), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, a(3, 1), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, a(3, 2), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(4.0, a(3, 3), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, a(3, 4), threshold);

  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, a(4, 0), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, a(4, 1), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, a(4, 2), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, a(4, 3), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, a(4, 4), threshold);
}



