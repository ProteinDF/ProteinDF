#include "TlPartialSymmetricMatrixTest.h"
#include <limits>
#include "TlFile.h"
#include "TlMatrix_RLHD.h"

const double TlPartialSymmetricMatrixTest::threshold =
    std::numeric_limits<double>::epsilon();

void TlPartialSymmetricMatrixTest::testConstructer() {
  TlPartialSymmetricMatrix psm(100, 20, 40, 10);

  CPPUNIT_ASSERT_EQUAL(100, psm.getNumOfRows());
  CPPUNIT_ASSERT_EQUAL(100, psm.getNumOfCols());
  CPPUNIT_ASSERT_EQUAL(20, psm.getStartRow());
  CPPUNIT_ASSERT_EQUAL(40, psm.getStartCol());
  CPPUNIT_ASSERT_EQUAL(10, psm.getRowRange());
  CPPUNIT_ASSERT_EQUAL(10, psm.getColRange());

  TlPartialSymmetricMatrix psm2(50, 30, 40, 20);

  CPPUNIT_ASSERT_EQUAL(50, psm2.getNumOfRows());
  CPPUNIT_ASSERT_EQUAL(50, psm2.getNumOfCols());
  CPPUNIT_ASSERT_EQUAL(30, psm2.getStartRow());
  CPPUNIT_ASSERT_EQUAL(40, psm2.getStartCol());
  CPPUNIT_ASSERT_EQUAL(20, psm2.getRowRange());
  CPPUNIT_ASSERT_EQUAL(10, psm2.getColRange());
}

void TlPartialSymmetricMatrixTest::testSetGet() {
  // partial matrix [30, 30] - (46, 46) of (100, 100)
  TlPartialSymmetricMatrix psm(100, 30, 30, 16);

  psm.set(32, 32, 12.0);
  psm.set(32, 37, 34.0);
  psm.set(42, 39, -3.0);

  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, psm.get(0, 0), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, psm.get(0, 99), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, psm.get(99, 0), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, psm.get(99, 99), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(12.0, psm.get(32, 32), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(34.0, psm.get(32, 37), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(34.0, psm.get(37, 32), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(-3.0, psm.get(42, 39), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(-3.0, psm.get(39, 42), threshold);
}

void TlPartialSymmetricMatrixTest::testAdd() {
  // partial matrix [30, 30] - (46, 46) of (100, 100)
  TlPartialSymmetricMatrix psm(100, 30, 30, 16);

  psm.set(32, 32, 12.0);
  psm.set(32, 37, 34.0);
  psm.set(42, 39, -3.0);

  psm.add(36, 36, -3.0);
  psm.add(32, 37, -4.0);
  psm.add(42, 39, 5.0);
  psm.add(34, 41, -4.0);

  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, psm.get(0, 0), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, psm.get(0, 99), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, psm.get(99, 0), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, psm.get(99, 99), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(12.0, psm.get(32, 32), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(-3.0, psm.get(36, 36), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(30.0, psm.get(32, 37), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(30.0, psm.get(37, 32), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(2.0, psm.get(42, 39), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(2.0, psm.get(39, 42), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(-4.0, psm.get(34, 41), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(-4.0, psm.get(41, 34), threshold);
}

void TlPartialSymmetricMatrixTest::testSetGet2() {
  // partial matrix [0, 16] - (16, 32) of (100, 100)
  TlPartialSymmetricMatrix psm(100, 0, 16, 16);

  psm.set(0, 20, 12.0);
  psm.set(10, 27, 34.0);

  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, psm.get(0, 0), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, psm.get(0, 99), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, psm.get(99, 0), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, psm.get(99, 99), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(12.0, psm.get(0, 20), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(34.0, psm.get(10, 27), threshold);
}

void TlPartialSymmetricMatrixTest::testAdd2() {
  // partial matrix [0, 16] - (16, 32) of (100, 100)
  TlPartialSymmetricMatrix psm(100, 0, 16, 16);

  psm.set(0, 20, 12.0);
  psm.set(10, 27, 34.0);
  psm.add(5, 19, -2.0);
  psm.add(10, 27, -4.0);

  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, psm.get(0, 0), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, psm.get(0, 99), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, psm.get(99, 0), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, psm.get(99, 99), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(12.0, psm.get(0, 20), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(30.0, psm.get(10, 27), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(-2.0, psm.get(5, 19), threshold);
}

void TlPartialSymmetricMatrixTest::testSetGet3() {
  TlPartialSymmetricMatrix psm(159, 144, 0, 144);

  psm.add(144, 15, -1.7);
  psm.add(16, 145, 3.5);
  psm.add(146, 22, -0.5);

  CPPUNIT_ASSERT_DOUBLES_EQUAL(-1.7, psm.get(144, 15), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(-1.7, psm.get(15, 144), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(3.5, psm.get(145, 16), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(3.5, psm.get(16, 145), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(-0.5, psm.get(146, 22), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(-0.5, psm.get(22, 146), threshold);
}

void TlPartialSymmetricMatrixTest::testGetRowVector() {
  TlPartialSymmetricMatrix psm1(100, 0, 20, 20);  // (0, 20) - (20, 40)

  psm1.set(14, 22, -1.7);
  psm1.set(14, 25, -0.5);
  psm1.set(14, 35, 3.5);

  TlVector v1 = psm1.getRowVector(14);

  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, v1.get(20), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, v1.get(21), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(-1.7, v1.get(22), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, v1.get(23), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, v1.get(24), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(-0.5, v1.get(25), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, v1.get(26), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, v1.get(27), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, v1.get(28), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, v1.get(29), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, v1.get(30), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, v1.get(31), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, v1.get(32), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, v1.get(33), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, v1.get(34), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(3.5, v1.get(35), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, v1.get(36), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, v1.get(37), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, v1.get(38), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, v1.get(39), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, v1.get(40), threshold);

  TlPartialSymmetricMatrix psm2(100, 60, 60, 20);  // (60, 60) - (80, 80)

  psm2.set(68, 61, -1.7);
  psm2.set(68, 65, 5.1);
  psm2.set(68, 68, -0.5);
  psm2.set(68, 72, 3.5);
  psm2.set(62, 68, -4.7);
  psm2.set(69, 68, 1.9);
  psm2.set(78, 68, 2.5);

  TlVector v2 = psm2.getRowVector(68);

  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, v2.get(60), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(-1.7, v2.get(61), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(-4.7, v2.get(62), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, v2.get(63), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, v2.get(64), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(5.1, v2.get(65), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, v2.get(66), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, v2.get(67), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(-0.5, v2.get(68), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.9, v2.get(69), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, v2.get(70), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, v2.get(71), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(3.5, v2.get(72), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, v2.get(73), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, v2.get(74), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, v2.get(75), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, v2.get(76), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, v2.get(77), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(2.5, v2.get(78), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, v2.get(79), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, v2.get(80), threshold);
}
