#include <iostream>
#include <limits>
#include "TlMatrix.h"
#include "TlRowVectorMatrix.h"
#include "gtest/gtest.h"

static const double EPS = 1.0E-10;  // std::numeric_limits<double>::epsilon();

TEST(TlRowVectorMatrix, constructer) {
  TlRowVectorMatrix A(100, 600);

  EXPECT_EQ(100, A.getNumOfRows());
  EXPECT_EQ(600, A.getNumOfCols());
  EXPECT_NEAR(0.0, A.get(0, 0), EPS);
}

TEST(TlRowVectorMatrix, constructer2) {
  TlRowVectorMatrix A0(100, 600, 10, 0);
  TlRowVectorMatrix A1(100, 600, 10, 1);

  EXPECT_EQ(100, A0.getNumOfRows());
  EXPECT_EQ(600, A0.getNumOfCols());
  EXPECT_EQ(10, A0.getNumOfSubunits());
  EXPECT_EQ(0, A0.getSubunitID());

  EXPECT_EQ(100, A1.getNumOfRows());
  EXPECT_EQ(600, A1.getNumOfCols());
  EXPECT_EQ(10, A1.getNumOfSubunits());
  EXPECT_EQ(1, A1.getSubunitID());
}

TEST(TlRowVectorMatrix, resize) {
  TlRowVectorMatrix A(100, 100);
  A.resize(200, 100);
  EXPECT_EQ(200, A.getNumOfRows());
  EXPECT_EQ(100, A.getNumOfCols());

  TlRowVectorMatrix B(100, 100);
  B.resize(100, 200);
  EXPECT_EQ(100, B.getNumOfRows());
  EXPECT_EQ(200, B.getNumOfCols());

  TlRowVectorMatrix C(50, 100);
  C.resize(50, 100);
  EXPECT_EQ(50, C.getNumOfRows());
  EXPECT_EQ(100, C.getNumOfCols());

  TlRowVectorMatrix D(100, 50);
  D.resize(100, 50);
  EXPECT_EQ(100, D.getNumOfRows());
  EXPECT_EQ(50, D.getNumOfCols());
}

TEST(TlRowVectorMatrix, contents) {
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
      EXPECT_NEAR(matA.get(r, c), vecA.get(r, c), EPS);
    }
  }
}

TEST(TlRowVectorMatrix, save_load) {
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
  vecB.load("/tmp/vecA.mat", 0);

  // test
  for (int r = 0; r < maxRow; ++r) {
    for (int c = 0; c < maxCol; ++c) {
      EXPECT_NEAR(matA.get(r, c), vecB.get(r, c), EPS);
    }
  }
}

TEST(TlRowVectorMatrix, toTlMatrix) {
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
      EXPECT_NEAR(matA.get(r, c), matB.get(r, c), EPS);
    }
  }
}
