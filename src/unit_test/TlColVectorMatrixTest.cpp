#include <iostream>
#include <limits>
#include "TlColVectorMatrix.h"
#include "TlMatrix.h"
#include "gtest/gtest.h"

static const double EPS = 1.0E-10;  // std::numeric_limits<double>::epsilon();

TEST(TlColVectorMatrix, constructer) {
  TlColVectorMatrix A(100, 600);

  EXPECT_EQ(100, A.getSizeOfVector());
  EXPECT_EQ(600, A.getNumOfVectors());
  EXPECT_NEAR(0.0, A.get(0, 0), EPS);
}

TEST(TlColVectorMatrix, constructer2) {
  TlColVectorMatrix A0(100, 600, 10, 0);
  TlColVectorMatrix A1(100, 600, 10, 1);

  EXPECT_EQ(100, A0.getSizeOfVector());
  EXPECT_EQ(600, A0.getNumOfVectors());
  EXPECT_EQ(10, A0.getNumOfSubunits());
  EXPECT_EQ(0, A0.getSubunitID());

  EXPECT_EQ(100, A1.getSizeOfVector());
  EXPECT_EQ(600, A1.getNumOfVectors());
  EXPECT_EQ(10, A1.getNumOfSubunits());
  EXPECT_EQ(1, A1.getSubunitID());
}

TEST(TlColVectorMatrix, resize) {
  TlColVectorMatrix A(100, 100);
  A.resize(200, 100);
  EXPECT_EQ(200, A.getSizeOfVector());
  EXPECT_EQ(100, A.getNumOfVectors());

  TlColVectorMatrix B(100, 100);
  B.resize(100, 200);
  EXPECT_EQ(100, B.getSizeOfVector());
  EXPECT_EQ(200, B.getNumOfVectors());

  TlColVectorMatrix C(50, 100);
  C.resize(50, 100);
  EXPECT_EQ(50, C.getSizeOfVector());
  EXPECT_EQ(100, C.getNumOfVectors());

  TlColVectorMatrix D(100, 50);
  D.resize(100, 50);
  EXPECT_EQ(100, D.getSizeOfVector());
  EXPECT_EQ(50, D.getNumOfVectors());
}

TEST(TlColVectorMatrix, contents) {
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
      EXPECT_NEAR(matA.get(r, c), vecA.get(r, c), EPS);
    }
  }
}

TEST(TlColVectorMatrix, save_load) {
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
      EXPECT_NEAR(matA.get(r, c), vecB.get(r, c), EPS);
    }
  }
}

TEST(TlColVectorMatrix, toTlMatrix) {
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
      EXPECT_NEAR(matA.get(r, c), matB.get(r, c), EPS);
    }
  }
}
