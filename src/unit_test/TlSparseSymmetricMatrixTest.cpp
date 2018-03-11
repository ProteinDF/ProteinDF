#include <limits>
#include "TlSparseSymmetricMatrix.h"
#include "gtest/gtest.h"

static const double EPS = 1.0E-10;  // std::numeric_limits<double>::epsilon();

TEST(TlSparseSymmetricMatrixTest, constructer) {
  TlSparseSymmetricMatrix a(3);

  EXPECT_NEAR(0.0, a(0, 0), EPS);
  EXPECT_NEAR(0.0, a(0, 1), EPS);
  EXPECT_NEAR(0.0, a(0, 2), EPS);
  EXPECT_NEAR(0.0, a(1, 0), EPS);
  EXPECT_NEAR(0.0, a(1, 1), EPS);
  EXPECT_NEAR(0.0, a(1, 2), EPS);
  EXPECT_NEAR(0.0, a(2, 0), EPS);
  EXPECT_NEAR(0.0, a(2, 1), EPS);
  EXPECT_NEAR(0.0, a(2, 2), EPS);
}

TEST(TlSparseSymmetricMatrix, merge) {
  TlSparseSymmetricMatrix a(5);
  TlSparseSymmetricMatrix b(5);

  a(0, 0) = 1.0;
  a(1, 0) = 2.0;
  b(2, 0) = 3.0;
  b(3, 3) = 4.0;

  a.merge(b);

  EXPECT_NEAR(1.0, a(0, 0), EPS);
  EXPECT_NEAR(2.0, a(0, 1), EPS);
  EXPECT_NEAR(3.0, a(0, 2), EPS);
  EXPECT_NEAR(0.0, a(0, 3), EPS);
  EXPECT_NEAR(0.0, a(0, 4), EPS);

  EXPECT_NEAR(2.0, a(1, 0), EPS);
  EXPECT_NEAR(0.0, a(1, 1), EPS);
  EXPECT_NEAR(0.0, a(1, 2), EPS);
  EXPECT_NEAR(0.0, a(1, 3), EPS);
  EXPECT_NEAR(0.0, a(1, 4), EPS);

  EXPECT_NEAR(3.0, a(2, 0), EPS);
  EXPECT_NEAR(0.0, a(2, 1), EPS);
  EXPECT_NEAR(0.0, a(2, 2), EPS);
  EXPECT_NEAR(0.0, a(2, 3), EPS);
  EXPECT_NEAR(0.0, a(2, 4), EPS);

  EXPECT_NEAR(0.0, a(3, 0), EPS);
  EXPECT_NEAR(0.0, a(3, 1), EPS);
  EXPECT_NEAR(0.0, a(3, 2), EPS);
  EXPECT_NEAR(4.0, a(3, 3), EPS);
  EXPECT_NEAR(0.0, a(3, 4), EPS);

  EXPECT_NEAR(0.0, a(4, 0), EPS);
  EXPECT_NEAR(0.0, a(4, 1), EPS);
  EXPECT_NEAR(0.0, a(4, 2), EPS);
  EXPECT_NEAR(0.0, a(4, 3), EPS);
  EXPECT_NEAR(0.0, a(4, 4), EPS);
}
