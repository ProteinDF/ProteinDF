#include "tl_sparse_vector.h"
#include <limits>
#include "gtest/gtest.h"

const double EPS = 1.0E-10;  // std::numeric_limits<double>::epsilon();

TEST(TlSparseVector, constructer) {
  TlSparseVector a(100);

  for (int i = 0; i < 100; ++i) {
    EXPECT_NEAR(0.0, a[i], EPS);
  }
}

TEST(TlSparseVector, add_double) {
  TlSparseVector a(100);

  a[10] = 10.0;
  a[10] += 20.0;
  a[51] += 3.0;

  EXPECT_NEAR(30.0, a[10], EPS);
  EXPECT_NEAR(3.0, a[51], EPS);
  EXPECT_NEAR(0.0, a[80], EPS);
}

TEST(TlSparseVector, add_TlSparseVector) {
  TlSparseVector a(100);
  TlSparseVector b(100);

  a[2] = 2.0;
  a[5] = 5.0;
  a[31] = 31.0;
  b[17] = 17.0;
  b[31] = 20.0;

  a += b;

  EXPECT_NEAR(2.0, a[2], EPS);
  EXPECT_NEAR(5.0, a[5], EPS);
  EXPECT_NEAR(17.0, a[17], EPS);
  EXPECT_NEAR(51.0, a[31], EPS);
  EXPECT_NEAR(0.0, a[99], EPS);
}
