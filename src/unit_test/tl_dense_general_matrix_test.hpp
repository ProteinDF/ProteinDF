#include <iostream>

#include "gtest/gtest.h"
#include "matrix_common.h"

static const double EPS = 1.0E-10;  // std::numeric_limits<double>::epsilon();
static const double EPS2 = 1.0E-2;
static const std::string mat_save_load_path = "temp.gen.blas.save_load.mat";

// -----------------------------------------------------------------------------
// test
// -----------------------------------------------------------------------------
template<typename T>
TEST(T, constructer) {
  T a(3, 3);

  EXPECT_EQ(3, a.getNumOfRows());
  EXPECT_EQ(3, a.getNumOfCols());
  EXPECT_DOUBLE_EQ(0.0, a.get(0, 0));
  EXPECT_DOUBLE_EQ(0.0, a.get(0, 1));
  EXPECT_DOUBLE_EQ(0.0, a.get(0, 2));
  EXPECT_DOUBLE_EQ(0.0, a.get(1, 0));
  EXPECT_DOUBLE_EQ(0.0, a.get(1, 1));
  EXPECT_DOUBLE_EQ(0.0, a.get(1, 2));
  EXPECT_DOUBLE_EQ(0.0, a.get(2, 0));
  EXPECT_DOUBLE_EQ(0.0, a.get(2, 1));
  EXPECT_DOUBLE_EQ(0.0, a.get(2, 2));
}
