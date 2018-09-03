#include "tl_sparse_symmetric_matrix_viennacl.h"
#include <limits>
#include "gtest/gtest.h"

static const double EPS = 1.0E-10;  // std::numeric_limits<double>::epsilon();

TEST(TlSparseSymmetricMatrix_ViennaCL, constructor) {
  TlSparseSymmetricMatrix_ViennaCL a(3);

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

TEST(TlSparseSymmetricMatrix_ViennaCL, copyConstructor) {
  TlSparseSymmetricMatrix_ViennaCL b(5);
  b.set(0, 0, 1.0);
  b.set(1, 0, 2.0);
  b.set(2, 0, 3.0);
  b.set(3, 3, 4.0);

  TlSparseSymmetricMatrix_ViennaCL a;
  a = b;

  EXPECT_EQ(5, a.getNumOfRows());
  EXPECT_EQ(5, a.getNumOfCols());
  EXPECT_DOUBLE_EQ(1.0, a.get(0, 0));
  EXPECT_DOUBLE_EQ(2.0, a.get(0, 1));
  EXPECT_DOUBLE_EQ(3.0, a.get(0, 2));
  EXPECT_DOUBLE_EQ(0.0, a.get(0, 3));
  EXPECT_DOUBLE_EQ(0.0, a.get(0, 4));

  EXPECT_DOUBLE_EQ(2.0, a.get(1, 0));
  EXPECT_DOUBLE_EQ(0.0, a.get(1, 1));
  EXPECT_DOUBLE_EQ(0.0, a.get(1, 2));
  EXPECT_DOUBLE_EQ(0.0, a.get(1, 3));
  EXPECT_DOUBLE_EQ(0.0, a.get(1, 4));

  EXPECT_DOUBLE_EQ(3.0, a.get(2, 0));
  EXPECT_DOUBLE_EQ(0.0, a.get(2, 1));
  EXPECT_DOUBLE_EQ(0.0, a.get(2, 2));
  EXPECT_DOUBLE_EQ(0.0, a.get(2, 3));
  EXPECT_DOUBLE_EQ(0.0, a.get(2, 4));

  EXPECT_DOUBLE_EQ(0.0, a.get(3, 0));
  EXPECT_DOUBLE_EQ(0.0, a.get(3, 1));
  EXPECT_DOUBLE_EQ(0.0, a.get(3, 2));
  EXPECT_DOUBLE_EQ(4.0, a.get(3, 3));
  EXPECT_DOUBLE_EQ(0.0, a.get(3, 4));

  EXPECT_DOUBLE_EQ(0.0, a.get(4, 0));
  EXPECT_DOUBLE_EQ(0.0, a.get(4, 1));
  EXPECT_DOUBLE_EQ(0.0, a.get(4, 2));
  EXPECT_DOUBLE_EQ(0.0, a.get(4, 3));
  EXPECT_DOUBLE_EQ(0.0, a.get(4, 4));
}

TEST(TlSparseSymmetricMatrix_ViennaCL, operator_eq) {
  TlSparseSymmetricMatrix_ViennaCL b(5);
  b.set(0, 0, 1.0);
  b.set(1, 0, 2.0);
  b.set(2, 0, 3.0);
  b.set(3, 3, 4.0);

  TlSparseSymmetricMatrix_ViennaCL a(10);
  a.set(0, 1, 1.0);
  a.set(5, 3, 1.0);
  a.set(9, 5, 1.0);
  
  a = b;

  EXPECT_EQ(5, a.getNumOfRows());
  EXPECT_EQ(5, a.getNumOfCols());

  EXPECT_DOUBLE_EQ(1.0, a.get(0, 0));
  EXPECT_DOUBLE_EQ(2.0, a.get(0, 1));
  EXPECT_DOUBLE_EQ(3.0, a.get(0, 2));
  EXPECT_DOUBLE_EQ(0.0, a.get(0, 3));
  EXPECT_DOUBLE_EQ(0.0, a.get(0, 4));

  EXPECT_DOUBLE_EQ(2.0, a.get(1, 0));
  EXPECT_DOUBLE_EQ(0.0, a.get(1, 1));
  EXPECT_DOUBLE_EQ(0.0, a.get(1, 2));
  EXPECT_DOUBLE_EQ(0.0, a.get(1, 3));
  EXPECT_DOUBLE_EQ(0.0, a.get(1, 4));

  EXPECT_DOUBLE_EQ(3.0, a.get(2, 0));
  EXPECT_DOUBLE_EQ(0.0, a.get(2, 1));
  EXPECT_DOUBLE_EQ(0.0, a.get(2, 2));
  EXPECT_DOUBLE_EQ(0.0, a.get(2, 3));
  EXPECT_DOUBLE_EQ(0.0, a.get(2, 4));

  EXPECT_DOUBLE_EQ(0.0, a.get(3, 0));
  EXPECT_DOUBLE_EQ(0.0, a.get(3, 1));
  EXPECT_DOUBLE_EQ(0.0, a.get(3, 2));
  EXPECT_DOUBLE_EQ(4.0, a.get(3, 3));
  EXPECT_DOUBLE_EQ(0.0, a.get(3, 4));

  EXPECT_DOUBLE_EQ(0.0, a.get(4, 0));
  EXPECT_DOUBLE_EQ(0.0, a.get(4, 1));
  EXPECT_DOUBLE_EQ(0.0, a.get(4, 2));
  EXPECT_DOUBLE_EQ(0.0, a.get(4, 3));
  EXPECT_DOUBLE_EQ(0.0, a.get(4, 4));
}
