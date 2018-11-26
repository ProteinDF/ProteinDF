#include <iostream>

#include "config.h"

#include "gtest/gtest.h"
#include "matrix_common.h"
#include "vector_common.h"
#include "tl_sparse_general_matrix_viennacl.h"
#include "tl_sparse_symmetric_matrix_viennacl.h"
#include "tl_dense_general_matrix_viennacl.h"
#include "tl_dense_vector_viennacl.h"

#ifdef HAVE_EIGEN
#include "tl_sparse_general_matrix_eigen.h"
#endif // HAVE_EIGEN

static const double EPS = 1.0E-10;  // std::numeric_limits<double>::epsilon();
static const double EPS2 = 1.0E-2;
static const std::string mat_save_path = "temp.sparse.gen.viennacl.save.mat";
static const std::string mat_load_path = "temp.sparse.gen.viennacl.load.mat";

TEST(TlSparseGeneralMatrix_ViennaCL, constructor) {
  TlSparseGeneralMatrix_ViennaCL a(3, 3);

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

TEST(TlSparseGeneralMatrix_ViennaCL, copyConstructor) {
  TlSparseGeneralMatrix_ViennaCL b(5, 5);
  b.set(0, 0, 1.0);
  b.set(1, 0, 2.0);
  b.set(2, 0, 3.0);
  b.set(3, 3, 4.0);

  TlSparseGeneralMatrix_ViennaCL a;
  a = b;

  EXPECT_EQ(5, a.getNumOfRows());
  EXPECT_EQ(5, a.getNumOfCols());
  EXPECT_DOUBLE_EQ(1.0, a.get(0, 0));
  EXPECT_DOUBLE_EQ(0.0, a.get(0, 1));
  EXPECT_DOUBLE_EQ(0.0, a.get(0, 2));
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

TEST(TlSparseGeneralMatrix_ViennaCL, constructorBySymmetricMatrix) {
  const int dim = 10;
  TlSparseSymmetricMatrix_ViennaCL SM1(dim);
  SM1.set(0, 0, 1.0);
  SM1.set(1, 0, 2.0);
  SM1.set(2, 0, 3.0);
  SM1.set(3, 3, 4.0);

  TlSparseGeneralMatrix_ViennaCL SM2(SM1);
  
  EXPECT_EQ(dim, SM2.getNumOfRows());
  EXPECT_EQ(dim, SM2.getNumOfCols());
  EXPECT_DOUBLE_EQ(1.0, SM2.get(0, 0));
  EXPECT_DOUBLE_EQ(2.0, SM2.get(0, 1));
  EXPECT_DOUBLE_EQ(3.0, SM2.get(0, 2));
  EXPECT_DOUBLE_EQ(0.0, SM2.get(0, 3));
  EXPECT_DOUBLE_EQ(0.0, SM2.get(0, 4));

  EXPECT_DOUBLE_EQ(2.0, SM2.get(1, 0));
  EXPECT_DOUBLE_EQ(0.0, SM2.get(1, 1));
  EXPECT_DOUBLE_EQ(0.0, SM2.get(1, 2));
  EXPECT_DOUBLE_EQ(0.0, SM2.get(1, 3));
  EXPECT_DOUBLE_EQ(0.0, SM2.get(1, 4));

  EXPECT_DOUBLE_EQ(3.0, SM2.get(2, 0));
  EXPECT_DOUBLE_EQ(0.0, SM2.get(2, 1));
  EXPECT_DOUBLE_EQ(0.0, SM2.get(2, 2));
  EXPECT_DOUBLE_EQ(0.0, SM2.get(2, 3));
  EXPECT_DOUBLE_EQ(0.0, SM2.get(2, 4));

  EXPECT_DOUBLE_EQ(0.0, SM2.get(3, 0));
  EXPECT_DOUBLE_EQ(0.0, SM2.get(3, 1));
  EXPECT_DOUBLE_EQ(0.0, SM2.get(3, 2));
  EXPECT_DOUBLE_EQ(4.0, SM2.get(3, 3));
  EXPECT_DOUBLE_EQ(0.0, SM2.get(3, 4));

  EXPECT_DOUBLE_EQ(0.0, SM2.get(4, 0));
  EXPECT_DOUBLE_EQ(0.0, SM2.get(4, 1));
  EXPECT_DOUBLE_EQ(0.0, SM2.get(4, 2));
  EXPECT_DOUBLE_EQ(0.0, SM2.get(4, 3));
  EXPECT_DOUBLE_EQ(0.0, SM2.get(4, 4));
}

TEST(TlSparseGeneralMatrix_ViennaCL, constructByDenseMat) {
  TlDenseGeneralMatrix_ViennaCL DM(10, 20);
  DM.set(1, 0, 1.0);
  DM.set(3, 3, 3.0);
  DM.set(2, 0, -0.5);
  DM.set(4, 5, 1.0E-20);
  DM.set(6, 2, -1.0E-17);
  DM.set(3, 1, 1.0E-15);
  DM.set(7, 1, -1.0E-15);

  TlSparseGeneralMatrix_ViennaCL SM(DM);

  EXPECT_EQ(DM.getNumOfRows(), SM.getNumOfRows());
  EXPECT_EQ(DM.getNumOfCols(), SM.getNumOfCols());
  EXPECT_DOUBLE_EQ(1.0, SM.get(1, 0));
  EXPECT_DOUBLE_EQ(3.0, SM.get(3, 3));
  EXPECT_DOUBLE_EQ(0.0, SM.get(4, 5));
  EXPECT_DOUBLE_EQ(0.0, SM.get(6, 2));
  EXPECT_DOUBLE_EQ(1.0E-15, SM.get(3, 1));
  EXPECT_DOUBLE_EQ(-1.0E-15, SM.get(7, 1));
}

TEST(TlSparseGeneralMatrix_ViennaCL, constructorBySparseGeneralEigen) {
  const int row = 100;
  const int col = 120;
  TlSparseGeneralMatrix_Eigen SM1(row, col);
  SM1.set(0, 0, 1.0);
  SM1.set(1, 0, 2.0);
  SM1.set(2, 0, 3.0);
  SM1.set(3, 3, 4.0);

  // TlSparseGeneralMatrix_ViennaCL SM2(SM1);
  
  // EXPECT_EQ(dim, SM2.getNumOfRows());
  // EXPECT_EQ(dim, SM2.getNumOfCols());
  // EXPECT_DOUBLE_EQ(1.0, SM2.get(0, 0));
  // EXPECT_DOUBLE_EQ(0.0, SM2.get(0, 1));
  // EXPECT_DOUBLE_EQ(0.0, SM2.get(0, 2));
  // EXPECT_DOUBLE_EQ(0.0, SM2.get(0, 3));
  // EXPECT_DOUBLE_EQ(0.0, SM2.get(0, 4));

  // EXPECT_DOUBLE_EQ(2.0, SM2.get(1, 0));
  // EXPECT_DOUBLE_EQ(0.0, SM2.get(1, 1));
  // EXPECT_DOUBLE_EQ(0.0, SM2.get(1, 2));
  // EXPECT_DOUBLE_EQ(0.0, SM2.get(1, 3));
  // EXPECT_DOUBLE_EQ(0.0, SM2.get(1, 4));

  // EXPECT_DOUBLE_EQ(3.0, SM2.get(2, 0));
  // EXPECT_DOUBLE_EQ(0.0, SM2.get(2, 1));
  // EXPECT_DOUBLE_EQ(0.0, SM2.get(2, 2));
  // EXPECT_DOUBLE_EQ(0.0, SM2.get(2, 3));
  // EXPECT_DOUBLE_EQ(0.0, SM2.get(2, 4));

  // EXPECT_DOUBLE_EQ(0.0, SM2.get(3, 0));
  // EXPECT_DOUBLE_EQ(0.0, SM2.get(3, 1));
  // EXPECT_DOUBLE_EQ(0.0, SM2.get(3, 2));
  // EXPECT_DOUBLE_EQ(4.0, SM2.get(3, 3));
  // EXPECT_DOUBLE_EQ(0.0, SM2.get(3, 4));

  // EXPECT_DOUBLE_EQ(0.0, SM2.get(4, 0));
  // EXPECT_DOUBLE_EQ(0.0, SM2.get(4, 1));
  // EXPECT_DOUBLE_EQ(0.0, SM2.get(4, 2));
  // EXPECT_DOUBLE_EQ(0.0, SM2.get(4, 3));
  // EXPECT_DOUBLE_EQ(0.0, SM2.get(4, 4));
}

TEST(TlSparseGeneralMatrix_ViennaCL, operator_eq) {
  TlSparseGeneralMatrix_ViennaCL b(5, 8);
  b.set(0, 0, 1.0);
  b.set(1, 0, 2.0);
  b.set(2, 0, 3.0);
  b.set(3, 3, 4.0);

  TlSparseGeneralMatrix_ViennaCL a(10, 20);
  a.set(0, 1, 1.0);
  a.set(5, 3, 1.0);
  a.set(9, 5, 1.0);
  
  a = b;

  EXPECT_EQ(5, a.getNumOfRows());
  EXPECT_EQ(8, a.getNumOfCols());

  EXPECT_DOUBLE_EQ(1.0, a.get(0, 0));
  EXPECT_DOUBLE_EQ(0.0, a.get(0, 1));
  EXPECT_DOUBLE_EQ(0.0, a.get(0, 2));
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

TEST(TlSparseGeneralMatrix_ViennaCL, operator_iadd) {
  const int row1 = 10;
  const int col1 = 15;
  TlSparseGeneralMatrix_ViennaCL SM1(row1, col1);
  SM1.set(0, 0, 1.0);
  SM1.set(1, 0, 2.0);
  SM1.set(2, 0, 3.0);
  SM1.set(3, 3, 4.0);

  TlSparseGeneralMatrix_ViennaCL SM2(row1, col1);
  SM2.set(2, 0, 3.0);
  SM2.set(3, 3, 4.0);
  SM2.set(5, 4, 5.0);
  
  SM1 += SM2;

  EXPECT_EQ(row1, SM2.getNumOfRows());
  EXPECT_EQ(col1, SM2.getNumOfCols());
  EXPECT_DOUBLE_EQ(1.0, SM1.get(0, 0));
  EXPECT_DOUBLE_EQ(0.0, SM1.get(0, 1));
  EXPECT_DOUBLE_EQ(0.0, SM1.get(0, 2));
  EXPECT_DOUBLE_EQ(0.0, SM1.get(0, 3));
  EXPECT_DOUBLE_EQ(0.0, SM1.get(0, 4));

  EXPECT_DOUBLE_EQ(2.0, SM1.get(1, 0));
  EXPECT_DOUBLE_EQ(0.0, SM1.get(1, 1));
  EXPECT_DOUBLE_EQ(0.0, SM1.get(1, 2));
  EXPECT_DOUBLE_EQ(0.0, SM1.get(1, 3));
  EXPECT_DOUBLE_EQ(0.0, SM1.get(1, 4));

  EXPECT_DOUBLE_EQ(6.0, SM1.get(2, 0));
  EXPECT_DOUBLE_EQ(0.0, SM1.get(2, 1));
  EXPECT_DOUBLE_EQ(0.0, SM1.get(2, 2));
  EXPECT_DOUBLE_EQ(0.0, SM1.get(2, 3));
  EXPECT_DOUBLE_EQ(0.0, SM1.get(2, 4));

  EXPECT_DOUBLE_EQ(0.0, SM1.get(3, 0));
  EXPECT_DOUBLE_EQ(0.0, SM1.get(3, 1));
  EXPECT_DOUBLE_EQ(0.0, SM1.get(3, 2));
  EXPECT_DOUBLE_EQ(8.0, SM1.get(3, 3));
  EXPECT_DOUBLE_EQ(0.0, SM1.get(3, 4));

  EXPECT_DOUBLE_EQ(5.0, SM1.get(5, 4));
}

TEST(TlSparseGeneralMatrix_ViennaCL, multi_sparsemat_densevtr) {
  const int row1 = 4;
  const int col1 = 5;
  const int vlen = col1;

  TlSparseGeneralMatrix_ViennaCL SM(row1, col1);
  SM.set(0, 0, 1.0);
  SM.set(1, 0, 2.0);
  SM.set(2, 0, 3.0);
  SM.set(3, 3, 4.0);
  //std::cout << SM << std::endl;

  TlDenseGeneralMatrix_ViennaCL DM(row1, col1);
  DM.set(0, 0, 1.0);
  DM.set(1, 0, 2.0);
  DM.set(2, 0, 3.0);
  DM.set(3, 3, 4.0);

  TlDenseVector_ViennaCL x(vlen);
  for (int i = 0; i < vlen; ++i) {
    x.set(i, i);
  }
  //std::cout << x << std::endl;

  TlDenseVector_ViennaCL y2 = DM * x;
  //std::cout << y2 << std::endl;
  TlDenseVector_ViennaCL y1 = SM * x;
  //std::cout << y1 << std::endl;

  EXPECT_EQ(y1.getSize(), row1);
  for (int i = 0; i < row1; ++i) {
    EXPECT_NEAR(y2.get(i), y1.get(i), 1E-5);
  }
}
