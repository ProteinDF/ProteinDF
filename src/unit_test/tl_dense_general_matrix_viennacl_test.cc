#include <iostream>

#include "config.h"
#include "gtest/gtest.h"
#include "matrix_common.h"
#include "tl_dense_general_matrix_viennacl.h"
#include "tl_dense_vector_viennacl.h"
#include "vector_common.h"
#include "tl_sparse_general_matrix_viennacl.h"

#ifdef HAVE_EIGEN
#include "tl_dense_general_matrix_eigen.h"
#endif // HAVE_EIGEN

static const double EPS = 1.0E-10;  // std::numeric_limits<double>::epsilon();
static const double EPS2 = 1.0E-2;
static const std::string mat_save_path = "temp.gen.viennacl.save.mat";
static const std::string mat_load_path = "temp.gen.viennacl.load.mat";

// -----------------------------------------------------------------------------
// test
// -----------------------------------------------------------------------------
// TEST(TlDenseGeneralMatrix_ViennaCL, constructBySparseGeneralMatrix) {
//   const int row = 20;
//   const int col = 30;
//   TlSparseGeneralMatrix_ViennaCL SM(row, col);

//   SM.set(1, 0, 1.0);
//   SM.set(2, 3, 4.0);
//   SM.set(5, 8, -10.0);

//   TlDenseGeneralMatrix_ViennaCL DM = SM;
//   EXPECT_EQ(row, DM.getNumOfRows());
//   EXPECT_EQ(col, DM.getNumOfCols());
//   EXPECT_DOUBLE_EQ(1.0, DM.get(1, 0));
//   EXPECT_DOUBLE_EQ(4.0, DM.get(2, 3));
//   EXPECT_DOUBLE_EQ(-10.0, DM.get(5, 8));
// }

TEST(TlDenseGeneralMatrix_ViennaCL, copyFromEigen) {
  TlDenseGeneralMatrix_Eigen a = getMatrixA<TlDenseGeneralMatrix_Eigen>();
  TlDenseGeneralMatrix_ViennaCL b = a;

  EXPECT_EQ(a.getNumOfRows(), b.getNumOfRows());
  EXPECT_EQ(a.getNumOfCols(), b.getNumOfCols());
  for (int i = 0; i < a.getNumOfRows(); ++i) {
    for (int j = 0; j < a.getNumOfCols(); ++j) {
      EXPECT_DOUBLE_EQ(a.get(i, j), b.get(i, j));
    }
  }
}

#ifdef HAVE_EIGEN
TEST(TlDenseGeneralMatrix_ViennaCL, constructBy_DenseGeneralMatrix_Eigen) {
  const int row = 20;
  const int col = 30;
  TlDenseGeneralMatrix_Eigen DM_E(row, col);

  DM_E.set(1, 0, 1.0);
  DM_E.set(2, 3, 4.0);
  DM_E.set(5, 8, -10.0);

  TlDenseGeneralMatrix_ViennaCL DM = DM_E;
  EXPECT_EQ(row, DM.getNumOfRows());
  EXPECT_EQ(col, DM.getNumOfCols());
  EXPECT_DOUBLE_EQ(1.0, DM.get(1, 0));
  EXPECT_DOUBLE_EQ(4.0, DM.get(2, 3));
  EXPECT_DOUBLE_EQ(-10.0, DM.get(5, 8));
}
#endif // HAVE_EIGEN

TEST(TlDenseGeneralMatrix_ViennaCL, vtr2mat) {
  const int row = 3;
  const int col = 4;
  const int elements = row * col;
  std::vector<double> vtr(elements);
  for (int i = 0; i < elements; ++i) {
    vtr[i] = i;
  }

  TlDenseGeneralMatrix_ViennaCL a(row, col);
  a.vtr2mat(vtr);
  // std::cout << a << std::endl;

  EXPECT_EQ(row, a.getNumOfRows());
  EXPECT_EQ(col, a.getNumOfCols());
  int i = 0;
  for (int c = 0; c < col; ++c) { // col-major
    for (int r = 0; r < row; ++r) {
      EXPECT_DOUBLE_EQ(vtr[i], a.get(r, c));
      ++i;
    }
  }
}

// TEST(TlDenseGeneralMatrix_ViennaCL, constructByTlSerializedData) {
//   // TODO
// }
//
// TEST(TlDenseGeneralMatrix_ViennaCL, vtr2mat) {
//   // TODO
// }


#ifdef HAVE_HDF5
// TEST(TlDenseGeneralMatrix_ViennaCL, hdf5) {
//   TlDenseGeneralMatrix_ViennaCL m = getMatrix_A();
//   m.saveHdf5("temp.mat.h5", "matrix_A");
//
//   TlDenseGeneralMatrix_ViennaCL a;
//   a.loadHdf5("temp.mat.h5", "matrix_A");
//   EXPECT_EQ(TlMatrixObject::CSFD, a.getType());
//   EXPECT_EQ(3, a.getNumOfRows());
//   EXPECT_EQ(3, a.getNumOfCols());
//   EXPECT_DOUBLE_EQ(0.0, a.get(0, 0));
//   EXPECT_DOUBLE_EQ(1.0, a.get(0, 1));
//   EXPECT_DOUBLE_EQ(2.0, a.get(0, 2));
//   EXPECT_DOUBLE_EQ(3.0, a.get(1, 0));
//   EXPECT_DOUBLE_EQ(4.0, a.get(1, 1));
//   EXPECT_DOUBLE_EQ(5.0, a.get(1, 2));
//   EXPECT_DOUBLE_EQ(6.0, a.get(2, 0));
//   EXPECT_DOUBLE_EQ(7.0, a.get(2, 1));
//   EXPECT_DOUBLE_EQ(8.0, a.get(2, 2));
// }
#endif  // HAVE_HDF5

// TEST(TlDenseGeneralMatrix_ViennaCL, inverse)
//
// {
//   TlDenseGeneralMatrix_ViennaCL a = getMatrix_C();
//   TlDenseGeneralMatrix_ViennaCL b = a;
//
//   b.inverse();
//
//   TlDenseGeneralMatrix_ViennaCL c = a * b;
//   EXPECT_NEAR(1.0, c.get(0, 0), EPS);
//   EXPECT_NEAR(0.0, c.get(0, 1), EPS);
//   EXPECT_NEAR(0.0, c.get(0, 2), EPS);
//   EXPECT_NEAR(0.0, c.get(1, 0), EPS);
//   EXPECT_NEAR(1.0, c.get(1, 1), EPS);
//   EXPECT_NEAR(0.0, c.get(1, 2), EPS);
//   EXPECT_NEAR(0.0, c.get(2, 0), EPS);
//   EXPECT_NEAR(0.0, c.get(2, 1), EPS);
//   EXPECT_NEAR(1.0, c.get(2, 2), EPS);
// }
//
// TEST(TlDenseGeneralMatrix_ViennaCL, operator_mul_AB) {
//   TlDenseGeneralMatrix_ViennaCL a =
//   getMatrixA<TlDenseGeneralMatrix_ViennaCL>();
//   TlDenseGeneralMatrix_ViennaCL b =
//   getMatrixB<TlDenseGeneralMatrix_ViennaCL>();
//   TlDenseGeneralMatrix_ViennaCL c = a * b;
//
//   EXPECT_EQ(3, c.getNumOfRows());
//   EXPECT_EQ(3, c.getNumOfCols());
//   EXPECT_DOUBLE_EQ(5.0, c.get(0, 0));
//   EXPECT_DOUBLE_EQ(14.0, c.get(0, 1));
//   EXPECT_DOUBLE_EQ(23.0, c.get(0, 2));
//   EXPECT_DOUBLE_EQ(14.0, c.get(1, 0));
//   EXPECT_DOUBLE_EQ(50.0, c.get(1, 1));
//   EXPECT_DOUBLE_EQ(86.0, c.get(1, 2));
//   EXPECT_DOUBLE_EQ(23.0, c.get(2, 0));
//   EXPECT_DOUBLE_EQ(86.0, c.get(2, 1));
//   EXPECT_DOUBLE_EQ(149.0, c.get(2, 2));
// }

// [0 1 2]   [0]   [ 5]
// [3 4 5] x [1] = [14]
// [6 7 8]   [2]   [23]
TEST(TlDenseGeneralMatrix_ViennaCL, operator_mul_mat_vec) {
  TlDenseGeneralMatrix_ViennaCL a = getMatrixA<TlDenseGeneralMatrix_ViennaCL>();
  TlDenseVector_ViennaCL v = getVectorA<TlDenseVector_ViennaCL>();

  TlDenseVector_ViennaCL z = a * v;

  EXPECT_EQ(3, z.getSize());
  EXPECT_DOUBLE_EQ(5.0, z.get(0));
  EXPECT_DOUBLE_EQ(14.0, z.get(1));
  EXPECT_DOUBLE_EQ(23.0, z.get(2));
}

//           [ 1  2  3  4]
// [0 1 2] * [ 5  6  7  8] = [23 26 29 32]
//           [ 9 10 11 12]
TEST(TlDenseGeneralMatrix_ViennaCL, operator_mul_vec_mat) {
  TlDenseGeneralMatrix_ViennaCL a = getMatrixD<TlDenseGeneralMatrix_ViennaCL>();
  TlDenseVector_ViennaCL v = getVectorA<TlDenseVector_ViennaCL>();

  TlDenseVector_ViennaCL z = v * a;

  EXPECT_EQ(4, z.getSize());
  EXPECT_DOUBLE_EQ(23.0, z.get(0));
  EXPECT_DOUBLE_EQ(26.0, z.get(1));
  EXPECT_DOUBLE_EQ(29.0, z.get(2));
  EXPECT_DOUBLE_EQ(32.0, z.get(3));
}

// TEST(TlDenseGeneralMatrix_ViennaCL, dot) {
//   TlDenseGeneralMatrix_ViennaCL a = getMatrix_A();
//   TlDenseGeneralMatrix_ViennaCL b = getMatrix_B();
//   a.dot(b);
//
//   EXPECT_DOUBLE_EQ(0.0, a.get(0, 0));
//   EXPECT_DOUBLE_EQ(3.0, a.get(1, 0));
//   EXPECT_DOUBLE_EQ(12.0, a.get(2, 0));
//   EXPECT_DOUBLE_EQ(3.0, a.get(0, 1));
//   EXPECT_DOUBLE_EQ(16.0, a.get(1, 1));
//   EXPECT_DOUBLE_EQ(35.0, a.get(2, 1));
//   EXPECT_DOUBLE_EQ(12.0, a.get(0, 2));
//   EXPECT_DOUBLE_EQ(35.0, a.get(1, 2));
//   EXPECT_DOUBLE_EQ(64.0, a.get(2, 2));
// }
//
// TEST(TlDenseGeneralMatrix_ViennaCL, sum) {
//   TlDenseGeneralMatrix_ViennaCL a = getMatrix_A();
//   double s = a.sum();
//
//   EXPECT_DOUBLE_EQ(36.0, s);
// }
//
// TEST(TlDenseGeneralMatrix_ViennaCL, solveLinearLeastSquaresProblem) {
//   TlDenseGeneralMatrix_ViennaCL A(6, 5);
//   A(0, 0) = -0.09;
//   A(0, 1) = 0.14;
//   A(0, 2) = -0.46;
//   A(0, 3) = 0.68;
//   A(0, 4) = 1.29;
//   A(1, 0) = -1.56;
//   A(1, 1) = 0.20;
//   A(1, 2) = 0.29;
//   A(1, 3) = 1.09;
//   A(1, 4) = 0.51;
//   A(2, 0) = -1.48;
//   A(2, 1) = -0.43;
//   A(2, 2) = 0.89;
//   A(2, 3) = -0.71;
//   A(2, 4) = -0.96;
//   A(3, 0) = -1.09;
//   A(3, 1) = 0.84;
//   A(3, 2) = 0.77;
//   A(3, 3) = 2.11;
//   A(3, 4) = -1.27;
//   A(4, 0) = 0.08;
//   A(4, 1) = 0.55;
//   A(4, 2) = -1.13;
//   A(4, 3) = 0.14;
//   A(4, 4) = 1.74;
//   A(5, 0) = -1.59;
//   A(5, 1) = -0.72;
//   A(5, 2) = 1.06;
//   A(5, 3) = 1.24;
//   A(5, 4) = 0.34;
//
//   TlDenseGeneralMatrix_ViennaCL B(6, 1);
//   B(0, 0) = 7.4;
//   B(1, 0) = 4.2;
//   B(2, 0) = -8.3;
//   B(3, 0) = 1.8;
//   B(4, 0) = 8.6;
//   B(5, 0) = 2.1;
//
//   TlDenseGeneralMatrix_ViennaCL X = A.solveLinearLeastSquaresProblem(B);
//   // X.print(std::cout);
//
//   TlDenseGeneralMatrix_ViennaCL AX = A * X;
//   // AX.print(std::cout);
//
//   EXPECT_EQ(5, X.getNumOfRows());
//   EXPECT_EQ(1, X.getNumOfCols());
//   EXPECT_NEAR(B(0, 0), AX(0, 0), EPS2);
//   EXPECT_NEAR(B(1, 0), AX(1, 0), EPS2);
//   EXPECT_NEAR(B(2, 0), AX(2, 0), EPS2);
//   EXPECT_NEAR(B(3, 0), AX(3, 0), EPS2);
//   EXPECT_NEAR(B(4, 0), AX(4, 0), EPS2);
// }

TEST(TlDenseGeneralMatrix_ViennaCL, reverseColumns) {
  const int row = 20;
  const int col = 20;
  TlDenseGeneralMatrix_Eigen ref = getTlMatrix<TlDenseGeneralMatrix_Eigen>(row, col);
  TlDenseGeneralMatrix_ViennaCL a(ref);

  a.reverseColumns();
  // std::cout << ref << std::endl;
  // std::cout << a << std::endl;

  for (int i = 0; i < row; ++i) {
    for (int j = 0; j < col; ++j) {
      EXPECT_DOUBLE_EQ(ref.get(i, col -j -1), a.get(i, j));
    }
  }
}
