#include <iostream>

#include "config.h"
#include "gtest/gtest.h"
#include "matrix_common.h"
#include "tl_dense_general_matrix_blas_old.h"
#include "tl_dense_general_matrix_eigen_old.h"
#include "tl_dense_symmetric_matrix_eigen_old.h"

static const double EPS = 1.0E-10;  // std::numeric_limits<double>::epsilon();
static const double EPS2 = 1.0E-2;
static const std::string mat_path = "temp.mat";

// -----------------------------------------------------------------------------
// test
// -----------------------------------------------------------------------------
TEST(TlDenseGeneralMatrix_Eigen_Old, constructer) {
  TlDenseGeneralMatrix_Eigen_Old a(3, 3);

  EXPECT_EQ(TlMatrixObject::CSFD, a.getType());
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

// TEST(TlDenseGeneralMatrix_Eigen_Old, constructByTlSerializedData) {
//   // TODO
// }
//
// TEST(TlDenseGeneralMatrix_Eigen_Old, vtr2mat) {
//   // TODO
// }

TEST(TlDenseGeneralMatrix_Eigen_Old, copyConstructer) {
  TlDenseGeneralMatrix_Eigen_Old a = getMatrixA<TlDenseGeneralMatrix_Eigen_Old>();
  TlDenseGeneralMatrix_Eigen_Old c(a);

  EXPECT_EQ(3, c.getNumOfRows());
  EXPECT_EQ(3, c.getNumOfCols());
  EXPECT_DOUBLE_EQ(0.0, c.get(0, 0));
  EXPECT_DOUBLE_EQ(1.0, c.get(0, 1));
  EXPECT_DOUBLE_EQ(2.0, c.get(0, 2));
  EXPECT_DOUBLE_EQ(3.0, c.get(1, 0));
  EXPECT_DOUBLE_EQ(4.0, c.get(1, 1));
  EXPECT_DOUBLE_EQ(5.0, c.get(1, 2));
  EXPECT_DOUBLE_EQ(6.0, c.get(2, 0));
  EXPECT_DOUBLE_EQ(7.0, c.get(2, 1));
  EXPECT_DOUBLE_EQ(8.0, c.get(2, 2));
}

TEST(TlDenseGeneralMatrix_Eigen_Old, symmat2mat) {
  TlDenseSymmetricMatrix_Eigen_Old s(4);
  s.set(0, 0, 1.0);
  s.set(1, 2, 2.0);
  s.set(3, 2, 3.0);

  TlDenseGeneralMatrix_Eigen_Old m(s);

  EXPECT_EQ(TlMatrixObject::CSFD, m.getType());
  EXPECT_EQ(4, m.getNumOfRows());
  EXPECT_EQ(4, m.getNumOfCols());
  EXPECT_DOUBLE_EQ(1.0, m.get(0, 0));
  EXPECT_DOUBLE_EQ(0.0, m.get(0, 1));
  EXPECT_DOUBLE_EQ(0.0, m.get(0, 2));
  EXPECT_DOUBLE_EQ(0.0, m.get(0, 3));
  EXPECT_DOUBLE_EQ(0.0, m.get(1, 0));
  EXPECT_DOUBLE_EQ(0.0, m.get(1, 1));
  EXPECT_DOUBLE_EQ(2.0, m.get(1, 2));
  EXPECT_DOUBLE_EQ(0.0, m.get(1, 3));
  EXPECT_DOUBLE_EQ(0.0, m.get(2, 0));
  EXPECT_DOUBLE_EQ(2.0, m.get(2, 1));
  EXPECT_DOUBLE_EQ(0.0, m.get(2, 2));
  EXPECT_DOUBLE_EQ(3.0, m.get(2, 3));
  EXPECT_DOUBLE_EQ(0.0, m.get(3, 0));
  EXPECT_DOUBLE_EQ(0.0, m.get(3, 1));
  EXPECT_DOUBLE_EQ(3.0, m.get(3, 2));
  EXPECT_DOUBLE_EQ(0.0, m.get(3, 3));
}

TEST(TlDenseGeneralMatrix_Eigen_Old, resize) {
  TlDenseGeneralMatrix_Eigen_Old a(2, 3);
  a.set(0, 1, 2.0);
  a.set(1, 2, -4.0);

  EXPECT_EQ(2, a.getNumOfRows());
  EXPECT_EQ(3, a.getNumOfCols());
  EXPECT_DOUBLE_EQ(0.0, a.get(0, 0));
  EXPECT_DOUBLE_EQ(2.0, a.get(0, 1));
  EXPECT_DOUBLE_EQ(0.0, a.get(0, 2));
  EXPECT_DOUBLE_EQ(0.0, a.get(1, 0));
  EXPECT_DOUBLE_EQ(0.0, a.get(1, 1));
  EXPECT_DOUBLE_EQ(-4.0, a.get(1, 2));

  a.resize(4, 4);
  EXPECT_EQ(4, a.getNumOfRows());
  EXPECT_EQ(4, a.getNumOfCols());
  EXPECT_DOUBLE_EQ(0.0, a.get(0, 0));
  EXPECT_DOUBLE_EQ(2.0, a.get(0, 1));
  EXPECT_DOUBLE_EQ(0.0, a.get(0, 2));
  EXPECT_DOUBLE_EQ(0.0, a.get(0, 3));
  EXPECT_DOUBLE_EQ(0.0, a.get(1, 0));
  EXPECT_DOUBLE_EQ(0.0, a.get(1, 1));
  EXPECT_DOUBLE_EQ(-4.0, a.get(1, 2));
  EXPECT_DOUBLE_EQ(0.0, a.get(1, 3));
  EXPECT_DOUBLE_EQ(0.0, a.get(2, 0));
  EXPECT_DOUBLE_EQ(0.0, a.get(2, 1));
  EXPECT_DOUBLE_EQ(0.0, a.get(2, 2));
  EXPECT_DOUBLE_EQ(0.0, a.get(2, 3));
  EXPECT_DOUBLE_EQ(0.0, a.get(3, 0));
  EXPECT_DOUBLE_EQ(0.0, a.get(3, 1));
  EXPECT_DOUBLE_EQ(0.0, a.get(3, 2));
  EXPECT_DOUBLE_EQ(0.0, a.get(3, 3));

  a.resize(2, 2);
  EXPECT_EQ(2, a.getNumOfRows());
  EXPECT_EQ(2, a.getNumOfCols());
  EXPECT_DOUBLE_EQ(0.0, a.get(0, 0));
  EXPECT_DOUBLE_EQ(2.0, a.get(0, 1));
  EXPECT_DOUBLE_EQ(0.0, a.get(1, 0));
  EXPECT_DOUBLE_EQ(0.0, a.get(1, 1));
}

TEST(TlDenseGeneralMatrix_Eigen_Old, operator_eq) {
  TlDenseGeneralMatrix_Eigen_Old a = getMatrixA<TlDenseGeneralMatrix_Eigen_Old>();
  TlDenseGeneralMatrix_Eigen_Old c;

  c = a;

  EXPECT_EQ(3, c.getNumOfRows());
  EXPECT_EQ(3, c.getNumOfCols());
  EXPECT_DOUBLE_EQ(0.0, c.get(0, 0));
  EXPECT_DOUBLE_EQ(1.0, c.get(0, 1));
  EXPECT_DOUBLE_EQ(2.0, c.get(0, 2));
  EXPECT_DOUBLE_EQ(3.0, c.get(1, 0));
  EXPECT_DOUBLE_EQ(4.0, c.get(1, 1));
  EXPECT_DOUBLE_EQ(5.0, c.get(1, 2));
  EXPECT_DOUBLE_EQ(6.0, c.get(2, 0));
  EXPECT_DOUBLE_EQ(7.0, c.get(2, 1));
  EXPECT_DOUBLE_EQ(8.0, c.get(2, 2));
}

TEST(TlDenseGeneralMatrix_Eigen_Old, operator_add) {
  TlDenseGeneralMatrix_Eigen_Old a = getMatrixA<TlDenseGeneralMatrix_Eigen_Old>();
  TlDenseGeneralMatrix_Eigen_Old b = getMatrixB<TlDenseGeneralMatrix_Eigen_Old>();

  TlDenseGeneralMatrix_Eigen_Old c = a + b;

  EXPECT_EQ(3, c.getNumOfRows());
  EXPECT_EQ(3, c.getNumOfCols());
  EXPECT_DOUBLE_EQ(0.0, c.get(0, 0));
  EXPECT_DOUBLE_EQ(4.0, c.get(0, 1));
  EXPECT_DOUBLE_EQ(8.0, c.get(0, 2));
  EXPECT_DOUBLE_EQ(4.0, c.get(1, 0));
  EXPECT_DOUBLE_EQ(8.0, c.get(1, 1));
  EXPECT_DOUBLE_EQ(12.0, c.get(1, 2));
  EXPECT_DOUBLE_EQ(8.0, c.get(2, 0));
  EXPECT_DOUBLE_EQ(12.0, c.get(2, 1));
  EXPECT_DOUBLE_EQ(16.0, c.get(2, 2));
}

TEST(TlDenseGeneralMatrix_Eigen_Old, operator_iadd) {
  TlDenseGeneralMatrix_Eigen_Old a = getMatrixA<TlDenseGeneralMatrix_Eigen_Old>();
  TlDenseGeneralMatrix_Eigen_Old b = getMatrixB<TlDenseGeneralMatrix_Eigen_Old>();

  b += a;

  EXPECT_EQ(3, b.getNumOfRows());
  EXPECT_EQ(3, b.getNumOfCols());
  EXPECT_DOUBLE_EQ(0.0, b.get(0, 0));
  EXPECT_DOUBLE_EQ(4.0, b.get(0, 1));
  EXPECT_DOUBLE_EQ(8.0, b.get(0, 2));
  EXPECT_DOUBLE_EQ(4.0, b.get(1, 0));
  EXPECT_DOUBLE_EQ(8.0, b.get(1, 1));
  EXPECT_DOUBLE_EQ(12.0, b.get(1, 2));
  EXPECT_DOUBLE_EQ(8.0, b.get(2, 0));
  EXPECT_DOUBLE_EQ(12.0, b.get(2, 1));
  EXPECT_DOUBLE_EQ(16.0, b.get(2, 2));
}

TEST(TlDenseGeneralMatrix_Eigen_Old, save) {
  TlDenseGeneralMatrix_Eigen_Old m = getMatrixA<TlDenseGeneralMatrix_Eigen_Old>();
  m.save(mat_path);

  TlDenseGeneralMatrix_BLAS_old ref;
  ref.load(mat_path);
  EXPECT_EQ(m.getNumOfRows(), ref.getNumOfRows());
  EXPECT_EQ(m.getNumOfCols(), ref.getNumOfCols());
  for (int r = 0; r < m.getNumOfRows(); ++r) {
    for (int c = 0; c < m.getNumOfCols(); ++c) {
      EXPECT_DOUBLE_EQ(m.get(r, c), ref.get(r, c));
    }
  }
}

TEST(TlDenseGeneralMatrix_Eigen_Old, load) {
  TlDenseGeneralMatrix_BLAS_old ref = getMatrixA<TlDenseGeneralMatrix_BLAS_old>();
  ref.save(mat_path);

  TlDenseGeneralMatrix_Eigen_Old a;
  a.load(mat_path);

  EXPECT_EQ(TlMatrixObject::CSFD, a.getType());
  EXPECT_EQ(3, a.getNumOfRows());
  EXPECT_EQ(3, a.getNumOfCols());
  EXPECT_DOUBLE_EQ(0.0, a.get(0, 0));
  EXPECT_DOUBLE_EQ(1.0, a.get(0, 1));
  EXPECT_DOUBLE_EQ(2.0, a.get(0, 2));
  EXPECT_DOUBLE_EQ(3.0, a.get(1, 0));
  EXPECT_DOUBLE_EQ(4.0, a.get(1, 1));
  EXPECT_DOUBLE_EQ(5.0, a.get(1, 2));
  EXPECT_DOUBLE_EQ(6.0, a.get(2, 0));
  EXPECT_DOUBLE_EQ(7.0, a.get(2, 1));
  EXPECT_DOUBLE_EQ(8.0, a.get(2, 2));
}

#ifdef HAVE_HDF5
// TEST(TlDenseGeneralMatrix_Eigen_Old, hdf5) {
//   TlDenseGeneralMatrix_Eigen_Old m = getMatrix_A();
//   m.saveHdf5("temp.mat.h5", "matrix_A");
//
//   TlDenseGeneralMatrix_Eigen_Old a;
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

// TEST(TlDenseGeneralMatrix_Eigen_Old, inverse)
//
// {
//   TlDenseGeneralMatrix_Eigen_Old a = getMatrix_C();
//   TlDenseGeneralMatrix_Eigen_Old b = a;
//
//   b.inverse();
//
//   TlDenseGeneralMatrix_Eigen_Old c = a * b;
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
TEST(TlDenseGeneralMatrix_Eigen_Old, operator_mul_AB) {
  TlDenseGeneralMatrix_Eigen_Old a = getMatrixA<TlDenseGeneralMatrix_Eigen_Old>();
  TlDenseGeneralMatrix_Eigen_Old b = getMatrixB<TlDenseGeneralMatrix_Eigen_Old>();
  TlDenseGeneralMatrix_Eigen_Old c = a * b;

  EXPECT_EQ(3, c.getNumOfRows());
  EXPECT_EQ(3, c.getNumOfCols());
  EXPECT_DOUBLE_EQ(5.0, c.get(0, 0));
  EXPECT_DOUBLE_EQ(14.0, c.get(0, 1));
  EXPECT_DOUBLE_EQ(23.0, c.get(0, 2));
  EXPECT_DOUBLE_EQ(14.0, c.get(1, 0));
  EXPECT_DOUBLE_EQ(50.0, c.get(1, 1));
  EXPECT_DOUBLE_EQ(86.0, c.get(1, 2));
  EXPECT_DOUBLE_EQ(23.0, c.get(2, 0));
  EXPECT_DOUBLE_EQ(86.0, c.get(2, 1));
  EXPECT_DOUBLE_EQ(149.0, c.get(2, 2));
}

// TEST(TlDenseGeneralMatrix_Eigen_Old, operator_mul_AX) {
//   TlDenseGeneralMatrix_Eigen_Old a = getMatrixA<TlDenseGeneralMatrix_Eigen_Old>();
//   TlVector x(3);
//   x[0] = 1.0;
//   x[1] = 2.0;
//   x[2] = 3.0;
//   TlVector z = a * x;
//
//   EXPECT_EQ(3, (int)z.getSize());
//   EXPECT_DOUBLE_EQ(8.0, z[0]);
//   EXPECT_DOUBLE_EQ(26.0, z[1]);
//   EXPECT_DOUBLE_EQ(44.0, z[2]);
// }

// TEST(TlDenseGeneralMatrix_Eigen_Old, operator_mul_X) {
//   TlDenseGeneralMatrix_Eigen_Old a = getMatrix_A();
//   TlVector x(3);
//   x[0] = 1.0;
//   x[1] = 2.0;
//   x[2] = 3.0;
//   TlVector z = x * a;
//
//   EXPECT_EQ(3, (int)z.getSize());
//   EXPECT_DOUBLE_EQ(24.0, z[0]);
//   EXPECT_DOUBLE_EQ(30.0, z[1]);
//   EXPECT_DOUBLE_EQ(36.0, z[2]);
// }
//
// TEST(TlDenseGeneralMatrix_Eigen_Old, dot) {
//   TlDenseGeneralMatrix_Eigen_Old a = getMatrix_A();
//   TlDenseGeneralMatrix_Eigen_Old b = getMatrix_B();
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
// TEST(TlDenseGeneralMatrix_Eigen_Old, sum) {
//   TlDenseGeneralMatrix_Eigen_Old a = getMatrix_A();
//   double s = a.sum();
//
//   EXPECT_DOUBLE_EQ(36.0, s);
// }
//
// TEST(TlDenseGeneralMatrix_Eigen_Old, solveLinearLeastSquaresProblem) {
//   TlDenseGeneralMatrix_Eigen_Old A(6, 5);
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
//   TlDenseGeneralMatrix_Eigen_Old B(6, 1);
//   B(0, 0) = 7.4;
//   B(1, 0) = 4.2;
//   B(2, 0) = -8.3;
//   B(3, 0) = 1.8;
//   B(4, 0) = 8.6;
//   B(5, 0) = 2.1;
//
//   TlDenseGeneralMatrix_Eigen_Old X = A.solveLinearLeastSquaresProblem(B);
//   // X.print(std::cout);
//
//   TlDenseGeneralMatrix_Eigen_Old AX = A * X;
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
