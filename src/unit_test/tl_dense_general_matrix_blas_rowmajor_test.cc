#include <iostream>

#include "config.h"
#include "gtest/gtest.h"
#include "matrix_common.h"
#include "tl_dense_symmetric_matrix_blas_old.h"
#include "tl_matrix_blas_row_major.h"

static const double EPS = 1.0E-10;  // std::numeric_limits<double>::epsilon();
static const double EPS2 = 1.0E-2;
static const std::string mat_path = "temp.mat";

// -----------------------------------------------------------------------------
// test
// -----------------------------------------------------------------------------
TEST(TlMatrixBlasRowMajor, constructer) {
  TlMatrixBlasRowMajor a(3, 3);

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

// TEST(TlMatrixBlasRowMajor, constructByTlSerializedData) {
//   // TODO
// }
//
// TEST(TlMatrixBlasRowMajor, vtr2mat) {
//   // TODO
// }
//
// TEST(TlMatrixBlasRowMajor, copyConstructer) {
//   TlMatrixBlasRowMajor a = getMatrixA();
//   TlMatrixBlasRowMajor c(a);
//
//   EXPECT_EQ(3, c.getNumOfRows());
//   EXPECT_EQ(3, c.getNumOfCols());
//   EXPECT_DOUBLE_EQ(0.0, c(0, 0));
//   EXPECT_DOUBLE_EQ(1.0, c(0, 1));
//   EXPECT_DOUBLE_EQ(2.0, c(0, 2));
//   EXPECT_DOUBLE_EQ(3.0, c(1, 0));
//   EXPECT_DOUBLE_EQ(4.0, c(1, 1));
//   EXPECT_DOUBLE_EQ(5.0, c(1, 2));
//   EXPECT_DOUBLE_EQ(6.0, c(2, 0));
//   EXPECT_DOUBLE_EQ(7.0, c(2, 1));
//   EXPECT_DOUBLE_EQ(8.0, c(2, 2));
// }
//
// TEST(TlMatrixBlasRowMajor, symmat2mat) {
//   TlDenseSymmetricMatrix_BLAS_Old s(4);
//   s(0, 0) = 1.0;
//   s(1, 2) = 2.0;
//   s(3, 2) = 3.0;
//
//   TlMatrixBlasRowMajor m(s);
//
//   EXPECT_EQ(TlMatrixObject::CSFD, m.getType());
//   EXPECT_EQ(4, m.getNumOfRows());
//   EXPECT_EQ(4, m.getNumOfCols());
//   EXPECT_DOUBLE_EQ(1.0, m(0, 0));
//   EXPECT_DOUBLE_EQ(0.0, m(0, 1));
//   EXPECT_DOUBLE_EQ(0.0, m(0, 2));
//   EXPECT_DOUBLE_EQ(0.0, m(0, 3));
//   EXPECT_DOUBLE_EQ(0.0, m(1, 0));
//   EXPECT_DOUBLE_EQ(0.0, m(1, 1));
//   EXPECT_DOUBLE_EQ(2.0, m(1, 2));
//   EXPECT_DOUBLE_EQ(0.0, m(1, 3));
//   EXPECT_DOUBLE_EQ(0.0, m(2, 0));
//   EXPECT_DOUBLE_EQ(2.0, m(2, 1));
//   EXPECT_DOUBLE_EQ(0.0, m(2, 2));
//   EXPECT_DOUBLE_EQ(3.0, m(2, 3));
//   EXPECT_DOUBLE_EQ(0.0, m(3, 0));
//   EXPECT_DOUBLE_EQ(0.0, m(3, 1));
//   EXPECT_DOUBLE_EQ(3.0, m(3, 2));
//   EXPECT_DOUBLE_EQ(0.0, m(3, 3));
// }
//
// TEST(TlMatrixBlasRowMajor, resize) {
//   TlMatrixBlasRowMajor a(2, 3);
//   a(0, 1) = 2.0;
//   a(1, 2) = -4.0;
//
//   EXPECT_EQ(2, a.getNumOfRows());
//   EXPECT_EQ(3, a.getNumOfCols());
//   EXPECT_DOUBLE_EQ(0.0, a(0, 0));
//   EXPECT_DOUBLE_EQ(2.0, a(0, 1));
//   EXPECT_DOUBLE_EQ(0.0, a(0, 2));
//   EXPECT_DOUBLE_EQ(0.0, a(1, 0));
//   EXPECT_DOUBLE_EQ(0.0, a(1, 1));
//   EXPECT_DOUBLE_EQ(-4.0, a(1, 2));
//
//   a.resize(4, 4);
//   EXPECT_EQ(4, a.getNumOfRows());
//   EXPECT_EQ(4, a.getNumOfCols());
//   EXPECT_DOUBLE_EQ(0.0, a(0, 0));
//   EXPECT_DOUBLE_EQ(2.0, a(0, 1));
//   EXPECT_DOUBLE_EQ(0.0, a(0, 2));
//   EXPECT_DOUBLE_EQ(0.0, a(0, 3));
//   EXPECT_DOUBLE_EQ(0.0, a(1, 0));
//   EXPECT_DOUBLE_EQ(0.0, a(1, 1));
//   EXPECT_DOUBLE_EQ(-4.0, a(1, 2));
//   EXPECT_DOUBLE_EQ(0.0, a(1, 3));
//   EXPECT_DOUBLE_EQ(0.0, a(2, 0));
//   EXPECT_DOUBLE_EQ(0.0, a(2, 1));
//   EXPECT_DOUBLE_EQ(0.0, a(2, 2));
//   EXPECT_DOUBLE_EQ(0.0, a(2, 3));
//   EXPECT_DOUBLE_EQ(0.0, a(3, 0));
//   EXPECT_DOUBLE_EQ(0.0, a(3, 1));
//   EXPECT_DOUBLE_EQ(0.0, a(3, 2));
//   EXPECT_DOUBLE_EQ(0.0, a(3, 3));
//
//   a.resize(2, 2);
//   EXPECT_EQ(2, a.getNumOfRows());
//   EXPECT_EQ(2, a.getNumOfCols());
//   EXPECT_DOUBLE_EQ(0.0, a(0, 0));
//   EXPECT_DOUBLE_EQ(2.0, a(0, 1));
//   EXPECT_DOUBLE_EQ(0.0, a(1, 0));
//   EXPECT_DOUBLE_EQ(0.0, a(1, 1));
// }
//
// TEST(TlMatrixBlasRowMajor, operator_eq) {
//   TlMatrixBlasRowMajor a = getMatrixA();
//   TlMatrixBlasRowMajor c;
//
//   c = a;
//
//   EXPECT_EQ(3, c.getNumOfRows());
//   EXPECT_EQ(3, c.getNumOfCols());
//   EXPECT_DOUBLE_EQ(0.0, c(0, 0));
//   EXPECT_DOUBLE_EQ(1.0, c(0, 1));
//   EXPECT_DOUBLE_EQ(2.0, c(0, 2));
//   EXPECT_DOUBLE_EQ(3.0, c(1, 0));
//   EXPECT_DOUBLE_EQ(4.0, c(1, 1));
//   EXPECT_DOUBLE_EQ(5.0, c(1, 2));
//   EXPECT_DOUBLE_EQ(6.0, c(2, 0));
//   EXPECT_DOUBLE_EQ(7.0, c(2, 1));
//   EXPECT_DOUBLE_EQ(8.0, c(2, 2));
// }
//
// TEST(TlMatrixBlasRowMajor, operator_add) {
//   TlMatrixBlasRowMajor a = getMatrixA();
//   TlMatrixBlasRowMajor b = getMatrixB();
//
//   TlMatrixBlasRowMajor c = a + b;
//
//   EXPECT_EQ(3, c.getNumOfRows());
//   EXPECT_EQ(3, c.getNumOfCols());
//   EXPECT_DOUBLE_EQ(0.0, c(0, 0));
//   EXPECT_DOUBLE_EQ(4.0, c(0, 1));
//   EXPECT_DOUBLE_EQ(8.0, c(0, 2));
//   EXPECT_DOUBLE_EQ(4.0, c(1, 0));
//   EXPECT_DOUBLE_EQ(8.0, c(1, 1));
//   EXPECT_DOUBLE_EQ(12.0, c(1, 2));
//   EXPECT_DOUBLE_EQ(8.0, c(2, 0));
//   EXPECT_DOUBLE_EQ(12.0, c(2, 1));
//   EXPECT_DOUBLE_EQ(16.0, c(2, 2));
// }
//
// TEST(TlMatrixBlasRowMajor, operator_iadd) {
//   TlMatrixBlasRowMajor a = getMatrixA();
//   TlMatrixBlasRowMajor b = getMatrixB();
//
//   b += a;
//
//   EXPECT_EQ(3, b.getNumOfRows());
//   EXPECT_EQ(3, b.getNumOfCols());
//   EXPECT_DOUBLE_EQ(0.0, b(0, 0));
//   EXPECT_DOUBLE_EQ(4.0, b(0, 1));
//   EXPECT_DOUBLE_EQ(8.0, b(0, 2));
//   EXPECT_DOUBLE_EQ(4.0, b(1, 0));
//   EXPECT_DOUBLE_EQ(8.0, b(1, 1));
//   EXPECT_DOUBLE_EQ(12.0, b(1, 2));
//   EXPECT_DOUBLE_EQ(8.0, b(2, 0));
//   EXPECT_DOUBLE_EQ(12.0, b(2, 1));
//   EXPECT_DOUBLE_EQ(16.0, b(2, 2));
// }

// TEST(TlMatrixBlasRowMajor, save) {
//   TlMatrixBlasRowMajor m = getMatrixA();
//   m.save(mat_path);
// }
//
// TEST(TlMatrixBlasRowMajor, load) {
//   TlMatrixBlasRowMajor a;
//   a.load(mat_path);
//
//   EXPECT_EQ(TlMatrixObject::CSFD, a.getType());
//   EXPECT_EQ(3, a.getNumOfRows());
//   EXPECT_EQ(3, a.getNumOfCols());
//   EXPECT_DOUBLE_EQ(0.0, a(0, 0));
//   EXPECT_DOUBLE_EQ(1.0, a(0, 1));
//   EXPECT_DOUBLE_EQ(2.0, a(0, 2));
//   EXPECT_DOUBLE_EQ(3.0, a(1, 0));
//   EXPECT_DOUBLE_EQ(4.0, a(1, 1));
//   EXPECT_DOUBLE_EQ(5.0, a(1, 2));
//   EXPECT_DOUBLE_EQ(6.0, a(2, 0));
//   EXPECT_DOUBLE_EQ(7.0, a(2, 1));
//   EXPECT_DOUBLE_EQ(8.0, a(2, 2));
// }

// #ifdef HAVE_HDF5
// TEST(TlMatrixBlasRowMajor, hdf5) {
//   TlMatrixBlasRowMajor m = getMatrixA();
//   m.saveHdf5("temp.mat.h5", "matrix_A");
//
//   TlMatrixBlasRowMajor a;
//   a.loadHdf5("temp.mat.h5", "matrix_A");
//   EXPECT_EQ(TlMatrixObject::CSFD, a.getType());
//   EXPECT_EQ(3, a.getNumOfRows());
//   EXPECT_EQ(3, a.getNumOfCols());
//   EXPECT_DOUBLE_EQ(0.0, a(0, 0));
//   EXPECT_DOUBLE_EQ(1.0, a(0, 1));
//   EXPECT_DOUBLE_EQ(2.0, a(0, 2));
//   EXPECT_DOUBLE_EQ(3.0, a(1, 0));
//   EXPECT_DOUBLE_EQ(4.0, a(1, 1));
//   EXPECT_DOUBLE_EQ(5.0, a(1, 2));
//   EXPECT_DOUBLE_EQ(6.0, a(2, 0));
//   EXPECT_DOUBLE_EQ(7.0, a(2, 1));
//   EXPECT_DOUBLE_EQ(8.0, a(2, 2));
// }
// #endif  // HAVE_HDF5
//
// TEST(TlMatrixBlasRowMajor, inverse)
//
// {
//   TlMatrixBlasRowMajor a = getMatrixC();
//   TlMatrixBlasRowMajor b = a;
//
//   b.inverse();
//
//   TlMatrixBlasRowMajor c = a * b;
//   EXPECT_NEAR(1.0, c(0, 0), EPS);
//   EXPECT_NEAR(0.0, c(0, 1), EPS);
//   EXPECT_NEAR(0.0, c(0, 2), EPS);
//   EXPECT_NEAR(0.0, c(1, 0), EPS);
//   EXPECT_NEAR(1.0, c(1, 1), EPS);
//   EXPECT_NEAR(0.0, c(1, 2), EPS);
//   EXPECT_NEAR(0.0, c(2, 0), EPS);
//   EXPECT_NEAR(0.0, c(2, 1), EPS);
//   EXPECT_NEAR(1.0, c(2, 2), EPS);
// }
//
// TEST(TlMatrixBlasRowMajor, operator_mul_AB) {
//   TlMatrixBlasRowMajor a = getMatrixA();
//   TlMatrixBlasRowMajor b = getMatrixB();
//   TlMatrixBlasRowMajor c = a * b;
//
//   EXPECT_EQ(3, c.getNumOfRows());
//   EXPECT_EQ(3, c.getNumOfCols());
//   EXPECT_DOUBLE_EQ(5.0, c(0, 0));
//   EXPECT_DOUBLE_EQ(14.0, c(0, 1));
//   EXPECT_DOUBLE_EQ(23.0, c(0, 2));
//   EXPECT_DOUBLE_EQ(14.0, c(1, 0));
//   EXPECT_DOUBLE_EQ(50.0, c(1, 1));
//   EXPECT_DOUBLE_EQ(86.0, c(1, 2));
//   EXPECT_DOUBLE_EQ(23.0, c(2, 0));
//   EXPECT_DOUBLE_EQ(86.0, c(2, 1));
//   EXPECT_DOUBLE_EQ(149.0, c(2, 2));
// }
//
// TEST(TlMatrixBlasRowMajor, operator_mul_AX) {
//   TlMatrixBlasRowMajor a = getMatrixA();
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
//
// TEST(TlMatrixBlasRowMajor, operator_mul_X) {
//   TlMatrixBlasRowMajor a = getMatrixA();
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
// TEST(TlMatrixBlasRowMajor, dot) {
//   TlMatrixBlasRowMajor a = getMatrixA();
//   TlMatrixBlasRowMajor b = getMatrixB();
//   a.dot(b);
//
//   EXPECT_DOUBLE_EQ(0.0, a(0, 0));
//   EXPECT_DOUBLE_EQ(3.0, a(1, 0));
//   EXPECT_DOUBLE_EQ(12.0, a(2, 0));
//   EXPECT_DOUBLE_EQ(3.0, a(0, 1));
//   EXPECT_DOUBLE_EQ(16.0, a(1, 1));
//   EXPECT_DOUBLE_EQ(35.0, a(2, 1));
//   EXPECT_DOUBLE_EQ(12.0, a(0, 2));
//   EXPECT_DOUBLE_EQ(35.0, a(1, 2));
//   EXPECT_DOUBLE_EQ(64.0, a(2, 2));
// }
//
// TEST(TlMatrixBlasRowMajor, sum) {
//   TlMatrixBlasRowMajor a = getMatrixA();
//   double s = a.sum();
//
//   EXPECT_DOUBLE_EQ(36.0, s);
// }
//
// TEST(TlMatrixBlasRowMajor, solveLinearLeastSquaresProblem) {
//   TlMatrixBlasRowMajor A(6, 5);
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
//   TlMatrixBlasRowMajor B(6, 1);
//   B(0, 0) = 7.4;
//   B(1, 0) = 4.2;
//   B(2, 0) = -8.3;
//   B(3, 0) = 1.8;
//   B(4, 0) = 8.6;
//   B(5, 0) = 2.1;
//
//   TlMatrixBlasRowMajor X = A.solveLinearLeastSquaresProblem(B);
//   // X.print(std::cout);
//
//   TlMatrixBlasRowMajor AX = A * X;
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
