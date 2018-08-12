#include <limits>

#include "config.h"
#include "gtest/gtest.h"
#include "matrix_common.h"
#include "tl_dense_general_matrix_viennacl.h"
#include "tl_dense_symmetric_matrix_eigen.h"
#include "tl_dense_symmetric_matrix_viennacl.h"

static const double EPS = 1.0E-10;  // std::numeric_limits<double>::epsilon();
static const std::string mat_save_load_path = "temp.sym.viennacl.save_load.mat";
static const std::string mat_h5 = "temp.sym.h5";

TEST(TlDenseSymmetricMatrix_ViennaCL, copyFromEigen) {
  TlDenseSymmetricMatrix_Eigen a =
      getSymMatrixA<TlDenseSymmetricMatrix_Eigen>();
  TlDenseSymmetricMatrix_ViennaCL b = a;

  EXPECT_EQ(a.getNumOfRows(), b.getNumOfRows());
  EXPECT_EQ(a.getNumOfCols(), b.getNumOfCols());
  for (int i = 0; i < a.getNumOfRows(); ++i) {
    for (int j = 0; j < a.getNumOfCols(); ++j) {
      EXPECT_DOUBLE_EQ(a.get(i, j), b.get(i, j));
    }
  }
}

TEST(TlDenseSymmetricMatrix_ViennaCL, doesSym2gen) {
  TlDenseSymmetricMatrix_ViennaCL a =
      getSymMatrixA<TlDenseSymmetricMatrix_ViennaCL>();
  TlDenseGeneralMatrix_ViennaCL b = a;

  ASSERT_EQ(3, b.getNumOfRows());
  ASSERT_EQ(3, b.getNumOfCols());
  EXPECT_DOUBLE_EQ(0.0, b.get(0, 0));
  EXPECT_DOUBLE_EQ(1.0, b.get(0, 1));
  EXPECT_DOUBLE_EQ(3.0, b.get(0, 2));
  EXPECT_DOUBLE_EQ(1.0, b.get(1, 0));
  EXPECT_DOUBLE_EQ(2.0, b.get(1, 1));
  EXPECT_DOUBLE_EQ(4.0, b.get(1, 2));
  EXPECT_DOUBLE_EQ(3.0, b.get(2, 0));
  EXPECT_DOUBLE_EQ(4.0, b.get(2, 1));
  EXPECT_DOUBLE_EQ(5.0, b.get(2, 2));
}

TEST(TlDenseSymmetricMatrix_ViennaCL, doesGen2sym) {
  TlDenseGeneralMatrix_ViennaCL a = getMatrixA<TlDenseGeneralMatrix_ViennaCL>();
  TlDenseSymmetricMatrix_ViennaCL b = a;

  ASSERT_EQ(3, b.getNumOfRows());
  ASSERT_EQ(3, b.getNumOfCols());
  EXPECT_DOUBLE_EQ(0.0, b.get(0, 0));
  EXPECT_DOUBLE_EQ(1.0, b.get(0, 1));
  EXPECT_DOUBLE_EQ(2.0, b.get(0, 2));
  EXPECT_DOUBLE_EQ(1.0, b.get(1, 0));
  EXPECT_DOUBLE_EQ(4.0, b.get(1, 1));
  EXPECT_DOUBLE_EQ(5.0, b.get(1, 2));
  EXPECT_DOUBLE_EQ(2.0, b.get(2, 0));
  EXPECT_DOUBLE_EQ(5.0, b.get(2, 1));
  EXPECT_DOUBLE_EQ(8.0, b.get(2, 2));
}

// TEST(TlDenseSymmetricMatrix_ViennaCL, convertFromTlVector1) {
//   // b =
//   // { 0 1 2 3 4 5 }
//   TlVector b(6);
//   for (int i = 0; i < 6; ++i) {
//     b[i] = i;
//   }
//
//   TlDenseSymmetricMatrix_ViennaCL A(b, 3);
//
//   // column oriented
//   EXPECT_DOUBLE_EQ(0.0, A(0, 0));
//   EXPECT_DOUBLE_EQ(1.0, A(0, 1));
//   EXPECT_DOUBLE_EQ(2.0, A(1, 1));
//   EXPECT_DOUBLE_EQ(3.0, A(0, 2));
//   EXPECT_DOUBLE_EQ(4.0, A(1, 2));
//   EXPECT_DOUBLE_EQ(5.0, A(2, 2));
// }

// TEST(TlDenseSymmetricMatrix_ViennaCL, convertFromTlVector2) {
//   TlDenseSymmetricMatrix_ViennaCL a = getSymMatrixA();
//
//   TlVector v = a.getVector();
//
//   EXPECT_DOUBLE_EQ(0.0, v[0]);
//   EXPECT_DOUBLE_EQ(1.0, v[1]);
//   EXPECT_DOUBLE_EQ(2.0, v[2]);
//   EXPECT_DOUBLE_EQ(3.0, v[3]);
//   EXPECT_DOUBLE_EQ(4.0, v[4]);
//   EXPECT_DOUBLE_EQ(5.0, v[5]);
//
//   TlDenseSymmetricMatrix_ViennaCL c(v, 3);
//
//   EXPECT_DOUBLE_EQ(0.0, c(0, 0));
//   EXPECT_DOUBLE_EQ(1.0, c(0, 1));
//   EXPECT_DOUBLE_EQ(3.0, c(0, 2));
//   EXPECT_DOUBLE_EQ(1.0, c(1, 0));
//   EXPECT_DOUBLE_EQ(2.0, c(1, 1));
//   EXPECT_DOUBLE_EQ(4.0, c(1, 2));
//   EXPECT_DOUBLE_EQ(3.0, c(2, 0));
//   EXPECT_DOUBLE_EQ(4.0, c(2, 1));
//   EXPECT_DOUBLE_EQ(5.0, c(2, 2));
// }

TEST(TlDenseSymmetricMatrix_ViennaCL, doesOperatorMulSymSym) {
  TlDenseSymmetricMatrix_ViennaCL a =
      getSymMatrixA<TlDenseSymmetricMatrix_ViennaCL>();
  TlDenseSymmetricMatrix_ViennaCL b =
      getSymMatrixB<TlDenseSymmetricMatrix_ViennaCL>();

  TlDenseGeneralMatrix_ViennaCL c = a * b;

  ASSERT_EQ(3, c.getNumOfRows());
  ASSERT_EQ(3, c.getNumOfCols());
  EXPECT_DOUBLE_EQ(7.0, c.get(0, 0));
  EXPECT_DOUBLE_EQ(15.0, c.get(0, 1));
  EXPECT_DOUBLE_EQ(19.0, c.get(0, 2));
  EXPECT_DOUBLE_EQ(10.0, c.get(1, 0));
  EXPECT_DOUBLE_EQ(23.0, c.get(1, 1));
  EXPECT_DOUBLE_EQ(30.0, c.get(1, 2));
  EXPECT_DOUBLE_EQ(14.0, c.get(2, 0));
  EXPECT_DOUBLE_EQ(35.0, c.get(2, 1));
  EXPECT_DOUBLE_EQ(47.0, c.get(2, 2));
}

// TEST(TlDenseSymmetricMatrix_ViennaCL, operator_mul) {
//   TlDenseSymmetricMatrix_ViennaCL a = getSymMatrixA();
//   TlDenseSymmetricMatrix_ViennaCL b = getSymMatrixB();
//
//   TlDenseGeneralMatrix_BLAS_old c = a * b;
//   // c.print(std::cout);
//
//   ASSERT_EQ(3, c.getNumOfRows());
//   ASSERT_EQ(3, c.getNumOfCols());
//   EXPECT_DOUBLE_EQ(7.0, c(0, 0));
//   EXPECT_DOUBLE_EQ(15.0, c(0, 1));
//   EXPECT_DOUBLE_EQ(19.0, c(0, 2));
//   EXPECT_DOUBLE_EQ(10.0, c(1, 0));
//   EXPECT_DOUBLE_EQ(23.0, c(1, 1));
//   EXPECT_DOUBLE_EQ(30.0, c(1, 2));
//   EXPECT_DOUBLE_EQ(14.0, c(2, 0));
//   EXPECT_DOUBLE_EQ(35.0, c(2, 1));
//   EXPECT_DOUBLE_EQ(47.0, c(2, 2));
// }

TEST(TlDenseSymmetricMatrix_ViennaCL, save_and_load) {
  TlDenseSymmetricMatrix_ViennaCL a =
      getSymMatrixA<TlDenseSymmetricMatrix_ViennaCL>();
  bool isSaved = a.save(mat_save_load_path);
  EXPECT_EQ(isSaved, true);

  TlDenseSymmetricMatrix_ViennaCL b;
  bool isLoaded = b.load(mat_save_load_path);
  EXPECT_EQ(isLoaded, true);

  EXPECT_DOUBLE_EQ(0.0, b.get(0, 0));
  EXPECT_DOUBLE_EQ(1.0, b.get(0, 1));
  EXPECT_DOUBLE_EQ(3.0, b.get(0, 2));
  EXPECT_DOUBLE_EQ(1.0, b.get(1, 0));
  EXPECT_DOUBLE_EQ(2.0, b.get(1, 1));
  EXPECT_DOUBLE_EQ(4.0, b.get(1, 2));
  EXPECT_DOUBLE_EQ(3.0, b.get(2, 0));
  EXPECT_DOUBLE_EQ(4.0, b.get(2, 1));
  EXPECT_DOUBLE_EQ(5.0, b.get(2, 2));
}

// #ifdef HAVE_HDF5
// TEST(TlDenseSymmetricMatrix_ViennaCL, save_hdf5) {
//   TlDenseSymmetricMatrix_ViennaCL a = getSymMatrixA();
//   a.saveHdf5(mat_h5, "matrix_A");
// }
//
// TEST(TlDenseSymmetricMatrix_ViennaCL, load_hdf5) {
//   TlDenseSymmetricMatrix_ViennaCL a;
//   a.loadHdf5(mat_h5, "matrix_A");
//
//   EXPECT_DOUBLE_EQ(0.0, a(0, 0));
//   EXPECT_DOUBLE_EQ(1.0, a(0, 1));
//   EXPECT_DOUBLE_EQ(3.0, a(0, 2));
//   EXPECT_DOUBLE_EQ(1.0, a(1, 0));
//   EXPECT_DOUBLE_EQ(2.0, a(1, 1));
//   EXPECT_DOUBLE_EQ(4.0, a(1, 2));
//   EXPECT_DOUBLE_EQ(3.0, a(2, 0));
//   EXPECT_DOUBLE_EQ(4.0, a(2, 1));
//   EXPECT_DOUBLE_EQ(5.0, a(2, 2));
// }
// #endif  // HAVE_HDF5
//
// TEST(TlDenseSymmetricMatrix_ViennaCL, inverse) {
//   TlDenseSymmetricMatrix_ViennaCL a = getSymMatrixA();
//   TlDenseSymmetricMatrix_ViennaCL b = a;
//
//   b.inverse();
//
//   TlDenseGeneralMatrix_BLAS_old c = a * b;
//
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
// TEST(TlDenseSymmetricMatrix_ViennaCL, operator_mul1) {
//   TlDenseSymmetricMatrix_ViennaCL A = getSymMatrixA();
//   // [ 0  -  - ]
//   // [ 1  2  - ]
//   // [ 3  4  5 ]
//
//   TlDenseGeneralMatrix_BLAS_old B(3, 3);
//   B(0, 0) = 0.0;
//   B(0, 1) = 1.0;
//   B(0, 2) = 2.0;
//   B(1, 0) = 3.0;
//   B(1, 1) = 4.0;
//   B(1, 2) = 5.0;
//   B(2, 0) = 6.0;
//   B(2, 1) = 7.0;
//   B(2, 2) = 8.0;
//
//   TlDenseGeneralMatrix_BLAS_old C = A * B;
//
//   EXPECT_DOUBLE_EQ(21.0, C(0, 0));
//   EXPECT_DOUBLE_EQ(25.0, C(0, 1));
//   EXPECT_DOUBLE_EQ(29.0, C(0, 2));
//   EXPECT_DOUBLE_EQ(30.0, C(1, 0));
//   EXPECT_DOUBLE_EQ(37.0, C(1, 1));
//   EXPECT_DOUBLE_EQ(44.0, C(1, 2));
//   EXPECT_DOUBLE_EQ(42.0, C(2, 0));
//   EXPECT_DOUBLE_EQ(54.0, C(2, 1));
//   EXPECT_DOUBLE_EQ(66.0, C(2, 2));
// }
//
// TEST(TlDenseSymmetricMatrix_ViennaCL, operator_multi2) {
//   TlDenseSymmetricMatrix_ViennaCL A = getSymMatrixA();
//   // [ 0  -  - ]
//   // [ 1  2  - ]
//   // [ 3  4  5 ]
//
//   TlDenseGeneralMatrix_BLAS_old B(3, 3);
//   B(0, 0) = 0.0;
//   B(0, 1) = 1.0;
//   B(0, 2) = 2.0;
//   B(1, 0) = 3.0;
//   B(1, 1) = 4.0;
//   B(1, 2) = 5.0;
//   B(2, 0) = 6.0;
//   B(2, 1) = 7.0;
//   B(2, 2) = 8.0;
//
//   TlDenseGeneralMatrix_BLAS_old C = B * A;
//
//   EXPECT_DOUBLE_EQ(7.0, C(0, 0));
//   EXPECT_DOUBLE_EQ(10.0, C(0, 1));
//   EXPECT_DOUBLE_EQ(14.0, C(0, 2));
//   EXPECT_DOUBLE_EQ(19.0, C(1, 0));
//   EXPECT_DOUBLE_EQ(31.0, C(1, 1));
//   EXPECT_DOUBLE_EQ(50.0, C(1, 2));
//   EXPECT_DOUBLE_EQ(31.0, C(2, 0));
//   EXPECT_DOUBLE_EQ(52.0, C(2, 1));
//   EXPECT_DOUBLE_EQ(86.0, C(2, 2));
// }
//
// TEST(TlDenseSymmetricMatrix_ViennaCL, imul1) {
//   TlDenseSymmetricMatrix_ViennaCL A = getSymMatrixA();
//   // [ 0  -  - ]
//   // [ 1  2  - ]
//   // [ 3  4  5 ]
//
//   TlDenseGeneralMatrix_BLAS_old B(3, 3);
//   B(0, 0) = 0.0;
//   B(0, 1) = 1.0;
//   B(0, 2) = 2.0;
//   B(1, 0) = 3.0;
//   B(1, 1) = 4.0;
//   B(1, 2) = 5.0;
//   B(2, 0) = 6.0;
//   B(2, 1) = 7.0;
//   B(2, 2) = 8.0;
//
//   TlDenseGeneralMatrix_BLAS_old C = A;
//   C *= B;
//
//   EXPECT_DOUBLE_EQ(21.0, C(0, 0));
//   EXPECT_DOUBLE_EQ(25.0, C(0, 1));
//   EXPECT_DOUBLE_EQ(29.0, C(0, 2));
//   EXPECT_DOUBLE_EQ(30.0, C(1, 0));
//   EXPECT_DOUBLE_EQ(37.0, C(1, 1));
//   EXPECT_DOUBLE_EQ(44.0, C(1, 2));
//   EXPECT_DOUBLE_EQ(42.0, C(2, 0));
//   EXPECT_DOUBLE_EQ(54.0, C(2, 1));
//   EXPECT_DOUBLE_EQ(66.0, C(2, 2));
// }
//
// TEST(TlDenseSymmetricMatrix_ViennaCL, testMultiEqual2) {
//   TlDenseSymmetricMatrix_ViennaCL A = getSymMatrixA();
//   // [ 0  -  - ]
//   // [ 1  2  - ]
//   // [ 3  4  5 ]
//
//   TlDenseGeneralMatrix_BLAS_old B(3, 3);
//   B(0, 0) = 0.0;
//   B(0, 1) = 1.0;
//   B(0, 2) = 2.0;
//   B(1, 0) = 3.0;
//   B(1, 1) = 4.0;
//   B(1, 2) = 5.0;
//   B(2, 0) = 6.0;
//   B(2, 1) = 7.0;
//   B(2, 2) = 8.0;
//
//   TlDenseGeneralMatrix_BLAS_old C = B;
//   C *= A;
//
//   EXPECT_DOUBLE_EQ(7.0, C(0, 0));
//   EXPECT_DOUBLE_EQ(10.0, C(0, 1));
//   EXPECT_DOUBLE_EQ(14.0, C(0, 2));
//   EXPECT_DOUBLE_EQ(19.0, C(1, 0));
//   EXPECT_DOUBLE_EQ(31.0, C(1, 1));
//   EXPECT_DOUBLE_EQ(50.0, C(1, 2));
//   EXPECT_DOUBLE_EQ(31.0, C(2, 0));
//   EXPECT_DOUBLE_EQ(52.0, C(2, 1));
//   EXPECT_DOUBLE_EQ(86.0, C(2, 2));
// }
//
// TEST(TlDenseSymmetricMatrix_ViennaCL, dot) {
//   TlDenseSymmetricMatrix_ViennaCL A = getSymMatrixA();
//   TlDenseSymmetricMatrix_ViennaCL B = getSymMatrixB();
//   A.dot(B);
//
//   EXPECT_DOUBLE_EQ(0.0, A(0, 0));
//   EXPECT_DOUBLE_EQ(1.0, A(0, 1));
//   EXPECT_DOUBLE_EQ(6.0, A(0, 2));
//   EXPECT_DOUBLE_EQ(1.0, A(1, 0));
//   EXPECT_DOUBLE_EQ(6.0, A(1, 1));
//   EXPECT_DOUBLE_EQ(16.0, A(1, 2));
//   EXPECT_DOUBLE_EQ(6.0, A(2, 0));
//   EXPECT_DOUBLE_EQ(16.0, A(2, 1));
//   EXPECT_DOUBLE_EQ(25.0, A(2, 2));
// }
//
// TEST(TlDenseSymmetricMatrix_ViennaCLTest, sum) {
//   TlDenseSymmetricMatrix_ViennaCL A = getSymMatrixA();
//   double s = A.sum();
//
//   EXPECT_DOUBLE_EQ(23.0, s);
// }
//
// TEST(TlDenseSymmetricMatrix_ViennaCL, choleskyDecomposition) {
//   TlDenseSymmetricMatrix_ViennaCL A = getSymMatrixC();
//   // A.print(std::cout);
//
//   // TlDenseGeneralMatrix_BLAS_old L = A.choleskyFactorization();
//   TlDenseGeneralMatrix_BLAS_old L = A.choleskyFactorization2(1.0E-16);
//   // L.print(std::cout);
//
//   TlDenseGeneralMatrix_BLAS_old Lt = L;
//   Lt.transpose();
//
//   TlDenseGeneralMatrix_BLAS_old LL = L * Lt;
//   // LL.print(std::cout);
//
//   EXPECT_DOUBLE_EQ(A(0, 0), LL(0, 0));
//   EXPECT_DOUBLE_EQ(A(0, 1), LL(0, 1));
//   EXPECT_DOUBLE_EQ(A(0, 2), LL(0, 2));
//   EXPECT_DOUBLE_EQ(A(0, 3), LL(0, 3));
//   EXPECT_DOUBLE_EQ(A(1, 0), LL(1, 0));
//   EXPECT_DOUBLE_EQ(A(1, 1), LL(1, 1));
//   EXPECT_DOUBLE_EQ(A(1, 2), LL(1, 2));
//   EXPECT_DOUBLE_EQ(A(1, 3), LL(1, 3));
//   EXPECT_DOUBLE_EQ(A(2, 0), LL(2, 0));
//   EXPECT_DOUBLE_EQ(A(2, 1), LL(2, 1));
//   EXPECT_DOUBLE_EQ(A(2, 2), LL(2, 2));
//   EXPECT_DOUBLE_EQ(A(2, 3), LL(2, 3));
//   EXPECT_DOUBLE_EQ(A(3, 0), LL(3, 0));
//   EXPECT_DOUBLE_EQ(A(3, 1), LL(3, 1));
//   EXPECT_DOUBLE_EQ(A(3, 2), LL(3, 2));
//   EXPECT_DOUBLE_EQ(A(3, 3), LL(3, 3));
// }
