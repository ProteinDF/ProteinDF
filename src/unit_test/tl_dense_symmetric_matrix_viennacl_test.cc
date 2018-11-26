#include <limits>

#include "config.h"
#include "gtest/gtest.h"
#include "matrix_common.h"
#include "tl_dense_general_matrix_viennacl.h"
#include "tl_dense_symmetric_matrix_eigen.h"
#include "tl_dense_symmetric_matrix_viennacl.h"
#include "tl_dense_vector_viennacl.h"

static const double EPS = 1.0E-10;  // std::numeric_limits<double>::epsilon();
static const std::string mat_save_load_path = "temp.sym.viennacl.save_load.mat";
static const std::string mat_h5 = "temp.sym.h5";

#ifdef HAVE_EIGEN
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
#endif // HAVE_EIGEN

TEST(TlDenseSymmetricMatrix_ViennaCL, vtr2mat) {
  const int dim = 4;
  const int elements = dim * (dim +1) / 2;
  std::vector<double> vtr(elements);
  for (int i = 0; i < elements; ++i) {
    vtr[i] = i;
  }

  TlDenseSymmetricMatrix_ViennaCL a(dim);
  a.vtr2mat(vtr);

  EXPECT_EQ(dim, a.getNumOfRows());
  EXPECT_EQ(dim, a.getNumOfCols());
  int i = 0;
  for (int c = 0; c < dim; ++c) { // col-major
    for (int r = 0; r <= c; ++r) {
      EXPECT_DOUBLE_EQ(vtr[i], a.get(r, c));
      ++i;
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

TEST(TlDenseSymmetricMatrix_ViennaCL, operator_mul1) {
  TlDenseSymmetricMatrix_ViennaCL A = getSymMatrixA<TlDenseSymmetricMatrix_ViennaCL>();
  // [ 0  -  - ]
  // [ 1  2  - ]
  // [ 3  4  5 ]

  TlDenseGeneralMatrix_ViennaCL B(3, 3);
  B.set(0, 0, 0.0);
  B.set(0, 1, 1.0);
  B.set(0, 2, 2.0);
  B.set(1, 0, 3.0);
  B.set(1, 1, 4.0);
  B.set(1, 2, 5.0);
  B.set(2, 0, 6.0);
  B.set(2, 1, 7.0);
  B.set(2, 2, 8.0);

  TlDenseGeneralMatrix_ViennaCL C = A * B;

  EXPECT_DOUBLE_EQ(21.0, C.get(0, 0));
  EXPECT_DOUBLE_EQ(25.0, C.get(0, 1));
  EXPECT_DOUBLE_EQ(29.0, C.get(0, 2));
  EXPECT_DOUBLE_EQ(30.0, C.get(1, 0));
  EXPECT_DOUBLE_EQ(37.0, C.get(1, 1));
  EXPECT_DOUBLE_EQ(44.0, C.get(1, 2));
  EXPECT_DOUBLE_EQ(42.0, C.get(2, 0));
  EXPECT_DOUBLE_EQ(54.0, C.get(2, 1));
  EXPECT_DOUBLE_EQ(66.0, C.get(2, 2));
}


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

// TEST(TlDenseSymmetricMatrix_ViennaCL, eig) {
//   TlDenseSymmetricMatrix_ViennaCL A = getSymMatrixD<TlDenseSymmetricMatrix_ViennaCL>();

//   TlDenseVector_ViennaCL eigVal;
//   TlDenseGeneralMatrix_ViennaCL eigVec;
//   A.eig(&eigVal, &eigVec);

//   eigVal.save("eigval_vcl.vtr");
//   eigVec.save("eigvec_vcl.mat");

//   // check
//   // -2.0531 -0.5146 -0.2943 12.8621
//   // 
//   //   0.7003 -0.5144 -0.2767  0.4103
//   //   0.3592  0.4851  0.6634  0.4422
//   //  -0.1569  0.5420 -0.6504  0.5085
//   //  -0.5965 -0.4543  0.2457  0.6144
//   EXPECT_NEAR(-2.0531, eigVal.get(0), 1.E-4);
//   EXPECT_NEAR(-0.5146, eigVal.get(1), 1.E-4);
//   EXPECT_NEAR(-0.2943, eigVal.get(2), 1.E-4);
//   EXPECT_NEAR(12.8621, eigVal.get(3), 1.E-4);

//   TlDenseSymmetricMatrix_ViennaCL d(eigVal.getSize());
//   for (int i = 0; i < eigVal.getSize(); ++i) {
//     d.set(i, i, eigVal.get(i));
//   }
//   TlDenseGeneralMatrix_ViennaCL lhs = A * eigVec;
//   TlDenseGeneralMatrix_ViennaCL rhs = eigVec * d;

//   for (int i = 0; i < A.getNumOfRows(); ++i) {
//     for (int j = 0; j < A.getNumOfCols(); ++j) {
//       EXPECT_NEAR(lhs.get(i, j), rhs.get(i, j), 1.0E-2);
//     }
//   }
// }

TEST(TlDenseSymmetricMatrix_ViennaCL, eig) {
  int dim = 100;
  TlDenseSymmetricMatrix_ViennaCL A = getSymmetricMatrix<TlDenseSymmetricMatrix_ViennaCL>(dim);
  A.save("eig_A.mat");

  TlDenseVector_ViennaCL eigVal;
  TlDenseGeneralMatrix_ViennaCL eigVec;
  A.eig(&eigVal, &eigVec);

  eigVal.save("eigval_vcl.vtr");
  eigVec.save("eigvec_vcl.mat");

  // check
  TlDenseSymmetricMatrix_ViennaCL d(eigVal.getSize());
  for (int i = 0; i < dim; ++i) {
    d.set(i, i, eigVal.get(i));
  }
  TlDenseGeneralMatrix_ViennaCL lhs = A * eigVec;
  TlDenseGeneralMatrix_ViennaCL rhs = eigVec * d;

  for (int i = 0; i < dim; ++i) {
    for (int j = 0; j < dim; ++j) {
      EXPECT_NEAR(lhs.get(i, j), rhs.get(i, j), 1.0E-2);
    }
  }
}
