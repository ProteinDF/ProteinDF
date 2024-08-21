#ifndef DENSE_SYMMETRIC_MATRIX_TEST_TEMPLATE_H
#define DENSE_SYMMETRIC_MATRIX_TEST_TEMPLATE_H

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif  // HAVE_CONFIG_H

#include <string>

#include "gtest/gtest.h"
#include "matrix_common.h"

static const double EPS = 1.0E-5;  // std::numeric_limits<double>::epsilon();

template <typename T>
class DenseSymmetricMatrixTest : public ::testing::Test {};

TYPED_TEST_SUITE_P(DenseSymmetricMatrixTest);

TYPED_TEST_P(DenseSymmetricMatrixTest, doesConstructor) {
    TypeParam a(3);

    ASSERT_EQ(3, a.getNumOfRows());
    ASSERT_EQ(3, a.getNumOfCols());
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

TYPED_TEST_P(DenseSymmetricMatrixTest, doesSetterGetter) {
    TypeParam a(3);

    ASSERT_EQ(3, a.getNumOfRows());
    ASSERT_EQ(3, a.getNumOfCols());
    a.set(0, 0, 1.0);
    a.set(1, 2, 3.0);
    a.set(2, 0, 5.0);
    EXPECT_DOUBLE_EQ(1.0, a.get(0, 0));
    EXPECT_DOUBLE_EQ(0.0, a.get(0, 1));
    EXPECT_DOUBLE_EQ(5.0, a.get(0, 2));
    EXPECT_DOUBLE_EQ(0.0, a.get(1, 0));
    EXPECT_DOUBLE_EQ(0.0, a.get(1, 1));
    EXPECT_DOUBLE_EQ(3.0, a.get(1, 2));
    EXPECT_DOUBLE_EQ(5.0, a.get(2, 0));
    EXPECT_DOUBLE_EQ(3.0, a.get(2, 1));
    EXPECT_DOUBLE_EQ(0.0, a.get(2, 2));
}

TYPED_TEST_P(DenseSymmetricMatrixTest, doesCopyConstructor) {
    TypeParam a = getSymMatrixA<TypeParam>();
    TypeParam c(a);

    ASSERT_EQ(3, c.getNumOfRows());
    ASSERT_EQ(3, c.getNumOfCols());
    EXPECT_DOUBLE_EQ(0.0, c.get(0, 0));
    EXPECT_DOUBLE_EQ(1.0, c.get(0, 1));
    EXPECT_DOUBLE_EQ(3.0, c.get(0, 2));
    EXPECT_DOUBLE_EQ(1.0, c.get(1, 0));
    EXPECT_DOUBLE_EQ(2.0, c.get(1, 1));
    EXPECT_DOUBLE_EQ(4.0, c.get(1, 2));
    EXPECT_DOUBLE_EQ(3.0, c.get(2, 0));
    EXPECT_DOUBLE_EQ(4.0, c.get(2, 1));
    EXPECT_DOUBLE_EQ(5.0, c.get(2, 2));
}

// TEST(TlDenseSymmetricMatrix_Lapack, convertFromTlVector1) {
//   // b =
//   // { 0 1 2 3 4 5 }
//   TlVector b(6);
//   for (int i = 0; i < 6; ++i) {
//     b[i] = i;
//   }
//
//   TlDenseSymmetricMatrix_Lapack A(b, 3);
//
//   // column oriented
//   EXPECT_DOUBLE_EQ(0.0, A(0, 0));
//   EXPECT_DOUBLE_EQ(1.0, A(0, 1));
//   EXPECT_DOUBLE_EQ(2.0, A(1, 1));
//   EXPECT_DOUBLE_EQ(3.0, A(0, 2));
//   EXPECT_DOUBLE_EQ(4.0, A(1, 2));
//   EXPECT_DOUBLE_EQ(5.0, A(2, 2));
// }
//
// TEST(TlDenseSymmetricMatrix_Lapack, convertFromTlVector2) {
//   TlDenseSymmetricMatrix_Lapack a = getSymMatrixA();
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
//   TlDenseSymmetricMatrix_Lapack c(v, 3);
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

TYPED_TEST_P(DenseSymmetricMatrixTest, doesOperatorEq) {
    TypeParam a = getSymMatrixA<TypeParam>();
    TypeParam c;

    c = a;

    ASSERT_EQ(3, c.getNumOfRows());
    ASSERT_EQ(3, c.getNumOfCols());
    EXPECT_DOUBLE_EQ(0.0, c.get(0, 0));
    EXPECT_DOUBLE_EQ(1.0, c.get(0, 1));
    EXPECT_DOUBLE_EQ(3.0, c.get(0, 2));
    EXPECT_DOUBLE_EQ(1.0, c.get(1, 0));
    EXPECT_DOUBLE_EQ(2.0, c.get(1, 1));
    EXPECT_DOUBLE_EQ(4.0, c.get(1, 2));
    EXPECT_DOUBLE_EQ(3.0, c.get(2, 0));
    EXPECT_DOUBLE_EQ(4.0, c.get(2, 1));
    EXPECT_DOUBLE_EQ(5.0, c.get(2, 2));
}

TYPED_TEST_P(DenseSymmetricMatrixTest, doesOperatorAdd) {
    TypeParam a = getSymMatrixA<TypeParam>();
    TypeParam b = getSymMatrixB<TypeParam>();

    TypeParam c = a + b;

    //   0 1 3
    //   1 2 4
    //   3 4 5

    //   0 1 2
    //   1 3 4
    //   2 4 5

    ASSERT_EQ(3, c.getNumOfRows());
    ASSERT_EQ(3, c.getNumOfCols());
    EXPECT_DOUBLE_EQ(0.0, c.get(0, 0));
    EXPECT_DOUBLE_EQ(2.0, c.get(0, 1));
    EXPECT_DOUBLE_EQ(5.0, c.get(0, 2));
    EXPECT_DOUBLE_EQ(2.0, c.get(1, 0));
    EXPECT_DOUBLE_EQ(5.0, c.get(1, 1));
    EXPECT_DOUBLE_EQ(8.0, c.get(1, 2));
    EXPECT_DOUBLE_EQ(5.0, c.get(2, 0));
    EXPECT_DOUBLE_EQ(8.0, c.get(2, 1));
    EXPECT_DOUBLE_EQ(10.0, c.get(2, 2));
}

TYPED_TEST_P(DenseSymmetricMatrixTest, doesOperatorIAdd) {
    TypeParam a = getSymMatrixA<TypeParam>();
    TypeParam b = getSymMatrixB<TypeParam>();

    b += a;

    ASSERT_EQ(3, b.getNumOfRows());
    ASSERT_EQ(3, b.getNumOfCols());
    EXPECT_DOUBLE_EQ(0.0, b.get(0, 0));
    EXPECT_DOUBLE_EQ(2.0, b.get(0, 1));
    EXPECT_DOUBLE_EQ(5.0, b.get(0, 2));
    EXPECT_DOUBLE_EQ(2.0, b.get(1, 0));
    EXPECT_DOUBLE_EQ(5.0, b.get(1, 1));
    EXPECT_DOUBLE_EQ(8.0, b.get(1, 2));
    EXPECT_DOUBLE_EQ(5.0, b.get(2, 0));
    EXPECT_DOUBLE_EQ(8.0, b.get(2, 1));
    EXPECT_DOUBLE_EQ(10.0, b.get(2, 2));
}

TYPED_TEST_P(DenseSymmetricMatrixTest, doesSaveAndLoad) {
    static const std::string mat_save_load_path = "temp.sym.save_load.mat";

    TypeParam a = getSymMatrixA<TypeParam>();
    bool isSaved = a.save(mat_save_load_path);
    EXPECT_EQ(isSaved, true);

    TypeParam b;
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

TYPED_TEST_P(DenseSymmetricMatrixTest, doesSaveAndLoadToHdf5) {
    static const std::string mat_h5_path = "temp.sym.save_load.h5";

#ifdef HAVE_HDF5
    TypeParam m = getSymMatrixA<TypeParam>();
    m.saveHdf5(mat_h5_path, "matrix_A");

    TypeParam a;
    a.loadHdf5(mat_h5_path, "matrix_A");
    EXPECT_EQ(3, a.getNumOfRows());
    EXPECT_EQ(3, a.getNumOfCols());
    EXPECT_DOUBLE_EQ(0.0, a.get(0, 0));
    EXPECT_DOUBLE_EQ(1.0, a.get(0, 1));
    EXPECT_DOUBLE_EQ(3.0, a.get(0, 2));
    EXPECT_DOUBLE_EQ(1.0, a.get(1, 0));
    EXPECT_DOUBLE_EQ(2.0, a.get(1, 1));
    EXPECT_DOUBLE_EQ(4.0, a.get(1, 2));
    EXPECT_DOUBLE_EQ(3.0, a.get(2, 0));
    EXPECT_DOUBLE_EQ(4.0, a.get(2, 1));
    EXPECT_DOUBLE_EQ(5.0, a.get(2, 2));
#else
    std::cerr << "HDF5 is not supported in this build." << std::endl;
#endif  // HAVE_HDF5
}

TYPED_TEST_P(DenseSymmetricMatrixTest, doesInverse) {
    const TypeParam a = getSymMatrixA<TypeParam>();
    const TypeParam b = a.inverse();

    const TypeParam c = a * b;

    EXPECT_NEAR(1.0, c.get(0, 0), EPS);
    EXPECT_NEAR(0.0, c.get(0, 1), EPS);
    EXPECT_NEAR(0.0, c.get(0, 2), EPS);
    EXPECT_NEAR(0.0, c.get(1, 0), EPS);
    EXPECT_NEAR(1.0, c.get(1, 1), EPS);
    EXPECT_NEAR(0.0, c.get(1, 2), EPS);
    EXPECT_NEAR(0.0, c.get(2, 0), EPS);
    EXPECT_NEAR(0.0, c.get(2, 1), EPS);
    EXPECT_NEAR(1.0, c.get(2, 2), EPS);
}

TYPED_TEST_P(DenseSymmetricMatrixTest, doesDotInPlace) {
    TypeParam A = getSymMatrixA<TypeParam>();
    TypeParam B = getSymMatrixB<TypeParam>();
    A.dotInPlace(B);

    EXPECT_DOUBLE_EQ(0.0, A.get(0, 0));
    EXPECT_DOUBLE_EQ(1.0, A.get(0, 1));
    EXPECT_DOUBLE_EQ(6.0, A.get(0, 2));
    EXPECT_DOUBLE_EQ(1.0, A.get(1, 0));
    EXPECT_DOUBLE_EQ(6.0, A.get(1, 1));
    EXPECT_DOUBLE_EQ(16.0, A.get(1, 2));
    EXPECT_DOUBLE_EQ(6.0, A.get(2, 0));
    EXPECT_DOUBLE_EQ(16.0, A.get(2, 1));
    EXPECT_DOUBLE_EQ(25.0, A.get(2, 2));
}

TYPED_TEST_P(DenseSymmetricMatrixTest, doesSum) {
    TypeParam A = getSymMatrixA<TypeParam>();
    const double s = A.sum();

    EXPECT_DOUBLE_EQ(23.0, s);
}

// TEST(TlDenseSymmetricMatrix_Lapack, choleskyDecomposition) {
//   TlDenseSymmetricMatrix_Lapack A = getSymMatrixC();
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

REGISTER_TYPED_TEST_SUITE_P(DenseSymmetricMatrixTest,
                            doesConstructor,
                            doesSetterGetter,
                            doesCopyConstructor, 
                            doesOperatorEq, doesOperatorAdd,
                            doesOperatorIAdd, doesSaveAndLoad,
                            doesSaveAndLoadToHdf5, doesInverse, doesDotInPlace,
                            doesSum);

#endif  // DENSE_SYMMETRIC_MATRIX_TEST_TEMPLATE_H
