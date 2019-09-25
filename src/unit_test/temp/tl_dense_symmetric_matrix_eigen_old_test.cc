#include <limits>

#include "config.h"
#include "gtest/gtest.h"
#include "matrix_common.h"
#include "tl_dense_general_matrix_eigen_old.h"
#include "tl_dense_symmetric_matrix_eigen_old.h"

static const double EPS = 1.0E-10;  // std::numeric_limits<double>::epsilon();
static const std::string mat_path = "temp.sym.mat";
static const std::string mat_h5 = "temp.sym.h5";

TEST(TlDenseSymmetricMatrix_Eigen_Old, constructer) {
    TlDenseSymmetricMatrix_Eigen_Old a(3);

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

TEST(TlDenseSymmetricMatrix_Eigen_Old, operaterRoundBracket) {
    TlDenseSymmetricMatrix_Eigen_Old a(3);

    // [ 0  -  - ]
    // [ 1  2  - ]
    // [ 3  4  5 ]
    int t = 0;
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j <= i; ++j) {
            a.set(i, j, t);
            ++t;
        }
    }

    EXPECT_DOUBLE_EQ(0.0, a.get(0, 0));
    EXPECT_DOUBLE_EQ(1.0, a.get(0, 1));
    EXPECT_DOUBLE_EQ(3.0, a.get(0, 2));
    EXPECT_DOUBLE_EQ(1.0, a.get(1, 0));
    EXPECT_DOUBLE_EQ(2.0, a.get(1, 1));
    EXPECT_DOUBLE_EQ(4.0, a.get(1, 2));
    EXPECT_DOUBLE_EQ(3.0, a.get(2, 0));
    EXPECT_DOUBLE_EQ(4.0, a.get(2, 1));
    EXPECT_DOUBLE_EQ(5.0, a.get(2, 2));
}

TEST(TlDenseSymmetricMatrix_Eigen_Old, copyConstructer) {
    TlDenseSymmetricMatrix_Eigen_Old a =
        getSymMatrixA<TlDenseSymmetricMatrix_Eigen_Old>();
    TlDenseSymmetricMatrix_Eigen_Old c(a);

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

// TEST(TlDenseSymmetricMatrix_Eigen_Old, convertFromTlVector1) {
//   // b =
//   // { 0 1 2 3 4 5 }
//   TlVector b(6);
//   for (int i = 0; i < 6; ++i) {
//     b[i] = i;
//   }
//
//   TlDenseSymmetricMatrix_Eigen_Old A(b, 3);
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
// TEST(TlDenseSymmetricMatrix_Eigen_Old, convertFromTlVector2) {
//   TlDenseSymmetricMatrix_Eigen_Old a =
//   getSymMatrixA<TlDenseSymmetricMatrix_Eigen_Old>();
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
//   TlDenseSymmetricMatrix_Eigen_Old c(v, 3);
//
//   EXPECT_DOUBLE_EQ(0.0, c.get(0, 0));
//   EXPECT_DOUBLE_EQ(1.0, c.get(0, 1));
//   EXPECT_DOUBLE_EQ(3.0, c.get(0, 2));
//   EXPECT_DOUBLE_EQ(1.0, c.get(1, 0));
//   EXPECT_DOUBLE_EQ(2.0, c.get(1, 1));
//   EXPECT_DOUBLE_EQ(4.0, c.get(1, 2));
//   EXPECT_DOUBLE_EQ(3.0, c.get(2, 0));
//   EXPECT_DOUBLE_EQ(4.0, c.get(2, 1));
//   EXPECT_DOUBLE_EQ(5.0, c.get(2, 2));
// }
//
TEST(TlDenseSymmetricMatrix_Eigen_Old, operator_eq) {
    TlDenseSymmetricMatrix_Eigen_Old a =
        getSymMatrixA<TlDenseSymmetricMatrix_Eigen_Old>();
    TlDenseSymmetricMatrix_Eigen_Old c;

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

TEST(TlDenseSymmetricMatrix_Eigen_Old, operator_add) {
    TlDenseSymmetricMatrix_Eigen_Old a =
        getSymMatrixA<TlDenseSymmetricMatrix_Eigen_Old>();
    TlDenseSymmetricMatrix_Eigen_Old b =
        getSymMatrixB<TlDenseSymmetricMatrix_Eigen_Old>();

    TlDenseSymmetricMatrix_Eigen_Old c = a + b;

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

TEST(TlDenseSymmetricMatrix_Eigen_Old, operator_iadd) {
    TlDenseSymmetricMatrix_Eigen_Old a =
        getSymMatrixA<TlDenseSymmetricMatrix_Eigen_Old>();
    TlDenseSymmetricMatrix_Eigen_Old b =
        getSymMatrixB<TlDenseSymmetricMatrix_Eigen_Old>();

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

TEST(TlDenseSymmetricMatrix_Eigen_Old, operator_mul) {
    TlDenseSymmetricMatrix_Eigen_Old a =
        getSymMatrixA<TlDenseSymmetricMatrix_Eigen_Old>();
    TlDenseSymmetricMatrix_Eigen_Old b =
        getSymMatrixB<TlDenseSymmetricMatrix_Eigen_Old>();

    TlDenseGeneralMatrix_Eigen_Old c = a * b;
    // c.print(std::cout);

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

TEST(TlDenseSymmetricMatrix_Eigen_Old, save) {
    TlDenseSymmetricMatrix_Eigen_Old a =
        getSymMatrixA<TlDenseSymmetricMatrix_Eigen_Old>();
    a.save(mat_path);
}

TEST(TlDenseSymmetricMatrix_Eigen_Old, load) {
    TlDenseSymmetricMatrix_Eigen_Old a;
    a.load(mat_path);

    EXPECT_DOUBLE_EQ(0.0, a.get(0, 0));
    EXPECT_DOUBLE_EQ(1.0, a.get(0, 1));
    EXPECT_DOUBLE_EQ(3.0, a.get(0, 2));
    EXPECT_DOUBLE_EQ(1.0, a.get(1, 0));
    EXPECT_DOUBLE_EQ(2.0, a.get(1, 1));
    EXPECT_DOUBLE_EQ(4.0, a.get(1, 2));
    EXPECT_DOUBLE_EQ(3.0, a.get(2, 0));
    EXPECT_DOUBLE_EQ(4.0, a.get(2, 1));
    EXPECT_DOUBLE_EQ(5.0, a.get(2, 2));
}

// #ifdef HAVE_HDF5
// TEST(TlDenseSymmetricMatrix_Eigen_Old, save_hdf5) {
//   TlDenseSymmetricMatrix_Eigen_Old a =
//   getSymMatrixA<TlDenseSymmetricMatrix_Eigen_Old>();
//   a.saveHdf5(mat_h5, "matrix_A");
// }
//
// TEST(TlDenseSymmetricMatrix_Eigen_Old, load_hdf5) {
//   TlDenseSymmetricMatrix_Eigen_Old a;
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
TEST(TlDenseSymmetricMatrix_Eigen_Old, inverse) {
    TlDenseSymmetricMatrix_Eigen_Old a =
        getSymMatrixA<TlDenseSymmetricMatrix_Eigen_Old>();
    ASSERT_EQ(3, a.getNumOfRows());
    ASSERT_EQ(3, a.getNumOfCols());
    EXPECT_DOUBLE_EQ(0.0, a.get(0, 0));
    EXPECT_DOUBLE_EQ(1.0, a.get(0, 1));
    EXPECT_DOUBLE_EQ(3.0, a.get(0, 2));
    EXPECT_DOUBLE_EQ(1.0, a.get(1, 0));
    EXPECT_DOUBLE_EQ(2.0, a.get(1, 1));
    EXPECT_DOUBLE_EQ(4.0, a.get(1, 2));
    EXPECT_DOUBLE_EQ(3.0, a.get(2, 0));
    EXPECT_DOUBLE_EQ(4.0, a.get(2, 1));
    EXPECT_DOUBLE_EQ(5.0, a.get(2, 2));

    TlDenseSymmetricMatrix_Eigen_Old b = a.inverse();
    TlDenseGeneralMatrix_Eigen_Old c = a * b;

    for (int i = 0; i < c.getNumOfRows(); ++i) {
        for (int j = 0; j < c.getNumOfCols(); ++j) {
            if (i != j) {
                ASSERT_NEAR(0.0, c.get(i, j), EPS);
            } else {
                ASSERT_NEAR(1.0, c.get(i, i), EPS);
            }
        }
    }
}

TEST(TlDenseSymmetricMatrix_Eigen_Old, operator_mul1) {
    TlDenseSymmetricMatrix_Eigen_Old A =
        getSymMatrixA<TlDenseSymmetricMatrix_Eigen_Old>();
    // [ 0  -  - ]
    // [ 1  2  - ]
    // [ 3  4  5 ]

    TlDenseGeneralMatrix_Eigen_Old B(3, 3);
    B.set(0, 0, 0.0);
    B.set(0, 1, 1.0);
    B.set(0, 2, 2.0);
    B.set(1, 0, 3.0);
    B.set(1, 1, 4.0);
    B.set(1, 2, 5.0);
    B.set(2, 0, 6.0);
    B.set(2, 1, 7.0);
    B.set(2, 2, 8.0);

    TlDenseGeneralMatrix_Eigen_Old C = A * B;

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

TEST(TlDenseSymmetricMatrix_Eigen_Old, operator_multi2) {
    TlDenseSymmetricMatrix_Eigen_Old A =
        getSymMatrixA<TlDenseSymmetricMatrix_Eigen_Old>();
    // [ 0  -  - ]
    // [ 1  2  - ]
    // [ 3  4  5 ]

    TlDenseGeneralMatrix_Eigen_Old B(3, 3);
    B.set(0, 0, 0.0);
    B.set(0, 1, 1.0);
    B.set(0, 2, 2.0);
    B.set(1, 0, 3.0);
    B.set(1, 1, 4.0);
    B.set(1, 2, 5.0);
    B.set(2, 0, 6.0);
    B.set(2, 1, 7.0);
    B.set(2, 2, 8.0);

    TlDenseGeneralMatrix_Eigen_Old C = B * A;

    EXPECT_DOUBLE_EQ(7.0, C.get(0, 0));
    EXPECT_DOUBLE_EQ(10.0, C.get(0, 1));
    EXPECT_DOUBLE_EQ(14.0, C.get(0, 2));
    EXPECT_DOUBLE_EQ(19.0, C.get(1, 0));
    EXPECT_DOUBLE_EQ(31.0, C.get(1, 1));
    EXPECT_DOUBLE_EQ(50.0, C.get(1, 2));
    EXPECT_DOUBLE_EQ(31.0, C.get(2, 0));
    EXPECT_DOUBLE_EQ(52.0, C.get(2, 1));
    EXPECT_DOUBLE_EQ(86.0, C.get(2, 2));
}

TEST(TlDenseSymmetricMatrix_Eigen_Old, imul1) {
    TlDenseSymmetricMatrix_Eigen_Old A =
        getSymMatrixA<TlDenseSymmetricMatrix_Eigen_Old>();
    // [ 0  -  - ]
    // [ 1  2  - ]
    // [ 3  4  5 ]

    TlDenseGeneralMatrix_Eigen_Old B(3, 3);
    B.set(0, 0, 0.0);
    B.set(0, 1, 1.0);
    B.set(0, 2, 2.0);
    B.set(1, 0, 3.0);
    B.set(1, 1, 4.0);
    B.set(1, 2, 5.0);
    B.set(2, 0, 6.0);
    B.set(2, 1, 7.0);
    B.set(2, 2, 8.0);

    TlDenseGeneralMatrix_Eigen_Old C = A;
    C *= B;

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

TEST(TlDenseSymmetricMatrix_Eigen_Old, testMultiEqual2) {
    TlDenseSymmetricMatrix_Eigen_Old A =
        getSymMatrixA<TlDenseSymmetricMatrix_Eigen_Old>();
    // [ 0  -  - ]
    // [ 1  2  - ]
    // [ 3  4  5 ]

    TlDenseGeneralMatrix_Eigen_Old B(3, 3);
    B.set(0, 0, 0.0);
    B.set(0, 1, 1.0);
    B.set(0, 2, 2.0);
    B.set(1, 0, 3.0);
    B.set(1, 1, 4.0);
    B.set(1, 2, 5.0);
    B.set(2, 0, 6.0);
    B.set(2, 1, 7.0);
    B.set(2, 2, 8.0);

    TlDenseGeneralMatrix_Eigen_Old C = B;
    C *= A;

    EXPECT_DOUBLE_EQ(7.0, C.get(0, 0));
    EXPECT_DOUBLE_EQ(10.0, C.get(0, 1));
    EXPECT_DOUBLE_EQ(14.0, C.get(0, 2));
    EXPECT_DOUBLE_EQ(19.0, C.get(1, 0));
    EXPECT_DOUBLE_EQ(31.0, C.get(1, 1));
    EXPECT_DOUBLE_EQ(50.0, C.get(1, 2));
    EXPECT_DOUBLE_EQ(31.0, C.get(2, 0));
    EXPECT_DOUBLE_EQ(52.0, C.get(2, 1));
    EXPECT_DOUBLE_EQ(86.0, C.get(2, 2));
}

TEST(TlDenseSymmetricMatrix_Eigen_Old, dotInPlace) {
    TlDenseSymmetricMatrix_Eigen_Old A =
        getSymMatrixA<TlDenseSymmetricMatrix_Eigen_Old>();
    TlDenseSymmetricMatrix_Eigen_Old B =
        getSymMatrixB<TlDenseSymmetricMatrix_Eigen_Old>();
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

TEST(TlDenseSymmetricMatrix_Eigen_OldTest, sum) {
    TlDenseSymmetricMatrix_Eigen_Old A =
        getSymMatrixA<TlDenseSymmetricMatrix_Eigen_Old>();
    double s = A.sum();

    EXPECT_DOUBLE_EQ(23.0, s);
}

// TEST(TlDenseSymmetricMatrix_Eigen_Old, choleskyDecomposition) {
//   TlDenseSymmetricMatrix_Eigen_Old A =
//   getSymMatrixC<TlDenseSymmetricMatrix_Eigen_Old>();
//   // A.print(std::cout);
//
//   // TlDenseGeneralMatrix_Eigen_Old L = A.choleskyFactorization();
//   TlDenseGeneralMatrix_Eigen_Old L = A.choleskyFactorization2(1.0E-16);
//   // L.print(std::cout);
//
//   TlDenseGeneralMatrix_Eigen_Old Lt = L;
//   Lt.transpose();
//
//   TlDenseGeneralMatrix_Eigen_Old LL = L * Lt;
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
