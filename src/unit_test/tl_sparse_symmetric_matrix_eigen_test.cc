#include "tl_sparse_symmetric_matrix_eigen.h"
#include <limits>
#include "gtest/gtest.h"
#include "tl_dense_general_matrix_eigen.h"
#include "tl_dense_symmetric_matrix_eigen.h"
#include "tl_dense_vector_eigen.h"

static const double EPS = 1.0E-10;  // std::numeric_limits<double>::epsilon();

TEST(TlSparseSymmetricMatrix_Eigen, constructor) {
    TlSparseSymmetricMatrix_Eigen a(3);

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

TEST(TlSparseSymmetricMatrix_Eigen, copyConstructor) {
    TlSparseSymmetricMatrix_Eigen b(5);
    b.set(0, 0, 1.0);
    b.set(1, 0, 2.0);
    b.set(2, 0, 3.0);
    b.set(3, 3, 4.0);

    TlSparseSymmetricMatrix_Eigen a;
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

TEST(TlSparseSymmetricMatrix_Eigen, constructByDenseMat) {
    TlDenseSymmetricMatrix_Eigen DM(20);
    DM.set(1, 0, 1.0);
    DM.set(3, 3, 3.0);
    DM.set(2, 0, -0.5);
    DM.set(4, 5, 1.0E-20);
    DM.set(6, 2, -1.0E-17);
    DM.set(3, 1, 1.0E-15);
    DM.set(7, 1, -1.0E-15);

    TlSparseSymmetricMatrix_Eigen SM(DM);

    EXPECT_EQ(DM.getNumOfRows(), SM.getNumOfRows());
    EXPECT_EQ(DM.getNumOfCols(), SM.getNumOfCols());
    EXPECT_DOUBLE_EQ(1.0, SM.get(1, 0));
    EXPECT_DOUBLE_EQ(3.0, SM.get(3, 3));
    EXPECT_DOUBLE_EQ(0.0, SM.get(4, 5));
    EXPECT_DOUBLE_EQ(0.0, SM.get(6, 2));
    EXPECT_DOUBLE_EQ(1.0E-15, SM.get(3, 1));
    EXPECT_DOUBLE_EQ(-1.0E-15, SM.get(7, 1));

    EXPECT_DOUBLE_EQ(1.0, SM.get(0, 1));
    EXPECT_DOUBLE_EQ(0.0, SM.get(5, 4));
    EXPECT_DOUBLE_EQ(0.0, SM.get(2, 6));
    EXPECT_DOUBLE_EQ(1.0E-15, SM.get(1, 3));
    EXPECT_DOUBLE_EQ(-1.0E-15, SM.get(1, 7));
}

TEST(TlSparseSymmetricMatrix_Eigen, operator_eq) {
    TlSparseSymmetricMatrix_Eigen b(5);
    b.set(0, 0, 1.0);
    b.set(1, 0, 2.0);
    b.set(2, 0, 3.0);
    b.set(3, 3, 4.0);

    TlSparseSymmetricMatrix_Eigen a(10);
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

TEST(TlSparseSymmetricMatrix_Eigen, mul_densegen_sparsesym) {
    const int row1 = 3;
    const int col1 = 5;
    const int row2 = col1;
    const int col2 = row2;

    TlDenseGeneralMatrix_Eigen a(row1, col1);
    {
        int count = 1;
        for (int i = 0; i < row1; ++i) {
            for (int j = 0; j < col1; ++j) {
                a.set(i, j, count);
                ++count;
            }
        }
    }

    TlSparseSymmetricMatrix_Eigen b(row2);
    b.set(0, 0, 1.0);
    b.set(1, 0, 2.0);
    b.set(2, 0, 3.0);
    b.set(3, 3, 4.0);

    TlDenseSymmetricMatrix_Eigen B(row2);
    B.set(0, 0, 1.0);
    B.set(1, 0, 2.0);
    B.set(2, 0, 3.0);
    B.set(3, 3, 4.0);

    TlDenseGeneralMatrix_Eigen c = a * b;
    TlDenseGeneralMatrix_Eigen C = a * B;

    EXPECT_EQ(row1, c.getNumOfRows());
    EXPECT_EQ(col2, c.getNumOfCols());
    for (int i = 0; i < row1; ++i) {
        for (int j = 0; j < col2; ++j) {
            EXPECT_DOUBLE_EQ(C.get(i, j), c.get(i, j));
        }
    }
    // EXPECT_DOUBLE_EQ(14.0, c.get(0, 0));
    // EXPECT_DOUBLE_EQ( 2.0, c.get(0, 1));
    // EXPECT_DOUBLE_EQ( 3.0, c.get(0, 2));
    // EXPECT_DOUBLE_EQ(16.0, c.get(0, 3));
    // EXPECT_DOUBLE_EQ( 0.0, c.get(0, 4));
    // EXPECT_DOUBLE_EQ(44.0, c.get(1, 0));
    // EXPECT_DOUBLE_EQ(12.0, c.get(1, 1));
    // EXPECT_DOUBLE_EQ(18.0, c.get(1, 2));
    // EXPECT_DOUBLE_EQ(36.0, c.get(1, 3));
    // EXPECT_DOUBLE_EQ( 0.0, c.get(1, 4));
    // EXPECT_DOUBLE_EQ(74.0, c.get(2, 0));
    // EXPECT_DOUBLE_EQ(22.0, c.get(2, 1));
    // EXPECT_DOUBLE_EQ(33.0, c.get(2, 2));
    // EXPECT_DOUBLE_EQ(56.0, c.get(2, 3));
    // EXPECT_DOUBLE_EQ( 0.0, c.get(2, 4));
}

TEST(TlSparseSymmetricMatrix_Eigen, mul_sparsesym_densegen) {
    const int row1 = 4;
    const int col1 = row1;
    const int row2 = col1;
    const int col2 = 8;

    TlSparseSymmetricMatrix_Eigen a(row1);
    a.set(0, 0, 1.0);
    a.set(1, 0, 2.0);
    a.set(2, 0, 3.0);
    a.set(3, 3, 4.0);

    TlDenseSymmetricMatrix_Eigen A(row1);
    A.set(0, 0, 1.0);
    A.set(1, 0, 2.0);
    A.set(2, 0, 3.0);
    A.set(3, 3, 4.0);

    TlDenseGeneralMatrix_Eigen b(row2, col2);
    {
        int count = 1;
        for (int i = 0; i < row2; ++i) {
            for (int j = 0; j < col2; ++j) {
                b.set(i, j, count);
                ++count;
            }
        }
    }

    TlDenseGeneralMatrix_Eigen c = a * b;
    TlDenseGeneralMatrix_Eigen C = A * b;

    EXPECT_EQ(row1, c.getNumOfRows());
    EXPECT_EQ(col2, c.getNumOfCols());
    for (int i = 0; i < row1; ++i) {
        for (int j = 0; j < col2; ++j) {
            EXPECT_DOUBLE_EQ(C.get(i, j), c.get(i, j));
        }
    }
}

TEST(TlSparseSymmetricMatrix_Eigen, multi_sparsemat_densevtr) {
    const int row1 = 4;
    const int col1 = row1;
    const int vlen = col1;

    TlSparseSymmetricMatrix_Eigen SM(row1);
    SM.set(0, 0, 1.0);
    SM.set(1, 0, 2.0);
    SM.set(2, 0, 3.0);
    SM.set(3, 3, 4.0);

    TlDenseSymmetricMatrix_Eigen DM(row1);
    DM.set(0, 0, 1.0);
    DM.set(1, 0, 2.0);
    DM.set(2, 0, 3.0);
    DM.set(3, 3, 4.0);

    TlDenseVector_Eigen x(vlen);
    for (int i = 0; i < vlen; ++i) {
        x.set(i, i);
    }

    TlDenseVector_Eigen y1 = SM * x;
    TlDenseVector_Eigen y2 = DM * x;

    EXPECT_EQ(y1.getSize(), row1);
    for (int i = 0; i < row1; ++i) {
        EXPECT_NEAR(y2.get(i), y1.get(i), 1E-5);
    }
}
