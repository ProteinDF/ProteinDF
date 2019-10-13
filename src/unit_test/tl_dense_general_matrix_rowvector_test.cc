#include <iostream>
#include <limits>
#include "gtest/gtest.h"
#include "tl_dense_general_matrix_arrays_roworiented.h"
#include "tl_dense_general_matrix_lapack.h"

static const double EPS = 1.0E-10;  // std::numeric_limits<double>::epsilon();

TEST(TlDenseGeneralMatrix_arrays_RowOriented, constructer) {
    TlDenseGeneralMatrix_arrays_RowOriented A(100, 600);

    EXPECT_EQ(100, A.getNumOfRows());
    EXPECT_EQ(600, A.getNumOfCols());
    EXPECT_NEAR(0.0, A.get(0, 0), EPS);
}

TEST(TlDenseGeneralMatrix_arrays_RowOriented, constructer2) {
    TlDenseGeneralMatrix_arrays_RowOriented A0(100, 600, 10, 0);
    TlDenseGeneralMatrix_arrays_RowOriented A1(100, 600, 10, 1);

    EXPECT_EQ(100, A0.getNumOfRows());
    EXPECT_EQ(600, A0.getNumOfCols());
    EXPECT_EQ(10, A0.getNumOfSubunits());
    EXPECT_EQ(0, A0.getSubunitID());

    EXPECT_EQ(100, A1.getNumOfRows());
    EXPECT_EQ(600, A1.getNumOfCols());
    EXPECT_EQ(10, A1.getNumOfSubunits());
    EXPECT_EQ(1, A1.getSubunitID());
}

TEST(TlDenseGeneralMatrix_arrays_RowOriented, resize) {
    TlDenseGeneralMatrix_arrays_RowOriented A(100, 100);
    A.resize(200, 100);
    EXPECT_EQ(200, A.getNumOfRows());
    EXPECT_EQ(100, A.getNumOfCols());

    TlDenseGeneralMatrix_arrays_RowOriented B(100, 100);
    B.resize(100, 200);
    EXPECT_EQ(100, B.getNumOfRows());
    EXPECT_EQ(200, B.getNumOfCols());

    TlDenseGeneralMatrix_arrays_RowOriented C(50, 100);
    C.resize(50, 100);
    EXPECT_EQ(50, C.getNumOfRows());
    EXPECT_EQ(100, C.getNumOfCols());

    TlDenseGeneralMatrix_arrays_RowOriented D(100, 50);
    D.resize(100, 50);
    EXPECT_EQ(100, D.getNumOfRows());
    EXPECT_EQ(50, D.getNumOfCols());
}

TEST(TlDenseGeneralMatrix_arrays_RowOriented, contents) {
    const int maxRow = 100;
    const int maxCol = 80;
    TlDenseGeneralMatrix_arrays_RowOriented vecA(maxRow, maxCol);
    TlDenseGeneralMatrix_Lapack matA(maxRow, maxCol);

    // setup
    int count = 0;
    for (int r = 0; r < maxRow; ++r) {
        for (int c = 0; c < maxCol; ++c) {
            double v = double(count);
            matA.set(r, c, v);
            vecA.set(r, c, v);

            ++count;
        }
    }

    // test
    for (int r = 0; r < maxRow; ++r) {
        for (int c = 0; c < maxCol; ++c) {
            EXPECT_NEAR(matA.get(r, c), vecA.get(r, c), EPS);
        }
    }
}

TEST(TlDenseGeneralMatrix_arrays_RowOriented, save_load) {
    const int maxRow = 100;
    const int maxCol = 80;
    TlDenseGeneralMatrix_arrays_RowOriented vecA(maxRow, maxCol);
    TlDenseGeneralMatrix_Lapack matA(maxRow, maxCol);

    // setup
    int count = 0;
    for (int r = 0; r < maxRow; ++r) {
        for (int c = 0; c < maxCol; ++c) {
            double v = double(count);
            matA.set(r, c, v);
            vecA.set(r, c, v);

            ++count;
        }
    }

    // save
    vecA.save("/tmp/vecA.mat");

    // load
    TlDenseGeneralMatrix_arrays_RowOriented vecB;
    vecB.load("/tmp/vecA.mat", 0);

    // test
    for (int r = 0; r < maxRow; ++r) {
        for (int c = 0; c < maxCol; ++c) {
            EXPECT_NEAR(matA.get(r, c), vecB.get(r, c), EPS);
        }
    }
}

TEST(TlDenseGeneralMatrix_arrays_RowOriented, toTlMatrix) {
    const int maxRow = 100;
    const int maxCol = 80;
    TlDenseGeneralMatrix_arrays_RowOriented vecA(maxRow, maxCol);
    TlDenseGeneralMatrix_Lapack matA(maxRow, maxCol);

    // setup
    int count = 0;
    for (int r = 0; r < maxRow; ++r) {
        for (int c = 0; c < maxCol; ++c) {
            double v = double(count);
            matA.set(r, c, v);
            vecA.set(r, c, v);

            ++count;
        }
    }

    TlDenseGeneralMatrix_Lapack matB = vecA.getTlMatrixObject();

    // test
    for (int r = 0; r < maxRow; ++r) {
        for (int c = 0; c < maxCol; ++c) {
            EXPECT_NEAR(matA.get(r, c), matB.get(r, c), EPS);
        }
    }
}
