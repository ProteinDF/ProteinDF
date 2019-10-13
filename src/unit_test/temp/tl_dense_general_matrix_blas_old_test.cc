#include <iostream>

#include "config.h"
#include "gtest/gtest.h"
#include "matrix_common.h"
#include "tl_dense_general_matrix_blas_old.h"
#include "tl_dense_symmetric_matrix_blas_old.h"

static const double EPS = 1.0E-10;  // std::numeric_limits<double>::epsilon();
static const double EPS2 = 1.0E-2;
static const std::string mat_path = "temp.mat";

// -----------------------------------------------------------------------------
// test
// -----------------------------------------------------------------------------
TEST(TlDenseGeneralMatrix_BLAS_old, constructer) {
    TlDenseGeneralMatrix_BLAS_old a(3, 3);

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

TEST(TlDenseGeneralMatrix_BLAS_old, constructByTlSerializedData) {
    // TODO
}

TEST(TlDenseGeneralMatrix_BLAS_old, vtr2mat) {
    // TODO
}

TEST(TlDenseGeneralMatrix_BLAS_old, copyConstructer) {
    TlDenseGeneralMatrix_BLAS_old a =
        getMatrixA<TlDenseGeneralMatrix_BLAS_old>();
    TlDenseGeneralMatrix_BLAS_old c(a);

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

TEST(TlDenseGeneralMatrix_BLAS_old, symmat2mat) {
    TlDenseSymmetricMatrix_BLAS_Old s(4);
    s.set(0, 0, 1.0);
    s.set(1, 2, 2.0);
    s.set(3, 2, 3.0);

    TlDenseGeneralMatrix_BLAS_old m(s);

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

TEST(TlDenseGeneralMatrix_BLAS_old, resize) {
    TlDenseGeneralMatrix_BLAS_old a(2, 3);
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

TEST(TlDenseGeneralMatrix_BLAS_old, operator_eq) {
    TlDenseGeneralMatrix_BLAS_old a =
        getMatrixA<TlDenseGeneralMatrix_BLAS_old>();
    TlDenseGeneralMatrix_BLAS_old c;

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

TEST(TlDenseGeneralMatrix_BLAS_old, operator_add) {
    TlDenseGeneralMatrix_BLAS_old a =
        getMatrixA<TlDenseGeneralMatrix_BLAS_old>();
    TlDenseGeneralMatrix_BLAS_old b =
        getMatrixB<TlDenseGeneralMatrix_BLAS_old>();

    TlDenseGeneralMatrix_BLAS_old c = a + b;

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

TEST(TlDenseGeneralMatrix_BLAS_old, operator_iadd) {
    TlDenseGeneralMatrix_BLAS_old a =
        getMatrixA<TlDenseGeneralMatrix_BLAS_old>();
    TlDenseGeneralMatrix_BLAS_old b =
        getMatrixB<TlDenseGeneralMatrix_BLAS_old>();

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

TEST(TlDenseGeneralMatrix_BLAS_old, save) {
    TlDenseGeneralMatrix_BLAS_old m =
        getMatrixA<TlDenseGeneralMatrix_BLAS_old>();
    bool ans = m.save(mat_path);
    m.save("dens-gen.mat");
    EXPECT_EQ(ans, true);
}

TEST(TlDenseGeneralMatrix_BLAS_old, load) {
    TlDenseGeneralMatrix_BLAS_old a;
    bool ans = a.load(mat_path);

    EXPECT_EQ(ans, true);
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
TEST(TlDenseGeneralMatrix_BLAS_old, hdf5) {
    TlDenseGeneralMatrix_BLAS_old m =
        getMatrixA<TlDenseGeneralMatrix_BLAS_old>();
    m.saveHdf5("temp.mat.h5", "matrix_A");

    TlDenseGeneralMatrix_BLAS_old a;
    a.loadHdf5("temp.mat.h5", "matrix_A");
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
#endif  // HAVE_HDF5

TEST(TlDenseGeneralMatrix_BLAS_old, inverse)

{
    TlDenseGeneralMatrix_BLAS_old a =
        getMatrixC<TlDenseGeneralMatrix_BLAS_old>();
    TlDenseGeneralMatrix_BLAS_old b = a;

    b.inverse();

    TlDenseGeneralMatrix_BLAS_old c = a * b;
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

TEST(TlDenseGeneralMatrix_BLAS_old, operator_mul_AB) {
    TlDenseGeneralMatrix_BLAS_old a =
        getMatrixA<TlDenseGeneralMatrix_BLAS_old>();
    TlDenseGeneralMatrix_BLAS_old b =
        getMatrixB<TlDenseGeneralMatrix_BLAS_old>();
    TlDenseGeneralMatrix_BLAS_old c = a * b;

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

TEST(TlDenseGeneralMatrix_BLAS_old, operator_mul_AX) {
    TlDenseGeneralMatrix_BLAS_old a =
        getMatrixA<TlDenseGeneralMatrix_BLAS_old>();
    TlVector_BLAS x(3);
    x[0] = 1.0;
    x[1] = 2.0;
    x[2] = 3.0;
    TlVector_BLAS z = a * x;

    EXPECT_EQ(3, (int)z.getSize());
    EXPECT_DOUBLE_EQ(8.0, z[0]);
    EXPECT_DOUBLE_EQ(26.0, z[1]);
    EXPECT_DOUBLE_EQ(44.0, z[2]);
}

TEST(TlDenseGeneralMatrix_BLAS_old, operator_mul_X) {
    TlDenseGeneralMatrix_BLAS_old a =
        getMatrixA<TlDenseGeneralMatrix_BLAS_old>();
    TlVector_BLAS x(3);
    x[0] = 1.0;
    x[1] = 2.0;
    x[2] = 3.0;
    TlVector_BLAS z = x * a;

    EXPECT_EQ(3, (int)z.getSize());
    EXPECT_DOUBLE_EQ(24.0, z[0]);
    EXPECT_DOUBLE_EQ(30.0, z[1]);
    EXPECT_DOUBLE_EQ(36.0, z[2]);
}

TEST(TlDenseGeneralMatrix_BLAS_old, dotInPlace) {
    TlDenseGeneralMatrix_BLAS_old a =
        getMatrixA<TlDenseGeneralMatrix_BLAS_old>();
    TlDenseGeneralMatrix_BLAS_old b =
        getMatrixB<TlDenseGeneralMatrix_BLAS_old>();
    a.dotInPlace(b);

    EXPECT_DOUBLE_EQ(0.0, a.get(0, 0));
    EXPECT_DOUBLE_EQ(3.0, a.get(1, 0));
    EXPECT_DOUBLE_EQ(12.0, a.get(2, 0));
    EXPECT_DOUBLE_EQ(3.0, a.get(0, 1));
    EXPECT_DOUBLE_EQ(16.0, a.get(1, 1));
    EXPECT_DOUBLE_EQ(35.0, a.get(2, 1));
    EXPECT_DOUBLE_EQ(12.0, a.get(0, 2));
    EXPECT_DOUBLE_EQ(35.0, a.get(1, 2));
    EXPECT_DOUBLE_EQ(64.0, a.get(2, 2));
}

TEST(TlDenseGeneralMatrix_BLAS_old, sum) {
    TlDenseGeneralMatrix_BLAS_old a =
        getMatrixA<TlDenseGeneralMatrix_BLAS_old>();
    double s = a.sum();

    EXPECT_DOUBLE_EQ(36.0, s);
}

TEST(TlDenseGeneralMatrix_BLAS_old, solveLinearLeastSquaresProblem) {
    TlDenseGeneralMatrix_BLAS_old A(6, 5);
    A.set(0, 0, -0.09);
    A.set(0, 1, 0.14);
    A.set(0, 2, -0.46);
    A.set(0, 3, 0.68);
    A.set(0, 4, 1.29);
    A.set(1, 0, -1.56);
    A.set(1, 1, 0.20);
    A.set(1, 2, 0.29);
    A.set(1, 3, 1.09);
    A.set(1, 4, 0.51);
    A.set(2, 0, -1.48);
    A.set(2, 1, -0.43);
    A.set(2, 2, 0.89);
    A.set(2, 3, -0.71);
    A.set(2, 4, -0.96);
    A.set(3, 0, -1.09);
    A.set(3, 1, 0.84);
    A.set(3, 2, 0.77);
    A.set(3, 3, 2.11);
    A.set(3, 4, -1.27);
    A.set(4, 0, 0.08);
    A.set(4, 1, 0.55);
    A.set(4, 2, -1.13);
    A.set(4, 3, 0.14);
    A.set(4, 4, 1.74);
    A.set(5, 0, -1.59);
    A.set(5, 1, -0.72);
    A.set(5, 2, 1.06);
    A.set(5, 3, 1.24);
    A.set(5, 4, 0.34);

    TlDenseGeneralMatrix_BLAS_old B(6, 1);
    B.set(0, 0, 7.4);
    B.set(1, 0, 4.2);
    B.set(2, 0, -8.3);
    B.set(3, 0, 1.8);
    B.set(4, 0, 8.6);
    B.set(5, 0, 2.1);

    TlDenseGeneralMatrix_BLAS_old X = A.solveLinearLeastSquaresProblem(B);
    // X.print(std::cout);

    TlDenseGeneralMatrix_BLAS_old AX = A * X;
    // AX.print(std::cout);

    EXPECT_EQ(5, X.getNumOfRows());
    EXPECT_EQ(1, X.getNumOfCols());
    EXPECT_NEAR(B.get(0, 0), AX.get(0, 0), EPS2);
    EXPECT_NEAR(B.get(1, 0), AX.get(1, 0), EPS2);
    EXPECT_NEAR(B.get(2, 0), AX.get(2, 0), EPS2);
    EXPECT_NEAR(B.get(3, 0), AX.get(3, 0), EPS2);
    EXPECT_NEAR(B.get(4, 0), AX.get(4, 0), EPS2);
}
