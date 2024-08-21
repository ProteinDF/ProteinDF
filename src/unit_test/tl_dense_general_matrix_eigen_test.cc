#include "tl_dense_general_matrix_eigen.h"

#include <iostream>

#include "config.h"
#include "gtest/gtest.h"
#include "matrix_common.h"
#include "tl_dense_general_matrix_eigen_float.h"
#include "tl_dense_vector_eigen.h"
#include "tl_sparse_general_matrix_eigen.h"
#include "vector_common.h"

#ifdef HAVE_VIENNACL
#include "tl_dense_general_matrix_viennacl.h"
#endif  // HAVE_VIENNACL

static const double EPS = 1.0E-10;  // std::numeric_limits<double>::epsilon();
static const double EPS2 = 1.0E-2;
static const std::string mat_save_path = "temp.gen.eigen.save.mat";
static const std::string mat_load_path = "temp.gen.eigen.load.mat";

// -----------------------------------------------------------------------------
// test
// -----------------------------------------------------------------------------
TEST(TlDenseGeneralMatrix_Eigen, constructByGeneralMatrixFloat) {
    const int row = 20;
    const int col = 30;
    TlDenseGeneralMatrix_EigenFloat Mf(row, col);
    Mf.set(1, 0, 1.0);
    Mf.set(2, 3, 4.0);
    Mf.set(5, 8, -10.0);

    TlDenseGeneralMatrix_Eigen M = Mf;
    EXPECT_EQ(row, M.getNumOfRows());
    EXPECT_EQ(col, M.getNumOfCols());
    EXPECT_DOUBLE_EQ(1.0, M.get(1, 0));
    EXPECT_DOUBLE_EQ(4.0, M.get(2, 3));
    EXPECT_DOUBLE_EQ(-10.0, M.get(5, 8));
}

TEST(TlDenseGeneralMatrix_Eigen, constructBySparseGeneralMatrix) {
    const int row = 20;
    const int col = 30;
    TlSparseGeneralMatrix_Eigen SM(row, col);

    SM.set(1, 0, 1.0);
    SM.set(2, 3, 4.0);
    SM.set(5, 8, -10.0);

    TlDenseGeneralMatrix_Eigen DM = SM;
    EXPECT_EQ(row, DM.getNumOfRows());
    EXPECT_EQ(col, DM.getNumOfCols());
    EXPECT_DOUBLE_EQ(1.0, DM.get(1, 0));
    EXPECT_DOUBLE_EQ(4.0, DM.get(2, 3));
    EXPECT_DOUBLE_EQ(-10.0, DM.get(5, 8));
}

#ifdef HAVE_VIENNACL
TEST(TlDenseGeneralMatrix_Eigen, constructBy_DenseGeneralMatrix_ViennaCL) {
    const int row = 20;
    const int col = 30;
    TlDenseGeneralMatrix_ViennaCL DM_VCL(row, col);

    DM_VCL.set(1, 0, 1.0);
    DM_VCL.set(2, 3, 4.0);
    DM_VCL.set(5, 8, -10.0);

    TlDenseGeneralMatrix_Eigen DM = DM_VCL;
    EXPECT_EQ(row, DM.getNumOfRows());
    EXPECT_EQ(col, DM.getNumOfCols());
    EXPECT_DOUBLE_EQ(1.0, DM.get(1, 0));
    EXPECT_DOUBLE_EQ(4.0, DM.get(2, 3));
    EXPECT_DOUBLE_EQ(-10.0, DM.get(5, 8));
}
#endif  // HAVE_VIENNACL

TEST(TlDenseGeneralMatrix_Eigen, vtr2mat) {
    const int row = 3;
    const int col = 4;
    const int elements = row * col;
    std::vector<double> vtr(elements);
    for (int i = 0; i < elements; ++i) {
        vtr[i] = i;
    }

    TlDenseGeneralMatrix_Eigen a(row, col, vtr.data());
    // std::cout << a << std::endl;

    EXPECT_EQ(row, a.getNumOfRows());
    EXPECT_EQ(col, a.getNumOfCols());
    int i = 0;
    for (int c = 0; c < col; ++c) {  // col-major
        for (int r = 0; r < row; ++r) {
            EXPECT_DOUBLE_EQ(vtr[i], a.get(r, c));
            ++i;
        }
    }
}

// TEST(TlDenseGeneralMatrix_Eigen, constructByTlSerializedData) {
//   // TODO
// }
//
// TEST(TlDenseGeneralMatrix_Eigen, vtr2mat) {
//   // TODO
// }

// TEST(TlDenseGeneralMatrix_Eigen, symmat2mat) {
//   TlDenseSymmetricMatrix_Eigen_Old s(4);
//   s.set(0, 0, 1.0);
//   s.set(1, 2, 2.0);
//   s.set(3, 2, 3.0);
//
//   TlDenseGeneralMatrix_Eigen m(s);
//
//   EXPECT_EQ(TlMatrixObject::CSFD, m.getType());
//   EXPECT_EQ(4, m.getNumOfRows());
//   EXPECT_EQ(4, m.getNumOfCols());
//   EXPECT_DOUBLE_EQ(1.0, m.get(0, 0));
//   EXPECT_DOUBLE_EQ(0.0, m.get(0, 1));
//   EXPECT_DOUBLE_EQ(0.0, m.get(0, 2));
//   EXPECT_DOUBLE_EQ(0.0, m.get(0, 3));
//   EXPECT_DOUBLE_EQ(0.0, m.get(1, 0));
//   EXPECT_DOUBLE_EQ(0.0, m.get(1, 1));
//   EXPECT_DOUBLE_EQ(2.0, m.get(1, 2));
//   EXPECT_DOUBLE_EQ(0.0, m.get(1, 3));
//   EXPECT_DOUBLE_EQ(0.0, m.get(2, 0));
//   EXPECT_DOUBLE_EQ(2.0, m.get(2, 1));
//   EXPECT_DOUBLE_EQ(0.0, m.get(2, 2));
//   EXPECT_DOUBLE_EQ(3.0, m.get(2, 3));
//   EXPECT_DOUBLE_EQ(0.0, m.get(3, 0));
//   EXPECT_DOUBLE_EQ(0.0, m.get(3, 1));
//   EXPECT_DOUBLE_EQ(3.0, m.get(3, 2));
//   EXPECT_DOUBLE_EQ(0.0, m.get(3, 3));
// }

// TEST(TlDenseGeneralMatrix_Eigen, save) {
//   TlDenseGeneralMatrix_Eigen m =
//   getMatrixA<TlDenseGeneralMatrix_Eigen>();
//   EXPECT_EQ(3, m.getNumOfRows());
//   EXPECT_EQ(3, m.getNumOfCols());
//   const bool isSaved = m.save(mat_save_path);
//   EXPECT_EQ(true, isSaved);
//
//   TlDenseGeneralMatrix_BLAS_old ref;
//   const bool isLoaded = ref.load(mat_save_path);
//   EXPECT_EQ(true, isLoaded);
//
//   EXPECT_EQ(m.getNumOfRows(), ref.getNumOfRows());
//   EXPECT_EQ(m.getNumOfCols(), ref.getNumOfCols());
//   for (int r = 0; r < m.getNumOfRows(); ++r) {
//     for (int c = 0; c < m.getNumOfCols(); ++c) {
//       EXPECT_DOUBLE_EQ(m.get(r, c), ref.get(r, c));
//     }
//   }
// }
//
// TEST(TlDenseGeneralMatrix_Eigen, load) {
//   TlDenseGeneralMatrix_BLAS_old ref =
//   getMatrixA<TlDenseGeneralMatrix_BLAS_old>();
//   const bool isSaved = ref.save(mat_load_path);
//   EXPECT_EQ(true, isSaved);
//
//   TlDenseGeneralMatrix_Eigen a;
//   const bool isLoaded = a.load(mat_load_path);
//   EXPECT_EQ(true, isLoaded);
//
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

#ifdef HAVE_HDF5
// TEST(TlDenseGeneralMatrix_Eigen, hdf5) {
//   TlDenseGeneralMatrix_Eigen m = getMatrix_A();
//   m.saveHdf5("temp.mat.h5", "matrix_A");
//
//   TlDenseGeneralMatrix_Eigen a;
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

// TEST(TlDenseGeneralMatrix_Eigen, inverse)
//
// {
//   TlDenseGeneralMatrix_Eigen a = getMatrix_C();
//   TlDenseGeneralMatrix_Eigen b = a;
//
//   b.inverse();
//
//   TlDenseGeneralMatrix_Eigen c = a * b;
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

// TEST(TlDenseGeneralMatrix_Eigen, operator_mul_AB) {
//   TlDenseGeneralMatrix_Eigen a =
//   getMatrixA<TlDenseGeneralMatrix_Eigen>();
//   TlDenseGeneralMatrix_Eigen b =
//   getMatrixB<TlDenseGeneralMatrix_Eigen>();
//   TlDenseGeneralMatrix_Eigen c = a * b;
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
TEST(TlDenseGeneralMatrix_Eigen, operator_mul_mat_vec) {
    TlDenseGeneralMatrix_Eigen a = getMatrixA<TlDenseGeneralMatrix_Eigen>();
    TlDenseVector_Eigen v = getVectorA<TlDenseVector_Eigen>();

    TlDenseVector_Eigen z = a * v;

    EXPECT_EQ(3, z.getSize());
    EXPECT_DOUBLE_EQ(5.0, z.get(0));
    EXPECT_DOUBLE_EQ(14.0, z.get(1));
    EXPECT_DOUBLE_EQ(23.0, z.get(2));
}

//           [ 1  2  3  4]
// [0 1 2] * [ 5  6  7  8] = [23 26 29 32]
//           [ 9 10 11 12]
TEST(TlDenseGeneralMatrix_Eigen, operator_mul_vec_mat) {
    TlDenseGeneralMatrix_Eigen a = getMatrixD<TlDenseGeneralMatrix_Eigen>();
    TlDenseVector_Eigen v = getVectorA<TlDenseVector_Eigen>();
    // std::cout << a << std::endl;
    // std::cout << v << std::endl;

    TlDenseVector_Eigen z = v * a;
    // std::cout << z << std::endl;

    EXPECT_EQ(4, z.getSize());
    EXPECT_DOUBLE_EQ(23.0, z.get(0));
    EXPECT_DOUBLE_EQ(26.0, z.get(1));
    EXPECT_DOUBLE_EQ(29.0, z.get(2));
    EXPECT_DOUBLE_EQ(32.0, z.get(3));
}

// TEST(TlDenseGeneralMatrix_Eigen, dot) {
//   TlDenseGeneralMatrix_Eigen a = getMatrix_A();
//   TlDenseGeneralMatrix_Eigen b = getMatrix_B();
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
// TEST(TlDenseGeneralMatrix_Eigen, sum) {
//   TlDenseGeneralMatrix_Eigen a = getMatrix_A();
//   double s = a.sum();
//
//   EXPECT_DOUBLE_EQ(36.0, s);
// }
//
// TEST(TlDenseGeneralMatrix_Eigen, solveLinearLeastSquaresProblem) {
//   TlDenseGeneralMatrix_Eigen A(6, 5);
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
//   TlDenseGeneralMatrix_Eigen B(6, 1);
//   B(0, 0) = 7.4;
//   B(1, 0) = 4.2;
//   B(2, 0) = -8.3;
//   B(3, 0) = 1.8;
//   B(4, 0) = 8.6;
//   B(5, 0) = 2.1;
//
//   TlDenseGeneralMatrix_Eigen X = A.solveLinearLeastSquaresProblem(B);
//   // X.print(std::cout);
//
//   TlDenseGeneralMatrix_Eigen AX = A * X;
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
