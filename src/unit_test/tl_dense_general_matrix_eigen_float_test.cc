#include "tl_dense_general_matrix_eigen_float.h"

#include <iostream>

#include "config.h"
#include "gtest/gtest.h"
#include "matrix_common.h"
#include "tl_dense_general_matrix_eigen.h"
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
TEST(TlDenseGeneralMatrix_EigenFloat, constructByGeneralMatrix) {
    const int row = 20;
    const int col = 30;
    TlDenseGeneralMatrix_Eigen M(row, col);
    M.set(1, 0, 1.0);
    M.set(2, 3, 4.0);
    M.set(5, 8, -10.0);

    TlDenseGeneralMatrix_EigenFloat Mf = M;
    EXPECT_EQ(row, Mf.getNumOfRows());
    EXPECT_EQ(col, Mf.getNumOfCols());
    EXPECT_DOUBLE_EQ(1.0, Mf.get(1, 0));
    EXPECT_DOUBLE_EQ(4.0, Mf.get(2, 3));
    EXPECT_DOUBLE_EQ(-10.0, Mf.get(5, 8));
}

// TEST(TlDenseGeneralMatrix_EigenFloat, constructBySparseGeneralMatrix) {
//     const int row = 20;
//     const int col = 30;
//     TlSparseGeneralMatrix_Eigen SM(row, col);

//     SM.set(1, 0, 1.0);
//     SM.set(2, 3, 4.0);
//     SM.set(5, 8, -10.0);

//     TlDenseGeneralMatrix_EigenFloat DM = SM;
//     EXPECT_EQ(row, DM.getNumOfRows());
//     EXPECT_EQ(col, DM.getNumOfCols());
//     EXPECT_DOUBLE_EQ(1.0, DM.get(1, 0));
//     EXPECT_DOUBLE_EQ(4.0, DM.get(2, 3));
//     EXPECT_DOUBLE_EQ(-10.0, DM.get(5, 8));
// }

// #ifdef HAVE_VIENNACL
// TEST(TlDenseGeneralMatrix_EigenFloat, constructBy_DenseGeneralMatrix_ViennaCL) {
//     const int row = 20;
//     const int col = 30;
//     TlDenseGeneralMatrix_ViennaCL DM_VCL(row, col);

//     DM_VCL.set(1, 0, 1.0);
//     DM_VCL.set(2, 3, 4.0);
//     DM_VCL.set(5, 8, -10.0);

//     TlDenseGeneralMatrix_EigenFloat DM = DM_VCL;
//     EXPECT_EQ(row, DM.getNumOfRows());
//     EXPECT_EQ(col, DM.getNumOfCols());
//     EXPECT_DOUBLE_EQ(1.0, DM.get(1, 0));
//     EXPECT_DOUBLE_EQ(4.0, DM.get(2, 3));
//     EXPECT_DOUBLE_EQ(-10.0, DM.get(5, 8));
// }
// #endif  // HAVE_VIENNACL

// TEST(TlDenseGeneralMatrix_EigenFloat, vtr2mat) {
//     const int row = 3;
//     const int col = 4;
//     const int elements = row * col;
//     std::vector<double> vtr(elements);
//     for (int i = 0; i < elements; ++i) {
//         vtr[i] = i;
//     }

//     TlDenseGeneralMatrix_EigenFloat a(row, col, vtr.data());
//     // std::cout << a << std::endl;

//     EXPECT_EQ(row, a.getNumOfRows());
//     EXPECT_EQ(col, a.getNumOfCols());
//     int i = 0;
//     for (int c = 0; c < col; ++c) {  // col-major
//         for (int r = 0; r < row; ++r) {
//             EXPECT_DOUBLE_EQ(vtr[i], a.get(r, c));
//             ++i;
//         }
//     }
// }

// // [0 1 2]   [0]   [ 5]
// // [3 4 5] x [1] = [14]
// // [6 7 8]   [2]   [23]
// TEST(TlDenseGeneralMatrix_EigenFloat, operator_mul_mat_vec) {
//     TlDenseGeneralMatrix_EigenFloat a = getMatrixA<TlDenseGeneralMatrix_EigenFloat>();
//     TlDenseVector_Eigen v = getVectorA<TlDenseVector_Eigen>();

//     TlDenseVector_Eigen z = a * v;

//     EXPECT_EQ(3, z.getSize());
//     EXPECT_DOUBLE_EQ(5.0, z.get(0));
//     EXPECT_DOUBLE_EQ(14.0, z.get(1));
//     EXPECT_DOUBLE_EQ(23.0, z.get(2));
// }

// //           [ 1  2  3  4]
// // [0 1 2] * [ 5  6  7  8] = [23 26 29 32]
// //           [ 9 10 11 12]
// TEST(TlDenseGeneralMatrix_EigenFloat, operator_mul_vec_mat) {
//     TlDenseGeneralMatrix_EigenFloat a = getMatrixD<TlDenseGeneralMatrix_EigenFloat>();
//     TlDenseVector_Eigen v = getVectorA<TlDenseVector_Eigen>();
//     // std::cout << a << std::endl;
//     // std::cout << v << std::endl;

//     TlDenseVector_Eigen z = v * a;
//     // std::cout << z << std::endl;

//     EXPECT_EQ(4, z.getSize());
//     EXPECT_DOUBLE_EQ(23.0, z.get(0));
//     EXPECT_DOUBLE_EQ(26.0, z.get(1));
//     EXPECT_DOUBLE_EQ(29.0, z.get(2));
//     EXPECT_DOUBLE_EQ(32.0, z.get(3));
// }
