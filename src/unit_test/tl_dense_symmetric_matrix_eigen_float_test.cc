#include "tl_dense_symmetric_matrix_eigen_float.h"

#include <limits>

#include "config.h"
#include "gtest/gtest.h"
#include "matrix_common.h"
#include "tl_dense_symmetric_matrix_eigen.h"
#include "tl_dense_vector_eigen.h"

static const double EPS = 1.0E-10;  // std::numeric_limits<double>::epsilon();
static const std::string mat_save_load_path = "temp.sym.eigen.save_load.mat";
static const std::string mat_h5 = "temp.sym.h5";

TEST(TlDenseSymmetricMatrix_EigenFloat, constructBySymmetricMatrix) {
    const int dim = 30;
    TlDenseSymmetricMatrix_Eigen M(dim);
    M.set(1, 0, 1.0);
    M.set(2, 3, 4.0);
    M.set(5, 8, -10.0);

    TlDenseSymmetricMatrix_EigenFloat Mf = M;
    EXPECT_EQ(dim, Mf.getNumOfRows());
    EXPECT_EQ(dim, Mf.getNumOfCols());
    EXPECT_DOUBLE_EQ(1.0, Mf.get(1, 0));
    EXPECT_DOUBLE_EQ(1.0, Mf.get(0, 1));
    EXPECT_DOUBLE_EQ(4.0, Mf.get(2, 3));
    EXPECT_DOUBLE_EQ(4.0, Mf.get(3, 2));
    EXPECT_DOUBLE_EQ(-10.0, Mf.get(5, 8));
    EXPECT_DOUBLE_EQ(-10.0, Mf.get(8, 5));
}

// TEST(TlDenseSymmetricMatrix_EigenFloat, resize) {
//     const int dim1 = 100;
//     const int dim2 = 200;

//     TlDenseSymmetricMatrix_EigenFloat mat(dim1);

//     {
//         int index = 0;
//         for (int r = 0; r < dim1; ++r) {
//             for (int c = 0; c <= r; ++c) {
//                 mat.set(r, c, index);
//                 ++index;
//             }
//         }
//     }

//     mat.resize(dim2);

//     EXPECT_EQ(dim2, mat.getNumOfRows());
//     EXPECT_EQ(dim2, mat.getNumOfCols());
//     {
//         int index = 0;
//         for (int r = 0; r < dim1; ++r) {
//             for (int c = 0; c <= r; ++c) {
//                 EXPECT_DOUBLE_EQ(double(index), mat.get(r, c));
//                 ++index;
//             }
//         }
//     }
// }

// TEST(TlDenseSymmetricMatrix_EigenFloat, vtr2mat) {
//     const int dim = 4;
//     const int elements = dim * (dim + 1) / 2;
//     std::vector<double> vtr(elements);
//     for (int i = 0; i < elements; ++i) {
//         vtr[i] = i;
//     }

//     TlDenseSymmetricMatrix_EigenFloat a(dim);
//     a.vtr2mat(vtr);

//     EXPECT_EQ(dim, a.getNumOfRows());
//     EXPECT_EQ(dim, a.getNumOfCols());
//     int i = 0;
//     for (int c = 0; c < dim; ++c) {  // col-major
//         for (int r = 0; r <= c; ++r) {
//             EXPECT_DOUBLE_EQ(vtr[i], a.get(r, c));
//             ++i;
//         }
//     }
// }

// TEST(TlDenseSymmetricMatrix_EigenFloat, sym2gen) {
//     TlDenseSymmetricMatrix_EigenFloat a =
//         getSymMatrixA<TlDenseSymmetricMatrix_EigenFloat>();
//     TlDenseGeneralMatrix_Eigen b = a;

//     ASSERT_EQ(3, b.getNumOfRows());
//     ASSERT_EQ(3, b.getNumOfCols());
//     EXPECT_DOUBLE_EQ(0.0, b.get(0, 0));
//     EXPECT_DOUBLE_EQ(1.0, b.get(0, 1));
//     EXPECT_DOUBLE_EQ(3.0, b.get(0, 2));
//     EXPECT_DOUBLE_EQ(1.0, b.get(1, 0));
//     EXPECT_DOUBLE_EQ(2.0, b.get(1, 1));
//     EXPECT_DOUBLE_EQ(4.0, b.get(1, 2));
//     EXPECT_DOUBLE_EQ(3.0, b.get(2, 0));
//     EXPECT_DOUBLE_EQ(4.0, b.get(2, 1));
//     EXPECT_DOUBLE_EQ(5.0, b.get(2, 2));
// }

// TEST(TlDenseSymmetricMatrix_EigenFloat, gen2sym) {
//     TlDenseGeneralMatrix_Eigen a = getMatrixA<TlDenseGeneralMatrix_Eigen>();
//     TlDenseSymmetricMatrix_EigenFloat b = a;

//     ASSERT_EQ(3, b.getNumOfRows());
//     ASSERT_EQ(3, b.getNumOfCols());
//     EXPECT_DOUBLE_EQ(0.0, b.get(0, 0));
//     EXPECT_DOUBLE_EQ(1.0, b.get(0, 1));
//     EXPECT_DOUBLE_EQ(2.0, b.get(0, 2));
//     EXPECT_DOUBLE_EQ(1.0, b.get(1, 0));
//     EXPECT_DOUBLE_EQ(4.0, b.get(1, 1));
//     EXPECT_DOUBLE_EQ(5.0, b.get(1, 2));
//     EXPECT_DOUBLE_EQ(2.0, b.get(2, 0));
//     EXPECT_DOUBLE_EQ(5.0, b.get(2, 1));
//     EXPECT_DOUBLE_EQ(8.0, b.get(2, 2));
// }

// TEST(TlDenseSymmetricMatrix_EigenFloat, save_and_load) {
//     TlDenseSymmetricMatrix_EigenFloat a =
//         getSymMatrixA<TlDenseSymmetricMatrix_EigenFloat>();
//     bool isSaved = a.save(mat_save_load_path);
//     EXPECT_EQ(isSaved, true);

//     TlDenseSymmetricMatrix_EigenFloat b;
//     bool isLoaded = b.load(mat_save_load_path);
//     EXPECT_EQ(isLoaded, true);

//     EXPECT_DOUBLE_EQ(0.0, b.get(0, 0));
//     EXPECT_DOUBLE_EQ(1.0, b.get(0, 1));
//     EXPECT_DOUBLE_EQ(3.0, b.get(0, 2));
//     EXPECT_DOUBLE_EQ(1.0, b.get(1, 0));
//     EXPECT_DOUBLE_EQ(2.0, b.get(1, 1));
//     EXPECT_DOUBLE_EQ(4.0, b.get(1, 2));
//     EXPECT_DOUBLE_EQ(3.0, b.get(2, 0));
//     EXPECT_DOUBLE_EQ(4.0, b.get(2, 1));
//     EXPECT_DOUBLE_EQ(5.0, b.get(2, 2));
// }

// TEST(TlDenseSymmetricMatrix_EigenFloat, mul_sym_sym) {
//     TlDenseSymmetricMatrix_EigenFloat a =
//         getSymMatrixA<TlDenseSymmetricMatrix_EigenFloat>();
//     TlDenseSymmetricMatrix_EigenFloat b =
//         getSymMatrixB<TlDenseSymmetricMatrix_EigenFloat>();

//     TlDenseGeneralMatrix_Eigen c = a * b;

//     ASSERT_EQ(3, c.getNumOfRows());
//     ASSERT_EQ(3, c.getNumOfCols());
//     EXPECT_DOUBLE_EQ(7.0, c.get(0, 0));
//     EXPECT_DOUBLE_EQ(15.0, c.get(0, 1));
//     EXPECT_DOUBLE_EQ(19.0, c.get(0, 2));
//     EXPECT_DOUBLE_EQ(10.0, c.get(1, 0));
//     EXPECT_DOUBLE_EQ(23.0, c.get(1, 1));
//     EXPECT_DOUBLE_EQ(30.0, c.get(1, 2));
//     EXPECT_DOUBLE_EQ(14.0, c.get(2, 0));
//     EXPECT_DOUBLE_EQ(35.0, c.get(2, 1));
//     EXPECT_DOUBLE_EQ(47.0, c.get(2, 2));
// }

// TEST(TlDenseSymmetricMatrix_EigenFloat, multi_sym_gen) {
//     TlDenseSymmetricMatrix_EigenFloat A =
//         getSymMatrixA<TlDenseSymmetricMatrix_EigenFloat>();
//     // [ 0  -  - ]
//     // [ 1  2  - ]
//     // [ 3  4  5 ]

//     TlDenseGeneralMatrix_Eigen B(3, 3);
//     B.set(0, 0, 0.0);
//     B.set(0, 1, 1.0);
//     B.set(0, 2, 2.0);
//     B.set(1, 0, 3.0);
//     B.set(1, 1, 4.0);
//     B.set(1, 2, 5.0);
//     B.set(2, 0, 6.0);
//     B.set(2, 1, 7.0);
//     B.set(2, 2, 8.0);

//     TlDenseGeneralMatrix_Eigen C = A * B;

//     EXPECT_DOUBLE_EQ(21.0, C.get(0, 0));
//     EXPECT_DOUBLE_EQ(25.0, C.get(0, 1));
//     EXPECT_DOUBLE_EQ(29.0, C.get(0, 2));
//     EXPECT_DOUBLE_EQ(30.0, C.get(1, 0));
//     EXPECT_DOUBLE_EQ(37.0, C.get(1, 1));
//     EXPECT_DOUBLE_EQ(44.0, C.get(1, 2));
//     EXPECT_DOUBLE_EQ(42.0, C.get(2, 0));
//     EXPECT_DOUBLE_EQ(54.0, C.get(2, 1));
//     EXPECT_DOUBLE_EQ(66.0, C.get(2, 2));
// }

// TEST(TlDenseSymmetricMatrix_EigenFloat, mul_gen_sym) {
//     TlDenseSymmetricMatrix_EigenFloat A =
//         getSymMatrixA<TlDenseSymmetricMatrix_EigenFloat>();
//     // [ 0  -  - ]
//     // [ 1  2  - ]
//     // [ 3  4  5 ]

//     TlDenseGeneralMatrix_Eigen B(3, 3);
//     B.set(0, 0, 0.0);
//     B.set(0, 1, 1.0);
//     B.set(0, 2, 2.0);
//     B.set(1, 0, 3.0);
//     B.set(1, 1, 4.0);
//     B.set(1, 2, 5.0);
//     B.set(2, 0, 6.0);
//     B.set(2, 1, 7.0);
//     B.set(2, 2, 8.0);

//     TlDenseGeneralMatrix_Eigen C = B * A;

//     EXPECT_DOUBLE_EQ(7.0, C.get(0, 0));
//     EXPECT_DOUBLE_EQ(10.0, C.get(0, 1));
//     EXPECT_DOUBLE_EQ(14.0, C.get(0, 2));
//     EXPECT_DOUBLE_EQ(19.0, C.get(1, 0));
//     EXPECT_DOUBLE_EQ(31.0, C.get(1, 1));
//     EXPECT_DOUBLE_EQ(50.0, C.get(1, 2));
//     EXPECT_DOUBLE_EQ(31.0, C.get(2, 0));
//     EXPECT_DOUBLE_EQ(52.0, C.get(2, 1));
//     EXPECT_DOUBLE_EQ(86.0, C.get(2, 2));
// }

// TEST(TlDenseSymmetricMatrix_EigenFloat, eig) {
//     int dim = 100;
//     TlDenseSymmetricMatrix_EigenFloat A =
//         getSymmetricMatrix<TlDenseSymmetricMatrix_EigenFloat>(dim);
//     A.save("eig_A.mat");

//     TlDenseVector_Eigen eigVal;
//     TlDenseGeneralMatrix_Eigen eigVec;
//     A.eig(&eigVal, &eigVec);

//     eigVal.save("eigval_vcl.vtr");
//     eigVec.save("eigvec_vcl.mat");

//     // check
//     TlDenseSymmetricMatrix_EigenFloat d(eigVal.getSize());
//     for (int i = 0; i < dim; ++i) {
//         d.set(i, i, eigVal.get(i));
//     }
//     TlDenseGeneralMatrix_Eigen lhs = A * eigVec;
//     TlDenseGeneralMatrix_Eigen rhs = eigVec * d;

//     for (int i = 0; i < dim; ++i) {
//         for (int j = 0; j < dim; ++j) {
//             EXPECT_NEAR(lhs.get(i, j), rhs.get(i, j), 1.0E-2);
//         }
//     }
// }
