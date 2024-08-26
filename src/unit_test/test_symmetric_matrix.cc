#include "dense_symmetric_matrix_test_template.h"

#include "tl_dense_general_matrix_lapack.h"
#include "tl_dense_symmetric_matrix_lapack.h"
INSTANTIATE_TYPED_TEST_SUITE_P(Lapack, DenseSymmetricMatrixTest, TlDenseSymmetricMatrix_Lapack);

#ifdef HAVE_EIGEN
#include "tl_dense_general_matrix_eigen.h"
#include "tl_dense_general_matrix_eigen_float.h"
#include "tl_dense_symmetric_matrix_eigen.h"
#include "tl_dense_symmetric_matrix_eigen_float.h"
INSTANTIATE_TYPED_TEST_SUITE_P(Eigen, DenseSymmetricMatrixTest, TlDenseSymmetricMatrix_Eigen);
INSTANTIATE_TYPED_TEST_SUITE_P(Eigen_FP32, DenseSymmetricMatrixTest, TlDenseSymmetricMatrix_EigenFloat);
#endif  // HAVE_EIGEN

#ifdef HAVE_VIENNACL
#include "tl_dense_general_matrix_viennacl.h"
#include "tl_dense_general_matrix_viennacl_float.h"
#include "tl_dense_symmetric_matrix_viennacl.h"
#include "tl_dense_symmetric_matrix_viennacl_float.h"
INSTANTIATE_TYPED_TEST_SUITE_P(ViennaCL, DenseSymmetricMatrixTest, TlDenseSymmetricMatrix_ViennaCL);
INSTANTIATE_TYPED_TEST_SUITE_P(ViennaCL_FP32, DenseSymmetricMatrixTest, TlDenseSymmetricMatrix_ViennaCLFloat);
#endif  // HAVE_VIENNACL
