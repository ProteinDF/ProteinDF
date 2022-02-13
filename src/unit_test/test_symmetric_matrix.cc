#include "dense_symmetric_matrix_test_template.h"
#include "tl_dense_general_matrix_eigen.h"
#include "tl_dense_general_matrix_lapack.h"
#include "tl_dense_general_matrix_viennacl.h"
#include "tl_dense_symmetric_matrix_eigen.h"
#include "tl_dense_symmetric_matrix_lapack.h"
#include "tl_dense_symmetric_matrix_viennacl.h"

INSTANTIATE_TYPED_TEST_SUITE_P(Lapack, DenseSymmetricMatrixTest, TlDenseSymmetricMatrix_Lapack);

#ifdef HAVE_EIGEN
INSTANTIATE_TYPED_TEST_SUITE_P(Eigen, DenseSymmetricMatrixTest, TlDenseSymmetricMatrix_Eigen);
#endif  // HAVE_EIGEN

#ifdef HAVE_VIENNACL
INSTANTIATE_TYPED_TEST_SUITE_P(ViennaCL, DenseSymmetricMatrixTest, TlDenseSymmetricMatrix_ViennaCL);
#endif  // HAVE_VIENNACL
