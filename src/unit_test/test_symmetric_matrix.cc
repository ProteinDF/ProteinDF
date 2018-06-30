#include "dense_symmetric_matrix_test_template.h"
#include "tl_dense_general_matrix_lapack.h"
#include "tl_dense_general_matrix_eigen.h"
#include "tl_dense_general_matrix_viennacl.h"
#include "tl_dense_symmetric_matrix_lapack.h"
#include "tl_dense_symmetric_matrix_eigen.h"
#include "tl_dense_symmetric_matrix_viennacl.h"

INSTANTIATE_TYPED_TEST_CASE_P(Lapack, DenseSymmetricMatrixTest, TlDenseSymmetricMatrix_Lapack);
INSTANTIATE_TYPED_TEST_CASE_P(Eigen, DenseSymmetricMatrixTest, TlDenseSymmetricMatrix_Eigen);

#ifdef HAVE_VIENNACL
INSTANTIATE_TYPED_TEST_CASE_P(ViennaCL, DenseSymmetricMatrixTest, TlDenseSymmetricMatrix_ViennaCL);
#endif // HAVE_VIENNACL
