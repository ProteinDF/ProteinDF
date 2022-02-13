#include "dense_general_matrix_test_template.h"
#include "tl_dense_general_matrix_eigen.h"
#include "tl_dense_general_matrix_lapack.h"
#include "tl_dense_general_matrix_viennacl.h"

INSTANTIATE_TYPED_TEST_SUITE_P(Lapack, DenseGeneralMatrixTest, TlDenseGeneralMatrix_Lapack);

#ifdef HAVE_EIGEN
INSTANTIATE_TYPED_TEST_SUITE_P(Eigen, DenseGeneralMatrixTest, TlDenseGeneralMatrix_Eigen);
#endif  // HAVE_EIGEN

#ifdef HAVE_VIENNACL
INSTANTIATE_TYPED_TEST_SUITE_P(ViennaCL, DenseGeneralMatrixTest, TlDenseGeneralMatrix_ViennaCL);
#endif  // HAVE_VIENNACL
