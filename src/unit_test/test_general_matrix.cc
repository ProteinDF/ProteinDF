#include "dense_general_matrix_test_template.h"
#include "tl_dense_general_matrix_lapack.h"
#include "tl_dense_general_matrix_eigen.h"
#include "tl_dense_general_matrix_viennacl.h"

INSTANTIATE_TYPED_TEST_CASE_P(Lapack, DenseGeneralMatrixTest, TlDenseGeneralMatrix_Lapack);
INSTANTIATE_TYPED_TEST_CASE_P(Eigen, DenseGeneralMatrixTest, TlDenseGeneralMatrix_Eigen);

#ifdef HAVE_VIENNACL
INSTANTIATE_TYPED_TEST_CASE_P(ViennaCL, DenseGeneralMatrixTest, TlDenseGeneralMatrix_ViennaCL);
#endif // HAVE_VIENNACL
