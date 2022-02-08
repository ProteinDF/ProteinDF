#include "dense_vector_test_template.h"
#include "tl_dense_vector_eigen.h"
#include "tl_dense_vector_lapack.h"
#include "tl_dense_vector_viennacl.h"

INSTANTIATE_TYPED_TEST_SUITE_P(Lapack, DenseVectorTest, TlDenseVector_Lapack);

#ifdef HAVE_EIGEN
INSTANTIATE_TYPED_TEST_SUITE_P(Eigen, DenseVectorTest, TlDenseVector_Eigen);
#endif  // HAVE_EIGEN

#ifdef HAVE_VIENNACL
INSTANTIATE_TYPED_TEST_SUITE_P(ViennaCL, DenseVectorTest, TlDenseVector_ViennaCL);
#endif  // HAVE_VIENNACL
