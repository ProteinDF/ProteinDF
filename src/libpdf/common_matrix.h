#ifdef HAVE_EIGEN3
#include "tl_dense_general_matrix_eigen.h"
#include "tl_dense_symmetric_matrix_eigen.h"
#endif // HAVE_EIGEN3

#ifdef HAVE_LAPACK
#include "tl_dense_general_matrix_lapack.h"
#include "tl_dense_symmetric_matrix_lapack.h"
#endif // HAVE_LAPACK

#ifdef HAVE_VIENNACL
#include "tl_dense_general_matrix_viennacl.h"
#include "tl_dense_symmetric_matrix_viennacl.h"
#endif // HAVE_VIENNACL

