#ifndef COMMON_H
#define COMMON_H

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif  // HAvE_CONFIG_H

#ifdef HAVE_LAPACK
#include "tl_dense_general_matrix_lapack.h"
#include "tl_dense_symmetric_matrix_lapack.h"
#include "tl_dense_vector_lapack.h"
// #include "tl_sparse_general_matrix_lapack.h"
#endif  // HAVE_LAPACK

#ifdef HAVE_EIGEN
#include "tl_dense_general_matrix_eigen.h"
#include "tl_dense_symmetric_matrix_eigen.h"
#include "tl_dense_vector_eigen.h"
#include "tl_sparse_general_matrix_eigen.h"
#include "tl_sparse_symmetric_matrix_eigen.h"
#endif  // HAVE_EIGEN

#ifdef HAVE_VIENNACL
#include "tl_dense_general_matrix_viennacl.h"
#include "tl_dense_symmetric_matrix_viennacl.h"
#include "tl_dense_vector_viennacl.h"
#include "tl_sparse_general_matrix_viennacl.h"
#include "tl_sparse_symmetric_matrix_viennacl.h"
#endif  // HAVE_VIENNACL

#endif  // COMMON_H
