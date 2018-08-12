#include "tl_sparse_general_matrix_viennacl.h"
#include "tl_dense_general_matrix_impl_viennacl.h"
#include "tl_dense_general_matrix_viennacl.h"
#include "tl_dense_vector_impl_viennacl.h"
#include "tl_dense_vector_viennacl.h"
#include "tl_sparse_general_matrix_impl_viennacl.h"
#include "tl_sparse_symmetric_matrix_impl_viennacl.h"
#include "tl_sparse_symmetric_matrix_viennacl.h"

#ifdef HAVE_EIGEN
#include "tl_sparse_general_matrix_eigen.h"
#include "tl_sparse_general_matrix_impl_eigen.h"
#endif  // HAVE_EIGEN

// ---------------------------------------------------------------------------
// constructor & destructor
// ---------------------------------------------------------------------------
TlSparseGeneralMatrix_ViennaCL::TlSparseGeneralMatrix_ViennaCL(
    const TlMatrixObject::index_type row,
    const TlMatrixObject::index_type col) {
  this->pImpl_ = new TlSparseGeneralMatrix_ImplViennaCL(row, col);
}

TlSparseGeneralMatrix_ViennaCL::TlSparseGeneralMatrix_ViennaCL(
    const TlSparseGeneralMatrix_ViennaCL& rhs) {
  this->pImpl_ = new TlSparseGeneralMatrix_ImplViennaCL(
      *(dynamic_cast<const TlSparseGeneralMatrix_ImplViennaCL*>(rhs.pImpl_)));
}

TlSparseGeneralMatrix_ViennaCL::TlSparseGeneralMatrix_ViennaCL(
    const TlSparseSymmetricMatrix_ViennaCL& rhs) {
  this->pImpl_ = new TlSparseGeneralMatrix_ImplViennaCL(
      *(dynamic_cast<const TlSparseSymmetricMatrix_ImplViennaCL*>(rhs.pImpl_)));
}

TlSparseGeneralMatrix_ViennaCL::TlSparseGeneralMatrix_ViennaCL(
    const TlSparseGeneralMatrix_ImplViennaCL& rhs) {
  this->pImpl_ = new TlSparseGeneralMatrix_ImplViennaCL(rhs);
}

TlSparseGeneralMatrix_ViennaCL::TlSparseGeneralMatrix_ViennaCL(
    const TlDenseGeneralMatrix_ViennaCL& rhs) {
  this->pImpl_ = new TlSparseGeneralMatrix_ImplViennaCL(
      *(dynamic_cast<TlDenseGeneralMatrix_ImplViennaCL*>(rhs.pImpl_)));
}

#ifdef HAVE_VIENNACL
TlSparseGeneralMatrix_ViennaCL::TlSparseGeneralMatrix_ViennaCL(
    const TlSparseGeneralMatrix_Eigen& rhs) {
  this->pImpl_ = new TlSparseGeneralMatrix_ImplViennaCL(
      *(dynamic_cast<TlSparseGeneralMatrix_ImplEigen*>(rhs.pImpl_)));
}
#endif  // HAVE_VIENNACL

TlSparseGeneralMatrix_ViennaCL::~TlSparseGeneralMatrix_ViennaCL() {
  delete this->pImpl_;
  this->pImpl_ = NULL;
}

// ---------------------------------------------------------------------------
// operators
// ---------------------------------------------------------------------------
TlSparseGeneralMatrix_ViennaCL& TlSparseGeneralMatrix_ViennaCL::operator=(
    const TlSparseGeneralMatrix_ViennaCL& rhs) {
  if (this != &rhs) {
    delete this->pImpl_;
    this->pImpl_ = new TlSparseGeneralMatrix_ImplViennaCL(
        *(dynamic_cast<TlSparseGeneralMatrix_ImplViennaCL*>(rhs.pImpl_)));
  }

  return *this;
}

TlSparseGeneralMatrix_ViennaCL& TlSparseGeneralMatrix_ViennaCL::operator+=(
    const TlSparseGeneralMatrix_ViennaCL& rhs) {
  *(dynamic_cast<TlSparseGeneralMatrix_ImplViennaCL*>(this->pImpl_)) +=
      *(dynamic_cast<TlSparseGeneralMatrix_ImplViennaCL*>(rhs.pImpl_));
  return *this;
}

TlSparseGeneralMatrix_ViennaCL& TlSparseGeneralMatrix_ViennaCL::operator-=(
    const TlSparseGeneralMatrix_ViennaCL& rhs) {
  *(dynamic_cast<TlSparseGeneralMatrix_ImplViennaCL*>(this->pImpl_)) -=
      *(dynamic_cast<TlSparseGeneralMatrix_ImplViennaCL*>(rhs.pImpl_));
  return *this;
}

TlSparseGeneralMatrix_ViennaCL& TlSparseGeneralMatrix_ViennaCL::operator*=(
    const double coef) {
  *(dynamic_cast<TlSparseGeneralMatrix_ImplViennaCL*>(this->pImpl_)) *= coef;
  return *this;
}

// ----------------------------------------------------------------------------
// others
// ----------------------------------------------------------------------------
TlSparseGeneralMatrix_ViennaCL operator*(const TlSparseGeneralMatrix_ViennaCL& sm1, const TlSparseGeneralMatrix_ViennaCL& sm2) {
    return TlSparseGeneralMatrix_ViennaCL(
        *(dynamic_cast<TlSparseGeneralMatrix_ImplViennaCL*>(sm1.pImpl_)) * 
    *(dynamic_cast<TlSparseGeneralMatrix_ImplViennaCL*>(sm2.pImpl_)));
}

TlDenseVector_ViennaCL operator*(const TlSparseGeneralMatrix_ViennaCL& mat,
                                 const TlDenseVector_ViennaCL& vtr) {
  return TlDenseVector_ViennaCL(
      *(dynamic_cast<TlSparseGeneralMatrix_ImplViennaCL*>(mat.pImpl_)) *
      *(dynamic_cast<TlDenseVector_ImplViennaCL*>(vtr.pImpl_)));
}
