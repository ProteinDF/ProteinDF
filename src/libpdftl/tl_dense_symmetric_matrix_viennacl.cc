#include "tl_dense_symmetric_matrix_viennacl.h"
#include "tl_dense_general_matrix_impl_viennacl.h"
#include "tl_dense_general_matrix_viennacl.h"
#include "tl_dense_symmetric_matrix_impl_viennacl.h"
#include "tl_dense_vector_impl_viennacl.h"
#include "tl_dense_vector_viennacl.h"

TlDenseSymmetricMatrix_ViennaCL::TlDenseSymmetricMatrix_ViennaCL(
    const TlMatrixObject::index_type dim) {
  this->pImpl_ = new TlDenseSymmetricMatrix_ImplViennaCL(dim);
}

TlDenseSymmetricMatrix_ViennaCL::TlDenseSymmetricMatrix_ViennaCL(
    const TlDenseSymmetricMatrix_ViennaCL& rhs) {
  this->pImpl_ = new TlDenseSymmetricMatrix_ImplViennaCL(
      *(dynamic_cast<const TlDenseSymmetricMatrix_ImplViennaCL*>(rhs.pImpl_)));
}

TlDenseSymmetricMatrix_ViennaCL::TlDenseSymmetricMatrix_ViennaCL(
    const TlDenseGeneralMatrix_ViennaCL& rhs) {
  this->pImpl_ = new TlDenseSymmetricMatrix_ImplViennaCL(
      *(dynamic_cast<const TlDenseGeneralMatrix_ImplViennaCL*>(rhs.pImpl_)));
}

TlDenseSymmetricMatrix_ViennaCL::~TlDenseSymmetricMatrix_ViennaCL() {
  delete this->pImpl_;
  this->pImpl_ = NULL;
}

// ---------------------------------------------------------------------------
// operators
// ---------------------------------------------------------------------------
TlDenseSymmetricMatrix_ViennaCL& TlDenseSymmetricMatrix_ViennaCL::operator=(
    const TlDenseSymmetricMatrix_ViennaCL& rhs) {
  delete this->pImpl_;
  this->pImpl_ = new TlDenseSymmetricMatrix_ImplViennaCL(
      *(dynamic_cast<TlDenseSymmetricMatrix_ImplViennaCL*>(rhs.pImpl_)));
}

const TlDenseSymmetricMatrix_ViennaCL TlDenseSymmetricMatrix_ViennaCL::
operator+(const TlDenseSymmetricMatrix_ViennaCL& rhs) const {
  TlDenseSymmetricMatrix_ViennaCL answer = *this;
  answer += rhs;
  return answer;
}

const TlDenseSymmetricMatrix_ViennaCL TlDenseSymmetricMatrix_ViennaCL::
operator-(const TlDenseSymmetricMatrix_ViennaCL& rhs) const {
  TlDenseSymmetricMatrix_ViennaCL answer = *this;
  answer -= rhs;
  return answer;
}

const TlDenseSymmetricMatrix_ViennaCL TlDenseSymmetricMatrix_ViennaCL::
operator*(const TlDenseSymmetricMatrix_ViennaCL& rhs) const {
  TlDenseSymmetricMatrix_ViennaCL answer = *this;
  answer *= rhs;
  return answer;
}

TlDenseSymmetricMatrix_ViennaCL& TlDenseSymmetricMatrix_ViennaCL::operator+=(
    const TlDenseSymmetricMatrix_ViennaCL& rhs) {
  *(dynamic_cast<TlDenseSymmetricMatrix_ImplViennaCL*>(this->pImpl_)) +=
      *(dynamic_cast<TlDenseSymmetricMatrix_ImplViennaCL*>(rhs.pImpl_));
}

TlDenseSymmetricMatrix_ViennaCL& TlDenseSymmetricMatrix_ViennaCL::operator-=(
    const TlDenseSymmetricMatrix_ViennaCL& rhs) {
  *(dynamic_cast<TlDenseSymmetricMatrix_ImplViennaCL*>(this->pImpl_)) -=
      *(dynamic_cast<TlDenseSymmetricMatrix_ImplViennaCL*>(rhs.pImpl_));
}

TlDenseSymmetricMatrix_ViennaCL& TlDenseSymmetricMatrix_ViennaCL::operator*=(
    const double coef) {
  *(dynamic_cast<TlDenseSymmetricMatrix_ImplViennaCL*>(this->pImpl_)) *= coef;
}

TlDenseSymmetricMatrix_ViennaCL& TlDenseSymmetricMatrix_ViennaCL::operator/=(
    const double coef) {
  *(dynamic_cast<TlDenseSymmetricMatrix_ImplViennaCL*>(this->pImpl_)) /= coef;
}

TlDenseSymmetricMatrix_ViennaCL& TlDenseSymmetricMatrix_ViennaCL::operator*=(
    const TlDenseSymmetricMatrix_ViennaCL& rhs) {
  *(dynamic_cast<TlDenseSymmetricMatrix_ImplViennaCL*>(this->pImpl_)) *=
      *(dynamic_cast<TlDenseSymmetricMatrix_ImplViennaCL*>(rhs.pImpl_));
}

// ---------------------------------------------------------------------------
// operations
// ---------------------------------------------------------------------------
const TlDenseSymmetricMatrix_ViennaCL&
TlDenseSymmetricMatrix_ViennaCL::dotInPlace(
    const TlDenseSymmetricMatrix_ViennaCL& rhs) {
  dynamic_cast<TlDenseSymmetricMatrix_ImplViennaCL*>(this->pImpl_)
      ->dotInPlace(
          *dynamic_cast<TlDenseSymmetricMatrix_ImplViennaCL*>(rhs.pImpl_));

  return *this;
}

bool TlDenseSymmetricMatrix_ViennaCL::eig(
    TlDenseVector_ViennaCL* pEigVal,
    TlDenseGeneralMatrix_ViennaCL* pEigVec) const {
  TlDenseVector_ImplViennaCL* pImpl_eigval =
      dynamic_cast<TlDenseVector_ImplViennaCL*>(pEigVal->pImpl_);
  TlDenseGeneralMatrix_ImplViennaCL* pImpl_eigvec =
      dynamic_cast<TlDenseGeneralMatrix_ImplViennaCL*>(pEigVec->pImpl_);

  const bool answer =
      dynamic_cast<TlDenseSymmetricMatrix_ImplViennaCL*>(this->pImpl_)
          ->eig(pImpl_eigval, pImpl_eigvec);
  return answer;
}

TlDenseSymmetricMatrix_ViennaCL TlDenseSymmetricMatrix_ViennaCL::inverse()
    const {
  TlDenseSymmetricMatrix_ViennaCL answer;
  answer.pImpl_ = new TlDenseSymmetricMatrix_ImplViennaCL(
      dynamic_cast<const TlDenseSymmetricMatrix_ImplViennaCL*>(this->pImpl_)
          ->inverse());
  return answer;
}

// ---------------------------------------------------------------------------
// friend functions
// ---------------------------------------------------------------------------
