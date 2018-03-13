#include "tl_dense_symmetric_matrix_eigen.h"
#include "tl_dense_general_matrix_eigen.h"
#include "tl_dense_symmetric_matrix_impl_eigen.h"
#include "tl_dense_vector_eigen.h"
#include "tl_dense_vector_impl_eigen.h"

TlDenseSymmetricMatrix_Eigen::TlDenseSymmetricMatrix_Eigen(
    const TlMatrixObject::index_type dim) {
  this->pImpl_ = new TlDenseSymmetricMatrix_ImplEigen(dim);
}

TlDenseSymmetricMatrix_Eigen::TlDenseSymmetricMatrix_Eigen(
    const TlDenseSymmetricMatrix_Eigen& rhs) {
  this->pImpl_ = new TlDenseSymmetricMatrix_ImplEigen(
      *(dynamic_cast<const TlDenseSymmetricMatrix_ImplEigen*>(rhs.pImpl_)));
}

TlDenseSymmetricMatrix_Eigen::TlDenseSymmetricMatrix_Eigen(
    const TlDenseGeneralMatrix_Eigen& rhs) {
  this->pImpl_ = new TlDenseSymmetricMatrix_ImplEigen(
      *(dynamic_cast<const TlDenseGeneralMatrix_ImplEigen*>(rhs.pImpl_)));
}

TlDenseSymmetricMatrix_Eigen::~TlDenseSymmetricMatrix_Eigen() {
  delete this->pImpl_;
  this->pImpl_ = NULL;
}

// ---------------------------------------------------------------------------
// operators
// ---------------------------------------------------------------------------
TlDenseSymmetricMatrix_Eigen& TlDenseSymmetricMatrix_Eigen::operator=(
    const TlDenseSymmetricMatrix_Eigen& rhs) {
  delete this->pImpl_;
  this->pImpl_ = new TlDenseSymmetricMatrix_ImplEigen(
      *(dynamic_cast<TlDenseSymmetricMatrix_ImplEigen*>(rhs.pImpl_)));
}

const TlDenseSymmetricMatrix_Eigen TlDenseSymmetricMatrix_Eigen::operator+(
    const TlDenseSymmetricMatrix_Eigen& rhs) const {
  TlDenseSymmetricMatrix_Eigen answer = *this;
  answer += rhs;
  return answer;
}

const TlDenseSymmetricMatrix_Eigen TlDenseSymmetricMatrix_Eigen::operator-(
    const TlDenseSymmetricMatrix_Eigen& rhs) const {
  TlDenseSymmetricMatrix_Eigen answer = *this;
  answer -= rhs;
  return answer;
}

const TlDenseSymmetricMatrix_Eigen TlDenseSymmetricMatrix_Eigen::operator*(
    const TlDenseSymmetricMatrix_Eigen& rhs) const {
  TlDenseSymmetricMatrix_Eigen answer = *this;
  answer *= rhs;
  return answer;
}

TlDenseSymmetricMatrix_Eigen& TlDenseSymmetricMatrix_Eigen::operator+=(
    const TlDenseSymmetricMatrix_Eigen& rhs) {
  *(dynamic_cast<TlDenseSymmetricMatrix_ImplEigen*>(this->pImpl_)) +=
      *(dynamic_cast<TlDenseSymmetricMatrix_ImplEigen*>(rhs.pImpl_));
}

TlDenseSymmetricMatrix_Eigen& TlDenseSymmetricMatrix_Eigen::operator-=(
    const TlDenseSymmetricMatrix_Eigen& rhs) {
  *(dynamic_cast<TlDenseSymmetricMatrix_ImplEigen*>(this->pImpl_)) -=
      *(dynamic_cast<TlDenseSymmetricMatrix_ImplEigen*>(rhs.pImpl_));
}

TlDenseSymmetricMatrix_Eigen& TlDenseSymmetricMatrix_Eigen::operator*=(
    const double coef) {
  *(dynamic_cast<TlDenseSymmetricMatrix_ImplEigen*>(this->pImpl_)) *= coef;
}

TlDenseSymmetricMatrix_Eigen& TlDenseSymmetricMatrix_Eigen::operator/=(
    const double coef) {
  *(dynamic_cast<TlDenseSymmetricMatrix_ImplEigen*>(this->pImpl_)) /= coef;
}

TlDenseSymmetricMatrix_Eigen& TlDenseSymmetricMatrix_Eigen::operator*=(
    const TlDenseSymmetricMatrix_Eigen& rhs) {
  *(dynamic_cast<TlDenseSymmetricMatrix_ImplEigen*>(this->pImpl_)) *=
      *(dynamic_cast<TlDenseSymmetricMatrix_ImplEigen*>(rhs.pImpl_));
}

// ---------------------------------------------------------------------------
// operations
// ---------------------------------------------------------------------------
const TlDenseSymmetricMatrix_Eigen& TlDenseSymmetricMatrix_Eigen::dotInPlace(
    const TlDenseSymmetricMatrix_Eigen& rhs) {
  dynamic_cast<TlDenseSymmetricMatrix_ImplEigen*>(this->pImpl_)
      ->dotInPlace(
          *dynamic_cast<TlDenseSymmetricMatrix_ImplEigen*>(rhs.pImpl_));

  return *this;
}

bool TlDenseSymmetricMatrix_Eigen::eig(
    TlDenseVector_Eigen* pEigVal, TlDenseGeneralMatrix_Eigen* pEigVec) const {
  TlDenseVector_ImplEigen* pImpl_eigval =
      dynamic_cast<TlDenseVector_ImplEigen*>(pEigVal->pImpl_);
  TlDenseGeneralMatrix_ImplEigen* pImpl_eigvec =
      dynamic_cast<TlDenseGeneralMatrix_ImplEigen*>(pEigVec->pImpl_);

  const bool answer =
      dynamic_cast<TlDenseSymmetricMatrix_ImplEigen*>(this->pImpl_)
          ->eig(pImpl_eigval, pImpl_eigvec);
  return answer;
}

TlDenseSymmetricMatrix_Eigen TlDenseSymmetricMatrix_Eigen::inverse() const {
  TlDenseSymmetricMatrix_Eigen answer;
  answer.pImpl_ = new TlDenseSymmetricMatrix_ImplEigen(
      dynamic_cast<const TlDenseSymmetricMatrix_ImplEigen*>(this->pImpl_)
          ->inverse());
  return answer;
}

// ---------------------------------------------------------------------------
// friend functions
// ---------------------------------------------------------------------------
