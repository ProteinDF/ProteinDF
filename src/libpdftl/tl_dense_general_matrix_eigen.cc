#include "tl_dense_general_matrix_eigen.h"
#include "tl_dense_general_matrix_impl_eigen.h"
#include "tl_dense_symmetric_matrix_eigen.h"
#include "tl_dense_symmetric_matrix_impl_eigen.h"
#include "tl_dense_vector_eigen.h"
#include "tl_dense_vector_impl_eigen.h"

TlDenseGeneralMatrix_Eigen::TlDenseGeneralMatrix_Eigen(
    const TlMatrixObject::index_type row,
    const TlMatrixObject::index_type col) {
  this->pImpl_ = new TlDenseGeneralMatrix_ImplEigen(row, col);
}

TlDenseGeneralMatrix_Eigen::TlDenseGeneralMatrix_Eigen(
    const TlDenseGeneralMatrix_Eigen& rhs) {
  this->pImpl_ = new TlDenseGeneralMatrix_ImplEigen(
      *(dynamic_cast<const TlDenseGeneralMatrix_ImplEigen*>(rhs.pImpl_)));
}

TlDenseGeneralMatrix_Eigen::TlDenseGeneralMatrix_Eigen(
    const TlDenseSymmetricMatrix_Eigen& rhs) {
  this->pImpl_ = new TlDenseGeneralMatrix_ImplEigen(
      *(dynamic_cast<const TlDenseSymmetricMatrix_ImplEigen*>(rhs.pImpl_)));
}

TlDenseGeneralMatrix_Eigen::TlDenseGeneralMatrix_Eigen(
    const TlDenseGeneralMatrix_ImplEigen& rhs) {
  this->pImpl_ = new TlDenseGeneralMatrix_ImplEigen(rhs);
}

TlDenseGeneralMatrix_Eigen::~TlDenseGeneralMatrix_Eigen() {
  delete this->pImpl_;
  this->pImpl_ = NULL;
}

// ---------------------------------------------------------------------------
// operators
// ---------------------------------------------------------------------------
TlDenseGeneralMatrix_Eigen& TlDenseGeneralMatrix_Eigen::operator=(
    const TlDenseGeneralMatrix_Eigen& rhs) {
  delete this->pImpl_;
  this->pImpl_ = new TlDenseGeneralMatrix_ImplEigen(
      *(dynamic_cast<TlDenseGeneralMatrix_ImplEigen*>(rhs.pImpl_)));
}

const TlDenseGeneralMatrix_Eigen TlDenseGeneralMatrix_Eigen::operator+(
    const TlDenseGeneralMatrix_Eigen& rhs) const {
  TlDenseGeneralMatrix_Eigen answer = *this;
  answer += rhs;
  return answer;
}

const TlDenseGeneralMatrix_Eigen TlDenseGeneralMatrix_Eigen::operator-(
    const TlDenseGeneralMatrix_Eigen& rhs) const {
  TlDenseGeneralMatrix_Eigen answer = *this;
  answer -= rhs;
  return answer;
}

const TlDenseGeneralMatrix_Eigen TlDenseGeneralMatrix_Eigen::operator*(
    const TlDenseGeneralMatrix_Eigen& rhs) const {
  TlDenseGeneralMatrix_Eigen answer = *this;
  answer *= rhs;
  return answer;
}

TlDenseGeneralMatrix_Eigen& TlDenseGeneralMatrix_Eigen::operator+=(
    const TlDenseGeneralMatrix_Eigen& rhs) {
  *(dynamic_cast<TlDenseGeneralMatrix_ImplEigen*>(this->pImpl_)) +=
      *(dynamic_cast<TlDenseGeneralMatrix_ImplEigen*>(rhs.pImpl_));
}

TlDenseGeneralMatrix_Eigen& TlDenseGeneralMatrix_Eigen::operator-=(
    const TlDenseGeneralMatrix_Eigen& rhs) {
  *(dynamic_cast<TlDenseGeneralMatrix_ImplEigen*>(this->pImpl_)) -=
      *(dynamic_cast<TlDenseGeneralMatrix_ImplEigen*>(rhs.pImpl_));
}

TlDenseGeneralMatrix_Eigen& TlDenseGeneralMatrix_Eigen::operator*=(
    const double coef) {
  *(dynamic_cast<TlDenseGeneralMatrix_ImplEigen*>(this->pImpl_)) *= coef;
}

TlDenseGeneralMatrix_Eigen& TlDenseGeneralMatrix_Eigen::operator/=(
    const double coef) {
  *(dynamic_cast<TlDenseGeneralMatrix_ImplEigen*>(this->pImpl_)) /= coef;
}

TlDenseGeneralMatrix_Eigen& TlDenseGeneralMatrix_Eigen::operator*=(
    const TlDenseGeneralMatrix_Eigen& rhs) {
  *(dynamic_cast<TlDenseGeneralMatrix_ImplEigen*>(this->pImpl_)) *=
      *(dynamic_cast<TlDenseGeneralMatrix_ImplEigen*>(rhs.pImpl_));
}

// ---------------------------------------------------------------------------
// operations
// ---------------------------------------------------------------------------
double TlDenseGeneralMatrix_Eigen::sum() const {
  return this->pImpl_->sum();
}

double TlDenseGeneralMatrix_Eigen::getRMS() const {
  return this->pImpl_->getRMS();
}

double TlDenseGeneralMatrix_Eigen::getMaxAbsoluteElement(
    TlMatrixObject::index_type* outRow,
    TlMatrixObject::index_type* outCol) const {
  return this->pImpl_->getMaxAbsoluteElement(outRow, outCol);
}

const TlDenseGeneralMatrix_Eigen&
TlDenseGeneralMatrix_Eigen::dotInPlace(
    const TlDenseGeneralMatrix_Eigen& rhs) {
  dynamic_cast<TlDenseGeneralMatrix_ImplEigen*>(this->pImpl_)
      ->dotInPlace(
          *(dynamic_cast<const TlDenseGeneralMatrix_ImplEigen*>(rhs.pImpl_)));
  return *this;
}

TlDenseGeneralMatrix_Eigen TlDenseGeneralMatrix_Eigen::transpose()
    const {
  return TlDenseGeneralMatrix_Eigen(
      dynamic_cast<TlDenseGeneralMatrix_ImplEigen*>(this->pImpl_)->transpose());
}

TlDenseGeneralMatrix_Eigen TlDenseGeneralMatrix_Eigen::inverse() const {
  return TlDenseGeneralMatrix_Eigen(
      dynamic_cast<TlDenseGeneralMatrix_ImplEigen*>(this->pImpl_)->inverse());
}

// ---------------------------------------------------------------------------
// friend functions
// ---------------------------------------------------------------------------
TlDenseVector_Eigen operator*(const TlDenseGeneralMatrix_Eigen& rhs1,
                              const TlDenseVector_Eigen& rhs2) {
  TlDenseVector_Eigen answer;
  *(dynamic_cast<TlDenseVector_ImplEigen*>(answer.pImpl_)) =
      *(dynamic_cast<TlDenseGeneralMatrix_ImplEigen*>(rhs1.pImpl_)) *
      *(dynamic_cast<TlDenseVector_ImplEigen*>(rhs2.pImpl_));

  return answer;
}

TlDenseVector_Eigen operator*(const TlDenseVector_Eigen& rhs1,
                              const TlDenseGeneralMatrix_Eigen& rhs2) {
  TlDenseVector_Eigen answer;
  *(dynamic_cast<TlDenseVector_ImplEigen*>(answer.pImpl_)) =
      *(dynamic_cast<TlDenseVector_ImplEigen*>(rhs1.pImpl_)) *
      *(dynamic_cast<TlDenseGeneralMatrix_ImplEigen*>(rhs2.pImpl_));

  return answer;
}
