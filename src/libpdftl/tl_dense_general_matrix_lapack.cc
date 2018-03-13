#include "tl_dense_general_matrix_lapack.h"
#include "tl_dense_symmetric_matrix_lapack.h"
#include "tl_dense_symmetric_matrix_impl_lapack.h"
#include "tl_dense_vector_lapack.h"
#include "tl_dense_vector_impl_lapack.h"

TlDenseGeneralMatrix_Lapack::TlDenseGeneralMatrix_Lapack(
    const TlMatrixObject::index_type row,
    const TlMatrixObject::index_type col) {
  this->pImpl_ = new TlDenseGeneralMatrix_ImplLapack(row, col);
}

TlDenseGeneralMatrix_Lapack::TlDenseGeneralMatrix_Lapack(
    const TlDenseGeneralMatrix_Lapack& rhs) {
  this->pImpl_ = new TlDenseGeneralMatrix_ImplLapack(
      *(dynamic_cast<TlDenseGeneralMatrix_ImplLapack*>(rhs.pImpl_)));
}

TlDenseGeneralMatrix_Lapack::TlDenseGeneralMatrix_Lapack(
    const TlDenseSymmetricMatrix_Lapack& rhs) {
  this->pImpl_ = new TlDenseGeneralMatrix_ImplLapack(
      *(dynamic_cast<const TlDenseSymmetricMatrix_ImplLapack*>(rhs.pImpl_)));
}

TlDenseGeneralMatrix_Lapack::TlDenseGeneralMatrix_Lapack(
    const TlDenseGeneralMatrix_ImplLapack& rhs) {
  this->pImpl_ = new TlDenseGeneralMatrix_ImplLapack(rhs);
}

TlDenseGeneralMatrix_Lapack::~TlDenseGeneralMatrix_Lapack() {
  delete this->pImpl_;
  this->pImpl_ = NULL;
}

// ---------------------------------------------------------------------------
// operators
// ---------------------------------------------------------------------------
TlDenseGeneralMatrix_Lapack& TlDenseGeneralMatrix_Lapack::operator=(
    const TlDenseGeneralMatrix_Lapack& rhs) {
  delete this->pImpl_;
  this->pImpl_ = new TlDenseGeneralMatrix_ImplLapack(
      *(dynamic_cast<TlDenseGeneralMatrix_ImplLapack*>(rhs.pImpl_)));
}

const TlDenseGeneralMatrix_Lapack TlDenseGeneralMatrix_Lapack::operator+(
    const TlDenseGeneralMatrix_Lapack& rhs) const {
  TlDenseGeneralMatrix_Lapack answer = *this;
  answer += rhs;
  return answer;
}

const TlDenseGeneralMatrix_Lapack TlDenseGeneralMatrix_Lapack::operator-(
    const TlDenseGeneralMatrix_Lapack& rhs) const {
  TlDenseGeneralMatrix_Lapack answer = *this;
  answer -= rhs;
  return answer;
}

const TlDenseGeneralMatrix_Lapack TlDenseGeneralMatrix_Lapack::operator*(
    const TlDenseGeneralMatrix_Lapack& rhs) const {
  TlDenseGeneralMatrix_Lapack answer = *this;
  answer *= rhs;
  return answer;
}

TlDenseGeneralMatrix_Lapack& TlDenseGeneralMatrix_Lapack::operator+=(
    const TlDenseGeneralMatrix_Lapack& rhs) {
  *(dynamic_cast<TlDenseGeneralMatrix_ImplLapack*>(this->pImpl_)) +=
      *(dynamic_cast<TlDenseGeneralMatrix_ImplLapack*>(rhs.pImpl_));
}

TlDenseGeneralMatrix_Lapack& TlDenseGeneralMatrix_Lapack::operator-=(
    const TlDenseGeneralMatrix_Lapack& rhs) {
  *(dynamic_cast<TlDenseGeneralMatrix_ImplLapack*>(this->pImpl_)) -=
      *(dynamic_cast<TlDenseGeneralMatrix_ImplLapack*>(rhs.pImpl_));
}

TlDenseGeneralMatrix_Lapack& TlDenseGeneralMatrix_Lapack::operator*=(
    const double coef) {
  *(dynamic_cast<TlDenseGeneralMatrix_ImplLapack*>(this->pImpl_)) *= coef;
}

TlDenseGeneralMatrix_Lapack& TlDenseGeneralMatrix_Lapack::operator/=(
    const double coef) {
  *(dynamic_cast<TlDenseGeneralMatrix_ImplLapack*>(this->pImpl_)) /= coef;
}

TlDenseGeneralMatrix_Lapack& TlDenseGeneralMatrix_Lapack::operator*=(
    const TlDenseGeneralMatrix_Lapack& rhs) {
  *(dynamic_cast<TlDenseGeneralMatrix_ImplLapack*>(this->pImpl_)) *=
      *(dynamic_cast<TlDenseGeneralMatrix_ImplLapack*>(rhs.pImpl_));
}

// ---------------------------------------------------------------------------
// operations
// ---------------------------------------------------------------------------
double TlDenseGeneralMatrix_Lapack::sum() const {
  return this->pImpl_->sum();
}

double TlDenseGeneralMatrix_Lapack::getRMS() const {
  return this->pImpl_->getRMS();
}

double TlDenseGeneralMatrix_Lapack::getMaxAbsoluteElement(
    TlMatrixObject::index_type* outRow,
    TlMatrixObject::index_type* outCol) const {
  return this->pImpl_->getMaxAbsoluteElement(outRow, outCol);
}

const TlDenseGeneralMatrix_Lapack& TlDenseGeneralMatrix_Lapack::dotInPlace(
    const TlDenseGeneralMatrix_Lapack& rhs) {
  dynamic_cast<TlDenseGeneralMatrix_ImplLapack*>(this->pImpl_)
      ->dotInPlace(*(dynamic_cast<TlDenseGeneralMatrix_ImplLapack*>(rhs.pImpl_)));

  return *this;
}

TlDenseGeneralMatrix_Lapack TlDenseGeneralMatrix_Lapack::transpose() const {
  return TlDenseGeneralMatrix_Lapack(
      dynamic_cast<const TlDenseGeneralMatrix_ImplLapack*>(this->pImpl_)
          ->transpose());
}

TlDenseGeneralMatrix_Lapack TlDenseGeneralMatrix_Lapack::inverse() const {
  return TlDenseGeneralMatrix_Lapack(
      dynamic_cast<const TlDenseGeneralMatrix_ImplLapack*>(this->pImpl_)
          ->inverse());
}

// ---------------------------------------------------------------------------
// friend functions
// ---------------------------------------------------------------------------
TlDenseVector_Lapack operator*(const TlDenseGeneralMatrix_Lapack& rhs1,
                               const TlDenseVector_Lapack& rhs2) {
  TlDenseVector_Lapack answer;
  *(dynamic_cast<TlDenseVector_ImplLapack*>(answer.pImpl_)) =
      *(dynamic_cast<TlDenseGeneralMatrix_ImplLapack*>(rhs1.pImpl_)) *
      *(dynamic_cast<TlDenseVector_ImplLapack*>(rhs2.pImpl_));

  return answer;
}

TlDenseVector_Lapack operator*(const TlDenseVector_Lapack& rhs1,
                               const TlDenseGeneralMatrix_Lapack& rhs2) {
  TlDenseVector_Lapack answer;
  *(dynamic_cast<TlDenseVector_ImplLapack*>(answer.pImpl_)) =
      *(dynamic_cast<TlDenseVector_ImplLapack*>(rhs1.pImpl_)) *
      *(dynamic_cast<TlDenseGeneralMatrix_ImplLapack*>(rhs2.pImpl_));

  return answer;
}
