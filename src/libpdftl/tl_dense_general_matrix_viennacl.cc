#include "tl_dense_general_matrix_viennacl.h"
#include "tl_dense_general_matrix_eigen.h"
#include "tl_dense_general_matrix_impl_eigen.h"
#include "tl_dense_general_matrix_impl_viennacl.h"
#include "tl_dense_symmetric_matrix_impl_viennacl.h"
#include "tl_dense_symmetric_matrix_viennacl.h"
#include "tl_dense_vector_viennacl.h"
#include "tl_dense_vector_impl_viennacl.h"

TlDenseGeneralMatrix_ViennaCL::TlDenseGeneralMatrix_ViennaCL(
    const TlMatrixObject::index_type row,
    const TlMatrixObject::index_type col) {
  this->pImpl_ = new TlDenseGeneralMatrix_ImplViennaCL(row, col);
}

TlDenseGeneralMatrix_ViennaCL::TlDenseGeneralMatrix_ViennaCL(
    const TlDenseGeneralMatrix_ViennaCL& rhs) {
  this->pImpl_ = new TlDenseGeneralMatrix_ImplViennaCL(
      *(dynamic_cast<const TlDenseGeneralMatrix_ImplViennaCL*>(rhs.pImpl_)));
}

TlDenseGeneralMatrix_ViennaCL::TlDenseGeneralMatrix_ViennaCL(
    const TlDenseSymmetricMatrix_ViennaCL& rhs) {
  this->pImpl_ = new TlDenseGeneralMatrix_ImplViennaCL(
      *(dynamic_cast<const TlDenseSymmetricMatrix_ImplViennaCL*>(rhs.pImpl_)));
}

TlDenseGeneralMatrix_ViennaCL::TlDenseGeneralMatrix_ViennaCL(
    const TlDenseGeneralMatrix_ImplViennaCL& rhs) {
  this->pImpl_ = new TlDenseGeneralMatrix_ImplViennaCL(rhs);
}

TlDenseGeneralMatrix_ViennaCL::TlDenseGeneralMatrix_ViennaCL(
    const TlDenseGeneralMatrix_Eigen& rhs) {
  this->pImpl_ = new TlDenseGeneralMatrix_ImplViennaCL(
      *dynamic_cast<const TlDenseGeneralMatrix_ImplEigen*>(rhs.pImpl_));
}

TlDenseGeneralMatrix_ViennaCL::~TlDenseGeneralMatrix_ViennaCL() {
  delete this->pImpl_;
  this->pImpl_ = NULL;
}

// ---------------------------------------------------------------------------
// operators
// ---------------------------------------------------------------------------
TlDenseGeneralMatrix_ViennaCL& TlDenseGeneralMatrix_ViennaCL::operator=(
    const TlDenseGeneralMatrix_ViennaCL& rhs) {
  delete this->pImpl_;
  this->pImpl_ = new TlDenseGeneralMatrix_ImplViennaCL(
      *(dynamic_cast<TlDenseGeneralMatrix_ImplViennaCL*>(rhs.pImpl_)));

  return *this;
}

// TlDenseGeneralMatrix_ViennaCL& TlDenseGeneralMatrix_ViennaCL::operator=(const
// TlDenseGeneralMatrix_Eigen_Old& rhs) {
//     *(this->pImpl_) = *(dynamic_cast<TlDenselMatrix_ImplEigen>(rhs.pImpl_));
// }

const TlDenseGeneralMatrix_ViennaCL TlDenseGeneralMatrix_ViennaCL::operator+(
    const TlDenseGeneralMatrix_ViennaCL& rhs) const {
  TlDenseGeneralMatrix_ViennaCL answer = *this;
  answer += rhs;
  return answer;
}

const TlDenseGeneralMatrix_ViennaCL TlDenseGeneralMatrix_ViennaCL::operator-(
    const TlDenseGeneralMatrix_ViennaCL& rhs) const {
  TlDenseGeneralMatrix_ViennaCL answer = *this;
  answer -= rhs;
  return answer;
}

const TlDenseGeneralMatrix_ViennaCL TlDenseGeneralMatrix_ViennaCL::operator*(
    const TlDenseGeneralMatrix_ViennaCL& rhs) const {
  TlDenseGeneralMatrix_ViennaCL answer = *this;
  answer *= rhs;
  return answer;
}

TlDenseGeneralMatrix_ViennaCL& TlDenseGeneralMatrix_ViennaCL::operator+=(
    const TlDenseGeneralMatrix_ViennaCL& rhs) {
  *(dynamic_cast<TlDenseGeneralMatrix_ImplViennaCL*>(this->pImpl_)) +=
      *(dynamic_cast<TlDenseGeneralMatrix_ImplViennaCL*>(rhs.pImpl_));

  return *this;
}

TlDenseGeneralMatrix_ViennaCL& TlDenseGeneralMatrix_ViennaCL::operator-=(
    const TlDenseGeneralMatrix_ViennaCL& rhs) {
  *(dynamic_cast<TlDenseGeneralMatrix_ImplViennaCL*>(this->pImpl_)) -=
      *(dynamic_cast<TlDenseGeneralMatrix_ImplViennaCL*>(rhs.pImpl_));

  return *this;
}

TlDenseGeneralMatrix_ViennaCL& TlDenseGeneralMatrix_ViennaCL::operator*=(
    const double coef) {
  *(dynamic_cast<TlDenseGeneralMatrix_ImplViennaCL*>(this->pImpl_)) *= coef;

  return *this;
}

TlDenseGeneralMatrix_ViennaCL& TlDenseGeneralMatrix_ViennaCL::operator/=(
    const double coef) {
  *(dynamic_cast<TlDenseGeneralMatrix_ImplViennaCL*>(this->pImpl_)) /= coef;

  return *this;
}

TlDenseGeneralMatrix_ViennaCL& TlDenseGeneralMatrix_ViennaCL::operator*=(
    const TlDenseGeneralMatrix_ViennaCL& rhs) {
  *(dynamic_cast<TlDenseGeneralMatrix_ImplViennaCL*>(this->pImpl_)) *=
      *(dynamic_cast<TlDenseGeneralMatrix_ImplViennaCL*>(rhs.pImpl_));

  return *this;
}

// ---------------------------------------------------------------------------
// operations
// ---------------------------------------------------------------------------
double TlDenseGeneralMatrix_ViennaCL::sum() const {
  return this->pImpl_->sum();
}

double TlDenseGeneralMatrix_ViennaCL::getRMS() const {
  return this->pImpl_->getRMS();
}

const TlDenseGeneralMatrix_ViennaCL& TlDenseGeneralMatrix_ViennaCL::dotInPlace(
    const TlDenseGeneralMatrix_ViennaCL& rhs) {
  dynamic_cast<TlDenseGeneralMatrix_ImplViennaCL*>(this->pImpl_)
      ->dotInPlace(
          *(dynamic_cast<TlDenseGeneralMatrix_ImplViennaCL*>(rhs.pImpl_)));
  return *this;
}

TlDenseGeneralMatrix_ViennaCL TlDenseGeneralMatrix_ViennaCL::transpose() const {
  return TlDenseGeneralMatrix_ViennaCL(
      dynamic_cast<const TlDenseGeneralMatrix_ImplViennaCL*>(this->pImpl_)
          ->transpose());
}

TlDenseGeneralMatrix_ViennaCL TlDenseGeneralMatrix_ViennaCL::inverse() const {
  return TlDenseGeneralMatrix_ViennaCL(
      dynamic_cast<const TlDenseGeneralMatrix_ImplViennaCL*>(this->pImpl_)
          ->inverse());
}

// ---------------------------------------------------------------------------
// friend functions
// ---------------------------------------------------------------------------
TlDenseVector_ViennaCL operator*(const TlDenseGeneralMatrix_ViennaCL& rhs1,
                                 const TlDenseVector_ViennaCL& rhs2) {
  TlDenseVector_ViennaCL answer;
  *(dynamic_cast<TlDenseVector_ImplViennaCL*>(answer.pImpl_)) =
      *(dynamic_cast<TlDenseGeneralMatrix_ImplViennaCL*>(rhs1.pImpl_)) *
      *(dynamic_cast<TlDenseVector_ImplViennaCL*>(rhs2.pImpl_));

  return answer;
}

TlDenseVector_ViennaCL operator*(const TlDenseVector_ViennaCL& rhs1,
                                 const TlDenseGeneralMatrix_ViennaCL& rhs2) {
  TlDenseVector_ViennaCL answer;
  *(dynamic_cast<TlDenseVector_ImplViennaCL*>(answer.pImpl_)) =
      *(dynamic_cast<TlDenseVector_ImplViennaCL*>(rhs1.pImpl_)) *
      *(dynamic_cast<TlDenseGeneralMatrix_ImplViennaCL*>(rhs2.pImpl_));

  return answer;
}
