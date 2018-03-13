#include "tl_dense_symmetric_matrix_lapack.h"
#include "tl_dense_general_matrix_lapack.h"
#include "tl_dense_general_matrix_impl_lapack.h"
#include "tl_dense_symmetric_matrix_impl_lapack.h"
#include "tl_dense_vector_impl_lapack.h"
#include "tl_dense_vector_lapack.h"

TlDenseSymmetricMatrix_Lapack::TlDenseSymmetricMatrix_Lapack(
    const TlMatrixObject::index_type dim) {
  this->pImpl_ = new TlDenseSymmetricMatrix_ImplLapack(dim);
}

TlDenseSymmetricMatrix_Lapack::TlDenseSymmetricMatrix_Lapack(
    const TlDenseSymmetricMatrix_Lapack& rhs) {
  this->pImpl_ = new TlDenseSymmetricMatrix_ImplLapack(
      *(dynamic_cast<const TlDenseSymmetricMatrix_ImplLapack*>(rhs.pImpl_)));
}

TlDenseSymmetricMatrix_Lapack::TlDenseSymmetricMatrix_Lapack(
    const TlDenseGeneralMatrix_Lapack& rhs) {
  this->pImpl_ = new TlDenseSymmetricMatrix_ImplLapack(
      *(dynamic_cast<const TlDenseGeneralMatrix_ImplLapack*>(rhs.pImpl_)));
}

TlDenseSymmetricMatrix_Lapack::~TlDenseSymmetricMatrix_Lapack() {
  delete this->pImpl_;
  this->pImpl_ = NULL;
}

// ---------------------------------------------------------------------------
// operators
// ---------------------------------------------------------------------------
// TlDenseSymmetricMatrix_Lapack::operator TlDenseGeneralMatrix_Lapack()
//     const {
//   TlDenseGeneralMatrix_Lapack answer(
//       *(dynamic_cast<TlDenseSymmetricMatrix_ImplLapack*>(this->pImpl_)));
//
//   return answer;
// }

TlDenseSymmetricMatrix_Lapack& TlDenseSymmetricMatrix_Lapack::operator=(
    const TlDenseSymmetricMatrix_Lapack& rhs) {
  delete this->pImpl_;
  this->pImpl_ = new TlDenseSymmetricMatrix_ImplLapack(
      *(dynamic_cast<TlDenseSymmetricMatrix_ImplLapack*>(rhs.pImpl_)));
}

const TlDenseSymmetricMatrix_Lapack TlDenseSymmetricMatrix_Lapack::
operator+(const TlDenseSymmetricMatrix_Lapack& rhs) const {
  TlDenseSymmetricMatrix_Lapack answer = *this;
  answer += rhs;
  return answer;
}

const TlDenseSymmetricMatrix_Lapack TlDenseSymmetricMatrix_Lapack::
operator-(const TlDenseSymmetricMatrix_Lapack& rhs) const {
  TlDenseSymmetricMatrix_Lapack answer = *this;
  answer -= rhs;
  return answer;
}

const TlDenseGeneralMatrix_Lapack TlDenseSymmetricMatrix_Lapack::operator*(
    const TlDenseSymmetricMatrix_Lapack& rhs) const {
  TlDenseGeneralMatrix_Lapack answer = *this;
  answer *= rhs;
  return answer;
}

TlDenseSymmetricMatrix_Lapack& TlDenseSymmetricMatrix_Lapack::operator+=(
    const TlDenseSymmetricMatrix_Lapack& rhs) {
  *(dynamic_cast<TlDenseSymmetricMatrix_ImplLapack*>(this->pImpl_)) +=
      *(dynamic_cast<TlDenseSymmetricMatrix_ImplLapack*>(rhs.pImpl_));
}

TlDenseSymmetricMatrix_Lapack& TlDenseSymmetricMatrix_Lapack::operator-=(
    const TlDenseSymmetricMatrix_Lapack& rhs) {
  *(dynamic_cast<TlDenseSymmetricMatrix_ImplLapack*>(this->pImpl_)) -=
      *(dynamic_cast<TlDenseSymmetricMatrix_ImplLapack*>(rhs.pImpl_));
}

TlDenseSymmetricMatrix_Lapack& TlDenseSymmetricMatrix_Lapack::operator*=(
    const double coef) {
  *(dynamic_cast<TlDenseSymmetricMatrix_ImplLapack*>(this->pImpl_)) *= coef;
}

TlDenseSymmetricMatrix_Lapack& TlDenseSymmetricMatrix_Lapack::operator/=(
    const double coef) {
  *(dynamic_cast<TlDenseSymmetricMatrix_ImplLapack*>(this->pImpl_)) /= coef;
}

TlDenseSymmetricMatrix_Lapack& TlDenseSymmetricMatrix_Lapack::operator*=(
    const TlDenseSymmetricMatrix_Lapack& rhs) {
  *(dynamic_cast<TlDenseSymmetricMatrix_ImplLapack*>(this->pImpl_)) *=
      *(dynamic_cast<TlDenseSymmetricMatrix_ImplLapack*>(rhs.pImpl_));
}

// ---------------------------------------------------------------------------
// operations
// ---------------------------------------------------------------------------
const TlDenseSymmetricMatrix_Lapack&
TlDenseSymmetricMatrix_Lapack::dotInPlace(
    const TlDenseSymmetricMatrix_Lapack& rhs) {
  dynamic_cast<TlDenseSymmetricMatrix_ImplLapack*>(this->pImpl_)
      ->dotInPlace(*dynamic_cast<TlDenseSymmetricMatrix_ImplLapack*>(rhs.pImpl_));

  return *this;
}

bool TlDenseSymmetricMatrix_Lapack::eig(
    TlDenseVector_Lapack* pEigVal,
    TlDenseGeneralMatrix_Lapack* pEigVec) const {
  TlDenseVector_ImplLapack* pImpl_eigval =
      dynamic_cast<TlDenseVector_ImplLapack*>(pEigVal->pImpl_);
  TlDenseGeneralMatrix_ImplLapack* pImpl_eigvec =
      dynamic_cast<TlDenseGeneralMatrix_ImplLapack*>(pEigVec->pImpl_);

  const bool answer =
      dynamic_cast<TlDenseSymmetricMatrix_ImplLapack*>(this->pImpl_)
          ->eig(pImpl_eigval, pImpl_eigvec);
  return answer;
}

TlDenseSymmetricMatrix_Lapack TlDenseSymmetricMatrix_Lapack::inverse()
    const {
  TlDenseSymmetricMatrix_Lapack answer;
  answer.pImpl_ = new TlDenseSymmetricMatrix_ImplLapack(
      dynamic_cast<const TlDenseSymmetricMatrix_ImplLapack*>(this->pImpl_)
          ->inverse());
  return answer;
}

// ---------------------------------------------------------------------------
// friend functions
// ---------------------------------------------------------------------------
