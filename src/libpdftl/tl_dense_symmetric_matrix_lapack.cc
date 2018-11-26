#include "tl_dense_symmetric_matrix_lapack.h"
#include "tl_dense_general_matrix_impl_lapack.h"
#include "tl_dense_general_matrix_lapack.h"
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

void TlDenseSymmetricMatrix_Lapack::vtr2mat(const std::vector<double>& vtr){
    dynamic_cast<TlDenseSymmetricMatrix_ImplLapack*>(this->pImpl_)->vtr2mat(vtr);
}

TlDenseSymmetricMatrix_Lapack::~TlDenseSymmetricMatrix_Lapack() {
  delete this->pImpl_;
  this->pImpl_ = NULL;
}

// ---------------------------------------------------------------------------
// properties
// ---------------------------------------------------------------------------

// ---------------------------------------------------------------------------
// operators
// ---------------------------------------------------------------------------
TlDenseSymmetricMatrix_Lapack& TlDenseSymmetricMatrix_Lapack::operator=(
    const TlDenseSymmetricMatrix_Lapack& rhs) {
  delete this->pImpl_;
  this->pImpl_ = new TlDenseSymmetricMatrix_ImplLapack(
      *(dynamic_cast<TlDenseSymmetricMatrix_ImplLapack*>(rhs.pImpl_)));

  return *this;
}

const TlDenseSymmetricMatrix_Lapack TlDenseSymmetricMatrix_Lapack::operator+(
    const TlDenseSymmetricMatrix_Lapack& rhs) const {
  TlDenseSymmetricMatrix_Lapack answer = *this;
  answer += rhs;
  return answer;
}

const TlDenseSymmetricMatrix_Lapack TlDenseSymmetricMatrix_Lapack::operator-(
    const TlDenseSymmetricMatrix_Lapack& rhs) const {
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

  return *this;
}

TlDenseSymmetricMatrix_Lapack& TlDenseSymmetricMatrix_Lapack::operator-=(
    const TlDenseSymmetricMatrix_Lapack& rhs) {
  *(dynamic_cast<TlDenseSymmetricMatrix_ImplLapack*>(this->pImpl_)) -=
      *(dynamic_cast<TlDenseSymmetricMatrix_ImplLapack*>(rhs.pImpl_));

  return *this;
}

TlDenseSymmetricMatrix_Lapack& TlDenseSymmetricMatrix_Lapack::operator*=(
    const double coef) {
  *(dynamic_cast<TlDenseSymmetricMatrix_ImplLapack*>(this->pImpl_)) *= coef;

  return *this;
}

TlDenseSymmetricMatrix_Lapack& TlDenseSymmetricMatrix_Lapack::operator/=(
    const double coef) {
  *(dynamic_cast<TlDenseSymmetricMatrix_ImplLapack*>(this->pImpl_)) /= coef;

  return *this;
}

// TlDenseSymmetricMatrix_Lapack& TlDenseSymmetricMatrix_Lapack::operator*=(
//     const TlDenseSymmetricMatrix_Lapack& rhs) {
//   *(dynamic_cast<TlDenseSymmetricMatrix_ImplLapack*>(this->pImpl_)) *=
//       *(dynamic_cast<TlDenseSymmetricMatrix_ImplLapack*>(rhs.pImpl_));
// }

// ---------------------------------------------------------------------------
// operations
// ---------------------------------------------------------------------------
const TlDenseSymmetricMatrix_Lapack& TlDenseSymmetricMatrix_Lapack::dotInPlace(
    const TlDenseSymmetricMatrix_Lapack& rhs) {
  dynamic_cast<TlDenseSymmetricMatrix_ImplLapack*>(this->pImpl_)
      ->dotInPlace(
          *dynamic_cast<TlDenseSymmetricMatrix_ImplLapack*>(rhs.pImpl_));

  return *this;
}

bool TlDenseSymmetricMatrix_Lapack::eig(
    TlDenseVector_Lapack* pEigVal, TlDenseGeneralMatrix_Lapack* pEigVec) const {
  TlDenseVector_ImplLapack* pImpl_eigval =
      dynamic_cast<TlDenseVector_ImplLapack*>(pEigVal->pImpl_);
  TlDenseGeneralMatrix_ImplLapack* pImpl_eigvec =
      dynamic_cast<TlDenseGeneralMatrix_ImplLapack*>(pEigVec->pImpl_);

  const bool answer =
      dynamic_cast<TlDenseSymmetricMatrix_ImplLapack*>(this->pImpl_)
          ->eig(pImpl_eigval, pImpl_eigvec);
  return answer;
}

TlDenseSymmetricMatrix_Lapack TlDenseSymmetricMatrix_Lapack::inverse() const {
  TlDenseSymmetricMatrix_Lapack answer;
  answer.pImpl_ = new TlDenseSymmetricMatrix_ImplLapack(
      dynamic_cast<const TlDenseSymmetricMatrix_ImplLapack*>(this->pImpl_)
          ->inverse());
  return answer;
}

double* TlDenseSymmetricMatrix_Lapack::data() {
  return dynamic_cast<TlDenseSymmetricMatrix_ImplLapack*>(this->pImpl_)->data();
}

const double* TlDenseSymmetricMatrix_Lapack::data() const {
  return dynamic_cast<TlDenseSymmetricMatrix_ImplLapack*>(this->pImpl_)->data();
}

// ---------------------------------------------------------------------------
// I/O
// ---------------------------------------------------------------------------
void TlDenseSymmetricMatrix_Lapack::dump(TlDenseVector_Lapack* v) const {
  dynamic_cast<TlDenseSymmetricMatrix_ImplLapack*>(this->pImpl_)
      ->dump(dynamic_cast<TlDenseVector_ImplLapack*>(v->pImpl_));
}

void TlDenseSymmetricMatrix_Lapack::restore(const TlDenseVector_Lapack& v) {
  dynamic_cast<TlDenseSymmetricMatrix_ImplLapack*>(this->pImpl_)
      ->restore(*(dynamic_cast<TlDenseVector_ImplLapack*>(v.pImpl_)));
}

// ---------------------------------------------------------------------------
// protected
// ---------------------------------------------------------------------------
TlMatrixObject::size_type TlDenseSymmetricMatrix_Lapack::getNumOfElements()
    const {
  return dynamic_cast<TlDenseSymmetricMatrix_ImplLapack*>(this->pImpl_)
      ->getNumOfElements();
}

// ---------------------------------------------------------------------------
// friend functions
// ---------------------------------------------------------------------------
TlDenseGeneralMatrix_Lapack operator*(const TlDenseSymmetricMatrix_Lapack& rhs1,
                                      const TlDenseGeneralMatrix_Lapack& rhs2) {
  TlDenseGeneralMatrix_Lapack answer;
  *(dynamic_cast<TlDenseGeneralMatrix_ImplLapack*>(answer.pImpl_)) =
      *(dynamic_cast<TlDenseSymmetricMatrix_ImplLapack*>(rhs1.pImpl_)) *
      *(dynamic_cast<TlDenseGeneralMatrix_ImplLapack*>(rhs2.pImpl_));

  return answer;
}

TlDenseGeneralMatrix_Lapack operator*(
    const TlDenseGeneralMatrix_Lapack& rhs1,
    const TlDenseSymmetricMatrix_Lapack& rhs2) {
  TlDenseGeneralMatrix_Lapack answer;
  *(dynamic_cast<TlDenseGeneralMatrix_ImplLapack*>(answer.pImpl_)) =
      *(dynamic_cast<TlDenseGeneralMatrix_ImplLapack*>(rhs1.pImpl_)) *
      *(dynamic_cast<TlDenseSymmetricMatrix_ImplLapack*>(rhs2.pImpl_));

  return answer;
}

TlDenseVector_Lapack operator*(const TlDenseSymmetricMatrix_Lapack& rhs1,
                               const TlDenseVector_Lapack& rhs2) {
  TlDenseVector_Lapack answer;
  *(dynamic_cast<TlDenseVector_ImplLapack*>(answer.pImpl_)) =
      *(dynamic_cast<TlDenseSymmetricMatrix_ImplLapack*>(rhs1.pImpl_)) *
      *(dynamic_cast<TlDenseVector_ImplLapack*>(rhs2.pImpl_));

  return answer;
}

TlDenseVector_Lapack operator*(const TlDenseVector_Lapack& rhs1,
                               const TlDenseSymmetricMatrix_Lapack& rhs2) {
  TlDenseVector_Lapack answer;
  *(dynamic_cast<TlDenseVector_ImplLapack*>(answer.pImpl_)) =
      *(dynamic_cast<TlDenseSymmetricMatrix_ImplLapack*>(rhs2.pImpl_)) *
      *(dynamic_cast<TlDenseVector_ImplLapack*>(rhs1.pImpl_));

  return answer;
}

// ---------------------------------------------------------------------------
// arithmetic
// ---------------------------------------------------------------------------
TlDenseSymmetricMatrix_Lapack operator*(
    const double coef, const TlDenseSymmetricMatrix_Lapack& matrix) {
  TlDenseSymmetricMatrix_Lapack answer = matrix;
  answer *= coef;

  return answer;
}

TlDenseSymmetricMatrix_Lapack operator*(
    const TlDenseSymmetricMatrix_Lapack& matrix, const double coef) {
  return coef * matrix;
}
