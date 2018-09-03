#include "tl_dense_symmetric_matrix_scalapack.h"
#include "tl_dense_general_matrix_impl_scalapack.h"
#include "tl_dense_general_matrix_scalapack.h"
#include "tl_dense_symmetric_matrix_impl_scalapack.h"
#include "tl_dense_vector_impl_lapack.h"
#include "tl_dense_vector_impl_scalapack.h"
#include "tl_dense_vector_lapack.h"
#include "tl_dense_vector_scalapack.h"

// ---------------------------------------------------------------------------
// constructor & destructor
// ---------------------------------------------------------------------------
TlDenseSymmetricMatrix_Scalapack::TlDenseSymmetricMatrix_Scalapack(
    const TlMatrixObject::index_type dim) {
  this->pImpl_ = new TlDenseSymmetricMatrix_ImplScalapack(dim);
}

TlDenseSymmetricMatrix_Scalapack::TlDenseSymmetricMatrix_Scalapack(
    const TlDenseSymmetricMatrix_Scalapack& rhs) {
  this->pImpl_ = new TlDenseSymmetricMatrix_ImplScalapack(
      *(dynamic_cast<const TlDenseSymmetricMatrix_ImplScalapack*>(rhs.pImpl_)));
}

TlDenseSymmetricMatrix_Scalapack::TlDenseSymmetricMatrix_Scalapack(
    const TlDenseGeneralMatrix_Scalapack& rhs) {
  this->pImpl_ = new TlDenseSymmetricMatrix_ImplScalapack(
      *(dynamic_cast<const TlDenseGeneralMatrix_ImplScalapack*>(rhs.pImpl_)));
}

// TlDenseSymmetricMatrix_Scalapack::TlDenseSymmetricMatrix_Scalapack(
//     const TlDenseVector_Scalapack& v, const TlMatrixObject::index_type dim) {
//   this->pImpl_ = new TlDenseSymmetricMatrix_ImplScalapack(
//       *(dynamic_cast<const TlDenseVector_ImplScalapack*>(v.pImpl_)), dim);
// }

TlDenseSymmetricMatrix_Scalapack::~TlDenseSymmetricMatrix_Scalapack() {
  delete this->pImpl_;
  this->pImpl_ = NULL;
}

// ---------------------------------------------------------------------------
// properties
// ---------------------------------------------------------------------------
TlMatrixObject::size_type TlDenseSymmetricMatrix_Scalapack::getNumOfElements()
    const {
  return dynamic_cast<TlDenseSymmetricMatrix_ImplScalapack*>(this->pImpl_)
      ->getNumOfElements();
}

double TlDenseSymmetricMatrix_Scalapack::getLocal(
    TlMatrixObject::index_type row, TlMatrixObject::index_type col) const {
  return dynamic_cast<TlDenseSymmetricMatrix_ImplScalapack*>(this->pImpl_)
      ->getLocal(row, col);
}

// ---------------------------------------------------------------------------
// operators
// ---------------------------------------------------------------------------
TlDenseSymmetricMatrix_Scalapack& TlDenseSymmetricMatrix_Scalapack::operator=(
    const TlDenseSymmetricMatrix_Scalapack& rhs) {
  delete this->pImpl_;
  this->pImpl_ = new TlDenseSymmetricMatrix_ImplScalapack(
      *(dynamic_cast<TlDenseSymmetricMatrix_ImplScalapack*>(rhs.pImpl_)));

  return *this;
}

const TlDenseSymmetricMatrix_Scalapack TlDenseSymmetricMatrix_Scalapack::
operator+(const TlDenseSymmetricMatrix_Scalapack& rhs) const {
  TlDenseSymmetricMatrix_Scalapack answer = *this;
  answer += rhs;
  return answer;
}

const TlDenseSymmetricMatrix_Scalapack TlDenseSymmetricMatrix_Scalapack::
operator-(const TlDenseSymmetricMatrix_Scalapack& rhs) const {
  TlDenseSymmetricMatrix_Scalapack answer = *this;
  answer -= rhs;
  return answer;
}

const TlDenseGeneralMatrix_Scalapack TlDenseSymmetricMatrix_Scalapack::
operator*(const TlDenseSymmetricMatrix_Scalapack& rhs) const {
  TlDenseGeneralMatrix_Scalapack answer = *this;
  answer *= rhs;
  return answer;
}

TlDenseSymmetricMatrix_Scalapack& TlDenseSymmetricMatrix_Scalapack::operator+=(
    const TlDenseSymmetricMatrix_Scalapack& rhs) {
  *(dynamic_cast<TlDenseSymmetricMatrix_ImplScalapack*>(this->pImpl_)) +=
      *(dynamic_cast<TlDenseSymmetricMatrix_ImplScalapack*>(rhs.pImpl_));

  return *this;
}

TlDenseSymmetricMatrix_Scalapack& TlDenseSymmetricMatrix_Scalapack::operator-=(
    const TlDenseSymmetricMatrix_Scalapack& rhs) {
  *(dynamic_cast<TlDenseSymmetricMatrix_ImplScalapack*>(this->pImpl_)) -=
      *(dynamic_cast<TlDenseSymmetricMatrix_ImplScalapack*>(rhs.pImpl_));

  return *this;
}

TlDenseSymmetricMatrix_Scalapack& TlDenseSymmetricMatrix_Scalapack::operator*=(
    const double coef) {
  *(dynamic_cast<TlDenseSymmetricMatrix_ImplScalapack*>(this->pImpl_)) *= coef;

  return *this;
}

TlDenseSymmetricMatrix_Scalapack& TlDenseSymmetricMatrix_Scalapack::operator/=(
    const double coef) {
  *(dynamic_cast<TlDenseSymmetricMatrix_ImplScalapack*>(this->pImpl_)) /= coef;

  return *this;
}

TlDenseSymmetricMatrix_Scalapack& TlDenseSymmetricMatrix_Scalapack::operator*=(
    const TlDenseSymmetricMatrix_Scalapack& rhs) {
  *(dynamic_cast<TlDenseSymmetricMatrix_ImplScalapack*>(this->pImpl_)) *=
      *(dynamic_cast<TlDenseSymmetricMatrix_ImplScalapack*>(rhs.pImpl_));

  return *this;
}

// ---------------------------------------------------------------------------
// operations
// ---------------------------------------------------------------------------
const TlDenseSymmetricMatrix_Scalapack&
TlDenseSymmetricMatrix_Scalapack::dotInPlace(
    const TlDenseSymmetricMatrix_Scalapack& rhs) {
  dynamic_cast<TlDenseSymmetricMatrix_ImplScalapack*>(this->pImpl_)
      ->dotInPlace(
          *dynamic_cast<TlDenseSymmetricMatrix_ImplScalapack*>(rhs.pImpl_));

  return *this;
}

bool TlDenseSymmetricMatrix_Scalapack::eig(
    TlDenseVector_Lapack* pEigVal,
    TlDenseGeneralMatrix_Scalapack* pEigVec) const {
  TlDenseVector_ImplLapack* pImpl_eigval =
      dynamic_cast<TlDenseVector_ImplLapack*>(pEigVal->pImpl_);
  TlDenseGeneralMatrix_ImplScalapack* pImpl_eigvec =
      dynamic_cast<TlDenseGeneralMatrix_ImplScalapack*>(pEigVec->pImpl_);

  const bool answer =
      dynamic_cast<TlDenseSymmetricMatrix_ImplScalapack*>(this->pImpl_)
          ->eig(pImpl_eigval, pImpl_eigvec);
  return answer;
}

TlDenseSymmetricMatrix_Scalapack TlDenseSymmetricMatrix_Scalapack::inverse()
    const {
  TlDenseSymmetricMatrix_Scalapack answer;
  answer.pImpl_ = new TlDenseSymmetricMatrix_ImplScalapack(
      dynamic_cast<const TlDenseSymmetricMatrix_ImplScalapack*>(this->pImpl_)
          ->inverse());
  return answer;
}

bool TlDenseSymmetricMatrix_Scalapack::getSparseMatrix(TlSparseMatrix* pMatrix,
                                                       bool isFinalize) const {
  return dynamic_cast<TlDenseSymmetricMatrix_ImplScalapack*>(this->pImpl_)
      ->getSparseMatrix(pMatrix, isFinalize);
}

void TlDenseSymmetricMatrix_Scalapack::mergeSparseMatrix(
    const TlSparseMatrix& M) {
  dynamic_cast<TlDenseSymmetricMatrix_ImplScalapack*>(this->pImpl_)
      ->mergeSparseMatrix(M);
}

std::vector<TlMatrixObject::index_type>
TlDenseSymmetricMatrix_Scalapack::getRowIndexTable() const {
  return dynamic_cast<TlDenseSymmetricMatrix_ImplScalapack*>(this->pImpl_)
      ->getRowIndexTable();
}

std::vector<TlMatrixObject::index_type>
TlDenseSymmetricMatrix_Scalapack::getColIndexTable() const {
  return dynamic_cast<TlDenseSymmetricMatrix_ImplScalapack*>(this->pImpl_)
      ->getColIndexTable();
}

void TlDenseSymmetricMatrix_Scalapack::getLocalMatrix(
    TlDenseGeneralMatrixObject* pOutputMatrix) const {
  dynamic_cast<TlDenseSymmetricMatrix_ImplScalapack*>(this->pImpl_)
      ->getLocalMatrix(pOutputMatrix);
}

// ---------------------------------------------------------------------------
// I/O
// ---------------------------------------------------------------------------
bool TlDenseSymmetricMatrix_Scalapack::load(const std::string& filePath) {
  return dynamic_cast<TlDenseSymmetricMatrix_ImplScalapack*>(this->pImpl_)
      ->load(filePath);
}

bool TlDenseSymmetricMatrix_Scalapack::save(const std::string& filePath) const {
  return dynamic_cast<TlDenseSymmetricMatrix_ImplScalapack*>(this->pImpl_)
      ->save(filePath);
}

void TlDenseSymmetricMatrix_Scalapack::dump(TlDenseVector_Scalapack* v) const {
    // TODO: implement
//   dynamic_cast<TlDenseSymmetricMatrix_ImplScalapack*>(this->pImpl_)
//       ->dump(dynamic_cast<TlDenseVector_ImplScalapack*>(v->pImpl_));
}

void TlDenseSymmetricMatrix_Scalapack::restore(
    const TlDenseVector_Scalapack& v) {
    // TODO: implement
//   dynamic_cast<TlDenseSymmetricMatrix_ImplScalapack*>(this->pImpl_)
//       ->restore(*(dynamic_cast<TlDenseVector_ImplScalapack*>(v.pImpl_)));
}

// ---------------------------------------------------------------------------
// friend functions
// ---------------------------------------------------------------------------
TlDenseGeneralMatrix_Scalapack operator*(
    const TlDenseSymmetricMatrix_Scalapack& rhs1,
    const TlDenseGeneralMatrix_Scalapack& rhs2) {
  TlDenseGeneralMatrix_Scalapack answer;
  *(dynamic_cast<TlDenseGeneralMatrix_ImplScalapack*>(answer.pImpl_)) =
      *(dynamic_cast<TlDenseSymmetricMatrix_ImplScalapack*>(rhs1.pImpl_)) *
      *(dynamic_cast<TlDenseGeneralMatrix_ImplScalapack*>(rhs2.pImpl_));

  return answer;
}

TlDenseGeneralMatrix_Scalapack operator*(
    const TlDenseGeneralMatrix_Scalapack& rhs1,
    const TlDenseSymmetricMatrix_Scalapack& rhs2) {
  TlDenseGeneralMatrix_Scalapack answer;
  *(dynamic_cast<TlDenseGeneralMatrix_ImplScalapack*>(answer.pImpl_)) =
      *(dynamic_cast<TlDenseGeneralMatrix_ImplScalapack*>(rhs1.pImpl_)) *
      *(dynamic_cast<TlDenseSymmetricMatrix_ImplScalapack*>(rhs2.pImpl_));

  return answer;
}

TlDenseVector_Scalapack operator*(const TlDenseSymmetricMatrix_Scalapack& rhs1,
                                  const TlDenseVector_Scalapack& rhs2) {
  TlDenseVector_Scalapack answer;
  *(dynamic_cast<TlDenseVector_ImplScalapack*>(answer.pImpl_)) =
      *(dynamic_cast<TlDenseSymmetricMatrix_ImplScalapack*>(rhs1.pImpl_)) *
      *(dynamic_cast<TlDenseVector_ImplScalapack*>(rhs2.pImpl_));

  return answer;
}

TlDenseVector_Scalapack operator*(
    const TlDenseVector_Scalapack& rhs1,
    const TlDenseSymmetricMatrix_Scalapack& rhs2) {
  TlDenseVector_Scalapack answer;
  *(dynamic_cast<TlDenseVector_ImplScalapack*>(answer.pImpl_)) =
      *(dynamic_cast<TlDenseSymmetricMatrix_ImplScalapack*>(rhs2.pImpl_)) *
      *(dynamic_cast<TlDenseVector_ImplScalapack*>(rhs1.pImpl_));

  return answer;
}

// ---------------------------------------------------------------------------
// arithmetic
// ---------------------------------------------------------------------------
TlDenseSymmetricMatrix_Scalapack operator*(
    const double coef, const TlDenseSymmetricMatrix_Scalapack& matrix) {
  TlDenseSymmetricMatrix_Scalapack answer = matrix;
  answer *= coef;

  return answer;
}

TlDenseSymmetricMatrix_Scalapack operator*(
    const TlDenseSymmetricMatrix_Scalapack& matrix, const double coef) {
  return coef * matrix;
}
