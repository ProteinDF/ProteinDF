#include "tl_dense_vector_lapack.h"
#include "tl_dense_vector_impl_lapack.h"

// ---------------------------------------------------------------------------
// constructor & destructor
// ---------------------------------------------------------------------------
TlDenseVector_Lapack::TlDenseVector_Lapack(
    TlDenseVectorObject::index_type size) {
  this->pImpl_ = new TlDenseVector_ImplLapack(size);
}

TlDenseVector_Lapack::TlDenseVector_Lapack(const TlDenseVector_Lapack& rhs) {
  this->pImpl_ = new TlDenseVector_ImplLapack(
      *dynamic_cast<const TlDenseVector_ImplLapack*>(rhs.pImpl_));
}

TlDenseVector_Lapack::TlDenseVector_Lapack(const std::vector<double>& rhs) {
  this->pImpl_ = new TlDenseVector_ImplLapack(rhs);
}

// TlDenseVector_Lapack::TlDenseVector_Lapack(const double* p, size_type size) {
//   this->pImpl_ = new TlDenseVector_ImplLapack(p, size);
// }

TlDenseVector_Lapack::~TlDenseVector_Lapack() {
  delete this->pImpl_;
  this->pImpl_ = NULL;
}

// ---------------------------------------------------------------------------
// operators
// ---------------------------------------------------------------------------
TlDenseVector_Lapack& TlDenseVector_Lapack::operator=(
    const TlDenseVector_Lapack& rhs) {
  if (this != &rhs) {
    delete this->pImpl_;
    this->pImpl_ = NULL;

    this->pImpl_ = new TlDenseVector_ImplLapack(
        *dynamic_cast<const TlDenseVector_ImplLapack*>(rhs.pImpl_));
  }

  return *this;
}

TlDenseVector_Lapack& TlDenseVector_Lapack::operator+=(
    const TlDenseVector_Lapack& rhs) {
  dynamic_cast<TlDenseVector_ImplLapack*>(this->pImpl_)
      ->
      operator+=(*dynamic_cast<TlDenseVector_ImplLapack*>(rhs.pImpl_));

  return *this;
}

TlDenseVector_Lapack& TlDenseVector_Lapack::operator-=(
    const TlDenseVector_Lapack& rhs) {
  dynamic_cast<TlDenseVector_ImplLapack*>(this->pImpl_)
      ->
      operator-=(*dynamic_cast<TlDenseVector_ImplLapack*>(rhs.pImpl_));

  return *this;
}

TlDenseVector_Lapack& TlDenseVector_Lapack::operator*=(const double rhs) {
  dynamic_cast<TlDenseVector_ImplLapack*>(this->pImpl_)->operator*=(rhs);

  return *this;
}

TlDenseVector_Lapack& TlDenseVector_Lapack::operator/=(const double rhs) {
  dynamic_cast<TlDenseVector_ImplLapack*>(this->pImpl_)->operator/=(rhs);

  return *this;
}

double TlDenseVector_Lapack::operator*(const TlDenseVector_Lapack& rhs) const {
  return dynamic_cast<const TlDenseVector_ImplLapack*>(this->pImpl_)
      ->
      operator*(*dynamic_cast<const TlDenseVector_ImplLapack*>(rhs.pImpl_));
}

// ---------------------------------------------------------------------------
// operations
// ---------------------------------------------------------------------------
TlDenseVector_Lapack& TlDenseVector_Lapack::dotInPlace(
    const TlDenseVector_Lapack& rhs) {
  dynamic_cast<TlDenseVector_ImplLapack*>(this->pImpl_)
      ->dotInPlace(*dynamic_cast<TlDenseVector_ImplLapack*>(rhs.pImpl_));

  return *this;
}

double* TlDenseVector_Lapack::data() {
  return dynamic_cast<TlDenseVector_ImplLapack*>(this->pImpl_)->data();
}

const double* TlDenseVector_Lapack::data() const {
  return dynamic_cast<TlDenseVector_ImplLapack*>(this->pImpl_)->data();
}

// ---------------------------------------------------------------------------
// others
// ---------------------------------------------------------------------------
TlDenseVector_Lapack operator+(const TlDenseVector_Lapack& rhs1,
                               const TlDenseVector_Lapack& rhs2) {
  TlDenseVector_Lapack answer = rhs1;
  answer += rhs2;

  return answer;
}

TlDenseVector_Lapack operator-(const TlDenseVector_Lapack& rhs1,
                               const TlDenseVector_Lapack& rhs2) {
  TlDenseVector_Lapack answer = rhs1;
  answer -= rhs2;

  return answer;
}

TlDenseVector_Lapack operator*(const TlDenseVector_Lapack& rhs1,
                               const double rhs2) {
  TlDenseVector_Lapack answer = rhs1;
  answer *= rhs2;

  return answer;
}

TlDenseVector_Lapack operator*(const double rhs1,
                               const TlDenseVector_Lapack& rhs2) {
  return rhs2 * rhs1;
}
