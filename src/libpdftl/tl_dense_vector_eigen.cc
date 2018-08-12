#include "tl_dense_vector_eigen.h"
#include "tl_dense_vector_impl_eigen.h"

// ---------------------------------------------------------------------------
// constructor & destructor
// ---------------------------------------------------------------------------
TlDenseVector_Eigen::TlDenseVector_Eigen(TlDenseVectorObject::index_type size) {
  this->pImpl_ = new TlDenseVector_ImplEigen(size);
}

TlDenseVector_Eigen::TlDenseVector_Eigen(const TlDenseVector_Eigen& rhs) {
  this->pImpl_ = new TlDenseVector_ImplEigen(
      *dynamic_cast<const TlDenseVector_ImplEigen*>(rhs.pImpl_));
}

TlDenseVector_Eigen::TlDenseVector_Eigen(const double* p, size_type size) {
  this->pImpl_ = new TlDenseVector_ImplEigen(p, size);
}

TlDenseVector_Eigen::~TlDenseVector_Eigen() {
  delete this->pImpl_;
  this->pImpl_ = NULL;
}

// ---------------------------------------------------------------------------
// operators
// ---------------------------------------------------------------------------
TlDenseVector_Eigen& TlDenseVector_Eigen::operator=(
    const TlDenseVector_Eigen& rhs) {
  if (this != &rhs) {
    delete this->pImpl_;
    this->pImpl_ = NULL;

    this->pImpl_ = new TlDenseVector_ImplEigen(
        *dynamic_cast<const TlDenseVector_ImplEigen*>(rhs.pImpl_));
  }

  return *this;
}

TlDenseVector_Eigen& TlDenseVector_Eigen::operator+=(
    const TlDenseVector_Eigen& rhs) {
  dynamic_cast<TlDenseVector_ImplEigen*>(this->pImpl_)
      ->
      operator+=(*dynamic_cast<TlDenseVector_ImplEigen*>(rhs.pImpl_));

  return *this;
}

TlDenseVector_Eigen& TlDenseVector_Eigen::operator-=(
    const TlDenseVector_Eigen& rhs) {
  dynamic_cast<TlDenseVector_ImplEigen*>(this->pImpl_)
      ->
      operator-=(*dynamic_cast<TlDenseVector_ImplEigen*>(rhs.pImpl_));

  return *this;
}

TlDenseVector_Eigen& TlDenseVector_Eigen::operator*=(const double rhs) {
  dynamic_cast<TlDenseVector_ImplEigen*>(this->pImpl_)->operator*=(rhs);

  return *this;
}

TlDenseVector_Eigen& TlDenseVector_Eigen::operator/=(const double rhs) {
  dynamic_cast<TlDenseVector_ImplEigen*>(this->pImpl_)->operator/=(rhs);

  return *this;
}

// ---------------------------------------------------------------------------
// operations
// ---------------------------------------------------------------------------
TlDenseVector_Eigen& TlDenseVector_Eigen::dotInPlace(
    const TlDenseVector_Eigen& rhs) {
  dynamic_cast<TlDenseVector_ImplEigen*>(this->pImpl_)
      ->dotInPlace(*dynamic_cast<TlDenseVector_ImplEigen*>(rhs.pImpl_));

  return *this;
}

// ---------------------------------------------------------------------------
// others
// ---------------------------------------------------------------------------
TlDenseVector_Eigen operator+(const TlDenseVector_Eigen& rhs1,
                              const TlDenseVector_Eigen& rhs2) {
  TlDenseVector_Eigen answer = rhs1;
  answer += rhs2;

  return answer;
}

TlDenseVector_Eigen operator-(const TlDenseVector_Eigen& rhs1,
                              const TlDenseVector_Eigen& rhs2) {
  TlDenseVector_Eigen answer = rhs1;
  answer -= rhs2;

  return answer;
}

TlDenseVector_Eigen operator*(const TlDenseVector_Eigen& rhs1,
                              const double rhs2) {
  TlDenseVector_Eigen answer = rhs1;
  answer *= rhs2;

  return answer;
}

TlDenseVector_Eigen operator*(const double rhs1,
                              const TlDenseVector_Eigen& rhs2) {
  return rhs2 * rhs1;
}
