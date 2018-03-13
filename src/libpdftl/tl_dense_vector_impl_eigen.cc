#include "tl_dense_vector_impl_eigen.h"

// ---------------------------------------------------------------------------
// constructor & destructor
// ---------------------------------------------------------------------------
TlDenseVector_ImplEigen::TlDenseVector_ImplEigen(
    TlDenseVectorObject::index_type dim)
    : vector_(VectorDataType::Zero(dim)) {}

TlDenseVector_ImplEigen::TlDenseVector_ImplEigen(
    const TlDenseVector_ImplEigen& rhs) {
  this->vector_ = rhs.vector_;
}

TlDenseVector_ImplEigen::TlDenseVector_ImplEigen(
    const double* p, const TlDenseVectorObject::size_type size)
    : vector_(MapTypeConst(p, size)) {}

TlDenseVector_ImplEigen::TlDenseVector_ImplEigen(const VectorDataType& rhs) {
  this->vector_ = rhs;
}

TlDenseVector_ImplEigen::~TlDenseVector_ImplEigen() {}

// ---------------------------------------------------------------------------
// properties
// ---------------------------------------------------------------------------
TlDenseVectorObject::size_type TlDenseVector_ImplEigen::getSize() const {
  return this->vector_.rows();
}

void TlDenseVector_ImplEigen::resize(
    const TlDenseVectorObject::index_type newSize) {
  this->vector_.conservativeResizeLike(VectorDataType::Zero(newSize, 1));
}

double TlDenseVector_ImplEigen::get(
    const TlDenseVectorObject::index_type i) const {
  return this->vector_.coeff(i);
}

void TlDenseVector_ImplEigen::set(const TlDenseVectorObject::index_type i,
                                  const double value) {
  this->vector_.coeffRef(i) = value;
}

void TlDenseVector_ImplEigen::add(const TlDenseVectorObject::index_type i,
                                  const double value) {
  this->vector_.coeffRef(i) += value;
}

// ---------------------------------------------------------------------------
// operators
// ---------------------------------------------------------------------------
TlDenseVector_ImplEigen& TlDenseVector_ImplEigen::operator=(
    const TlDenseVector_ImplEigen& rhs) {
  if (this != &rhs) {
    this->vector_ = rhs.vector_;
  }

  return (*this);
}

TlDenseVector_ImplEigen& TlDenseVector_ImplEigen::operator+=(
    const TlDenseVector_ImplEigen& rhs) {
  this->vector_ += rhs.vector_;

  return *this;
}

TlDenseVector_ImplEigen& TlDenseVector_ImplEigen::operator-=(
    const TlDenseVector_ImplEigen& rhs) {
  this->vector_ -= rhs.vector_;

  return *this;
}

TlDenseVector_ImplEigen& TlDenseVector_ImplEigen::operator*=(const double rhs) {
  const TlDenseVectorObject::size_type size = this->getSize();
  this->vector_ *= rhs;

  return *this;
}

TlDenseVector_ImplEigen& TlDenseVector_ImplEigen::operator/=(const double rhs) {
  return this->operator*=(1.0 / rhs);
}

// ---------------------------------------------------------------------------
// operations
// ---------------------------------------------------------------------------
double TlDenseVector_ImplEigen::sum() const {
  return this->vector_.array().sum();
}

void TlDenseVector_ImplEigen::sortByGreater() {
  std::sort(this->vector_.data(), this->vector_.data() + this->getSize(),
            std::greater<double>());
}

TlDenseVector_ImplEigen& TlDenseVector_ImplEigen::dotInPlace(
    const TlDenseVector_ImplEigen& rhs) {
  assert(this->getSize() == rhs.getSize());
  this->vector_.array() *= rhs.vector_.array();

  return *this;
}

// ---------------------------------------------------------------------------
// others
// ---------------------------------------------------------------------------
TlDenseVector_ImplEigen operator+(const TlDenseVector_ImplEigen& rhs1,
                                  const TlDenseVector_ImplEigen& rhs2) {
  TlDenseVector_ImplEigen answer = rhs1;
  answer += rhs2;

  return answer;
}

TlDenseVector_ImplEigen operator-(const TlDenseVector_ImplEigen& rhs1,
                                  const TlDenseVector_ImplEigen& rhs2) {
  TlDenseVector_ImplEigen answer = rhs1;
  answer -= rhs2;

  return answer;
}

TlDenseVector_ImplEigen operator*(const TlDenseVector_ImplEigen& rhs1,
                                  const double rhs2) {
  TlDenseVector_ImplEigen answer = rhs1;

  answer *= rhs2;

  return answer;
}

TlDenseVector_ImplEigen operator*(const double rhs1,
                                  const TlDenseVector_ImplEigen& rhs2) {
  return (rhs2 * rhs1);
}
