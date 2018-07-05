#include "tl_dense_vector_impl_viennacl.h"
#include <algorithm>
#include <cassert>
#include <functional>
#include <vector>
#include "viennacl/linalg/sum.hpp"
#include "viennacl/vector.hpp"

// ---------------------------------------------------------------------------
// constructor & destructor
// ---------------------------------------------------------------------------
TlDenseVector_ImplViennaCL::TlDenseVector_ImplViennaCL(
    TlDenseVectorObject::index_type dim)
    : vector_(dim) {}

TlDenseVector_ImplViennaCL::TlDenseVector_ImplViennaCL(
    const TlDenseVector_ImplViennaCL& rhs) {
  this->vector_ = rhs.vector_;
}

TlDenseVector_ImplViennaCL::~TlDenseVector_ImplViennaCL() {}

// ---------------------------------------------------------------------------
// properties
// ---------------------------------------------------------------------------
TlDenseVectorObject::size_type TlDenseVector_ImplViennaCL::getSize() const {
  return this->vector_.size();
}

void TlDenseVector_ImplViennaCL::resize(
    TlDenseVectorObject::index_type newSize) {
  this->vector_.resize(newSize, true);
}

double TlDenseVector_ImplViennaCL::get(
    const TlDenseVectorObject::index_type i) const {
  return this->vector_[i];
}

void TlDenseVector_ImplViennaCL::set(const TlDenseVectorObject::index_type i,
                                     const double value) {
  this->vector_[i] = value;
}

void TlDenseVector_ImplViennaCL::add(const TlDenseVectorObject::index_type i,
                                     const double value) {
  this->vector_[i] += value;
}

void TlDenseVector_ImplViennaCL::mul(const TlDenseVectorObject::index_type i,
                                     const double value) {
  this->vector_[i] *= value;
}

// ---------------------------------------------------------------------------
// operators
// ---------------------------------------------------------------------------
TlDenseVector_ImplViennaCL& TlDenseVector_ImplViennaCL::operator=(
    const TlDenseVector_ImplViennaCL& rhs) {
  if (this != &rhs) {
    this->vector_ = rhs.vector_;
  }

  return *this;
}

TlDenseVector_ImplViennaCL& TlDenseVector_ImplViennaCL::operator+=(
    const TlDenseVector_ImplViennaCL& rhs) {
  this->vector_ += rhs.vector_;

  return *this;
}

TlDenseVector_ImplViennaCL& TlDenseVector_ImplViennaCL::operator-=(
    const TlDenseVector_ImplViennaCL& rhs) {
  this->vector_ -= rhs.vector_;

  return *this;
}

TlDenseVector_ImplViennaCL& TlDenseVector_ImplViennaCL::operator*=(
    const double rhs) {
  this->vector_ *= rhs;

  return *this;
}
TlDenseVector_ImplViennaCL& TlDenseVector_ImplViennaCL::operator/=(
    const double rhs) {
  return this->operator*=(1.0 / rhs);
}

// ---------------------------------------------------------------------------
// operations
// ---------------------------------------------------------------------------
double TlDenseVector_ImplViennaCL::sum() const {
  const double sum = viennacl::linalg::sum(this->vector_);
  return sum;
}

void TlDenseVector_ImplViennaCL::sortByGreater() {
  // calculate on CPU
  std::vector<double> vector_cpu(this->getSize());
  viennacl::copy(this->vector_, vector_cpu);

  std::sort(&(vector_cpu[0]), &(vector_cpu[0]) + this->getSize(),
            std::greater<double>());
  viennacl::copy(vector_cpu, this->vector_);
}

TlDenseVector_ImplViennaCL& TlDenseVector_ImplViennaCL::dotInPlace(
    const TlDenseVector_ImplViennaCL& rhs) {
  const VectorDataType tmp =
      viennacl::linalg::element_prod(this->vector_, rhs.vector_);
  this->vector_ = tmp;

  return *this;
}

// ---------------------------------------------------------------------------
// others
// ---------------------------------------------------------------------------
TlDenseVector_ImplViennaCL operator+(const TlDenseVector_ImplViennaCL& rhs1,
                                     const TlDenseVector_ImplViennaCL& rhs2) {
  TlDenseVector_ImplViennaCL answer = rhs1;
  answer += rhs2;

  return answer;
}

TlDenseVector_ImplViennaCL operator-(const TlDenseVector_ImplViennaCL& rhs1,
                                     const TlDenseVector_ImplViennaCL& rhs2) {
  TlDenseVector_ImplViennaCL answer = rhs1;
  answer -= rhs2;

  return answer;
}

TlDenseVector_ImplViennaCL operator*(const TlDenseVector_ImplViennaCL& rhs1,
                                     const double rhs2) {
  TlDenseVector_ImplViennaCL answer = rhs1;

  answer *= rhs2;

  return answer;
}

TlDenseVector_ImplViennaCL operator*(const double rhs1,
                                     const TlDenseVector_ImplViennaCL& rhs2) {
  return (rhs2 * rhs1);
}
