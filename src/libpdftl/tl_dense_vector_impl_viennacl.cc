#include "tl_dense_vector_impl_viennacl.h"
#include <algorithm>
#include <cassert>
#include <functional>
#include <vector>
#include "viennacl/linalg/sum.hpp"
#include "viennacl/vector.hpp"
#include "viennacl/vector_proxy.hpp"
#include "tl_dense_vector_impl_eigen.h"

#ifdef HAVE_EIGEN
#include "tl_dense_vector_eigen.h"
#endif  // HAVE_EIGEN

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

TlDenseVector_ImplViennaCL::TlDenseVector_ImplViennaCL(
    const std::vector<double>& rhs)
    : vector_(rhs.size()) {
  viennacl::copy(&(rhs[0]), &(rhs[0]) + rhs.size(), this->vector_.begin());
}

#ifdef HAVE_EIGEN
TlDenseVector_ImplViennaCL::TlDenseVector_ImplViennaCL(
    const TlDenseVector_ImplEigen& rhs) : vector_(rhs.getSize()) {
  viennacl::copy(rhs.vector_, this->vector_);
}
#endif  // HAVE_EIGEN

TlDenseVector_ImplViennaCL::operator std::vector<double>() const {
  std::vector<double> answer(this->getSize());
  viennacl::copy(this->vector_, answer);
  return answer;
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

TlDenseVector_ImplViennaCL& TlDenseVector_ImplViennaCL::reverse() {
  const TlVectorObject::index_type size = this->getSize();
  viennacl::slice s(size - 1, -1, size);

  viennacl::vector_slice<VectorDataType> vs(this->vector_, s);
  const VectorDataType tmp = vs;
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
