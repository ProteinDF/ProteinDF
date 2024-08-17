#ifdef HAVE_CONFIG_H
#include "config.h"
#endif  // HAVE_CONFIG_H

#include <algorithm>
#include <cassert>
#include <functional>
#include <vector>

#define VIENNACL_HAVE_EIGEN 1
#include <viennacl/linalg/sum.hpp>
#include <viennacl/vector.hpp>
#include <viennacl/vector_proxy.hpp>

#include "tl_dense_vector_impl_eigen_float.h"
#include "tl_dense_vector_impl_viennacl_float.h"

#ifdef HAVE_EIGEN
#include "tl_dense_vector_eigen_float.h"
#endif  // HAVE_EIGEN

// ---------------------------------------------------------------------------
// constructor & destructor
// ---------------------------------------------------------------------------
TlDenseVector_ImplViennaCLFloat::TlDenseVector_ImplViennaCLFloat(TlDenseVectorObject::index_type dim)
    : vector_(dim) {}

TlDenseVector_ImplViennaCLFloat::TlDenseVector_ImplViennaCLFloat(const TlDenseVector_ImplViennaCLFloat& rhs) {
    this->vector_ = rhs.vector_;
}

TlDenseVector_ImplViennaCLFloat::TlDenseVector_ImplViennaCLFloat(const std::vector<double>& rhs)
    : vector_(rhs.size()) {
    viennacl::copy(&(rhs[0]), &(rhs[0]) + rhs.size(), this->vector_.begin());
}

#ifdef HAVE_EIGEN
TlDenseVector_ImplViennaCLFloat::TlDenseVector_ImplViennaCLFloat(const TlDenseVector_ImplEigenFloat& rhs)
    : vector_(rhs.getSize()) {
    viennacl::copy(rhs.vector_, this->vector_);
}
#endif  // HAVE_EIGEN

TlDenseVector_ImplViennaCLFloat::operator std::vector<double>() const {
    const TlDenseVectorObject::index_type size = this->getSize();
    std::vector<double> answer(size);
    viennacl::copy(this->vector_.begin(), this->vector_.end(), answer.begin());
    return answer;
}

TlDenseVector_ImplViennaCLFloat::~TlDenseVector_ImplViennaCLFloat() {}

// ---------------------------------------------------------------------------
// properties
// ---------------------------------------------------------------------------
TlDenseVectorObject::size_type TlDenseVector_ImplViennaCLFloat::getSize() const {
    return this->vector_.size();
}

void TlDenseVector_ImplViennaCLFloat::resize(
    TlDenseVectorObject::index_type newSize) {
    this->vector_.resize(newSize, true);
}

double TlDenseVector_ImplViennaCLFloat::get(const TlDenseVectorObject::index_type i) const {
    return static_cast<double>(this->vector_[i]);
}

void TlDenseVector_ImplViennaCLFloat::set(const TlDenseVectorObject::index_type i,
                                          const double value) {
    this->vector_[i] = static_cast<float>(value);
}

void TlDenseVector_ImplViennaCLFloat::add(const TlDenseVectorObject::index_type i,
                                          const double value) {
    this->vector_[i] += static_cast<double>(value);
}

void TlDenseVector_ImplViennaCLFloat::mul(const TlDenseVectorObject::index_type i,
                                          const double value) {
    this->vector_[i] *= static_cast<double>(value);
}

// ---------------------------------------------------------------------------
// operators
// ---------------------------------------------------------------------------
TlDenseVector_ImplViennaCLFloat& TlDenseVector_ImplViennaCLFloat::operator=(
    const TlDenseVector_ImplViennaCLFloat& rhs) {
    if (this != &rhs) {
        this->vector_ = rhs.vector_;
    }

    return *this;
}

TlDenseVector_ImplViennaCLFloat& TlDenseVector_ImplViennaCLFloat::operator+=(
    const TlDenseVector_ImplViennaCLFloat& rhs) {
    this->vector_ += rhs.vector_;

    return *this;
}

TlDenseVector_ImplViennaCLFloat& TlDenseVector_ImplViennaCLFloat::operator-=(
    const TlDenseVector_ImplViennaCLFloat& rhs) {
    this->vector_ -= rhs.vector_;

    return *this;
}

TlDenseVector_ImplViennaCLFloat& TlDenseVector_ImplViennaCLFloat::operator*=(
    const float rhs) {
    this->vector_ *= rhs;

    return *this;
}
TlDenseVector_ImplViennaCLFloat& TlDenseVector_ImplViennaCLFloat::operator/=(
    const float rhs) {
    return this->operator*=(1.0 / rhs);
}

// ---------------------------------------------------------------------------
// operations
// ---------------------------------------------------------------------------
double TlDenseVector_ImplViennaCLFloat::sum() const {
    const float sum = viennacl::linalg::sum(this->vector_);
    return static_cast<double>(sum);
}

void TlDenseVector_ImplViennaCLFloat::sortByGreater() {
    // calculate on CPU
    std::vector<float> vector_cpu(this->getSize());
    viennacl::copy(this->vector_, vector_cpu);

    std::sort(&(vector_cpu[0]), &(vector_cpu[0]) + this->getSize(),
              std::greater<float>());
    viennacl::copy(vector_cpu, this->vector_);
}

TlDenseVector_ImplViennaCLFloat& TlDenseVector_ImplViennaCLFloat::dotInPlace(
    const TlDenseVector_ImplViennaCLFloat& rhs) {
    const VectorDataType tmp =
        viennacl::linalg::element_prod(this->vector_, rhs.vector_);
    this->vector_ = tmp;

    return *this;
}

TlDenseVector_ImplViennaCLFloat& TlDenseVector_ImplViennaCLFloat::reverse() {
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
TlDenseVector_ImplViennaCLFloat operator+(const TlDenseVector_ImplViennaCLFloat& rhs1,
                                          const TlDenseVector_ImplViennaCLFloat& rhs2) {
    TlDenseVector_ImplViennaCLFloat answer = rhs1;
    answer += rhs2;

    return answer;
}

TlDenseVector_ImplViennaCLFloat operator-(const TlDenseVector_ImplViennaCLFloat& rhs1,
                                          const TlDenseVector_ImplViennaCLFloat& rhs2) {
    TlDenseVector_ImplViennaCLFloat answer = rhs1;
    answer -= rhs2;

    return answer;
}

TlDenseVector_ImplViennaCLFloat operator*(const TlDenseVector_ImplViennaCLFloat& rhs1,
                                          const float rhs2) {
    TlDenseVector_ImplViennaCLFloat answer = rhs1;

    answer *= rhs2;

    return answer;
}

TlDenseVector_ImplViennaCLFloat operator*(const float rhs1,
                                          const TlDenseVector_ImplViennaCLFloat& rhs2) {
    return (rhs2 * rhs1);
}
