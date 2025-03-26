#include "tl_dense_vector_impl_eigen_float.h"

#include <iostream>

#ifdef HAVE_VIENNACL
#define VIENNACL_HAVE_EIGEN
#include <viennacl/vector.hpp>

#include "tl_dense_vector_impl_viennacl_float.h"
#endif  // HAVE_VIENNACL

// ---------------------------------------------------------------------------
// constructor & destructor
// ---------------------------------------------------------------------------
TlDenseVector_ImplEigenFloat::TlDenseVector_ImplEigenFloat(TlDenseVectorObject::index_type dim)
    : vector_(VectorDataType::Zero(dim)) {}

TlDenseVector_ImplEigenFloat::TlDenseVector_ImplEigenFloat(const TlDenseVector_ImplEigenFloat& rhs) {
    this->vector_ = rhs.vector_;
}

TlDenseVector_ImplEigenFloat::TlDenseVector_ImplEigenFloat(const std::vector<double>& rhs)
    : vector_() {
    this->vector_ = Eigen::Map<const Eigen::VectorXd>(rhs.data(), rhs.size()).cast<float>();
}

TlDenseVector_ImplEigenFloat::TlDenseVector_ImplEigenFloat(const double* p, const TlDenseVectorObject::size_type size)
    : vector_() {
    this->vector_ = Eigen::Map<const Eigen::VectorXd>(p, size).cast<float>();
}

TlDenseVector_ImplEigenFloat::TlDenseVector_ImplEigenFloat(const VectorDataType& rhs) {
    this->vector_ = rhs;
}

#ifdef HAVE_VIENNACL
TlDenseVector_ImplEigenFloat::TlDenseVector_ImplEigenFloat(const TlDenseVector_ImplViennaCLFloat& rhs)
    : vector_(VectorDataType::Zero(rhs.getSize())) {
    viennacl::copy(rhs.vector_, this->vector_);
}
#endif  // HAVE_VIENNACL

TlDenseVector_ImplEigenFloat::operator std::vector<double>() const {
    const std::size_t size = this->getSize();
    std::vector<float> tmp(size);
    MapType(&(tmp[0]), size) = this->vector_;

    return std::vector<double>(tmp.begin(), tmp.end());
}

TlDenseVector_ImplEigenFloat::~TlDenseVector_ImplEigenFloat() {}

// ---------------------------------------------------------------------------
// properties
// ---------------------------------------------------------------------------
TlDenseVectorObject::size_type TlDenseVector_ImplEigenFloat::getSize() const {
    return this->vector_.rows();
}

void TlDenseVector_ImplEigenFloat::resize(const TlDenseVectorObject::index_type newSize) {
    this->vector_.conservativeResizeLike(VectorDataType::Zero(newSize, 1));
}

double TlDenseVector_ImplEigenFloat::get(const TlDenseVectorObject::index_type i) const {
    return static_cast<double>(this->vector_.coeff(i));
}

void TlDenseVector_ImplEigenFloat::set(const TlDenseVectorObject::index_type i, const double value) {
    this->vector_.coeffRef(i) = static_cast<float>(value);
}

void TlDenseVector_ImplEigenFloat::add(const TlDenseVectorObject::index_type i, const double value) {
    this->vector_.coeffRef(i) += static_cast<float>(value);
}

void TlDenseVector_ImplEigenFloat::mul(const TlDenseVectorObject::index_type i, const double value) {
    this->vector_.coeffRef(i) *= static_cast<float>(value);
}

// ---------------------------------------------------------------------------
// operators
// ---------------------------------------------------------------------------
TlDenseVector_ImplEigenFloat& TlDenseVector_ImplEigenFloat::operator=(
    const TlDenseVector_ImplEigenFloat& rhs) {
    if (this != &rhs) {
        this->vector_ = rhs.vector_;
    }

    return (*this);
}

TlDenseVector_ImplEigenFloat& TlDenseVector_ImplEigenFloat::operator+=(
    const TlDenseVector_ImplEigenFloat& rhs) {
    this->vector_ += rhs.vector_;

    return *this;
}

TlDenseVector_ImplEigenFloat& TlDenseVector_ImplEigenFloat::operator-=(
    const TlDenseVector_ImplEigenFloat& rhs) {
    this->vector_ -= rhs.vector_;

    return *this;
}

TlDenseVector_ImplEigenFloat& TlDenseVector_ImplEigenFloat::operator*=(const double rhs) {
    this->vector_ *= rhs;

    return *this;
}

TlDenseVector_ImplEigenFloat& TlDenseVector_ImplEigenFloat::operator/=(const double rhs) {
    return this->operator*=(1.0 / rhs);
}

double TlDenseVector_ImplEigenFloat::operator*(
    const TlDenseVector_ImplEigenFloat& rhs) const {
    return this->dot(rhs);
}

// ---------------------------------------------------------------------------
// operations
// ---------------------------------------------------------------------------
double TlDenseVector_ImplEigenFloat::sum() const {
    return this->vector_.array().sum();
}

void TlDenseVector_ImplEigenFloat::sortByGreater() {
    std::sort(this->vector_.data(), this->vector_.data() + this->getSize(),
              std::greater<float>());
}

double TlDenseVector_ImplEigenFloat::dot(const TlDenseVector_ImplEigenFloat& rhs) const {
    assert(this->getSize() == rhs.getSize());
    const double answer = this->vector_.dot(rhs.vector_);

    return answer;
}

TlDenseVector_ImplEigenFloat& TlDenseVector_ImplEigenFloat::dotInPlace(
    const TlDenseVector_ImplEigenFloat& rhs) {
    assert(this->getSize() == rhs.getSize());
    this->vector_.array() *= rhs.vector_.array();

    return *this;
}

// ---------------------------------------------------------------------------
// others
// ---------------------------------------------------------------------------
TlDenseVector_ImplEigenFloat operator+(const TlDenseVector_ImplEigenFloat& rhs1,
                                       const TlDenseVector_ImplEigenFloat& rhs2) {
    TlDenseVector_ImplEigenFloat answer = rhs1;
    answer += rhs2;

    return answer;
}

TlDenseVector_ImplEigenFloat operator-(const TlDenseVector_ImplEigenFloat& rhs1,
                                       const TlDenseVector_ImplEigenFloat& rhs2) {
    TlDenseVector_ImplEigenFloat answer = rhs1;
    answer -= rhs2;

    return answer;
}

TlDenseVector_ImplEigenFloat operator*(const TlDenseVector_ImplEigenFloat& rhs1,
                                       const double rhs2) {
    TlDenseVector_ImplEigenFloat answer = rhs1;

    answer *= rhs2;

    return answer;
}

TlDenseVector_ImplEigenFloat operator*(const double rhs1,
                                       const TlDenseVector_ImplEigenFloat& rhs2) {
    return (rhs2 * rhs1);
}
