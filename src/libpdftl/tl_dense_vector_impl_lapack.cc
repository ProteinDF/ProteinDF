#include <algorithm>
#include <cassert>
#include <functional>
#include <numeric>

#include "TlUtils.h"
#include "lapack.h"
#include "tl_dense_vector_impl_lapack.h"

TlDenseVector_ImplLapack::TlDenseVector_ImplLapack(
    const TlDenseVectorObject::index_type size)
    : size_(size), vector_(NULL) {
  this->initialize(true);
}

TlDenseVector_ImplLapack::TlDenseVector_ImplLapack(
    const TlDenseVector_ImplLapack& rhs)
    : size_(rhs.getSize()), vector_(NULL) {
  this->initialize(false);
  std::copy(rhs.vector_, rhs.vector_ + rhs.getSize(), this->vector_);
}

TlDenseVector_ImplLapack::TlDenseVector_ImplLapack(
    const std::vector<double>& rhs)
    : size_(rhs.size()), vector_(NULL) {
  this->initialize(false);
  std::copy(rhs.begin(), rhs.end(), this->vector_);
}

TlDenseVector_ImplLapack::TlDenseVector_ImplLapack(
    const std::valarray<double>& rhs)
    : size_(rhs.size()), vector_(NULL) {
  this->initialize(false);
  std::copy(&(rhs[0]), &(rhs[0]) + rhs.size(), this->vector_);
}

TlDenseVector_ImplLapack::operator std::vector<double>() const {
  std::vector<double> answer(this->vector_, this->vector_ + this->getSize());
  return answer;
}

TlDenseVector_ImplLapack::~TlDenseVector_ImplLapack() {
  delete[] this->vector_;
  this->vector_ = NULL;
}

// ---------------------------------------------------------------------------
// properties
// ---------------------------------------------------------------------------
TlDenseVectorObject::index_type TlDenseVector_ImplLapack::getSize() const {
  return this->size_;
}

void TlDenseVector_ImplLapack::resize(
    const TlDenseVectorObject::index_type newSize) {
  assert(newSize > 0);
  const TlDenseVector_ImplLapack oldVector(*this);

  // destroy object
  delete[] this->vector_;
  this->vector_ = NULL;

  // initialize
  this->size_ = newSize;
  this->initialize(true);

  // copy old data
  const TlDenseVectorObject::index_type maxSizeForCopy =
      std::min(oldVector.getSize(), newSize);
#pragma omp parallel for
  for (TlDenseVectorObject::index_type i = 0; i < maxSizeForCopy; ++i) {
    this->set(i, oldVector.get(i));
  }
}

double TlDenseVector_ImplLapack::get(
    const TlDenseVectorObject::index_type i) const {
  return this->vector_[i];
}

void TlDenseVector_ImplLapack::set(const TlDenseVectorObject::index_type i,
                                   const double value) {
#pragma omp critical(TlDenseGeneralVector_ImplLapack__set)
  { this->vector_[i] = value; }
}

void TlDenseVector_ImplLapack::add(const TlDenseVectorObject::index_type i,
                                   const double value) {
#pragma omp atomic
  this->vector_[i] += value;
}

void TlDenseVector_ImplLapack::mul(const TlDenseVectorObject::index_type i,
                                   const double value) {
#pragma omp atomic
  this->vector_[i] *= value;
}

// ---------------------------------------------------------------------------
// operators
// ---------------------------------------------------------------------------
TlDenseVector_ImplLapack& TlDenseVector_ImplLapack::operator=(
    const TlDenseVector_ImplLapack& rhs) {
  if (this != &rhs) {
    delete[] this->vector_;
    this->vector_ = NULL;

    this->size_ = rhs.getSize();
    this->initialize(false);
    std::copy(rhs.vector_, rhs.vector_ + rhs.getSize(), this->vector_);
  }

  return (*this);
}

TlDenseVector_ImplLapack& TlDenseVector_ImplLapack::operator+=(
    const TlDenseVector_ImplLapack& rhs) {
  const TlDenseVectorObject::size_type size = this->getSize();
#pragma omp parallel for
  for (TlDenseVectorObject::size_type i = 0; i < size; ++i) {
    this->vector_[i] += rhs.get(i);
  }

  return *this;
}

TlDenseVector_ImplLapack& TlDenseVector_ImplLapack::operator-=(
    const TlDenseVector_ImplLapack& rhs) {
  const TlDenseVectorObject::size_type size = this->getSize();
#pragma omp parallel for
  for (TlDenseVectorObject::size_type i = 0; i < size; ++i) {
    this->vector_[i] -= rhs.get(i);
  }

  return *this;
}

TlDenseVector_ImplLapack& TlDenseVector_ImplLapack::operator*=(
    const double rhs) {
  const TlDenseVectorObject::size_type size = this->getSize();
#pragma omp parallel for
  for (TlDenseVectorObject::size_type i = 0; i < size; ++i) {
    this->vector_[i] *= rhs;
  }

  return *this;
}

TlDenseVector_ImplLapack& TlDenseVector_ImplLapack::operator/=(
    const double rhs) {
  return this->operator*=(1.0 / rhs);
}

double TlDenseVector_ImplLapack::operator*(
    const TlDenseVector_ImplLapack& rhs) const {
      return this->dot(rhs);
}

// ---------------------------------------------------------------------------
// operations
// ---------------------------------------------------------------------------
double* TlDenseVector_ImplLapack::data() { return this->vector_; }

const double* TlDenseVector_ImplLapack::data() const { return this->vector_; }

double TlDenseVector_ImplLapack::sum() const {
  return std::accumulate(this->vector_, this->vector_ + this->getSize(), 0.0);
}

void TlDenseVector_ImplLapack::sortByGreater() {
  std::sort(this->vector_, this->vector_ + this->getSize(),
            std::greater<double>());
}

double TlDenseVector_ImplLapack::dot(const TlDenseVector_ImplLapack& rhs) const {
  assert(this->getSize() == rhs.getSize());
  const int N = this->getSize();
  const int incX = 1;
  const int incY = 1;

  const double v = ddot_(&N, this->vector_, &incX, rhs.vector_, &incY);
  return v;
}

TlDenseVector_ImplLapack& TlDenseVector_ImplLapack::dotInPlace(
    const TlDenseVector_ImplLapack& rhs) {
  assert(this->getSize() == rhs.getSize());
  std::transform(rhs.vector_, rhs.vector_ + rhs.getSize(), this->vector_,
                 this->vector_, std::multiplies<double>());

  return *this;
}

// ---------------------------------------------------------------------------
// protected
// ---------------------------------------------------------------------------
void TlDenseVector_ImplLapack::initialize(bool clearIfNeeded) {
  const TlDenseVectorObject::size_type size = this->getSize();
  if (size > 0) {
    try {
      this->vector_ = new double[size];
    } catch (std::bad_alloc& ba) {
      this->log_.critical(
          TlUtils::format("bad_alloc caught: %s:  size=%ld", ba.what(), size));
    } catch (...) {
      this->log_.critical("unknown error.");
      throw;
    }
    assert(this->vector_ != NULL);

    if (clearIfNeeded) {
      std::fill(this->vector_, this->vector_ + size, 0.0);
    }
  }
}

// ---------------------------------------------------------------------------
// others
// ---------------------------------------------------------------------------
TlDenseVector_ImplLapack operator+(const TlDenseVector_ImplLapack& rhs1,
                                   const TlDenseVector_ImplLapack& rhs2) {
  TlDenseVector_ImplLapack answer = rhs1;
  answer += rhs2;

  return answer;
}

TlDenseVector_ImplLapack operator-(const TlDenseVector_ImplLapack& rhs1,
                                   const TlDenseVector_ImplLapack& rhs2) {
  TlDenseVector_ImplLapack answer = rhs1;
  answer -= rhs2;

  return answer;
}

TlDenseVector_ImplLapack operator*(const TlDenseVector_ImplLapack& rhs1,
                                   const double rhs2) {
  TlDenseVector_ImplLapack answer = rhs1;

  answer *= rhs2;

  return answer;
}

TlDenseVector_ImplLapack operator*(const double rhs1,
                                   const TlDenseVector_ImplLapack& rhs2) {
  return (rhs2 * rhs1);
}
