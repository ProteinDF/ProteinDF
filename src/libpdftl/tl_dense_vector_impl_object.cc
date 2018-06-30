#include <cassert>
#include <cmath>
#include "tl_dense_vector_impl_object.h"
#include "TlUtils.h"
#include "tl_vector_utils.h"

// ---------------------------------------------------------------------------
// constructor & destructor
// ---------------------------------------------------------------------------
TlDenseVector_ImplObject::TlDenseVector_ImplObject()
    : log_(TlLogging::getInstance()) {}

TlDenseVector_ImplObject::~TlDenseVector_ImplObject() {}

// ---------------------------------------------------------------------------
// operations
// ---------------------------------------------------------------------------
double TlDenseVector_ImplObject::sum() const {
    double answer = 0.0;

    const TlDenseVectorObject::size_type dim = this->getSize();
    for (TlDenseVectorObject::size_type i = 0; i < dim; ++i) {
        answer += this->get(i);
    }

    return answer;
}

double TlDenseVector_ImplObject::norm() const {
  return std::sqrt(this->norm2());
}

double TlDenseVector_ImplObject::norm2() const {
  double answer = 0.0;
  const TlDenseVectorObject::size_type size = this->getSize();
#pragma omp parallel for reduction(+ : answer)
  for (TlDenseVectorObject::size_type i = 0; i < size; ++i) {
    const double v = this->get(i);
    answer += v * v;
  }

  return answer;
}

TlDenseVectorObject::index_type TlDenseVector_ImplObject::argmax(
    const TlDenseVectorObject::index_type& begin,
    const TlDenseVectorObject::index_type& end) const {
    assert((0 <= begin) && (begin < this->getSize()));
    assert((0 <= end) && (end <= this->getSize()));
  TlDenseVectorObject::index_type answer = begin;
  double maxValue = this->get(begin);

  for (TlDenseVectorObject::index_type i = begin; i < end; ++i) {
    const double value = this->get(i);
    if (maxValue < value) {
      answer = i;
      maxValue = value;
    }
  }

  return answer;
}
