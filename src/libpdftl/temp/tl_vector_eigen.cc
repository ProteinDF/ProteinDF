#include <iostream>
#include "TlUtils.h"
#include "tl_dense_vector_eigen.h"
#include "tl_vector_utils.h"

TlEigenVector::TlEigenVector(TlVectorAbstract::index_type size)
    : vector_(VectorType::Zero(size)) {}

TlEigenVector::TlEigenVector(const TlEigenVector& rhs) : vector_(rhs.vector_) {}

TlEigenVector::TlEigenVector(const double* p, size_type size)
    : vector_(MapTypeConst(p, size)) {}

TlEigenVector::TlEigenVector(const std::vector<double>& rhs)
    : vector_(MapTypeConst(rhs.data(), rhs.size())) {}

TlEigenVector::TlEigenVector(const TlEigenVector::VectorType& rhs)
    : vector_(rhs) {}

TlEigenVector::~TlEigenVector() {}

TlVectorAbstract::size_type TlEigenVector::getSize() const {
  return this->vector_.rows();
}

void TlEigenVector::resize(TlVectorAbstract::index_type newSize) {
  this->vector_.conservativeResizeLike(VectorType::Zero(newSize, 1));
}

double TlEigenVector::get(const TlVectorAbstract::index_type i) const {
  return this->vector_.coeff(i);
}

void TlEigenVector::set(const TlEigenVector::index_type i, const double value) {
  this->vector_.coeffRef(i) = value;
}

void TlEigenVector::add(const TlEigenVector::index_type i, const double value) {
  this->vector_.coeffRef(i) += value;
}

double TlEigenVector::operator[](const TlVectorAbstract::index_type i) const {
  return this->get(i);
}

double& TlEigenVector::operator[](const TlEigenVector::index_type i) {
  return this->vector_.coeffRef(i);
}

// -----------------------------------------------------------------------------
// I/O
// -----------------------------------------------------------------------------
bool TlEigenVector::load(const std::string& filePath) {
  bool answer = false;

  TlVectorAbstract::index_type length = 0;
  const TlVectorUtils::FileSize headerSize =
      TlVectorUtils::getHeaderInfo(filePath, &length);
  if (headerSize > 0) {
    this->resize(length);
    const TlVectorAbstract::size_type copied =
        TlVectorUtils::load(filePath, this->vector_.data(), length);
    answer = (copied == length);
  } else {
    std::cerr << TlUtils::format("load failed.: %d@%s", __FILE__, __LINE__)
              << std::endl;
    answer = false;
  }

  return answer;
}

bool TlEigenVector::save(const std::string& filePath) const {
  return TlVectorUtils::save(filePath, this->getSize(), this->vector_.data(),
                             this->getSize());
}
