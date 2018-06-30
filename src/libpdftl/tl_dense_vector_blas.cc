// Copyright (C) 2002-2014 The ProteinDF project
// see also AUTHORS and README.
//
// This file is part of ProteinDF.
//
// ProteinDF is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// ProteinDF is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with ProteinDF.  If not, see <http://www.gnu.org/licenses/>.

#ifdef HAVE_CONFIG_H
#include "config.h"  // this file created by autotools
#endif               // HAVE_CONFIG_H

#include <algorithm>
#include <cmath>
#include <cstring>
#include <fstream>
#include <functional>
#include <iostream>
#include <numeric>
#include "tl_dense_vector_blas.h"
#include "tl_vector_utils.h"

#ifdef HAVE_HDF5
#include "TlHdf5Utils.h"
#endif  // HAVE_HDF5

TlVector_BLAS::TlVector_BLAS(TlVectorAbstract::index_type size)
    : size_(size), data_(NULL) {
  this->initialize();
}

TlVector_BLAS::TlVector_BLAS(const double* p, TlVectorAbstract::size_type size)
    : size_(size), data_(NULL) {
  this->initialize(false);
  std::copy(p, p + size, this->data_);
}

TlVector_BLAS::TlVector_BLAS(const TlVector_BLAS& rhs)
    : size_(rhs.size_), data_(NULL) {
  this->initialize(false);
  std::copy(rhs.data_, rhs.data_ + rhs.getSize(), this->data_);
}

TlVector_BLAS::TlVector_BLAS(const std::vector<double>& rhs)
    : size_(rhs.size()), data_(NULL) {
  this->initialize(false);
  std::copy(rhs.begin(), rhs.end(), this->data_);
}

TlVector_BLAS::~TlVector_BLAS() { this->destroy(); }

void TlVector_BLAS::initialize(bool isZeroClear) {
  const size_type size = this->getSize();
  this->data_ = new double[size];

  if (isZeroClear == true) {
    this->zeroClear();
  }
}

void TlVector_BLAS::destroy() {
  this->size_ = 0;
  delete[] this->data_;
  this->data_ = NULL;
}

void TlVector_BLAS::resize(const TlVectorAbstract::index_type newSize) {
  double* pOld = this->data_;
  const TlVectorAbstract::index_type oldSize = this->size_;

  this->size_ = newSize;
  this->initialize(true);

  const size_type fillSize = std::min(newSize, oldSize);
  std::copy(pOld, pOld + fillSize, this->data_);

  delete[] pOld;
  pOld = NULL;
}

void TlVector_BLAS::push_back(const double value) {
  const size_type size = this->getSize();
  this->resize(size + 1);
  this->data_[size] = value;
}

double TlVector_BLAS::getMaxAbsoluteElement() const {
  double answer = 0.0;
  const size_type size = this->getSize();
  for (size_type i = 0; i < size; ++i) {
    answer = std::max(answer, std::fabs(this->get(i)));
  }

  return answer;
}

TlVector_BLAS& TlVector_BLAS::dotInPlace(const TlVector_BLAS& rhs) {
  assert(this->getSize() == rhs.getSize());
  std::transform(rhs.data_, rhs.data_ + rhs.getSize(), this->data_, this->data_,
                 std::multiplies<double>());

  return *this;
}

double TlVector_BLAS::sum() const {
  return std::accumulate(this->data_, this->data_ + this->getSize(), 0.0);
}

void TlVector_BLAS::sortByGrater() {
  // std::cout << "call TlVector_BLAS::sortByGraterEqual()" << std::endl;
  std::sort(this->data_, this->data_ + this->getSize(), std::greater<double>());
}

std::vector<TlVectorAbstract::size_type>::const_iterator TlVector_BLAS::argmax(
    const std::vector<TlVectorAbstract::size_type>::const_iterator& begin,
    const std::vector<TlVectorAbstract::size_type>::const_iterator& end) const {
  std::vector<TlVectorAbstract::size_type>::const_iterator answer = begin;
  double maxValue = this->get(*begin);
  for (std::vector<size_type>::const_iterator it = begin; it != end; ++it) {
    const size_type index = *it;
    const double value = this->get(index);
    if (maxValue < value) {
      answer = it;
      maxValue = value;
    }
  }

  return answer;
}

double TlVector_BLAS::norm2() const {
  double answer = 0.0;
  const std::size_t size = this->getSize();
#pragma omp parallel for reduction(+ : answer)
  for (std::size_t i = 0; i < size; ++i) {
    const double v = this->get(i);
    answer += v * v;
  }

  return answer;
}

double TlVector_BLAS::norm() const { return std::sqrt(this->norm2()); }

TlVector_BLAS& TlVector_BLAS::operator=(const TlVector_BLAS& rhs) {
  if (this != &rhs) {
    // this->clear();
    delete[] this->data_;
    this->data_ = NULL;

    this->size_ = rhs.size_;
    this->initialize(false);

    std::copy(rhs.data_, rhs.data_ + this->getSize(), this->data_);
  }

  return (*this);
}

TlVector_BLAS& TlVector_BLAS::operator+=(const TlVector_BLAS& rhs) {
  assert(this->getSize() == rhs.getSize());

  const size_type size = this->getSize();
#pragma omp parallel for
  for (size_type i = 0; i < size; ++i) {
    this->data_[i] += rhs.get(i);
  }

  return (*this);
}

TlVector_BLAS& TlVector_BLAS::operator+=(const TlSparseVector& rhs) {
  assert(this->getSize() == rhs.getSize());

  TlSparseVector::const_iterator itEnd = rhs.end();
  for (TlSparseVector::const_iterator it = rhs.begin(); it != itEnd; ++it) {
    this->data_[it->first] += it->second;
  }

  return (*this);
}

TlVector_BLAS& TlVector_BLAS::operator-=(const TlVector_BLAS& rhs) {
  assert(this->getSize() == rhs.getSize());

  const size_type size = this->getSize();
#pragma omp parallel for
  for (size_type i = 0; i < size; ++i) {
    this->data_[i] -= rhs.get(i);
  }

  return (*this);
}

TlVector_BLAS& TlVector_BLAS::operator*=(const double& rhs) {
  const size_type size = this->getSize();
#pragma omp parallel for
  for (size_type i = 0; i < size; ++i) {
    this->data_[i] *= rhs;
  }

  return (*this);
}

TlVector_BLAS& TlVector_BLAS::operator/=(const double& rhs) {
  return (this->operator*=(rhs));
}

TlVector_BLAS operator+(const TlVector_BLAS& X, const TlVector_BLAS& Y) {
  assert(X.getSize() == Y.getSize());

  TlVector_BLAS answer = X;
  answer += Y;
  return answer;
}

TlVector_BLAS operator-(const TlVector_BLAS& X, const TlVector_BLAS& Y) {
  assert(X.getSize() == Y.getSize());

  TlVector_BLAS answer = X;
  answer -= Y;
  return answer;
}

TlVector_BLAS operator*(const TlVector_BLAS& X, const double& Y) {
  TlVector_BLAS answer = X;
  answer *= Y;

  return answer;
}

double operator*(const TlVector_BLAS& X, const TlVector_BLAS& Y) {
  assert(X.getSize() == Y.getSize());

  return std::inner_product(X.data_, X.data_ + X.getSize(), Y.data_, 0.0);
}

bool TlVector_BLAS::load(const std::string& filePath) {
  bool answer = false;

  TlVectorAbstract::index_type length = 0;
  const TlVectorUtils::FileSize headerSize =
      TlVectorUtils::getHeaderInfo(filePath, &length);
  if (headerSize > 0) {
    this->resize(length);
    TlVectorUtils::load(filePath, this->data_, length);
  } else {
    std::cerr << TlUtils::format("load failed.: %d@%s", __FILE__, __LINE__)
              << std::endl;
    answer = false;
  }

  return answer;
}

bool TlVector_BLAS::loadText(const std::string& filePath) {
  bool answer = false;

  std::ifstream ifs;
  ifs.open(filePath.c_str(), std::ios::in);
  if (ifs.fail()) {
#ifdef DEBUG
    std::cerr << "[error] TlVector_BLAS::loadText(): could not open file. "
              << filePath << std::endl;
#endif  // DEBUG
    return false;
  }

  // load contents
  std::string line = "";
  std::getline(ifs, line);  // read 1st line

  if (line == "TEXT") {
    std::string tmp = "";
    ifs >> tmp;
    const int size = std::atoi(tmp.c_str());
    ifs >> tmp;  // equal to 'size'
    ifs >> tmp;  // equal to '0'
    this->resize(size);
    for (int i = 0; i < size; ++i) {
      ifs >> tmp;
      const double v = std::atof(tmp.c_str());
      (*this)[i] = v;
    }

    answer = true;
  }

  ifs.close();

  return answer;
}

bool TlVector_BLAS::save(const std::string& filePath) const {
  bool answer = true;

  std::ofstream ofs;
  ofs.open(filePath.c_str(), std::ofstream::out | std::ofstream::binary);

  const TlVectorAbstract::size_type nSize = this->getSize();
  ofs.write(reinterpret_cast<const char*>(&nSize),
            sizeof(TlVectorAbstract::size_type));

  const char* p =
      reinterpret_cast<const char*>(const_cast<TlVector_BLAS*>(this)->data_);
  ofs.write(p, sizeof(double) * nSize);

  ofs.close();

  return answer;
}

#ifdef HAVE_HDF5
bool TlVector_BLAS::saveHdf5(const std::string& filepath,
                             const std::string& h5path) const {
  TlHdf5Utils h5(filepath);

  const index_type size = this->getSize();
  h5.write(h5path, this->data_, size);
  h5.setAttr(h5path, "size", size);

  return true;
}

bool TlVector_BLAS::loadHdf5(const std::string& filepath,
                             const std::string& h5path) {
  TlHdf5Utils h5(filepath);

  index_type size;
  h5.getAttr(h5path, "size", &size);

  this->destroy();
  this->size_ = size;
  this->initialize();
  h5.get(h5path, reinterpret_cast<double*>(this->data_), size);

  return true;
}

#endif  // HAVE_HDF5

void TlVector_BLAS::outputText(std::ostream& os) const {
  const size_type nSize = this->getSize();

  os << "TEXT\n";
  os << nSize << "\n";
  os << nSize << "\n";
  os << "0\n";

  for (size_type j = 0; j < nSize; j += 10) {
    for (size_type i = j; ((i < j + 10) && (i < nSize)); ++i) {
      os << TlUtils::format("  %10.4lf", this->get(i));
    }
    os << std::endl;
  }
  os << std::endl;
}
