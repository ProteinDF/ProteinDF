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

#include "tl_sparse_vector.h"
#include <cassert>

TlSparseVector::TlSparseVector(TlSparseVector::size_type size) {
  this->resize(size);
}

TlSparseVector::~TlSparseVector() {}

TlSparseVector::size_type TlSparseVector::getSize() const {
  return this->size_;
}

void TlSparseVector::resize(TlSparseVector::size_type newSize) {
  assert(newSize > 0);
  this->size_ = newSize;
}

double TlSparseVector::get(const TlSparseVector::size_type index) const {
  double answer = 0.0;
  SparseVectorData::const_iterator it = this->data_.find(index);
  if (it != this->data_.end()) {
    answer = it->second;
  }
  return answer;
}

void TlSparseVector::set(const TlSparseVector::size_type index,
                         const double value) {
  SparseVectorData::iterator it = this->data_.find(index);
  if (it != this->data_.end()) {
    this->data_[index] = value;
  } else {
    it->second = value;
  }
}

void TlSparseVector::add(const TlSparseVector::size_type index,
                         const double value) {
  // this->data_[index] += value;
  SparseVectorData::iterator it = this->data_.find(index);
  if (it != this->data_.end()) {
    it->second += value;
  } else {
    this->data_[index] = value;
  }
}

double TlSparseVector::operator[](const TlSparseVector::size_type index) const {
  return this->get(index);
}

double& TlSparseVector::operator[](const TlSparseVector::size_type index) {
  SparseVectorData::iterator it = this->data_.find(index);
  if (it != this->data_.end()) {
    return it->second;
  } else {
    this->data_[index] = 0.0;
    return this->data_[index];
  }
}

TlSparseVector& TlSparseVector::operator+=(const TlSparseVector& rhs) {
  assert(this->getSize() == rhs.getSize());

  SparseVectorData::const_iterator itEnd = rhs.end();
  for (SparseVectorData::const_iterator it = rhs.begin(); it != itEnd; ++it) {
    KeyType key = it->first;
    double value = it->second;
    SparseVectorData::iterator p = this->data_.find(key);
    if (p != this->data_.end()) {
      p->second += value;
    } else {
      this->data_.insert(std::pair<KeyType, double>(key, value));
    }
  }

  return *this;
}
