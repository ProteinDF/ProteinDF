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

#ifndef TLSIMPLEVECTOR_H
#define TLSIMPLEVECTOR_H

#include <algorithm>
#include <cassert>

template <typename T>
class TlSimpleVector {
 public:
  explicit TlSimpleVector(size_t nSize = 0, const T& value = T())
      : m_nSize(nSize), m_pData(NULL) {
    this->m_pData = new T[m_nSize];
    //     const size_t nEnd = m_nSize;
    //     for (size_t i = 0; i < nEnd; ++i){
    //       this->m_pData[i] = value;
    //     }
    std::fill(this->m_pData, (this->m_pData + this->m_nSize), value);
  }

  TlSimpleVector(const TlSimpleVector<T>& rhs)
      : m_nSize(rhs.m_nSize), m_pData(NULL) {
    this->m_pData = new T[m_nSize];
    //     const size_t nEnd = m_nSize;
    //     for (size_t i = 0; i < nEnd; ++i){
    //       pData[i] = rhs.m_pData[i];
    //     }
    std::copy(rhs.m_pData, (rhs.m_pData + rhs.m_nSize), this->m_pData);
  }

  ~TlSimpleVector() {
    delete[] this->m_pData;
    this->m_pData = NULL;
  }

  TlSimpleVector& operator=(const TlSimpleVector<T>& rhs) {
    if (this != &rhs) {
      if (this->m_pData != NULL) {
        delete[] this->m_pData;
      }

      m_nSize = rhs.m_nSize;
      this->m_pData = new T[m_nSize];
      //       const size_t nEnd = m_nSize;
      //       for (size_t i = 0; i < nEnd; ++i){
      //    pData[i] = rhs.pData[i];
      //       }
      std::copy(rhs.m_pData, (rhs.m_pData + rhs.m_nSize), this->m_pData);
    }

    return *this;
  }

  TlSimpleVector& operator+=(const TlSimpleVector<T>& rhs) {
    assert(m_nSize == rhs.m_nSize);

    const size_t nEnd = m_nSize;
    for (size_t i = 0; i < nEnd; ++i) {
      this->m_pData[i] += rhs.pData[i];
    }

    return *this;
  }

  T& operator[](size_t index) {
    assert(index < m_nSize);
    return this->m_pData[index];
  }

  const T& operator[](size_t index) const {
    assert(index < m_nSize);
    return this->m_pData[index];
  }

  size_t size() const { return this->m_nSize; }

  T max() const {
    assert(m_nSize > 0);
    T answer = this->m_pData[0];
    const size_t nEnd = m_nSize;
    for (size_t i = 1; i < nEnd; ++i) {
      answer = std::max(answer, this->m_pData[i]);
    }

    return answer;
  }

 public:
  friend TlSimpleVector<T>& operator+(const TlSimpleVector<T>& x,
                                      TlSimpleVector<T>& y) {
    TlSimpleVector<T> answer = x;
    answer += y;

    return answer;
  }

 private:
  size_t m_nSize;
  T* m_pData;
};

#endif  // TLSIMPLEVECTOR_H
