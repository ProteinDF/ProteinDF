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

#ifndef TLSHAREDPOINTER_H
#define TLSHAREDPOINTER_H

#include <cassert>
#include <cstddef>

/** リファレンスカウント式Smart Pointer
 *
 *  - newで作成したポインタのみ保持させること。new
 * 以外で生成したポインタは解放できない。
 *  -
 * newで作成したポインタを複数のTlSharedPoinnterに渡さないこと(2重解放になる)。
 *  - 配列を保持させないこと(delete []を行わないため)。
 */
template <typename T>
class TlSharedPointer {
 public:
  TlSharedPointer(T* pObject);
  TlSharedPointer(const TlSharedPointer<T>& psp);

  ~TlSharedPointer();

  TlSharedPointer<T>& operator=(const TlSharedPointer<T>& sp);

  bool operator==(const T* p);
  bool operator!=(const T* p);

  T* operator->() const;
  T& operator*() const;

 private:
  void set(const TlSharedPointer<T>& psp);
  void release();

 private:
  T* m_pObject;
  int* m_pCount;

  const static T* m_pNullObject;  // NULLポインタ
};

////////////////////////////////////////////////////////////////////////////////

template <typename T>
const T* TlSharedPointer<T>::m_pNullObject = NULL;

template <typename T>
TlSharedPointer<T>::TlSharedPointer(T* pObject) {
  this->m_pObject = pObject;

  this->m_pCount = new int();
  *(this->m_pCount) = 1;
}

template <typename T>
TlSharedPointer<T>::TlSharedPointer(const TlSharedPointer<T>& psp) {
  this->m_pObject = NULL;
  this->m_pCount = NULL;

  this->set(psp);
}

template <typename T>
TlSharedPointer<T>::~TlSharedPointer() {
  this->release();
}

template <typename T>
TlSharedPointer<T>& TlSharedPointer<T>::operator=(
    const TlSharedPointer<T>& psp) {
  this->set(psp);
  return (*this);
}

template <typename T>
bool TlSharedPointer<T>::operator==(const T* p) {
  return (this->m_pObject == p);
}

template <typename T>
bool TlSharedPointer<T>::operator!=(const T* p) {
  return !(this->operator==(p));
}

template <typename T>
T* TlSharedPointer<T>::operator->() const {
  return this->m_pObject;
}

template <typename T>
T& TlSharedPointer<T>::operator*() const {
  return *(this->m_pObject);
}

template <typename T>
void TlSharedPointer<T>::set(const TlSharedPointer<T>& psp) {
  if (this != &psp) {
    this->release();

    this->m_pObject = psp.m_pObject;
    this->m_pCount = psp.m_pCount;
    ++(*(this->m_pCount));
  }
}

template <typename T>
void TlSharedPointer<T>::release() {
  if (this->m_pCount != NULL) {
    --(*(this->m_pCount));

    if (*(this->m_pCount) == 0) {
      delete this->m_pObject;
      delete this->m_pCount;

      this->m_pObject = NULL;
      this->m_pCount = NULL;
    }
  }
}

#endif  // TLSHAREDPOINTER_H
