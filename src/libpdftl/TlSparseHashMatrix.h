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

#ifndef TLSPARSEHASHMATRIX_H
#define TLSPARSEHASHMATRIX_H

#include <cassert>
#include <ext/hash_map> // -> std::unordered_map (C++11)
//#include "TlHashMap.h"

class TlSparseHashMatrix {
public:
  typedef __gnu_cxx::hash_map<unsigned long, double> ContainerType;
  typedef ContainerType::const_iterator const_iterator;
  typedef ContainerType::iterator iterator;

public:	 
  TlSparseHashMatrix(int numOfRows =1, int numOfCols =1);
  TlSparseHashMatrix(const TlSparseHashMatrix& rhs);
  ~TlSparseHashMatrix();

public:
  /// 最初のイテレータを返す
  const_iterator begin() const {
    return this->container.begin();
  }

  /// 最初のイテレータを返す
  iterator begin() {
    return this->container.begin();
  }

  /// 最後のイテレータを返す
  const_iterator end() const {
    return this->container.end();
  }

  /// 最後のイテレータを返す
  iterator end() {
    return this->container.end();
  }

public:
  TlSparseHashMatrix& operator=(const TlSparseHashMatrix& rhs);

  int getNumOfRows() const;
  int getNumOfCols() const;
  int getSize() const;
  void resize(int row, int col);

  double& operator()(int row, int col);

public:
  unsigned long index(const int row, const int col) const {
    const unsigned int r = static_cast<unsigned int>(row);
    const unsigned int c = static_cast<unsigned int>(col);

    return ((static_cast<unsigned long>(r) << INT_BITS) + c);
  }

  void index(const unsigned long i, int* pRow, int* pCol) const {
    assert(pRow != NULL);
    assert(pCol != NULL);
    
    const unsigned int r = static_cast<unsigned int>(i >> INT_BITS);
    const unsigned int c = static_cast<unsigned int>(i & TlSparseHashMatrix::MAX_INT);
    
    *pRow = static_cast<int>(r);
    *pCol = static_cast<int>(c);
  }

protected:
  int numOfRows;
  int numOfCols;

  //TlHashMap container;
  //__gnu_cxx::hash_map<unsigned long, int, std::tr1::hash<unsigned long> > container;
  ContainerType container;

  static const int INT_BITS;
  static const int MAX_INT;

  friend class TlCommunicate;
};

#endif // TLSPARSEHASHMATRIX_H
