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

#include <cassert>
#include <limits>

#include "TlSparseHashMatrix.h"

const int TlSparseHashMatrix::INT_BITS = sizeof(int) * 8;
const int TlSparseHashMatrix::MAX_INT =
    std::numeric_limits<unsigned int>::max();

TlSparseHashMatrix::TlSparseHashMatrix(int rows, int cols)
    : numOfRows(rows), numOfCols(cols) {}

TlSparseHashMatrix::TlSparseHashMatrix(const TlSparseHashMatrix& rhs)
    : numOfRows(rhs.numOfRows),
      numOfCols(rhs.numOfCols),
      container(rhs.container) {}

TlSparseHashMatrix::~TlSparseHashMatrix() {}

int TlSparseHashMatrix::getNumOfRows() const { return this->numOfRows; }

int TlSparseHashMatrix::getNumOfCols() const { return this->numOfCols; }

void TlSparseHashMatrix::resize(const int row, const int col) {
  assert(0 < row);
  assert(0 < col);

  if ((row > this->getNumOfRows()) || (col > this->getNumOfCols())) {
    ContainerType::iterator pEnd = this->container.end();
    for (ContainerType::iterator p = this->container.begin(); p != pEnd; ++p) {
      int r;
      int c;
      this->index(p->first, &r, &c);

      if ((r >= row) || (c > col)) {
        this->container.erase(p);
      }
    }
  }

  this->numOfRows = row;
  this->numOfCols = col;
}

int TlSparseHashMatrix::getSize() const { return this->container.size(); }

TlSparseHashMatrix& TlSparseHashMatrix::operator=(
    const TlSparseHashMatrix& rhs) {
  this->numOfRows = rhs.numOfRows;
  this->numOfCols = rhs.numOfCols;
  this->container.clear();
  this->container = rhs.container;

  return (*this);
}

double& TlSparseHashMatrix::operator()(const int row, const int col) {
  const unsigned long i = this->index(row, col);

  return (this->container[i]);
}
