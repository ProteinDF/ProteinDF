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

#include <cmath>
#include <iostream>
#include <limits>

#include "TlMatrixObject.h"

TlMatrixObject::TlMatrixObject(const MatrixType matrixType)
    : matrixType_(matrixType), log_(TlLogging::getInstance()) {}

TlMatrixObject::MatrixType TlMatrixObject::getType() const {
  return this->matrixType_;
}

TlMatrixObject::size_type TlMatrixObject::getNumOfElements_RSFD() const {
  const size_type row = static_cast<size_type>(this->getNumOfRows());
  const size_type col = static_cast<size_type>(this->getNumOfCols());

  const size_type elements = row * col;
  return elements;
}

TlMatrixObject::size_type TlMatrixObject::getNumOfElements_RLHD() const {
  const size_type dim = static_cast<size_type>(this->getNumOfRows());
  assert(dim == static_cast<size_type>(this->getNumOfCols()));

  const size_type elements = dim * (dim + 1) / 2;
  return elements;
}

TlMatrixObject::size_type TlMatrixObject::getIndex_RSFD(index_type row,
                                                        index_type col) const {
  assert(0 <= row);
  assert(row < this->getNumOfRows());
  assert(0 <= col);
  assert(col < this->getNumOfCols());

  const size_type r = static_cast<size_type>(row);
  const size_type c = static_cast<size_type>(col);
  const size_type maxC = static_cast<size_type>(this->getNumOfCols());

  const size_type addr = r * maxC + c;
  return addr;
}

TlMatrixObject::size_type TlMatrixObject::getIndex_CSFD(index_type row,
                                                        index_type col) const {
  assert(0 <= row);
  assert(row < this->getNumOfRows());
  assert(0 <= col);
  assert(col < this->getNumOfCols());

  const size_type r = static_cast<size_type>(row);
  const size_type c = static_cast<size_type>(col);
  const size_type maxR = static_cast<size_type>(this->getNumOfRows());

  const size_type addr = c * maxR + r;
  return addr;
}

TlMatrixObject::size_type TlMatrixObject::getIndex_RLHD(index_type row,
                                                        index_type col) const {
  if (row < col) {
    std::swap(row, col);
  }

  const size_type r = static_cast<size_type>(row);
  const size_type c = static_cast<size_type>(col);

  const size_type addr = r * (r + 1) / 2 + c;
  return addr;
}

void TlMatrixObject::addByList(const index_type* pIndexPairs,
                               const double* pValues, const std::size_t size) {
  for (std::size_t i = 0; i < size; ++i) {
    const index_type globalRow = pIndexPairs[i * 2];
    const index_type globalCol = pIndexPairs[i * 2 + 1];
    const double value = pValues[i];

    this->add(globalRow, globalCol, value);
  }
}

TlVector TlMatrixObject::getRowVector(const int row) const {
  assert((0 <= row) && (row < this->getNumOfRows()));

  const int numOfCols = this->getNumOfCols();
  TlVector answer(numOfCols);

  for (int i = 0; i < numOfCols; ++i) {
    answer[i] = this->get(row, i);
  }

  return answer;
}

TlVector TlMatrixObject::getColVector(const int col) const {
  assert((0 <= col) && (col < this->getNumOfCols()));

  const int numOfRows = this->getNumOfRows();
  TlVector answer(numOfRows);

  for (int i = 0; i < numOfRows; ++i) {
    answer[i] = this->get(i, col);
  }

  return answer;
}

double TlMatrixObject::getMaxAbsoluteElement(int* pOutRow, int* pOutCol) const {
  const int maxRow = this->getNumOfRows();
  const int maxCol = this->getNumOfCols();
  int outRow = -1;
  int outCol = -1;
  double value = 0.0;
  for (int r = 0; r < maxRow; ++r) {
    for (int c = 0; c < maxCol; ++c) {
      const double tmp = std::fabs(this->get(r, c));
      if (value < tmp) {
        outRow = r;
        outCol = c;
        value = tmp;
      }
    }
  }

  if (pOutRow != NULL) {
    *pOutRow = outRow;
  }
  if (pOutCol != NULL) {
    *pOutCol = outCol;
  }

  return value;
}

// TlMatrixObject& TlMatrixObject::operator+=(const TlMatrixObject& rhs)
// {
//     const index_type numOfRows = this->getNumOfRows();
//     const index_type numOfCols = this->getNumOfCols();
//     assert(numOfRows == rhs.getNumOfRows());
//     assert(numOfCols == rhs.getNumOfCols());

//     for (index_type r = 0; r < numOfRows; ++r) {
//         for (index_type c = 0; c < numOfCols; ++c) {
//             this->add(r, c, rhs.get(r, c));
//         }
//     }

//     return *this;
// }
