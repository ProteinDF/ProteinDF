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

#include "tl_dense_symmetric_matrix_io.h"
#include <iomanip>
#include <iostream>
#include "tl_matrix_utils.h"

TlFileSymmetricMatrix::TlFileSymmetricMatrix(
    const std::string& filePath, const TlMatrixObject::index_type dim,
    const std::size_t cacheSize)
    : TlDenseMatrix_IO_object(TlMatrixObject::RLHD, filePath, dim, dim,
                              cacheSize) {
  this->createNewFile();
  this->open();
}

TlFileSymmetricMatrix::TlFileSymmetricMatrix(const std::string& filePath,
                                             const std::size_t cacheSize)
    : TlDenseMatrix_IO_object(TlMatrixObject::RLHD, filePath, cacheSize) {
  this->open();
}

TlFileSymmetricMatrix::~TlFileSymmetricMatrix() {}

void TlFileSymmetricMatrix::resize(const index_type newDim) {
  TlDenseMatrix_IO_object::resize<TlFileSymmetricMatrix>(newDim, newDim);
}

TlMatrixObject::size_type TlFileSymmetricMatrix::getIndex(
    const TlMatrixObject::index_type row,
    const TlMatrixObject::index_type col) const {
  return TlMatrixObject::getIndex_RLHD(row, col);
}

TlMatrixObject::size_type TlFileSymmetricMatrix::getNumOfElements() const {
  return TlMatrixObject::getNumOfElements_RLHD();
}

TlFileSymmetricMatrix& TlFileSymmetricMatrix::operator*=(const double coef) {
  // ToDo: 高速化
  const size_type numOfRows = this->getNumOfRows();
  const size_type numOfCols = this->getNumOfCols();
  for (size_type r = 0; r < numOfRows; ++r) {
    for (size_type c = 0; c < numOfCols; ++c) {
      const double value = this->get(r, c) * coef;
      this->set(r, c, value);
    }
  }

  return *this;
}
