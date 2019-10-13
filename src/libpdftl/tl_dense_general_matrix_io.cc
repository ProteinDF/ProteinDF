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

#include "tl_dense_general_matrix_io.h"
#include <iostream>

TlFileGenericMatrix::TlFileGenericMatrix(const std::string& filePath,
                                         const TlMatrixObject::index_type row,
                                         const TlMatrixObject::index_type col,
                                         const size_t cacheSize)
    : TlDenseMatrix_IO_object(TlMatrixObject::CSFD, filePath, row, col,
                              cacheSize) {
    this->createNewFile();
    this->open();
}

TlFileGenericMatrix::TlFileGenericMatrix(const std::string& filePath,
                                         const size_t cacheSize)
    : TlDenseMatrix_IO_object(TlMatrixObject::CSFD, filePath, cacheSize) {
    this->open();
}

TlFileGenericMatrix::~TlFileGenericMatrix() {}

void TlFileGenericMatrix::resize(const TlMatrixObject::index_type row,
                                 const TlMatrixObject::index_type col) {
    TlDenseMatrix_IO_object::resize<TlFileGenericMatrix>(row, col);
}

TlMatrixObject::size_type TlFileGenericMatrix::getIndex(
    const index_type row, const index_type col) const {
    return TlMatrixObject::getIndex_CSFD(row, col);
}

TlMatrixObject::size_type TlFileGenericMatrix::getNumOfElements() const {
    return TlMatrixObject::getNumOfElements_CSFD();
}
