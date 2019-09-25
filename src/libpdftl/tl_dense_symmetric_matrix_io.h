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

#ifndef TL_DENSE_SYMMETRIC_MATRIX_IO_H
#define TL_DENSE_SYMMETRIC_MATRIX_IO_H

#include "tl_dense_matrix_io_object.h"

class TlFileSymmetricMatrix : public TlDenseMatrix_IO_object {
   public:
    explicit TlFileSymmetricMatrix(
        const std::string& filePath, const TlMatrixObject::index_type dim,
        const std::size_t cacheSize = DEFAULT_CACHE_SIZE);
    explicit TlFileSymmetricMatrix(
        const std::string& filePath,
        const std::size_t cacheSize = DEFAULT_CACHE_SIZE);
    virtual ~TlFileSymmetricMatrix();

   public:
    void resize(const TlMatrixObject::index_type dim);

   private:
    // void resize(index_type newRow, index_type newCol){};  // access forbidden

   public:
    virtual TlFileSymmetricMatrix& operator*=(const double coef);
    virtual TlFileSymmetricMatrix& operator/=(const double coef) {
        return this->operator*=(1.0 / coef);
    }

   protected:
    virtual TlMatrixObject::size_type getIndex(
        const TlMatrixObject::index_type row,
        const TlMatrixObject::index_type col) const;
    virtual TlMatrixObject::size_type getNumOfElements() const;
};

#endif  // TL_DENSE_SYMMETRIC_MATRIX_IO_H
