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

#ifndef TL_DENSE_GENERAL_MATRIX_IO_H
#define TL_DENSE_GENERAL_MATRIX_IO_H

#include "tl_dense_matrix_io_object.h"

// #define CACHE_GROUP_BIT (20) //  1048576 (2^20)個のローカルインデックス
// (=8MB)
// #define CACHE_GROUP_BIT (24) // 16777216 (2^24)個のローカルインデックス
// (=128MB for double)
// #define DEFAULT_CACHE_SIZE (1 * 1024 * 1024 * 1024) // 1 GB
#define DEFAULT_CACHE_SIZE (2147483648)  // 2 GB = 2 * 1073741824

class TlFileGenericMatrix : public TlDenseMatrix_IO_object {
 public:
  explicit TlFileGenericMatrix(const std::string& filePath,
                               TlMatrixObject::index_type row,
                               TlMatrixObject::index_type col,
                               std::size_t cacheSize = DEFAULT_CACHE_SIZE);
  explicit TlFileGenericMatrix(const std::string& filePath,
                               std::size_t cacheSize = DEFAULT_CACHE_SIZE);
  virtual ~TlFileGenericMatrix();

 public:
  void resize(const TlMatrixObject::index_type row,
              const TlMatrixObject::index_type col);

 protected:
  virtual TlMatrixObject::size_type getIndex(
      const TlMatrixObject::index_type row,
      const TlMatrixObject::index_type col) const;
  virtual TlMatrixObject::size_type getNumOfElements() const;
};

#endif  // TL_DENSE_GENERAL_MATRIX_IO_H
