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

#include "TlFileSymmetricMatrix.h"
#include <iomanip>
#include <iostream>
#include "TlFile.h"
#include "TlSymmetricMatrix.h"

TlFileSymmetricMatrix::TlFileSymmetricMatrix(const std::string& filePath,
                                             const int dim,
                                             const size_t cacheSize)
    : TlFileMatrix(filePath, dim, dim, false, cacheSize) {
  // TlFileMatrix::open() is not called, because base constructer was called
  // using (initialize = false)!
  this->open();
}

TlFileSymmetricMatrix::~TlFileSymmetricMatrix() {
  // this->fs_ is closed by parent class.
}

void TlFileSymmetricMatrix::resize(const index_type newDim) {
  assert(0 < newDim);

  const index_type oldDim = this->getNumOfRows();

  // rename(move) this object file
  std::string tempFilePath = TlFile::getTempFilePath();
  {
    this->finalize();
    TlFile::rename(this->filePath_, tempFilePath);
  }

  // create new file matrix
  this->numOfRows_ = newDim;
  this->numOfCols_ = newDim;
  this->initializeCache();
  this->open();

  // copy elements
  {
    const TlFileSymmetricMatrix refMatrix(tempFilePath);
    const index_type maxDim = std::min(oldDim, newDim);
    for (index_type r = 0; r < maxDim; ++r) {
      for (index_type c = 0; c <= r; ++c) {
        this->set(r, c, refMatrix.get(r, c));
      }
    }
  }

  // delete temp file
  TlFile::remove(tempFilePath);
}

void TlFileSymmetricMatrix::open() {
  if (TlFile::isExistFile(this->filePath_) == false) {
    // create new file
    this->fs_.open(this->filePath_.c_str(),
                   std::ios::binary | std::ios::trunc | std::ios::out);

    // バッファ機能を停止する
    this->fs_ << std::setiosflags(std::ios::unitbuf);

    const int nType = 2;  // means RLHD
    this->fs_.write(reinterpret_cast<const char*>(&nType), sizeof(int));
    this->fs_.write(reinterpret_cast<const char*>(&this->numOfRows_),
                    sizeof(index_type));
    this->fs_.write(reinterpret_cast<const char*>(&this->numOfCols_),
                    sizeof(index_type));

    const std::size_t size = std::size_t(this->numOfRows_) *
                             std::size_t(this->numOfCols_ + 1) / std::size_t(2);
    // size分ファイルにzeroを埋める
    const long unitSize = this->cacheSize_ / sizeof(double);
    ldiv_t ldivt = ldiv(size, unitSize);
    if (ldivt.quot > 0) {
      double* pBuf = new double[unitSize];
      for (int i = 0; i < unitSize; ++i) {
        pBuf[i] = 0.0;
      }

      for (long i = 0; i < ldivt.quot; ++i) {
        this->fs_.write(reinterpret_cast<const char*>(pBuf),
                        sizeof(double) * unitSize);
      }
      delete[] pBuf;
      pBuf = NULL;
    }
    if (ldivt.rem > 0) {
      double* pBuf = new double[ldivt.rem];
      for (int i = 0; i < ldivt.rem; ++i) {
        pBuf[i] = 0.0;
      }

      this->fs_.write(reinterpret_cast<const char*>(pBuf),
                      sizeof(double) * ldivt.rem);
      delete[] pBuf;
      pBuf = NULL;
    }

    // バッファ機能を再開する
    this->fs_ << std::resetiosflags(std::ios::unitbuf);
    this->fs_.flush();

    this->fs_.close();
  }

  this->fs_.open(this->filePath_.c_str(),
                 std::ios::binary | std::ios::in | std::ios::out);
  this->readHeader();
}

bool TlFileSymmetricMatrix::readHeader() {
  int matrixType = 0;
  index_type dim = 0;
  bool answer = TlSymmetricMatrix::getHeaderInfo(this->fs_, &matrixType, &dim);

  this->numOfRows_ = dim;
  this->numOfCols_ = dim;

  this->startPos_ = this->fs_.tellg();

  return answer;
}

TlMatrixObject::size_type TlFileSymmetricMatrix::getIndex(
    index_type row, index_type col) const {
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
