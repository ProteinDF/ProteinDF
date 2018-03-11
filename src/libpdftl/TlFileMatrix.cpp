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

#include "TlFileMatrix.h"
#include <cassert>
#include <iomanip>
#include <iostream>
#include "TlFile.h"
#include "TlMatrix.h"

TlFileMatrix::TlFileMatrix(const std::string& filePath, const int row,
                           const int col, const size_t cacheSize)
    : TlMatrixObject(RSFD),
      filePath_(filePath),
      numOfRows_(row),
      numOfCols_(col),
      cacheCount_(0),
      cacheSize_(cacheSize) {
  this->initializeCache();
  this->open();
}

// called from sub-class only!
TlFileMatrix::TlFileMatrix(const std::string& filePath, const int row,
                           const int col, bool doOpen, const size_t cacheSize)
    : TlMatrixObject(RSFD),
      filePath_(filePath),
      numOfRows_(row),
      numOfCols_(col),
      cacheCount_(0),
      cacheSize_(cacheSize) {
  this->initializeCache();
  if (doOpen == true) {
    this->open();
  }
}

TlFileMatrix::~TlFileMatrix() { this->finalize(); }

void TlFileMatrix::initializeCache() {
  this->cache_.clear();
  if (this->cacheSize_ < DEFAULT_CACHE_SIZE) {
    this->cacheSize_ = DEFAULT_CACHE_SIZE;
  }
}

void TlFileMatrix::finalize() const {
  std::list<CacheUnit>::const_iterator pEnd = this->cache_.end();
  for (std::list<CacheUnit>::const_iterator p = this->cache_.begin(); p != pEnd;
       ++p) {
    if (p->isUpdate == true) {
      this->writeDisk(*p);
    }
  }

  this->fs_.flush();
  this->fs_.close();
}

void TlFileMatrix::resize(const index_type newRow, const index_type newCol) {
  assert(0 < newRow);
  assert(0 < newCol);

  const index_type oldRow = this->getNumOfRows();
  const index_type oldCol = this->getNumOfCols();

  // rename(move) this object file
  std::string tempFilePath = TlFile::getTempFilePath();
  {
    this->finalize();
    TlFile::rename(this->filePath_, tempFilePath);
  }

  // create new file matrix
  this->numOfRows_ = newRow;
  this->numOfCols_ = newCol;
  this->initializeCache();
  this->open();

  // copy elements
  {
    const TlFileMatrix refMatrix(tempFilePath);
    const index_type maxRow = std::min(oldRow, newRow);
    const index_type maxCol = std::min(oldCol, newCol);
    for (index_type r = 0; r < maxRow; ++r) {
      for (index_type c = 0; c < maxCol; ++c) {
        this->set(r, c, refMatrix.get(r, c));
      }
    }
  }

  // delete temp file
  TlFile::remove(tempFilePath);
}

TlMatrixObject::index_type TlFileMatrix::getNumOfRows() const {
  return this->numOfRows_;
}

TlMatrixObject::index_type TlFileMatrix::getNumOfCols() const {
  return this->numOfCols_;
}

std::size_t TlFileMatrix::getMemSize() const {
  std::size_t answer = sizeof(TlFileMatrix) + this->cacheSize_;

  return answer;
}

void TlFileMatrix::open() {
  if (TlFile::isExistFile(this->filePath_) == false) {
    // create new file
    this->fs_.open(this->filePath_.c_str(),
                   std::ios::binary | std::ios::trunc | std::ios::out);

    // バッファ機能を停止する
    this->fs_ << std::setiosflags(std::ios::unitbuf);

    const int nType = 0;  // means RSFD
    this->fs_.write(reinterpret_cast<const char*>(&nType), sizeof(int));
    this->fs_.write(reinterpret_cast<const char*>(&this->numOfRows_),
                    sizeof(int));
    this->fs_.write(reinterpret_cast<const char*>(&this->numOfCols_),
                    sizeof(int));

    const std::size_t size =
        std::size_t(this->numOfRows_) * std::size_t(this->numOfCols_);
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

bool TlFileMatrix::readHeader() {
  int matrixType = 0;
  index_type row = 0;
  index_type col = 0;
  bool answer = TlMatrix::getHeaderInfo(this->fs_, &matrixType, &row, &col);

  this->numOfRows_ = row;
  this->numOfCols_ = col;

  this->startPos_ = this->fs_.tellg();

  return answer;
}

TlMatrixObject::size_type TlFileMatrix::getIndex(const index_type row,
                                                 const index_type col) const {
  return TlMatrixObject::getIndex_RSFD(row, col);
}

TlMatrixObject::size_type TlFileMatrix::getNumOfElements() const {
  return TlMatrixObject::getNumOfElements_RSFD();
}

void TlFileMatrix::set(const index_type row, const index_type col,
                       const double value) {
  double* const pBuf = this->getCachedData(row, col);
  *pBuf = value;
}

void TlFileMatrix::add(const index_type row, const index_type col,
                       const double value) {
  const double tmp = this->get(row, col) + value;
  double* const pBuf = this->getCachedData(row, col);
  *pBuf = tmp;
}

double TlFileMatrix::get(const index_type row, const index_type col) const {
  return this->getCachedData(row, col);
}

void TlFileMatrix::setRowVector(const index_type row, const TlVector& v) {
  const index_type numOfCols = this->getNumOfCols();
  assert(v.getSize() == numOfCols);

  for (index_type i = 0; i < numOfCols; ++i) {
    this->set(row, i, v[i]);
  }
}

void TlFileMatrix::setColVector(const index_type col, const TlVector& v) {
  const index_type numOfRows = this->getNumOfRows();
  assert(v.getSize() == numOfRows);

  for (index_type i = 0; i < numOfRows; ++i) {
    this->set(i, col, v[i]);
  }
}

TlVector TlFileMatrix::getRowVector(const index_type nRow) const {
  assert((0 <= nRow) && (nRow < this->getNumOfRows()));

  const int nNumOfCols = this->getNumOfCols();
  TlVector answer(nNumOfCols);

  for (int i = 0; i < nNumOfCols; ++i) {
    answer[i] = this->get(nRow, i);
  }

  return answer;
}

TlVector TlFileMatrix::getColumnVector(const index_type nCol) const {
  assert((0 <= nCol) && (nCol < this->getNumOfCols()));

  const int nNumOfRows = this->getNumOfRows();
  TlVector answer(nNumOfRows);

  for (int i = 0; i < nNumOfRows; ++i) {
    answer[i] = this->get(i, nCol);
  }

  return answer;
}

double* TlFileMatrix::getCachedData(const index_type row,
                                    const index_type col) {
  static const unsigned long localIndexBit =
      (static_cast<unsigned long>(1) << CACHE_GROUP_BIT) - 1;

  const size_t index = this->getIndex(row, col);
  const size_t localIndex = (index & localIndexBit);
  this->updateCache(index);

  this->cache_.front().isUpdate = true;
  return &(this->cache_.front().data[localIndex]);
}

double TlFileMatrix::getCachedData(const index_type row,
                                   const index_type col) const {
  static const unsigned long localIndexBit =
      (static_cast<unsigned long>(1) << CACHE_GROUP_BIT) - 1;

  const size_t index = this->getIndex(row, col);
  const size_t localIndex = (index & localIndexBit);
  this->updateCache(index);

  return this->cache_.front().data[localIndex];
}

void TlFileMatrix::updateCache(const size_t index) const {
  // groupの探索
  const size_t group = (static_cast<unsigned long>(index) >> CACHE_GROUP_BIT);

  std::list<CacheUnit>::iterator pCacheUnitEnd = this->cache_.end();
  std::list<CacheUnit>::iterator found =
      std::find_if(this->cache_.begin(), pCacheUnitEnd, CacheUnitComp(group));

  if (found != pCacheUnitEnd) {
    // 先頭に移動
    this->cache_.splice(this->cache_.begin(), this->cache_, found);
  } else {
    // groupの読み込み
    const size_t blockSize =
        sizeof(double) * (static_cast<unsigned long>(1) << CACHE_GROUP_BIT);
    const std::fstream::pos_type pos =
        std::fstream::pos_type(this->startPos_) +
        std::fstream::pos_type(group) * std::fstream::pos_type(blockSize);
    this->fs_.seekg(pos, std::ios_base::beg);

    const size_t groupCount = (1 << CACHE_GROUP_BIT);
    const size_t count =
        std::min(groupCount, (this->getNumOfElements() - group * groupCount));

    CacheUnit cu(group);
    cu.data.resize(count, 0.0);
    this->fs_.read((char*)&(cu.data[0]), sizeof(double) * count);

    this->cache_.push_front(cu);
    ++(this->cacheCount_);
  }

  // バッファサイズから溢れた分を削除
  static const size_t groupCacheSize =
      sizeof(double) * (static_cast<unsigned long>(1) << CACHE_GROUP_BIT);
  while ((this->cacheCount_ * groupCacheSize) > this->cacheSize_) {
    if (this->cache_.back().isUpdate == true) {
      this->writeDisk(this->cache_.back());
    }
    this->cache_.pop_back();
    --(this->cacheCount_);
  }
}

void TlFileMatrix::writeDisk(const TlFileMatrix::CacheUnit& cu) const {
  std::string stateStr = "";
  std::ios_base::iostate state = this->fs_.rdstate();
  if ((state & std::ios_base::eofbit) || (state & std::ios_base::failbit)) {
    this->fs_.clear();
  }

  static const size_t groupCacheSize =
      sizeof(double) * (static_cast<unsigned long>(1) << CACHE_GROUP_BIT);
  const std::fstream::pos_type pos =
      this->startPos_ +
      std::fstream::pos_type(cu.group) * std::fstream::pos_type(groupCacheSize);
  this->fs_.seekp(pos, std::ios_base::beg);

  const size_t bufferSize = sizeof(double) * cu.data.size();

  this->fs_.write((char*)&(cu.data[0]), bufferSize);
}

TlFileMatrix& TlFileMatrix::operator*=(const double coef) {
  // ToDo 高速化
  const index_type numOfRows = this->getNumOfRows();
  const index_type numOfCols = this->getNumOfCols();
  for (index_type r = 0; r < numOfRows; ++r) {
    for (index_type c = 0; c < numOfCols; ++c) {
      const double value = this->get(r, c) * coef;
      this->set(r, c, value);
    }
  }

  return *this;
}

TlMatrix TlFileMatrix::getBlockMatrix(const index_type row,
                                      const index_type col,
                                      const index_type rowDistance,
                                      const index_type colDistance) const {
  assert((0 <= row) && (row < this->getNumOfRows()));
  assert((0 <= col) && (col < this->getNumOfCols()));
  assert(0 < rowDistance);
  assert(0 < colDistance);

  assert(0 <= (row + rowDistance) &&
         (row + rowDistance) <= this->getNumOfRows());
  assert(0 <= (col + colDistance) &&
         (col + colDistance) <= this->getNumOfCols());

  TlMatrix answer(rowDistance, colDistance);
#pragma omp parallel for
  for (index_type dr = 0; dr < rowDistance; ++dr) {
    const index_type r = row + dr;
    for (index_type dc = 0; dc < colDistance; ++dc) {
      const index_type c = col + dc;

      answer.set(dr, dc, this->get(r, c));
    }
  }

  return answer;
}

void TlFileMatrix::setBlockMatrix(const index_type row, const index_type col,
                                  const TlMatrix& matrix) {
  const index_type row_distance = matrix.getNumOfRows();
  const index_type col_distance = matrix.getNumOfCols();

  assert(0 <= row && row < this->getNumOfRows());
  assert(0 <= col && col < this->getNumOfCols());
  assert(0 < (row + row_distance) &&
         (row + row_distance) <= this->getNumOfRows());
  assert(0 < (col + col_distance) &&
         (col + col_distance) <= this->getNumOfCols());

  // this->set()はスレッドセーフではないのでOpenMPでは注意すること
  for (index_type dr = 0; dr < row_distance; ++dr) {
    const index_type r = row + dr;
    for (index_type dc = 0; dc < col_distance; ++dc) {
      const index_type c = col + dc;

      this->set(r, c, matrix.get(dr, dc));
    }
  }
}
