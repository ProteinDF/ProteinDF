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

#ifndef TL_MATRIX_FILE_OBJECT_H
#define TL_MATRIX_FILE_OBJECT_H

#include <algorithm>
#include <cassert>
#include <fstream>
#include <list>
#include <string>
#include <vector>
#include "TlFile.h"
#include "tl_dense_general_matrix_blas_old.h"
#include "tl_matrix_object.h"
#include "tl_dense_vector_blas.h"

// #define CACHE_GROUP_BIT (20) //  1048576 (2^20)個のローカルインデックス
// (=8MB)
#define CACHE_GROUP_BIT \
  (24)  // 16777216 (2^24)個のローカルインデックス (=128MB for double)
// #define DEFAULT_CACHE_SIZE (1 * 1024 * 1024 * 1024) // 1 GB
#define DEFAULT_CACHE_SIZE (2147483648)  // 2 GB = 2 * 1073741824

class TlDenseMatrix_IO_object : public TlMatrixObject {
  friend class TlCommunicate;

 protected:
  // for cache
  struct CacheUnit {
   public:
    CacheUnit(unsigned long g) : group(g), isUpdate(false) {}

   public:
    unsigned long group;
    bool isUpdate;
    std::vector<double> data;

    friend bool operator==(const CacheUnit& x, const CacheUnit& y) {
      return (x.group == y.group);
    }
  };

  struct CacheUnitComp : public std::unary_function<const CacheUnit&, bool> {
   public:
    CacheUnitComp(unsigned long group) : group_(group) {}

    bool operator()(const CacheUnit& cu) { return (group_ == cu.group); }

   private:
    unsigned long group_;
  };

 public:
  // create new file
  explicit TlDenseMatrix_IO_object(const TlMatrixObject::MatrixType matrixType,
                                   const std::string& filePath, index_type row,
                                   index_type col,
                                   std::size_t cacheSize = DEFAULT_CACHE_SIZE);
  // load file
  explicit TlDenseMatrix_IO_object(const TlMatrixObject::MatrixType matrixType,
                                   const std::string& filePath,
                                   std::size_t cacheSize = DEFAULT_CACHE_SIZE);
  virtual ~TlDenseMatrix_IO_object();

 protected:
  /// subclass用コンストラクタ
  ///
  /// ファイルをsubclassで作成する場合には、`doOpen = false` にすること。
  // TlDenseMatrix_IO_object(const std::string& filePath, index_type row,
  //                    index_type col, bool doOpen, std::size_t cacheSize);

  void initializeCache();

  /// cacheに残ったデータの書き込みなど、終了処理を行う
  void finalize() const;

  // public:
  // index_type getNumOfRows() const;
  // index_type getNumOfCols() const;

  virtual std::size_t getMemSize() const;

 public:
  virtual void set(index_type row, index_type col, double value);
  virtual double get(index_type row, index_type col) const;
  virtual void add(index_type row, index_type col, double value);

  virtual TlDenseMatrix_IO_object& operator*=(double coef);
  virtual TlDenseMatrix_IO_object& operator/=(double coef) {
    return this->operator*=(1.0 / coef);
  }

  /// 指定した行の要素から構成されるベクトルを返す
  ///
  /// @param[in] nRow 指定する行
  virtual TlVector_BLAS getRowVector(index_type row) const;

  /// 指定した列の要素から構成されるベクトルを返す
  ///
  /// @param[in] nCol 指定する列
  virtual TlVector_BLAS getColumnVector(index_type col) const;

  virtual void setRowVector(const index_type row, const TlVector_BLAS& v);
  virtual void setColVector(const index_type col, const TlVector_BLAS& v);

  /// ブロック行列を返す
  ///
  /// @param[in] row 始点となる行
  /// @param[in] col 始点となる列
  /// @param[in] row_distance 取得する行数
  /// @param[in] col_distance 取得する列数
  /// @return row_distance × col_distance 次元のブロック行列
  virtual TlDenseGeneralMatrix_BLAS_old getBlockMatrix(
      index_type row, index_type col, index_type rowDistance,
      index_type colDistance) const;

  /// 行列要素を指定された位置に上書きする
  ///
  /// @param[in] row 始点となる行
  /// @param[in] col 始点となる列
  /// @param[in] matrix 行列要素
  virtual void setBlockMatrix(index_type row, index_type col,
                              const TlDenseGeneralMatrix_BLAS_old& matrix);

 protected:
  void createNewFile();
  void open();
  virtual bool readHeader();

  virtual TlMatrixObject::size_type getIndex(
      const TlMatrixObject::index_type row,
      const TlMatrixObject::index_type col) const = 0;
  virtual TlMatrixObject::size_type getNumOfElements() const = 0;

  double* getCachedData(const int row, const int col);
  double getCachedData(const int row, const int col) const;
  void updateCache(size_t index) const;

 protected:
  template <typename FileMatrixObject>
  void resize(const TlMatrixObject::index_type row,
              const TlMatrixObject::index_type col);

 protected:
  void writeDisk(const CacheUnit& cu) const;

 protected:
  virtual bool load(const std::string& path) { return false; }

  virtual bool save(const std::string& path) const { return false; }

 protected:
  std::string filePath_;
  mutable std::fstream fs_;
  std::fstream::pos_type startPos_;

  mutable std::list<CacheUnit> cache_;
  mutable size_t cacheCount_;  // == cache_.size()
  size_t cacheSize_;
};

template <typename FileMatrixObject>
void TlDenseMatrix_IO_object::resize(const TlMatrixObject::index_type newRow,
                                     const TlMatrixObject::index_type newCol) {
  assert(0 < newRow);
  assert(0 < newCol);

  const index_type oldRow = this->getNumOfRows();
  const index_type oldCol = this->getNumOfCols();

  // finalize this object
  this->finalize();

  // copy
  const std::string tempFilePath = TlFile::getTempFilePath();
  TlFile::copy(this->filePath_, tempFilePath);
  TlFile::remove(this->filePath_);
  assert(TlFile::isExistFile(this->filePath_) == false);

  // create new file matrix
  this->row_ = newRow;
  this->col_ = newCol;
  this->initializeCache();
  this->createNewFile();
  this->open();
  assert(this->getNumOfRows() == newRow);
  assert(this->getNumOfCols() == newCol);

  // copy elements
  {
    const FileMatrixObject refMatrix(tempFilePath);
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
#endif  // TL_MATRIX_FILE_OBJECT_H
