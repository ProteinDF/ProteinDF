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

#ifndef TL_MATRIX_OBJECT_H
#define TL_MATRIX_OBJECT_H

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif  // HAVE_CONFIG_H

#include "TlLogging.h"
#include "tl_dense_vector_lapack.h"

/// MAX_INDEX_BITSが指定するビット数まで行列の次数を確保できる。
/// 例えばMAX_INDEX_BITS=20であれば、1,048,576まで行列が確保できる。
#define MAX_INDEX_BITS (20)

/// 行列クラスのインターフェースを規定する
///
class TlMatrixObject {
  // ---------------------------------------------------------------------------
 public:
  /// 行数、列数を表す整数型(範囲外を示す値として-1を取る場合がある)
  typedef signed int index_type;
  /// 数値配列の総数を表す整数型(範囲外を示す値として-1を取る場合がある)
  typedef signed long size_type;

  // ---------------------------------------------------------------------------
 public:
  /// 行列要素を格納するための構造体
  ///
  /// TlCommunicate で通信可能
  struct MatrixElement {
    typedef TlMatrixObject::index_type index_type;

   public:
    MatrixElement(index_type r = 0, index_type c = 0, double v = 0.0)
        : row(r), col(c), value(v) {}

   public:
    index_type row;
    index_type col;
    double value;
  };

  // ---------------------------------------------------------------------------
 public:
  enum MatrixType {
    UNDEFINED = -1,
    RSFD = 0,  /// RSFD: Row-oriented Standard Full Dens-matrix
    CSFD = 1,  /// Coulmn-oriented Standard Full Dens-matrix
    RLHD = 2,  /// Row-oriented Lower Half Dens matrix
    RUHD = 3,  /// Row-oriented Upper Half Dens matrix
    CLHD = 4,  /// Coulmn-oriented Lower Half Dens matrix
    CUHD = 5,  /// Coulmn-oriented Upper Half Dens matrix
    COOF = 6,  /// list of tuple (row, col, value) full sparse matrix
    COOS = 7   /// list of tuple (row, col, value) half sparse matrix
  };

  // ---------------------------------------------------------------------------
 public:
  explicit TlMatrixObject(const MatrixType matrixType = UNDEFINED,
                          const index_type row = 1, index_type col = 1);

  virtual ~TlMatrixObject() {}

  // ---------------------------------------------------------------------------
 public:
  virtual MatrixType getType() const { return this->matrixType_; };
  virtual index_type getNumOfRows() const { return this->row_; };
  virtual index_type getNumOfCols() const { return this->col_; };

  // ---------------------------------------------------------------------------
  // I/O
  // ---------------------------------------------------------------------------
 public:
  /// 指定されたパスから内容を読み込む
  ///
  /// @retval true 成功
  /// @retval false 失敗
  virtual bool load(const std::string& path) = 0;

  /// 指定されたパスに内容を書き出す
  ///
  /// @retval true 成功
  /// @retval false 失敗
  virtual bool save(const std::string& path) const = 0;

 protected:
  size_type getNumOfElements_CSFD() const;
  size_type getNumOfElements_RSFD() const;
  size_type getNumOfElements_RLHD() const;

 protected:
  size_type getIndex_RSFD(index_type row, index_type col) const;
  size_type getIndex_CSFD(index_type row, index_type col) const;
  size_type getIndex_RLHD(index_type row, index_type col) const;

  // ---------------------------------------------------------------------------
  // get/set
  // ---------------------------------------------------------------------------
 public:
  /** 要素を返す(読み取り専用)
   *
   *  内部では、行列要素を(2次元配列ではなく)
   *  1次元配列として保持しているので、
   *  他のメンバ関数内でもこのメンバ関数が呼ばれる。
   *
   *  @param[in] row 行数
   *  @param[in] col 列数
   *  @return 要素
   */
  virtual double get(index_type row, index_type col) const = 0;

  /** 要素に値を代入する
   *
   *  内部では、行列要素を(2次元配列ではなく)
   *  1次元配列として保持しているので、
   *  他のメンバ関数内でもこのメンバ関数が呼ばれる。
   *
   *  @param[in] row 行数
   *  @param[in] col 列数
   *  @return 要素
   */
  virtual void set(index_type row, index_type col, double value) = 0;

  /** 要素に値を加える
   *
   *  内部では、行列要素を(2次元配列ではなく)
   *  1次元配列として保持しているので、
   *  他のメンバ関数内でもこのメンバ関数が呼ばれる。
   *
   *  @param[in] row 行数
   *  @param[in] col 列数
   *  @return 要素
   */
  virtual void add(index_type row, index_type col, double value) = 0;

  // ---------------------------------------------------------------------------
  // checked.
  // ---------------------------------------------------------------------------
 public:
  virtual double getLocal(index_type row, index_type col) const {
    return this->get(row, col);
  }

  // ---------------------------------------------------------------------------
  // deprecated
  // ---------------------------------------------------------------------------
 private:
  /// インスタンスのメモリサイズを返す
  //  virtual std::size_t getMemSize() const = 0;

 protected:
  MatrixType matrixType_;
  index_type row_;
  index_type col_;

 protected:
  TlLogging& log_;
};

std::ostream& operator<<(std::ostream& stream, const TlMatrixObject& mat);

#endif  // TL_MATRIX_OBJECT_H
