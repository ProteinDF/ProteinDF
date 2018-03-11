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

#ifndef TLSPARSESYMMETRICMATRIX_H
#define TLSPARSESYMMETRICMATRIX_H

#include <algorithm>
#include "TlSparseMatrix.h"

/// 対称疎行列クラス
///
/// 常に(row >= col)の要素にアクセスされる
class TlSparseSymmetricMatrix : public TlSparseMatrix {
 public:
  /// 行列オブジェクトを作成する
  ///
  /// @param[in] size 作成する行列の次元数
  explicit TlSparseSymmetricMatrix(int size = 1);

  /// コピーコンストラクタ
  TlSparseSymmetricMatrix(const TlSparseSymmetricMatrix& rhs);

  /// デストラクタ
  virtual ~TlSparseSymmetricMatrix();

 public:
  /// 行列のサイズを変更する
  ///
  /// 行列を大きくする場合、追加される要素は0で初期化される。
  /// 行列を小さくする場合、切り詰められる要素は破棄される。
  ///
  /// @param[in] row 行数
  /// @param[in] col 列数
  virtual void resize(index_type row, index_type col) {
    assert(row == col);
    this->resize(row);
  }

  /// 行列のサイズを変更する
  ///
  /// 行列を大きくする場合、追加される要素は0で初期化される。
  /// 行列を小さくする場合、切り詰められる要素は破棄される。
  /// @param[in] size 要素数
  virtual void resize(index_type size);

  /// 要素を返す(読み取り専用)
  ///
  /// 内部では、行列要素を(2次元配列ではなく)
  /// 1次元配列として保持しているので、
  /// 他のメンバ関数内でもこのメンバ関数が呼ばれる。
  ///
  /// @param[in] row 行数
  /// @param[in] col 列数
  /// @return 要素
  virtual double operator()(index_type row, index_type col) const;
  virtual double get(index_type row, index_type col) const;

  /// 要素を返す(代入可能)
  ///
  /// 内部では、行列要素を(2次元配列ではなく)
  /// 1次元配列として保持しているので、
  /// 他のメンバ関数内でもこのメンバ関数が呼ばれる。
  ///
  ///  @param[in] row 行数
  /// @param[in] col 列数
  /// @return 要素
  virtual double& operator()(index_type row, index_type col);
  virtual void set(index_type row, index_type col, double value);
  // virtual void set(const std::pair<unsigned long, double>& obj) {
  //     TlSparseMatrix::set(obj);
  // }

  virtual void add(index_type row, index_type col, double value);
  // virtual void add(const std::pair<unsigned long, double>& obj) {
  //     TlSparseMatrix::add(obj);
  // }

  /** 指定された要素が存在すればtrueを返す
   *
   *  @retval true 要素が存在する
   *  @retval false 要素が存在しない
   */
  virtual bool hasKey(index_type row, index_type col) {
    if (row < col) {
      std::swap(row, col);
    }
    return (this->m_aMatrix.find(KeyType(row, col)) != this->m_aMatrix.end());
  }

  /// 指定した行の要素から構成されるベクトルを返す
  ///
  /// @param[in] nRow 指定する行
  virtual TlVector getRowVector(int nRow) const;

  /// 指定した列の要素から構成されるベクトルを返す
  ///
  /// @param[in] nCol 指定する列
  virtual TlVector getColVector(int nCol) const;

  /// Hadamard product
  // const TlSparseSymmetricMatrix& dot(const TlSparseSymmetricMatrix& X);

 public:
  /// オブジェクトの内容をテキスト出力する
  ///
  /// @param[in] out 出力用の<< 演算子を持つオブジェクト
  template <typename T>
  void print(T& out) const;
};

////////////////////////////////////////////////////////////////////////
// inline functions
//
inline double TlSparseSymmetricMatrix::operator()(int row, int col) const {
  return this->get(row, col);
}

inline double TlSparseSymmetricMatrix::get(int row, int col) const {
  assert(0 <= row && row < this->m_nRows);
  assert(0 <= col && col < this->m_nCols);

  if (row < col) {
    std::swap(row, col);
  }

  return TlSparseMatrix::get(row, col);
}

inline double& TlSparseSymmetricMatrix::operator()(int row, int col) {
  assert(0 <= row && row < this->m_nRows);
  assert(0 <= col && col < this->m_nCols);

  if (row < col) {
    std::swap(row, col);
  }

  return TlSparseMatrix::operator()(row, col);
}

inline void TlSparseSymmetricMatrix::set(int row, int col, const double value) {
  assert(0 <= row && row < this->m_nRows);
  assert(0 <= col && col < this->m_nCols);

  if (row < col) {
    std::swap(row, col);
  }

  TlSparseMatrix::set(row, col, value);
}

inline void TlSparseSymmetricMatrix::add(int row, int col, const double value) {
  assert(0 <= row && row < this->m_nRows);
  assert(0 <= col && col < this->m_nCols);

  if (row < col) {
    std::swap(row, col);
  }

  TlSparseMatrix::add(row, col, value);
}

template <typename T>
void TlSparseSymmetricMatrix::print(T& out) const {
  const int nNumOfDim = this->getNumOfRows();  // == this->getNumOfCols()

  out << "\n\n";
  for (int ord = 0; ord < nNumOfDim; ord += 10) {
    out << "       ";
    for (int j = ord; ((j < ord + 10) && (j < nNumOfDim)); j++) {
      out << TlUtils::format("   %5d th", j + 1);
    }
    out << "\n"
        << " ----";

    for (int j = ord; ((j < ord + 10) && (j < nNumOfDim)); j++) {
      out << "-----------";
    }
    out << "----\n";

    for (int i = 0; i < nNumOfDim; i++) {
      out << TlUtils::format(" %5d  ", i + 1);

      for (int j = ord; ((j < ord + 10) && (j < nNumOfDim)); j++) {
        if (j > i) {
          out << "    ----   ";
        } else {
          out << TlUtils::format(" %10.6lf", (*this)(i, j));
        }
      }
      out << "\n";
    }
    out << "\n\n";
  }
  out.flush();
}

#endif  // TLSPARSESYMMETRICMATRIX_H
