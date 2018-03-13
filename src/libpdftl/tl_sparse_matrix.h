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

#ifndef TLSPARSEMATRIX_H
#define TLSPARSEMATRIX_H

#ifdef HAVE_CONFIG_H
#include "config.h"  // this file created by autotools
#endif               // HAVE_CONFIG_H

#include <cassert>
#include <functional>
#include <iostream>

// -------------------------------------
#ifdef HAVE_UNORDERED_MAP
#include <unordered_map>
#elifdef HAVE_TR1_UNORDERED_MAP
#include <tr1/unordered_map>
#elifdef HAVE_GOOGLE_SPARSE_HASH_MAP
#include <google/sparse_hash_map>
#else
#include <map>
#endif
// -------------------------------------

#include "TlUtils.h"
#include "tl_matrix_object.h"
#include "tl_dense_vector_blas.h"

struct TlMatrixElement {
  typedef TlMatrixObject::index_type index_type;

 public:
  index_type row;
  index_type col;
  double value;
};

/// 疎行列クラス
class TlSparseMatrix : public TlMatrixObject {
 public:
  /// 行列オブジェクトを作成する
  ///
  /// @param[in] row 作成する行列の行数
  /// @param[in] col 作成する行列の列数
  explicit TlSparseMatrix(index_type row = 1, index_type col = 1);

  /// コピーコンストラクタ
  TlSparseMatrix(const TlSparseMatrix& rhs);

  /// デストラクタ
  virtual ~TlSparseMatrix();

 public:
  struct Index2 {
   public:
    Index2(index_type r = 0, index_type c = 0) : row(r), col(c) {}

    bool operator<(const Index2& rhs) const {
      bool answer = false;
      if ((this->row < rhs.row) ||
          ((this->row == rhs.row) && (this->col < rhs.col))) {
        answer = true;
      }

      return answer;
    }

    bool operator>(const Index2& rhs) const {
      bool answer = false;
      if ((this->row > rhs.row) ||
          ((this->row == rhs.row) && (this->col > rhs.col))) {
        answer = true;
      }

      return answer;
    }

    bool operator==(const Index2& rhs) const {
      return ((this->row == rhs.row) && (this->col == rhs.col));
    }

    bool operator!=(const Index2& rhs) const {
      return !(this->operator==(rhs));
    }

   public:
    index_type row;
    index_type col;
  };

 public:
  // typedef std::size_t KeyType;
  typedef Index2 KeyType;

  // #ifdef HAVE_UNORDERED_MAP
  //     typedef std::unordered_map<KeyType, double> SparseMatrixData;
  // #elifdef HAVE_TR1_UNORDERED_MAP
  //     typedef std::tr1::unordered_map<KeyType, double> SparseMatrixData;
  // #elifdef HAVE_GOOGLE_SPARSE_HASH_MAP
  //     typedef google::sparse_hash_map<KeyType, double> SparseMatrixData;
  // #else
  //     typedef std::map<KeyType, double> SparseMatrixData;
  //     #define TSM_DATATYPE_BINTREE 1
  // #endif
  typedef std::map<KeyType, double> SparseMatrixData;

  typedef SparseMatrixData::const_iterator const_iterator;
  typedef SparseMatrixData::iterator iterator;

  // operator
 public:
  /// 最初のイテレータを返す
  const_iterator begin() const { return m_aMatrix.begin(); }

  /// 最初のイテレータを返す
  iterator begin() { return m_aMatrix.begin(); }

  /// 最後のイテレータを返す
  const_iterator end() const { return m_aMatrix.end(); }

  /// 最後のイテレータを返す
  iterator end() { return m_aMatrix.end(); }

  /// オブジェクトの内容を破棄する
  virtual void clear();

  virtual void zeroClear() { this->clear(); }

  /// 要素を削除する
  virtual void erase(index_type row, index_type col);
  virtual void erase(iterator p);

  /// 行数を返す
  ///
  /// @return 行数
  virtual index_type getNumOfRows() const;

  /// 列数を返す
  ///
  /// @return 列数
  virtual index_type getNumOfCols() const;

  /// 要素の数を返す
  ///
  /// @return 要素の数
  virtual size_type getSize() const;

  virtual std::size_t getMemSize() const;

  /// 行列のサイズを変更する
  ///
  /// 行列を大きくする場合、追加される要素は0で初期化される。
  /// 行列を小さくする場合、切り詰められる要素は破棄される。
  /// @param[in] row 行数
  /// @param[in] col 列数
  virtual void resize(index_type row, index_type col);

  double pop(index_type* pRow, index_type* pCol);

  /// 要素を返す(読み取り専用)
  ///
  /// 内部では、行列要素を(2次元配列ではなく)
  /// 1次元配列として保持しているので、
  /// 他のメンバ関数内でもこのメンバ関数が呼ばれる。
  ///
  /// @param[in] row 行数
  /// @param[in] col 列数
  /// @return 要素
  virtual double get(index_type row, index_type col) const;

  /** 要素を返す(代入可能)
   *
   *  内部では、行列要素を(2次元配列ではなく)
   *  1次元配列として保持しているので、
   *  他のメンバ関数内でもこのメンバ関数が呼ばれる。
   *
   *  @param[in] row 行数
   *  @param[in] col 列数
   *  @return 要素
   */
  virtual void set(const index_type row, const index_type col,
                   const double value);
  // virtual void set(const std::pair<unsigned long, double>& obj);

  virtual void add(const index_type row, const index_type col,
                   const double value);
  // virtual void add(const std::pair<unsigned long, double>& obj);

  void add(const std::vector<MatrixElement>& elements);

  /** 指定された要素が存在すればtrueを返す
   *
   *  @retval true 要素が存在する
   *  @retval false 要素が存在しない
   */
  virtual bool hasKey(index_type row, index_type col) {
    return (this->m_aMatrix.find(KeyType(row, col)) != this->m_aMatrix.end());
  }

  virtual void merge(const TlSparseMatrix& rhs);

  /// 代入演算子
  TlSparseMatrix& operator=(const TlSparseMatrix& rhs);

  /// 行列を定数倍する
  ///
  /// @param[in] rhs 定数倍の値
  /// @return 計算後のオブジェクト
  TlSparseMatrix& operator*=(const double& rhs);

  /// 行列を定数で割る
  ///
  /// @param[in] rhs 割る定数の値
  /// @return 計算後のオブジェクト
  TlSparseMatrix& operator/=(const double& rhs);

  /// 指定した行の要素から構成されるベクトルを返す
  ///
  /// @param[in] row 指定する行
  virtual TlVector_BLAS getRowVector(index_type row) const;

  /// 指定した列の要素から構成されるベクトルを返す
  ///
  /// @param[in] col 指定する列
  virtual TlVector_BLAS getColVector(index_type col) const;

  /// Hadamard product
  const TlSparseMatrix& dot(const TlSparseMatrix& X);

  double sum() const;

  std::vector<int> getRowIndexList() const;
  std::vector<int> getColIndexList() const;

  std::vector<MatrixElement> getMatrixElements() const;

 public:
  /// オブジェクトの内容をテキスト出力する
  ///
  /// @param[in] out 出力用の<< 演算子を持つオブジェクト
  template <typename T>
  void print(T& out) const;

 public:
  // unsigned long index(const index_type row, const index_type col) const;
  // void index(const KeyType& i, index_type* pRow, index_type* pCol) const;

 public:
  virtual bool load(const std::string& path);
  virtual bool save(const std::string& path) const;

 protected:
  bool load(std::ifstream& ifs);
  bool save(std::ofstream& ofs) const;

 protected:
  index_type m_nRows;                  /// 行数
  index_type m_nCols;                  /// 列数
  mutable SparseMatrixData m_aMatrix;  /// 行列要素

  friend class TlCommunicate;
};

////////////////////////////////////////////////////////////////////////
// inline functions

inline void TlSparseMatrix::set(const index_type row, const index_type col,
                                const double value) {
  assert((0 <= row) && (row < this->getNumOfRows()));
  assert((0 <= col) && (col < this->getNumOfCols()));

#pragma omp critical(TlSparseMatrix__set)
  { this->m_aMatrix[KeyType(row, col)] = value; }
}

inline void TlSparseMatrix::add(const index_type row, const index_type col,
                                const double value) {
  assert((0 <= row) && (row < this->getNumOfRows()));
  assert((0 <= col) && (col < this->getNumOfCols()));

//#pragma omp atomic
#pragma omp critical(TlSparseMatrix__add)
  { this->m_aMatrix[KeyType(row, col)] += value; }
}

inline double TlSparseMatrix::get(const index_type row,
                                  const index_type col) const {
  assert((0 <= row) && (row < this->m_nRows));
  assert((0 <= col) && (col < this->m_nCols));

  double answer = 0.0;
  const_iterator p = this->m_aMatrix.find(KeyType(row, col));
  if (p != this->m_aMatrix.end()) {
    answer = p->second;
  }

  return answer;
}

template <typename T>
void TlSparseMatrix::print(T& out) const {
  const int nNumOfRows = this->getNumOfRows();
  const int nNumOfCols = this->getNumOfCols();

  for (int ord = 0; ord < nNumOfCols; ord += 10) {
    out << "       ";
    for (int j = ord; j < ord + 10 && j < nNumOfCols; j++) {
      out << TlUtils::format("   %5d th", j + 1);
    }
    out << "\n ----";

    for (int j = ord; ((j < ord + 10) && (j < nNumOfCols)); j++) {
      out << "-----------";
    }
    out << "----\n";

    for (int i = 0; i < nNumOfRows; i++) {
      out << TlUtils::format(" %5d  ", i + 1);

      for (int j = ord; ((j < ord + 10) && (j < nNumOfCols)); j++) {
        out << TlUtils::format(" %10.6lf", this->get(i, j));
      }
      out << "\n";
    }
    out << "\n\n";
  }
  out.flush();
}

#endif  // TLMATRIX_SPARSE_H
