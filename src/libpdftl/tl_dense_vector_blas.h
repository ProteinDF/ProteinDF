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

#ifndef TLVECTOR_H
#define TLVECTOR_H

#include <cassert>
#include <valarray>

#include "TlUtils.h"
#include "tl_sparse_vector.h"
#include "tl_vector_abstract.h"

class TlDenseGeneralMatrix_BLAS_old;
class TlDenseSymmetricMatrix_BLAS_Old;
class TlDenseGeneralMatrix_blacs;
class TlDenseSymmetricMatrix_blacs;
class TlCommunicate;

/// double専用のベクトルクラス
class TlVector_BLAS : public TlVectorAbstract {
  friend class TlCommunicate;

 public:
  /// コンストラクタ
  ///
  /// @param[in] size ベクトルの要素数
  explicit TlVector_BLAS(TlVectorAbstract::index_type size = 0);

  /// コンストラクタ
  ///
  /// @param[in] rhs ベクトル
  TlVector_BLAS(const double* p, TlVectorAbstract::size_type size);

  /// コンストラクタ(変換)
  TlVector_BLAS(const std::vector<double>& rhs);

  /// コピーコンストラクタ
  TlVector_BLAS(const TlVector_BLAS& rhs);

  /// デストラクタ
  virtual ~TlVector_BLAS();

 public:
  /// 要素数を返す
  ///
  /// @return 要素数
  size_type getSize() const;

  /// サイズを変更する
  ///
  /// 現在の要素数よりも大きいサイズを指定された場合は、要素数が拡張される.
  /// 現在の要素数よりも小さいサイズを指定された場合は、切り詰められる.
  /// @param[in] size 要素数
  void resize(TlVectorAbstract::index_type size);

  /// 末尾に値を追加する
  ///
  /// @param[in] value 末尾に追加する値
  void push_back(double value);

  /// 要素のうち、要素の絶対値の最大を返す.
  ///
  /// @return 要素の絶対値の最大値
  double getMaxAbsoluteElement() const;

  TlVector_BLAS& dotInPlace(const TlVector_BLAS& rhs);

  /// 全要素の和を返す
  ///
  /// @return 全要素の和
  double sum() const;

  /// 要素を値の大きい順にソートする
  void sortByGrater();

  /// 大きな要素のインデックスを返す
  /// @param startIndex[in]: 探査を開始するインデックス
  std::vector<TlVectorAbstract::size_type>::const_iterator argmax(
      const std::vector<TlVectorAbstract::size_type>::const_iterator& begin,
      const std::vector<TlVectorAbstract::size_type>::const_iterator& end) const;

  /// ユークリッドノルムを返す
  double norm() const;

  /// ユークリッドノルムの2乗を返す
  double norm2() const;

 public:
  // static bool isLoadable(const std::string& sFilePath);
  // static bool isLoadable(std::ifstream& ifs);
  // static bool getHeaderInfo(std::ifstream& ifs,
  //                           TlVectorAbstract::index_type* pnSize);

  /// ファイルからデータを読み込む
  ///
  /// @param[in] sFilePath ファイルのパス名
  /// @retval true 読み込み成功
  /// @retval false 読み込み失敗
  virtual bool load(const std::string& sFilePath);

  /// std::ifstreamオブジェクトからデータを読み込む
  ///
  /// @param[in,out] ifs std::ifstreamオブジェクト
  /// @retval true 読み込み成功
  /// @retval false 読み込み失敗
  // virtual bool load(std::ifstream& ifs);

  virtual bool loadText(const std::string& filePath);

  /// ファイルへデータを出力する
  ///
  /// @param[in] sFilePath ファイルのパス名
  /// @retval true 読み込み成功
  /// @retval false 読み込み失敗
  virtual bool save(const std::string& sFilePath) const;

/// std::ofstreamオブジェクトへデータを出力する
///
/// @param[out] ofs std::ofstreamオブジェクト
/// @retval true 出力成功
/// @retval false 出力失敗
// virtual bool save(std::ofstream& ofs) const;

#ifdef HAVE_HDF5
  bool saveHdf5(const std::string& filepath, const std::string& h5path) const;
  bool loadHdf5(const std::string& filepath, const std::string& h5path);
#endif  // HAVE_HDF5

  virtual void outputText(std::ostream& os) const;

 public:
  /// オブジェクトの内容をテキスト出力する
  ///
  /// @param[in] out 出力用の<< 演算子を持つオブジェクト
  template <typename T>
  void print(T& out) const;

 public:
  /// 内容を0(zero)で埋める
  void zeroClear();

  double get(TlVectorAbstract::index_type index) const;

  void set(TlVectorAbstract::index_type index, double value);

  virtual void add(const TlVectorAbstract::index_type index, const double value);

  /// 要素を返す(読み取り専用)
  ///
  /// @param[in] index 要素番号
  /// @return 要素
  double operator[](const TlVectorAbstract::index_type index) const;

  /// 要素を返す(代入用)
  ///
  /// @param[in] index 要素番号
  /// @return 要素
  double& operator[](const TlVectorAbstract::index_type index);

  /// 代入演算子
  TlVector_BLAS& operator=(const TlVector_BLAS& rhs);

  /// 要素数は同じである必要がある
  TlVector_BLAS& operator+=(const TlVector_BLAS& rhs);

  TlVector_BLAS& operator+=(const TlSparseVector& rhs);

  /// 要素数は同じである必要がある
  TlVector_BLAS& operator-=(const TlVector_BLAS& rhs);

  /// 全要素に対し乗算する
  ///
  /// @param[in] rhs 乗算する値
  TlVector_BLAS& operator*=(const double& rhs);

  /// 全要素に対し除算する
  ///
  /// @param[in] rhs 除算する値
  TlVector_BLAS& operator/=(const double& rhs);

  /// ベクトル同士の加算
  friend TlVector_BLAS operator+(const TlVector_BLAS& X,
                                 const TlVector_BLAS& Y);

  /// ベクトル同士の減算
  friend TlVector_BLAS operator-(const TlVector_BLAS& X,
                                 const TlVector_BLAS& Y);

  /// ベクトルの定数倍
  friend TlVector_BLAS operator*(const TlVector_BLAS& X, const double& Y);

  /// ベクトルの除算
  friend TlVector_BLAS operator*(const double& X, const TlVector_BLAS& Y) {
    return (Y * X);
  }

  /// 行列と縦ベクトルの積
  friend TlVector_BLAS operator*(const TlDenseGeneralMatrix_BLAS_old& X,
                                 const TlVector_BLAS& Y);

  /// 横ベクトルと行列の積
  friend TlVector_BLAS operator*(const TlVector_BLAS& Y,
                                 const TlDenseGeneralMatrix_BLAS_old& X);

  /// ベクトルの内積
  friend double operator*(const TlVector_BLAS& X, const TlVector_BLAS& Y);

  /// ベクトルの除算
  ///
  /// @param[in] X TlVector_BLASオブジェクト
  /// @param[in] Y 除算する値
  friend TlVector_BLAS operator/(const TlVector_BLAS& X, const double& Y) {
    assert(std::fabs(Y) < 1.0E-20);
    return (X * (1.0 / Y));
  }

 protected:
  void initialize(bool isZeroClear = true);

  virtual void destroy();

 protected:
  size_type size_;

  // BufferType m_aVector; /// 要素格納用
  double* data_;  /// 要素格納用

#ifdef HAVE_LAPACK
  friend TlVector_BLAS multiplicationByLapack(
      const TlDenseGeneralMatrix_BLAS_old& A, const TlVector_BLAS& X);
  friend TlVector_BLAS multiplicationByLapack(
      const TlVector_BLAS& X, const TlDenseGeneralMatrix_BLAS_old& A);
  friend TlVector_BLAS multiplicationByLapack(
      const TlDenseSymmetricMatrix_BLAS_Old& A, const TlVector_BLAS& X);
  friend bool diagonalByLapack(const TlDenseSymmetricMatrix_BLAS_Old& inMatrix,
                               TlVector_BLAS* outEigVal,
                               TlDenseGeneralMatrix_BLAS_old* outEigVec);
#endif  // HAVE_LAPACK

#ifdef HAVE_SCALAPACK
  friend bool diagonalByScaLapack_QR(
      const TlDenseSymmetricMatrix_blacs& inMatrix, TlVector_BLAS* outEigVal,
      TlDenseGeneralMatrix_blacs* outEigVec);
  friend bool diagonalByScaLapack_DC(
      const TlDenseSymmetricMatrix_blacs& inMatrix, TlVector_BLAS* outEigVal,
      TlDenseGeneralMatrix_blacs* outEigVec);
#endif  // HAVE_SCALAPACK
};

////////////////////////////////////////////////////////////////////////
// inline functions
inline TlVector_BLAS::size_type TlVector_BLAS::getSize() const {
  return this->size_;
}

inline void TlVector_BLAS::zeroClear() {
  const size_type size = this->getSize();
  std::fill(this->data_, this->data_ + size, 0.0);
}

inline double TlVector_BLAS::get(const TlVectorAbstract::index_type index) const {
  assert(index < this->getSize());
  return this->data_[index];
}

inline void TlVector_BLAS::set(const TlVectorAbstract::index_type index,
                               const double value) {
  assert(index < this->getSize());

#pragma omp critical(TlVector_BLAS__set)
  { this->data_[index] = value; }
}

inline void TlVector_BLAS::add(const TlVectorAbstract::index_type index,
                               const double value) {
  assert(index < this->getSize());

#pragma omp atomic
  this->data_[index] += value;
}

inline double TlVector_BLAS::operator[](
    const TlVectorAbstract::index_type index) const {
  assert(index < this->getSize());

  return (this->data_[index]);
}

inline double& TlVector_BLAS::operator[](
    const TlVectorAbstract::index_type index) {
  assert(index < this->getSize());

  return (this->data_[index]);
}

template <typename T>
void TlVector_BLAS::print(T& out) const {
  const size_type size = this->getSize();
  for (size_type ord = 0; ord < size; ord += 10) {
    out << "\n";
    for (size_type j = ord; (j < ord + 10) && (j < size); ++j) {
      out << TlUtils::format("   %5d th", j + 1);
    }
    out << "\n";

    for (size_type j = ord; (j < ord + 10) && (j < size); ++j) {
      out << "-----------";
    }
    out << "----\n\n";

    for (size_type j = ord; (j < ord + 10) && (j < size); ++j) {
      out << TlUtils::format(" %10.6lf", this->get(j));
    }
    out << "\n\n";
  }

  out.flush();
}

#endif  // TLVECTOR_H
