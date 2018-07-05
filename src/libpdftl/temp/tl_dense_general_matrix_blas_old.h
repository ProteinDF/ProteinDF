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

#ifndef TL_DENSE_GENERAL_MATRIX_BLAS_H
#define TL_DENSE_GENERAL_MATRIX_BLAS_H

#ifdef HAVE_CONFIG_H
#include "config.h"  // this file created by autotools
#endif               // HAVE_CONFIG_H

#include <cassert>
#include <fstream>
#include <sstream>
#include <valarray>
#include "TlSerializeData.h"
#include "TlUtils.h"
#include "tl_dense_general_matrix_abstract.h"

class TlDenseSymmetricMatrix_BLAS_Old;

/// 密行列クラス
/// 最低でも32bit環境を想定している。
/// 要素の指定はlongではなく、intを用いている。
/// int(32bit)の範囲は-2,147,483,648〜2,147,483,647であり、
/// 1,000残基規模(約10万軌道)の計算には十分耐えられる。
class TlDenseGeneralMatrix_BLAS_old : public TlDenseGeneralMatrixAbstract {
  friend class TlCommunicate;
  friend class TlVector_BLAS;

 public:
  /// 行列オブジェクトを作成する
  ///
  /// 行列の要素は0で初期化される。
  /// @param[in] row 作成する行列の行数
  /// @param[in] col 作成する行列の列数
  explicit TlDenseGeneralMatrix_BLAS_old(TlMatrixObject::index_type row = 1,
                                     TlMatrixObject::index_type col = 1);

  TlDenseGeneralMatrix_BLAS_old(const TlSerializeData& data);

  /// 行列オブジェクトを作成する
  ///
  /// 行列の要素はベクトルから作成される。
  /// ベクトルサイズはrowとcolの積でなければならない。
  /// @param[in] rVector ベクトル
  /// @param[in] row 作成する行列の行数
  /// @param[in] col 作成する行列の列数
  TlDenseGeneralMatrix_BLAS_old(const TlVector_BLAS& rVector,
                            TlMatrixObject::index_type row,
                            TlMatrixObject::index_type col);

 public:
  /// コピーコンストラクタ
  TlDenseGeneralMatrix_BLAS_old(const TlDenseGeneralMatrix_BLAS_old& rhs);

  /// 対称行列オブジェクトからTlDenseGeneralMatrix_BLAS_oldオブジェクトを作成する
  TlDenseGeneralMatrix_BLAS_old(const TlDenseSymmetricMatrix_BLAS_Old& rhs);

  virtual ~TlDenseGeneralMatrix_BLAS_old();

 public:
  /// 行列のサイズを変更する
  ///
  /// 行列を大きくする場合、追加される要素は0で初期化される。
  /// 行列を小さくする場合、切り詰められる要素は破棄される。
  /// @param[in] row 行数
  /// @param[in] col 列数
  void resize(TlMatrixObject::index_type nRow, TlMatrixObject::index_type nCol);

  /// 行列要素をベクトルにして返す
  ///
  /// 3x3 の行列であれば以下の順に返す。
  /// [ 0 3 6 ]
  /// [ 1 4 7 ]
  /// [ 2 5 8 ]
  virtual TlVector_BLAS getVector() const;

  /// インスタンスのメモリサイズを返す
  virtual std::size_t getMemSize() const;

  // ---------------------------------------------------------------------------
  // operators
  // ---------------------------------------------------------------------------
  /// 要素を返す(読み取り専用)
  ///
  /// 内部では、行列要素を(2次元配列ではなく)
  /// 1次元配列として保持しているので、
  /// 他のメンバ関数内でもこのメンバ関数が呼ばれる。
  /// @param[in] row 行数
  /// @param[in] col 列数
  /// @return 要素
  // virtual double operator()(TlMatrixObject::index_type row,
  //                           TlMatrixObject::index_type col) const;

  virtual double get(const TlMatrixObject::index_type row,
                     const TlMatrixObject::index_type col) const;

  /// 要素を返す(代入可能)
  ///
  /// 内部では、行列要素を(2次元配列ではなく)
  /// 1次元配列として保持しているので、
  /// 他のメンバ関数内でもこのメンバ関数が呼ばれる。
  /// @param[in] row 行数
  /// @param[in] col 列数
  /// @return 要素
  // virtual double& operator()(TlMatrixObject::index_type row,
  //                            TlMatrixObject::index_type col);

  virtual void set(const TlMatrixObject::index_type row,
                   const TlMatrixObject::index_type col, const double value);
  virtual void add(const TlMatrixObject::index_type row,
                   const TlMatrixObject::index_type col, const double value);

  /// 代入演算子
  TlDenseGeneralMatrix_BLAS_old& operator=(const TlDenseGeneralMatrix_BLAS_old& rhs);

  /// 代入演算子
  /// @param[in] rhs 対称行列
  TlDenseGeneralMatrix_BLAS_old& operator=(const TlDenseSymmetricMatrix_BLAS_Old& rhs);

  TlDenseGeneralMatrix_BLAS_old& operator+=(const TlDenseGeneralMatrix_BLAS_old& rhs);
  TlDenseGeneralMatrix_BLAS_old& operator+=(const TlDenseSymmetricMatrix_BLAS_Old& rhs);
  TlDenseGeneralMatrix_BLAS_old& operator-=(const TlDenseGeneralMatrix_BLAS_old& rhs);
  TlDenseGeneralMatrix_BLAS_old& operator-=(const TlDenseSymmetricMatrix_BLAS_Old& rhs);
  TlDenseGeneralMatrix_BLAS_old& operator*=(const TlDenseGeneralMatrix_BLAS_old& rhs);
  TlDenseGeneralMatrix_BLAS_old& operator*=(const TlDenseSymmetricMatrix_BLAS_Old& rhs);

  /// 行列を定数倍する
  ///
  /// @param[in] rhs 定数倍の値
  /// @return 計算後のオブジェクト
  TlDenseGeneralMatrix_BLAS_old& operator*=(const double& rhs);

  /// 行列を定数で割る
  ///
  /// @param[in] rhs 割る定数の値
  /// @return 計算後のオブジェクト
  TlDenseGeneralMatrix_BLAS_old& operator/=(const double& rhs);

  friend TlDenseGeneralMatrix_BLAS_old operator+(
      const TlDenseGeneralMatrix_BLAS_old& X, const TlDenseGeneralMatrix_BLAS_old& Y);
  friend TlDenseGeneralMatrix_BLAS_old operator-(
      const TlDenseGeneralMatrix_BLAS_old& X, const TlDenseGeneralMatrix_BLAS_old& Y);
  friend TlDenseGeneralMatrix_BLAS_old operator*(
      const TlDenseGeneralMatrix_BLAS_old& X, const TlDenseGeneralMatrix_BLAS_old& Y);
  friend TlDenseGeneralMatrix_BLAS_old operator*(const TlDenseGeneralMatrix_BLAS_old& X,
                                             double Y);
  friend TlDenseGeneralMatrix_BLAS_old operator*(
      double X, const TlDenseGeneralMatrix_BLAS_old& Y) {
    return (Y * X);
  };

  friend TlVector_BLAS operator*(const TlDenseGeneralMatrix_BLAS_old& X,
                                 const TlVector_BLAS& Y);
  friend TlVector_BLAS operator*(const TlVector_BLAS& Y,
                                 const TlDenseGeneralMatrix_BLAS_old& X);

  // ---------------------------------------------------------------------------
  // operations
  // ---------------------------------------------------------------------------

  /// 指定した行の要素から構成されるベクトルを返す
  ///
  /// @param[in] row 指定する行
  virtual TlVector_BLAS getRowVector(TlMatrixObject::index_type row) const;

  /// 指定した列の要素から構成されるベクトルを返す
  ///
  /// @param[in] col 指定する行
  virtual TlVector_BLAS getColVector(TlMatrixObject::index_type col) const;

  /// ブロック行列を返す
  ///
  /// @param[in] row 始点となる行
  /// @param[in] col 始点となる列
  /// @param[in] row_distance 取得する行数
  /// @param[in] col_distance 取得する列数
  /// @return row_distance × col_distance 次元のブロック行列
  virtual TlDenseGeneralMatrix_BLAS_old getBlockMatrix(
      TlMatrixObject::index_type nRow, TlMatrixObject::index_type nCol,
      TlMatrixObject::index_type nRowDistance,
      TlMatrixObject::index_type nColDistance) const;

  /// 行列要素を指定された位置に上書きする
  ///
  /// @param[in] row 始点となる行
  /// @param[in] col 始点となる列
  /// @param[in] matrix 行列要素
  virtual void setBlockMatrix(TlMatrixObject::index_type nRow,
                              TlMatrixObject::index_type nCol,
                              const TlDenseGeneralMatrix_BLAS_old& matrix);

  /// 行列要素を指定された位置に加算する
  ///
  /// @param[in] row 始点となる行
  /// @param[in] col 始点となる列
  /// @param[in] matrix 行列要素
  virtual void addBlockMatrix(int nRow, int nCol,
                              const TlDenseGeneralMatrix_BLAS_old& matrix);

  /// Hadamard product
  const TlDenseGeneralMatrix_BLAS_old& dotInPlace(
      const TlDenseGeneralMatrix_BLAS_old& X);
  const TlDenseGeneralMatrix_BLAS_old& dotInPlace(
      const TlDenseSymmetricMatrix_BLAS_Old& X);

  virtual bool inverse();

  /// 全要素の和を返す
  virtual double sum() const;

  /// 対角要素の和を返す
  ///
  /// @return 対角要素の和
  virtual double trace() const;

  /// 転置行列にする
  ///
  /// @return 転置行列となったこのオブジェクト
  virtual const TlDenseGeneralMatrix_BLAS_old& transposeInPlace();

  /// 要素の絶対値の最大値を返す
  ///
  /// @param[out] outRow 該当要素の行
  /// @param[out] outCol 該当要素の列
  /// @return 要素の絶対値の最大値
  virtual double getMaxAbsoluteElement(
      TlMatrixObject::index_type* pOutRow = NULL,
      TlMatrixObject::index_type* pOutCol = NULL) const;

  /// calc RMS
  virtual double getRMS() const;

  // ---------------------------------------------------------------------------
  // I/O
  // ---------------------------------------------------------------------------
 public:
  /// 指定されたファイルパス名から行列要素を読み込む
  ///
  /// パスの区切り文字(UNIXなら'/'、windowsなら'\')は、
  /// この関数内では変換しない。
  /// 内部格納形式は、いわゆるRSFDとなる。
  /// @param[in] sFilePath ファイルパス名
  /// @retval true 読み取り成功
  /// @retval false 読み取り失敗
  virtual bool load(const std::string& filePath);

  /// 指定されたファイルパス名に行列要素を書き込む
  ///
  /// パスの区切り文字(UNIXなら'/'、windowsなら'\')は、
  /// この関数内では変換しない。
  /// 内部格納形式は、いわゆるRSFDとなる。
  /// @param[in] sFilePath ファイルパス名
  /// @retval true 書き込み成功
  /// @retval false 書き込み失敗
  virtual bool save(const std::string& sFilePath) const;
  // virtual bool save(const std::string& filePath) const =0;

  virtual bool saveText(const std::string& sFilePath) const;
  virtual bool saveText(std::ostream& os) const;

#ifdef HAVE_HDF5
  virtual bool saveHdf5(const std::string& filepath,
                        const std::string& h5path) const;
  virtual bool loadHdf5(const std::string& filepath, const std::string& h5path);
#endif  // HAVE_HDF5

  virtual TlSerializeData getSerialize() const;

 public:
  /// 対角成分をベクトルとして抽出する
  virtual TlVector_BLAS getDiagonalElements() const;

  /// Ax=bを求める
  TlDenseGeneralMatrix_BLAS_old solveLinearLeastSquaresProblem(
      const TlDenseGeneralMatrix_BLAS_old& inB) const;

 public:
  friend std::ostream& operator<<(std::ostream& out,
                                  const TlDenseGeneralMatrix_BLAS_old& rhs);
  template <typename T>
  void print(T& out) const;
  virtual std::string getCsv() const;

 protected:
  /// サブクラス用コンストラクタ
  /// このコンストラクタを用いた場合は、行列要素のメモリ確保を
  /// このクラスでは行いません。サブクラスで行う必要があります。
  TlDenseGeneralMatrix_BLAS_old(int row, int col, double* pData);

  /// 行列要素用メモリを確保する
  ///
  /// @param[in] isZeroClear ゼロで初期化する場合はtrueを指定する
  virtual void initialize(bool isZeroClear = true);

  /// 必要な行列要素の数を返す。
  /// initialize()などから呼び出される。
  virtual size_type getNumOfElements() const {
    return (this->getNumOfRows() * this->getNumOfCols());
  }

  /// オブジェクトの内容を破棄する
  virtual void clear();

  virtual size_type index(TlMatrixObject::index_type row,
                          TlMatrixObject::index_type col) const;

 protected:
  // variables
  static const size_type MAX_LOOP;

  double* data_;  /// 行列要素

#ifdef HAVE_LAPACK
 protected:
  friend TlVector_BLAS multiplicationByLapack(
      const TlDenseGeneralMatrix_BLAS_old& A, const TlVector_BLAS& X);
  friend TlVector_BLAS multiplicationByLapack(
      const TlVector_BLAS& X, const TlDenseGeneralMatrix_BLAS_old& A);

  /// LAPACKを使用した行列積(X x Y)を計算する
  ///
  /// X の列数とY の行数が一致しないといけない。
  /// @param[in] X 乗算の左側
  /// @param[in] Y 乗算の右側
  /// @return 行列積
  friend TlDenseGeneralMatrix_BLAS_old multiplicationByLapack(
      const TlDenseGeneralMatrix_BLAS_old& X, const TlDenseGeneralMatrix_BLAS_old& Y);

  /// 対称行列の積を求める
  friend TlDenseGeneralMatrix_BLAS_old multiplicationByLapack(
      const TlDenseSymmetricMatrix_BLAS_Old& X, const TlDenseGeneralMatrix_BLAS_old& Y);
  friend TlDenseGeneralMatrix_BLAS_old multiplicationByLapack(
      const TlDenseGeneralMatrix_BLAS_old& X, const TlDenseSymmetricMatrix_BLAS_Old& Y);

  /// 対称行列の固有値を求める
  /// @param[in] inMatrix 対称行列
  /// @param[out] outEigvVal 固有値が格納されたベクトル
  /// @param[out] outEigVec 固有値ベクトルが格納された行列
  /// @retval true 固有値が求められた
  /// @retval false エラーが発生した
  friend bool diagonalByLapack(const TlDenseSymmetricMatrix_BLAS_Old& inMatrix,
                               TlVector_BLAS* outEigVal,
                               TlDenseGeneralMatrix_BLAS_old* outEigVec);

  friend bool inverseByLapack(TlDenseGeneralMatrix_BLAS_old& inoutMatrix);

  /// the minimum norm solution to a real linear least
  /// squares problem:
  ///  Minimize 2-norm(| b - A*x |).
  friend bool solveLinearLeastSquaresProblemByLapack(
      const TlDenseGeneralMatrix_BLAS_old& inA,
      const TlDenseGeneralMatrix_BLAS_old& inB, TlDenseGeneralMatrix_BLAS_old* pX);
#else
// cause compile error
#error NOT found algebra package: need BLAS/LAPACK library
#endif  // BLAS_FOUND
};

////////////////////////////////////////////////////////////////////////
// inline functions
// inline double TlDenseGeneralMatrix_BLAS_old::operator()(
//     const TlMatrixObject::index_type nRow,
//     const TlMatrixObject::index_type nCol) const {
//   return this->get(nRow, nCol);
// }
//
// inline double& TlDenseGeneralMatrix_BLAS_old::operator()(
//     const TlMatrixObject::index_type nRow,
//     const TlMatrixObject::index_type nCol) {
//   //   const int index = nRow + (this->m_nRows * nCol);
//   return this->data_[this->index(nRow, nCol)];
// }

template <typename T>
void TlDenseGeneralMatrix_BLAS_old::print(T& out) const {
  const int nNumOfRows = this->getNumOfRows();
  const int nNumOfCols = this->getNumOfCols();

  for (int ord = 0; ord < nNumOfCols; ord += 10) {
    out << "       ";
    for (int j = ord; ((j < ord + 10) && (j < nNumOfCols)); ++j) {
      out << TlUtils::format("   %5d th", j + 1);
    }
    out << "\n ----";

    for (int j = ord; ((j < ord + 10) && (j < nNumOfCols)); ++j) {
      out << "-----------";
    }
    out << "----\n";

    for (int i = 0; i < nNumOfRows; ++i) {
      out << TlUtils::format(" %5d  ", i + 1);

      for (int j = ord; ((j < ord + 10) && (j < nNumOfCols)); ++j) {
        out << TlUtils::format(" %10.6lf", this->get(i, j));
      }
      out << "\n";
    }
    out << "\n\n";
  }
  out.flush();
}

#endif  // TL_DENSE_GENERAL_MATRIX_BLAS_H
