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

#ifndef TLSYMMETRICMATRIX_H
#define TLSYMMETRICMATRIX_H

#ifdef HAVE_CONFIG_H
#include "config.h"  // this file created by autotools
#endif               // HAVE_CONFIG_H

#include <algorithm>
#include <fstream>
#include <iostream>
#include <vector>

#include "TlMatrix.h"
#include "TlSparseSymmetricMatrix.h"
#include "TlVector.h"

class TlCommunicate;  // 通信クラス

/// 対称密行列クラス
///
/// 行列は正方行列
/// 要素の格納方法はLapackでいうところの'U'
/// AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;
///
// Packed Storage
// Symmetric, Hermitian or triangular matrices may be stored more compactly,
// if the relevant triangle (again as specified by UPLO) is packed by columns in
// a one-dimensional array. In LAPACK, arrays that hold matrices in packed
// storage, have names ending in `P'. So:
//    if UPLO = `U', aij is stored in AP(i+j(j-1)/2) for $i \leq j$;
//    if UPLO = `L', aij is stored in AP( i+(2n-j)(j-1)/2) for  $j \leq i$.
class TlSymmetricMatrix : public TlMatrix {
  // friend
  friend class TlCommunicate;
  friend class TlVector;

 public:
  /// 正方対称行列オブジェクトを作成する
  ///
  /// 行列の要素は0で初期化される。
  /// @param[in] dim 作成する行列の次元数
  explicit TlSymmetricMatrix(int dim = 1);

  /// 正方対称行列オブジェクトを作成する
  ///
  /// 行列の要素はバッファから作成される。
  /// バッファサイズは dim + dim * (dim -1) /2 でなければならない。
  /// @param[in] dim 作成する行列の次元数
  TlSymmetricMatrix(const TlVector& rVector, int nDim);

  /// コピーコンストラクタ
  TlSymmetricMatrix(const TlSymmetricMatrix& rhs);

  /// 行列オブジェクトから正方対称行列オブジェクトを作成する
  ///
  /// 対称行列になっているかのチェックは行わず、
  /// 左下の要素を格納する。
  TlSymmetricMatrix(const TlMatrix& rhs);

  TlSymmetricMatrix(const TlSerializeData& data);

  virtual ~TlSymmetricMatrix();

 public:
  /** 行列のサイズを変更する
   *
   *  行列を大きくする場合、追加される要素は0で初期化される。
   *  行列を小さくする場合、切り詰められる要素は破棄される。
   *
   *  @param[in] dim 次元数
   */
  virtual void resize(index_type dim);

  /// 行列要素をベクトルにして返す
  ///
  /// 3x3 の行列であれば以下の順に返す。
  /// [ 0 - - ]
  /// [ 1 3 - ]
  /// [ 2 4 5 ]
  // virtual TlVector getTlVector() const;

  //   /// ベクトルから行列に変換する
  //   void convertFromVector(int dim, const TlVector& vector);

  // static std::size_t vtr_index(index_type row,
  //                              index_type col,
  //                              index_type dim);

  /// 全要素の和を返す
  virtual double sum() const;

  /// 要素を返す(読み取り専用)
  ///
  /// 内部では、行列要素を(2次元配列ではなく)
  /// 1次元配列として保持しているので、
  /// 他のメンバ関数内でもこのメンバ関数が呼ばれる。
  /// @param[in] row 行数
  /// @param[in] col 列数
  /// @return 要素
  virtual double operator()(int row, int col) const;
  virtual double get(int row, int col) const;

  /// 要素を返す(代入可能)
  ///
  /// 内部では、行列要素を(2次元配列ではなく)
  /// 1次元配列として保持しているので、
  /// 他のメンバ関数内でもこのメンバ関数が呼ばれる。
  /// @param[in] row 行数
  /// @param[in] col 列数
  /// @return 要素
  virtual double& operator()(int row, int col);

  // virtual void set(index_type row, index_type col, double value);
  // virtual void add(index_type row, index_type col, double value);

  /// 代入演算子
  TlSymmetricMatrix& operator=(const TlSymmetricMatrix& rhs);

  TlSymmetricMatrix& operator+=(const TlSymmetricMatrix& rhs);
  // virtual TlSymmetricMatrix& operator+=(const TlSparseSymmetricMatrix& rhs);

  // TlSymmetricMatrix& operator+=(const TlMatrix& rhs);
  TlSymmetricMatrix& operator-=(const TlSymmetricMatrix& rhs);

  /// 行列を定数倍する
  ///
  /// @param[in] rhs 定数倍の値
  /// @return 計算後のオブジェクト
  TlSymmetricMatrix& operator*=(const double& rhs);

  /// 行列を定数で割る
  /// @param[in] rhs 割る定数の値
  /// @return 計算後のオブジェクト
  TlSymmetricMatrix& operator/=(const double& rhs);

  friend TlSymmetricMatrix operator+(const TlSymmetricMatrix& X,
                                     const TlSymmetricMatrix& Y);
  friend TlMatrix operator+(const TlMatrix& X, const TlSymmetricMatrix& Y);
  friend TlMatrix operator+(const TlSymmetricMatrix& X, const TlMatrix& Y) {
    return (Y + X);
  }
  friend TlSymmetricMatrix operator-(const TlSymmetricMatrix& X,
                                     const TlSymmetricMatrix& Y);
  friend TlMatrix operator-(const TlMatrix& X, const TlSymmetricMatrix& Y);
  friend TlMatrix operator-(const TlSymmetricMatrix& X, const TlMatrix& Y);
  friend TlMatrix operator*(const TlSymmetricMatrix& X,
                            const TlSymmetricMatrix& Y);
  friend TlMatrix operator*(const TlMatrix& X, const TlSymmetricMatrix& Y);
  friend TlMatrix operator*(const TlSymmetricMatrix& X, const TlMatrix& Y);
  friend TlSymmetricMatrix operator*(const TlSymmetricMatrix& X, double Y);
  friend TlSymmetricMatrix operator*(double X, const TlSymmetricMatrix& Y);
  friend TlVector operator*(const TlSymmetricMatrix& A, const TlVector& X);
  friend TlVector operator*(const TlVector& X, const TlSymmetricMatrix& A);

  const TlSymmetricMatrix& dot(const TlSymmetricMatrix& X);

 public:
  /// 指定された入力ストリームが読み込み可能かどうかを返す
  ///
  /// @param[in] ifs 入力ファイルストリーム
  /// @retval true 読み取り可能
  /// @retval false 読み取り不可能
  static bool isLoadable(const std::string& rFilePath);
  static bool isLoadable(std::ifstream& ifs);

  /// ヘッダ情報を読み取る
  static bool getHeaderInfo(const std::string& filePath, int* pType = NULL,
                            index_type* pNumOfRows = NULL,
                            index_type* pNumOfCols = NULL);

  /// ヘッダ情報を読み取る
  /// 正常に読み取れた場合、ファイルポインタはデータ領域の先頭(ヘッダのすぐ後)を指す
  static bool getHeaderInfo(std::fstream& fs, int* pType = NULL,
                            index_type* pNumOfRows = NULL,
                            index_type* pNumOfCols = NULL);
  static bool getHeaderInfo(std::ifstream& ifs, int* pType = NULL,
                            index_type* pNumOfRows = NULL,
                            index_type* pNumOfCols = NULL);

 protected:
  template <typename StreamType>
  static bool getHeaderInfo_tmpl(StreamType& s, int* pMatrixType = NULL,
                                 index_type* pNumOfRows = NULL,
                                 index_type* pNumOfCols = NULL);

 public:
  /// 指定されたファイルパス名から行列要素を読み込む
  ///
  /// パスの区切り文字(UNIXなら'/'、windowsなら'\')は、
  /// この関数内では変換しない。
  /// 内部格納形式は、いわゆるRLHDとなる。
  /// @param[in] sFile ファイルパス名
  /// @retval true 読み取り成功
  /// @retval false 読み取り失敗
  virtual bool load(const std::string& sFile);

  /// 入力ストリームから行列要素を読み込む
  ///
  /// パスの区切り文字(UNIXなら'/'、windowsなら'\')は、
  /// この関数内では変換しない。
  /// 内部格納形式は、いわゆるRLHDとなる。
  /// @param[in,out] ifs 入力ストリーム
  /// @retval true 読み取り成功
  /// @retval false 読み取り失敗
  virtual bool load(std::ifstream& ifs);

  /// 指定されたファイルパス名に行列要素を書き込む
  ///
  /// パスの区切り文字(UNIXなら'/'、windowsなら'\')は、
  /// この関数内では変換しない。
  /// 内部格納形式は、いわゆるRLHDとなる。
  /// @param[in] sFile ファイルパス名
  /// @retval true 書き込み成功
  /// @retval false 書き込み失敗
  virtual bool save(const std::string& sFile) const;

  /// 出力ストリームに行列要素を書き込む
  ///
  ///  パスの区切り文字(UNIXなら'/'、windowsなら'\')は、
  /// この関数内では変換しない。
  /// 内部格納形式は、いわゆるRLHDとなる。
  /// @param[in,out] ofs 出力ストリーム
  /// @retval true 書き込み成功
  /// @retval false 書き込み失敗
  virtual bool save(std::ofstream& ofs) const;

  virtual TlSerializeData getSerialize() const;

#ifdef HAVE_HDF5
  virtual bool saveHdf5(const std::string& filepath,
                        const std::string& h5path) const;
  virtual bool loadHdf5(const std::string& filepath, const std::string& h5path);
#endif  // HAVE_HDF5

 public:
  /// 転置行列にする
  ///
  ///  転置行列となったこのオブジェクト
  virtual TlSymmetricMatrix& transpose();

  /// 要素の絶対値の最大値を返す
  ///
  /// @param[out] outRow 該当要素の行
  /// @param[out] outCol 該当要素の列
  /// @return 要素の絶対値の最大値
  virtual double getMaxAbsoluteElement(int* outRow = NULL,
                                       int* outCol = NULL) const;

  /// 各行の要素の絶対値の最大値をベクトルとして返す
  TlVector getMaxAbsoluteVectorOnEachRow() const;

  /// 指定された値の行または列の要素の絶対値の最大値を返す
  ///
  /// @return 要素の絶対値の最大値
  double getMaxAbsoluteElementByIndex(int index) const;

  /// calc RMS
  virtual double getRMS() const;

  /// 固有値を求める
  ///
  /// @param[out] pEigVal 固有値が格納されたベクトル
  /// @param[out] pEigVec 固有値ベクトルが格納された行列
  /// @retval true 固有値が求められた
  /// @retval false エラーが発生した
  virtual bool diagonal(TlVector* pEigVal, TlMatrix* pEigVec) const;

  virtual bool inverse();

  TlMatrix choleskyFactorization();
  TlMatrix choleskyFactorization2(const double threshold = 1.0E-16) const;
  TlMatrix choleskyFactorization2omp(const double threshold = 1.0E-16) const;

 public:
  friend std::ostream& operator<<(std::ostream& out,
                                  const TlSymmetricMatrix& rhs);
  template <typename T>
  void print(T& out) const;

 protected:
  /// サブクラス用コンストラクタ
  /// データ領域の確保(this->data_)はサブクラスに任せる
  // TlSymmetricMatrix(int dim, double* pData);

  /// 行列要素用メモリを確保する
  ///
  /// @param[in] isZeroClear ゼロで初期化する場合はtrueを指定する
  // virtual void initialize(bool isZeroClear = true);

  /// 必要な行列要素の数を返す。
  /// initialize()などから呼び出される。
  virtual size_type getNumOfElements() const;

  /// オブジェクトの内容を破棄する
  // virtual void clear();

  virtual size_type index(index_type row, index_type col) const;

  bool load_RLHD(std::ifstream& ifs);
  bool load_CLHD(std::ifstream& ifs);

#ifdef HAVE_LAPACK
  /// 対称行列の積を求める
  friend TlMatrix multiplicationByLapack(const TlSymmetricMatrix& X,
                                         const TlMatrix& Y);
  friend TlMatrix multiplicationByLapack(const TlMatrix& X,
                                         const TlSymmetricMatrix& Y);

  friend TlVector multiplicationByLapack(const TlSymmetricMatrix& A,
                                         const TlVector& X);

  /// 対称行列の固有値を求める
  ///
  /// @param[in] inMatrix 対称行列
  /// @param[out] outEigVal 固有値が格納されたベクトル
  /// @param[out] outEigVec 固有値ベクトルが格納された行列
  /// @retval true 固有値が求められた
  /// @retval false エラーが発生した
  friend bool diagonalByLapack(const TlSymmetricMatrix& inMatrix,
                               TlVector* outEigVal, TlMatrix* outEigVec);

  friend bool inverseByLapack(TlSymmetricMatrix& inoutMatrix);

  // friend int choleskyFactorization(TlSymmetricMatrix* A,
  //                                  std::vector<int>* pPivot);

#else
  // cause compile error
#error NOT found algebra package: need LAPACK library
#endif  // HAVE_LAPACK

 private:
  /// 行列のサイズを変更する
  ///
  /// 行列を大きくする場合、追加される要素は0で初期化される。
  /// 行列を小さくする場合、切り詰められる要素は破棄される。
  /// このメンバ関数は下位互換性のために用意されている。
  /// @param[in] row 行数
  /// @param[in] col 列数
  virtual void resize(index_type row, index_type col) {
    assert(row == col);
    this->resize(row);
  }
};

////////////////////////////////////////////////////////////////////////
// inline functions
//
inline double TlSymmetricMatrix::get(int row, int col) const {
  return (this->data_[this->index(row, col)]);
}

inline double TlSymmetricMatrix::operator()(int row, int col) const {
  return this->get(row, col);
}

inline double& TlSymmetricMatrix::operator()(int row, int col) {
  return (this->data_[this->index(row, col)]);
}

template <typename T>
void TlSymmetricMatrix::print(T& out) const {
  const int nNumOfDim = this->getNumOfRows();  // == this->getNumOfCols()

  out << "\n\n";
  for (int ord = 0; ord < nNumOfDim; ord += 10) {
    out << "       ";
    for (int j = ord; ((j < ord + 10) && (j < nNumOfDim)); ++j) {
      out << TlUtils::format("   %5d th", j + 1);
    }
    out << "\n"
        << " ----";

    for (int j = ord; ((j < ord + 10) && (j < nNumOfDim)); ++j) {
      out << "-----------";
    }
    out << "----\n";

    for (int i = 0; i < nNumOfDim; ++i) {
      out << TlUtils::format(" %5d  ", i + 1);

      for (int j = ord; ((j < ord + 10) && (j < nNumOfDim)); ++j) {
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

#endif  // TLSYMMETRICMATRIX_H
