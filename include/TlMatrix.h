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

#ifndef TLMATRIX_H
#define TLMATRIX_H

#ifdef HAVE_CONFIG_H
#include "config.h"    // this file created by autotools
#endif // HAVE_CONFIG_H

#include <cassert>
#include <fstream>
#include <sstream>
#include <valarray>

#include "TlMatrixObject.h"
#include "TlUtils.h"
#include "TlSerializeData.h"

/// 密行列クラス
/// 最低でも32bit環境を想定している。
/// 要素の指定はlongではなく、intを用いている。
/// int(32bit)の範囲は-2,147,483,648〜2,147,483,647であり、
/// 1,000残基規模(約10万軌道)の計算には十分耐えられる。
class TlMatrix : public TlMatrixObject {
    friend class TlCommunicate;
    friend class TlVector;

public:
    /// 今後作成される行列オブジェクトに対し、
    /// メモリ管理オブジェクトとしてTlMemManagerを使用する(true)か否か(false)を決定する。
    ///
    /// @see isUsingMemManager_
    /// @see isUsingMemManagerDefault_
    static void useMemManager(bool isUsingMemManager);

    /// 行列オブジェクトを作成する
    ///
    /// 行列の要素は0で初期化される。
    /// @param[in] row 作成する行列の行数
    /// @param[in] col 作成する行列の列数
    explicit TlMatrix(index_type row = 1, index_type col = 1);

    TlMatrix(const TlSerializeData& data);
    
    /// 行列オブジェクトを作成する
    ///
    /// 行列の要素はベクトルから作成される。
    /// ベクトルサイズはrowとcolの積でなければならない。
    /// @param[in] rVector ベクトル
    /// @param[in] row 作成する行列の行数
    /// @param[in] col 作成する行列の列数
    TlMatrix(const TlVector& rVector,
             index_type row, index_type col);
    
public:
    /// コピーコンストラクタ
    TlMatrix(const TlMatrix& rhs);

    /// 対称行列オブジェクトからTlMatrixオブジェクトを作成する

    TlMatrix(const TlSymmetricMatrix& rhs);

    virtual ~TlMatrix();

public:
    /// 行数を返す
    ///
    /// @return 行数
    virtual index_type getNumOfRows() const;

    /// 列数を返す
    ///
    /// @return 列数
    virtual index_type getNumOfCols() const;

    /// 行列のサイズを変更する
    ///
    /// 行列を大きくする場合、追加される要素は0で初期化される。
    /// 行列を小さくする場合、切り詰められる要素は破棄される。
    /// @param[in] row 行数
    /// @param[in] col 列数
    void resize(index_type nRow, index_type nCol);

    /// 行列要素をベクトルにして返す
    ///
    /// 3x3 の行列であれば以下の順に返す。
    /// [ 0 3 6 ]
    /// [ 1 4 7 ]
    /// [ 2 5 8 ]
    virtual TlVector getVector() const;

    /// 指定した行の要素から構成されるベクトルを返す
    ///
    /// @param[in] nRow 指定する行
    //virtual TlVector getRowVector(int nRow) const;

    /// 指定した列の要素から構成されるベクトルを返す
    ///
    /// @param[in] nCol 指定する列
    //virtual TlVector getColumnVector(int nCol) const;

    // ===================================================================
    /// ブロック行列を返す
    ///
    /// @param[in] row 始点となる行
    /// @param[in] col 始点となる列
    /// @param[in] row_distance 取得する行数
    /// @param[in] col_distance 取得する列数
    /// @return row_distance × col_distance 次元のブロック行列
    virtual TlMatrix getBlockMatrix(index_type nRow, index_type nCol,
                                    index_type nRowDistance, index_type nColDistance) const;

    /// 行列要素を指定された位置に上書きする
    ///
    /// @param[in] row 始点となる行
    /// @param[in] col 始点となる列
    /// @param[in] matrix 行列要素
    virtual void setBlockMatrix(index_type nRow, index_type nCol,
                                const TlMatrix& matrix);

    /// 行列要素を指定された位置に加算する
    ///
    /// @param[in] row 始点となる行
    /// @param[in] col 始点となる列
    /// @param[in] matrix 行列要素
    virtual void addBlockMatrix(int nRow, int nCol, const TlMatrix& matrix);

    /// 対角要素の和を返す
    ///
    /// @return 対角要素の和
    virtual double trace() const;

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
    virtual double operator()(int nRow, int nCol) const;

    virtual double get(const index_type row,
                       const index_type col) const;

    /// 要素を返す(代入可能)
    ///
    /// 内部では、行列要素を(2次元配列ではなく)
    /// 1次元配列として保持しているので、
    /// 他のメンバ関数内でもこのメンバ関数が呼ばれる。
    /// @param[in] row 行数
    /// @param[in] col 列数
    /// @return 要素
    virtual double& operator()(int nRow, int nCol);

    virtual void set(const index_type row,
                     const index_type col,
                     const double value);
    virtual void add(const index_type row,
                     const index_type col,
                     const double value);

    /// 代入演算子
    TlMatrix& operator =(const TlMatrix& rhs);

    /// 代入演算子
    /// @param[in] rhs 対称行列
    TlMatrix& operator =(const TlSymmetricMatrix& rhs);

    TlMatrix& operator+=(const TlMatrix& rhs);
    TlMatrix& operator+=(const TlSymmetricMatrix& rhs);
    TlMatrix& operator-=(const TlMatrix& rhs);
    TlMatrix& operator-=(const TlSymmetricMatrix& rhs);
    TlMatrix& operator*=(const TlMatrix& rhs);
    TlMatrix& operator*=(const TlSymmetricMatrix& rhs);

    /// 行列を定数倍する
    ///
    /// @param[in] rhs 定数倍の値
    /// @return 計算後のオブジェクト
    TlMatrix& operator*=(const double& rhs);

    /// 行列を定数で割る
    ///
    /// @param[in] rhs 割る定数の値
    /// @return 計算後のオブジェクト
    TlMatrix& operator/=(const double& rhs);

    friend TlMatrix operator+(const TlMatrix& X, const TlMatrix& Y);
    friend TlMatrix operator-(const TlMatrix& X, const TlMatrix& Y);
    friend TlMatrix operator*(const TlMatrix& X, const TlMatrix& Y);
    friend TlMatrix operator*(const TlMatrix& X, double Y);
    friend TlMatrix operator*(double X, const TlMatrix& Y) {
        return (Y * X);
    };

    friend TlVector operator*(const TlMatrix& X, const TlVector& Y);
    friend TlVector operator*(const TlVector& Y, const TlMatrix& X);

    /// Hadamard product
    const TlMatrix& dot(const TlMatrix& X);
    const TlMatrix& dot(const TlSymmetricMatrix& X);

    /// インスタンスのメモリサイズを返す
    virtual std::size_t getMemSize() const;
    
    /// 指定された入力ストリームが読み込み可能かどうかを返す
    ///
    /// @param[in,out] ifs 入力ファイルストリーム
    /// @retval true 読み取り可能
    /// @retval false 読み取り不可能
    static bool isLoadable(std::ifstream& ifs);
    static bool isLoadable(const std::string& rFilePath);

    /// ヘッダ情報を読み取る
    static bool getHeaderInfo(const std::string& filePath,
                              int* pType = NULL,
                              int* pNumOfRows = NULL, int* pNumOfCols = NULL);

    /// 指定されたファイルパス名から行列要素を読み込む
    ///
    /// パスの区切り文字(UNIXなら'/'、windowsなら'\')は、
    /// この関数内では変換しない。
    /// 内部格納形式は、いわゆるRSFDとなる。
    /// @param[in] sFilePath ファイルパス名
    /// @retval true 読み取り成功
    /// @retval false 読み取り失敗
    virtual bool load(const std::string& sFilePath);

    /// 入力ストリームから行列要素を読み込む
    ///
    /// パスの区切り文字(UNIXなら'/'、windowsなら'\')は、
    /// この関数内では変換しない。
    /// 内部格納形式は、いわゆるRSFDとなる。
    /// @param[in,out] ifs 入力ストリーム
    /// @retval true 読み取り成功
    /// @retval false 読み取り失敗
    virtual bool load(std::ifstream& ifs);

    /// 指定されたファイルパス名に行列要素を書き込む
    ///
    /// パスの区切り文字(UNIXなら'/'、windowsなら'\')は、
    /// この関数内では変換しない。
    /// 内部格納形式は、いわゆるRSFDとなる。
    /// @param[in] sFilePath ファイルパス名
    /// @retval true 書き込み成功
    /// @retval false 書き込み失敗
    virtual bool save(const std::string& sFilePath) const;

    /// 出力ストリームに行列要素を書き込む
    ///
    /// パスの区切り文字(UNIXなら'/'、windowsなら'\')は、
    /// この関数内では変換しない。
    /// 内部格納形式は、いわゆるRSFDとなる。
    /// @param[in,out] ofs 出力ストリーム
    /// @retval true 書き込み成功
    /// @retval false 書き込み失敗
    virtual bool save(std::ofstream& ofs) const;

    virtual bool saveText(const std::string& sFilePath) const;
    virtual bool saveText(std::ostream& os) const;

    virtual TlSerializeData getSerialize() const;
   
    /// 転置行列にする
    ///
    /// @return 転置行列となったこのオブジェクト
    virtual const TlMatrix& transpose();

    /// 要素の絶対値の最大値を返す
    ///
    /// @param[out] outRow 該当要素の行
    /// @param[out] outCol 該当要素の列
    /// @return 要素の絶対値の最大値
    virtual double getMaxAbsoluteElement(int* pOutRow =NULL, int* pOutCol =NULL) const;

    virtual bool inverse();

    /// 対角成分をベクトルとして抽出する
    virtual TlVector getDiagonalElements() const;

    /// Ax=bを求める
    TlMatrix solveLinearLeastSquaresProblem(const TlMatrix& inB) const;

public:
    friend std::ostream& operator <<(std::ostream& out, const TlMatrix& rhs);
    template <typename T> void print(T& out) const;
    virtual std::string getCsv() const;

protected:
    /// サブクラス用コンストラクタ
    /// このコンストラクタを用いた場合は、行列要素のメモリ確保を
    /// このクラスでは行いません。サブクラスで行う必要があります。
    TlMatrix(int row, int col, double* pData);

    /// 行列要素用メモリを確保する
    ///
    /// @param[in] isZeroClear ゼロで初期化する場合はtrueを指定する
    virtual void initialize(bool isZeroClear = true);

    void initialize_usingStandard(bool isZeroClear);
    void initialize_usingMemManager(bool isZeroClear);
    
    /// 必要な行列要素の数を返す。
    /// initialize()などから呼び出される。
    virtual std::size_t getNumOfElements() const {
        return (this->getNumOfRows() * this->getNumOfCols());
    }

    /// オブジェクトの内容を破棄する
    virtual void clear();

    void clear_usingStandard();
    void clear_usingMemManager();
    
    virtual std::size_t index(index_type row,
                              index_type col) const;

    static bool getHeaderInfo(std::ifstream& ifs, int* pType = NULL,
                              int* pNumOfRows = NULL, int* pNumOfCols = NULL);

protected:
    // variables
    static const size_type MAX_LOOP;

    int m_nRows;    /// 行数
    int m_nCols;    /// 列数
    double* data_;  /// 行列要素

    bool isUsingMemManager_; /// TlMemManagerオブジェクトを用いてメモリ管理を行う(true)/行わない(false)

    /// TlMemManagerオブジェクトを使うか否か、デフォルトを規定する
    /// @see useMemManager()
    static bool isUsingMemManagerDefault_; 
    
#ifdef HAVE_LAPACK
protected:
    friend TlVector multiplicationByLapack(const TlMatrix& A, const TlVector& X);
    friend TlVector multiplicationByLapack(const TlVector& X, const TlMatrix& A);

    /// LAPACKを使用した行列積(X x Y)を計算する
    ///
    /// X の列数とY の行数が一致しないといけない。
    /// @param[in] X 乗算の左側
    /// @param[in] Y 乗算の右側
    /// @return 行列積
    friend TlMatrix multiplicationByLapack(const TlMatrix& X, const TlMatrix& Y);

    /// 対称行列の積を求める
    friend TlMatrix multiplicationByLapack(const TlSymmetricMatrix& X, const TlMatrix& Y);
    friend TlMatrix multiplicationByLapack(const TlMatrix& X, const TlSymmetricMatrix& Y);

    /// 対称行列の固有値を求める
    /// @param[in] inMatrix 対称行列
    /// @param[out] outEigvVal 固有値が格納されたベクトル
    /// @param[out] outEigVec 固有値ベクトルが格納された行列
    /// @retval true 固有値が求められた
    /// @retval false エラーが発生した
    friend bool diagonalByLapack(const TlSymmetricMatrix& inMatrix, TlVector* outEigVal, TlMatrix* outEigVec);

    friend bool inverseByLapack(TlMatrix& inoutMatrix);

    /// the minimum norm solution to a real linear least
    /// squares problem:
    ///  Minimize 2-norm(| b - A*x |).
    friend bool solveLinearLeastSquaresProblemByLapack(const TlMatrix& inA,
                                                       const TlMatrix& inB,
                                                       TlMatrix* pX);
#else
    // cause compile error
#error NOT found algebra package: need LAPACK library
#endif // HAVE_LAPACK
};


////////////////////////////////////////////////////////////////////////
// inline functions
inline double TlMatrix::operator()(const int nRow, const int nCol) const
{
    return this->get(nRow, nCol);
}


inline double& TlMatrix::operator()(const int nRow, const int nCol)
{
//   assert((0 <= nRow) && (nRow < this->m_nRows));
//   assert((0 <= nCol) && (nCol < this->m_nCols));

//   const int index = nRow + (this->m_nRows * nCol);
    return this->data_[this->index(nRow, nCol)];
}


template <typename T> void TlMatrix::print(T& out) const
{
    const int nNumOfRows = this->getNumOfRows();
    const int nNumOfCols = this->getNumOfCols();

    for (int ord = 0; ord < nNumOfCols; ord += 10) {
        out << "       ";
        for (int j = ord; ((j < ord+10) && (j < nNumOfCols)); ++j) {
            out << TlUtils::format("   %5d th", j+1);
        }
        out << "\n ----";

        for (int j = ord; ((j < ord+10) && (j < nNumOfCols)); ++j) {
            out << "-----------";
        }
        out << "----\n";

        for (int i = 0; i < nNumOfRows; ++i) {
            out << TlUtils::format(" %5d  ", i+1);

            for (int j = ord; ((j < ord+10) && (j < nNumOfCols)); ++j) {
                out << TlUtils::format(" %10.6lf", (*this)(i, j));
            }
            out << "\n";
        }
        out << "\n\n";
    }
    out.flush();
}

#endif // TLMATRIX_H

