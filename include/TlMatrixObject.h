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

#ifndef TLMATRIXOBJECT_H
#define TLMATRIXOBJECT_H

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif // HAVE_CONFIG_H

#include "TlVector.h"

/// MAX_INDEX_BITSが指定するビット数まで行列の次数を確保できる。
/// 例えばMAX_INDEX_BITS=20であれば、1,048,576まで行列が確保できる。
#define MAX_INDEX_BITS (20) 

/// 行列クラスのインターフェースを規定する
class TlMatrixObject {
public:
    typedef signed int index_type; /// 行数、列数を表す整数型(範囲外を示す値として-1を取る場合がある)
    typedef signed long size_type; /// 数値配列の総数を表す整数型(範囲外を示す値として-1を取る場合がある)
    
public:
    /// 仮想デストラクタ
    virtual ~TlMatrixObject() {
    }

public:
    virtual index_type getNumOfRows() const =0;
    virtual index_type getNumOfCols() const =0;

    /// インスタンスのメモリサイズを返す
    virtual std::size_t getMemSize() const =0;
    
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
    virtual double get(index_type row, index_type col) const =0;

    virtual double getLocal(index_type row, index_type col) const {
        return this->get(row, col);
    }
    
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
    virtual void set(index_type row, index_type col, double value) =0;

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
    virtual void add(index_type row, index_type col, double value) =0;

    /// 指定した行の要素から構成されるベクトルを返す
    ///
    /// @param[in] row 指定する行
    virtual TlVector getRowVector(index_type row) const;

    /// 指定した列の要素から構成されるベクトルを返す
    ///
    /// @param[in] col 指定する行
    virtual TlVector getColVector(index_type col) const;

    /// リストからまとめて要素を加算する
    ///
    /// @param[in] pIndexPairs インデックスのリスト(row1, col1, row2, col2, ...)
    /// @param[in] pValues 値のリスト
    /// @param[in] size 要素数
    virtual void addByList(const index_type* pIndexPairs,
                           const double* pValues,
                           const std::size_t size);

public:
    //virtual TlMatrixObject& operator+=(const TlMatrixObject& rhs);

public:
    /// 指定されたパスから内容を読み込む
    ///
    /// @retval true 成功
    /// @retval false 失敗
    virtual bool load(const std::string& path) =0;

    /// 指定されたパスに内容を書き出す
    ///
    /// @retval true 成功
    /// @retval false 失敗
    virtual bool save(const std::string& path) const =0;
    
public:
    virtual double getMaxAbsoluteElement(index_type* pOutRow =NULL, index_type* pOutCol =NULL) const;
};

#endif // TLMATRIXOBJECT_H

