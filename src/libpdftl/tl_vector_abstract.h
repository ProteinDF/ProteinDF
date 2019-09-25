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

#ifndef TLVECTOROBJECT_H
#define TLVECTOROBJECT_H

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif  // HAVE_CONFIG_H

#include <fstream>

/** 行列クラスのインターフェースを規定する
 */
class TlVectorAbstract {
   public:
    /// 行数、列数を表す整数型(範囲外を示す値として-1を取る場合がある)
    typedef signed int index_type;

    /// 数値配列の総数を表す整数型
    typedef signed int size_type;

   public:
    virtual ~TlVectorAbstract() {}

   public:
    virtual void resize(index_type nSize) = 0;

   public:
    virtual void add(const index_type index, const double value) = 0;

    /// 要素を返す(読み取り専用)
    ///
    /// @param[in] index 要素番号
    /// @return 要素
    virtual double operator[](index_type index) const = 0;

    /// 要素を返す(代入用)
    ///
    /// @param[in] index 要素番号
    /// @return 要素
    virtual double& operator[](index_type index) = 0;

   public:
    virtual bool save(const std::string& str) const = 0;
    virtual bool load(const std::string& str) = 0;
};

#endif  // TLVECTOROBJECT_H
