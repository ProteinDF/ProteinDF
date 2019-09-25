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

#ifndef TLSTRINGTOKENIZER_H
#define TLSTRINGTOKENIZER_H

#include <iostream>
#include <string>

/// std::string オブジェクトをトークンに分割するクラス
class TlStringTokenizer {
   public:
    /// コンストラクタ
    ///
    /// @param[in] str 分割される文字列
    /// @param[in] delimiter デリミタとなる文字列.
    /// 空文字列の場合はホワイトスペースが使われる.
    explicit TlStringTokenizer(const std::string& str,
                               const std::string& delimiter = "");

    /// デストラクタ
    ~TlStringTokenizer();

   public:
    /// トークンの数を返す
    ///
    /// @return トークンの数
    std::size_t countTokens();

    /// 次のトークンの有無を返す
    ///
    /// @retval true 次のトークンあり
    /// @retval false 次のトークンなし
    bool hasMoreTokens();

    /// 次のトークンを返す
    ///
    /// @return 次のトークン. 無い場合は空文字列を返す.
    std::string nextToken();

   private:
    std::string m_delimiter;  /// デリミタを保持する
    std::string m_str;        /// 分割される文字列
    int m_count;              /// トークンの数
    std::string::size_type m_begin;  /// トークン処理を行う始めの文字カウント
    std::string::size_type m_end;  /// トークン処理を行う最後の文字カウント
};

#endif  // TLSTRINGTOKENIZER_H
