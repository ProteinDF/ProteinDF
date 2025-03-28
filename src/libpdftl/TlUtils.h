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

#ifndef TLUTILS_H
#define TLUTILS_H

#include <stdarg.h>

#include <algorithm>
#include <cctype>
#include <cwctype>
#include <set>
#include <sstream>
#include <string>
#include <typeinfo>
#include <vector>

/// 文字列処理クラス
class TlUtils {
public:
    /// 安全な sprintf() (std::string 版)
    static std::string format(const char* psFormat, ...);

    /** class Tをstringに変換
     *
     */
    template <typename T>
    static std::string xtos(const T& t);

    /** 文字列のpadding
     *
     *  指定された文字s が長さn になるまで、文字c でpaddingする
     */
    template <typename T>
    static void pad(std::basic_string<T>& s,
                    typename std::basic_string<T>::size_type n, T c);

    /// 先頭から同一の文字を削除する
    template <typename T>
    static void trim(std::basic_string<T>& s, T c);

    /// 末尾から同一の文字を削除する
    template <typename T>
    static void rtrim(std::basic_string<T>& s, T c);

    /** 先頭の空白を除去する
     *
     *  @param[in,out] str 文字列
     */
    static void trim_ws(std::string& str);

    /** 先頭の空白を除去
     *
     *  @param[in,out] str 文字列
     */
    static void trim_ws(std::wstring& str);

    static void rtrim_ws(std::string& str);
    static void rtrim_ws(std::wstring& str);

    /** 文字列を区切り文字で区切る
     *
     */
    static std::vector<std::string> split(const std::string& str, const std::string& delimiter);

    static std::string join(const std::vector<std::string>& strings, const std::string& delimiter);

    /** 文字列を置換する
     *
     */
    static std::string& replace(std::string& str, const std::string& sFrom, const std::string& sTo);

    /// 大文字に変換する
    static std::string toUpper(const std::string& str);

    /// 小文字に変換する
    static std::string toLower(const std::string& str);

    /// 文字列の繰り返し
    ///
    /// 文字列strをtimes回繰り返した文字列を返す
    /// @param str 繰り返す文字列
    /// @param times 回数
    /// @return 指定回数分繰り返した文字列
    static std::string repeat(const std::string& str, int times);

    /**
     * ホワイトスペースで挟まれた文字列を返す
     */
    static std::string getWord(std::string& str);

    /** テキストを折り返す
     *
     */
    static std::string textWrap(const std::string& str, size_t width);

    /**
     *  ホワイトスペースもしくは"[]"で囲まれた文字列を返す。
     *  入力文字列は該当部分が除去される。
     *
     * @param str [in/out] 入力文字列。該当文字列は除去される。
     * @return ホワイトスペースもしくは"[]"で囲まれた文字列
     */
    static std::string getPdfParam(std::string& str);

    /**
     * "/" で区切られている最初の文字列を返す
     * 抜き出す文字がない場合は、ヌル文字列を返す。
     * {}で括られている部分は／で区切られていても一つのシンタックスと認識する。
     * 抜き出した文字とその次の "/" は切り詰められる。
     *
     * @param str [in/out] 入力文字列。該当文字列は除去される。
     * @return "/" で区切られている最初の文字列
     */
    static std::string getPdfSlash(std::string& str);

    //   // 末尾から同一の文字を削除する
    //   template<typename T>
    //   void rTrim(std::basic_string<T>& s, T c);

    //   void rTrimWs(std::string& s){
    //     TlUtils::rTrimWs(s, std::isspace);
    //   }

    //   void rTrimWs(std::wstring& s){
    //     TlUtils::rTrimWs(s, std::iswspace);
    //   }

    //   template<typename T, typename F>
    //   void rTrimWs(std::basic_string<T>& s, F f);

    template <typename T>
    static T changeEndian(T value) {
        char* p = reinterpret_cast<char*>(&value);
        std::reverse(p, p + sizeof(T));

        return value;
    }

    static void changeEndian(char* p, std::size_t size) {
        std::reverse(p, p + size);
    }

    static bool isLittleEndian() {
        int i = 1;
        return *(reinterpret_cast<char*>(&i)) == 1;
    }

    static bool isBigEndian() {
        return !(TlUtils::isLittleEndian());
    }

    template <typename T>
    static T toBigEndian(T value) {
        if (TlUtils::isLittleEndian() == true) {
            value = TlUtils::changeEndian(value);
        }
        return value;
    }

    template <typename T>
    static T toLittleEndian(T value) {
        if (TlUtils::isBigEndian() == true) {
            value = TlUtils::changeEndian(value);
        }
        return value;
    }

    template <typename T>
    static T bigEndianToLocal(T value) {
        if (TlUtils::isLittleEndian() == true) {
            value = TlUtils::changeEndian(value);
        }
        return value;
    }

    static void progressbar(const float progress);

    /// change memory layout
    ///
    /// a11, a12, a13, a21, a22, a23,...
    /// ->
    /// a11, a21, a31, ..., a12, a22, a32, ..., a13, a23, a33, ...
    template <typename T>
    static void changeMemoryLayout(const T* pInput, const std::size_t axis1,
                                   const std::size_t axis2, T* pOutput);

public:
    /// vector notation
    /// "[1-3, 5]" -> [1, 2, 3, 5]
    static std::vector<int> vector_notation(const std::string& inputStr);

    /// [1, 2, 3, 4, 5, 2, 3, 4] -> [1, 2, 3, 4, 5]
    template <typename T>
    static std::vector<T> nonreduntant_vector(const std::vector<T>& vtr);

public:
    template <typename T>
    static std::string vector2str(const std::vector<T>& vtr);

public:
    template <typename T>
    static std::string getTypeName(const T& t);

private:
    template <typename T, typename F>
    static void trim_ws_t(std::basic_string<T>& s, F f);
};

////////////////////////////////////////////////////////////////////////
// テンプレートの定義

template <typename T>
std::string TlUtils::xtos(const T& t) {
    std::stringstream s;
    s << t;
    return s.str();
}

template <typename T>
void TlUtils::pad(std::basic_string<T>& s,
                  typename std::basic_string<T>::size_type n, T c) {
    if (n > s.length()) {
        s.append(n - s.length(), c);
    }
}

template <typename T>
void TlUtils::trim(std::basic_string<T>& s, T c) {
    if (s.empty()) {
        return;
    }

    typename std::basic_string<T>::iterator p;
    for (p = s.begin(); p != s.end() && *p == c; ++p) {
    }

    if (p != s.begin()) {
        s.erase(s.begin(), p);
    }
}

// template<typename T, typename F>
// void TlUtils::trim_ws_t(std::basic_string<T>& s, F f){
//   if (s.empty()){
//     return;
//   }

//   typename std::basic_string<T>::iterator p;
//   for (p = s.begin(); p != s.end() && f(*++p); ){
//   }

//   if (!f(*p)){
//     p++;
//   }

//   s.erase(s.begin(), p);
// }

template <typename T>
void TlUtils::rtrim(std::basic_string<T>& s, T c) {
    if (s.empty()) {
        return;
    }

    typename std::basic_string<T>::iterator p;
    for (p = s.end(); p != s.begin() && *--p == c;) {
    }

    if (*p != c) {
        ++p;
    }

    s.erase(p, s.end());
}

// template<typename T, typename F>
// void rTrimWs(std::basic_string<T>& s, F f){
//   if (s.empty()){
//     return;
//   }

//   typename std::basic_string<T>::iterator p;
//   for (p = s.end(); p != s.begin() && f(*--p); ){
//   }

//   if (!f(*p)){
//     p++;
//   }

//   s.erase(p, s.end());
// }

template <typename T>
std::vector<T> TlUtils::nonreduntant_vector(const std::vector<T>& vtr) {
    std::set<T> s;
    {
        typename std::vector<T>::const_iterator itEnd = vtr.end();
        for (typename std::vector<T>::const_iterator it = vtr.begin();
             it != itEnd; ++it) {
            s.insert(*it);
        }
    }

    std::vector<T> v;
    v.reserve(s.size());
    {
        typename std::set<T>::const_iterator itEnd = s.end();
        for (typename std::set<T>::const_iterator it = s.begin(); it != itEnd;
             ++it) {
            v.push_back(*it);
        }
    }

    return v;
}

template <typename T>
std::string TlUtils::vector2str(const std::vector<T>& vtr) {
    std::string output = "";

    typename std::vector<T>::const_iterator itEnd = vtr.end();
    for (typename std::vector<T>::const_iterator it = vtr.begin(); it != itEnd;
         ++it) {
        output += TlUtils::xtos(*it) + ",";
    }

    return output;
}

template <typename T>
void TlUtils::changeMemoryLayout(const T* pInput, const std::size_t axis1,
                                 const std::size_t axis2, T* pOutput) {
    for (std::size_t i = 0; i < axis1; ++i) {
        for (std::size_t j = 0; j < axis2; ++j) {
            pOutput[j * axis1 + i] = pInput[axis2 * i + j];
        }
    }
}

template <typename T>
std::string TlUtils::getTypeName(const T& t) {
    return typeid(t).name();
}

#endif  // TLUTILS_H
