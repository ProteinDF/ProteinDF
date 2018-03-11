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

#include <cctype>
#include <cstdlib>
#include <iostream>
#include <sstream>

#include <stdarg.h>
#include <stdio.h>
#include <string.h>

#include "TlUtils.h"

std::string TlUtils::format(const char* psFormat, ...) {
  va_list ap;
  va_start(ap, psFormat);

  size_t nLength = 256;
  char* pBuf = NULL;
  int nResult;

  do {
    va_list aq;
    // va_copy(aq, ap); // !! this function is defined by C99 !!
    memcpy(&aq, &ap, sizeof(va_list));

    nLength *= 2;
    char* pNewBuf = new char[nLength];

    if (pBuf != NULL) {
      delete[] pBuf;
    }
    pBuf = pNewBuf;
    pNewBuf = NULL;

    nResult = vsnprintf(pBuf, nLength, psFormat, aq);
    // nResult  = snprintf(pBuf, nLength, psFormat, aq);
    if (nResult < 0) {
      std::abort();
    }

    va_end(aq);
  } while (static_cast<size_t>(nResult) >= nLength);

  va_end(ap);
  std::string str(pBuf);

  delete[] pBuf;

  return str;
}

template <>
std::string TlUtils::xtos<bool>(const bool& t) {
  std::stringstream s;
  s << std::boolalpha << t;
  return s.str();
}

void TlUtils::trim_ws(std::string& s) {
  if (s.empty()) {
    return;
  }

  std::string::iterator p;
  for (p = s.begin(); p != s.end() && std::isspace(*p); ++p) {
  }

  if (p != s.begin()) {
    s.erase(s.begin(), p);
  }
}

void TlUtils::trim_ws(std::wstring& s) {
  if (s.empty()) {
    return;
  }

  std::wstring::iterator p;
  for (p = s.begin(); p != s.end() && std::iswspace(*++p);) {
  }

  if (!std::iswspace(*p)) {
    p++;
  }

  s.erase(s.begin(), p);
}

void TlUtils::rtrim_ws(std::string& s) {
  if (s.empty()) {
    return;
  }

  std::string::iterator p;
  for (p = s.end(); p != s.begin() && std::isspace(*--p);) {
  }

  if (!std::iswspace(*p)) {
    p++;
  }

  if (p != s.end()) {
    s.erase(p, s.end());
  }
}

void TlUtils::rtrim_ws(std::wstring& s) {
  if (s.empty()) {
    return;
  }

  std::wstring::iterator p;
  for (p = s.end(); p != s.begin() && std::iswspace(*--p);) {
  }

  if (!std::iswspace(*p)) {
    p++;
  }

  if (p != s.end()) {
    s.erase(p, s.end());
  }
}

/**
 *  XMLのテキストのうち、余分な改行などを取り除く
 */
// std::string TlUtils::trimXmlText(const std::string& str){
//   std::istringstream in(str);
//   char c;
//   std::string sWhiteSpace = " \t\n";
//   std::string sAnswer = "";

//   // remove head white space
//   while (in.get(c)){
//     if (sWhiteSpace.find(c) == std::string::npos){
//       sAnswer += c;
//       break;;
//     }
//   }
//   // remove tail white space
//   std::string sTempWhiteSpace = "";
//   while (in.get(c)){
//     if (sWhiteSpace.find(c) == std::string::npos){
//       sAnswer += sTempWhiteSpace + c;
//       sTempWhiteSpace = "";
//     } else {
//       sTempWhiteSpace += c;
//     }
//   }

//   return sAnswer;
// }

/**
 *  文字列をトークンによって分割
 *  @param str [in] 変換文字列
 *  @param token [in] トークン
 *  @return 分割した文字列のリスト
 */
std::vector<std::string> TlUtils::split(const std::string& str,
                                        const std::string& token) {
  std::vector<std::string> aList;
  aList.clear();

  std::istringstream in(str);
  char c;

  std::string temp_str = "";
  while (in.get(c)) {
    if (token.find(c) == std::string::npos) {
      temp_str += c;
    } else {
      aList.push_back(temp_str);
      temp_str = "";
    }
  }
  if (!temp_str.empty()) {
    aList.push_back(temp_str);
  }

  return aList;
}

std::string& TlUtils::replace(std::string& str, const std::string& sFrom,
                              const std::string& sTo) {
  std::string::size_type n = 0;
  std::string::size_type next = 0;
  while ((n = str.find(sFrom, next)) != std::string::npos) {
    str.replace(n, sFrom.size(), sTo);
    next = n + sTo.size();
  }

  return str;
}

std::string TlUtils::toUpper(const std::string& str) {
  std::ostringstream out;
  std::istringstream in(str);
  char c;

  while (in.get(c)) {
    char tc = std::toupper(c);
    out << tc;
  }

  return out.str();
}

std::string TlUtils::toLower(const std::string& str) {
  std::ostringstream out;
  std::istringstream in(str);
  char c;

  while (in.get(c)) {
    char tc = std::tolower(c);
    out << tc;
  }

  return out.str();
}

// 文字列strをtimes回繰り返した文字列を返す
std::string TlUtils::repeat(const std::string& str, const int times) {
  std::ostringstream out;

  for (int i = 0; i < times; ++i) {
    out << str;
  }

  return out.str();
}

std::string TlUtils::getWord(std::string& str) {
  TlUtils::trim_ws(str);
  std::string sAnswer = "";

  std::string::size_type nPosition = str.find_first_of(" \f\n\r\t\v");

  if (nPosition != std::string::npos) {
    sAnswer = str.substr(0, nPosition);
    str = str.substr(nPosition);
  } else {
    sAnswer = str;
    str = "";
  }

  return sAnswer;
}

std::string TlUtils::textWrap(const std::string& str, size_t width) {
  std::string tmp = "";
  char cur = '\0';
  char last = '\0';
  size_t i = 0;

  std::istringstream in(str);
  std::ostringstream out;
  while (in.get(cur)) {
    if (++i == width) {
      TlUtils::trim_ws(tmp);
      out << '\n' << tmp;
      tmp.clear();
      i = 0;
    } else if ((std::isspace(cur) == true) &&
               (!std::isspace(last))) {  // 単語の終わり
      out << tmp;
      tmp.clear();
    }

    if (cur == '\n') {
      // flush
      out << tmp << '\n';
      tmp.clear();
      i = 0;
      last = '\0';
    } else {
      tmp += cur;
      last = cur;
    }
  }
  out << tmp;

  return out.str();
}

/**
 *  ホワイトスペースもしくは"[]"で囲まれた文字列を返す。
 *  入力文字列は該当部分が除去される。
 *
 *  例) "[hoge] [foo]" => 返値: hoge, 入力: " [foo]"
 *      "hoge foo"     => 返値: hoge, 入力: " foo"
 *
 * @param str [in/out] 入力文字列。該当文字列は除去される。
 * @return ホワイトスペースもしくは"[]"で囲まれた文字列
 */
std::string TlUtils::getPdfParam(std::string& str) {
  // 先頭の空白を除去
  TlUtils::trim_ws(str);

  std::string sAnswer = "";

  bool bBraketMode = false;
  if (str[0] == '[') {
    std::string::size_type nBraket = str.find(']');

    if (nBraket != std::string::npos) {
      // [] で囲まれた領域を返す
      bBraketMode = true;
      sAnswer = str.substr(1, nBraket - 1);
      str = str.substr(nBraket + 1);
    }
  }

  if (bBraketMode == false) {
    std::string::size_type nWhiteSpace = str.find_first_of(" \t");
    if (nWhiteSpace != std::string::npos) {
      // ホワイトスペースで囲まれた領域を返す

      // 注:本来なら、以下のコードで足りる。
      // >>>> code start:
      sAnswer = str.substr(0, nWhiteSpace);
      // <<<< code end:
      //
      // しかし、このコードでは何故か並列時に問題を生じる。
      // 従って、以下のコードで代用する。
      //       std::string tmp = "";
      //       std::cout << "ws = " << nWhiteSpace << std::endl;
      //       for (std::string::size_type i=0; i<nWhiteSpace; ++i){
      // 	std::cout << i << ": " << str[i] << std::endl;
      //   	tmp.append(&str[i]);
      //       }
      //       std::cout << sAnswer << "=?" << tmp << "." << std::endl;
      //       {
      // 	int nSize = nWhiteSpace +1; // +1 for '\0'
      // 	char* pBuf = new char[nSize];
      // 	str.copy(pBuf, nSize -1);
      // 	pBuf[nSize -1] = '\0';

      // 	sAnswer = std::string(pBuf);
      // 	delete[] pBuf;
      // 	pBuf = NULL;
      //       }

      str = str.substr(nWhiteSpace + 1);
    } else {
      // 文字列全体を返す
      sAnswer = str;
      str = "";
    }
  }

  return sAnswer;
}

std::string TlUtils::getPdfSlash(std::string& str) {
  std::string sAnswer = "";
  TlUtils::trim_ws(str);

  bool bFoundLeftBraket = false;
  while (str.length() > 0) {
    if ((str[0] == '/') && (bFoundLeftBraket == false)) {
      str = str.substr(1);  // 先頭の"/"を除去
      break;
    }

    if (str[0] == '{') {
      bFoundLeftBraket = true;
    } else if ((bFoundLeftBraket == true) && (str[0] == '}')) {
      bFoundLeftBraket = false;
    }

    sAnswer += str[0];
    if (str.length() > 1) {
      str = str.substr(1);
    } else {
      str = "";
    }
  }

  return sAnswer;
}

void TlUtils::progressbar(const float progress) {
  static const int barWidth = 70;

  std::cout << "[";
  const int pos = barWidth * progress;
  for (int i = 0; i < barWidth; ++i) {
    if (i < pos) {
      std::cout << "=";
    } else if (i == pos) {
      std::cout << ">";
    } else {
      std::cout << " ";
    }
  }
  std::cout << "] " << std::min<int>(int(progress * 100.0), 100) << " %\r";
  std::cout.flush();
}
