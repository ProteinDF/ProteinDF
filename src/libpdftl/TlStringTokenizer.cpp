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

#include "TlStringTokenizer.h"

TlStringTokenizer::TlStringTokenizer(const std::string& str,
                                     const std::string& delimiter)
    : m_str(str), m_count(-1), m_begin(0), m_end(0) {
  if (delimiter == "") {
    this->m_delimiter = " \f\n\r\t\v";  // デフォルトはホワイトスペース
  } else {
    this->m_delimiter = delimiter;
  }

  // 1つめのトークンをポイントする
  this->m_begin = str.find_first_not_of(this->m_delimiter);
  this->m_end = str.find_first_of(this->m_delimiter, this->m_begin);
}

TlStringTokenizer::~TlStringTokenizer() {}

std::size_t TlStringTokenizer::countTokens() {
  if (this->m_count >= 0) {
    // 既に計算されているので制御を戻す
    return m_count;
  }

  std::string::size_type n = 0;
  std::string::size_type i = 0;

  for (;;) {
    // 1つめのトークンに進む
    if ((i = this->m_str.find_first_not_of(this->m_delimiter, i)) ==
        std::string::npos) {
      break;
    }

    // 次の区切り文字に進む
    i = this->m_str.find_first_of(this->m_delimiter, i + 1);
    n++;

    if (i == std::string::npos) {
      break;
    }
  }

  return (this->m_count = n);
}

bool TlStringTokenizer::hasMoreTokens() {
  return (this->m_begin != this->m_end);
}

std::string TlStringTokenizer::nextToken() {
  std::string s = "";

  if (this->m_begin != std::string::npos) {
    if (this->m_end != std::string::npos) {
      s = this->m_str.substr(this->m_begin, (this->m_end - this->m_begin));
      this->m_begin =
          this->m_str.find_first_not_of(this->m_delimiter, this->m_end);
      this->m_end = this->m_str.find_first_of(this->m_delimiter, this->m_begin);
    } else {
      s = this->m_str.substr(this->m_begin,
                             (this->m_str.length() - this->m_begin));
      this->m_begin =
          this->m_str.find_first_not_of(this->m_delimiter, this->m_end);
    }
  }

  return s;
}
