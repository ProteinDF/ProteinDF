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

#include <iostream>
#include <sstream>

#include "TlGetopt.h"
#include "TlUtils.h"

/**
 *  @param argc [in] main()のargc
 *  @param argv [in] main()のargv
 *  @param list [in] 有効にするオプション(引数を取るものはコロンをつける)
 */
TlGetopt::TlGetopt(int argc, char* argv[], const char* list) {
    this->initialize();

    // argv[0]は通常コマンド名が入る
    this->args_.push_back(argv[0]);
    // this->m_Data["0"] = std::string(argv[0]);
    argc--;
    argv++;
    // this->m_nCount++;

    this->parseArgv(argc, argv, list);
}

TlGetopt::~TlGetopt() {}

/**
 *  メンバ変数の初期化
 *
 */
void TlGetopt::initialize() {
    // this->m_nCount = 0;
    this->m_Data.clear();
    this->args_.clear();
    this->m_sError = "";
}

/**
 *  []演算子(参照用)
 *
 */
const std::string TlGetopt::operator[](const std::string& sKey) const {
    std::string answer = "";

    std::map<std::string, std::string>::const_iterator p =
        this->m_Data.find(sKey);
    if (p != this->m_Data.end()) {
        answer = p->second;
    }

    return answer;
}

/**
 *  []演算子(参照用)
 *
 */
const std::string TlGetopt::operator[](unsigned int n) const {
    std::string answer = "";

    if (n < this->args_.size()) {
        answer = this->args_[n];
    }

    return answer;
}

/**
 *  オプション解析
 *
 *  オプションを解析し、結果をメンバ変数m_Dataに代入する。
 *  @param argc main()のargc-1
 *  @param argv[] main()のargv++
 *  @param list 有効にするオプション(引数を取るものはコロンをつける)
 */
//   aOptListのkeyにはオプション文字列を入れる
//   そのvalueには""(null; オプションのみ), ":"(引数が必須)
void TlGetopt::parseArgv(int argc, char* argv[], const char* list) {
    // parse list
    //
    std::string sList(list);
    std::map<std::string, std::string> aOptList;
    aOptList.clear();
    std::string prev_s = "";
    while (!sList.empty()) {
        if (sList.at(0) == ':') {
            aOptList[prev_s] = ":";
        } else {
            std::string tmp(sList, 0, 1);
            aOptList[tmp] = "";
            prev_s = tmp;
        }
        sList = sList.substr(1);
    }

    // optionのparseをしている時はtrue, 引数の時はfalseを表すフラグ
    bool bOptionSection = true;
    // optionの引数用フラグ兼ハッシュのkey。""の場合flagなし
    std::string sOptArg = "";
    for (int i = 0; i < argc; i++) {
        std::string str = std::string(argv[i]);
        if (str.size() == 0) {
            continue;
        }

        // strの先頭が'-'でないものはオプションではない
        if (str.at(0) != '-') {
            bOptionSection = false;
        }

        // option, 引数ごとの処理
        if (bOptionSection) {
            // optionの処理
            str = str.substr(1);  // str[0]は'-'のハズのはずなので、取り除く

            if (str.at(0) == '-') {  // 先頭から２つ目が'-'の場合も取り除く
                str = str.substr(1);
            }

            if (str.empty()) {  // strが'-'または'--'はoptionの終わり。次は引数
                bOptionSection = false;
            }

            int nLen = str.length();  // length()は'\0'を含める
            while (!str.empty()) {
                std::string check_char(str, 0, 1);
                for (std::map<std::string, std::string>::const_iterator p =
                         aOptList.begin();
                     p != aOptList.end(); p++) {
                    if (check_char == p->first) {
                        if (nLen == 1 &&
                            p->second == ":") {  // optionは引数を取る
                            if (i + 1 < argc) {
                                this->m_Data[check_char] =
                                    std::string(argv[i + 1]);
                                i++;
                                break;
                            } else {  // error! optionが引数を取らなかった
                                std::stringstream s;
                                s << "illegal option: " << check_char;
                                m_sError = s.str();
                            }
                        } else if (p->second == "") {  // optionは引数を取らない
                            // ''の場合
                            this->m_Data[check_char] = "defined";
                            break;
                        } else {  // 引数リストに無いoption
                            std::stringstream s;
                            s << "illegal option: " << check_char;
                            m_sError = s.str();
                        }
                    }
                }
                str = str.substr(1);
            }

        } else {
            // 引数の処理
            this->args_.push_back(str);
            // std::cout << "m_nCount=" << this->m_nCount << std::endl;
            // this->m_Data[TlUtils::xtos<int>(this->m_nCount)] = str;
            // this->m_nCount++;
        }
    }
}

// EOF
