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

#ifndef TLGETOPT_H
#define TLGETOPT_H

#include <map>
#include <string>
#include <vector>

// コマンドラインオプションを処理するクラス
//
// コンストラクタにargc, argvを設定し、listに所定の書式で
// 読み込むオプションを指定すると、このオブジェクトに対し、
// []メソッドでオプションの値にアクセスできる。
// オプションの書式は、いわゆるgetoptの書式に準ずる。
// long optionは指定できない。
//
// 引数の値には、0から始まるインデックスで、
// 同じく[]メソッドでアクセスできる。
// argv[0]の値はコマンド名になる前提で、
// オプション解析をargv[1]から処理しているが、
// argv[0]がコマンド名になるかは、実装定義であるので注意が必要である。
//
// オプションフラグが指定された場合は"defined"が格納される。
// 引数を受け取るオプションフラグ(文字の後に":"をつける)が指定された場合は、
// 文字列が格納される。
// それ以外は空文字が格納されている。

class TlGetopt {
public:
    TlGetopt(int argc, char* argv[], const char* list);
    ~TlGetopt();

    // accession
public:
    int getCount() const {
        return this->args_.size();
    };
    const std::string operator[](const std::string& sKey) const;
    const std::string operator[](unsigned int n) const;

protected:
    void initialize();
    void parseArgv(int argc, char* argv[], const char* list);

protected:
    /**
     *  引数の数
     */
    int m_nCount;

    /**
     *  データ保持用
     */
    std::map<std::string, std::string> m_Data;
    std::vector<std::string> args_;

    /**
     *  エラー文字列
     */
    std::string m_sError;
};

#endif  // TLGETOPT_H
