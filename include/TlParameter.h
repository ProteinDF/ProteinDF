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

#ifndef TLPARAMETER_H
#define TLPARAMETER_H

#include <string>
#include <map>
#include <fstream>
#include "TlUtils.h"

/// パラメータコンテナクラス
///
/// パラメータのグループ化と値の格納を目的としたクラス
/// param[group][keyword] = value のような使い方を意図している。
class TlParameter {
public:
    /// パラメータクラス
    ///
    /// パラメータ部分の実装
    class Param {
    public:
        /// デフォルトコンストラクタ
        Param();

        /// コピーコンストラクタ
        Param(const Param& rhs);

        /// デストラクタ
        ~Param();

    public:
        /// 代入演算子
        Param& operator=(const Param& rhs);

        /// [] 演算子(参照用)
        ///
        /// @param[in] sKey キーワード
        /// @return キーワードの対応する文字列。該当する値が無い場合は空文字列を返す。
        std::string  operator[](const std::string& sKey) const;

        /// [] 演算子(代入用)
        ///
        /// @param[in] sKey キーワード
        std::string& operator[](const std::string& sKey);

        /// 格納されているキーワードの数を返す
        ///
        /// @return キーワード数
        int size() const;

        std::vector<std::string> getKeywords() const;

        /// TlParameter::Paramオブジュクトをマージする
        ///
        /// 新たなキーワードは追加される。
        /// 既に存在するキーワードは上書きされる。
        ///
        /// @param[in] rhs 上書きするTlParameter::Paramオブジュクト
        void merge(const Param& rhs);

        /// TlParameter::Paramオブジュクトをマージする
        ///
        /// 新たなキーワードは追加される。
        /// 既に存在するキーワードは上書きされる。
        ///
        /// @param[in] rhs 上書きするTlParameter::Paramオブジュクト
        Param& operator+=(const Param& rhs);

        /// コンテナの内容を出力する
        ///
        /// 出力先のオブジェクトにはstd::coutなどを想定している。
        /// 出力先のオブジェクトは << 演算子を定義する必要がある。
        ///
        /// @param[in,out] out 出力先オブジェクト
        template<typename T>
        void print(T& out) const;

    private:
        typedef std::map<std::string, std::string> ParamData;
        ParamData m_data; /// パラメータ保持用
    };

    //////////////////////////////////////////////////////////////////////
public:
    /// デフォルトコンストラクタ
    TlParameter();

    /// コピーコンストラクタ
    TlParameter(const TlParameter& rhs);

    /// デストラクタ
    virtual ~TlParameter();

public:
    typedef std::map<std::string, Param> Group;
    typedef Group::const_iterator const_iterator;
    typedef Group::iterator iterator;

public:
    /// 要素(キーワードグループ)数を返す
    int size() const;

//   TlParameter::const_iterator begin() const;
//   TlParameter::const_iterator end() const;
//   TlParameter::iterator begin();
//   TlParameter::iterator end();

    /// キーワードグループのリストを返す
    std::vector<std::string> getGroups() const;


public:
    /// 入出力に用いるファイルパスを設定する
    ///
    /// @param[in] sFilePath ファイルパス
    void setFilePath(const std::string& sFilePath);

    /// 入出力に用いるファイルパスを返す
    ///
    /// @return ファイルパス
    std::string getFilePath() const;

public:
    /// 代入演算子
    TlParameter& operator=(const TlParameter& rhs);

    /// [] 演算子(参照用)
    ///
    /// @param[in] sGroupKey グループ名
    /// @return Paramオブジェクト
    Param  operator[](const std::string& sGroupKey) const;

    /// [] 演算子(代入用)
    ///
    /// @param[in] sGroupKey グループ名
    /// @return Paramオブジェクト
    Param& operator[](const std::string& sGroupKey);

    /// TlParameterオブジェクトをマージする
    ///
    /// 既にあるグループ・キーワードは上書きされる
    /// @param[in] rhs 上書きするTlParameterオブジェクト
    void merge(const TlParameter& rhs);

    /// TlParameterオブジェクトをマージする
    ///
    /// 既にあるグループ・キーワードは上書きされる
    /// @param[in] rhs 上書きするTlParameterオブジェクト
    TlParameter& operator+=(const TlParameter& rhs);

public:
    /// 指定されたファイルからデータを読み込む
    ///
    /// @param[in] sFilePath ファイルパス名
    /// @retval true 読み込み成功
    /// @retval false 読み込み失敗
    bool load(const std::string& sFilePath);

    /// 指定されたstd::ifstreamオブジェクトからデータを読み込む
    ///
    /// @param[in,out] ifs std::ifstreamオブジェクト
    /// @retval true 読み込み成功
    /// @retval false 読み込み失敗
    bool load(std::ifstream& ifs);

    /// ファイルにデータを出力する
    ///
    /// @retval true 読み込み成功
    /// @retval false 読み込み失敗
    virtual bool save() const;

    /// 指定されたファイルパスにデータを出力する
    ///
    /// @param[in] sFilePath ファイルパス
    /// @retval true 読み込み成功
    /// @retval false 読み込み失敗
    bool save(const std::string& sFilePath) const;

    /// std::ofstreamオブジェクトにデータを出力する
    ///
    /// @param[out] ofs 出力用std::ofstreamオブジェクト
    /// @retval true 読み込み成功
    /// @retval false 読み込み失敗
    bool save(std::ofstream& ofs) const;

public:
    /// 全内容をテキスト出力する
    ///
    /// outオブジェクトは <<演算子を定義してある必要がある.
    template<typename T>
    void print(T& out) const;

    /// 指定されたブロックの内容をテキスト出力する
    ///
    /// outオブジェクトは <<演算子を定義してある必要がある.
    template<typename T>
    void print(const std::string& sBlock, T& out) const;

protected:
    std::string m_sFilePath; /// ファイルパスを保持する
    Group m_data;            /// データ保持
};

////////////////////////////////////////////////////////////////////////
// template
template<typename T>
void TlParameter::Param::print(T& out) const
{
    for (ParamData::const_iterator q = this->m_data.begin(); q != this->m_data.end(); q++) {
        if (q->first == "") {
            continue;
        }

        const std::string key   = "[" + q->first  + "]";
        const std::string value = "[" + q->second + "]";

        std::string line = key;
        if ((line.length() + value.length()) < 72) {
            TlUtils::pad(line, (72 - value.length()), ' ');
        }
        line += value;

        out << line << "\n";
    }
}

template<typename T>
void TlParameter::print(T& out) const
{
//   out << "//-----------------------------------------------------------------------\n";
//   out << "// keywords and values for the calculation                             //\n";
//   out << "//-----------------------------------------------------------------------\n";

    for (Group::const_iterator p = this->m_data.begin(); p != this->m_data.end(); p++) {
        out << TlUtils::format(">>>>%s // Number of stored keywords is %2d\n",
                               p->first.c_str(), p->second.size());

        p->second.print(out);
        out << "\n";
    }
    out << "//-----------------------------------------------------------------------\n";

    out.flush();
}

template<typename T>
void TlParameter::print(const std::string& sBlock, T& out) const
{
    Group::const_iterator p = this->m_data.find(sBlock);
    if (p != this->m_data.end()) {
        out << TlUtils::format(">>>>%s // Number of stored keywords is %2d\n",
                               sBlock.c_str(), p->second.size());
        p->second.print(out);
        out << "\n";
    }

    out.flush();
}


#endif // TLPARAMETER_H
