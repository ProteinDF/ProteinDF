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

#ifndef TLVECTORMATRIXOBJECT_H
#define TLVECTORMATRIXOBJECT_H

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif // HAVE_CONFIG_H

#include "TlMatrixObject.h"

/// 配列の配列として行列を扱うためのコンテナ
class TlVectorMatrixObject {
public:
    typedef TlMatrixObject::index_type index_type;

public:
    explicit TlVectorMatrixObject(index_type numOfVectors = 1,
                                  index_type sizeOfVector = 1, 
                                  int numOfMembers = 1, int ID = 0,
                                  bool isUsingMemManager = false);
    virtual ~TlVectorMatrixObject();
        
public:
    void resize(index_type newNumOfVectors,
                index_type newVectorSize);

    /// 値を代入する
    void set(index_type vectorIndex,
             index_type index,
             double value);

    double get(index_type vectorIndex,
               index_type index) const;

    std::vector<double> getVector(index_type vectorIndex) const;

    void save(const std::string& basename) const;
    void load(const std::string& basename);

protected:
    /// 前もってvector sizeを設定する
    void reserveVectorSize(index_type vectorSize);

public:
    index_type getSizeOfVector() const {
        return this->sizeOfVector_;
    };

    index_type getNumOfVectors() const {
        return this->numOfVectors_;
    };

    int getNumOfSubunits() const {
        return this->numOfSubunits_;
    };

    /// 自分のIDを返す
    int getSubunitID() const {
        return this->subunitID_;
    };

    /// 該当するベクトルを持っているPE番号を返す
    int getSubunitID(const index_type vectorIndex) const;


private:
    /// 暗黙なコピーコンストラクタの作成を防止する
    TlVectorMatrixObject(const TlVectorMatrixObject& rhs);

    /// 暗黙な代入演算子の作成を防止する
    TlVectorMatrixObject& operator=(const TlVectorMatrixObject& rhs);

private:
    index_type numOfVectors_;      /// ベクトルの総数(global)
    index_type sizeOfVector_;        /// 1ベクトルの大きさ

    /// 行列全体を構成するオブジェクトの総数
    /// 通常はプロセスの総数
    int numOfSubunits_; 

    /// このオブジェクトのID
    /// 通常はプロセスID
    int subunitID_;

    index_type numOfLocalVectors_; // ベクトルの総数(local)

    index_type reservedVectorSize_; /// あらかじめ保持している1ベクトルの大きさ

    bool isUsingMemManager_;     /// 独自のメモリマネージャを使う(true)
    std::vector<double* > data_; /// データ

};

#endif // TLCOLVECTORMATRIX2_H
