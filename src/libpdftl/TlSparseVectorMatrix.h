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

#ifndef TLSPARSEVECTORMATRIX_H
#define TLSPARSEVECTORMATRIX_H

#include <algorithm>
#include <list>
#include <map>
#include "TlMatrixObject.h"

#define SVM_DEFAULT_CACHE_SIZE (100UL * 1024UL * 1024UL) // 100 MB

/// 疎行列クラス
/// 行列要素は行をインデックスとする列ベクトルで保持する
class TlSparseVectorMatrix : public TlMatrixObject {
public:
    /// 行列オブジェクトを作成する
    ///
    /// @param[in] size 作成する行列の次元数
    explicit TlSparseVectorMatrix(index_type row =1, index_type col =1,
                                  std::size_t cacheMemSize = SVM_DEFAULT_CACHE_SIZE);

    /// コピーコンストラクタ
    TlSparseVectorMatrix(const TlSparseVectorMatrix& rhs);
    
    /// デストラクタ
    virtual ~TlSparseVectorMatrix();

public:
    virtual int getNumOfRows() const;
    virtual int getNumOfCols() const;

    virtual std::size_t getMemSize() const;

    virtual double get(index_type row, index_type col) const;
    virtual void set(index_type row, index_type col, double value);
    virtual void add(index_type row, index_type col, double value);

    virtual TlVector getRowVector(index_type row) const;
    virtual TlVector getColVector(index_type col) const;

    /// 行ベクトルを設定する
    ///
    /// isKeeped指定は上書きされます。すなわち、ある行インデックスxに対しisKeeped=true指定した後、
    /// 再度行インデックスxに対しisKeeped=falseを指定した場合は、行xはキャッシュ削除対象になります。
    /// @param [in] isKeeped 規定のメモリを越えた場合でも保持すべき(true)かどうかを指定します。
    /// 
    void setRowVector(index_type row, const TlVector& cols, bool isKeeped = false);
    bool findRow(index_type row) const;
    
public:
    virtual bool load(const std::string& path) {
        return false;
    }

    virtual bool save(const std::string& path) const {
        return false;
    }

protected:
    void updateCache(const std::size_t row);
    
protected:
    index_type numOfRows_;
    index_type numOfCols_;
    typedef std::map<index_type, TlVector> DataType;
    DataType data_;

    std::size_t cacheMemSize_;
    std::list<int> cacheIndexList_;
};

#endif // TLSPARSEVECTORSYMMETRICMATRIX_H

