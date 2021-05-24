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

#ifndef TL_DENSE_MATRIX_ARRAYS_OBJECT_H
#define TL_DENSE_MATRIX_ARRAYS_OBJECT_H

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif  // HAVE_CONFIG_H

#include <vector>

#include "tl_dense_vector_lapack.h"
#include "tl_matrix_object.h"

/// 配列の配列として行列を扱うためのコンテナ
class TlDenseMatrix_arrays_Object : public TlMatrixObject {
   public:
    explicit TlDenseMatrix_arrays_Object(index_type numOfVectors = 1, index_type sizeOfVector = 1, int numOfMembers = 1,
                                         int ID = 0, bool isUsingMemManager = false);
    TlDenseMatrix_arrays_Object(const TlDenseMatrix_arrays_Object& rhs);

    virtual ~TlDenseMatrix_arrays_Object();

    TlDenseMatrix_arrays_Object& operator=(const TlDenseMatrix_arrays_Object& rhs);

   public:
    virtual index_type getNumOfRows() const = 0;
    virtual index_type getNumOfCols() const = 0;
    int getSizeOfChunk() const;

    /// インスタンスのメモリサイズを返す
    virtual std::size_t getMemSize() const;

   private:
    virtual double get(index_type row, index_type col) const = 0;
    virtual void set(index_type row, index_type col, double value) = 0;
    virtual void add(index_type row, index_type col, double value) = 0;

   public:
    virtual bool load(const std::string& path);
    virtual bool save(const std::string& path) const;

   public:
    virtual bool loadSubunitFile(const std::string& path);
    virtual bool saveSubunitFileWithReservedSizeOfVector(const std::string& basename) const;

   public:
    /// @param (subunitID) >= 0 add suffix, < 0 open the basename as path
    virtual bool load(const std::string& basename, int subunitID);

   public:
    void resize(index_type newNumOfVectors, index_type newVectorSize);

    /// 値を代入する
    void set_to_vm(index_type vectorIndex, index_type index, double value);

    void add_to_vm(index_type vectorIndex, index_type index, double value);

    double get_from_vm(index_type vectorIndex, index_type index) const;

    std::vector<double> getVector(index_type vectorIndex) const;
    void getVector(const index_type vectorIndex, double* pBuf, const index_type length) const;
    void setVector(index_type vectorIndex, const std::vector<double>& v);

    /// copy chunk data
    ///
    //  retval: size of copied buffer
    std::size_t getChunk(const index_type vectorIndex, double* pBuf, const std::size_t length) const;

   public:
    static std::string getFileName(const std::string& basename, const int subunitID);

    /// 指定したファイルがこのクラスで読み込めるかどうかをチェックする。
    /// 読み込める場合はtrueを返す。
    /// このとき、pNumOfSubunitsが指定されていれば、全サブユニット数が代入される。
    /// また、pSubunitIDが指定されていれば、サブユニットIDが代入される。
    static bool isLoadable(const std::string& filepath, index_type* pNumOfVectors = NULL,
                           index_type* pSizeOfVector = NULL, int* pNumOfSubunits = NULL, int* pSubunitID = NULL,
                           int* pSizeOfChunk = NULL);

   public:
    index_type getSizeOfVector() const {
        return this->sizeOfVector_;
    };

    index_type getNumOfVectors() const {
        return this->numOfVectors_;
    };

    index_type getNumOfLocalVectors() const {
        // return this->numOfLocalVectors_;
        const index_type numOfLocalVectors =
            this->getNumOfLocalVectors(this->numOfVectors_, this->numOfSubunits_, this->sizeOfChunk_);
        return numOfLocalVectors;
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

   protected:
    /// 前もってvector sizeを設定する
    void reserveVectorSize(index_type vectorSize);

    /// data_ メンバ変数を破棄する
    void destroy();

   public:
    static index_type getNumOfLocalChunks(const index_type numOfVectors, const int numOfSubunits,
                                          const int sizeOfChunk);

   protected:
    // defunct
    static index_type getNumOfLocalVectors(const index_type numOfVectors, const int numOfSubunits,
                                           const int sizeOfChunk);

    index_type getLocalVectorIndex(const index_type vectorIndex, int* pSubunitId, int* pLocalChunkId = NULL,
                                   int* pLocalChunkVectorIndex = NULL) const;
    bool saveByTheOtherType(const std::string& basename) const;

   private:
    index_type numOfVectors_;  /// ベクトルの総数(global)
    index_type sizeOfVector_;  /// 1ベクトルの大きさ
    index_type sizeOfChunk_;   ///

    /// 行列全体を構成するオブジェクトの総数
    /// 通常はプロセスの総数
    int numOfSubunits_;

    /// このオブジェクトのID
    /// 通常はプロセスID
    int subunitID_;

    index_type reservedVectorSize_;  /// あらかじめ保持している1ベクトルの大きさ

    bool isUsingMemManager_;  /// 独自のメモリマネージャを使う(true)

    int numOfLocalChunks_;
    std::vector<double*> chunks_;
};

#endif  // TL_DENSE_MATRIX_ARRAYS_OBJECT_H
