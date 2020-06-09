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

#ifndef TLROWVECTORMATRIX_H
#define TLROWVECTORMATRIX_H

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif  // HAVE_CONFIG_H

#include "tl_dense_general_matrix_lapack.h"
#include "tl_dense_matrix_arrays_object.h"
#include "tl_matrix_object.h"

/// 行ベクトルで表される行列
/// 行個数の、列個サイズのベクトルが集合したコンテナクラス
class TlDenseGeneralMatrix_arrays_RowOriented
    : public TlDenseMatrix_arrays_Object {
   public:
    explicit TlDenseGeneralMatrix_arrays_RowOriented(
        index_type row = 1, index_type col = 1, int numOfSubunits = 1,
        int subunitID = 0, bool isUsingMemManager = false);
    explicit TlDenseGeneralMatrix_arrays_RowOriented(
        const TlDenseGeneralMatrix_Lapack& rhs, int numOfSubunits = 1,
        int subunitID = 0, bool isUsingMemManager = false);
    TlDenseGeneralMatrix_arrays_RowOriented(
        const TlDenseGeneralMatrix_arrays_RowOriented& rhs);

    virtual ~TlDenseGeneralMatrix_arrays_RowOriented();

   public:
    virtual void resize(index_type row, index_type col);
    void reserveColSize(const index_type reserveColSize);

    index_type getNumOfRows() const { return this->getNumOfVectors(); };

    index_type getNumOfCols() const { return this->getSizeOfVector(); };

    virtual void set(index_type row, index_type col, double value);
    virtual void add(index_type row, index_type col, double value);
    virtual double get(index_type row, index_type col) const;

    virtual std::vector<double> getRowVector(const index_type row) const;
    virtual void getRowVector(const index_type row, double* pBuf,
                              const index_type length) const;

   public:
    /// TlDenseGeneralMatrix_Lapackオブジェクトを返す(for debug)
    TlDenseGeneralMatrix_Lapack getTlMatrixObject() const;

    /// TlDenseGeneralMatrix_arrays_ColOriented形式で保存
    void saveByTlDenseGeneralMatrix_arrays_ColOriented(
        const std::string& basename) const;
};

std::ostream& operator<<(std::ostream& stream,
                         const TlDenseGeneralMatrix_arrays_RowOriented& mat);

bool RowVectorMatrix2CSFD(const std::string& rvmBasePath,
                          const std::string& csfdPath, bool verbose = false,
                          bool showProgress = false);

#endif  // TLROWVECTORMATRIX_H
