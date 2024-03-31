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

#ifndef DFKMATRIX_H
#define DFKMATRIX_H

#include "DfObject.h"
#include "tl_dense_symmetric_matrix_lapack.h"

class DfKMatrix : public DfObject {
public:
    DfKMatrix(TlSerializeData* pPdfParam);
    virtual ~DfKMatrix();

public:
    void buildK();

protected:
    // virtual void getK_RI();
    virtual void getK_CD();
    virtual void getK_conventional();

    // void getK_RI_local(const RUN_TYPE runType, TlDenseSymmetricMatrix_Lapack
    // *pK);
    virtual void getK_CD_replica(const RUN_TYPE runType);

    void getK_conventional_local(const RUN_TYPE runType,
                                 TlDenseSymmetricMatrix_Lapack* pK);

protected:
    virtual TlDenseSymmetricMatrix_Lapack getKMatrix(const RUN_TYPE runType,
                                                     const int iteration);
    virtual void saveKMatrix(const RUN_TYPE runType,
                             const TlDenseSymmetricMatrix_Lapack& K);

    template <class SymmetricMatrixType>
    SymmetricMatrixType getDiffDensityMatrix(RUN_TYPE runType);

    template <class SymmetricMatrixType>
    SymmetricMatrixType getDensityMatrix(RUN_TYPE runType);
};

template <class SymmetricMatrixType>
SymmetricMatrixType DfKMatrix::getDiffDensityMatrix(const RUN_TYPE runType) {
    SymmetricMatrixType diffP;
    switch (runType) {
        case RUN_RKS:
            diffP = DfObject::getDiffDensityMatrix<SymmetricMatrixType>(
                RUN_RKS, this->m_nIteration);
            diffP *= 0.5;
            break;

        default:
            // RUN_UKS_ALPHA, RUN_UKS_BETA, RUN_ROKS_CLOSE, RUN_ROKS_OPEN
            diffP = DfObject::getDiffDensityMatrix<SymmetricMatrixType>(
                runType, this->m_nIteration);
            break;
    }

    return diffP;
}

template <class SymmetricMatrixType>
SymmetricMatrixType DfKMatrix::getDensityMatrix(const RUN_TYPE runType) {
    SymmetricMatrixType P;
    switch (runType) {
        case RUN_RKS:
            P = DfObject::getPInMatrix<SymmetricMatrixType>(RUN_RKS, this->m_nIteration);
            P *= 0.5;
            break;

        default:
            // RUN_UKS_ALPHA, RUN_UKS_BETA, RUN_ROKS_CLOSE, RUN_ROKS_OPEN
            P = DfObject::getPInMatrix<SymmetricMatrixType>(runType, this->m_nIteration);
            break;
    }

    return P;
}

#endif  // DFKMATRIX_H
