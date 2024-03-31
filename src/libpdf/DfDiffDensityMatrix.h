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

#ifndef DFDIFFDENSITYMATRIX_H
#define DFDIFFDENSITYMATRIX_H

#include "DfObject.h"

class DfDiffDensityMatrix : public DfObject {
public:
    DfDiffDensityMatrix(TlSerializeData* pPdfParam);
    virtual ~DfDiffDensityMatrix();

public:
    /// 差電子密度行列を求める
    ///
    virtual void exec();

protected:
    template <class SymmetricMatrixType>
    void calc(DfObject::RUN_TYPE runType, int iteration);

    // template<class SymmetricMatrixType>
    // void calc_ROKS();

protected:
    /// 差電子密度行列をディスクに保存する(true)かどうか
    bool isSaveDiffMatrix_;
};

template <class SymmetricMatrixType>
void DfDiffDensityMatrix::calc(const DfObject::RUN_TYPE runType, const int iteration) {
    SymmetricMatrixType P = DfObject::getPInMatrix<SymmetricMatrixType>(runType, iteration);
    if (iteration > 1) {
        P -= (DfObject::getPInMatrix<SymmetricMatrixType>(runType, iteration - 1));
    }

    this->saveDiffDensityMatrix(runType, iteration, P);
}

// template<class SymmetricMatrixType>
// void DfDiffDensityMatrix::calc_ROKS()
// {
//     SymmetricMatrixType P_close;
//     P_close.load(DfObject::getP1pqMatrixPath(this->m_nIteration -1));
//     SymmetricMatrixType P_open;
//     P_open.load(DfObject::getP2pqMatrixPath(this->m_nIteration -1));
//     P_open += P_close;

//     {
//         SymmetricMatrixType prevP_open;
//         prevP_open.load(DfObject::getP2pqMatrixPath(this->m_nIteration
//         -2)); P_open -= prevP_open;
//     }

//     {
//         SymmetricMatrixType prevP_close;
//         prevP_close.load(DfObject::getP1pqMatrixPath(this->m_nIteration
//         -2)); P_open -= prevP_close; P_close -= prevP_close;
//     }

//     P_open.save(this->getDiffDensityMatrixPath(DfObject::RUN_ROKS_OPEN,
//     this->m_nIteration));
//     P_close.save(this->getDiffDensityMatrixPath(DfObject::RUN_ROKS_CLOSED,
//     this->m_nIteration));
// }

#endif  // DFDIFFDENSITYMATRIX_H
