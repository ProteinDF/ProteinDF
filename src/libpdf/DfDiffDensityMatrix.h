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

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif  // HAVE_CONFIG_H

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
#ifdef HAVE_LAPACK
    void exec_lapack();
#endif  // HAVE_LAPACK

#ifdef HAVE_EIGEN
    void exec_eigen();
#endif  // HAVE_EIGEN

protected:
    template <class SymmetricMatrixType>
    void calc(DfObject::RUN_TYPE runType, int iteration);
};

template <class SymmetricMatrixType>
void DfDiffDensityMatrix::calc(const DfObject::RUN_TYPE runType, const int iteration) {
    SymmetricMatrixType P = DfObject::getPInMatrix<SymmetricMatrixType>(runType, iteration);
    if (iteration > 1) {
        P -= (DfObject::getPInMatrix<SymmetricMatrixType>(runType, iteration - 1));
    }

    // Maximum value of the difference electron density matrix
    const double maxDiff = P.getMaxAbsoluteElement();
    this->log_.info(TlUtils::format("max abs. of diff-density mat.: %e", maxDiff));
    (*this->pPdfParam_)["DfDiffDensityMatrix"]["max_abs"][iteration] = maxDiff;

    this->saveDiffDensityMatrix(runType, iteration, P);
}

#endif  // DFDIFFDENSITYMATRIX_H
