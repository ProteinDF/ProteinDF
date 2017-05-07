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

#ifndef DFDIAGONAL_H
#define DFDIAGONAL_H

#include <string>
#include "DfObject.h"
#include "TlVector.h"

/// KS行列の対角化を行う
class DfDiagonal : public DfObject {
public:
    DfDiagonal(TlSerializeData* pPdfData);
    virtual ~DfDiagonal();

public:
    virtual void DfDiagMain();
    virtual void DfDiagQclo(DfObject::RUN_TYPE runType, const std::string& fragname, int norbcut);

protected:
    template<typename MatrixType, typename SymmetricMatrixType>
    void main(DfObject::RUN_TYPE runType, const std::string& fragname ="", bool bPdfQcloMode = false);

};

// template
template<typename MatrixType, typename SymmetricMatrixType>
void DfDiagonal::main(const DfObject::RUN_TYPE runType, const std::string& fragname, bool bPdfQcloMode)
{
    // input file
    SymmetricMatrixType Fprime = this->getFprimeMatrix<SymmetricMatrixType>(runType, this->m_nIteration, fragname);

    // diagonal
    TlVector* pEigVal = new TlVector();
    MatrixType* pEigVec = new MatrixType();
    Fprime.diagonal(pEigVal, pEigVec);

    // save eigval
    pEigVal->save(this->getEigenvaluesPath(runType, this->m_nIteration));

    // save eigvec
    this->saveCprimeMatrix(runType, this->m_nIteration, fragname, *pEigVec);

    // finalize
    delete pEigVal;
    pEigVal = NULL;

    delete pEigVec;
    pEigVec = NULL;
}

#endif // DFDIAGONAL_H

