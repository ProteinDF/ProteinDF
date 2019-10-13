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

#include "CnError.h"
#include "DfObject.h"
#include "common.h"
#include "tl_dense_vector_lapack.h"

/// KS行列の対角化を行う
class DfDiagonal : public DfObject {
   public:
    DfDiagonal(TlSerializeData* pPdfData);
    virtual ~DfDiagonal();

   public:
    virtual void run();
    virtual void runQclo(DfObject::RUN_TYPE runType,
                         const std::string& fragname, int norbcut);

   protected:
    template <typename GeneralMatrix, typename SymmetricMatrix, typename Vector>
    inline void run_impl();

    template <typename MatrixType, typename SymmetricMatrixType,
              typename Vector>
    void main(DfObject::RUN_TYPE runType, const std::string& fragname = "",
              bool bPdfQcloMode = false);
};

// template
template <typename GeneralMatrix, typename SymmetricMatrix, typename Vector>
inline void DfDiagonal::run_impl() {
    // output informations
    switch (this->m_nMethodType) {
        case METHOD_RKS:
            this->main<GeneralMatrix, SymmetricMatrix, Vector>(RUN_RKS);
            break;

        case METHOD_UKS:
            this->main<GeneralMatrix, SymmetricMatrix, Vector>(RUN_UKS_ALPHA);
            this->main<GeneralMatrix, SymmetricMatrix, Vector>(RUN_UKS_BETA);
            break;

        case METHOD_ROKS:
            this->main<GeneralMatrix, SymmetricMatrix, Vector>(RUN_ROKS);
            break;

        default:
            CnErr.abort("DfDiagonal", "", "DfDiagMain",
                        "the value of scftype is illegal");
            break;
    }
}

template <typename GeneralMatrix, typename SymmetricMatrix, typename Vector>
inline void DfDiagonal::main(const DfObject::RUN_TYPE runType,
                             const std::string& fragname, bool bPdfQcloMode) {
    // input file
    SymmetricMatrix Fprime = this->getFprimeMatrix<SymmetricMatrix>(
        runType, this->m_nIteration, fragname);

    // diagonal
    Vector eigVal;
    GeneralMatrix eigVec;
    Fprime.eig(&eigVal, &eigVec);

    // save eigvec
    DfObject::saveCprimeMatrix(runType, this->m_nIteration, fragname, eigVec);
    // save eigval
    eigVal.save(this->getEigenvaluesPath(runType, this->m_nIteration));
}

#endif  // DFDIAGONAL_H
