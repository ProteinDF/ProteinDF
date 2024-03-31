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

#include "DfDiffDensityMatrix.h"

#include "TlFile.h"
#include "tl_dense_symmetric_matrix_lapack.h"

DfDiffDensityMatrix::DfDiffDensityMatrix(TlSerializeData* pPdfData)
    : DfObject(pPdfData) {
    this->isSaveDiffMatrix_ = (TlUtils::toUpper((*pPdfData)["save_diff_density_matrix"].getStr()) == "YES");
}

DfDiffDensityMatrix::~DfDiffDensityMatrix() {}

void DfDiffDensityMatrix::exec() {
    // check memory
    // const std::size_t needMem =
    //     this->m_nNumOfAOs * (this->m_nNumOfAOs + 1) * sizeof(double);
    // if ((this->isWorkOnDisk_ == true) || (this->procMaxMemSize_ < needMem)) {
    //   this->log_.info(" The differencial density matrix is build on disk.");
    //   TlMatrix::useMemManager(true);
    // } else {
    //   this->log_.info(" The differencial density matrix is build on
    //   memory."); TlMatrix::useMemManager(false);
    // }

    this->log_.info("make diff matrix by using LAPACK");
    switch (this->m_nMethodType) {
        case METHOD_RKS:
            this->calc<TlDenseSymmetricMatrix_Lapack>(RUN_RKS, this->m_nIteration);
            break;

        case METHOD_UKS:
            this->calc<TlDenseSymmetricMatrix_Lapack>(RUN_UKS_ALPHA, this->m_nIteration);
            this->calc<TlDenseSymmetricMatrix_Lapack>(RUN_UKS_BETA, this->m_nIteration);
            break;

        case METHOD_ROKS:
            this->calc<TlDenseSymmetricMatrix_Lapack>(RUN_ROKS_CLOSED, this->m_nIteration);
            this->calc<TlDenseSymmetricMatrix_Lapack>(RUN_ROKS_OPEN, this->m_nIteration);
            break;

        default:
            std::cerr << "program error. " << __FILE__ << ":" << __LINE__
                      << std::endl;
            abort();
            break;
    }
}
