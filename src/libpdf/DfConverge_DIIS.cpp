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

#include "DfConverge_DIIS.h"
#include <cassert>
#include "TlTime.h"

DfConverge_DIIS::DfConverge_DIIS(TlSerializeData* pPdfParam)
    : DfConverge_Damping(pPdfParam) {
    const TlSerializeData& pdfParam = *pPdfParam;
    this->startIterationOfDIIS_ =
        std::max(pdfParam["scf_acceleration/diis/start"].getInt(), 1);
    this->numOfLastItems_ =
        std::max(pdfParam["scf_acceleration/diis/last_items"].getInt(), 3);

    if (this->pPdfParam_->hasKey("DIIS_e_max_pass_itr") == false) {
        (*this->pPdfParam_)["DIIS_e_max_pass_itr"] = -1;
    }
}

DfConverge_DIIS::~DfConverge_DIIS() {}

void DfConverge_DIIS::convergeRhoTilde() {
    this->log_.critical("NOT support: DIIS for rho.");

    // switch (this->m_nMethodType) {
    // case METHOD_RKS:
    //     this->convergeRhoTilde<TlDenseGeneralMatrix_Lapack,
    //     TlDenseSymmetricMatrix_Lapack,
    //     TlDenseVector_Lapack>(DfObject::RUN_RKS); break;
    // case METHOD_UKS:
    //     this->convergeRhoTilde<TlDenseGeneralMatrix_Lapack,
    //     TlDenseSymmetricMatrix_Lapack,
    //     TlDenseVector_Lapack>(DfObject::RUN_UKS_ALPHA);
    //     this->convergeRhoTilde<TlDenseGeneralMatrix_Lapack,
    //     TlDenseSymmetricMatrix_Lapack,
    //     TlDenseVector_Lapack>(DfObject::RUN_UKS_BETA);
    //     break;
    // case METHOD_ROKS:
    //     this->convergeRhoTilde<TlDenseGeneralMatrix_Lapack,
    //     TlDenseSymmetricMatrix_Lapack,
    //     TlDenseVector_Lapack>(DfObject::RUN_UKS_ALPHA);
    //     this->convergeRhoTilde<TlDenseGeneralMatrix_Lapack,
    //     TlDenseSymmetricMatrix_Lapack,
    //     TlDenseVector_Lapack>(DfObject::RUN_UKS_BETA);
    //     break;
    // default:
    //     std::cerr << "program error. @DfConverge_Damping::convergeRhoTilde()"
    //     << std::endl; break;
    // }
}

void DfConverge_DIIS::convergeKSMatrix() {
    switch (this->m_nMethodType) {
        case METHOD_RKS:
            this->convergeKSMatrix<TlDenseGeneralMatrix_Lapack,
                                   TlDenseSymmetricMatrix_Lapack>(
                DfObject::RUN_RKS);
            break;
        case METHOD_UKS:
            this->convergeKSMatrix<TlDenseGeneralMatrix_Lapack,
                                   TlDenseSymmetricMatrix_Lapack>(
                DfObject::RUN_UKS_ALPHA);
            this->convergeKSMatrix<TlDenseGeneralMatrix_Lapack,
                                   TlDenseSymmetricMatrix_Lapack>(
                DfObject::RUN_UKS_BETA);
            break;
        case METHOD_ROKS:
            this->convergeKSMatrix<TlDenseGeneralMatrix_Lapack,
                                   TlDenseSymmetricMatrix_Lapack>(
                DfObject::RUN_ROKS_CLOSED);
            this->convergeKSMatrix<TlDenseGeneralMatrix_Lapack,
                                   TlDenseSymmetricMatrix_Lapack>(
                DfObject::RUN_ROKS_OPEN);
            break;
        default:
            std::cerr
                << "program error. @DfConverge_Damping::convergeKSMatrix()"
                << std::endl;
            break;
    }
}

void DfConverge_DIIS::convergePMatrix() {
    switch (this->m_nMethodType) {
        case METHOD_RKS:
            this->convergePMatrix<TlDenseGeneralMatrix_Lapack,
                                  TlDenseSymmetricMatrix_Lapack>(
                DfObject::RUN_RKS);
            break;
        case METHOD_UKS:
            this->convergePMatrix<TlDenseGeneralMatrix_Lapack,
                                  TlDenseSymmetricMatrix_Lapack>(
                DfObject::RUN_UKS_ALPHA);
            this->convergePMatrix<TlDenseGeneralMatrix_Lapack,
                                  TlDenseSymmetricMatrix_Lapack>(
                DfObject::RUN_UKS_BETA);
            break;
        case METHOD_ROKS:
            this->convergePMatrix<TlDenseGeneralMatrix_Lapack,
                                  TlDenseSymmetricMatrix_Lapack>(
                DfObject::RUN_ROKS_CLOSED);
            this->convergePMatrix<TlDenseGeneralMatrix_Lapack,
                                  TlDenseSymmetricMatrix_Lapack>(
                DfObject::RUN_ROKS_OPEN);
            break;
        default:
            std::cerr << "program error. @DfConverge_Damping::convergePMatrix()"
                      << std::endl;
            break;
    }
}
