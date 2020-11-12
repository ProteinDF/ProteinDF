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

#include "DfConverge_Damping.h"

#include <cassert>

#include "TlTime.h"
#include "tl_dense_symmetric_matrix_lapack.h"
#include "tl_dense_vector_lapack.h"

DfConverge_Damping::DfConverge_Damping(TlSerializeData* pPdfParam)
    : DfConverge(pPdfParam) {
    const TlSerializeData& pdfParam = *pPdfParam;
    this->m_nStartIteration = std::max(
        pdfParam["scf_acceleration/damping/number_of_damping"].getInt(), 2);
    this->m_dDampingFactor =
        pdfParam["scf_acceleration/damping/damping_factor"].getDouble();
}

DfConverge_Damping::~DfConverge_Damping() {}

void DfConverge_Damping::convergeRhoTilde() {
    switch (this->m_nMethodType) {
        case METHOD_RKS:
            this->convergeRhoTilde<TlDenseVector_Lapack>(DfObject::RUN_RKS);
            break;
        case METHOD_UKS:
            this->convergeRhoTilde<TlDenseVector_Lapack>(
                DfObject::RUN_UKS_ALPHA);
            this->convergeRhoTilde<TlDenseVector_Lapack>(
                DfObject::RUN_UKS_BETA);
            break;
        case METHOD_ROKS:
            this->convergeRhoTilde<TlDenseVector_Lapack>(
                DfObject::RUN_ROKS_CLOSED);
            this->convergeRhoTilde<TlDenseVector_Lapack>(
                DfObject::RUN_ROKS_OPEN);
            break;
        default:
            std::cerr
                << "program error. @DfConverge_Damping::convergeRhoTilde()"
                << std::endl;
            break;
    }
}

void DfConverge_Damping::convergeKSMatrix() {
    switch (this->m_nMethodType) {
        case METHOD_RKS:
            this->convergeKSMatrix<TlDenseSymmetricMatrix_Lapack>(
                DfObject::RUN_RKS);
            break;
        case METHOD_UKS:
            this->convergeKSMatrix<TlDenseSymmetricMatrix_Lapack>(
                DfObject::RUN_UKS_ALPHA);
            this->convergeKSMatrix<TlDenseSymmetricMatrix_Lapack>(
                DfObject::RUN_UKS_BETA);
            break;
        case METHOD_ROKS:
            this->convergeKSMatrix<TlDenseSymmetricMatrix_Lapack>(
                DfObject::RUN_ROKS_CLOSED);
            this->convergeKSMatrix<TlDenseSymmetricMatrix_Lapack>(
                DfObject::RUN_ROKS_OPEN);
            break;
        default:
            std::cerr
                << "program error. @DfConverge_Damping::convergeKSMatrix()"
                << std::endl;
            break;
    }
}

void DfConverge_Damping::convergePMatrix() {
    switch (this->m_nMethodType) {
        case METHOD_RKS:
            this->convergePMatrix<TlDenseSymmetricMatrix_Lapack>(
                DfObject::RUN_RKS);
            break;
        case METHOD_UKS:
            this->convergePMatrix<TlDenseSymmetricMatrix_Lapack>(
                DfObject::RUN_UKS_ALPHA);
            this->convergePMatrix<TlDenseSymmetricMatrix_Lapack>(
                DfObject::RUN_UKS_BETA);
            break;
        case METHOD_ROKS:
            this->convergePMatrix<TlDenseSymmetricMatrix_Lapack>(
                DfObject::RUN_ROKS_CLOSED);
            this->convergePMatrix<TlDenseSymmetricMatrix_Lapack>(
                DfObject::RUN_ROKS_OPEN);
            break;
        default:
            std::cerr << "program error. @DfConverge_Damping::convergePMatrix()"
                      << std::endl;
            break;
    }
}

double DfConverge_Damping::getDampingFactor() const {
    const double dampingFactor = (*this->pPdfParam_)["scf_acceleration/damping/damping_factor"].getDouble();

    return dampingFactor;
}
