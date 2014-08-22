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

#include "DfConverge_Anderson.h"
#include "TlSymmetricMatrix.h"

DfConverge_Anderson::DfConverge_Anderson(TlSerializeData* pPdfParam)
    : DfConverge_Damping(pPdfParam)
{
    const TlSerializeData& pdfParam = *pPdfParam;
    this->m_nStartIterationOfAnderson =
        std::max(pdfParam["scf_acceleration/anderson/start_number"].getInt(), 3);

    this->m_dDampingFactorOfAnderson =
        pdfParam["scf_acceleration/anderson/damping_factor"].getDouble();

}

DfConverge_Anderson::~DfConverge_Anderson()
{
}


void DfConverge_Anderson::convergeRhoTilde()
{
    switch (this->m_nMethodType) {
    case METHOD_RKS:
        this->convergeRhoTilde<TlVector>(DfObject::RUN_RKS);
        break;
    case METHOD_UKS:
        this->convergeRhoTilde<TlVector>(DfObject::RUN_UKS_ALPHA);
        this->convergeRhoTilde<TlVector>(DfObject::RUN_UKS_BETA);
        break;
    case METHOD_ROKS:
        this->convergeRhoTilde<TlVector>(DfObject::RUN_ROKS_CLOSED);
        this->convergeRhoTilde<TlVector>(DfObject::RUN_ROKS_OPEN);
        break;
    default:
        std::cerr << "program error. @DfConverge_Anderson::convergeRhoTilde()" << std::endl;
        break;
    }
}

void DfConverge_Anderson::convergeKSMatrix()
{
    switch (this->m_nMethodType) {
    case METHOD_RKS:
        this->convergeKSMatrix<TlSymmetricMatrix, TlVector>(DfObject::RUN_RKS);
        break;
    case METHOD_UKS:
        this->convergeKSMatrix<TlSymmetricMatrix, TlVector>(DfObject::RUN_UKS_ALPHA);
        this->convergeKSMatrix<TlSymmetricMatrix, TlVector>(DfObject::RUN_UKS_BETA);
        break;
    case METHOD_ROKS:
        this->convergeKSMatrix<TlSymmetricMatrix, TlVector>(DfObject::RUN_ROKS_CLOSED);
        this->convergeKSMatrix<TlSymmetricMatrix, TlVector>(DfObject::RUN_ROKS_OPEN);
        break;
    default:
        std::cerr << "program error. @DfConverge_Anderson::convergeKSMatrix()" << std::endl;
        break;
    }
}

void DfConverge_Anderson::convergePMatrix()
{
    switch (this->m_nMethodType) {
    case METHOD_RKS:
        this->convergePMatrix<TlSymmetricMatrix, TlVector>(DfObject::RUN_RKS);
        break;
    case METHOD_UKS:
        this->convergePMatrix<TlSymmetricMatrix, TlVector>(DfObject::RUN_UKS_ALPHA);
        this->convergePMatrix<TlSymmetricMatrix, TlVector>(DfObject::RUN_UKS_BETA);
        break;
    case METHOD_ROKS:
        this->convergePMatrix<TlSymmetricMatrix, TlVector>(DfObject::RUN_ROKS_CLOSED);
        this->convergePMatrix<TlSymmetricMatrix, TlVector>(DfObject::RUN_ROKS_OPEN);
        break;
    default:
        std::cerr << "program error. @DfConverge_Anderson::convergePMatrix()" << std::endl;
        break;
    }
}


