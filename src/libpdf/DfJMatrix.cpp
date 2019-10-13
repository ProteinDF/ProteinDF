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

#include "DfJMatrix.h"
#include "CnError.h"
#include "DfCD.h"
#include "DfEriX.h"
#include "TlUtils.h"
#include "tl_dense_symmetric_matrix_lapack.h"

DfJMatrix::DfJMatrix(TlSerializeData* pPdfParam) : DfObject(pPdfParam) {}

DfJMatrix::~DfJMatrix() {}

void DfJMatrix::buildJ() {
    const index_type numOfAOs = this->m_nNumOfAOs;
    TlDenseSymmetricMatrix_Lapack J(numOfAOs);

    switch (this->J_engine_) {
        case J_ENGINE_RI_J:
            this->log_.info("use RI_J engine");
            this->getJ_RI();
            break;
        case J_ENGINE_CD:
            this->log_.info("use CD engine");
            this->getJ_CD();
            break;
        default:
            this->log_.info("use conventional engine");
            this->getJ_conventional();
            break;
    }
}

void DfJMatrix::getJ_RI() {
    TlDenseSymmetricMatrix_Lapack J(this->m_nNumOfAOs);
    this->getJ_RI_local(&J);
    this->saveJMatrix(J);
}

void DfJMatrix::getJ_CD() {
    TlDenseSymmetricMatrix_Lapack J(this->m_nNumOfAOs);
    this->getJ_CD_local(&J);
    this->saveJMatrix(J);
}

void DfJMatrix::getJ_conventional() {
    TlDenseSymmetricMatrix_Lapack J(this->m_nNumOfAOs);
    this->getJ_conventional_local(&J);
    this->saveJMatrix(J);
}

void DfJMatrix::saveJMatrix(const TlDenseSymmetricMatrix_Lapack& J) {
    DfObject::saveJMatrix(this->m_nIteration, J);
}

void DfJMatrix::getJ_RI_local(TlDenseSymmetricMatrix_Lapack* pJ) {
    TlDenseVector_Lapack Rho;
    switch (this->m_nMethodType) {
        case METHOD_RKS:
            Rho = DfObject::getRho<TlDenseVector_Lapack>(RUN_RKS,
                                                         this->m_nIteration);
            break;
        case METHOD_UKS:
            Rho = DfObject::getRho<TlDenseVector_Lapack>(RUN_UKS_ALPHA,
                                                         this->m_nIteration);
            Rho += DfObject::getRho<TlDenseVector_Lapack>(RUN_UKS_BETA,
                                                          this->m_nIteration);
            break;
        default:
            this->log_.critical("this method is unsopported. sorry.");
            CnErr.abort();
            break;
    }

    if (this->isUpdateMethod_ == true) {
        if (this->m_nIteration > 1) {
            TlDenseVector_Lapack prevRho;
            switch (this->m_nMethodType) {
                case METHOD_RKS:
                    prevRho = DfObject::getRho<TlDenseVector_Lapack>(
                        RUN_RKS, this->m_nIteration - 1);
                    break;
                case METHOD_UKS:
                    prevRho = DfObject::getRho<TlDenseVector_Lapack>(
                        RUN_UKS_ALPHA, this->m_nIteration - 1);
                    prevRho += DfObject::getRho<TlDenseVector_Lapack>(
                        RUN_UKS_BETA, this->m_nIteration - 1);
                    break;
                default:
                    this->log_.critical("this method is unsopported. sorry.");
                    CnErr.abort();
                    break;
            }

            Rho -= prevRho;
        }
    }

    DfEriX dfEri(this->pPdfParam_);
    dfEri.getJ(Rho, pJ);

    if (this->isUpdateMethod_ == true) {
        if (this->m_nIteration > 1) {
            const TlDenseSymmetricMatrix_Lapack prevJ =
                DfObject::getJMatrix<TlDenseSymmetricMatrix_Lapack>(
                    this->m_nIteration - 1);
            *pJ += prevJ;
        }
    }
}

TlDenseVector_Lapack DfJMatrix::getRho(const RUN_TYPE runType,
                                       const int iteration) {
    TlDenseVector_Lapack rho =
        DfObject::getRho<TlDenseVector_Lapack>(runType, iteration);
    return rho;
}

void DfJMatrix::getJ_CD_local(TlDenseSymmetricMatrix_Lapack* pJ) {
    DfCD dfCD(this->pPdfParam_);
    dfCD.getJ(pJ);
}

void DfJMatrix::getJ_conventional_local(TlDenseSymmetricMatrix_Lapack* pJ) {
    TlDenseSymmetricMatrix_Lapack P;
    if (this->isUpdateMethod_ == true) {
        P = this->getDiffDensityMatrix<TlDenseSymmetricMatrix_Lapack>();
    } else {
        P = this->getDensityMatrix<TlDenseSymmetricMatrix_Lapack>();
    }
    // P.save("P.mat");

    DfEriX dfEri(this->pPdfParam_);
    dfEri.getJpq(P, pJ);

    if (this->isUpdateMethod_ == true) {
        if (this->m_nIteration > 1) {
            const TlDenseSymmetricMatrix_Lapack prevJ =
                this->getJMatrix(this->m_nIteration - 1);
            *pJ += prevJ;
        }
    }
}

TlDenseSymmetricMatrix_Lapack DfJMatrix::getJMatrix(const int iteration) {
    const TlDenseSymmetricMatrix_Lapack J =
        DfObject::getJMatrix<TlDenseSymmetricMatrix_Lapack>(iteration);
    return J;
}
