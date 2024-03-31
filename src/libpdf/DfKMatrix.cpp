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

#include "DfKMatrix.h"

#include "CnError.h"
#include "DfCD.h"
#include "DfEriX.h"
#include "common.h"
#include "df_cdk_matrix.h"

DfKMatrix::DfKMatrix(TlSerializeData* pPdfParam)
    : DfObject(pPdfParam) {
    this->updateLinearAlgebraPackageParam(
        (*(this->pPdfParam_))["linear_algebra_package/K"].getStr());
}

DfKMatrix::~DfKMatrix() {}

void DfKMatrix::buildK() {
    switch (this->K_engine_) {
        // case K_ENGINE_RI_K:
        //     this->log_.info("use RI_K engine");
        //     this->getK_RI();
        //     break;
        case K_ENGINE_CD:
        case K_ENGINE_FASTCDK:
            this->log_.info("use CD method");
            this->getK_CD();
            break;
        default:
            this->log_.info("use conventional method");
            this->getK_conventional();
            break;
    }
}

// void DfKMatrix::getK_RI()
// {
//     TlDenseSymmetricMatrix_Lapack K(this->m_nNumOfAOs);
//     this->getK_RI_local(RUN_RKS, &K);
//     this->saveKMatrix(RUN_RKS, K);
// }

void DfKMatrix::getK_CD() {
    DfCdkMatrix dfCDK(this->pPdfParam_);
    dfCDK.getK();
}

void DfKMatrix::getK_conventional() {
    switch (this->m_nMethodType) {
        case METHOD_RKS: {
            TlDenseSymmetricMatrix_Lapack K(this->m_nNumOfAOs);
            this->getK_conventional_local(RUN_RKS, &K);
            this->saveKMatrix(RUN_RKS, K);
        } break;

        case METHOD_UKS: {
            TlDenseSymmetricMatrix_Lapack K(this->m_nNumOfAOs);
            this->getK_conventional_local(RUN_UKS_ALPHA, &K);
            this->saveKMatrix(RUN_UKS_ALPHA, K);
        }
            {
                TlDenseSymmetricMatrix_Lapack K(this->m_nNumOfAOs);
                this->getK_conventional_local(RUN_UKS_BETA, &K);
                this->saveKMatrix(RUN_UKS_BETA, K);
            }
            break;

        case METHOD_ROKS: {
            TlDenseSymmetricMatrix_Lapack K(this->m_nNumOfAOs);
            this->getK_conventional_local(RUN_ROKS_ALPHA, &K);
            this->saveKMatrix(RUN_ROKS_ALPHA, K);
        }
            {
                TlDenseSymmetricMatrix_Lapack K(this->m_nNumOfAOs);
                this->getK_conventional_local(RUN_ROKS_BETA, &K);
                this->saveKMatrix(RUN_ROKS_BETA, K);
            }
            break;

        default:
            this->log_.critical("program error.");
            CnErr.abort();
            break;
    }
}

// void DfKMatrix::getK_RI_local(const RUN_TYPE runType,
//                               TlDenseSymmetricMatrix_Lapack *pK)
// {
//     TlDenseSymmetricMatrix_Lapack P;
//     if (this->isUpdateMethod_ == true) {
//         P = this->getDiffDensityMatrix(runType);
//     } else {
//         P = this->getPMatrix(runType, this->m_nIteration -1);
//     }

//     DfEri2 dfEri2(this->pPdfParam_);
//     *pK = dfEri2.getKMatrix(P);

//     if (this->isUpdateMethod_ == true) {
//         if (this->m_nIteration > 1) {
//             const TlDenseSymmetricMatrix_Lapack prevK =
//             DfObject::getHFxMatrix<TlDenseSymmetricMatrix_Lapack>(runType,
//             this->m_nIteration -1); *pK += prevK;
//         }
//     }
// }

void DfKMatrix::getK_CD_replica(const RUN_TYPE runType) {
    DfCD dfCD(this->pPdfParam_);
    dfCD.getK(runType);
}

void DfKMatrix::getK_conventional_local(const RUN_TYPE runType,
                                        TlDenseSymmetricMatrix_Lapack* pK) {
    TlDenseSymmetricMatrix_Lapack P;
    if (this->isUpdateMethod_ == true) {
        P = this->getDiffDensityMatrix<TlDenseSymmetricMatrix_Lapack>(runType);
    } else {
        P = this->getDensityMatrix<TlDenseSymmetricMatrix_Lapack>(runType);
    }

    DfEriX dfEri(this->pPdfParam_);
    dfEri.getK(P, pK);

    if (this->isUpdateMethod_ == true) {
        if (this->m_nIteration > 1) {
            const TlDenseSymmetricMatrix_Lapack prevK =
                this->getKMatrix(runType, this->m_nIteration - 1);
            *pK += prevK;
        }
    }
}

TlDenseSymmetricMatrix_Lapack DfKMatrix::getKMatrix(const RUN_TYPE runType,
                                                    const int iteration) {
    const TlDenseSymmetricMatrix_Lapack K =
        DfObject::getHFxMatrix<TlDenseSymmetricMatrix_Lapack>(runType,
                                                              iteration);
    return K;
}

void DfKMatrix::saveKMatrix(const RUN_TYPE runType,
                            const TlDenseSymmetricMatrix_Lapack& K) {
    DfObject::saveHFxMatrix(runType, this->m_nIteration, K);
}
