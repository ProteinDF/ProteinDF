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

#include "CnError.h"
#include "DfJMatrix.h"
#include "DfEriX.h"
#include "DfCD.h"
#include "TlSymmetricMatrix.h"
#include "TlUtils.h"

DfJMatrix::DfJMatrix(TlSerializeData* pPdfParam)
    : DfObject(pPdfParam) {
}


DfJMatrix::~DfJMatrix()
{
}


void DfJMatrix::buildJ()
{
    const index_type numOfAOs = this->m_nNumOfAOs;
    TlSymmetricMatrix J(numOfAOs);

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


void DfJMatrix::getJ_RI()
{
    TlSymmetricMatrix J(this->m_nNumOfAOs);
    this->getJ_RI_local(&J);
    this->saveJMatrix(J);
}


void DfJMatrix::getJ_CD()
{
    TlSymmetricMatrix J(this->m_nNumOfAOs);
    this->getJ_CD_local(&J);
    this->saveJMatrix(J);
}


void DfJMatrix::getJ_conventional()
{
    TlSymmetricMatrix J(this->m_nNumOfAOs);
    this->getJ_conventional_local(&J);
    this->saveJMatrix(J);
}


void DfJMatrix::saveJMatrix(const TlSymmetricMatrix& J)
{
    DfObject::saveJMatrix(this->m_nIteration, J);
}


void DfJMatrix::getJ_RI_local(TlSymmetricMatrix* pJ)
{
    TlVector Rho;
    switch (this->m_nMethodType) {
    case METHOD_RKS:
        Rho = DfObject::getRho<TlVector>(RUN_RKS, this->m_nIteration);
        break;
    case METHOD_UKS:
        Rho = DfObject::getRho<TlVector>(RUN_UKS_ALPHA, this->m_nIteration);
        Rho += DfObject::getRho<TlVector>(RUN_UKS_BETA, this->m_nIteration);
        break;
    default:
        this->log_.critical("this method is unsopported. sorry.");
        CnErr.abort();
        break;
    }

    if (this->isUpdateMethod_ == true) {
        if (this->m_nIteration > 1) {
            TlVector prevRho;
            switch (this->m_nMethodType) {
            case METHOD_RKS:
                prevRho = DfObject::getRho<TlVector>(RUN_RKS, this->m_nIteration -1);
                break;
            case METHOD_UKS:
                prevRho = DfObject::getRho<TlVector>(RUN_UKS_ALPHA, this->m_nIteration -1);
                prevRho += DfObject::getRho<TlVector>(RUN_UKS_BETA, this->m_nIteration -1);
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
            const TlSymmetricMatrix prevJ = DfObject::getJMatrix<TlSymmetricMatrix>(this->m_nIteration -1);
            *pJ += prevJ;
        }
    }
}


TlVector DfJMatrix::getRho(const RUN_TYPE runType, const int iteration)
{
    TlVector rho = DfObject::getRho<TlVector>(runType, iteration);
    return rho;
}


void DfJMatrix::getJ_CD_local(TlSymmetricMatrix* pJ)
{
    DfCD dfCD(this->pPdfParam_);
    dfCD.getJ(pJ);
}


void DfJMatrix::getJ_conventional_local(TlSymmetricMatrix* pJ)
{
    TlSymmetricMatrix P;
    if (this->isUpdateMethod_ == true) {
        P = this->getDiffDensityMatrix<TlSymmetricMatrix>();
    } else {
        P = this->getDensityMatrix<TlSymmetricMatrix>();
    }
    // P.save("P.mat");

    DfEriX dfEri(this->pPdfParam_);
    dfEri.getJpq(P, pJ);

    if (this->isUpdateMethod_ == true) {
        if (this->m_nIteration > 1) {
            const TlSymmetricMatrix prevJ = this->getJMatrix(this->m_nIteration -1);
            *pJ += prevJ;
        }
    }
}

TlSymmetricMatrix DfJMatrix::getJMatrix(const int iteration)
{
    const TlSymmetricMatrix J = DfObject::getJMatrix<TlSymmetricMatrix>(iteration);
    return J;
}


