#include "DfKMatrix.h"
#include "DfEriX.h"
#include "DfEri2.h"
#include "DfCD.h"

DfKMatrix::DfKMatrix(TlSerializeData* pPdfParam)
    : DfObject(pPdfParam) {
}


DfKMatrix::~DfKMatrix()
{
}


void DfKMatrix::buildK()
{
    const index_type numOfAOs = this->m_nNumOfAOs;
    TlSymmetricMatrix K(numOfAOs);

    switch (this->K_engine_) {
    case K_ENGINE_RI_K:
        this->log_.info("use RI_K engine");
        this->getK_byRI_K(RUN_RKS, &K);
        break;
    case K_ENGINE_CD:
        this->log_.info("use CD engine");
        this->getK_byCD(RUN_RKS, &K);
        break;
    default:
        this->log_.info("use conventional engine");
        this->getK_conventional(RUN_RKS, &K);
        break;
    }

    DfObject::saveHFxMatrix(RUN_RKS, this->m_nIteration, K);
}

void DfKMatrix::getK_byRI_K(const RUN_TYPE runType,
                            TlSymmetricMatrix *pK)
{
    TlSymmetricMatrix P;
    if (this->isUpdateMethod_ == true) {
        P = this->getDiffDensityMatrix<TlSymmetricMatrix>(runType, this->m_nIteration);
    } else {
        P = this->getPpqMatrix<TlSymmetricMatrix>(runType, this->m_nIteration -1);
    }

    DfEri2 dfEri2(this->pPdfParam_);
    *pK = dfEri2.getKMatrix(P);

    if (this->isUpdateMethod_ == true) {
        if (this->m_nIteration > 1) {
            const TlSymmetricMatrix prevK = DfObject::getHFxMatrix<TlSymmetricMatrix>(runType, this->m_nIteration -1);
            *pK += prevK;
        }
    }
}

void DfKMatrix::getK_conventional(const RUN_TYPE runType,
                                  TlSymmetricMatrix *pK)
{
    TlSymmetricMatrix P;
    if (this->isUpdateMethod_ == true) {
        P = this->getDiffDensityMatrix<TlSymmetricMatrix>(runType, this->m_nIteration);
    } else {
        P = this->getPpqMatrix<TlSymmetricMatrix>(runType, this->m_nIteration -1);
    }

    DfEriX dfEri(this->pPdfParam_);
    dfEri.getK(P, pK);

    if (this->isUpdateMethod_ == true) {
        if (this->m_nIteration > 1) {
            const TlSymmetricMatrix prevK = DfObject::getHFxMatrix<TlSymmetricMatrix>(runType, this->m_nIteration -1);
            *pK += prevK;
        }
    }
}


void DfKMatrix::getK_byCD(const RUN_TYPE runType,
                          TlSymmetricMatrix *pK)
{
    DfCD dfCD(this->pPdfParam_);
    dfCD.getK(runType, pK);
}


