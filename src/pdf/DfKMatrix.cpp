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
    switch (this->K_engine_) {
    // case K_ENGINE_RI_K:
    //     this->log_.info("use RI_K engine");
    //     this->getK_RI();
    //     break;
    case K_ENGINE_CD:
        this->log_.info("use CD engine");
        this->getK_CD();
        break;
    default:
        this->log_.info("use conventional engine");
        this->getK_conventional();
        break;
    }
}


// void DfKMatrix::getK_RI()
// {
//     TlSymmetricMatrix K(this->m_nNumOfAOs);
//     this->getK_RI_local(RUN_RKS, &K);
//     this->saveKMatrix(RUN_RKS, K);
// }


void DfKMatrix::getK_CD()
{
    TlSymmetricMatrix K(this->m_nNumOfAOs);
    this->getK_CD_local(RUN_RKS, &K);
    this->saveKMatrix(RUN_RKS, K);
}


void DfKMatrix::getK_conventional()
{
    TlSymmetricMatrix K(this->m_nNumOfAOs);
    this->getK_conventional_local(RUN_RKS, &K);
    this->saveKMatrix(RUN_RKS, K);
}


// void DfKMatrix::getK_RI_local(const RUN_TYPE runType,
//                               TlSymmetricMatrix *pK)
// {
//     TlSymmetricMatrix P;
//     if (this->isUpdateMethod_ == true) {
//         P = this->getDiffDensityMatrix(runType);
//     } else {
//         P = this->getPMatrix(runType, this->m_nIteration -1);
//     }

//     DfEri2 dfEri2(this->pPdfParam_);
//     *pK = dfEri2.getKMatrix(P);

//     if (this->isUpdateMethod_ == true) {
//         if (this->m_nIteration > 1) {
//             const TlSymmetricMatrix prevK = DfObject::getHFxMatrix<TlSymmetricMatrix>(runType, this->m_nIteration -1);
//             *pK += prevK;
//         }
//     }
// }

void DfKMatrix::getK_CD_local(const RUN_TYPE runType,
                              TlSymmetricMatrix *pK)
{
    DfCD dfCD(this->pPdfParam_);
    dfCD.getK(runType, pK);
}


void DfKMatrix::getK_conventional_local(const RUN_TYPE runType,
                                        TlSymmetricMatrix *pK)
{
    TlSymmetricMatrix P;
    if (this->isUpdateMethod_ == true) {
        P = this->getDiffDensityMatrix(runType);
    } else {
        P = this->getPMatrix(runType, this->m_nIteration -1);
    }

    DfEriX dfEri(this->pPdfParam_);
    dfEri.getK(P, pK);

    if (this->isUpdateMethod_ == true) {
        if (this->m_nIteration > 1) {
            const TlSymmetricMatrix prevK = this->getKMatrix(runType, this->m_nIteration -1);
            *pK += prevK;
        }
    }
}


TlSymmetricMatrix DfKMatrix::getPMatrix(const RUN_TYPE runType,
                                        const int iteration)
{
    const TlSymmetricMatrix P = DfObject::getPpqMatrix<TlSymmetricMatrix>(runType, iteration);
    return P;
}


TlSymmetricMatrix DfKMatrix::getDiffDensityMatrix(const RUN_TYPE runType)
{
    const TlSymmetricMatrix diffP = 
        DfObject::getDiffDensityMatrix<TlSymmetricMatrix>(runType, this->m_nIteration);
    return diffP;
}


TlSymmetricMatrix DfKMatrix::getKMatrix(const RUN_TYPE runType,
                                        const int iteration)
{
    const TlSymmetricMatrix K = DfObject::getHFxMatrix<TlSymmetricMatrix>(runType, iteration);
    return K;
}


void DfKMatrix::saveKMatrix(const RUN_TYPE runType,
                            const TlSymmetricMatrix& K)
{
    DfObject::saveHFxMatrix(runType, this->m_nIteration, K);
}


