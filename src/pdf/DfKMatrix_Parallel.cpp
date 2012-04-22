#include <cassert>
#include "TlCommunicate.h"
#include "DfKMatrix_Parallel.h"
#include "DfEriX_Parallel.h"
#include "DfCD_Parallel.h"
#include "CnError.h"

DfKMatrix_Parallel::DfKMatrix_Parallel(TlSerializeData* pPdfParam)
    : DfKMatrix(pPdfParam) {
}


DfKMatrix_Parallel::~DfKMatrix_Parallel()
{
}


void DfKMatrix_Parallel::getK_CD()
{
#ifdef HAVE_SCALAPACK
    if (this->m_bUsingSCALAPACK == true) {
        TlDistributeSymmetricMatrix K(this->m_nNumOfAOs);
        this->getK_CD_distributed(RUN_RKS, &K);
        DfObject::saveHFxMatrix(RUN_RKS, this->m_nIteration, K);
    } else {
        TlSymmetricMatrix K(this->m_nNumOfAOs);
        this->getK_CD_local(RUN_RKS, &K);
        this->saveKMatrix(RUN_RKS, K);
    }
#else
    {
        TlSymmetricMatrix K(this->m_nNumOfAOs);
        this->getK_CD_local(&K);
        this->saveKMatrix(K);
    }
#endif // HAVE_SCALAPACK
}


void DfKMatrix_Parallel::getK_conventional()
{
#ifdef HAVE_SCALAPACK
    if (this->m_bUsingSCALAPACK == true) {
        TlDistributeSymmetricMatrix K(this->m_nNumOfAOs);
        this->getK_conventional_distributed(RUN_RKS, &K);
        DfObject::saveHFxMatrix(RUN_RKS, this->m_nIteration, K);
    } else {
        TlSymmetricMatrix K(this->m_nNumOfAOs);
        this->getK_conventional_local(RUN_RKS, &K);
        this->saveKMatrix(RUN_RKS, K);
    }
#else
    {
        TlSymmetricMatrix K(this->m_nNumOfAOs);
        this->getK_conventional_local(&K);
        this->saveKMatrix(K);
    }
#endif // HAVE_SCALAPACK
}


void DfKMatrix_Parallel::saveKMatrix(const RUN_TYPE runType,
                                     const TlSymmetricMatrix& K)
{
    TlCommunicate& rComm = TlCommunicate::getInstance();
    if (rComm.isMaster() == true) {
        DfKMatrix::saveKMatrix(runType, K);
    }
}


void DfKMatrix_Parallel::getK_CD_local(const RUN_TYPE runType, 
                                       TlSymmetricMatrix *pK)
{
    DfCD_Parallel dfCD(this->pPdfParam_);
    dfCD.getK(runType, pK);
}


void DfKMatrix_Parallel::getK_CD_distributed(const RUN_TYPE runType,
                                             TlDistributeSymmetricMatrix *pK)
{
    DfCD_Parallel dfCD(this->pPdfParam_);
    dfCD.getK_distributed(runType, pK);
}


void DfKMatrix_Parallel::getK_conventional_local(const RUN_TYPE runType,
                                                 TlSymmetricMatrix *pK)
{
    TlCommunicate& rComm = TlCommunicate::getInstance();
    
    TlSymmetricMatrix P;
    if (this->isUpdateMethod_ == true) {
        if (rComm.isMaster() == true) {
            P = this->getDiffDensityMatrix(runType);
        }
    } else {
        if (rComm.isMaster() == true) {
            P = this->getPMatrix(runType, this->m_nIteration -1);
        }
    }
    rComm.broadcast(P);
    assert(P.getNumOfRows() == this->m_nNumOfAOs);
    
    DfEriX_Parallel dfEri(this->pPdfParam_);
    dfEri.getK(P, pK);

    if (this->isUpdateMethod_ == true) {
        if (this->m_nIteration > 1) {
            if (rComm.isMaster() == true) {
                const TlSymmetricMatrix prevK = this->getKMatrix(runType, this->m_nIteration -1);
                *pK += prevK;
            }
        }
    }
}


void DfKMatrix_Parallel::getK_conventional_distributed(const RUN_TYPE runType,
                                                       TlDistributeSymmetricMatrix *pK)
{
    TlCommunicate& rComm = TlCommunicate::getInstance();
    TlDistributeSymmetricMatrix P;
    if (this->isUpdateMethod_ == true) {
        P = DfObject::getDiffDensityMatrix<TlDistributeSymmetricMatrix>(runType, this->m_nIteration);
    } else {
        P = DfObject::getPpqMatrix<TlDistributeSymmetricMatrix>(runType, this->m_nIteration -1);
    }
    assert(P.getNumOfRows() == this->m_nNumOfAOs);

    DfEriX_Parallel dfEri(this->pPdfParam_);
    dfEri.getK_D(P, pK);
    
    if (this->isUpdateMethod_ == true) {
        if (this->m_nIteration > 1) {
            if (rComm.isMaster() == true) {
                const TlDistributeSymmetricMatrix prevK = 
                    DfObject::getHFxMatrix<TlDistributeSymmetricMatrix>(runType,
                                                                        this->m_nIteration -1);
                *pK += prevK;
            }
        }
    }
}


TlSymmetricMatrix DfKMatrix_Parallel::getPMatrix(const RUN_TYPE runType,
                                                 const int iteration)
{
    TlCommunicate& rComm = TlCommunicate::getInstance();

    TlSymmetricMatrix P;
    if (rComm.isMaster() == true) {
        P = DfKMatrix::getPMatrix(runType, iteration);
    }
    rComm.broadcast(P);

    return P;
}


TlSymmetricMatrix DfKMatrix_Parallel::getDiffDensityMatrix(const RUN_TYPE runType)
{
    TlCommunicate& rComm = TlCommunicate::getInstance();

    TlSymmetricMatrix diffP;
    if (rComm.isMaster() == true) {
        diffP = DfKMatrix::getDiffDensityMatrix(runType);
    }
    rComm.broadcast(diffP);

    return diffP;
}


TlSymmetricMatrix DfKMatrix_Parallel::getKMatrix(const RUN_TYPE runType,
                                                 const int iteration)
{
    TlCommunicate& rComm = TlCommunicate::getInstance();

    TlSymmetricMatrix K;
    if (rComm.isMaster() == true) {
        K = this->getKMatrix(runType, iteration);
    }
    rComm.broadcast(K);

    return K;
}

