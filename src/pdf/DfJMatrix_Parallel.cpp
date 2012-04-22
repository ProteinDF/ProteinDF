#include <cassert>
#include "TlCommunicate.h"
#include "DfJMatrix_Parallel.h"
#include "DfEriX_Parallel.h"
#include "DfCD_Parallel.h"
#include "CnError.h"

DfJMatrix_Parallel::DfJMatrix_Parallel(TlSerializeData* pPdfParam)
    : DfJMatrix(pPdfParam) {
}


DfJMatrix_Parallel::~DfJMatrix_Parallel()
{
}


void DfJMatrix_Parallel::getJ_RI()
{
#ifdef HAVE_SCALAPACK
    if (this->m_bUsingSCALAPACK == true) {
        TlDistributeSymmetricMatrix J(this->m_nNumOfAOs);
        this->getJ_RI_distributed(&J);
        DfObject::saveJMatrix(this->m_nIteration, J);
    } else {
        TlSymmetricMatrix J(this->m_nNumOfAOs);
        this->getJ_RI_local(&J);
        this->saveJMatrix(J);
    }
#else
    {
        TlSymmetricMatrix J(this->m_nNumOfAOs);
        this->getJ_RI_local(&J);
        this->saveJMatrix(J);
    }
#endif // HAVE_SCALAPACK
}


void DfJMatrix_Parallel::getJ_CD()
{
#ifdef HAVE_SCALAPACK
    if (this->m_bUsingSCALAPACK == true) {
        TlDistributeSymmetricMatrix J(this->m_nNumOfAOs);
        this->getJ_CD_distributed(&J);
        DfObject::saveJMatrix(this->m_nIteration, J);
    } else {
        TlSymmetricMatrix J(this->m_nNumOfAOs);
        this->getJ_CD_local(&J);
        this->saveJMatrix(J);
    }
#else
    {
        TlSymmetricMatrix J(this->m_nNumOfAOs);
        this->getJ_CD_local(&J);
        this->saveJMatrix(J);
    }
#endif // HAVE_SCALAPACK
}


void DfJMatrix_Parallel::getJ_conventional()
{
#ifdef HAVE_SCALAPACK
    if (this->m_bUsingSCALAPACK == true) {
        TlDistributeSymmetricMatrix J(this->m_nNumOfAOs);
        this->getJ_conventional_distributed(&J);
        DfObject::saveJMatrix(this->m_nIteration, J);
    } else {
        TlSymmetricMatrix J(this->m_nNumOfAOs);
        this->getJ_conventional_local(&J);
        this->saveJMatrix(J);
    }
#else
    {
        TlSymmetricMatrix J(this->m_nNumOfAOs);
        this->getJ_conventional_local(&J);
        this->saveJMatrix(J);
    }
#endif // HAVE_SCALAPACK
}


void DfJMatrix_Parallel::saveJMatrix(const TlSymmetricMatrix& J)
{
    TlCommunicate& rComm = TlCommunicate::getInstance();
    if (rComm.isMaster() == true) {
        DfObject::saveJMatrix(this->m_nIteration, J);
    }
}


void DfJMatrix_Parallel::getJ_RI_local(TlSymmetricMatrix *pJ)
{
    TlVector rho;
    switch (this->m_nMethodType) {
    case METHOD_RKS:
        rho = this->getRho(RUN_RKS, this->m_nIteration);
        break;
    case METHOD_UKS:
        rho  = this->getRho(RUN_UKS_ALPHA, this->m_nIteration);
        rho += this->getRho(RUN_UKS_BETA,  this->m_nIteration);
        break;
    default:
        std::cerr << "unsopported. sorry." << std::endl;
        CnErr.abort();
        break;
    }

    if (this->isUpdateMethod_ == true) {
        if (this->m_nIteration > 1) {
            TlVector prevRho;
            switch (this->m_nMethodType) {
            case METHOD_RKS:
                prevRho = this->getRho(RUN_RKS, this->m_nIteration -1);
                break;
            case METHOD_UKS:
                prevRho = this->getRho(RUN_UKS_ALPHA, this->m_nIteration -1);
                prevRho += this->getRho(RUN_UKS_BETA, this->m_nIteration -1);
                break;
            default:
                std::cerr << "unsopported. sorry." << std::endl;
                CnErr.abort();
                break;
            }

            rho -= prevRho;
        }
    }

    DfEriX dfEri(this->pPdfParam_);
    dfEri.getdeltaHpqA(rho, *pJ);
}


TlVector DfJMatrix_Parallel::getRho(const RUN_TYPE runType, const int iteration)
{
    TlCommunicate& rComm = TlCommunicate::getInstance();
    TlVector rho;
    if (rComm.isMaster() == true) {
        rho = DfObject::getRho<TlVector>(runType, iteration);
    }
    rComm.broadcast(rho);

    return rho;
}


void DfJMatrix_Parallel::getJ_RI_distributed(TlDistributeSymmetricMatrix *pJ)
{
    TlVector rho;
    switch (this->m_nMethodType) {
    case METHOD_RKS:
        rho = this->getRho(RUN_RKS, this->m_nIteration);
        break;
    case METHOD_UKS:
        rho  = this->getRho(RUN_UKS_ALPHA, this->m_nIteration);
        rho += this->getRho(RUN_UKS_BETA,  this->m_nIteration);
        break;
    default:
        std::cerr << "unsopported. sorry." << std::endl;
        CnErr.abort();
        break;
    }

    if (this->isUpdateMethod_ == true) {
        if (this->m_nIteration > 1) {
            TlVector prevRho;
            switch (this->m_nMethodType) {
            case METHOD_RKS:
                prevRho = this->getRho(RUN_RKS, this->m_nIteration -1);
                break;
            case METHOD_UKS:
                prevRho  = this->getRho(RUN_UKS_ALPHA, this->m_nIteration -1);
                prevRho += this->getRho(RUN_UKS_BETA, this->m_nIteration -1);
                break;
            default:
                std::cerr << "unsopported. sorry." << std::endl;
                CnErr.abort();
                break;
            }

            rho -= prevRho;
        }
    }

    DfEriX_Parallel dfEri(this->pPdfParam_);
    dfEri.getJ_D(rho, pJ);
}


void DfJMatrix_Parallel::getJ_CD_local(TlSymmetricMatrix *pJ)
{
    DfCD_Parallel dfCD(this->pPdfParam_);
    dfCD.getJ(pJ);
}


void DfJMatrix_Parallel::getJ_CD_distributed(TlDistributeSymmetricMatrix *pJ)
{
    DfCD_Parallel dfCD(this->pPdfParam_);
    dfCD.getJ_distributed(pJ);
}


void DfJMatrix_Parallel::getJ_conventional_local(TlSymmetricMatrix *pJ)
{
    TlCommunicate& rComm = TlCommunicate::getInstance();
    
    TlSymmetricMatrix P;
    if (this->isUpdateMethod_ == true) {
        if (rComm.isMaster() == true) {
            P = this->getDiffDensityMatrix();
        }
    } else {
        if (rComm.isMaster() == true) {
            P = this->getPMatrix(this->m_nIteration -1);
        }
    }
    rComm.broadcast(P);
    assert(P.getNumOfRows() == this->m_nNumOfAOs);
    
    DfEriX_Parallel dfEri(this->pPdfParam_);
    dfEri.getJpq(P, pJ);

    if (this->isUpdateMethod_ == true) {
        if (this->m_nIteration > 1) {
            if (rComm.isMaster() == true) {
                const TlSymmetricMatrix prevJ = this->getJMatrix(this->m_nIteration -1);
                *pJ += prevJ;
            }
        }
    }
}


void DfJMatrix_Parallel::getJ_conventional_distributed(TlDistributeSymmetricMatrix *pJ)
{
    TlCommunicate& rComm = TlCommunicate::getInstance();
    TlDistributeSymmetricMatrix P;
    if (this->isUpdateMethod_ == true) {
        P = DfObject::getDiffDensityMatrix<TlDistributeSymmetricMatrix>(RUN_RKS, this->m_nIteration);
    } else {
        P = DfObject::getPpqMatrix<TlDistributeSymmetricMatrix>(RUN_RKS, this->m_nIteration -1);
    }
    assert(P.getNumOfRows() == this->m_nNumOfAOs);

    DfEriX_Parallel dfEri(this->pPdfParam_);
    dfEri.getJpq_D(P, pJ);
    
    if (this->isUpdateMethod_ == true) {
        if (this->m_nIteration > 1) {
            if (rComm.isMaster() == true) {
                const TlDistributeSymmetricMatrix prevJ = 
                    DfObject::getJMatrix<TlDistributeSymmetricMatrix>(this->m_nIteration -1);
                *pJ += prevJ;
            }
        }
    }
}


TlSymmetricMatrix DfJMatrix_Parallel::getPMatrix(const int iteration)
{
    TlCommunicate& rComm = TlCommunicate::getInstance();

    TlSymmetricMatrix P;
    if (rComm.isMaster() == true) {
        P = DfJMatrix::getPMatrix(iteration);
    }
    rComm.broadcast(P);

    return P;
}


TlSymmetricMatrix DfJMatrix_Parallel::getDiffDensityMatrix()
{
    TlCommunicate& rComm = TlCommunicate::getInstance();

    TlSymmetricMatrix diffP;
    if (rComm.isMaster() == true) {
        diffP = DfJMatrix::getDiffDensityMatrix();
    }
    rComm.broadcast(diffP);

    return diffP;
}


TlSymmetricMatrix DfJMatrix_Parallel::getJMatrix(const int iteration)
{
    TlCommunicate& rComm = TlCommunicate::getInstance();

    TlSymmetricMatrix J;
    if (rComm.isMaster() == true) {
        J = this->getJMatrix(iteration);
    }
    rComm.broadcast(J);

    return J;
}


