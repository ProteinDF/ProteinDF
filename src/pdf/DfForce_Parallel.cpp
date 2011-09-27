#include "DfForce_Parallel.h"
#include "TlCommunicate.h"
#include "TlSymmetricMatrix.h"
#include "TlVector.h"
#include "DfEriX_Parallel.h"
#include "DfXCFunctional_Parallel.h"


DfForce_Parallel::DfForce_Parallel(TlSerializeData* pPdfParam)
    : DfForce(pPdfParam) {
}


DfForce_Parallel::~DfForce_Parallel()
{
}


void DfForce_Parallel::logger(const std::string& str) const
{
    TlCommunicate& rComm = TlCommunicate::getInstance();
    if (rComm.isMaster() == true) {
        DfForce::logger(str);
    }
}


void DfForce_Parallel::calcForceFromCoulomb_RIJ(const RUN_TYPE runType)
{
    this->calcForceFromCoulomb_RIJ_DC(runType);
}


void DfForce_Parallel::calcForceFromCoulomb_RIJ_DC(const RUN_TYPE runType)
{
    this->loggerTime(" calc force from J (RIJ, DC)");

    TlCommunicate& rComm = TlCommunicate::getInstance();
    
    const int iteration = this->m_nIteration;
    const int numOfAtoms = this->m_nNumOfAtoms;
    
    TlVector rho;
    if (rComm.isMaster() == true) {
        rho = this->getRho<TlVector>(runType, iteration);
    }
    rComm.broadcast(rho);

    TlSymmetricMatrix P;
    if (rComm.isMaster() == true) {
        P = this->getPpqMatrix<TlSymmetricMatrix>(runType, iteration);
    }
    rComm.broadcast(P);

    DfEriX_Parallel dfEri(this->pPdfParam_);

    // ((pq)'|a)
    TlMatrix F_pqa(numOfAtoms, 3);
    dfEri.getForceJ(P, rho, &F_pqa);

    // (a'|b)
    TlMatrix F_ab(numOfAtoms, 3);
    dfEri.getForceJ(rho, &F_ab);

    if (this->isDebugOutMatrix_ == true) {
        if (rComm.isMaster() == true) {
            F_pqa.save("F_pqa.mtx");
            F_ab.save("F_ab.mtx");
        }
    }

    const TlMatrix F_J = (F_pqa - 0.5 * F_ab);
    this->force_ += F_J;
}


void DfForce_Parallel::calcForceFromK(const RUN_TYPE runType)
{
    this->calcForceFromK_DC(runType);
}


void DfForce_Parallel::calcForceFromK_DC(const RUN_TYPE runType)
{
    this->loggerTime(" calc force from K (RIJ, DC)");

    TlCommunicate& rComm = TlCommunicate::getInstance();
    const DfXCFunctional_Parallel dfXCFunctional(this->pPdfParam_);

    if (dfXCFunctional.isHybridFunctional() == true) {
        const int iteration = this->m_nIteration;
        const int numOfAtoms = this->m_nNumOfAtoms;
        
        TlSymmetricMatrix P;
        if (rComm.isMaster() == true) {
            P = this->getPpqMatrix<TlSymmetricMatrix>(runType, iteration);
        }
        rComm.broadcast(P);
        
        DfEriX_Parallel dfEri(this->pPdfParam_);
        
        TlMatrix F_K(numOfAtoms, 3);
        dfEri.getForceK(0.5 * P, &F_K);

        F_K *= -1.0;
        F_K *= dfXCFunctional.getFockExchangeCoefficient(); // for B3LYP

        if (this->isDebugOutMatrix_ == true) {
            if (rComm.isMaster() == true) {
                F_K.save("F_K.mtx");
            }
        }
        this->force_ += F_K;
    }
}

