#include <cassert>

#include "DfDensityFitting_Parallel.h"
#include "TlCommunicate.h"

#include "TlSymmetricMatrix.h"
#include "TlDistributeSymmetricMatrix.h"
#include "TlDistributeVector.h"

#include "DfDensityFitting_Parallel.h"
#include "TlCommunicate.h"
#include "TlFile.h"

// LAPACK ==============================================================
DfDensityFitting_Parallel::DfDensityFitting_Parallel(TlSerializeData* pPdfParam)
    : DfDensityFittingTmpl<TlSymmetricMatrix, TlVector, DfEri_Parallel>(pPdfParam)
{
}


DfDensityFitting_Parallel::~DfDensityFitting_Parallel()
{
}


void DfDensityFitting_Parallel::exec()
{
    this->calc();
}


void DfDensityFitting_Parallel::logger(const std::string& str) const
{
    TlCommunicate& rComm = TlCommunicate::getInstance();

    if (rComm.isMaster() == true) {
        DfObject::logger(str);
    }

    rComm.barrier();
}


TlVector DfDensityFitting_Parallel::getNalpha()
{
    TlCommunicate& rComm = TlCommunicate::getInstance();

    TlVector Nalpha;
    if (rComm.isMaster() == true) {
        Nalpha = DfDensityFittingTmpl<TlSymmetricMatrix, TlVector, DfEri_Parallel>::getNalpha();
    }
    rComm.broadcast(Nalpha);

    return Nalpha;
}

TlSymmetricMatrix DfDensityFitting_Parallel::getSabinv()
{
    TlCommunicate& rComm = TlCommunicate::getInstance();

    TlSymmetricMatrix Sabinv;
    if (rComm.isMaster() == true) {
        Sabinv = DfDensityFittingTmpl<TlSymmetricMatrix, TlVector, DfEri_Parallel>::getSabinv();
    }
    rComm.broadcast(Sabinv);

    return Sabinv;
}


TlVector DfDensityFitting_Parallel::calcTAlpha_DIRECT(const TlSymmetricMatrix& P)
{
    return DfDensityFittingTmpl<TlSymmetricMatrix, TlVector, DfEri_Parallel>::calcTAlpha_DIRECT(P);
}


TlVector DfDensityFitting_Parallel::getTalpha(const RUN_TYPE runType,
                                              const int iteration)
{
    TlCommunicate& rComm = TlCommunicate::getInstance();

    TlVector flVctTalpha;
    if (rComm.isMaster() == true) {
        flVctTalpha = DfDensityFittingTmpl<TlSymmetricMatrix, TlVector, DfEri_Parallel>::getTalpha(runType, iteration);
    }
    rComm.broadcast(flVctTalpha);

    return flVctTalpha;
}


void DfDensityFitting_Parallel::getTalpha_ROKS(TlVector* pT_alphaA, TlVector* pT_alphaB)
{
    TlCommunicate& rComm = TlCommunicate::getInstance();

    if (rComm.isMaster() == true) {
        DfDensityFittingTmpl<TlSymmetricMatrix, TlVector, DfEri_Parallel>::getTalpha_ROKS(pT_alphaA, pT_alphaB);
    }

    rComm.broadcast(*pT_alphaA);
    rComm.broadcast(*pT_alphaB);
}


TlSymmetricMatrix DfDensityFitting_Parallel::getDiffDensityMatrix(RUN_TYPE runType)
{
    TlCommunicate& rComm = TlCommunicate::getInstance();

    TlSymmetricMatrix diffP;
    if (rComm.isMaster() == true) {
        diffP = DfObject::getDiffDensityMatrix<TlSymmetricMatrix>(runType, this->m_nIteration);
    }
    rComm.broadcast(diffP);

    return diffP;
}


TlSymmetricMatrix DfDensityFitting_Parallel::getP1pq(const int nIteration)
{
    TlCommunicate& rComm = TlCommunicate::getInstance();

    TlSymmetricMatrix P;
    if (rComm.isMaster() == true) {
        P = DfDensityFittingTmpl<TlSymmetricMatrix, TlVector, DfEri_Parallel>::getP1pq(nIteration);
    }
    rComm.broadcast(P);

    return P;
}


TlSymmetricMatrix DfDensityFitting_Parallel::getP2pq(const int nIteration)
{
    TlCommunicate& rComm = TlCommunicate::getInstance();

    TlSymmetricMatrix P;
    if (rComm.isMaster() == true) {
        P = DfDensityFittingTmpl<TlSymmetricMatrix, TlVector, DfEri_Parallel>::getP2pq(nIteration);
    }
    rComm.broadcast(P);

    return P;
}


double DfDensityFitting_Parallel::getLamda(const TlVector& SabinvN, const TlVector& t_alpha,
                                           const TlVector& N_alpha, const double dNumOfElec)
{
    TlCommunicate& rComm = TlCommunicate::getInstance();

    double dAnswer;
    if (rComm.isMaster() == true) {
        dAnswer = DfDensityFittingTmpl<TlSymmetricMatrix, TlVector, DfEri_Parallel>::getLamda(SabinvN, t_alpha, N_alpha, dNumOfElec);
    }
    rComm.broadcast(dAnswer);

    return dAnswer;
}

void DfDensityFitting_Parallel::saveRho(const TlVector& rRho, const RUN_TYPE runType)
{
    TlCommunicate& rComm = TlCommunicate::getInstance();

    if (rComm.isMaster() == true) {
        DfDensityFittingTmpl<TlSymmetricMatrix, TlVector, DfEri_Parallel>::saveRho(rRho, runType);
    }
    rComm.barrier();
}

// SCALAPACK ===========================================================
DfDensityFitting_ScaLAPACK::DfDensityFitting_ScaLAPACK(TlSerializeData* pPdfParam)
    : DfDensityFittingTmpl<TlDistributeSymmetricMatrix, TlDistributeVector, DfEri_Parallel>(pPdfParam)
{
}


DfDensityFitting_ScaLAPACK::~DfDensityFitting_ScaLAPACK()
{
}


void DfDensityFitting_ScaLAPACK::exec()
{
    this->calc();
}


void DfDensityFitting_ScaLAPACK::logger(const std::string& str) const
{
    TlCommunicate& rComm = TlCommunicate::getInstance();

    if (rComm.isMaster() == true) {
        DfObject::logger(str);
    }
}


template<>
TlDistributeVector
DfDensityFittingTmpl<TlDistributeSymmetricMatrix, TlDistributeVector, DfEri_Parallel>::getTalpha(const RUN_TYPE runType)
{
    TlDistributeVector t_alpha;
    // std::string suffix = "";
    // if (runType == RUN_UKS_ALPHA) {
    //     suffix = "a";
    // } else if (runType == RUN_UKS_BETA) {
    //     suffix = "b";
    // }

    if (this->m_bDiskUtilization == false) {
        TlDistributeSymmetricMatrix P = DfObject::getDiffDensityMatrix<TlDistributeSymmetricMatrix>(runType, this->m_nIteration);
//         if (this->isSaveDistributedMatrixToLocalDisk_ == true) {
//             this->logger(" load distributed density matrix from local disk.\n");
//             P.loadLocal(this->getDiffDensityMatrixPath(runType));
//         } else {
//             P = this->getDiffDensityMatrix(runType);
//         }
        t_alpha = this->calcTAlpha_DIRECT(P);
        t_alpha += this->getTalpha(runType, this->m_nIteration -1);
    } else {
        abort();
    }

    // save
    if (this->m_bDiskUtilization == false) {
        t_alpha.save(this->getTalphaPath(runType, this->m_nIteration));
    }

    return t_alpha;
}

