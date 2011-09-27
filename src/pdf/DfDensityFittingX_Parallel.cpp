#include <cassert>

#include "DfDensityFittingX_Parallel.h"
#include "TlCommunicate.h"

#include "TlSymmetricMatrix.h"
#include "TlDistributeSymmetricMatrix.h"
#include "TlDistributeVector.h"
#include "TlFile.h"

// LAPACK ==============================================================
DfDensityFittingX_Parallel::DfDensityFittingX_Parallel(TlSerializeData* pPdfParam)
    : DfDensityFittingTmpl<TlSymmetricMatrix, TlVector, DfEriX_Parallel>(pPdfParam)
{
}


DfDensityFittingX_Parallel::~DfDensityFittingX_Parallel()
{
}


void DfDensityFittingX_Parallel::exec()
{
    this->calc();
}


void DfDensityFittingX_Parallel::logger(const std::string& str) const
{
    TlCommunicate& rComm = TlCommunicate::getInstance();

    if (rComm.isMaster() == true) {
        DfObject::logger(str);
    }

    rComm.barrier();
}


TlVector DfDensityFittingX_Parallel::getNalpha()
{
    TlCommunicate& rComm = TlCommunicate::getInstance();

    TlVector Nalpha;
    if (rComm.isMaster() == true) {
        Nalpha = DfDensityFittingTmpl<TlSymmetricMatrix, TlVector, DfEriX_Parallel>::getNalpha();
    }
    rComm.broadcast(Nalpha);

    return Nalpha;
}

TlSymmetricMatrix DfDensityFittingX_Parallel::getSabinv()
{
    TlCommunicate& rComm = TlCommunicate::getInstance();

    TlSymmetricMatrix Sabinv;
    if (rComm.isMaster() == true) {
        Sabinv = DfDensityFittingTmpl<TlSymmetricMatrix, TlVector, DfEriX_Parallel>::getSabinv();
    }
    rComm.broadcast(Sabinv);

    return Sabinv;
}


TlVector DfDensityFittingX_Parallel::calcTAlpha_DIRECT(const TlSymmetricMatrix& P)
{
    return DfDensityFittingTmpl<TlSymmetricMatrix, TlVector, DfEriX_Parallel>::calcTAlpha_DIRECT(P);
}


TlVector DfDensityFittingX_Parallel::get_flVctTalpha(const int nIteration, const std::string& type)
{
    TlCommunicate& rComm = TlCommunicate::getInstance();

    TlVector flVctTalpha;
    if (rComm.isMaster() == true) {
        flVctTalpha = DfDensityFittingTmpl<TlSymmetricMatrix, TlVector, DfEriX_Parallel>::get_flVctTalpha(nIteration, type);
    }
    rComm.broadcast(flVctTalpha);

    return flVctTalpha;
}


void DfDensityFittingX_Parallel::getTalpha_ROKS(TlVector* pT_alphaA, TlVector* pT_alphaB)
{
    TlCommunicate& rComm = TlCommunicate::getInstance();

    if (rComm.isMaster() == true) {
        DfDensityFittingTmpl<TlSymmetricMatrix, TlVector, DfEriX_Parallel>::getTalpha_ROKS(pT_alphaA, pT_alphaB);
    }

    rComm.broadcast(*pT_alphaA);
    rComm.broadcast(*pT_alphaB);
}


TlSymmetricMatrix DfDensityFittingX_Parallel::getDiffDensityMatrix(RUN_TYPE runType)
{
    TlCommunicate& rComm = TlCommunicate::getInstance();

    TlSymmetricMatrix diffP;
    if (rComm.isMaster() == true) {
        diffP = DfObject::getDiffDensityMatrix<TlSymmetricMatrix>(runType, this->m_nIteration);
    }
    rComm.broadcast(diffP);

    return diffP;
}


TlSymmetricMatrix DfDensityFittingX_Parallel::getP1pq(const int nIteration)
{
    TlCommunicate& rComm = TlCommunicate::getInstance();

    TlSymmetricMatrix P;
    if (rComm.isMaster() == true) {
        P = DfDensityFittingTmpl<TlSymmetricMatrix, TlVector, DfEriX_Parallel>::getP1pq(nIteration);
    }
    rComm.broadcast(P);

    return P;
}


TlSymmetricMatrix DfDensityFittingX_Parallel::getP2pq(const int nIteration)
{
    TlCommunicate& rComm = TlCommunicate::getInstance();

    TlSymmetricMatrix P;
    if (rComm.isMaster() == true) {
        P = DfDensityFittingTmpl<TlSymmetricMatrix, TlVector, DfEriX_Parallel>::getP2pq(nIteration);
    }
    rComm.broadcast(P);

    return P;
}


double DfDensityFittingX_Parallel::getLamda(const TlVector& SabinvN, const TlVector& t_alpha,
                                           const TlVector& N_alpha, const double dNumOfElec)
{
    TlCommunicate& rComm = TlCommunicate::getInstance();

    double dAnswer;
    if (rComm.isMaster() == true) {
        dAnswer = DfDensityFittingTmpl<TlSymmetricMatrix, TlVector, DfEriX_Parallel>::getLamda(SabinvN, t_alpha, N_alpha, dNumOfElec);
    }
    rComm.broadcast(dAnswer);

    return dAnswer;
}

void DfDensityFittingX_Parallel::saveRho(const TlVector& rRho, const RUN_TYPE runType)
{
    TlCommunicate& rComm = TlCommunicate::getInstance();

    if (rComm.isMaster() == true) {
        DfDensityFittingTmpl<TlSymmetricMatrix, TlVector, DfEriX_Parallel>::saveRho(rRho, runType);
    }
    rComm.barrier();
}

