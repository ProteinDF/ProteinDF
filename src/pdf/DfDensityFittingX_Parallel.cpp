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


TlVector DfDensityFittingX_Parallel::getTalpha(const RUN_TYPE runType, const int iteration)
{
    TlCommunicate& rComm = TlCommunicate::getInstance();

    TlVector flVctTalpha;
    if (rComm.isMaster() == true) {
        flVctTalpha = DfDensityFittingTmpl<TlSymmetricMatrix, TlVector, DfEriX_Parallel>::getTalpha(runType, iteration);
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

