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
#include "TlFile.h"
#include "tl_dense_vector_scalapack.h"

// LAPACK ==============================================================
DfDensityFittingX_Parallel::DfDensityFittingX_Parallel(
    TlSerializeData* pPdfParam)
    : DfDensityFittingTmpl<TlDenseSymmetricMatrix_Lapack, TlDenseVector_Lapack,
                           DfEriX_Parallel>(pPdfParam) {}

DfDensityFittingX_Parallel::~DfDensityFittingX_Parallel() {}

void DfDensityFittingX_Parallel::exec() { this->calc(); }

TlDenseVector_Lapack DfDensityFittingX_Parallel::getNalpha() {
    TlCommunicate& rComm = TlCommunicate::getInstance();

    TlDenseVector_Lapack Nalpha;
    if (rComm.isMaster() == true) {
        Nalpha = DfDensityFittingTmpl<TlDenseSymmetricMatrix_Lapack,
                                      TlDenseVector_Lapack,
                                      DfEriX_Parallel>::getNalpha();
    }
    rComm.broadcast(&Nalpha);

    return Nalpha;
}

TlDenseSymmetricMatrix_Lapack DfDensityFittingX_Parallel::getSabinv() {
    TlCommunicate& rComm = TlCommunicate::getInstance();

    TlDenseSymmetricMatrix_Lapack Sabinv;
    if (rComm.isMaster() == true) {
        Sabinv = DfDensityFittingTmpl<TlDenseSymmetricMatrix_Lapack,
                                      TlDenseVector_Lapack,
                                      DfEriX_Parallel>::getSabinv();
    }
    rComm.broadcast(&Sabinv);

    return Sabinv;
}

TlDenseVector_Lapack DfDensityFittingX_Parallel::calcTAlpha_DIRECT(
    const TlDenseSymmetricMatrix_Lapack& P) {
    return DfDensityFittingTmpl<TlDenseSymmetricMatrix_Lapack,
                                TlDenseVector_Lapack,
                                DfEriX_Parallel>::calcTAlpha_DIRECT(P);
}

TlDenseVector_Lapack DfDensityFittingX_Parallel::getTalpha(
    const RUN_TYPE runType, const int iteration) {
    TlCommunicate& rComm = TlCommunicate::getInstance();

    TlDenseVector_Lapack flVctTalpha;
    if (rComm.isMaster() == true) {
        flVctTalpha =
            DfDensityFittingTmpl<TlDenseSymmetricMatrix_Lapack,
                                 TlDenseVector_Lapack,
                                 DfEriX_Parallel>::getTalpha(runType,
                                                             iteration);
    }
    rComm.broadcast(&flVctTalpha);

    return flVctTalpha;
}

void DfDensityFittingX_Parallel::getTalpha_ROKS(
    TlDenseVector_Lapack* pT_alphaA, TlDenseVector_Lapack* pT_alphaB) {
    assert(pT_alphaA != NULL);
    assert(pT_alphaB != NULL);

    TlCommunicate& rComm = TlCommunicate::getInstance();

    if (this->m_bDiskUtilization == false) {
        if (this->isUpdateMethod_ == false) {
            TlDenseSymmetricMatrix_Lapack PA;
            TlDenseSymmetricMatrix_Lapack PB;
            if (rComm.isMaster() == true) {
                const TlDenseSymmetricMatrix_Lapack Pc =
                    this->getPInMatrix<TlDenseSymmetricMatrix_Lapack>(
                        RUN_ROKS_CLOSED, this->m_nIteration);
                const TlDenseSymmetricMatrix_Lapack Po =
                    this->getPInMatrix<TlDenseSymmetricMatrix_Lapack>(
                        RUN_ROKS_OPEN, this->m_nIteration);
                PB = 0.5 * Pc;
                PA = PB + Po;
            }
            rComm.broadcast(&PA);
            rComm.broadcast(&PB);

            *pT_alphaA = this->calcTAlpha_DIRECT(PA);
            *pT_alphaB = this->calcTAlpha_DIRECT(PB);
        } else {
            const TlDenseSymmetricMatrix_Lapack dPc =
                this->getDiffDensityMatrix(RUN_ROKS_CLOSED);
            const TlDenseSymmetricMatrix_Lapack dPo =
                this->getDiffDensityMatrix(RUN_ROKS_OPEN);
            const TlDenseSymmetricMatrix_Lapack dPB = 0.5 * dPc;
            const TlDenseSymmetricMatrix_Lapack dPA = dPB + dPo;

            *pT_alphaA = this->calcTAlpha_DIRECT(dPA);
            *pT_alphaB = this->calcTAlpha_DIRECT(dPB);
            if (this->m_nIteration > 1) {
                *pT_alphaA +=
                    this->getTalpha(RUN_UKS_ALPHA, this->m_nIteration - 1);
                *pT_alphaB +=
                    this->getTalpha(RUN_UKS_BETA, this->m_nIteration - 1);
            }
        }

        if (rComm.isMaster() == true) {
            pT_alphaA->save(
                this->getTalphaPath(RUN_UKS_ALPHA, this->m_nIteration));
            pT_alphaB->save(
                this->getTalphaPath(RUN_UKS_BETA, this->m_nIteration));
        }
        rComm.barrier();
    } else {
        abort();
    }
}

TlDenseSymmetricMatrix_Lapack DfDensityFittingX_Parallel::getDiffDensityMatrix(
    RUN_TYPE runType) {
    TlCommunicate& rComm = TlCommunicate::getInstance();

    TlDenseSymmetricMatrix_Lapack diffP;
    if (rComm.isMaster() == true) {
        diffP = DfObject::getDiffDensityMatrix<TlDenseSymmetricMatrix_Lapack>(
            runType, this->m_nIteration);
    }
    rComm.broadcast(&diffP);

    return diffP;
}

TlDenseSymmetricMatrix_Lapack DfDensityFittingX_Parallel::getP1pq(
    const int nIteration) {
    TlCommunicate& rComm = TlCommunicate::getInstance();

    TlDenseSymmetricMatrix_Lapack P;
    if (rComm.isMaster() == true) {
        P = DfDensityFittingTmpl<TlDenseSymmetricMatrix_Lapack,
                                 TlDenseVector_Lapack,
                                 DfEriX_Parallel>::getP1pq(nIteration);
    }
    rComm.broadcast(&P);

    return P;
}

TlDenseSymmetricMatrix_Lapack DfDensityFittingX_Parallel::getP2pq(
    const int nIteration) {
    TlCommunicate& rComm = TlCommunicate::getInstance();

    TlDenseSymmetricMatrix_Lapack P;
    if (rComm.isMaster() == true) {
        P = DfDensityFittingTmpl<TlDenseSymmetricMatrix_Lapack,
                                 TlDenseVector_Lapack,
                                 DfEriX_Parallel>::getP2pq(nIteration);
    }
    rComm.broadcast(&P);

    return P;
}

double DfDensityFittingX_Parallel::getLamda(const TlDenseVector_Lapack& SabinvN,
                                            const TlDenseVector_Lapack& t_alpha,
                                            const TlDenseVector_Lapack& N_alpha,
                                            const double dNumOfElec) {
    TlCommunicate& rComm = TlCommunicate::getInstance();

    double dAnswer;
    if (rComm.isMaster() == true) {
        dAnswer = DfDensityFittingTmpl<
            TlDenseSymmetricMatrix_Lapack, TlDenseVector_Lapack,
            DfEriX_Parallel>::getLamda(SabinvN, t_alpha, N_alpha, dNumOfElec);
    }
    rComm.broadcast(dAnswer);

    return dAnswer;
}

void DfDensityFittingX_Parallel::saveRho(const TlDenseVector_Lapack& rRho,
                                         const RUN_TYPE runType) {
    TlCommunicate& rComm = TlCommunicate::getInstance();

    if (rComm.isMaster() == true) {
        DfDensityFittingTmpl<TlDenseSymmetricMatrix_Lapack,
                             TlDenseVector_Lapack,
                             DfEriX_Parallel>::saveRho(rRho, runType);
    }
    rComm.barrier();
}
