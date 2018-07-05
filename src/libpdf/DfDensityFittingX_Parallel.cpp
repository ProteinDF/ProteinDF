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
                             DfEriX_Parallel>::getTalpha(runType, iteration);
  }
  rComm.broadcast(&flVctTalpha);

  return flVctTalpha;
}

void DfDensityFittingX_Parallel::getTalpha_ROKS(
    TlDenseVector_Lapack* pT_alphaA, TlDenseVector_Lapack* pT_alphaB) {
  TlCommunicate& rComm = TlCommunicate::getInstance();

  if (rComm.isMaster() == true) {
    DfDensityFittingTmpl<TlDenseSymmetricMatrix_Lapack, TlDenseVector_Lapack,
                         DfEriX_Parallel>::getTalpha_ROKS(pT_alphaA, pT_alphaB);
  }

  rComm.broadcast(pT_alphaA);
  rComm.broadcast(pT_alphaB);
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
    dAnswer =
        DfDensityFittingTmpl<TlDenseSymmetricMatrix_Lapack,
                             TlDenseVector_Lapack,
                             DfEriX_Parallel>::getLamda(SabinvN, t_alpha,
                                                        N_alpha, dNumOfElec);
  }
  rComm.broadcast(dAnswer);

  return dAnswer;
}

void DfDensityFittingX_Parallel::saveRho(const TlDenseVector_Lapack& rRho,
                                         const RUN_TYPE runType) {
  TlCommunicate& rComm = TlCommunicate::getInstance();

  if (rComm.isMaster() == true) {
    DfDensityFittingTmpl<TlDenseSymmetricMatrix_Lapack, TlDenseVector_Lapack,
                         DfEriX_Parallel>::saveRho(rRho, runType);
  }
  rComm.barrier();
}
