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
#include "tl_dense_symmetric_matrix_blas_old.h"
#include "tl_dense_vector_blacs.h"

// LAPACK ==============================================================
DfDensityFittingX_Parallel::DfDensityFittingX_Parallel(
    TlSerializeData* pPdfParam)
    : DfDensityFittingTmpl<TlDenseSymmetricMatrix_BLAS_Old, TlVector_BLAS,
                           DfEriX_Parallel>(pPdfParam) {}

DfDensityFittingX_Parallel::~DfDensityFittingX_Parallel() {}

void DfDensityFittingX_Parallel::exec() { this->calc(); }

TlVector_BLAS DfDensityFittingX_Parallel::getNalpha() {
  TlCommunicate& rComm = TlCommunicate::getInstance();

  TlVector_BLAS Nalpha;
  if (rComm.isMaster() == true) {
    Nalpha = DfDensityFittingTmpl<TlDenseSymmetricMatrix_BLAS_Old, TlVector_BLAS,
                                  DfEriX_Parallel>::getNalpha();
  }
  rComm.broadcast(Nalpha);

  return Nalpha;
}

TlDenseSymmetricMatrix_BLAS_Old DfDensityFittingX_Parallel::getSabinv() {
  TlCommunicate& rComm = TlCommunicate::getInstance();

  TlDenseSymmetricMatrix_BLAS_Old Sabinv;
  if (rComm.isMaster() == true) {
    Sabinv = DfDensityFittingTmpl<TlDenseSymmetricMatrix_BLAS_Old, TlVector_BLAS,
                                  DfEriX_Parallel>::getSabinv();
  }
  rComm.broadcast(Sabinv);

  return Sabinv;
}

TlVector_BLAS DfDensityFittingX_Parallel::calcTAlpha_DIRECT(
    const TlDenseSymmetricMatrix_BLAS_Old& P) {
  return DfDensityFittingTmpl<TlDenseSymmetricMatrix_BLAS_Old, TlVector_BLAS,
                              DfEriX_Parallel>::calcTAlpha_DIRECT(P);
}

TlVector_BLAS DfDensityFittingX_Parallel::getTalpha(const RUN_TYPE runType,
                                                    const int iteration) {
  TlCommunicate& rComm = TlCommunicate::getInstance();

  TlVector_BLAS flVctTalpha;
  if (rComm.isMaster() == true) {
    flVctTalpha =
        DfDensityFittingTmpl<TlDenseSymmetricMatrix_BLAS_Old, TlVector_BLAS,
                             DfEriX_Parallel>::getTalpha(runType, iteration);
  }
  rComm.broadcast(flVctTalpha);

  return flVctTalpha;
}

void DfDensityFittingX_Parallel::getTalpha_ROKS(TlVector_BLAS* pT_alphaA,
                                                TlVector_BLAS* pT_alphaB) {
  TlCommunicate& rComm = TlCommunicate::getInstance();

  if (rComm.isMaster() == true) {
    DfDensityFittingTmpl<TlDenseSymmetricMatrix_BLAS_Old, TlVector_BLAS,
                         DfEriX_Parallel>::getTalpha_ROKS(pT_alphaA, pT_alphaB);
  }

  rComm.broadcast(*pT_alphaA);
  rComm.broadcast(*pT_alphaB);
}

TlDenseSymmetricMatrix_BLAS_Old DfDensityFittingX_Parallel::getDiffDensityMatrix(
    RUN_TYPE runType) {
  TlCommunicate& rComm = TlCommunicate::getInstance();

  TlDenseSymmetricMatrix_BLAS_Old diffP;
  if (rComm.isMaster() == true) {
    diffP = DfObject::getDiffDensityMatrix<TlDenseSymmetricMatrix_BLAS_Old>(
        runType, this->m_nIteration);
  }
  rComm.broadcast(diffP);

  return diffP;
}

TlDenseSymmetricMatrix_BLAS_Old DfDensityFittingX_Parallel::getP1pq(
    const int nIteration) {
  TlCommunicate& rComm = TlCommunicate::getInstance();

  TlDenseSymmetricMatrix_BLAS_Old P;
  if (rComm.isMaster() == true) {
    P = DfDensityFittingTmpl<TlDenseSymmetricMatrix_BLAS_Old, TlVector_BLAS,
                             DfEriX_Parallel>::getP1pq(nIteration);
  }
  rComm.broadcast(P);

  return P;
}

TlDenseSymmetricMatrix_BLAS_Old DfDensityFittingX_Parallel::getP2pq(
    const int nIteration) {
  TlCommunicate& rComm = TlCommunicate::getInstance();

  TlDenseSymmetricMatrix_BLAS_Old P;
  if (rComm.isMaster() == true) {
    P = DfDensityFittingTmpl<TlDenseSymmetricMatrix_BLAS_Old, TlVector_BLAS,
                             DfEriX_Parallel>::getP2pq(nIteration);
  }
  rComm.broadcast(P);

  return P;
}

double DfDensityFittingX_Parallel::getLamda(const TlVector_BLAS& SabinvN,
                                            const TlVector_BLAS& t_alpha,
                                            const TlVector_BLAS& N_alpha,
                                            const double dNumOfElec) {
  TlCommunicate& rComm = TlCommunicate::getInstance();

  double dAnswer;
  if (rComm.isMaster() == true) {
    dAnswer =
        DfDensityFittingTmpl<TlDenseSymmetricMatrix_BLAS_Old, TlVector_BLAS,
                             DfEriX_Parallel>::getLamda(SabinvN, t_alpha,
                                                        N_alpha, dNumOfElec);
  }
  rComm.broadcast(dAnswer);

  return dAnswer;
}

void DfDensityFittingX_Parallel::saveRho(const TlVector_BLAS& rRho,
                                         const RUN_TYPE runType) {
  TlCommunicate& rComm = TlCommunicate::getInstance();

  if (rComm.isMaster() == true) {
    DfDensityFittingTmpl<TlDenseSymmetricMatrix_BLAS_Old, TlVector_BLAS,
                         DfEriX_Parallel>::saveRho(rRho, runType);
  }
  rComm.barrier();
}
