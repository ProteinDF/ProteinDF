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

#include "DfKMatrix_Parallel.h"
#include <cassert>
#include "CnError.h"
#include "DfCD_Parallel.h"
#include "DfEriX_Parallel.h"
#include "TlCommunicate.h"
#include "tl_dense_symmetric_matrix_scalapack.h"

DfKMatrix_Parallel::DfKMatrix_Parallel(TlSerializeData* pPdfParam)
    : DfKMatrix(pPdfParam) {}

DfKMatrix_Parallel::~DfKMatrix_Parallel() {}

void DfKMatrix_Parallel::getK_CD() {
#ifdef HAVE_SCALAPACK
  if (this->m_bUsingSCALAPACK == true) {
    TlDenseSymmetricMatrix_Scalapack K(this->m_nNumOfAOs);
    this->getK_CD_distributed(RUN_RKS, &K);
    DfObject::saveHFxMatrix(RUN_RKS, this->m_nIteration, K);
  } else {
    TlDenseSymmetricMatrix_Lapack K(this->m_nNumOfAOs);
    this->getK_CD_local(RUN_RKS, &K);
    this->saveKMatrix(RUN_RKS, K);
  }
#else
  {
    TlDenseSymmetricMatrix_Lapack K(this->m_nNumOfAOs);
    this->getK_CD_local(RUN_RKS, &K);
    this->saveKMatrix(RUN_RKS, K);
  }
#endif  // HAVE_SCALAPACK
}

void DfKMatrix_Parallel::getK_conventional() {
#ifdef HAVE_SCALAPACK
  if (this->m_bUsingSCALAPACK == true) {
    TlDenseSymmetricMatrix_Scalapack K(this->m_nNumOfAOs);
    this->getK_conventional_distributed(RUN_RKS, &K);
    DfObject::saveHFxMatrix(RUN_RKS, this->m_nIteration, K);
  } else {
    TlDenseSymmetricMatrix_Lapack K(this->m_nNumOfAOs);
    this->getK_conventional_local(RUN_RKS, &K);
    this->saveKMatrix(RUN_RKS, K);
  }
#else
  {
    TlDenseSymmetricMatrix_Lapack K(this->m_nNumOfAOs);
    this->getK_conventional_local(RUN_RKS, &K);
    this->saveKMatrix(RUN_RKS, K);
  }
#endif  // HAVE_SCALAPACK
}

void DfKMatrix_Parallel::saveKMatrix(const RUN_TYPE runType,
                                     const TlDenseSymmetricMatrix_Lapack& K) {
  TlCommunicate& rComm = TlCommunicate::getInstance();
  if (rComm.isMaster() == true) {
    DfKMatrix::saveKMatrix(runType, K);
  }
}

void DfKMatrix_Parallel::getK_CD_local(const RUN_TYPE runType,
                                       TlDenseSymmetricMatrix_Lapack* pK) {
  DfCD_Parallel dfCD(this->pPdfParam_);
  dfCD.getK(runType, pK);
}

void DfKMatrix_Parallel::getK_CD_distributed(
    const RUN_TYPE runType, TlDenseSymmetricMatrix_Scalapack* pK) {
  DfCD_Parallel dfCD(this->pPdfParam_);
  dfCD.getK_D(runType, pK);
}

void DfKMatrix_Parallel::getK_conventional_local(
    const RUN_TYPE runType, TlDenseSymmetricMatrix_Lapack* pK) {
  TlCommunicate& rComm = TlCommunicate::getInstance();
  assert(rComm.checkNonBlockingCommunications());
  rComm.barrier();
  this->log_.info("build K on replica parallel method");

  TlDenseSymmetricMatrix_Lapack P;
  if (rComm.isMaster() == true) {
    if (this->isUpdateMethod_ == true) {
      P = this->getDiffDensityMatrix<TlDenseSymmetricMatrix_Lapack>(runType);
    } else {
      P = this->getDensityMatrix<TlDenseSymmetricMatrix_Lapack>(runType);
    }
  }
  rComm.broadcast(&P);
  assert(rComm.checkNonBlockingCommunications());
  this->log_.info("broadcast density matrix");
  this->log_.info(TlUtils::format("density matrix size: %dx%d",
                                  P.getNumOfRows(), P.getNumOfCols()));
  this->log_.info(TlUtils::format("num of AOs: %d", this->m_nNumOfAOs));
  assert(P.getNumOfRows() == this->m_nNumOfAOs);

  this->log_.info("ERI operation: start");
  DfEriX_Parallel dfEri(this->pPdfParam_);
  dfEri.getK(P, pK);
  this->log_.info("ERI operation: end");

  if (this->isUpdateMethod_ == true) {
    if (this->m_nIteration > 1) {
      this->log_.info("update K matrix");
      const TlDenseSymmetricMatrix_Lapack prevK =
          this->getKMatrix(runType, this->m_nIteration - 1);
      *pK += prevK;
    }
  }
}

void DfKMatrix_Parallel::getK_conventional_distributed(
    const RUN_TYPE runType, TlDenseSymmetricMatrix_Scalapack* pK) {
  this->log_.info("build K on distributed parallel method");

  TlDenseSymmetricMatrix_Scalapack P;
  if (this->isUpdateMethod_ == true) {
    P = DfObject::getDiffDensityMatrix<TlDenseSymmetricMatrix_Scalapack>(
        runType, this->m_nIteration);
    if (runType == RUN_RKS) {
      P *= 0.5;
    }
  } else {
    P = DfObject::getPpqMatrix<TlDenseSymmetricMatrix_Scalapack>(
        runType, this->m_nIteration - 1);
    if (runType == RUN_RKS) {
      P *= 0.5;
    }
  }

  // TlCommunicate& rComm = TlCommunicate::getInstance();
  // assert(rComm.checkNonBlockingCommunications());
  assert(P.getNumOfRows() == this->m_nNumOfAOs);
  this->log_.info("density matrix loaded.");

  DfEriX_Parallel dfEri(this->pPdfParam_);
  dfEri.getK_D(P, pK);

  if (this->isUpdateMethod_ == true) {
    if (this->m_nIteration > 1) {
      this->log_.info("update K matrix: start");
      const TlDenseSymmetricMatrix_Scalapack prevK =
          DfObject::getHFxMatrix<TlDenseSymmetricMatrix_Scalapack>(
              runType, this->m_nIteration - 1);
      *pK += prevK;
      this->log_.info("update K matrix: end");
    }
  }
}

TlDenseSymmetricMatrix_Lapack DfKMatrix_Parallel::getKMatrix(
    const RUN_TYPE runType, const int iteration) {
  TlCommunicate& rComm = TlCommunicate::getInstance();

  TlDenseSymmetricMatrix_Lapack K;
  if (rComm.isMaster() == true) {
    K = DfKMatrix::getKMatrix(runType, iteration);
  }
  rComm.broadcast(&K);

  return K;
}
