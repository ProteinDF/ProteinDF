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

#include "DfJMatrix_Parallel.h"
#include <cassert>
#include "CnError.h"
#include "DfCD_Parallel.h"
#include "DfEriX_Parallel.h"
#include "TlCommunicate.h"
#include "tl_dense_symmetric_matrix_lapack.h"
#include "tl_dense_symmetric_matrix_scalapack.h"

DfJMatrix_Parallel::DfJMatrix_Parallel(TlSerializeData* pPdfParam)
    : DfJMatrix(pPdfParam) {}

DfJMatrix_Parallel::~DfJMatrix_Parallel() {}

void DfJMatrix_Parallel::getJ_RI() {
#ifdef HAVE_SCALAPACK
  if (this->m_bUsingSCALAPACK == true) {
    TlDenseSymmetricMatrix_Scalapack J(this->m_nNumOfAOs);
    this->getJ_RI_distributed(&J);
    DfObject::saveJMatrix(this->m_nIteration, J);
  } else {
    TlDenseSymmetricMatrix_Lapack J(this->m_nNumOfAOs);
    this->getJ_RI_local(&J);
    this->saveJMatrix(J);
  }
#else
  {
    TlDenseSymmetricMatrix_Lapack J(this->m_nNumOfAOs);
    this->getJ_RI_local(&J);
    this->saveJMatrix(J);
  }
#endif  // HAVE_SCALAPACK
}

void DfJMatrix_Parallel::getJ_CD() {
#ifdef HAVE_SCALAPACK
  if (this->m_bUsingSCALAPACK == true) {
    TlDenseSymmetricMatrix_Scalapack J(this->m_nNumOfAOs);
    this->getJ_CD_distributed(&J);
    DfObject::saveJMatrix(this->m_nIteration, J);
  } else {
    TlDenseSymmetricMatrix_Lapack J(this->m_nNumOfAOs);
    this->getJ_CD_local(&J);
    this->saveJMatrix(J);
  }
#else
  {
    TlDenseSymmetricMatrix_Lapack J(this->m_nNumOfAOs);
    this->getJ_CD_local(&J);
    this->saveJMatrix(J);
  }
#endif  // HAVE_SCALAPACK
}

void DfJMatrix_Parallel::getJ_conventional() {
#ifdef HAVE_SCALAPACK
  if (this->m_bUsingSCALAPACK == true) {
    TlDenseSymmetricMatrix_Scalapack J(this->m_nNumOfAOs);
    this->getJ_conventional_distributed(&J);
    DfObject::saveJMatrix(this->m_nIteration, J);
  } else {
    TlDenseSymmetricMatrix_Lapack J(this->m_nNumOfAOs);
    this->getJ_conventional_local(&J);
    this->saveJMatrix(J);
  }
#else
  {
    TlDenseSymmetricMatrix_Lapack J(this->m_nNumOfAOs);
    this->getJ_conventional_local(&J);
    this->saveJMatrix(J);
  }
#endif  // HAVE_SCALAPACK
}

void DfJMatrix_Parallel::saveJMatrix(const TlDenseSymmetricMatrix_Lapack& J) {
  TlCommunicate& rComm = TlCommunicate::getInstance();
  if (rComm.isMaster() == true) {
    DfObject::saveJMatrix(this->m_nIteration, J);
  }
}

void DfJMatrix_Parallel::getJ_RI_local(TlDenseSymmetricMatrix_Lapack* pJ) {
  TlDenseVector_Lapack rho;
  switch (this->m_nMethodType) {
    case METHOD_RKS:
      rho = this->getRho(RUN_RKS, this->m_nIteration);
      break;
    case METHOD_UKS:
      rho = this->getRho(RUN_UKS_ALPHA, this->m_nIteration);
      rho += this->getRho(RUN_UKS_BETA, this->m_nIteration);
      break;
    default:
      std::cerr << "unsopported. sorry." << std::endl;
      CnErr.abort();
      break;
  }

  if (this->isUpdateMethod_ == true) {
    if (this->m_nIteration > 1) {
      TlDenseVector_Lapack prevRho;
      switch (this->m_nMethodType) {
        case METHOD_RKS:
          prevRho = this->getRho(RUN_RKS, this->m_nIteration - 1);
          break;
        case METHOD_UKS:
          prevRho = this->getRho(RUN_UKS_ALPHA, this->m_nIteration - 1);
          prevRho += this->getRho(RUN_UKS_BETA, this->m_nIteration - 1);
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
  dfEri.getJ(rho, pJ);

  if (this->isUpdateMethod_ == true) {
    if (this->m_nIteration > 1) {
      TlDenseSymmetricMatrix_Lapack prevJ =
          this->getJMatrix(this->m_nIteration - 1);
      *pJ += prevJ;
    }
  }
}

TlDenseVector_Lapack DfJMatrix_Parallel::getRho(const RUN_TYPE runType,
                                                const int iteration) {
  TlCommunicate& rComm = TlCommunicate::getInstance();
  TlDenseVector_Lapack rho;
  if (rComm.isMaster() == true) {
    rho = DfObject::getRho<TlDenseVector_Lapack>(runType, iteration);
  }
  rComm.broadcast(&rho);

  return rho;
}

void DfJMatrix_Parallel::getJ_RI_distributed(
    TlDenseSymmetricMatrix_Scalapack* pJ) {
  TlDenseVector_Lapack rho;
  switch (this->m_nMethodType) {
    case METHOD_RKS:
      rho = this->getRho(RUN_RKS, this->m_nIteration);
      break;
    case METHOD_UKS:
      rho = this->getRho(RUN_UKS_ALPHA, this->m_nIteration);
      rho += this->getRho(RUN_UKS_BETA, this->m_nIteration);
      break;
    default:
      std::cerr << "unsopported. sorry." << std::endl;
      CnErr.abort();
      break;
  }

  if (this->isUpdateMethod_ == true) {
    if (this->m_nIteration > 1) {
      TlDenseVector_Lapack prevRho;
      switch (this->m_nMethodType) {
        case METHOD_RKS:
          prevRho = this->getRho(RUN_RKS, this->m_nIteration - 1);
          break;
        case METHOD_UKS:
          prevRho = this->getRho(RUN_UKS_ALPHA, this->m_nIteration - 1);
          prevRho += this->getRho(RUN_UKS_BETA, this->m_nIteration - 1);
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
  dfEri.getJ(rho, pJ);

  if (this->isUpdateMethod_ == true) {
    if (this->m_nIteration > 1) {
      TlDenseSymmetricMatrix_Scalapack prevJ =
          DfObject::getJMatrix<TlDenseSymmetricMatrix_Scalapack>(
              this->m_nIteration - 1);
      *pJ += prevJ;
    }
  }
}

void DfJMatrix_Parallel::getJ_CD_local(TlDenseSymmetricMatrix_Lapack* pJ) {
  DfCD_Parallel dfCD(this->pPdfParam_);
  dfCD.getJ(pJ);
}

void DfJMatrix_Parallel::getJ_CD_distributed(
    TlDenseSymmetricMatrix_Scalapack* pJ) {
  DfCD_Parallel dfCD(this->pPdfParam_);
  dfCD.getJ_D(pJ);
}

void DfJMatrix_Parallel::getJ_conventional_local(
    TlDenseSymmetricMatrix_Lapack* pJ) {
  TlCommunicate& rComm = TlCommunicate::getInstance();

  TlDenseSymmetricMatrix_Lapack P;
  if (rComm.isMaster() == true) {
    if (this->isUpdateMethod_ == true) {
      P = this->getDiffDensityMatrix<TlDenseSymmetricMatrix_Lapack>();
    } else {
      P = this->getDensityMatrix<TlDenseSymmetricMatrix_Lapack>();
    }
  }
  rComm.broadcast(&P);
  assert(P.getNumOfRows() == this->m_nNumOfAOs);

  DfEriX_Parallel dfEri(this->pPdfParam_);
  dfEri.getJpq(P, pJ);

  if (this->isUpdateMethod_ == true) {
    if (this->m_nIteration > 1) {
      const TlDenseSymmetricMatrix_Lapack prevJ =
          this->getJMatrix(this->m_nIteration - 1);
      *pJ += prevJ;
    }
  }
}

void DfJMatrix_Parallel::getJ_conventional_distributed(
    TlDenseSymmetricMatrix_Scalapack* pJ) {
  // TlCommunicate& rComm = TlCommunicate::getInstance();
  TlDenseSymmetricMatrix_Scalapack P;
  if (this->isUpdateMethod_ == true) {
    P = this->getDiffDensityMatrix<TlDenseSymmetricMatrix_Scalapack>();
  } else {
    P = DfObject::getPpqMatrix<TlDenseSymmetricMatrix_Scalapack>(
        RUN_RKS, this->m_nIteration - 1);
  }
  assert(P.getNumOfRows() == this->m_nNumOfAOs);

  DfEriX_Parallel dfEri(this->pPdfParam_);
  dfEri.getJpq_D(P, pJ);

  if (this->isUpdateMethod_ == true) {
    if (this->m_nIteration > 1) {
      this->log_.info("update J matrix: start");
      const TlDenseSymmetricMatrix_Scalapack prevJ =
          DfObject::getJMatrix<TlDenseSymmetricMatrix_Scalapack>(
              this->m_nIteration - 1);
      *pJ += prevJ;
      this->log_.info("update J matrix: end");
    }
  }
}

TlDenseSymmetricMatrix_Lapack DfJMatrix_Parallel::getJMatrix(
    const int iteration) {
  TlCommunicate& rComm = TlCommunicate::getInstance();

  TlDenseSymmetricMatrix_Lapack J;
  if (rComm.isMaster() == true) {
    J = DfJMatrix::getJMatrix(iteration);
  }
  rComm.broadcast(&J);

  return J;
}
