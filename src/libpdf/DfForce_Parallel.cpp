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

#include "DfForce_Parallel.h"
#include "DfEriX_Parallel.h"
#include "DfXCFunctional_Parallel.h"
#include "TlCommunicate.h"
#include "TlSymmetricMatrix.h"
#include "TlVector.h"
#include "config.h"

DfForce_Parallel::DfForce_Parallel(TlSerializeData* pPdfParam)
    : DfForce(pPdfParam) {
  TlCommunicate& rComm = TlCommunicate::getInstance();
  rComm.broadcast(this->pdfParamForForce_);
}

DfForce_Parallel::~DfForce_Parallel() {}

void DfForce_Parallel::logger(const std::string& str) const {
  TlCommunicate& rComm = TlCommunicate::getInstance();
  if (rComm.isMaster() == true) {
    DfForce::logger(str);
  }
}

void DfForce_Parallel::output() {
  TlCommunicate& rComm = TlCommunicate::getInstance();
  if (rComm.isMaster()) {
    DfForce::output();
  }
}

void DfForce_Parallel::calcForceFromNuclei() {
  TlCommunicate& rComm = TlCommunicate::getInstance();
  if (rComm.isMaster() == true) {
    // perform in master-node only
    DfForce::calcForceFromNuclei();
  }
}

void DfForce_Parallel::calcForceFromWS(const RUN_TYPE runType) {
  TlCommunicate& rComm = TlCommunicate::getInstance();
  if (rComm.isMaster() == true) {
    // perform in master-node only
    DfForce::calcForceFromWS(runType);
  }
}

void DfForce_Parallel::calcForceFromHpq(const TlSymmetricMatrix& P) {
  TlCommunicate& rComm = TlCommunicate::getInstance();
  if (rComm.isMaster() == true) {
    // perform in master-node only
    DfForce::calcForceFromHpq(P);
  }
}

void DfForce_Parallel::calcForceFromCoulomb_exact(const RUN_TYPE runType) {
  // #ifdef HAVE_SCALAPACK
  // #else
  // #endif // HAVE_SCALAPACK
  this->calcForceFromCoulomb_exact_replicated(runType);
}

void DfForce_Parallel::calcForceFromCoulomb_exact_replicated(
    const RUN_TYPE runType) {
  this->loggerTime("calc force from J(parallel, replicated)");

  TlCommunicate& rComm = TlCommunicate::getInstance();
  const int iteration = this->m_nIteration;
  const int numOfAtoms = this->m_nNumOfAtoms;

  TlSymmetricMatrix P;
  if (rComm.isMaster() == true) {
    P = this->getPpqMatrix<TlSymmetricMatrix>(runType, iteration);
  }
  rComm.broadcast(P);

  DfEriX_Parallel dfEri(&(this->pdfParamForForce_));

  // ((pq)'|(rs))
  TlMatrix F_J(numOfAtoms, 3);
  dfEri.getForceJ(P, &F_J);

  // F_J *= 0.5;
  if (this->isDebugOutMatrix_ == true) {
    if (rComm.isMaster() == true) {
      F_J.save("F_J.mtx");
    }
  }

  this->force_ += F_J;
}

void DfForce_Parallel::calcForceFromCoulomb_RIJ(const RUN_TYPE runType) {
  this->calcForceFromCoulomb_RIJ_DC(runType);
}

void DfForce_Parallel::calcForceFromCoulomb_RIJ_DC(const RUN_TYPE runType) {
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

  DfEriX_Parallel dfEri(&(this->pdfParamForForce_));

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

void DfForce_Parallel::calcForceFromK(const RUN_TYPE runType) {
  // #ifdef HAVE_SCALAPACK
  // #else
  // #endif // HAVE_SCALAPACK
  this->calcForceFromK_replicated(runType);
}

void DfForce_Parallel::calcForceFromK_replicated(const RUN_TYPE runType) {
  this->loggerTime("calc force from K (parallel; replicated)");

  TlCommunicate& rComm = TlCommunicate::getInstance();
  const DfXCFunctional_Parallel dfXCFunctional(&(this->pdfParamForForce_));

  if (dfXCFunctional.isHybridFunctional() == true) {
    const int iteration = this->m_nIteration;
    const int numOfAtoms = this->m_nNumOfAtoms;

    TlSymmetricMatrix P;
    if (rComm.isMaster() == true) {
      P = this->getPpqMatrix<TlSymmetricMatrix>(runType, iteration);
    }
    rComm.broadcast(P);

    DfEriX_Parallel dfEri(&(this->pdfParamForForce_));

    TlMatrix F_K(numOfAtoms, 3);
    // for RKS
    dfEri.getForceK(P, &F_K);
    if (runType == RUN_RKS) {
      F_K *= 0.5;
    }

    F_K *= -1.0;
    F_K *= dfXCFunctional.getFockExchangeCoefficient();  // for B3LYP

    if (this->isDebugOutMatrix_ == true) {
      if (rComm.isMaster() == true) {
        F_K.save("F_K.mtx");
      }
    }
    this->force_ += F_K;
  }
}

DfCalcGridX* DfForce_Parallel::getCalcGridObj() {
  this->logger("create calc-grid engine(parallel)");
  DfCalcGridX* pDfCalcGrid =
      new DfCalcGridX_Parallel(&(this->pdfParamForForce_));

  return pDfCalcGrid;
}
