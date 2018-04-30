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

#ifdef HAVE_CONFIG_H
#include "config.h"  // this file created by autotools
#endif               // HAVE_CONFIG_H

#include "DfXMatrix_Parallel.h"
#include "TlCommunicate.h"
#include "TlDistributeMatrix.h"
#include "TlDistributeSymmetricMatrix.h"

#define FAST_TRANCATE

DfXMatrix_Parallel::DfXMatrix_Parallel(TlSerializeData* pPdfParam)
    : DfXMatrix(pPdfParam) {}

DfXMatrix_Parallel::~DfXMatrix_Parallel() {}

void DfXMatrix_Parallel::buildX() {
#ifdef HAVE_SCALAPACK
  if (this->m_bUsingSCALAPACK == true) {
    this->buildX_ScaLAPACK();
  } else {
    this->buildX_LAPACK();
  }
#else
  { this->buildX_LAPACK(); }
#endif  // HAVE_SCALAPACK
}

void DfXMatrix_Parallel::buildX_LAPACK() {
  this->log_.info("build X using replica mode.");

  TlCommunicate& rComm = TlCommunicate::getInstance();
  if (rComm.isMaster() == true) {
    DfXMatrix::buildX();
  }

  index_type numOfMOs = 0;
  if (rComm.isMaster() == true) {
    numOfMOs = (*(this->pPdfParam_))["num_of_MOs"].getInt();
  }
  rComm.broadcast(numOfMOs);
  (*(this->pPdfParam_))["num_of_MOs"] = numOfMOs;
}

void DfXMatrix_Parallel::buildX_ScaLAPACK() {
  TlDistributeSymmetricMatrix S =
      this->getSpqMatrix<TlDistributeSymmetricMatrix>();
  TlDistributeMatrix X;
  TlDistributeMatrix Xinv;

  DfXMatrix::canonicalOrthogonalizeTmpl<TlDistributeSymmetricMatrix,
                                        TlDistributeMatrix>(
      S, &X, &Xinv, this->XEigvalFilePath_);

  (*(this->pPdfParam_))["num_of_MOs"] = X.getNumOfCols();

  DfObject::saveXMatrix(X);
  DfObject::saveXInvMatrix(Xinv);
}

void DfXMatrix_Parallel::canonicalOrthogonalize(
    const TlSymmetricMatrix& S, TlMatrix* pX, TlMatrix* pXinv,
    const std::string& eigvalFilePath) {
  TlCommunicate& rComm = TlCommunicate::getInstance();
  if (rComm.isMaster()) {
    DfXMatrix::canonicalOrthogonalize(S, pX, pXinv, eigvalFilePath);
  }
}

void DfXMatrix_Parallel::lowdinOrthogonalize(
    const TlSymmetricMatrix& S, TlMatrix* pX, TlMatrix* pXinv,
    const std::string& eigvalFilePath) {
  TlCommunicate& rComm = TlCommunicate::getInstance();
  if (rComm.isMaster()) {
    DfXMatrix::lowdinOrthogonalize(S, pX, pXinv, eigvalFilePath);
  }
}

void DfXMatrix_Parallel::canonicalOrthogonalize(
    const TlDistributeSymmetricMatrix& S, TlDistributeMatrix* pX,
    TlDistributeMatrix* pXinv, const std::string& eigvalFilePath) {
  DfXMatrix::canonicalOrthogonalizeTmpl(S, pX, pXinv, eigvalFilePath);
}

void DfXMatrix_Parallel::lowdinOrthogonalize(
    const TlDistributeSymmetricMatrix& S, TlDistributeMatrix* pX,
    TlDistributeMatrix* pXinv, const std::string& eigvalFilePath) {
  DfXMatrix::lowdinOrthogonalizeTmpl(S, pX, pXinv, eigvalFilePath);
}
