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
#include <cmath>
#include <fstream>
#include <string>

#include "DfXMatrix.h"
#include "TlMatrix.h"
#include "TlSymmetricMatrix.h"
#include "TlVector.h"

DfXMatrix::DfXMatrix(TlSerializeData* pPdfParam) : DfObject(pPdfParam) {
  assert(pPdfParam != NULL);
  const TlSerializeData& pdfParam = *pPdfParam;
  const double threshold_trancation =
      pdfParam["orbital_independence_threshold"].getDouble();

  this->threshold_trancation_canonical_ = threshold_trancation;
  if (pdfParam["orbital_independence_threshold/canonical"].getStr().empty() !=
      true) {
    this->threshold_trancation_canonical_ =
        pdfParam["orbital_independence_threshold/canonical"].getDouble();
  }

  this->threshold_trancation_lowdin_ = threshold_trancation;
  if (pdfParam["orbital_independence_threshold/lowdin"].getStr().empty() !=
      true) {
    this->threshold_trancation_lowdin_ =
        pdfParam["orbital_independence_threshold/lowdin"].getDouble();
  }

  this->XEigvalFilePath_ = "";
  if (pdfParam["XMatrix/save_eigval"].getBoolean()) {
    this->XEigvalFilePath_ = DfObject::getXEigvalVtrPath();
  }

  this->debug_save_mat_ = pdfParam["debug/DfXMatrix/save_mat"].getBoolean();
  this->debug_check_X_ = pdfParam["debug/DfXMatrix/check_X"].getBoolean();
}

DfXMatrix::~DfXMatrix() {}

void DfXMatrix::buildX() {
  TlSymmetricMatrix S = this->getSpqMatrix<TlSymmetricMatrix>();
  TlMatrix X;
  TlMatrix Xinv;

  this->canonicalOrthogonalize(S, &X, &Xinv, this->XEigvalFilePath_);

  DfObject::saveXMatrix(X);
  DfObject::saveXInvMatrix(Xinv);
  (*(this->pPdfParam_))["num_of_MOs"] = X.getNumOfCols();
}

void DfXMatrix::canonicalOrthogonalize(const TlSymmetricMatrix& S, TlMatrix* pX,
                                       TlMatrix* pXinv,
                                       const std::string& eigvalFilePath) {
  this->canonicalOrthogonalizeTmpl<TlSymmetricMatrix, TlMatrix>(S, pX, pXinv,
                                                                eigvalFilePath);
}

void DfXMatrix::lowdinOrthogonalize(const TlSymmetricMatrix& S, TlMatrix* pX,
                                    TlMatrix* pXinv,
                                    const std::string& eigvalFilePath) {
  this->lowdinOrthogonalizeTmpl<TlSymmetricMatrix, TlMatrix>(S, pX, pXinv,
                                                             eigvalFilePath);
}
