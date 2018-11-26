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
#include "tl_dense_general_matrix_lapack.h"
#include "tl_dense_symmetric_matrix_lapack.h"

DfXMatrix::DfXMatrix(TlSerializeData* pPdfParam) : DfObject(pPdfParam) {
  assert(pPdfParam != NULL);
  const TlSerializeData& pdfParam = *pPdfParam;
  const double threshold_trancation =
      pdfParam["orbital-independence-threshold"].getDouble();

  this->threshold_trancation_canonical_ = threshold_trancation;
  if (pdfParam["orbital-independence-threshold/canonical"].getStr().empty() !=
      true) {
    this->threshold_trancation_canonical_ =
        pdfParam["orbital-independence-threshold/canonical"].getDouble();
  }

  this->threshold_trancation_lowdin_ = threshold_trancation;
  if (pdfParam["orbital-independence-threshold/lowdin"].getStr().empty() !=
      true) {
    this->threshold_trancation_lowdin_ =
        pdfParam["orbital-independence-threshold/lowdin"].getDouble();
  }

  this->debugSaveEigval_ = pdfParam["debug/DfXMatrix/save-eigval"].getBoolean();
  this->debugSaveMatrix_ = pdfParam["debug/DfXMatrix/save-mat"].getBoolean();
  this->debugCheckX_ = pdfParam["debug/DfXMatrix/check-X"].getBoolean();
}

DfXMatrix::~DfXMatrix() {}

void DfXMatrix::buildX() {
  TlDenseSymmetricMatrix_Lapack S =
      this->getSpqMatrix<TlDenseSymmetricMatrix_Lapack>();
  TlDenseGeneralMatrix_Lapack X;
  TlDenseGeneralMatrix_Lapack Xinv;

  std::string eigvalFilePath = "";
  if (this->debugSaveEigval_) {
    eigvalFilePath = DfObject::getXEigvalVtrPath();
  }
  this->canonicalOrthogonalize(S, &X, &Xinv, eigvalFilePath);

  DfObject::saveXMatrix(X);
  DfObject::saveXInvMatrix(Xinv);
  (*(this->pPdfParam_))["num_of_MOs"] = X.getNumOfCols();
}

void DfXMatrix::canonicalOrthogonalize(const TlDenseSymmetricMatrix_Lapack& S,
                                       TlDenseGeneralMatrix_Lapack* pX,
                                       TlDenseGeneralMatrix_Lapack* pXinv,
                                       const std::string& eigvalFilePath) {
  this->canonicalOrthogonalizeTmpl<TlDenseSymmetricMatrix_Lapack,
                                   TlDenseGeneralMatrix_Lapack>(S, pX, pXinv,
                                                                eigvalFilePath);
}

void DfXMatrix::lowdinOrthogonalize(const TlDenseSymmetricMatrix_Lapack& S,
                                    TlDenseGeneralMatrix_Lapack* pX,
                                    TlDenseGeneralMatrix_Lapack* pXinv,
                                    const std::string& eigvalFilePath) {
  this->lowdinOrthogonalizeTmpl<TlDenseSymmetricMatrix_Lapack,
                                TlDenseGeneralMatrix_Lapack>(S, pX, pXinv,
                                                             eigvalFilePath);
}
