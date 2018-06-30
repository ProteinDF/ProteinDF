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

#include "CnError.h"
#include "DfDiagonal.h"
#include "TlUtils.h"
#include "tl_dense_symmetric_matrix_blas_old.h"

DfDiagonal::DfDiagonal(TlSerializeData* pPdfParam) : DfObject(pPdfParam) {}

DfDiagonal::~DfDiagonal() {}

void DfDiagonal::DfDiagMain() {
  // output informations
  switch (this->m_nMethodType) {
    case METHOD_RKS:
      this->main<TlDenseGeneralMatrix_BLAS_old, TlDenseSymmetricMatrix_BLAS_Old>(
          RUN_RKS);
      break;

    case METHOD_UKS:
      this->main<TlDenseGeneralMatrix_BLAS_old, TlDenseSymmetricMatrix_BLAS_Old>(
          RUN_UKS_ALPHA);
      this->main<TlDenseGeneralMatrix_BLAS_old, TlDenseSymmetricMatrix_BLAS_Old>(
          RUN_UKS_BETA);
      break;

    case METHOD_ROKS:
      this->main<TlDenseGeneralMatrix_BLAS_old, TlDenseSymmetricMatrix_BLAS_Old>(
          RUN_ROKS);
      break;

    default:
      CnErr.abort("DfDiagonal", "", "DfDiagMain",
                  "the value of scftype is illegal");
      break;
  }
}

// for extended QCLO method
void DfDiagonal::DfDiagQclo(DfObject::RUN_TYPE runType,
                            const std::string& fragname, int norbcut) {
  this->m_nNumOfMOs = norbcut;
  this->main<TlDenseGeneralMatrix_BLAS_old, TlDenseSymmetricMatrix_BLAS_Old>(
      runType, fragname, true);
}
