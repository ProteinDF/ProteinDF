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
#include "common.h"

DfDiagonal::DfDiagonal(TlSerializeData* pPdfParam) : DfObject(pPdfParam) {
    this->updateLinearAlgebraPackageParam(
      (*(this->pPdfParam_))["linear_algebra_package/diagonal"].getStr());
}

DfDiagonal::~DfDiagonal() {}

void DfDiagonal::run() {
  switch (this->linearAlgebraPackage_) {
    case LAP_LAPACK: {
      this->log_.info("Linear Algebra Package: LAPACK");
      this->run_impl<TlDenseGeneralMatrix_Lapack, TlDenseSymmetricMatrix_Lapack, TlDenseVector_Lapack>();
    } break;

#ifdef HAVE_EIGEN
    case LAP_EIGEN: {
      this->log_.info("Linear Algebra Package: Eigen");
      this->run_impl<TlDenseGeneralMatrix_Eigen, TlDenseSymmetricMatrix_Eigen, TlDenseVector_Eigen>();
    } break;
#endif // HAVE_EIGEN

#ifdef HAVE_VIENNACL
    case LAP_VIENNACL: {
      this->log_.info("Linear Algebra Package: ViennaCL");
      this->run_impl<TlDenseGeneralMatrix_ViennaCL, TlDenseSymmetricMatrix_ViennaCL, TlDenseVector_ViennaCL>();
    } break;
#endif // HAVE_VIENNACL

    default:
      CnErr.abort(TlUtils::format("program error: @%s,%d", __FILE__, __LINE__));
  }
}

// for extended QCLO method
void DfDiagonal::runQclo(DfObject::RUN_TYPE runType,
                            const std::string& fragname, int norbcut) {
  this->m_nNumOfMOs = norbcut;
  this->main<TlDenseGeneralMatrix_Lapack, TlDenseSymmetricMatrix_Lapack, TlDenseVector_Lapack>(
      runType, fragname, true);
}
