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

#include "DfTransatob.h"
#include "CnError.h"
#include "TlUtils.h"
#include "common.h"

DfTransatob::DfTransatob(TlSerializeData* pPdfParam) : DfObject(pPdfParam) {
    this->updateLinearAlgebraPackageParam(
        (*(this->pPdfParam_))["linear_algebra_package/trans_c"].getStr());
}

DfTransatob::~DfTransatob() {}

void DfTransatob::run() {
    switch (this->linearAlgebraPackage_) {
        case LAP_LAPACK: {
            this->log_.info("Linear Algebra Package: LAPACK");
            this->run_method<TlDenseGeneralMatrix_Lapack>(this->m_nNumOfMOs);
        } break;

#ifdef HAVE_EIGEN
        case LAP_EIGEN: {
            this->log_.info("Linear Algebra Package: Eigen");
            this->run_method<TlDenseGeneralMatrix_Eigen>(this->m_nNumOfMOs);
        } break;
#endif  // HAVE_EIGEN

#ifdef HAVE_VIENNACL
        case LAP_VIENNACL: {
            this->log_.info("Linear Algebra Package: ViennaCL");
            this->run_method<TlDenseGeneralMatrix_ViennaCL>(this->m_nNumOfMOs);
        } break;
#endif  // HAVE_VIENNACL

        default:
            CnErr.abort(
                TlUtils::format("program error: @%s,%d", __FILE__, __LINE__));
    }
}

// for extended QCLO method
void DfTransatob::runQclo(const std::string& fragname, int norbcut) {
    switch (this->linearAlgebraPackage_) {
        case LAP_LAPACK: {
            this->log_.info("Linear Algebra Package: LAPACK");
            this->run_method<TlDenseGeneralMatrix_Lapack>(norbcut, fragname);
        } break;

#ifdef HAVE_EIGEN
        case LAP_EIGEN: {
            this->log_.info("Linear Algebra Package: Eigen");
            this->run_method<TlDenseGeneralMatrix_Eigen>(norbcut, fragname);
        } break;
#endif  // HAVE_EIGEN

#ifdef HAVE_VIENNACL
        case LAP_VIENNACL: {
            this->log_.info("Linear Algebra Package: ViennaCL");
            this->run_method<TlDenseGeneralMatrix_ViennaCL>(norbcut, fragname);
        } break;
#endif  // HAVE_VIENNACL

        default:
            CnErr.abort(
                TlUtils::format("program error: @%s,%d", __FILE__, __LINE__));
    }
}
