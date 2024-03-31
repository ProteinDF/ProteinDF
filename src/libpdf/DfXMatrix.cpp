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

#include "DfXMatrix.h"

#include <cassert>
#include <cmath>
#include <fstream>
#include <string>

#include "CnError.h"
#include "tl_dense_general_matrix_lapack.h"
#include "tl_dense_symmetric_matrix_lapack.h"
#include "tl_dense_vector_lapack.h"

#ifdef HAVE_EIGEN
#include "tl_dense_general_matrix_eigen.h"
#include "tl_dense_symmetric_matrix_eigen.h"
#include "tl_dense_vector_eigen.h"
#endif  // HAVE_EIGEN

DfXMatrix::DfXMatrix(TlSerializeData* pPdfParam)
    : DfObject(pPdfParam) {
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

    this->debugSaveEigval_ =
        pdfParam["debug/DfXMatrix/save-eigval"].getBoolean();
    this->debugSaveMatrix_ = pdfParam["debug/DfXMatrix/save-mat"].getBoolean();
    this->debugCheckX_ = pdfParam["debug/DfXMatrix/check-X"].getBoolean();
}

DfXMatrix::~DfXMatrix() {}

void DfXMatrix::buildX() {
    switch (this->linearAlgebraPackage_) {
        case LAP_LAPACK: {
            this->log_.info("Linear Algebra Package: Lapack");
            this->buildX_templ<TlDenseSymmetricMatrix_Lapack,
                               TlDenseGeneralMatrix_Lapack, TlDenseVector_Lapack>();
        } break;

#ifdef HAVE_EIGEN
        case LAP_EIGEN: {
            this->log_.info("Linear Algebra Package: Eigen");
            this->buildX_templ<TlDenseSymmetricMatrix_Eigen,
                               TlDenseGeneralMatrix_Eigen, TlDenseVector_Eigen>();
        } break;
#endif  // HAVE_EIGEN

#ifdef HAVE_VIENNACL
        case LAP_VIENNACL: {
            this->log_.info("Linear Algebra Package: ViennaCL -> Eigen");
            this->buildX_templ<TlDenseSymmetricMatrix_Eigen,
                               TlDenseGeneralMatrix_Eigen, TlDenseVector_Eigen>();
        } break;
#endif  // HAVE_VIENNACL

        default:
            CnErr.abort(TlUtils::format("program error: @%s,%d", __FILE__, __LINE__));
    }
}

void DfXMatrix::canonicalOrthogonalize(const TlDenseSymmetricMatrix_Lapack& S,
                                       TlDenseGeneralMatrix_Lapack* pX,
                                       TlDenseGeneralMatrix_Lapack* pXinv,
                                       const std::string& eigvalFilePath) {
    this->canonicalOrthogonalizeTmpl<TlDenseSymmetricMatrix_Lapack,
                                     TlDenseGeneralMatrix_Lapack, TlDenseVector_Lapack>(
        S, pX, pXinv, eigvalFilePath);
}

void DfXMatrix::lowdinOrthogonalize(const TlDenseSymmetricMatrix_Lapack& S,
                                    TlDenseGeneralMatrix_Lapack* pX,
                                    TlDenseGeneralMatrix_Lapack* pXinv,
                                    const std::string& eigvalFilePath) {
    this->lowdinOrthogonalizeTmpl<TlDenseSymmetricMatrix_Lapack,
                                  TlDenseGeneralMatrix_Lapack>(S, pX, pXinv,
                                                               eigvalFilePath);
}
