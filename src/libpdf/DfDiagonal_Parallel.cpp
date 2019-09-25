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
#include "DfDiagonal_Parallel.h"
#include "TlCommunicate.h"
#include "TlUtils.h"
#include "tl_dense_general_matrix_scalapack.h"
#include "tl_dense_symmetric_matrix_scalapack.h"
#include "tl_dense_vector_scalapack.h"

DfDiagonal_Parallel::DfDiagonal_Parallel(TlSerializeData* pPdfParam)
    : DfDiagonal(pPdfParam) {}

DfDiagonal_Parallel::~DfDiagonal_Parallel() {}

void DfDiagonal_Parallel::run() {
#ifdef HAVE_SCALAPACK
    if (this->m_bUsingSCALAPACK == true) {
        this->log_.info("diagonal(parallel) using SCALAPACK.");
        this->run_SCALAPACK();
        return;
    }
#endif  // HAVE_SCALAPACK

    // LAPACK
    this->log_.info("diagonal(parallel) using LAPACK.");
    TlCommunicate& rComm = TlCommunicate::getInstance();
    if (rComm.isMaster() == true) {
        DfDiagonal::run();
    }
    rComm.barrier();
}

void DfDiagonal_Parallel::DfDiagQclo(const DfObject::RUN_TYPE runType,
                                     const std::string& fragname, int norbcut) {
#ifdef HAVE_SCALAPACK
    if (this->m_bUsingSCALAPACK == true) {
        this->log_.info("diagonalization using SCALAPACK.");
        this->DfDiagQclo_SCALAPACK(runType, fragname, norbcut);
        return;
    } else {
        this->log_.info("diagonalization using LAPACK.");
    }
#endif  // HAVE_SCALAPACK

    TlCommunicate& rComm = TlCommunicate::getInstance();

    if (rComm.isMaster() == true) {
        DfDiagonal::runQclo(runType, fragname, norbcut);
    }
    rComm.barrier();
}

void DfDiagonal_Parallel::run_SCALAPACK() {
    switch (this->m_nMethodType) {
        case METHOD_RKS:
            this->main<TlDenseGeneralMatrix_Scalapack,
                       TlDenseSymmetricMatrix_Scalapack, TlDenseVector_Lapack>(
                RUN_RKS);
            break;

        case METHOD_UKS:
            this->main<TlDenseGeneralMatrix_Scalapack,
                       TlDenseSymmetricMatrix_Scalapack, TlDenseVector_Lapack>(
                RUN_UKS_ALPHA);
            this->main<TlDenseGeneralMatrix_Scalapack,
                       TlDenseSymmetricMatrix_Scalapack, TlDenseVector_Lapack>(
                RUN_UKS_BETA);
            break;

        case METHOD_ROKS:
            this->main<TlDenseGeneralMatrix_Scalapack,
                       TlDenseSymmetricMatrix_Scalapack, TlDenseVector_Lapack>(
                RUN_ROKS);
            break;

        default:
            CnErr.abort("DfDiagonal", "", "DfDiagMain",
                        "the value of scftype is illegal");
            break;
    }
}

void DfDiagonal_Parallel::DfDiagQclo_SCALAPACK(const DfObject::RUN_TYPE runType,
                                               const std::string& fragname,
                                               int norbcut) {
    this->m_nNumOfMOs = norbcut;
    this->main<TlDenseGeneralMatrix_Scalapack, TlDenseSymmetricMatrix_Scalapack,
               TlDenseVector_Lapack>(runType, fragname, true);
}
