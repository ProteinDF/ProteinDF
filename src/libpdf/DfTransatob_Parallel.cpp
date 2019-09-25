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

#include "DfTransatob_Parallel.h"
#include "CnError.h"
#include "TlCommunicate.h"
#include "tl_dense_general_matrix_scalapack.h"

DfTransatob_Parallel::DfTransatob_Parallel(TlSerializeData* pPdfParam)
    : DfTransatob(pPdfParam) {}

DfTransatob_Parallel::~DfTransatob_Parallel() {}

void DfTransatob_Parallel::logger(const std::string& str) const {
    TlCommunicate& rComm = TlCommunicate::getInstance();
    if (rComm.isMaster() == true) {
        DfTransatob::logger(str);
    }
}

void DfTransatob_Parallel::run() {
#ifdef HAVE_SCALAPACK
    if (this->m_bUsingSCALAPACK == true) {
        this->logger("DfTransatob(parallel) using SCALAPACK.");
        this->run_Scalapack();
        return;
    } else {
        this->logger("DfTransatob(parallel) using LAPACK.");
    }
#endif  // HAVE_SCALAPACK

    TlCommunicate& rComm = TlCommunicate::getInstance();
    if (rComm.isMaster() == true) {
        DfTransatob::run();
    }
    rComm.barrier();
}

void DfTransatob_Parallel::run_Scalapack() {
    DfTransatob::run_method<TlDenseGeneralMatrix_Scalapack>(this->m_nNumOfMOs);
}

void DfTransatob_Parallel::runQclo(const std::string& fragname, int norbcut) {
#ifdef HAVE_SCALAPACK
    if (this->m_bUsingSCALAPACK == true) {
        this->logger("DfTransatob(parallel) using SCALAPACK.\n");
        this->runQclo_Scalapack(fragname, norbcut);
        return;
    } else {
        this->logger("DfTransatob(parallel) using LAPACK.\n");
    }
#endif  // HAVE_SCALAPACK

    TlCommunicate& rComm = TlCommunicate::getInstance();
    if (rComm.isMaster() == true) {
        DfTransatob::runQclo(fragname, norbcut);
    }
    rComm.barrier();
}

void DfTransatob_Parallel::runQclo_Scalapack(const std::string& fragname,
                                             int norbcut) {
    DfTransatob::run_method<TlDenseGeneralMatrix_Scalapack>(norbcut, fragname);
}
