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

#include "DfSummary_Parallel.h"
#include "TlCommunicate.h"
#include "tl_dense_general_matrix_scalapack.h"

DfSummary_Parallel::DfSummary_Parallel(TlSerializeData* pPdfParam)
    : DfSummary(pPdfParam) {}

DfSummary_Parallel::~DfSummary_Parallel() {
    TlCommunicate& rComm = TlCommunicate::getInstance();
    rComm.barrier();
}

void DfSummary_Parallel::logger(const std::string& str) const {
    TlCommunicate& rComm = TlCommunicate::getInstance();
    if (rComm.isMaster() == true) {
        DfObject::logger(str);
    }
}

void DfSummary_Parallel::exec() {
#ifdef HAVE_SCALAPACK
    if (this->m_bUsingSCALAPACK == true) {
        this->exec_ScaLAPACK();
        return;
    }
#endif  // HAVE_SCALAPACK

    this->exec_LAPACK();
}

void DfSummary_Parallel::exec_LAPACK() {
    TlCommunicate& rComm = TlCommunicate::getInstance();

    if (rComm.isMaster() == true) {
        DfSummary::exec();
    }
}

void DfSummary_Parallel::exec_ScaLAPACK() {
    DfSummary::exec<TlDenseGeneralMatrix_Scalapack>();
}
