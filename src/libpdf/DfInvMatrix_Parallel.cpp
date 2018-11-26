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

#include "DfInvMatrix_Parallel.h"
#include "TlCommunicate.h"
#include "tl_dense_symmetric_matrix_scalapack.h"
#include "tl_dense_symmetric_matrix_lapack.h"

DfInvMatrix_Parallel::DfInvMatrix_Parallel(TlSerializeData* pPdfParam)
    : DfInvMatrix(pPdfParam) {}

DfInvMatrix_Parallel::~DfInvMatrix_Parallel() {}

void DfInvMatrix_Parallel::DfInvMain() {
#ifdef HAVE_SCALAPACK
  if (this->m_bUsingSCALAPACK == true) {
    this->exec<TlDenseSymmetricMatrix_Scalapack>();
    return;
  }
#endif  // HAVE_SCALAPACK

  TlCommunicate& rComm = TlCommunicate::getInstance();
  if (rComm.isMaster() == true) {
    DfInvMatrix::DfInvMain();
  }
  rComm.barrier();
}
