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

#include "DfInitialGuessHarris_Parallel.h"
#include "DfOverlapX_Parallel.h"
#include "DfPopulation_Parallel.h"
#include "TlCommunicate.h"
#include "tl_dense_general_matrix_scalapack.h"
#include "tl_dense_symmetric_matrix_scalapack.h"

DfInitialGuessHarris_Parallel::DfInitialGuessHarris_Parallel(
    TlSerializeData* pPdfParam)
    : DfInitialGuessHarris(pPdfParam) {}

DfInitialGuessHarris_Parallel::~DfInitialGuessHarris_Parallel() {}

void DfInitialGuessHarris_Parallel::main() {
  TlCommunicate& rComm = TlCommunicate::getInstance();

#ifdef HAVE_SCALAPACK
  if (this->m_bUsingSCALAPACK) {
    this->distributeHarrisDB();

    DfInitialGuessHarris::calcInitialDensityMatrix<
        TlDenseGeneralMatrix_Scalapack, TlDenseSymmetricMatrix_Scalapack,
        DfOverlapX_Parallel, DfPopulation_Parallel>();
  } else {
    if (rComm.isMaster()) {
      DfInitialGuessHarris::main();
    }
  }
#else
  if (rComm.isMaster()) {
    DfInitialGuessHarris::main();
  }
#endif  // HAVE_SCALAPACK
}

void DfInitialGuessHarris_Parallel::distributeHarrisDB() {
  TlCommunicate& rComm = TlCommunicate::getInstance();

  if (rComm.isMaster()) {
    DfInitialGuessHarris::loadHarrisDB();
  }

  rComm.broadcast(this->pdfParam_harrisDB_);
}
