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

#include "DfPopulation_Parallel.h"
#include "TlCommunicate.h"
#include "tl_dense_general_matrix_scalapack.h"
#include "tl_dense_symmetric_matrix_scalapack.h"

DfPopulation_Parallel::DfPopulation_Parallel(TlSerializeData* pPdfParam)
    : DfPopulation(pPdfParam) {}

DfPopulation_Parallel::~DfPopulation_Parallel() {
  TlCommunicate& rComm = TlCommunicate::getInstance();
  rComm.barrier();
}

double DfPopulation_Parallel::getSumOfElectrons(
    const TlDenseSymmetricMatrix_Lapack& P) {
  double answer = 0.0;
  TlCommunicate& rComm = TlCommunicate::getInstance();
  if (rComm.isMaster() == true) {
    answer = DfPopulation::getSumOfElectrons(P);
  }
  rComm.broadcast(answer);
  return answer;
}

double DfPopulation_Parallel::getSumOfElectrons(
    const TlDenseSymmetricMatrix_Scalapack& P) {
  const TlDenseVector_Lapack trPS =
      this->getPS<TlDenseGeneralMatrix_Scalapack,
                  TlDenseSymmetricMatrix_Scalapack>(P);
  return trPS.sum();
}

void DfPopulation_Parallel::calcPop(const int iteration) {
#ifdef HAVE_SCALAPACK
  if (this->m_bUsingSCALAPACK == true) {
    DfPopulation::calcPop<TlDenseGeneralMatrix_Scalapack,
                          TlDenseSymmetricMatrix_Scalapack>(iteration);
    return;
  }
#endif  // HAVE_SCALAPACK

  TlCommunicate& rComm = TlCommunicate::getInstance();
  if (rComm.isMaster() == true) {
    DfPopulation::calcPop(iteration);
  }
  rComm.barrier();
}