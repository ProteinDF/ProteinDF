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

#include "DfConvcheck_Parallel.h"
#include <iostream>
#include "TlCommunicate.h"
#include "TlDistributeSymmetricMatrix.h"
#include "TlSymmetricMatrix.h"

DfConvcheck_Parallel::DfConvcheck_Parallel(TlSerializeData* pPdfParam,
                                           int num_iter)
    : DfConvcheck(pPdfParam, num_iter) {
  //   TlCommunicate& rComm = TlCommunicate::getInstance();
  //   std::cout << "DfConvcheck_Parallel::DfConvcheck_Parallel() at " <<
  //   rComm.getRank() << std::endl;
}

DfConvcheck_Parallel::~DfConvcheck_Parallel() {
  //   TlCommunicate& rComm = TlCommunicate::getInstance();
  //   std::cout << "DfConvcheck_Parallel::~DfConvcheck_Parallel() at " <<
  //   rComm.getRank() << std::endl;
}

void DfConvcheck_Parallel::check() {
  TlCommunicate& rComm = TlCommunicate::getInstance();

  // no judgement for the first iteration
  if (this->m_nIteration == 1) {
    this->isConverged_ = false;
    return;
  }

#ifdef HAVE_SCALAPACK
  if (this->m_bUsingSCALAPACK) {
    this->log_.info("convgergence check (parallel) using ScaLAPACK");
    DfConvcheck::check<TlDistributeSymmetricMatrix>(this->m_nIteration);
  } else {
    this->log_.info("convgergence check (parallel) using LAPACK");
    if (rComm.isMaster()) {
      DfConvcheck::check<TlSymmetricMatrix>(this->m_nIteration);
    }
  }
#else
  {
    this->log_.info("convgergence check (parallel) using LAPACK");
    if (rComm.isMaster()) {
      DfConvcheck::check<TlSymmetricMatrix>(this->m_nIteration);
    }
  }
#endif  // HAVE_SCALAPACK

  // check for convergence
  if (rComm.isMaster()) {
    this->isConverged_ =
        ((this->judgeRmsMatrixA_) && (this->judgeRmsMatrixB_) &&
         (this->judgeMaxMatrixA_) && (this->judgeMaxMatrixB_) &&
         (this->judgeTotalEnergy_));
  }
  rComm.broadcast(this->isConverged_);

  this->showResults();
}
