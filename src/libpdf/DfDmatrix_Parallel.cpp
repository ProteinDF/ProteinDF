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

#include "DfDmatrix_Parallel.h"
#include "TlCommunicate.h"
#include "TlDistributeMatrix.h"
#include "TlDistributeSymmetricMatrix.h"
#include "TlFile.h"
#include "TlTime.h"

DfDmatrix_Parallel::DfDmatrix_Parallel(TlSerializeData* pPdfParam)
    : DfDmatrix(pPdfParam) {}

DfDmatrix_Parallel::~DfDmatrix_Parallel() {}

void DfDmatrix_Parallel::main(const DfObject::RUN_TYPE runType) {
#ifdef HAVE_SCALAPACK
  if (this->m_bUsingSCALAPACK == true) {
    this->main_SCALAPACK(runType);
  } else {
    TlCommunicate& rComm = TlCommunicate::getInstance();
    if (rComm.isMaster() == true) {
      DfDmatrix::main(runType);
    }
    rComm.barrier();
  }
#else
  {
    TlCommunicate& rComm = TlCommunicate::getInstance();
    if (rComm.isMaster() == true) {
      DfDmatrix::main(runType);
    }
    rComm.barrier();
  }
#endif  // HAVE_SCALAPACK
}

void DfDmatrix_Parallel::main_SCALAPACK(const DfObject::RUN_TYPE runType) {
  TlCommunicate& rComm = TlCommunicate::getInstance();
  this->log_.info("build density matrix using ScaLAPACK.");

  // occupation
  TlVector currOcc;
  switch (this->orbitalCorrespondenceMethod_) {
    case OCM_OVERLAP:
      this->log_.info(" orbital correspondence method: MO-overlap");
      currOcc = this->getOccupationUsingOverlap<TlDistributeMatrix>(runType);
      if (rComm.isMaster() == true) {
        currOcc.save(this->getOccupationPath(runType));
      }
      break;

    case OCM_PROJECTION:
      this->log_.info(" orbital correspondence method: MO-projection");
      currOcc = this->getOccupationUsingProjection<TlDistributeMatrix,
                                                   TlDistributeSymmetricMatrix>(
          runType);
      if (rComm.isMaster() == true) {
        currOcc.save(this->getOccupationPath(runType));
      }
      break;

    default:
      if (rComm.isMaster() == true) {
        this->log_.info(" orbital correspondence method: none");
        currOcc.load(this->getOccupationPath(runType));
      }
      rComm.broadcast(currOcc);
      break;
  }

  rComm.barrier();
  this->generateDensityMatrix<TlDistributeMatrix, TlDistributeSymmetricMatrix>(
      runType, currOcc);
}

void DfDmatrix_Parallel::checkOccupation(const TlVector& prevOcc,
                                         const TlVector& currOcc) {
  TlCommunicate& rComm = TlCommunicate::getInstance();

  if (rComm.isMaster() == true) {
    DfDmatrix::checkOccupation(prevOcc, currOcc);
  }
}

void DfDmatrix_Parallel::printOccupation(const TlVector& occ) {
  TlCommunicate& rComm = TlCommunicate::getInstance();

  if (rComm.isMaster() == true) {
    DfDmatrix::printOccupation(occ);
  }
}

TlVector DfDmatrix_Parallel::getOccVtr(const DfObject::RUN_TYPE runType) {
  TlCommunicate& rComm = TlCommunicate::getInstance();
  TlVector occ;
  if (rComm.isMaster() == true) {
    occ = DfObject::getOccVtr(runType);
  }
  rComm.broadcast(occ);

  return occ;
}
