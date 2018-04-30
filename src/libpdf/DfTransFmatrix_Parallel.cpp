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

#include "DfTransFmatrix_Parallel.h"
#include "CnError.h"
#include "TlCommunicate.h"
#include "TlDistributeMatrix.h"
#include "TlDistributeSymmetricMatrix.h"

DfTransFmatrix_Parallel::DfTransFmatrix_Parallel(TlSerializeData* pPdfParam,
                                                 bool bExecDiis)
    : DfTransFmatrix(pPdfParam, bExecDiis) {}

DfTransFmatrix_Parallel::~DfTransFmatrix_Parallel() {}

void DfTransFmatrix_Parallel::logger(const std::string& str) const {
  TlCommunicate& rComm = TlCommunicate::getInstance();
  if (rComm.isMaster() == true) {
    DfTransFmatrix::logger(str);
  }
}

void DfTransFmatrix_Parallel::DfTrsFmatMain() {
#ifdef HAVE_SCALAPACK
  if (this->m_bUsingSCALAPACK == true) {
    this->logger("DfTransFmatrix(parallel) using SCALAPACK.\n");
    this->DfTrsFmatMain_SCALAPACK();
    return;
  } else {
    this->logger("DfTransFmatrix(parallel) using LAPACK.\n");
  }
#endif  // HAVE_SCALAPACK

  TlCommunicate& rComm = TlCommunicate::getInstance();
  if (rComm.isMaster() == true) {
    DfTransFmatrix::DfTrsFmatMain();
  }
  rComm.barrier();
}

void DfTransFmatrix_Parallel::DfTrsFmatMain_SCALAPACK() {
  switch (this->m_nMethodType) {
    case METHOD_RKS:
      DfTransFmatrix::main<TlDistributeMatrix, TlDistributeSymmetricMatrix>(
          RUN_RKS);
      break;

    case METHOD_UKS:
      DfTransFmatrix::main<TlDistributeMatrix, TlDistributeSymmetricMatrix>(
          RUN_UKS_ALPHA);
      DfTransFmatrix::main<TlDistributeMatrix, TlDistributeSymmetricMatrix>(
          RUN_UKS_BETA);
      break;

    case METHOD_ROKS:
      DfTransFmatrix::main<TlDistributeMatrix, TlDistributeSymmetricMatrix>(
          RUN_ROKS);
      break;

    default:
      CnErr.abort();
      break;
  }
}

void DfTransFmatrix_Parallel::DfTrsFmatQclo(const std::string& fragname,
                                            int norbcut) {
#ifdef HAVE_SCALAPACK
  if (this->m_bUsingSCALAPACK == true) {
    this->logger("DfTransFmatrix(parallel) using SCALAPACK.\n");
    this->DfTrsFmatQclo_SCALAPACK(fragname, norbcut);
    return;
  } else {
    this->logger("DfTransFmatrix(parallel) using LAPACK.\n");
  }
#endif  // HAVE_SCALAPACK

  TlCommunicate& rComm = TlCommunicate::getInstance();
  if (rComm.isMaster() == true) {
    DfTransFmatrix::DfTrsFmatQclo(fragname, norbcut);
  }
  rComm.barrier();
}

void DfTransFmatrix_Parallel::DfTrsFmatQclo_SCALAPACK(
    const std::string& fragname, int norbcut) {
  this->m_nNumOfMOs = norbcut;
  switch (this->m_nMethodType) {
    case METHOD_RKS:
      DfTransFmatrix::main<TlDistributeMatrix, TlDistributeSymmetricMatrix>(
          RUN_RKS, fragname, true);
      break;

    case METHOD_UKS:
      DfTransFmatrix::main<TlDistributeMatrix, TlDistributeSymmetricMatrix>(
          RUN_UKS_ALPHA, fragname, true);
      DfTransFmatrix::main<TlDistributeMatrix, TlDistributeSymmetricMatrix>(
          RUN_UKS_BETA, fragname, true);
      break;

    case METHOD_ROKS:
      DfTransFmatrix::main<TlDistributeMatrix, TlDistributeSymmetricMatrix>(
          RUN_ROKS, fragname, true);
      break;

    default:
      CnErr.abort();
      break;
  }
}
