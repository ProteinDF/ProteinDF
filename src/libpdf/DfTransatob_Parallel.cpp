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
#include "TlDistributeMatrix.h"

DfTransatob_Parallel::DfTransatob_Parallel(TlSerializeData* pPdfParam)
    : DfTransatob(pPdfParam) {}

DfTransatob_Parallel::~DfTransatob_Parallel() {}

void DfTransatob_Parallel::logger(const std::string& str) const {
  TlCommunicate& rComm = TlCommunicate::getInstance();
  if (rComm.isMaster() == true) {
    DfTransatob::logger(str);
  }
}

void DfTransatob_Parallel::DfTrsatobMain() {
#ifdef HAVE_SCALAPACK
  if (this->m_bUsingSCALAPACK == true) {
    this->logger("DfTransatob(parallel) using SCALAPACK.\n");
    this->DfTrsatobMain_SCALAPACK();
    return;
  } else {
    this->logger("DfTransatob(parallel) using LAPACK.\n");
  }
#endif  // HAVE_SCALAPACK

  TlCommunicate& rComm = TlCommunicate::getInstance();
  if (rComm.isMaster() == true) {
    DfTransatob::DfTrsatobMain();
  }
  rComm.barrier();
}

void DfTransatob_Parallel::DfTrsatobMain_SCALAPACK() {
  switch (this->m_nMethodType) {
    case METHOD_RKS:
      this->main<TlDistributeMatrix>(RUN_RKS);  // RKS
      break;

    case METHOD_UKS:
      this->main<TlDistributeMatrix>(RUN_UKS_ALPHA);  // UKS alpha spin
      this->main<TlDistributeMatrix>(RUN_UKS_BETA);   // UKS beta spin
      break;

    case METHOD_ROKS:
      this->main<TlDistributeMatrix>(RUN_ROKS);
      break;

    default:
      CnErr.abort();
      break;
  }
}

void DfTransatob_Parallel::DfTrsatobQclo(const std::string& fragname,
                                         int norbcut) {
#ifdef HAVE_SCALAPACK
  if (this->m_bUsingSCALAPACK == true) {
    this->logger("DfTransatob(parallel) using SCALAPACK.\n");
    this->DfTrsatobQclo_SCALAPACK(fragname, norbcut);
    return;
  } else {
    this->logger("DfTransatob(parallel) using LAPACK.\n");
  }
#endif  // HAVE_SCALAPACK

  TlCommunicate& rComm = TlCommunicate::getInstance();
  if (rComm.isMaster() == true) {
    DfTransatob::DfTrsatobQclo(fragname, norbcut);
  }
  rComm.barrier();
}

void DfTransatob_Parallel::DfTrsatobQclo_SCALAPACK(const std::string& fragname,
                                                   int norbcut) {
  this->m_nNumOfMOs = norbcut;

  switch (this->m_nMethodType) {
    case METHOD_RKS:
      this->main<TlDistributeMatrix>(RUN_RKS, fragname, true);  // RKS
      break;

    case METHOD_UKS:
      this->main<TlDistributeMatrix>(RUN_UKS_ALPHA, fragname,
                                     true);  // UKS alpha spin
      this->main<TlDistributeMatrix>(RUN_UKS_BETA, fragname,
                                     true);  // UKS beta spin
      break;

    case METHOD_ROKS:
      this->main<TlDistributeMatrix>(RUN_ROKS, fragname, true);
      break;

    default:
      CnErr.abort();
      break;
  }
}
