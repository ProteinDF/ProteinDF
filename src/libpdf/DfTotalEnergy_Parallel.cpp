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
#include "DfEriX_Parallel.h"
#include "DfOverlapX_Parallel.h"
#include "DfTotalEnergy_Parallel.h"
#include "DfXCFunctional_Parallel.h"
#include "Fl_Geometry.h"
#include "TlCommunicate.h"

DfTotalEnergy_Parallel::DfTotalEnergy_Parallel(TlSerializeData* pPdfParam)
    : DfTotalEnergy(pPdfParam) {
  this->m_bUseDistributeMatrix = false;
}

DfTotalEnergy_Parallel::~DfTotalEnergy_Parallel() {}

void DfTotalEnergy_Parallel::logger(const std::string& str) const {
  TlCommunicate& rComm = TlCommunicate::getInstance();

  if (rComm.isMaster() == true) {
    DfTotalEnergy::logger(str);
  }
}

void DfTotalEnergy_Parallel::output() {
  TlCommunicate& rComm = TlCommunicate::getInstance();
  if (rComm.isMaster() == true) {
    DfTotalEnergy::output();
  }
}

void DfTotalEnergy_Parallel::exec() {
#ifdef HAVE_SCALAPACK
  if (this->m_bUsingSCALAPACK == true) {
    this->exec_ScaLAPACK();
    return;
  }
#endif  // HAVE_SCALAPACK

  this->exec_LAPACK();
}

void DfTotalEnergy_Parallel::exec_LAPACK() {
  TlCommunicate& rComm = TlCommunicate::getInstance();
  if (rComm.isMaster() == true) {
    DfTotalEnergy::exec();
  }
  rComm.barrier();
}

void DfTotalEnergy_Parallel::exec_ScaLAPACK() {
  DfTotalEnergy::exec_template<DfOverlapX_Parallel, DfEriX_Parallel,
                               TlDistributeSymmetricMatrix,
                               TlDistributeVector>();
}

// total energy including dummy charge
void DfTotalEnergy_Parallel::calculate_real_energy() {
#ifdef HAVE_SCALAPACK
  if (this->m_bUsingSCALAPACK == true) {
    // ScaLAPACK
    this->calcRealEnergy<TlDistributeSymmetricMatrix>();
    return;
  }
#endif  // HAVE_SCALAPACK

  // LAPACK
  TlCommunicate& rComm = TlCommunicate::getInstance();
  if (rComm.isMaster() == true) {
    DfTotalEnergy::calculate_real_energy();
  }
}
