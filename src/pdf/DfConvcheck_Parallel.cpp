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

#include <iostream>
#include "DfConvcheck_Parallel.h"
#include "TlCommunicate.h"
#include "TlDistributeSymmetricMatrix.h"

DfConvcheck_Parallel::DfConvcheck_Parallel(TlSerializeData* pPdfParam, int num_iter)
    : DfConvcheck(pPdfParam, num_iter)
{
//   TlCommunicate& rComm = TlCommunicate::getInstance();
//   std::cout << "DfConvcheck_Parallel::DfConvcheck_Parallel() at " << rComm.getRank() << std::endl;
}

DfConvcheck_Parallel::~DfConvcheck_Parallel()
{
//   TlCommunicate& rComm = TlCommunicate::getInstance();
//   std::cout << "DfConvcheck_Parallel::~DfConvcheck_Parallel() at " << rComm.getRank() << std::endl;
}


void DfConvcheck_Parallel::DfConvcheckMain()
{
#ifdef HAVE_SCALAPACK
    if (this->m_bUsingSCALAPACK == true) {
        this->log_.info("convergence check using SCALAPACK.");
        this->main_ScaLAPACK();
        return;
    }
#endif // HAVE_SCALAPACK

    // LAPACK
    this->log_.info("convergence check using LAPACK.");
    TlCommunicate& rComm = TlCommunicate::getInstance();
    if (rComm.isMaster() == true) {
        DfConvcheck::DfConvcheckMain();
    }
    rComm.broadcast(this->converged_flag);
}


void DfConvcheck_Parallel::main_ScaLAPACK() {
    TlCommunicate& rComm = TlCommunicate::getInstance();
    
    this->converged_flag = 0;

    // no judgement for the first iteration
    if (this->m_nIteration == 1) {
        return;
    }
    this->main<TlDistributeSymmetricMatrix>(this->m_nIteration);

    // check for convergence
    double return_threshold = 0.0;
    if (this->convergence_type == "fock") {
        return_threshold = this->dev_ks;
    } else if (this->convergence_type == "density") {
        return_threshold = this->dev_dm;
    } else if (this->convergence_type == "energy") {
        return_threshold = this->dev_te;
    } else {
        return_threshold = this->dev_cd;
    }

    if ((return_threshold < this->threshold_cri) && (dev_te < this->threshold_cri_ene)) {
        this->converged_flag = 1;
    }

    if (rComm.isMaster() == true) {
        this->showResults();
    }
}
