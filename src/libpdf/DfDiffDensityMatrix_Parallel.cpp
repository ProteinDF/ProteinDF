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

#include "DfDiffDensityMatrix_Parallel.h"
#include "TlCommunicate.h"
#include "TlFile.h"
#include "tl_dense_symmetric_matrix_blacs.h"

DfDiffDensityMatrix_Parallel::DfDiffDensityMatrix_Parallel(
    TlSerializeData* pPdfParam)
    : DfDiffDensityMatrix(pPdfParam) {}

DfDiffDensityMatrix_Parallel::~DfDiffDensityMatrix_Parallel() {}

void DfDiffDensityMatrix_Parallel::exec() {
  TlCommunicate& rComm = TlCommunicate::getInstance();

#ifdef HAVE_SCALAPACK
  if (this->m_bUsingSCALAPACK == true) {
    // using ScaLAPACK
    switch (this->m_nMethodType) {
      case METHOD_RKS:
        DfDiffDensityMatrix::calc<TlDenseSymmetricMatrix_blacs>(
            RUN_RKS, this->m_nIteration);
        break;

      case METHOD_UKS:
        DfDiffDensityMatrix::calc<TlDenseSymmetricMatrix_blacs>(
            RUN_UKS_ALPHA, this->m_nIteration);
        DfDiffDensityMatrix::calc<TlDenseSymmetricMatrix_blacs>(
            RUN_UKS_BETA, this->m_nIteration);
        break;

      case METHOD_ROKS:
        DfDiffDensityMatrix::calc<TlDenseSymmetricMatrix_blacs>(
            RUN_ROKS_CLOSED, this->m_nIteration);
        DfDiffDensityMatrix::calc<TlDenseSymmetricMatrix_blacs>(
            RUN_ROKS_OPEN, this->m_nIteration);
        break;

      default:
        std::cerr << "program error. " << __FILE__ << ":" << __LINE__
                  << std::endl;
        abort();
    }

  } else {
    // using LAPACK
    if (rComm.isMaster() == true) {
      DfDiffDensityMatrix::exec();
    }
    rComm.barrier();
  }
#else
  {
    // using LAPACK
    if (rComm.isMaster() == true) {
      DfDiffDensityMatrix::exec();
    }
    rComm.barrier();
  }
#endif  // HAVE_SCALAPACK
}

// void DfDiffDensityMatrix_Parallel::calc_usingScaLAPACK(const
// DfObject::RUN_TYPE runType,
//                                                        const int iteration)
// {
//     this->log_.info("(delta P) is build using ScaLAPACK.");
//     TlDenseSymmetricMatrix_blacs P =
//     DfObject::getPpqMatrix<TlDenseSymmetricMatrix_blacs>(runType, iteration
//     -1); P.save(TlUtils::format("diffP_%d.mat", iteration)); if
//     (TlFile::isExist(this->getPpqMatrixPath(runType, iteration -2)) == true)
//     {
//         P -= (this->getPpqMatrix<TlDenseSymmetricMatrix_blacs>(runType,
//         iteration -2));
//     }

//     this->saveDiffDensityMatrix(runType, iteration, P);
// }
