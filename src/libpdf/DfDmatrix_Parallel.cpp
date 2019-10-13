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
#include "TlFile.h"
#include "TlTime.h"
#include "tl_dense_general_matrix_scalapack.h"
#include "tl_dense_symmetric_matrix_scalapack.h"
#include "tl_dense_vector_scalapack.h"

DfDmatrix_Parallel::DfDmatrix_Parallel(TlSerializeData* pPdfParam)
    : DfDmatrix(pPdfParam) {}

DfDmatrix_Parallel::~DfDmatrix_Parallel() {}

void DfDmatrix_Parallel::run() {
#ifdef HAVE_SCALAPACK
    if (this->m_bUsingSCALAPACK == true) {
        this->run_Scalapack();
    } else {
        TlCommunicate& rComm = TlCommunicate::getInstance();
        if (rComm.isMaster() == true) {
            DfDmatrix::run();
        }
        rComm.barrier();
    }
#else
    {
        TlCommunicate& rComm = TlCommunicate::getInstance();
        if (rComm.isMaster() == true) {
            DfDmatrix::run();
        }
        rComm.barrier();
    }
#endif  // HAVE_SCALAPACK
}

void DfDmatrix_Parallel::run_Scalapack() {
    switch (this->m_nMethodType) {
        case METHOD_RKS:
            this->makeOccupation<TlDenseGeneralMatrix_Scalapack,
                                 TlDenseSymmetricMatrix_Scalapack,
                                 TlDenseVector_Scalapack>(RUN_RKS);
            this->generateDensityMatrix<TlDenseGeneralMatrix_Scalapack,
                                        TlDenseSymmetricMatrix_Scalapack,
                                        TlDenseVector_Scalapack>(RUN_RKS);
            break;

        case METHOD_UKS:
            this->makeOccupation<TlDenseGeneralMatrix_Scalapack,
                                 TlDenseSymmetricMatrix_Scalapack,
                                 TlDenseVector_Scalapack>(RUN_UKS_ALPHA);
            this->generateDensityMatrix<TlDenseGeneralMatrix_Scalapack,
                                        TlDenseSymmetricMatrix_Scalapack,
                                        TlDenseVector_Scalapack>(RUN_UKS_ALPHA);

            this->makeOccupation<TlDenseGeneralMatrix_Scalapack,
                                 TlDenseSymmetricMatrix_Scalapack,
                                 TlDenseVector_Scalapack>(RUN_UKS_BETA);
            this->generateDensityMatrix<TlDenseGeneralMatrix_Scalapack,
                                        TlDenseSymmetricMatrix_Scalapack,
                                        TlDenseVector_Scalapack>(RUN_UKS_BETA);
            break;

        case METHOD_ROKS:
            this->makeOccupation<TlDenseGeneralMatrix_Scalapack,
                                 TlDenseSymmetricMatrix_Scalapack,
                                 TlDenseVector_Scalapack>(RUN_ROKS_CLOSED);
            this->generateDensityMatrix<TlDenseGeneralMatrix_Scalapack,
                                        TlDenseSymmetricMatrix_Scalapack,
                                        TlDenseVector_Scalapack>(
                RUN_ROKS_CLOSED);

            this->makeOccupation<TlDenseGeneralMatrix_Scalapack,
                                 TlDenseSymmetricMatrix_Scalapack,
                                 TlDenseVector_Scalapack>(RUN_ROKS_OPEN);
            this->generateDensityMatrix<TlDenseGeneralMatrix_Scalapack,
                                        TlDenseSymmetricMatrix_Scalapack,
                                        TlDenseVector_Scalapack>(RUN_ROKS_OPEN);

            // ROKS_alpha,beta
            {
                TlDenseSymmetricMatrix_Scalapack PC =
                    DfObject::getSpinDensityMatrix<
                        TlDenseSymmetricMatrix_Scalapack>(RUN_ROKS_CLOSED,
                                                          this->m_nIteration);
                TlDenseSymmetricMatrix_Scalapack PO =
                    DfObject::getSpinDensityMatrix<
                        TlDenseSymmetricMatrix_Scalapack>(RUN_ROKS_OPEN,
                                                          this->m_nIteration);

                TlDenseSymmetricMatrix_Scalapack PA = PC + PO;
                DfObject::saveSpinDensityMatrix(RUN_ROKS_ALPHA,
                                                this->m_nIteration, PA);

                TlDenseSymmetricMatrix_Scalapack PB = PC;
                DfObject::saveSpinDensityMatrix(RUN_ROKS_BETA,
                                                this->m_nIteration, PB);
            }
            break;

        default:
            CnErr.abort();
            break;
    }
}

// void DfDmatrix_Parallel::main_SCALAPACK(const DfObject::RUN_TYPE runType) {
//   TlCommunicate& rComm = TlCommunicate::getInstance();
//   this->log_.info("build density matrix using ScaLAPACK.");

//   // occupation
//   TlDenseVector_Lapack currOcc;
//   switch (this->orbitalCorrespondenceMethod_) {
//     case OCM_OVERLAP:
//       this->log_.info(" orbital correspondence method: MO-overlap");
//       currOcc =
//       this->getOccupationUsingOverlap<TlDenseGeneralMatrix_Scalapack>(
//           runType);
//       if (rComm.isMaster() == true) {
//         currOcc.save(this->getOccupationPath(runType));
//       }
//       break;

//     case OCM_PROJECTION:
//       this->log_.info(" orbital correspondence method: MO-projection");
//       currOcc =
//           this->getOccupationUsingProjection<TlDenseGeneralMatrix_Scalapack,
//                                              TlDenseSymmetricMatrix_Scalapack>(
//               runType);
//       if (rComm.isMaster() == true) {
//         currOcc.save(this->getOccupationPath(runType));
//       }
//       break;

//     default:
//       if (rComm.isMaster() == true) {
//         this->log_.info(" orbital correspondence method: none");
//         currOcc.load(this->getOccupationPath(runType));
//       }
//       rComm.broadcast(&currOcc);
//       break;
//   }

//   rComm.barrier();
//   this->generateDensityMatrix<TlDenseGeneralMatrix_Scalapack,
//                               TlDenseSymmetricMatrix_Scalapack>(runType,
//                                                                 currOcc);
// }

void DfDmatrix_Parallel::checkOccupation(const TlDenseVector_Lapack& prevOcc,
                                         const TlDenseVector_Lapack& currOcc) {
    TlCommunicate& rComm = TlCommunicate::getInstance();

    if (rComm.isMaster() == true) {
        DfDmatrix::checkOccupation(prevOcc, currOcc);
    }
}

void DfDmatrix_Parallel::printOccupation(const TlDenseVector_Lapack& occ) {
    TlCommunicate& rComm = TlCommunicate::getInstance();

    if (rComm.isMaster() == true) {
        DfDmatrix::printOccupation(occ);
    }
}

TlDenseVector_Lapack DfDmatrix_Parallel::getOccVtr(
    const DfObject::RUN_TYPE runType) {
    TlCommunicate& rComm = TlCommunicate::getInstance();
    TlDenseVector_Lapack occ;
    if (rComm.isMaster() == true) {
        occ = DfObject::getOccVtr<TlDenseVector_Lapack>(runType);
    }
    rComm.broadcast(&occ);

    return occ;
}
