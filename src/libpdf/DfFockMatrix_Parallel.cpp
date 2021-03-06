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

#include "DfFockMatrix_Parallel.h"
#include "DfEriX_Parallel.h"
#include "DfOverlapX_Parallel.h"
#include "TlCommunicate.h"
#include "tl_dense_symmetric_matrix_scalapack.h"

DfFockMatrix_Parallel::DfFockMatrix_Parallel(TlSerializeData* pPdfParam)
    : DfFockMatrix(pPdfParam) {
    // std::cerr << "DfFockMatrix_Parallel::DfFockMatrix_Parallel()" <<
    // std::endl;
}

DfFockMatrix_Parallel::~DfFockMatrix_Parallel() {
    //     std::cerr << "DfFockMatrix_Parallel::~DfFockMatrix_Parallel()" <<
    //     std::endl;
}

void DfFockMatrix_Parallel::logger(const std::string& str) const {
    TlCommunicate& rComm = TlCommunicate::getInstance();
    if (rComm.isMaster() == true) {
        DfFockMatrix::logger(str);
    }
}

void DfFockMatrix_Parallel::mainDIRECT_RKS() {
    TlCommunicate& rComm = TlCommunicate::getInstance();

    if (this->m_bUsingSCALAPACK == true) {
        // ScaLAPACK
        this->log_.info("buid KS matrix using ScaLAPACK.");
        DfFockMatrix::mainDIRECT_RKS<TlDenseSymmetricMatrix_Scalapack>();
    } else {
        // LAPACK
        this->log_.info("buid KS matrix using LAPACK.");
        if (rComm.isMaster() == true) {
            DfFockMatrix::mainDIRECT_RKS();
        }
    }
}

// void DfFockMatrix_Parallel::mainDIRECT_RKS_LAPACK()
// {
//     TlCommunicate& rComm = TlCommunicate::getInstance();

//     if (rComm.isMaster() == true) {
//         this->logger("Direct scheme method is employed\n");
//     }

//     TlDenseSymmetricMatrix_Lapack F(this->m_nNumOfAOs);
//     if (this->m_bMemorySave == true) {
//         // DfThreeindexintegrals を使わない
//         if (this->m_bIsXCFitting == true) {
//             this->setXC_RI(RUN_RKS, F);
//         } else {
//             this->setXC_DIRECT(RUN_RKS, F);
//         }

//         this->setCoulomb(METHOD_RKS, F);
//     } else {
//         // DfThreeindexintegrals を使う
//         assert(this->m_bIsXCFitting == true);
//         F = this->getFpqMatrix(RUN_RKS, this->m_nIteration);
//     }

//     if (rComm.isMaster() == true) {
//         this->setHpq(RUN_RKS, F);
//         DfObject::saveFpqMatrix(RUN_RKS, this->m_nIteration, F);
//     }
// }

// void DfFockMatrix_Parallel::mainDIRECT_RKS_ScaLAPACK()
// {
//     TlCommunicate& rComm = TlCommunicate::getInstance();

//     if (rComm.isMaster() == true) {
//         this->logger("Direct scheme method is employed\n");
//     }

//     TlDenseSymmetricMatrix_Scalapack F(this->m_nNumOfAOs);
//     if (this->m_bMemorySave == true) {
//         // DfThreeindexintegrals を使わない
//         if (this->m_bIsXCFitting == true) {
//             DfFockMatrix::setXC_RI<TlDenseSymmetricMatrix_Scalapack,
//             TlDenseVector_Lapack,
//             DfOverlap_Parallel>(RUN_RKS, F);
//         } else {
//             DfFockMatrix::setXC_DIRECT<TlDenseSymmetricMatrix_Scalapack>(RUN_RKS,
//             F);
//         }

//         this->setCoulomb(METHOD_RKS, F);
//     } else {
//         // DfThreeindexintegrals を使う
//         assert(this->m_bIsXCFitting == true);
//         F = DfObject::getFpqMatrix<TlDenseSymmetricMatrix_Scalapack>(RUN_RKS,
//         this->m_nIteration);
//     }

//     this->setHpq(RUN_RKS, F);

//     DfObject::saveFpqMatrix(RUN_RKS, this->m_nIteration, F);
// }

void DfFockMatrix_Parallel::mainDIRECT_UKS() {
    TlCommunicate& rComm = TlCommunicate::getInstance();

    if (this->m_bUsingSCALAPACK == true) {
        // ScaLAPACK
        if (rComm.isMaster() == true) {
            this->logger("buid KS matrix using ScaLAPACK.\n");
        }
        std::cout << "do implement!" << std::endl;
        std::abort();
    } else {
        // LAPACK
        if (rComm.isMaster() == true) {
            this->logger("buid KS matrix using LAPACK.\n");
        }
        DfFockMatrix::mainDIRECT_UKS();
    }
}

void DfFockMatrix_Parallel::mainDIRECT_ROKS() {
    TlCommunicate& rComm = TlCommunicate::getInstance();

    if (this->m_bUsingSCALAPACK == true) {
        // ScaLAPACK
        if (rComm.isMaster() == true) {
            this->logger("buid KS matrix using ScaLAPACK.\n");
        }
        std::cout << "do implement!" << std::endl;
        std::abort();
    } else {
        // LAPACK
        if (rComm.isMaster() == true) {
            this->logger("buid KS matrix using LAPACK.\n");
        }
        DfFockMatrix::mainDIRECT_ROKS();
    }
}

void DfFockMatrix_Parallel::setXC_RI(const RUN_TYPE nRunType,
                                     TlDenseSymmetricMatrix_Lapack& F) {
    assert(this->m_bUsingSCALAPACK == false);
    DfFockMatrix::setXC_RI<TlDenseSymmetricMatrix_Lapack, TlDenseVector_Lapack,
                           DfOverlapX_Parallel>(nRunType, F);
}

void DfFockMatrix_Parallel::setXC_DIRECT(const RUN_TYPE nRunType,
                                         TlDenseSymmetricMatrix_Lapack& F) {
    assert(this->m_bUsingSCALAPACK == false);
    TlCommunicate& rComm = TlCommunicate::getInstance();
    if (rComm.isMaster() == true) {
        DfFockMatrix::setXC_DIRECT<TlDenseSymmetricMatrix_Lapack>(nRunType, F);
    }
    rComm.broadcast(&F);
}

void DfFockMatrix_Parallel::setCoulomb(const METHOD_TYPE nMethodType,
                                       TlDenseSymmetricMatrix_Lapack& F) {
    assert(this->m_bUsingSCALAPACK == false);
    TlCommunicate& rComm = TlCommunicate::getInstance();

    TlDenseSymmetricMatrix_Lapack J(this->m_nNumOfAOs);
    if (this->J_engine_ == J_ENGINE_RI_J) {
        DfFockMatrix::setCoulomb<TlDenseSymmetricMatrix_Lapack,
                                 TlDenseVector_Lapack, DfEriX_Parallel>(
            nMethodType, J);
        // if (this->isUseNewEngine_ == true) {
        //     this->logger(" use new engine\n");
        //     DfFockMatrix::setCoulomb<TlDenseSymmetricMatrix_Lapack,
        //     TlDenseVector_Lapack,
        //     DfEriX_Parallel>(nMethodType, J);
        // } else {
        //     DfFockMatrix::setCoulomb<TlDenseSymmetricMatrix_Lapack,
        //     TlDenseVector_Lapack,
        //     DfEri_Parallel>(nMethodType, J);
        // }
        F += J;

        if (rComm.isMaster() == true) {
            // update method
            if (this->m_nIteration > 1) {
                TlDenseSymmetricMatrix_Lapack tmpJ;
                tmpJ = DfObject::getJMatrix<TlDenseSymmetricMatrix_Lapack>(
                    this->m_nIteration - 1);
                J += tmpJ;
            }
            DfObject::saveJMatrix(this->m_nIteration, J);
        }
    } else {
        if (rComm.isMaster() == true) {
            J = this->getJMatrix<TlDenseSymmetricMatrix_Lapack>(
                this->m_nIteration);
            // update method
            if (this->m_nIteration > 1) {
                const TlDenseSymmetricMatrix_Lapack prevJ =
                    DfObject::getJMatrix<TlDenseSymmetricMatrix_Lapack>(
                        this->m_nIteration - 1);
                J -= prevJ;
            }

            F += J;
        }
        rComm.broadcast(&F);
    }
}

void DfFockMatrix_Parallel::setCoulomb(const METHOD_TYPE nMethodType,
                                       TlDenseSymmetricMatrix_Scalapack& F) {
    TlDenseSymmetricMatrix_Scalapack J(this->m_nNumOfAOs);
    if (this->J_engine_ == J_ENGINE_RI_J) {
        DfFockMatrix::setCoulomb<TlDenseSymmetricMatrix_Scalapack,
                                 TlDenseVector_Lapack, DfEriX_Parallel>(
            nMethodType, J);
        // if (this->isUseNewEngine_ == true) {
        //     this->logger(" use new engine\n");
        //     DfFockMatrix::setCoulomb<TlDenseSymmetricMatrix_Scalapack,
        //     TlDenseVector_Lapack,
        //     DfEriX_Parallel>(nMethodType, J);
        // } else {
        //     DfFockMatrix::setCoulomb<TlDenseSymmetricMatrix_Scalapack,
        //     TlDenseVector_Lapack,
        //     DfEri_Parallel>(nMethodType, J);
        // }
        F += J;

        // update method
        if (this->m_nIteration > 1) {
            TlDenseSymmetricMatrix_Scalapack tmpJ;
            tmpJ = DfObject::getJMatrix<TlDenseSymmetricMatrix_Scalapack>(
                this->m_nIteration - 1);
            J += tmpJ;
        }
        DfObject::saveJMatrix(this->m_nIteration, J);
    } else {
        J = this->getJMatrix<TlDenseSymmetricMatrix_Scalapack>(
            this->m_nIteration);

        // update method
        if (this->m_nIteration > 1) {
            const TlDenseSymmetricMatrix_Scalapack prevJ =
                DfObject::getJMatrix<TlDenseSymmetricMatrix_Scalapack>(
                    this->m_nIteration - 1);
            J -= prevJ;
        }

        F += J;
    }
}

TlDenseSymmetricMatrix_Lapack DfFockMatrix_Parallel::getFpqMatrix(
    const RUN_TYPE nRunType, const int nIteration) {
    assert(this->m_bUsingSCALAPACK == false);
    TlDenseSymmetricMatrix_Lapack Fpq;
    TlCommunicate& rComm = TlCommunicate::getInstance();
    if (rComm.isMaster() == true) {
        Fpq = DfObject::getFpqMatrix<TlDenseSymmetricMatrix_Lapack>(nRunType,
                                                                    nIteration);
    }
    rComm.broadcast(&Fpq);
    return Fpq;
}

TlDenseVector_Lapack DfFockMatrix_Parallel::getRho(const RUN_TYPE nRunType,
                                                   const int nIteration) {
    TlCommunicate& rComm = TlCommunicate::getInstance();
    TlDenseVector_Lapack rho;
    if (rComm.isMaster() == true) {
        rho = DfFockMatrix::getRho(nRunType, nIteration);
    }
    rComm.broadcast(&rho);

    return rho;
}

TlDenseVector_Lapack DfFockMatrix_Parallel::getMyu(const RUN_TYPE nRunType,
                                                   const int nIteration) {
    TlCommunicate& rComm = TlCommunicate::getInstance();
    TlDenseVector_Lapack myu;
    if (rComm.isMaster() == true) {
        myu = DfFockMatrix::getMyu(nRunType, nIteration);
    }
    rComm.broadcast(&myu);

    return myu;
}

// void DfFockMatrix_Parallel::saveFpqMatrix(const RUN_TYPE nRunType, const
// TlDenseSymmetricMatrix_Lapack& F)
// {
//     assert(this->m_bUsingSCALAPACK == false);
//     const TlCommunicate& rComm = TlCommunicate::getInstance();
//     if (rComm.isMaster() == true) {
//         DfObject::saveFpqMatrix<TlDenseSymmetricMatrix_Lapack>(nRunType,
//         this->m_nIteration, F);
//     }
// }
