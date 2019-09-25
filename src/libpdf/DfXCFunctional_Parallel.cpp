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

#include "DfXCFunctional_Parallel.h"
#include "CnError.h"
#include "DfCalcGridX_Parallel.h"
#include "DfEriX_Parallel.h"
#include "TlCommunicate.h"
#include "TlMemManager.h"
#include "TlTime.h"
#include "tl_dense_symmetric_matrix_scalapack.h"

DfXCFunctional_Parallel::DfXCFunctional_Parallel(TlSerializeData* pPdfParam)
    : DfXCFunctional(pPdfParam) {}

DfXCFunctional_Parallel::~DfXCFunctional_Parallel() {
    //     std::cerr << "DfXCFunctional_Parallel::~DfXCFunctional_Parallel()" <<
    //     std::endl;
}

void DfXCFunctional_Parallel::logger(const std::string& str) const {
    TlCommunicate& rComm = TlCommunicate::getInstance();

    if (rComm.isMaster() == true) {
        DfXCFunctional::logger(str);
    }
}

void DfXCFunctional_Parallel::buildXcMatrix() {
    // const std::size_t needMem =
    //     this->m_nNumOfAOs * (this->m_nNumOfAOs + 1) * sizeof(double);
    // if ((isWorkOnDisk_ == true) || (this->procMaxMemSize_ < needMem)) {
    //   this->logger(" build XC matrix on disk.\n");
    //   TlMatrix::useMemManager(true);
    // } else {
    //   this->logger(" build XC matrix on memory.\n");
    //   TlMatrix::useMemManager(false);
    // }

    if (this->m_bUsingSCALAPACK == true) {
        this->buildXC_ScaLAPACK();
    } else {
        this->buildXC_LAPACK();
    }
}

void DfXCFunctional_Parallel::buildXC_LAPACK() {
    TlCommunicate& rComm = TlCommunicate::getInstance();
    DfCalcGridX_Parallel dfCalcGrid(this->pPdfParam_);

    switch (this->m_nMethodType) {
        case METHOD_RKS: {
            TlDenseSymmetricMatrix_Lapack Ppq;
            if (rComm.isMaster() == true) {
                if (this->m_bIsUpdateXC == true) {
                    Ppq = 0.5 * this->getDiffDensityMatrix<
                                    TlDenseSymmetricMatrix_Lapack>(
                                    RUN_RKS, this->m_nIteration);
                } else {
                    Ppq =
                        0.5 * this->getPpqMatrix<TlDenseSymmetricMatrix_Lapack>(
                                  RUN_RKS, this->m_nIteration - 1);
                }
            }
            rComm.broadcast(&Ppq);

            TlDenseSymmetricMatrix_Lapack Fxc(this->m_nNumOfAOs);
            DfXCFunctional::getFxc(Ppq, &dfCalcGrid, &Fxc);
            if (this->isSaveFxcPure_ == true) {
                if (rComm.isMaster() == true) {
                    this->saveFxcPureMatrix(RUN_RKS, this->m_nIteration, Fxc);
                }
            }

            // if (this->m_bIsHybrid == true) {
            //     Fxc += this->getFockExchange(this->m_dFockExchangeCoef * Ppq,
            //     RUN_RKS); this->XC_energy_ +=
            //     DfXCFunctional::m_dFockExchangeEnergyAlpha;
            // }

            if (rComm.isMaster() == true) {
                this->saveFxcMatrix<TlDenseSymmetricMatrix_Lapack>(
                    RUN_RKS, this->m_nIteration, Fxc);
            }
        } break;

        case METHOD_UKS: {
            TlDenseSymmetricMatrix_Lapack PApq, PBpq;
            if (rComm.isMaster() == true) {
                if (this->m_bIsUpdateXC == true) {
                    PApq = this->getDiffDensityMatrix<
                        TlDenseSymmetricMatrix_Lapack>(RUN_UKS_ALPHA,
                                                       this->m_nIteration);
                    PBpq = this->getDiffDensityMatrix<
                        TlDenseSymmetricMatrix_Lapack>(RUN_UKS_BETA,
                                                       this->m_nIteration);
                } else {
                    PApq = this->getPpqMatrix<TlDenseSymmetricMatrix_Lapack>(
                        RUN_UKS_ALPHA, this->m_nIteration - 1);
                    PBpq = this->getPpqMatrix<TlDenseSymmetricMatrix_Lapack>(
                        RUN_UKS_BETA, this->m_nIteration - 1);
                }
            }
            rComm.broadcast(&PApq);
            rComm.broadcast(&PBpq);

            TlDenseSymmetricMatrix_Lapack FxcA(this->m_nNumOfAOs);
            TlDenseSymmetricMatrix_Lapack FxcB(this->m_nNumOfAOs);
            DfXCFunctional::getFxc(PApq, PBpq, &dfCalcGrid, &FxcA, &FxcB);
            if (this->isSaveFxcPure_ == true) {
                if (rComm.isMaster() == true) {
                    this->saveFxcPureMatrix(RUN_UKS_ALPHA, this->m_nIteration,
                                            FxcA);
                    this->saveFxcPureMatrix(RUN_UKS_BETA, this->m_nIteration,
                                            FxcB);
                }
            }

            // if (this->m_bIsHybrid == true) {
            //     FxcA += this->getFockExchange(this->m_dFockExchangeCoef *
            //     PApq, RUN_UKS_ALPHA); FxcB +=
            //     this->getFockExchange(this->m_dFockExchangeCoef * PBpq,
            //     RUN_UKS_BETA); this->XC_energy_ +=
            //     (DfXCFunctional::m_dFockExchangeEnergyAlpha +
            //     DfXCFunctional::m_dFockExchangeEnergyBeta);
            // }

            if (rComm.isMaster() == true) {
                this->saveFxcMatrix<TlDenseSymmetricMatrix_Lapack>(
                    RUN_UKS_ALPHA, this->m_nIteration, FxcA);
                this->saveFxcMatrix<TlDenseSymmetricMatrix_Lapack>(
                    RUN_UKS_BETA, this->m_nIteration, FxcB);
            }
        } break;

        case METHOD_ROKS:
            CnErr.abort(
                "sorry. ROKS method is not implemented. @ "
                "DfXCFunctional_Parallel::buildXcMatrix()");
            break;
        default:
            break;
    }
}

void DfXCFunctional_Parallel::buildXC_ScaLAPACK() {
    DfCalcGridX_Parallel dfCalcGrid(this->pPdfParam_);

    switch (this->m_nMethodType) {
        case METHOD_RKS: {
            TlDenseSymmetricMatrix_Scalapack Ppq;
            if (this->m_bIsUpdateXC == true) {
                Ppq = 0.5 * this->getDiffDensityMatrix<
                                TlDenseSymmetricMatrix_Scalapack>(
                                RUN_RKS, this->m_nIteration);
            } else {
                Ppq =
                    0.5 * this->getPpqMatrix<TlDenseSymmetricMatrix_Scalapack>(
                              RUN_RKS, this->m_nIteration - 1);
            }

            TlDenseSymmetricMatrix_Scalapack Fxc(this->m_nNumOfAOs);
            DfXCFunctional::getFxc(Ppq, &dfCalcGrid, &Fxc);
            if (this->isSaveFxcPure_ == true) {
                this->saveFxcPureMatrix(RUN_RKS, this->m_nIteration, Fxc);
            }

            // if (this->m_bIsHybrid == true) {
            //     Fxc += this->getFockExchange(this->m_dFockExchangeCoef * Ppq,
            //     RUN_RKS); this->XC_energy_ +=
            //     DfXCFunctional::m_dFockExchangeEnergyAlpha;
            // }

            this->saveFxcMatrix<TlDenseSymmetricMatrix_Scalapack>(
                RUN_RKS, this->m_nIteration, Fxc);
        } break;

        case METHOD_UKS: {
            TlDenseSymmetricMatrix_Scalapack PApq, PBpq;
            if (this->m_bIsUpdateXC == true) {
                PApq = this->getDiffDensityMatrix<
                    TlDenseSymmetricMatrix_Scalapack>(RUN_UKS_ALPHA,
                                                      this->m_nIteration);
                PBpq = this->getDiffDensityMatrix<
                    TlDenseSymmetricMatrix_Scalapack>(RUN_UKS_BETA,
                                                      this->m_nIteration);
            } else {
                PApq = this->getPpqMatrix<TlDenseSymmetricMatrix_Scalapack>(
                    RUN_UKS_ALPHA, this->m_nIteration - 1);
                PBpq = this->getPpqMatrix<TlDenseSymmetricMatrix_Scalapack>(
                    RUN_UKS_BETA, this->m_nIteration - 1);
            }

            TlDenseSymmetricMatrix_Scalapack FxcA(this->m_nNumOfAOs);
            TlDenseSymmetricMatrix_Scalapack FxcB(this->m_nNumOfAOs);
            DfXCFunctional::getFxc(PApq, PBpq, &dfCalcGrid, &FxcA, &FxcB);
            if (this->isSaveFxcPure_ == true) {
                this->saveFxcPureMatrix(RUN_UKS_ALPHA, this->m_nIteration,
                                        FxcA);
                this->saveFxcPureMatrix(RUN_UKS_BETA, this->m_nIteration, FxcB);
            }

            // if (this->m_bIsHybrid == true) {
            //     FxcA += this->getFockExchange(this->m_dFockExchangeCoef *
            //     PApq, RUN_UKS_ALPHA); FxcB +=
            //     this->getFockExchange(this->m_dFockExchangeCoef * PBpq,
            //     RUN_UKS_BETA); this->XC_energy_ +=
            //     (DfXCFunctional::m_dFockExchangeEnergyAlpha +
            //     DfXCFunctional::m_dFockExchangeEnergyBeta);
            // }

            this->saveFxcMatrix<TlDenseSymmetricMatrix_Scalapack>(
                RUN_UKS_ALPHA, this->m_nIteration, FxcA);
            this->saveFxcMatrix<TlDenseSymmetricMatrix_Scalapack>(
                RUN_UKS_BETA, this->m_nIteration, FxcB);
        } break;

        case METHOD_ROKS:
            CnErr.abort(
                "sorry. ROKS method is not implemented. @ "
                "DfXCFunctional_Parallel::buildXcMatrix()");
            break;
        default:
            break;
    }
}

DfEriX* DfXCFunctional_Parallel::getDfEriXObject() {
    this->logger(" use new engine\n");
    DfEriX* pDfEriX = new DfEriX_Parallel(this->pPdfParam_);

    return pDfEriX;
}

double DfXCFunctional_Parallel::getGrimmeDispersionEnergy() {
    if (this->isCalcd_E_disp_ != true) {
        TlCommunicate& rComm = TlCommunicate::getInstance();
        if (rComm.isMaster() == true) {
            DfXCFunctional::getGrimmeDispersionEnergy();
        }
        rComm.broadcast(this->isCalcd_E_disp_);
        rComm.broadcast(this->E_disp_);
    }

    return this->E_disp_;
}
