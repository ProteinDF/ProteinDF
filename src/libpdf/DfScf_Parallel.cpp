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

#include "DfScf_Parallel.h"

#include "DfCleanup_Parallel.h"
#include "DfConvcheck_Parallel.h"
#include "DfConverge_Anderson_Parallel.h"
#include "DfConverge_Damping_Parallel.h"
#include "DfDensityFittingX_Parallel.h"
#include "DfDensityFittingX_ScaLAPACK.h"
#include "DfDiagonal_Parallel.h"
#include "DfDiffDensityMatrix_Parallel.h"
#include "DfDmatrix_Parallel.h"
#include "DfFockMatrix_Parallel.h"
#include "DfGridFreeXC_Parallel.h"
#include "DfJMatrix_Parallel.h"
#include "DfKMatrix_Parallel.h"
#include "DfPopulation_Parallel.h"
#include "DfSummary_Parallel.h"
#include "DfTotalEnergy_Parallel.h"
#include "DfTransFmatrix_Parallel.h"
#include "DfTransatob_Parallel.h"
#include "DfXCFunctional_Parallel.h"
#include "TlCommunicate.h"

#ifdef HAVE_SCALAPACK
#include "df_total_energy_scalapack.h"
#include "tl_dense_general_matrix_scalapack.h"
#include "tl_dense_symmetric_matrix_scalapack.h"
#include "tl_dense_vector_scalapack.h"
#endif  // HAVE_SCALAPACK

#define NUMBER_OF_CHECK 2

DfScf_Parallel::DfScf_Parallel(TlSerializeData* pPdfParam)
    : DfScf(pPdfParam) {
    TlCommunicate& rComm = TlCommunicate::getInstance();
    DfObject::rank_ = rComm.getRank();
}

DfScf_Parallel::~DfScf_Parallel() {}

void DfScf_Parallel::logger(const std::string& str) const {
    TlCommunicate& rComm = TlCommunicate::getInstance();
    if (rComm.isMaster() == true) {
        DfScf::logger(str);
    }
}

void DfScf_Parallel::saveParam() const {
    TlCommunicate& rComm = TlCommunicate::getInstance();
    if (rComm.isMaster() == true) {
        DfScf::saveParam();
    }
    rComm.broadcast(*(this->pPdfParam_));
}

void DfScf_Parallel::updateParam() {
    TlCommunicate& rComm = TlCommunicate::getInstance();
    if (rComm.isMaster() == true) {
        DfScf::setScfParam();
    }
    rComm.broadcast(*(this->pPdfParam_));
}

// void DfScf_Parallel::setScfParam() {
//     TlCommunicate& rComm = TlCommunicate::getInstance();
//     if (rComm.isMaster() == true) {
//         DfScf::setScfParam();
//     }
//     rComm.broadcast(*(this->pPdfParam_));

//     int dampObject = this->m_nDampObject;
//     rComm.broadcast(dampObject);
//     this->m_nDampObject = (DfScf::DampObjectType)dampObject;

//     int scfAcceleration = this->m_nScfAcceleration;
//     rComm.broadcast(scfAcceleration);
//     this->m_nScfAcceleration = (DfScf::ScfAccelerationType)scfAcceleration;

//     rComm.broadcast(this->m_nConvergenceCounter);
// }

// =====================================================================
// pre-SCF loop
// =====================================================================

// TODO: parallelization
// void DfScf_Parallel::executePreScf()
// {
//     DfPreScf_Parallel dfPreScf(this->pPdfParam_);
//     dfPreScf.prepareGuess();
// }

// =====================================================================
// SCF loop
// =====================================================================
void DfScf_Parallel::diffDensityMatrix() {
    if (this->m_nDampObject == DfScf::DAMP_DENSITY_MATRIX) {
        this->converge();
    }

    if (this->m_bDiskUtilization == false) {
        this->loggerStartTitle("diff density matrix (parallel)");
        DfDiffDensityMatrix_Parallel dfDiffDensityMatrix(this->pPdfParam_);
        dfDiffDensityMatrix.exec();

        this->loggerEndTitle();
    }
}

DfDensityFittingObject* DfScf_Parallel::getDfDensityFittingObject() {
    DfDensityFittingObject* pDfDensityFittingObj = NULL;

#ifdef HAVE_SCALAPACK
    if (this->m_bUsingSCALAPACK == true) {
        pDfDensityFittingObj =
            new DfDensityFittingX_ScaLAPACK(this->pPdfParam_);
        this->logger(" density fitting (parallel; ScaLAPACK)\n");
    } else {
        pDfDensityFittingObj = new DfDensityFittingX_Parallel(this->pPdfParam_);
        this->logger(" density fitting (parallel; LAPACK)");
    }
#else
    {
        pDfDensityFittingObj = new DfDensityFittingX_Parallel(this->pPdfParam_);
        this->logger(" density fitting (parallel; LAPACK)");
    }
#endif  // HAVE_SCALAPACK

    return pDfDensityFittingObj;
}

void DfScf_Parallel::doXCIntegral() {
    TlCommunicate& rComm = TlCommunicate::getInstance();

    if (rComm.isMaster() == true) {
        DfScf::doXCIntegral();
    }
    // rComm.barrier();
}

void DfScf_Parallel::doThreeIndexIntegral() {
    TlCommunicate& rComm = TlCommunicate::getInstance();

    if (rComm.isMaster() == true) {
        DfScf::doThreeIndexIntegral();
    }

    rComm.barrier();
}

DfGridFreeXC* DfScf_Parallel::getDfGridFreeXcObject() {
    DfGridFreeXC* pDfGridFreeXC = new DfGridFreeXC_Parallel(this->pPdfParam_);
    return pDfGridFreeXC;
}

DfXCFunctional* DfScf_Parallel::getDfXCFunctional() {
    DfXCFunctional* pDfXCFunctional =
        new DfXCFunctional_Parallel(this->pPdfParam_);
    return pDfXCFunctional;
}

DfJMatrix* DfScf_Parallel::getDfJMatrixObject() {
    DfJMatrix* pDfJMatrix = new DfJMatrix_Parallel(this->pPdfParam_);
    return pDfJMatrix;
}

DfKMatrix* DfScf_Parallel::getDfKMatrixObject() {
    DfKMatrix* pDfKMatrix = new DfKMatrix_Parallel(this->pPdfParam_);
    return pDfKMatrix;
}

DfFockMatrix* DfScf_Parallel::getDfFockMatrixObject() {
    DfFockMatrix* pDfFockMatrix = new DfFockMatrix_Parallel(this->pPdfParam_);
    return pDfFockMatrix;
}

DfTransFmatrix* DfScf_Parallel::getDfTransFmatrixObject(bool isExecDiis) {
    DfTransFmatrix_Parallel* pDfTransFmatrix =
        new DfTransFmatrix_Parallel(this->pPdfParam_, isExecDiis);
    return pDfTransFmatrix;
}

void DfScf_Parallel::doLevelShift() {
    TlCommunicate& rComm = TlCommunicate::getInstance();

    if (rComm.isMaster() == true) {
        DfScf::doLevelShift();
    }

    rComm.barrier();
}

DfDiagonal* DfScf_Parallel::getDfDiagonalObject() {
    DfDiagonal* pDfDiagonal = new DfDiagonal_Parallel(this->pPdfParam_);
    return pDfDiagonal;
}

DfTransatob* DfScf_Parallel::getDfTransatobObject() {
    DfTransatob* pDfTransAtoB = new DfTransatob_Parallel(this->pPdfParam_);
    return pDfTransAtoB;
}

void DfScf_Parallel::calcDensityMatrix() {
    bool done = false;
#ifdef HAVE_SCALAPACK
    if (this->m_bUsingSCALAPACK == true) {
        TlTime timer;
        this->loggerStartTitle("Density Matirx (ScaLAPACK)");

        DfDmatrix_Parallel dfDmatrix(this->pPdfParam_);
        dfDmatrix.run();

        this->loggerEndTitle();
        (*this->pPdfParam_)["stat"]["elapsed_time"]["density_matrix"]
                           [this->m_nIteration] = timer.getElapseTime();

        this->matrixCache_.flush();
        done = true;
    }
#endif  // HAVE_SCALAPACK

    if (!done) {
        TlCommunicate& rComm = TlCommunicate::getInstance();
        this->log_.info("Density Matirx @master");
        if (rComm.isMaster() == true) {
            DfScf::calcDensityMatrix();
        }
        rComm.barrier();
    }
}

DfDmatrix* DfScf_Parallel::getDfDmatrixObject() {
    DfDmatrix* pDfDmatrix = new DfDmatrix_Parallel(this->pPdfParam_);
    return pDfDmatrix;
}

void DfScf_Parallel::calcTotalEnergy() {
    if (this->linearAlgebraPackage_ == DfObject::LAP_SCALAPACK) {
        this->calcTotalEnergy_tmpl<DfTotalEnergy_Scalapack>();
    } else {
        TlCommunicate& rComm = TlCommunicate::getInstance();
        if (rComm.isMaster()) {
            DfScf::calcTotalEnergy();
        }
    }
}

// DfTotalEnergy* DfScf_Parallel::getDfTotalEnergyObject() {
//     DfTotalEnergy* pDfTotalEnergy =
//         new DfTotalEnergy_Parallel(this->pPdfParam_);
//     return pDfTotalEnergy;
// }

DfPopulation* DfScf_Parallel::getDfPopulationObject() {
    DfPopulation* pDfPopulation = new DfPopulation_Parallel(this->pPdfParam_);
    return pDfPopulation;
}

DfSummary* DfScf_Parallel::getDfSummaryObject() {
    DfSummary* pDfSummary = new DfSummary_Parallel(this->pPdfParam_);
    return pDfSummary;
}

bool DfScf_Parallel::judge() {
    TlCommunicate& rComm = TlCommunicate::getInstance();

    bool bAnswer = false;
    {
        if (rComm.isMaster() == true) {
            // TlLogX& Log = TlLogX::getInstance();
            // Log.outputStartTitle("DfConvcheck");
            this->loggerStartTitle("DfConvcheck(parallel)");
        }

        const bool bJudge = this->checkConverge();

        if (rComm.isMaster() == true) {
            // TlLogX& Log = TlLogX::getInstance();
            // Log.outputEndTitle();
            this->loggerEndTitle();
        }

        //  set convergence
        if (bJudge == true) {
            // 収束の閾値をみたすとcounterが１増える
            this->m_nConvergenceCounter++;

            if (rComm.isMaster() == true) {
                if (this->m_nConvergenceCounter == 1) {
                    std::cout << " *** Convergence conditions are satisfied: "
                                 "(1st) ***\n";
                } else if (this->m_nConvergenceCounter == 2) {
                    std::cout << " *** Convergence conditions are satisfied: "
                                 "(2nd) ***\n";
                    std::cout << " *** SCF is well converged ***" << std::endl;
                }
            }

            // conv_counterとNUMBER_OF_CHECKが一致したら収束とみなす
            if (this->m_nConvergenceCounter == NUMBER_OF_CHECK) {
                bAnswer = true;
            }
        } else {
            // conv_counterとNUMBER_OF_CHECKが一致する前に
            // 収束の閾値を満たさないことがあれば、
            // conv_counterは0に戻る
            this->m_nConvergenceCounter = 0;
        }
    }

    return bAnswer;
}

bool DfScf_Parallel::checkConverge() {
    DfConvcheck_Parallel dfConvcheck(this->pPdfParam_, this->m_nIteration);

    return dfConvcheck.isConverged();
}

DfConverge* DfScf_Parallel::getDfConverge() {
    DfConverge* pDfConverge = NULL;
    if (this->m_nScfAcceleration == SCF_ACCELERATION_SIMPLE) {
        pDfConverge = new DfConverge_Damping_Parallel(this->pPdfParam_);
    } else if (this->m_nScfAcceleration == SCF_ACCELERATION_ANDERSON) {
        pDfConverge = new DfConverge_Anderson_Parallel(this->pPdfParam_);
    } else {
        // diis 法の最初のdampingなど
        pDfConverge = new DfConverge_Damping_Parallel(this->pPdfParam_);
    }
    return pDfConverge;
}

void DfScf_Parallel::cleanup() {
    //     if (TlUtils::toUpper((*this->pPdfParam_)["cleanup"].getStr()) !=
    //     "NO")
    //     {
    //         this->loggerStartTitle("cleanup files");
    //         DfCleanup_Parallel dfCleanup(this->pPdfParam_);
    //         dfCleanup.cleanup();
    //         this->loggerEndTitle();
    //     }
}

bool DfScf_Parallel::checkMaxIteration() {
    TlCommunicate& rComm = TlCommunicate::getInstance();

    bool answer = false;
    if (rComm.isMaster() == true) {
        answer = DfScf::checkMaxIteration();
    }
    rComm.broadcast(answer);

    return answer;
}
