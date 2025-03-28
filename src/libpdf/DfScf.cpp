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

#ifdef __FUJITSU
#include <fjcoll.h>
#endif  // __FUJITSU

#include <cmath>

#include "DfCalcGrid.h"
#include "DfCleanup.h"
#include "DfConvcheck.h"
// #include "DfConverge_Anderson.h"
// #include "DfConverge_DIIS.h"
// #include "DfConverge_Damping.h"
#include "DfDensityFittingX.h"
#include "DfDiagonal.h"
#include "DfDiffDensityMatrix.h"
#include "DfDmatrix.h"
#include "DfFockMatrix.h"
#include "DfGridFreeXC.h"
#include "DfJMatrix.h"
#include "DfKMatrix.h"
#include "DfLevelshift.h"
#include "DfPopulation.h"
#include "DfScf.h"
#include "DfSummary.h"
#include "DfThreeindexintegrals.h"
#include "TlUtils.h"
#include "common.h"
#include "df_converge_damping_anderson_lapack.h"
#include "df_converge_damping_diis_lapack.h"
#include "df_converge_damping_lapack.h"
#include "df_converge_damping_oda_lapack.h"
#include "tl_matrix_utils.h"
// #include "DfTotalEnergy.h"

#ifdef HAVE_LAPACK
#include "df_total_energy_lapack.h"
#endif  // HAVE_LAPACK

#ifdef HAVE_EIGEN
#include "df_total_energy_eigen.h"
#endif  // HAVE_EIGEN

#include "DfTransFmatrix.h"
#include "DfTransatob.h"
#include "DfXCFunctional.h"

// FoR Extended QCLO
#include "DfCqclomatrix.h"
#include "DfPreScf.h"
#include "DfQclo.h"
#include "TlFile.h"
#include "TlMsgPack.h"
#include "tl_dense_general_matrix_eigen.h"
#include "tl_dense_general_matrix_lapack.h"
#include "tl_dense_symmetric_matrix_eigen.h"
#include "tl_dense_symmetric_matrix_lapack.h"
#include "tl_dense_vector_eigen.h"
#include "tl_dense_vector_lapack.h"

#define NUMBER_OF_CHECK 2

DfScf::DfScf(TlSerializeData* pPdfParam)
    : DfObject(pPdfParam), m_nDampObject(DAMP_NONE) {
    this->isUseNewEngine_ = (*pPdfParam)["new_engine"].getBoolean();
}

DfScf::~DfScf() {
}

void DfScf::saveParam() const {
    (*(this->pPdfParam_))["num_of_iterations"] = this->m_nIteration;

    // save PDF parameter
    const std::string pdfParamPath = (*this->pPdfParam_)["pdf_param_path"].getStr();
    TlMsgPack pdfParam_mpac(*(this->pPdfParam_));
    pdfParam_mpac.save(pdfParamPath);
}

// return  0 : not convergence
//         1 : convergence
int DfScf::run() {
    this->updateParam();
    this->setScfParam();
    const TlSerializeData& pdfParam = *(this->pPdfParam_);

    std::string sStepControl = pdfParam["step_control"].getStr();
    std::string group = "";
    do {
        group = TlUtils::toUpper(TlUtils::getWord(sStepControl));
    } while ((group != "SCF") && (group != "SCFQCLO"));

    // QCLO法用。PreSCFに持って行くべき
    if (group == "SCFQCLO") {
        if ((this->isRestart_ == false) || (pdfParam["DfScf"]["scf-restart-point"] == "startPDF")) {
            this->loggerStartTitle("DfCqclomatrix");

            DfCqclomatrix dfcqclomatrix(this->pPdfParam_);
            dfcqclomatrix.main();

            this->loggerEndTitle();
        }
    }

    // start SCF LOOP
    this->saveParam();
    return this->execScfLoop();
}

void DfScf::updateParam() {
    // numOfMOs
    DfObject::index_type X_cols = 0;
    {
        TlMatrixObject::HeaderInfo headerInfo;
        const std::string XmatPath = this->getXMatrixPath();
        const bool isLoadable = TlMatrixUtils::getHeaderInfo(XmatPath, &headerInfo);
        if (isLoadable) {
            X_cols = headerInfo.numOfCols;
        } else {
            CnErr.abort("cannot load X matrix.");
        }
    }

    if (X_cols != this->m_nNumOfMOs) {
        this->log_.warn("the number of MOs is not equal to the number of columns in the X matrix.");
        this->log_.warn(TlUtils::format("  [#MOs=%d] != [X_cols=%d]", this->m_nNumOfMOs, X_cols));

        this->log_.warn(TlUtils::format("force update of the number of MOs to %d.", X_cols));
        this->m_nNumOfMOs = X_cols;
        (*(this->pPdfParam_))["num_of_MOs"] = this->m_nNumOfMOs;
    }
}

void DfScf::setScfParam() {
    DfObject::setParam(*(this->pPdfParam_));

    const TlSerializeData& pdfParam = *(this->pPdfParam_);

    // iteration number
    this->m_nIteration = 1;
    if (this->isRestart_) {
        this->m_nIteration = std::max<int>(1, pdfParam["num_of_iterations"].getInt());
    } else {
        (*(this->pPdfParam_))["num_of_iterations"] = this->m_nIteration;
    }

    // damping switch
    this->m_nScfAcceleration = SCF_ACCELERATION_SIMPLE;
    {
        const std::string sScfAcceleration = TlUtils::toUpper(pdfParam["scf_acceleration"].getStr());

        if (sScfAcceleration == "DAMPING") {
            this->m_nScfAcceleration = SCF_ACCELERATION_SIMPLE;
        } else if (sScfAcceleration == "ODA") {
            this->m_nScfAcceleration = SCF_ACCELERATION_ODA;
        } else if (sScfAcceleration == "ANDERSON") {
            this->m_nScfAcceleration = SCF_ACCELERATION_ANDERSON;
        } else if (sScfAcceleration == "DIIS") {
            this->m_nScfAcceleration = SCF_ACCELERATION_DIIS;
        } else {
            this->log_.warn(TlUtils::format("unknown acceleration method: %s", sScfAcceleration.c_str()));
        }
    }

    // damping object
    this->m_nDampObject = DAMP_DENSITY;
    {
        const std::string sDampObject = TlUtils::toUpper(pdfParam["scf_acceleration/damping/damping_type"].getStr());

        if (sDampObject == "DENSITY_MATRIX") {
            this->m_nDampObject = DAMP_DENSITY_MATRIX;
        } else if (sDampObject == "DENSITY") {
            this->m_nDampObject = DAMP_DENSITY;
        } else if (sDampObject == "FOCK") {
            this->m_nDampObject = DAMP_FOCK;
        } else {
            this->log_.warn(TlUtils::format("unknown acceleration method: %s", sDampObject.c_str()));
        }
    }
    if (this->m_nDampObject == DAMP_DENSITY) {
        if (this->J_engine_ != J_ENGINE_RI_J) {
            this->log_.warn("damping \"density\" is not supported except for the RI method.");
            this->log_.warn("\"density_matrix\" is used.");
            (*this->pPdfParam_)["scf_acceleration/damping/damping_type"] = "density_matrix";
            this->m_nDampObject = DAMP_DENSITY_MATRIX;
        }
    }
}

int DfScf::execScfLoop() {
    const TlSerializeData& pdfParam = *(this->pPdfParam_);
    // int OUTSCF_FLAG = 1;
    std::string sStepControl = pdfParam["step_control"].getStr();

    std::string group = "";
    do {
        group = TlUtils::getWord(sStepControl);
    } while (group != "scf" && group != "scfqclo");

    // prepare restart
    enum SCF_STATE {
        UNDEFINED,
        DIFF_DENSITY_MATRIX,
        DENSITY_FITTING,
        XC_INTEGRAL,
        BEGIN,
        THREE_INDEX_INTEGRAL,
        XC_MATRIX,
        K_MATRIX,
        J_MATRIX,
        FOCK,
        ENDFOCK,
        TRANSFORM_FOCK,
        LEVEL_SHIFT,
        DIAGONAL,
        ENDFOCK_TRANSC,
        DENSITY_MATRIX,
        TOTAL_ENERGY,
        POPULATION,
        SUMMARY,
        JUDGE,
        JUDGE_TAIL,
        END_OF_SCF_LOOP,
        END_SCF_LOOP
    };

    SCF_STATE nScfState = BEGIN;
    (*this->pPdfParam_)["control"]["scf_converged"] = false;
    if (this->isRestart_ == true) {
        const std::string respoint = TlUtils::toUpper(pdfParam["control"]["scf_state"].getStr());

        if (respoint == "prepareGuess") {
            nScfState = BEGIN;
        } else if (respoint == "DIFF_DENSITY_MATRIX") {
            nScfState = DENSITY_FITTING;
        } else if ("DENSITY_FITTING" == respoint) {
            nScfState = XC_INTEGRAL;
        } else if ("XC_INTEGRAL" == respoint) {
            nScfState = THREE_INDEX_INTEGRAL;
        } else if ("THREE_INDEX_INTEGRAL" == respoint) {
            nScfState = XC_MATRIX;
        } else if ("XC_MATRIX" == respoint) {
            nScfState = K_MATRIX;
        } else if ("K_MATRIX" == respoint) {
            nScfState = J_MATRIX;
        } else if ("J_MATRIX" == respoint) {
            nScfState = FOCK;
        } else if ("FOCK" == respoint) {
            nScfState = ENDFOCK;
        } else if ("TRANSFORM_FOCK" == respoint) {
            nScfState = LEVEL_SHIFT;
        } else if ("LEVEL_SHIFT" == respoint) {
            nScfState = DIAGONAL;
        } else if ("DIAGONAL" == respoint) {
            nScfState = ENDFOCK_TRANSC;
        } else if ("TRANSC" == respoint) {
            nScfState = DENSITY_MATRIX;
        } else if ("QCLO" == respoint) {
            nScfState = DENSITY_MATRIX;
        } else if ("DENSITY_MATRIX" == respoint) {
            nScfState = TOTAL_ENERGY;
        } else if ("TOTAL_ENERGY" == respoint) {
            nScfState = POPULATION;
        } else if ("POPULATION" == respoint) {
            nScfState = SUMMARY;
        } else if ("SUMMARY" == respoint) {
            nScfState = JUDGE;
        } else {
            this->logger(
                " restart calculation is indicated, but state is not "
                "defined.\n");
            nScfState = BEGIN;
        }
    }

    while (nScfState != END_SCF_LOOP) {
        TlTime timer(true);

        switch (nScfState) {
            case BEGIN:
                this->logger("// =====================================\n");
                this->logger(TlUtils::format("number of iteration = %d\n", this->m_nIteration));
                this->logger("// =====================================\n");
                nScfState = DIFF_DENSITY_MATRIX;
                break;

            case DIFF_DENSITY_MATRIX:
                // this->log_.info("diff density matrix");
                this->diffDensityMatrix();
                // this->log_.info("diff density matrix: done");
                // this->setScfRestartPoint("DIFF_DENSITY_MATRIX");
                // this->log_.info("diff density matrix: done2");
                nScfState = DENSITY_FITTING;
                break;

            case DENSITY_FITTING:
                this->doDensityFitting();
                this->setScfRestartPoint("DENSITY_FITTING");
                nScfState = XC_INTEGRAL;
                break;

            case XC_INTEGRAL:
                this->doXCIntegral();
                this->setScfRestartPoint("XC_INTEGRAL");
                nScfState = THREE_INDEX_INTEGRAL;
                break;

            case THREE_INDEX_INTEGRAL:
                this->doThreeIndexIntegral();
                this->setScfRestartPoint("THREE_INDEX_INTEGRAL");
                nScfState = XC_MATRIX;
                break;

            case XC_MATRIX:
                this->buildXcMatrix();
                this->setScfRestartPoint("XC_MATRIX");
                nScfState = K_MATRIX;
                break;

            case K_MATRIX:
                this->buildKMatrix();
                this->setScfRestartPoint("K_MATRIX");
                nScfState = J_MATRIX;
                break;

            case J_MATRIX:
                this->buildJMatrix();
                this->setScfRestartPoint("J_MATRIX");
                nScfState = FOCK;
                break;

            case FOCK:
                this->buildFock();
                this->setScfRestartPoint("FOCK");
                nScfState = ENDFOCK;
                break;

            case ENDFOCK:
                if (group == "scf") {
                    nScfState = TRANSFORM_FOCK;
                } else if (group == "scfqclo") {
                    this->loggerStartTitle("DfQclo");

                    DfQclo dfQclo(this->pPdfParam_, this->m_nIteration, false);
                    dfQclo.DfQcloMain();

                    this->loggerEndTitle();

                    this->setScfRestartPoint("QCLO");
                    nScfState = DENSITY_MATRIX;
                }
                break;

            case TRANSFORM_FOCK:
                assert(group == "scf");
                this->transformFock();
                this->setScfRestartPoint("TRANSFORM_FOCK");
                nScfState = LEVEL_SHIFT;
                break;

            case LEVEL_SHIFT:
                assert(group == "scf");
                this->doLevelShift();
                this->setScfRestartPoint("LEVEL_SHIFT");
                nScfState = DIAGONAL;
                break;

            case DIAGONAL:
                assert(group == "scf");
                this->diagonal();
                this->setScfRestartPoint("DIAGONAL");
                nScfState = ENDFOCK_TRANSC;
                break;

            case ENDFOCK_TRANSC:
                assert(group == "scf");
                this->execScfLoop_EndFock_TransC();
                this->setScfRestartPoint("TRANSC");
                nScfState = DENSITY_MATRIX;
                break;

            case DENSITY_MATRIX:
                this->calcDensityMatrix();
                this->setScfRestartPoint("DENSITY_MATRIX");
                nScfState = TOTAL_ENERGY;
                break;

            case TOTAL_ENERGY:
                this->calcTotalEnergy();
                this->setScfRestartPoint("TOTAL_ENERGY");
                nScfState = POPULATION;
                break;

            case POPULATION:
                if (TlUtils::toUpper((*(this->pPdfParam_))["analyze_population"].getStr()) == "EVERY-SCF") {
                    this->calcPopulation();
                }
                this->setScfRestartPoint("POPULATION");
                nScfState = SUMMARY;
                break;

            case SUMMARY:
                if (TlUtils::toUpper(pdfParam["summary"].getStr()) == "EVERY-SCF") {
                    this->summarize();
                }
                this->setScfRestartPoint("SUMMARY");
                nScfState = JUDGE;
                break;

            case JUDGE:
                if (this->judge() == false) {
                    nScfState = JUDGE_TAIL;
                } else {
                    (*this->pPdfParam_)["control"]["scf_converged"] = true;
                    nScfState = END_OF_SCF_LOOP;
                }

                break;

            case JUDGE_TAIL: {
                this->cleanup();
                const bool isExitScf = this->checkMaxIteration();

                if (isExitScf == false) {
                    ++(this->m_nIteration);
                    nScfState = BEGIN;
                } else {
                    nScfState = END_OF_SCF_LOOP;
                }
            } break;

            case END_OF_SCF_LOOP:
                nScfState = END_SCF_LOOP;
                break;

            default:
                std::cerr << "unknown ScfState(" << (int)nScfState << "). stop." << std::endl;
                exit(1);
                break;
        }

        (*this->pPdfParam_)["stat"]["elapsed_time"]["scf"][this->m_nIteration] = timer.getElapseTime();
        this->saveParam();
    }

    // pupulation analysis
    const std::string sAnalizePopulation = TlUtils::toUpper(pdfParam["analyze_population"].getStr());
    if ((sAnalizePopulation == "EVERY-SCF") || (sAnalizePopulation == "CONVERGENCE")) {
        this->calcPopulation();
    }

    // DfSummary
    const std::string sDfSummary = TlUtils::toUpper(pdfParam["summary"].getStr());
    if ((sDfSummary == "EVERY-SCF") || (sDfSummary == "CONVERGENCE")) {
        this->summarize();
    }

    return 0;
}

void DfScf::setScfRestartPoint(const std::string& str) {
    (*(this->pPdfParam_))["control"]["scf_state"] = str;
    this->saveParam();
}

void DfScf::diffDensityMatrix() {
    if (this->m_nDampObject == DAMP_DENSITY_MATRIX) {
        this->converge();
    }

    if (this->m_bDiskUtilization == false) {
        TlTime timer(true);
        this->loggerStartTitle("diff density matrix");
        DfDiffDensityMatrix dddm(this->pPdfParam_);
        dddm.exec();

        this->loggerEndTitle();
        (*this->pPdfParam_)["stat"]["elapsed_time"]["diff_density_matrix"][this->m_nIteration] = timer.getElapseTime();
    }
}

void DfScf::doDensityFitting() {
    // #ifdef __FUJITSU
    //     start_collection("density_fitting");
    // #endif  // __FUJITSU

    if (this->J_engine_ == J_ENGINE_RI_J) {
        if (this->m_nIteration == 1) {
            if ((this->initialGuessType_ == GUESS_RHO) || (this->initialGuessType_ == GUESS_FILE_RHO)) {
                this->logger(
                    "Initial rho is provided. Density fitting is "
                    "unnecessary.\n");
                return;
            }
        }

        TlTime timer(true);
        this->loggerStartTitle("Density Fitting");

        DfDensityFittingObject* pDfDensityFitting = this->getDfDensityFittingObject();
        pDfDensityFitting->exec();
        delete pDfDensityFitting;
        pDfDensityFitting = NULL;

        this->loggerEndTitle();
        (*this->pPdfParam_)["stat"]["elapsed_time"]["density_fitting"][this->m_nIteration] = timer.getElapseTime();

        if (this->m_nDampObject == DAMP_DENSITY) {
            this->converge();
        }
    }

    // #ifdef __FUJITSU
    //     stop_collection("density_fitting");
    // #endif  // __FUJITSU
}

DfDensityFittingObject* DfScf::getDfDensityFittingObject() {
    DfDensityFittingObject* pDfDensityFittingObj = NULL;
    pDfDensityFittingObj = new DfDensityFittingX(this->pPdfParam_);

    return pDfDensityFittingObj;
}

void DfScf::doXCIntegral() {
    if ((this->isDFT_ == true) && (this->m_bIsXCFitting == true)) {
        this->loggerStartTitle("Grid fitting Myu");

        DfCalcGrid dg(this->pPdfParam_, this->m_nIteration);
        dg.dfGrdMain();

        this->loggerEndTitle();
    }
}

void DfScf::doThreeIndexIntegral() {
    if ((this->isDFT_ == true) && (this->m_bMemorySave != true) && (this->m_bIsXCFitting == true)) {
        this->loggerStartTitle("DfThreeindexintegrals");
        DfThreeindexintegrals dfThreeindexintegrals(this->pPdfParam_);
        dfThreeindexintegrals.DfThreeindexintegralsMain();
        this->loggerEndTitle();
    }
}

void DfScf::buildXcMatrix() {
#ifdef __FUJITSU
    fapp_start("XC_matrix", 1, 0);
#endif  // __FUJITSU

    if ((this->isDFT_ == true) && (this->m_bIsXCFitting == false)) {
        TlTime timer(true);
        this->loggerStartTitle("generate XC matrix");

        if (this->XC_engine_ != XC_ENGINE_GRID) {
            this->log_.info("using grid-free method");
            DfGridFreeXC* pDfGridFreeXC = this->getDfGridFreeXcObject();
            pDfGridFreeXC->buildFxc();

            delete pDfGridFreeXC;
            pDfGridFreeXC = NULL;
        } else {
            this->log_.info("using grid method");

            // for restart
            if (this->isRestart_ == true) {
                const std::string prevGridDataFilePath =
                    TlUtils::format("%s.itr%d", this->getGridDataFilePath().c_str(), this->m_nIteration - 1);
                TlFile::copy(prevGridDataFilePath, this->getGridDataFilePath());
            }

            DfXCFunctional* pDfXCFunctional = this->getDfXCFunctional();
            pDfXCFunctional->buildXcMatrix();

            delete pDfXCFunctional;
            pDfXCFunctional = NULL;
        }

        this->loggerEndTitle();
        (*this->pPdfParam_)["stat"]["elapsed_time"]["xc_matrix"][this->m_nIteration] = timer.getElapseTime();

        // flush
        this->matrixCache_.flush();
    }

#ifdef __FUJITSU
    fapp_stop("XC_matrix", 1, 0);
#endif  // __FUJITSU
}

DfGridFreeXC* DfScf::getDfGridFreeXcObject() {
    DfGridFreeXC* pDfGridFreeXC = new DfGridFreeXC(this->pPdfParam_);
    return pDfGridFreeXC;
}

DfXCFunctional* DfScf::getDfXCFunctional() {
    DfXCFunctional* pDfXCFunctional = new DfXCFunctional(this->pPdfParam_);
    return pDfXCFunctional;
}

void DfScf::buildJMatrix() {
#ifdef __FUJITSU
    fapp_start("J_matrix", 1, 0);
#endif  // __FUJITSU

    TlTime timer(true);

    this->loggerStartTitle("J matrix");
    DfJMatrix* pDfJMatrix = this->getDfJMatrixObject();
    pDfJMatrix->buildJ();
    this->loggerEndTitle();

    (*this->pPdfParam_)["stat"]["elapsed_time"]["j_matrix"][this->m_nIteration] = timer.getElapseTime();

#ifdef __FUJITSU
    fapp_stop("J_matrix", 1, 0);
#endif  // __FUJITSU
}

DfJMatrix* DfScf::getDfJMatrixObject() {
    DfJMatrix* pDfJMatrix = new DfJMatrix(this->pPdfParam_);
    return pDfJMatrix;
}

void DfScf::buildKMatrix() {
#ifdef __FUJITSU
    fapp_start("K_matrix", 1, 0);
#endif  // __FUJITSU

    const DfXCFunctional dfXCFunctional(this->pPdfParam_);
    if (dfXCFunctional.isHybridFunctional() == true) {
        TlTime timer(true);

        this->loggerStartTitle("K matrix");
        DfKMatrix* pDfKMatrix = this->getDfKMatrixObject();
        pDfKMatrix->buildK();
        this->loggerEndTitle();

        (*this->pPdfParam_)["stat"]["elapsed_time"]["k_matrix"][this->m_nIteration] = timer.getElapseTime();
    }

#ifdef __FUJITSU
    fapp_stop("K_matrix", 1, 0);
#endif  // __FUJITSU
}

DfKMatrix* DfScf::getDfKMatrixObject() {
    DfKMatrix* pDfKMatrix = new DfKMatrix(this->pPdfParam_);
    return pDfKMatrix;
}

void DfScf::buildFock() {
#ifdef __FUJITSU
    fapp_start("Fock", 1, 0);
#endif  // __FUJITSU

    TlTime timer(true);

    this->loggerStartTitle("Fock matrix");
    DfFockMatrix* pDfFockMatrix = this->getDfFockMatrixObject();
    pDfFockMatrix->DfFockMatrixMain();
    this->loggerEndTitle();

    (*this->pPdfParam_)["stat"]["elapsed_time"]["fock_matrix"][this->m_nIteration] = timer.getElapseTime();

    if (this->m_nDampObject == DAMP_FOCK) {
        this->converge();
    }

    delete pDfFockMatrix;
    pDfFockMatrix = NULL;

    // flush
    this->matrixCache_.flush();

#ifdef __FUJITSU
    fapp_stop("Fock", 1, 0);
#endif  // __FUJITSU
}

DfFockMatrix* DfScf::getDfFockMatrixObject() {
    DfFockMatrix* pDfFockMatrix = new DfFockMatrix(this->pPdfParam_);
    return pDfFockMatrix;
}

void DfScf::transformFock() {
#ifdef __FUJITSU
    fapp_start("Transform_matrix", 1, 0);
#endif  // __FUJITSU

    // transformed to orth. A.O. based Fock matrix
    TlTime timer(true);
    this->loggerStartTitle("Transform KS matrix");

    DfTransFmatrix* pDfTransFmatrix = this->getDfTransFmatrixObject(false);
    pDfTransFmatrix->DfTrsFmatMain();
    delete pDfTransFmatrix;
    pDfTransFmatrix = NULL;

    this->loggerEndTitle();
    (*this->pPdfParam_)["stat"]["elapsed_time"]["transform_F_matrix"][this->m_nIteration] = timer.getElapseTime();

#ifdef __FUJITSU
    fapp_stop("Transform_matrix", 1, 0);
#endif  // __FUJITSU
}

DfTransFmatrix* DfScf::getDfTransFmatrixObject(bool isExecDiis) {
    DfTransFmatrix* pDfTransFmatrix = new DfTransFmatrix(this->pPdfParam_, isExecDiis);
    return pDfTransFmatrix;
}

void DfScf::doLevelShift() {
    // add level shift to Kohn-Sham matrix
    const int start_iter = (*(this->pPdfParam_))["level_shift/start_iteration"].getInt();
    const bool levelShift = (*(this->pPdfParam_))["level_shift"].getBoolean();
    if ((levelShift == true) && (this->m_nIteration >= start_iter)) {
        TlTime timer(true);
        this->loggerStartTitle("Level shift");

        DfLevelshift LS(this->pPdfParam_, this->m_nIteration);
        LS.DfLshiftMain();

        this->loggerEndTitle();
        (*this->pPdfParam_)["stat"]["elapsed_time"]["level_shift"][this->m_nIteration] = timer.getElapseTime();
    }
}

void DfScf::diagonal() {
#ifdef __FUJITSU
    fapp_start("Diagonal", 1, 0);
#endif  // __FUJITSU

    // Diagonarize Fock matrix
    TlTime timer(true);
    this->loggerStartTitle("Diagonal");

    DfDiagonal* pDfDiagonal = this->getDfDiagonalObject();
    pDfDiagonal->run();
    delete pDfDiagonal;
    pDfDiagonal = NULL;

    // switch (this->linearAlgebraPackage_) {
    //   case LAP_LAPACK: {
    //     this->log_.info("Linear Algebra Package: LAPACK");
    //     DfDiagonalTempl<TlDenseGeneralMatrix_Lapack,
    //                     TlDenseSymmetricMatrix_Lapack, TlDenseVector_Lapack>
    //         dfDiagonal(this->pPdfParam_);
    //     dfDiagonal.run();
    //   } break;

    //   case LAP_EIGEN:
    //   case LAP_VIENNACL:
    //   {
    //     this->log_.info("Linear Algebra Package: Eigen");
    //     DfDiagonalTempl<TlDenseGeneralMatrix_Eigen,
    //     TlDenseSymmetricMatrix_Eigen,
    //                     TlDenseVector_Eigen>
    //         dfDiagonal(this->pPdfParam_);
    //     dfDiagonal.run();
    //   } break;

    //   // case LAP_VIENNACL: {
    //   //   this->log_.info("Linear Algebra Package: ViennaCL");
    //   //   DfDiagonalTempl<TlDenseGeneralMatrix_ViennaCL,
    //   //                   TlDenseSymmetricMatrix_ViennaCL,
    //   TlDenseVector_ViennaCL>
    //   //       dfDiagonal(this->pPdfParam_);
    //   //   dfDiagonal.run();
    //   // } break;

    //   default:
    //     CnErr.abort(TlUtils::format("program error: @%s,%d", __FILE__,
    //     __LINE__));
    // }

    this->loggerEndTitle();
    (*this->pPdfParam_)["stat"]["elapsed_time"]["diagonal"][this->m_nIteration] = timer.getElapseTime();

    // flush
    this->matrixCache_.flush();

#ifdef __FUJITSU
    fapp_stop("Diagonal", 1, 0);
#endif  // __FUJITSU
}

DfDiagonal* DfScf::getDfDiagonalObject() {
    DfDiagonal* pDfDiagonal = new DfDiagonal(this->pPdfParam_);
    return pDfDiagonal;
}

void DfScf::execScfLoop_EndFock_TransC() {
#ifdef __FUJITSU
    fapp_start("Transform_C", 1, 0);
#endif  // __FUJITSU

    // transformed to original nonorth. A.O.based space
    TlTime timer(true);

    this->loggerStartTitle("Transform Matrix");
    DfTransatob* pDfTransAtoB = this->getDfTransatobObject();
    pDfTransAtoB->run();
    delete pDfTransAtoB;
    pDfTransAtoB = NULL;
    this->loggerEndTitle();

    (*this->pPdfParam_)["stat"]["elapsed_time"]["transform_C_matrix"][this->m_nIteration] = timer.getElapseTime();

#ifdef __FUJITSU
    fapp_stop("Transform_C", 1, 0);
#endif  // __FUJITSU
}

DfTransatob* DfScf::getDfTransatobObject() {
    DfTransatob* pDfTransAtoB = new DfTransatob(this->pPdfParam_);
    return pDfTransAtoB;
}

void DfScf::calcDensityMatrix() {
#ifdef __FUJITSU
    fapp_start("P_matrix", 1, 0);
#endif  // __FUJITSU

    // density matrix generation
    TlTime timer(true);
    this->loggerStartTitle("Density Matirx");

    DfDmatrix dfDmatrix(this->pPdfParam_);
    dfDmatrix.run();

    this->loggerEndTitle();
    (*this->pPdfParam_)["stat"]["elapsed_time"]["density_matrix"][this->m_nIteration] = timer.getElapseTime();

    // flush
    this->matrixCache_.flush();

#ifdef __FUJITSU
    fapp_stop("P_matrix", 1, 0);
#endif  // __FUJITSU
}

DfDmatrix* DfScf::getDfDmatrixObject() {
    DfDmatrix* pDfDmatrix = new DfDmatrix(this->pPdfParam_);
    return pDfDmatrix;
}

// void DfScf::calcTotalEnergy() {
//     // calculate total energy
//     TlTime timer(true);
//     this->loggerStartTitle("Total Energy");
//     DfTotalEnergy* pDfTotalEnergy = this->getDfTotalEnergyObject();
//     pDfTotalEnergy->exec();
//     delete pDfTotalEnergy;
//     pDfTotalEnergy = NULL;
//     this->loggerEndTitle();
//     (*this->pPdfParam_)["stat"]["elapsed_time"]["total_energy"]
//                        [this->m_nIteration] = timer.getElapseTime();
// }

// DfTotalEnergy* DfScf::getDfTotalEnergyObject() {
//     DfTotalEnergy* pDfTotalEnergy = new DfTotalEnergy(this->pPdfParam_);
//     return pDfTotalEnergy;
// }

// calculate total energy
void DfScf::calcTotalEnergy() {
    this->calcTotalEnergy_tmpl<DfTotalEnergy_Lapack>();
}

DfPopulation* DfScf::getDfPopulationObject() {
    DfPopulation* pDfPopulation = new DfPopulation(this->pPdfParam_);
    return pDfPopulation;
}

void DfScf::calcPopulation() {
    // pupulation analysis
    this->loggerStartTitle("Population analysis");

    DfPopulation* pDfPopulation = this->getDfPopulationObject();
    {
        std::stringstream ss;
        pDfPopulation->getReport(this->m_nIteration, ss);
        this->log_.info(ss.str());
    }

    delete pDfPopulation;
    pDfPopulation = NULL;
    this->loggerEndTitle();
}

DfSummary* DfScf::getDfSummaryObject() {
    DfSummary* pDfSummary = new DfSummary(this->pPdfParam_);
    return pDfSummary;
}

void DfScf::summarize() {
    // std::cout << "DfScf::summarize() called." << std::endl;
    this->loggerStartTitle("Summary");
    DfSummary* pDfSummary = this->getDfSummaryObject();
    pDfSummary->exec();
    delete pDfSummary;
    pDfSummary = NULL;

    this->loggerEndTitle();
    // std::cout << "DfScf::summarize() exit." << std::endl;
}

bool DfScf::judge() {
    // std::cerr << "enter JUDGE" << std::endl;
    // OUTSCF_FLAG=1;

    bool bAnswer = false;
    {
        this->loggerStartTitle("Convergence Check");

        const bool bJudge = this->checkConverge();

        this->loggerEndTitle();

        //  set convergence
        if (bJudge == true) {
            // 収束の閾値をみたすとconv_counterが１増える
            // conv_counter++;
            this->m_nConvergenceCounter++;

            // if (conv_counter == 1){
            if (this->m_nConvergenceCounter == 1) {
                const std::string str = " *** Convergence conditions are satisfied: (1st) ***\n";
                this->logger(str);
                std::cout << str;
                //} else if (conv_counter == 2){
            } else if (this->m_nConvergenceCounter == 2) {
                std::string str = " *** Convergence conditions are satisfied: (2nd) ***\n";
                str += "*** SCF is well converged ***\n";
                this->logger(str);
                std::cout << str;
            }

            // conv_counterとNUMBER_OF_CHECKが一致したら収束とみなす
            if (this->m_nConvergenceCounter == NUMBER_OF_CHECK) {
                // OUTSCF_FLAG=0;
                bAnswer = true;
            }
            //}
        } else {
            // conv_counterとNUMBER_OF_CHECKが一致する前に
            // 収束の閾値を満たさないことがあれば、
            // conv_counterは0に戻る
            // conv_counter = 0;
            this->m_nConvergenceCounter = 0;
        }
    }

    return bAnswer;
}

bool DfScf::checkConverge() {
    // std::cout << ">>conv" << std::endl;
    // std::cout << (*this->pPdfParam_)["TEs"].str() << std::endl;
    // std::cout << "<<" << std::endl;

    DfConvcheck dfConvcheck(this->pPdfParam_, this->m_nIteration);

    return dfConvcheck.isConverged();
}

void DfScf::converge() {
    TlTime timer(true);
    this->loggerStartTitle("Converge");

    DfConverge* pDfConverge = this->getDfConverge();
    pDfConverge->doConverge();

    delete pDfConverge;
    pDfConverge = NULL;

    this->loggerEndTitle();
    (*this->pPdfParam_)["stat"]["elapsed_time"]["converge"][this->m_nIteration] = timer.getElapseTime();
}

DfConverge* DfScf::getDfConverge() {
    DfConverge* pDfConverge = NULL;
    if (this->m_nScfAcceleration == SCF_ACCELERATION_SIMPLE) {
        pDfConverge = new DfConverge_Damping_Lapack(this->pPdfParam_);
    } else if (this->m_nScfAcceleration == SCF_ACCELERATION_ODA) {
        pDfConverge = new DfConverge_Damping_Oda_Lapack(this->pPdfParam_);
    } else if (this->m_nScfAcceleration == SCF_ACCELERATION_ANDERSON) {
        pDfConverge = new DfConverge_Damping_Anderson_Lapack(this->pPdfParam_);
    } else if (this->m_nScfAcceleration == SCF_ACCELERATION_DIIS) {
        pDfConverge = new DfConverge_Damping_Diis_Lapack(this->pPdfParam_);
    } else {
        this->log_.info("unknown acceleration method. use damping method.");
        pDfConverge = new DfConverge_Damping_Lapack(this->pPdfParam_);
    }
    return pDfConverge;
}

void DfScf::cleanup() {
    if (TlUtils::toUpper((*this->pPdfParam_)["cleanup"].getStr()) != "NO") {
        this->loggerStartTitle("cleanup files");
        DfCleanup dfCleanup(this->pPdfParam_);
        dfCleanup.cleanup();
        this->loggerEndTitle();
    }
}

bool DfScf::checkMaxIteration() {
    bool answer = false;
    if (this->m_nIteration >= (*this->pPdfParam_)["max_iteration"].getInt()) {
        const std::string str =
            TlUtils::format(" max_iteration %d is reached.\n", (*this->pPdfParam_)["max_iteration"].getInt());
        this->logger(str);
        std::cout << str;

        answer = true;
    }

    return answer;
}
