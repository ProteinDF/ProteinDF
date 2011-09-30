#include <cmath>

#include "DfScf.h"
#include "TlUtils.h"

#include "DfDiffDensityMatrix.h"
#include "DfDensityFitting.h"
#include "DfCalcGrid.h"
#include "DfThreeindexintegrals.h"
#include "DfXCFunctional.h"
#include "DfJMatrix.h"
#include "DfFockMatrix.h"
#include "DfTransFmatrix.h"
#include "DfDiagonal.h"
#include "DfTransatob.h"
#include "DfDmatrix.h"
#include "DfTotalEnergy.h"
#include "DfPopulation.h"
#include "DfSummary.h"
#include "DfConvcheck.h"
#include "DfXcpotfitting.h"

#include "DfConverge_Damping.h"
#include "DfConverge_Anderson.h"

#include "DfConverge2.h"

#include "DfXcenefitting.h"
#include "DfLevelshift.h"

#include "DfCleanup.h"

// FoR Extended QCLO
#include "DfQclo.h"
#include "DfCqclomatrix.h"

#include "DfTwoElectronIntegral.h"
#include "DfPreScf.h"

#include "TlFile.h"
#include "TlLogX.h"
#include "TlMsgPack.h"

#include "Fl_GlobalinputX.h"

#define NUMBER_OF_CHECK 2

DfScf::DfScf(TlSerializeData* pPdfParam)
    : DfObject(pPdfParam), m_nDampObject(DAMP_NONE)
{
    this->isUseNewEngine_ = (*pPdfParam)["new_engine"].getBoolean();
}


DfScf::~DfScf()
{
}


void DfScf::saveParam() const
{
    (*(this->pPdfParam_))["iterations"] = this->m_nIteration;

    // save PDF parameter
    const std::string pdfParamPath = (*this->pPdfParam_)["pdf_param_path"].getStr();
    TlMsgPack pdfParam_mpac(*(this->pPdfParam_));
    pdfParam_mpac.save(pdfParamPath);
}


// return  0 : not convergence
//         1 : convergence
int DfScf::dfScfMain()
{
    const TlSerializeData& pdfParam = *(this->pPdfParam_);
    this->setScfParam();

    this->logger(" restart calculation is " + pdfParam["restart"].getStr() + "\n");

    std::string sStepControl = pdfParam["step_control"].getStr();
    std::string group = "";
    do {
        group = TlUtils::toUpper(TlUtils::getWord(sStepControl));
    } while ((group != "SCF") && (group != "SCFQCLO"));

    // QCLO法用。PreSCFに持って行くべき
    if (group == "SCFQCLO") {
        if ((this->isRestart_ == false) ||
            (pdfParam["DfScf"]["scf-restart-point"] == "startPDF")) {
            this->loggerStartTitle("DfCqclomatrix");

            DfCqclomatrix dfcqclomatrix(this->pPdfParam_);
            dfcqclomatrix.main();

            this->loggerEndTitle();
        }
    }

    // guess の作成
    // preSCFに持って行くべき。
    if ((this->m_nIteration == 0) || (this->isRestart_ == false)) {
        this->m_nIteration = 1;
    }

    // start SCF LOOP
    this->saveParam();
    return this->execScfLoop();
}

void DfScf::setScfParam()
{
    const TlSerializeData& pdfParam = *(this->pPdfParam_);
    
    // iteration number
    this->m_nIteration = 0;
    if (this->isRestart_ == true) {
        this->m_nIteration = std::max<int>(1, pdfParam["iterations"].getInt());
    }

    // damping switch
    this->m_nScfAcceleration = SCF_ACCELERATION_SIMPLE;
    {
        const std::string sScfAcceleration =
            TlUtils::toUpper(pdfParam["scf-acceleration"].getStr());

        if (sScfAcceleration == "DAMPING") {
            this->m_nScfAcceleration = SCF_ACCELERATION_SIMPLE;
        } else if (sScfAcceleration == "ANDERSON") {
            this->m_nScfAcceleration = SCF_ACCELERATION_ANDERSON;
        }
    }

    // DIIS
    this->diisflg = false;
    if ((pdfParam["scf-acceleration"] == "diis") ||
        (pdfParam["scf-acceleration"] == "mix")) {
        this->diisflg = true;
    }
    this->diisworkflg = false;
    this->diiscycle = 0; // diis stream specifier -> off

    // Damp Object Type
    this->m_nDampObject = DAMP_DENSITY;
    if (this->isRI_J_ == false) {
        this->m_nDampObject = DAMP_DENSITY_MATRIX;
    }
    {
        const std::string sDampObject =
            TlUtils::toUpper(pdfParam["scf-acceleration/damping/damping-type"].getStr());
        
        if (sDampObject == "DENSITY_MATRIX") {
            this->m_nDampObject = DAMP_DENSITY_MATRIX;
        } else if (sDampObject == "DENSITY") {
            this->m_nDampObject = DAMP_DENSITY;
        } else if (sDampObject == "FOCK") {
            this->m_nDampObject = DAMP_FOCK;
        }
    }
}


int DfScf::execScfLoop()
{
    const TlSerializeData& pdfParam = *(this->pPdfParam_);
    //int OUTSCF_FLAG = 1;
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
        XCENEFIT,
        THREE_INDEX_INTEGRAL,
        XC_MATRIX,
        J_MATRIX,
        FOCK,
        ENDFOCK,
        TRANSFORM_FOCK,
        LEVEL_SHIFT,
        DIAGONAL,
        ENDFOCK_TRANSC,
        DENSITY_MATRIX,
        DIIS_TOTALENE,
        TOTAL_ENERGY,
        POPULATION,
        SUMMARY,
        JUDGE,
        JUDGE_CONV2,
        JUDGE_TAIL,
        END_OF_SCF_LOOP,
        END_SCF_LOOP
    };

    SCF_STATE nScfState = BEGIN;
    if (this->isRestart_ == true) {
        const std::string respoint = TlUtils::toUpper(pdfParam["control"]["scf_state"].getStr());

        if (respoint == "prepareGuess") {
            nScfState = BEGIN;
        } else if (respoint == "DIFF_DENSITY_MATRIX") {
            nScfState = DENSITY_FITTING;
        } else if ("DENSITY_FITTING" == respoint) {
            nScfState = XC_INTEGRAL;
        } else if ("XC_INTEGRAL" == respoint) {
            nScfState = XCENEFIT;
        } else if ("XCENEFIT" == respoint) {
            nScfState = THREE_INDEX_INTEGRAL;
        } else if ("THREE_INDEX_INTEGRAL" == respoint) {
            nScfState = XC_MATRIX;
        } else if ("XC_MATRIX" == respoint) {
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
            nScfState = DIIS_TOTALENE;
        } else if ("TOTAL_ENERGY" == respoint) {
            nScfState = POPULATION;
        } else if ("POPULATION" == respoint) {
            nScfState = SUMMARY;
        } else if ("SUMMARY" == respoint) {
            nScfState = JUDGE;
        } else {
            this->logger(" restart calculation is indicated, but state is not defined.\n");
            nScfState = BEGIN;
        }
    }

    while (nScfState != END_SCF_LOOP) {
        switch (nScfState) {
        case BEGIN:
            this->logger("// =====================================\n");
            this->logger(TlUtils::format("number of iteration = %d\n", this->m_nIteration));
            this->logger("// =====================================\n");
            nScfState = DIFF_DENSITY_MATRIX;
            break;

        case DIFF_DENSITY_MATRIX:
            this->diffDensityMatrix();
            this->setScfRestartPoint("DIFF_DENSITY_MATRIX");
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
            nScfState = XCENEFIT;
            break;

        case XCENEFIT:
            this->execScfLoop_XcEneFit();
            this->setScfRestartPoint("XCENEFIT");
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
            nScfState = J_MATRIX;
            break;

        case J_MATRIX:
            this->setScfRestartPoint("J_MATRIX");
            this->buildJMatrix();
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

                // if (diiscycle == 1) {
                //     this->logger("DIIS内挿Fpqを処理する\n");
                // }
                bool bExecDiis = (this->diiscycle == 1) ? true : false;
                DfQclo dfQclo(this->pPdfParam_, this->m_nIteration, bExecDiis);
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
            nScfState = DIIS_TOTALENE;
            break;

        case DIIS_TOTALENE:
            if (diiscycle == 1) {
                // this->logger("DIIS内挿Fpqを処理した\n");
                diiscycle = 0;

                nScfState = JUDGE_TAIL;
            } else {
                nScfState = TOTAL_ENERGY;
            }
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
                nScfState = JUDGE_CONV2;
            } else {
                nScfState = END_OF_SCF_LOOP;
            }

            break;

        case JUDGE_CONV2:
            if (this->diisflg == true) {
                this->loggerStartTitle("DfConverge2");

                DfConverge2 dfconverge2(this->pPdfParam_, this->m_nIteration);
                this->diisworkflg = dfconverge2.DfConv2Main();

                this->loggerEndTitle();

                // if DIIS not work well, changed to damping
                if (this->diisworkflg == false) {
                    nScfState = JUDGE_TAIL;
                } else {
                    diiscycle = 1;
                    nScfState = TRANSFORM_FOCK;
                }
            } else {
                nScfState = JUDGE_TAIL;
            }
            break;

        case JUDGE_TAIL:
            {
                this->cleanup();
                const bool isExitScf = this->checkMaxIteration();
                
                if (isExitScf == false) {
                    ++(this->m_nIteration);
                    (*(this->pPdfParam_))["iterations"] = this->m_nIteration;
                    this->saveParam();
                    nScfState = BEGIN;
                } else {
                    nScfState = END_OF_SCF_LOOP;
                }
            }
            break;

        case END_OF_SCF_LOOP:
            this->saveParam();
            nScfState = END_SCF_LOOP;
            break;

        default:
            std::cerr << "unknown ScfState(" << (int)nScfState << "). stop." << std::endl;
            exit(1);
            break;
        }

    }

    // pupulation analysis
    const std::string sAnalizePopulation = TlUtils::toUpper(pdfParam["analyze_population"].getStr());
    if ((sAnalizePopulation == "EVERY-SCF") ||
        (sAnalizePopulation == "CONVERGENCE")) {
        this->calcPopulation();
    }

    //calculation of Total energy, energy part which comes from dummy charge
    this->calcTotalRealEnergy();

    // DfSummary
    const std::string sDfSummary = TlUtils::toUpper(pdfParam["summary"].getStr());
    if ((sDfSummary == "EVERY-SCF") ||
        (sDfSummary == "CONVERGENCE")) {
        this->summarize();
    }

    return 0;
}

void DfScf::setScfRestartPoint(const std::string& str)
{
    (*(this->pPdfParam_))["control"]["scf_state"] = str;
    this->saveParam();
}


void DfScf::diffDensityMatrix()
{
    if (this->m_nDampObject == DAMP_DENSITY_MATRIX) {
        this->converge();
    }

    if (this->m_bDiskUtilization == false) {
        TlTime timer;
        this->loggerStartTitle("diff density matrix");
        DfDiffDensityMatrix dddm(this->pPdfParam_);
        dddm.exec();

        this->loggerEndTitle();
        (*this->pPdfParam_)["stat"]["elapse_time"]["diff_density_matrix"][this->m_nIteration] = timer.getElapseTime();
    }
}


void DfScf::doDensityFitting()
{
    if (this->isRI_J_ == true) {
        if (this->m_nIteration == 1) {
            if ((this->initialGuessType_ == GUESS_RHO) ||
                (this->initialGuessType_ == GUESS_FILE_RHO)) {
                this->logger("Initial rho is provided. Density fitting is unnecessary.\n");
                return;
            }
        }
        
        TlTime timer;
        this->loggerStartTitle("Density Fitting");

        DfDensityFittingObject* pDfDensityFitting = this->getDfDensityFittingObject();
        pDfDensityFitting->exec();
        delete pDfDensityFitting;
        pDfDensityFitting = NULL;
        
        this->loggerEndTitle();
        (*this->pPdfParam_)["stat"]["elapse_time"]["density_fitting"][this->m_nIteration] = timer.getElapseTime();
        
        if (this->m_nDampObject == DAMP_DENSITY) {
            this->converge();
        }
    }
}


DfDensityFittingObject* DfScf::getDfDensityFittingObject()
{
    DfDensityFittingObject* pDfDensityFittingObj = NULL;    
    if (this->isUseNewEngine_ == true) {
        pDfDensityFittingObj = new DfDensityFittingX(this->pPdfParam_);
    } else {
        pDfDensityFittingObj = new DfDensityFitting(this->pPdfParam_);
    }

    return pDfDensityFittingObj;
}


void DfScf::doXCIntegral()
{
    if (this->m_bIsXCFitting == true) {
        if (this->m_sXCFunctional != "xalpha") {
            // grid-method

            this->loggerStartTitle("DfGrid fitting Myu");

            DfCalcGrid dg(this->pPdfParam_, this->m_nIteration);
            dg.dfGrdMain();

            this->loggerEndTitle();
        } else {
            // xalpha-method
            if ((*(this->pPdfParam_))["SCF"]["xc-potential/xalpha/method"] == "newton-raphson") {
                // newton-raphson
                this->loggerStartTitle("DfXcpotfitting");

                DfXcpotfitting dfXcpotfitting(this->pPdfParam_, this->m_nIteration);
                dfXcpotfitting.dfXcpMain();

                this->loggerEndTitle();
            } else {
                // broyden
                this->loggerStartTitle("DfBroyden");
                this->loggerEndTitle();
            }
        }
    }
}


void DfScf::execScfLoop_XcEneFit()
{
    // x.c. energy fitting
    if (this->m_sXCFunctional == "xalpha" ||
        this->m_sXCFunctional == "gxalpha") {
        this->loggerStartTitle("Xc ene fitting");

        DfXcenefitting dfXcenefitting(this->pPdfParam_, this->m_nIteration);
        dfXcenefitting.dfXceMain();

        this->loggerEndTitle();
    }
}


void DfScf::doThreeIndexIntegral()
{
    if ((this->m_bMemorySave != true) &&
        (this->m_bIsXCFitting == true)) {
        this->loggerStartTitle("DfThreeindexintegrals");
        DfThreeindexintegrals dfThreeindexintegrals(this->pPdfParam_);
        dfThreeindexintegrals.DfThreeindexintegralsMain();
        this->loggerEndTitle();
    }
}


void DfScf::buildXcMatrix()
{
    if (this->m_bIsXCFitting == false) {
        // for restart
        if (this->isRestart_ == true) {
            const std::string prevGridDataFilePath = TlUtils::format("%s.itr%d",
                                                                     this->getGridDataFilePath().c_str(),
                                                                     this->m_nIteration -1);
            TlFile::copy(prevGridDataFilePath, this->getGridDataFilePath());
        }

        TlTime timer;
        this->loggerStartTitle("generate XC matrix");
        DfXCFunctional* pDfXCFunctional = this->getDfXCFunctional();
        pDfXCFunctional->buildXcMatrix();
        this->loggerEndTitle();
        (*this->pPdfParam_)["stat"]["elapse_time"]["xc_matrix"][this->m_nIteration] = timer.getElapseTime();

        delete pDfXCFunctional;
        pDfXCFunctional = NULL;

        this->saveParam();

        // flush
        this->matrixCache_.flush();
    }
}


DfXCFunctional* DfScf::getDfXCFunctional()
{
    DfXCFunctional* pDfXCFunctional = new DfXCFunctional(this->pPdfParam_);
    return pDfXCFunctional;
}


void DfScf::buildJMatrix()
{
    if (this->isRI_J_ == false) {
        this->loggerStartTitle("Coulomb term");
        DfJMatrix* pDfJMatrix = this->getDfJMatrixObject();
        pDfJMatrix->buildJMatrix();
        this->loggerEndTitle();
    }
}


DfJMatrix* DfScf::getDfJMatrixObject()
{
    DfJMatrix* pDfJMatrix = new DfJMatrix(this->pPdfParam_);
    return pDfJMatrix;
}


void DfScf::buildFock()
{
    TlTime timer;
    this->loggerStartTitle("Fock matrix");
    DfFockMatrix* pDfFockMatrix = this->getDfFockMatrixObject();
    pDfFockMatrix->DfFockMatrixMain();
    this->loggerEndTitle();
    (*this->pPdfParam_)["staat"]["elapse_time"]["fock_matrix"][this->m_nIteration] = timer.getElapseTime();

    if (this->m_nDampObject == DAMP_FOCK) {
        this->converge();
    }
 
    delete pDfFockMatrix;
    pDfFockMatrix = NULL;

    // flush
    this->matrixCache_.flush();
}


DfFockMatrix* DfScf::getDfFockMatrixObject()
{
    DfFockMatrix* pDfFockMatrix = new DfFockMatrix(this->pPdfParam_);
    return pDfFockMatrix;
}


void DfScf::transformFock()
{
    // transformed to orth. A.O. based Fock matrix
    TlTime timer;
    this->loggerStartTitle("Transform KS matrix");
    bool isExecDiis = (this->diiscycle == 1) ? true : false;

    DfTransFmatrix* pDfTransFmatrix = this->getDfTransFmatrixObject(isExecDiis);
    pDfTransFmatrix->DfTrsFmatMain();
    delete pDfTransFmatrix;
    pDfTransFmatrix = NULL;

    this->loggerEndTitle();
    (*this->pPdfParam_)["stat"]["elapse_time"]["transform_F_matrix"][this->m_nIteration] = timer.getElapseTime();
}


DfTransFmatrix* DfScf::getDfTransFmatrixObject(bool isExecDiis)
{
    DfTransFmatrix* pDfTransFmatrix = new DfTransFmatrix(this->pPdfParam_, isExecDiis);
    return pDfTransFmatrix;
}


void DfScf::doLevelShift()
{
    // add level shift to Kohn-Sham matrix
    const int start_iter = (*(this->pPdfParam_))["level-shift/start-iteration"].getInt();
    const bool levelShift = (*(this->pPdfParam_))["level-shift"].getBoolean();
    if ((levelShift == true) &&
        (this->m_nIteration >= start_iter)) {
        TlTime timer;
        this->loggerStartTitle("Level shift");

        DfLevelshift LS(this->pPdfParam_, this->m_nIteration);
        LS.DfLshiftMain();

        this->loggerEndTitle();
        (*this->pPdfParam_)["stat"]["elapse_time"]["level_shift"][this->m_nIteration] = timer.getElapseTime();
    }
}


void DfScf::diagonal()
{
    // Diagonarize Fock matrix
    TlTime timer;
    this->loggerStartTitle("Diagonal");
    DfDiagonal* pDfDiagonal = this->getDfDiagonalObject();
    pDfDiagonal->DfDiagMain();
    delete pDfDiagonal;
    pDfDiagonal = NULL;
    this->loggerEndTitle();
    (*this->pPdfParam_)["stat"]["elapse_time"]["diagonal"][this->m_nIteration] = timer.getElapseTime();

    // flush
    this->matrixCache_.flush();
}


DfDiagonal* DfScf::getDfDiagonalObject()
{
    DfDiagonal* pDfDiagonal = new DfDiagonal(this->pPdfParam_);
    return pDfDiagonal;
}


void DfScf::execScfLoop_EndFock_TransC()
{
    // transformed to original nonorth. A.O.based space
    TlTime timer;
    this->loggerStartTitle("Transform Matrix");
    DfTransatob* pDfTransAtoB = this->getDfTransatobObject();
    pDfTransAtoB->DfTrsatobMain();
    delete pDfTransAtoB;
    pDfTransAtoB = NULL;
    this->loggerEndTitle();
    (*this->pPdfParam_)["stat"]["elapse_time"]["transform_C_matrix"][this->m_nIteration] = timer.getElapseTime();
}


DfTransatob* DfScf::getDfTransatobObject()
{
    DfTransatob* pDfTransAtoB = new DfTransatob(this->pPdfParam_);
    return pDfTransAtoB;
}


void DfScf::calcDensityMatrix()
{
    // density matrix generation
    TlTime timer;
    this->loggerStartTitle("Density Matirx");
    DfDmatrix* pDfDmatrix = this->getDfDmatrixObject();
    pDfDmatrix->DfDmatrixMain();
    delete pDfDmatrix;
    pDfDmatrix = NULL;
    this->loggerEndTitle();
    (*this->pPdfParam_)["stat"]["elapse_time"]["density_matrix"][this->m_nIteration] = timer.getElapseTime();

    // flush
    this->matrixCache_.flush();
}


DfDmatrix* DfScf::getDfDmatrixObject()
{
    DfDmatrix* pDfDmatrix = new DfDmatrix(this->pPdfParam_);
    return pDfDmatrix;
}


void DfScf::calcTotalEnergy()
{
    // calculate total energy
    TlTime timer;
    this->loggerStartTitle("Total Energy");
    DfTotalEnergy* pDfTotalEnergy = this->getDfTotalEnergyObject();
    pDfTotalEnergy->exec();
    delete pDfTotalEnergy;
    pDfTotalEnergy = NULL;
    this->loggerEndTitle();
    (*this->pPdfParam_)["stat"]["elapse_time"]["total_energy"][this->m_nIteration] = timer.getElapseTime();
}


DfTotalEnergy* DfScf::getDfTotalEnergyObject()
{
    DfTotalEnergy* pDfTotalEnergy = new DfTotalEnergy(this->pPdfParam_);
    return pDfTotalEnergy;
}


void DfScf::calcTotalRealEnergy()
{
    // calculate total energy
    this->loggerStartTitle("Total Energy derived from point charges");
    DfTotalEnergy* pDfTotalEnergy = this->getDfTotalEnergyObject();
    pDfTotalEnergy->calculate_real_energy();
    delete pDfTotalEnergy;
    pDfTotalEnergy = NULL;
    this->loggerEndTitle();
}


DfPopulation* DfScf::getDfPopulationObject()
{
    DfPopulation* pDfPopulation = new DfPopulation(this->pPdfParam_);
    return pDfPopulation;
}


void DfScf::calcPopulation()
{
    TlLogX& log = TlLogX::getInstance();

    // pupulation analysis
    this->loggerStartTitle("Population analysis");

    DfPopulation* pDfPopulation = this->getDfPopulationObject();
    pDfPopulation->getReport(this->m_nIteration, log);
    
    delete pDfPopulation;
    pDfPopulation = NULL;
    this->loggerEndTitle();
}


DfSummary* DfScf::getDfSummaryObject()
{
    DfSummary* pDfSummary = new DfSummary(this->pPdfParam_);
    return pDfSummary;
}


void DfScf::summarize()
{
    //std::cout << "DfScf::summarize() called." << std::endl;
    this->loggerStartTitle("Summary");
    DfSummary* pDfSummary = this->getDfSummaryObject();
    pDfSummary->exec();
    delete pDfSummary;
    pDfSummary = NULL;
    
    this->loggerEndTitle();
    //std::cout << "DfScf::summarize() exit." << std::endl;
}


bool DfScf::judge()
{
    //std::cerr << "enter JUDGE" << std::endl;
    //OUTSCF_FLAG=1;

    bool bAnswer = false;
    {
        this->loggerStartTitle("Convergence Check");

        int bJudge = this->checkConverge();
//     {
//       DfConvcheck dfConvcheck(this->m_flGbi, niter);
//       dfConvcheck.DfConvcheckMain();
//     }

        this->loggerEndTitle();

        //  set convergence
        //if ((convergence = dfConvcheck.judge()) != 0){
//     if (dfConvcheck.judge() != 0){
        if (bJudge != 0) {
            // 収束の閾値をみたすとconv_counterが１増える
            //conv_counter++;
            this->m_nConvergenceCounter++;

            //if (conv_counter == 1){
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
            //if(conv_counter == NUMBER_OF_CHECK) {
            if (this->m_nConvergenceCounter == NUMBER_OF_CHECK) {
                //OUTSCF_FLAG=0;
                bAnswer = true;
            }
            //}
        } else {
            // conv_counterとNUMBER_OF_CHECKが一致する前に
            // 収束の閾値を満たさないことがあれば、
            // conv_counterは0に戻る
            //conv_counter = 0;
            this->m_nConvergenceCounter = 0;
        }
    }

    return bAnswer;
}

int DfScf::checkConverge()
{
    DfConvcheck dfConvcheck(this->pPdfParam_, this->m_nIteration);
    dfConvcheck.DfConvcheckMain();

    return dfConvcheck.judge();
}


void DfScf::converge()
{
    if (this->m_nIteration == 1) {
        return;
    }

    if ((this->m_nScfAcceleration == SCF_ACCELERATION_SIMPLE) ||
        (this->m_nScfAcceleration == SCF_ACCELERATION_ANDERSON) ||
        ((this->diisflg == true) && (this->diisworkflg == false))) {

        TlTime timer;
        this->loggerStartTitle("Converge");

        DfConverge* pDfConverge = this->getDfConverge();
        pDfConverge->doConverge();

        delete pDfConverge;
        pDfConverge = NULL;

        this->loggerEndTitle();
        (*this->pPdfParam_)["stat"]["elapse_time"]["converge"][this->m_nIteration] = timer.getElapseTime();
    }
}


DfConverge* DfScf::getDfConverge()
{
    DfConverge* pDfConverge = NULL;
    if (this->m_nScfAcceleration == SCF_ACCELERATION_SIMPLE) {
        pDfConverge = new DfConverge_Damping(this->pPdfParam_);
    } else if (this->m_nScfAcceleration == SCF_ACCELERATION_ANDERSON) {
        pDfConverge = new DfConverge_Anderson(this->pPdfParam_);
    } else {
        // diis 法の最初のdampingなど
        pDfConverge = new DfConverge_Damping(this->pPdfParam_);
    }
    return pDfConverge;
}


void DfScf::cleanup()
{
    if (TlUtils::toUpper((*this->pPdfParam_)["cleanup"].getStr()) != "NO") {
        this->loggerStartTitle("cleanup files");
        DfCleanup dfCleanup(this->pPdfParam_);
        dfCleanup.cleanup();
        this->loggerEndTitle();
    }
}


bool DfScf::checkMaxIteration()
{
    bool answer = false;
    if (this->m_nIteration >= (*this->pPdfParam_)["max-iteration"].getInt()) {
        const std::string str = TlUtils::format(" max-iteration %d is reached.\n",
                                                (*this->pPdfParam_)["max-iteration"].getInt());
        this->logger(str);
        std::cout << str;

        answer = true;
    }

    return answer;
}


