#ifdef HAVE_CONFIG_H
#include "config.h"    // this file created by autotools
#endif // HAVE_CONFIG_H

#ifdef _OPENMP
#include <omp.h>
#endif // _OPENMP

#include <iostream>
#include <cstdlib>
#include "DfObject.h"
#include "TlLogging.h"
#include "TlUtils.h"
#include "TlSystem.h"
#include "TlTime.h"

const std::string DfObject::m_sWorkDirPath = "fl_Work";
const std::string DfObject::m_sTempDirPath = "fl_Work";

int DfObject::objectCount_ = 0;
TlMatrixCache DfObject::matrixCache_;
int DfObject::rank_ = 0;

const std::string DfObject::m_sRunTypeSuffix[DfObject::RUN_MAXINDEX] = {
    "undefined",
    "rks",
    "uks-alpha",
    "uks-beta",
    "roks",
    "roks_close",
    "roks_open"
};


DfObject::DfObject(TlSerializeData* pPdfParam)
    : pPdfParam_(pPdfParam), log_(TlLogging::getInstance())
{
    this->setParam(*pPdfParam);
    ++DfObject::objectCount_;
}


DfObject::~DfObject()
{
    --DfObject::objectCount_;
    if (DfObject::objectCount_ == 0) {
        // DfObjectの個数が0になったので、行列cacheを空にする
        this->matrixCache_.flush();
    }
}


void DfObject::setParam(const TlSerializeData& data)
{
    int numOfThreads = 1;
#ifdef _OPENMP
    numOfThreads = omp_get_max_threads();
#endif // _OPENMP

    // computational resource
    this->procMaxMemSize_ = 1024UL * 1024UL * 1024UL;
    {
        const std::string memSizeStr = TlUtils::toUpper(data["memory_size"].getStr());
        if (memSizeStr.empty() == false) {
            std::size_t value = std::atol(memSizeStr.c_str());
            if (memSizeStr.rfind("MB") != std::string::npos) {
                value *= (1024UL * 1024UL);
            } else if (memSizeStr.rfind("GB") != std::string::npos) {
                value *= (1024UL * 1024UL * 1024UL);
            }

            this->procMaxMemSize_ = value;
        }
    }

    this->isWorkOnDisk_ = (TlUtils::toUpper(data["work_on_disk"].getStr()) == "YES");
    this->localTempDir_ = data["local_temp_dir"].getStr();
    if (this->localTempDir_ == "") {
        this->localTempDir_ = "/tmp/";
    }

    this->isRestart_ = false;
    if (TlUtils::toUpper(data["restart"].getStr()) == "YES") {
        this->isRestart_ = true;
    }

    this->isUseNewEngine_ = data["new_engine"].getBoolean();

    // SCF type
    const std::string sMethodType = TlUtils::toUpper(data["method"].getStr());
    if (sMethodType == "NSP") {
        this->m_nMethodType = METHOD_RKS;
    } else if (sMethodType == "SP") {
        this->m_nMethodType = METHOD_UKS;
    } else if (sMethodType == "ROKS") {
        this->m_nMethodType = METHOD_ROKS;
    } else {
        //this->logger("no method is specified. use RKS method.\n");
        this->m_nMethodType = METHOD_RKS;
    }

    // model
    this->m_nNumOfAtoms = data["num_of_atoms"].getInt();
    this->m_nNumOfDummyAtoms = data["num_of_dummy_atoms"].getInt();
    this->numOfRealAtoms_ = this->m_nNumOfAtoms - this->m_nNumOfDummyAtoms;
    
    this->m_nIteration = data["num_of_iterations"].getInt();
    this->m_nNumOfAOs = data["num_of_AOs"].getInt();
    this->m_nNumOfMOs = data["num_of_MOs"].getInt();
    this->m_nNumOfAux = data["num_of_auxCDs"].getInt();
    this->numOfAuxXC_ = data["num_of_auxXCs"].getInt();

    this->m_nNumOfElectrons = data["RKS/electrons"].getInt();
    this->m_nNumOfAlphaElectrons = data["UKS/alphaElectrons"].getInt();
    this->m_nNumOfBetaElectrons = data["UKS/betaElectrons"].getInt();

    // guess
    this->initialGuessType_ = GUESS_UNKNOWN;
    {
        const std::string guess = TlUtils::toUpper(data["scf-start-guess"].getStr());
        if ((guess == "RHO") || (guess == "ATOM_RHO")) {
            this->initialGuessType_ = GUESS_RHO;
        } else if (guess == "FILE_RHO") {
            this->initialGuessType_ = GUESS_FILE_RHO;
        } else if ((guess == "LCAO") || (guess == "FILE_LCAO")) {
            this->initialGuessType_ = GUESS_LCAO;
        } else if (guess == "DENSITY") {
            this->initialGuessType_ = GUESS_DENSITY;
        } else if (guess == "HUCKEL") {
            this->initialGuessType_ = GUESS_HUCKEL;
        } else if (guess == "CORE") {
            this->initialGuessType_ = GUESS_CORE;
        } else if (guess == "HARRIS") {
            this->initialGuessType_ = GUESS_HARRIS;
        } else {
            if (guess.empty() == true) {
                std::cerr << "initial guess parameter is not configured." << std::endl;
            } else {
                std::cerr << "unknown initial guess parameter: " << guess << std::endl;
            }
        }
    }

    // calculaton properties ===================================================
    this->chargeExtrapolateNumber_ = data["charge_extrapolate_number"].getInt(); 
    // disk utilization(no == DIRECT, yes == DISK)
    this->m_bDiskUtilization = (TlUtils::toUpper(data["disk-utilization"].getStr()) == "YES") ? true : false;
    this->m_bMemorySave = (TlUtils::toUpper(data["scf-memory-saving"].getStr()) == "YES") ? true : false;

    // XC potential
    this->m_sXCFunctional = TlUtils::toUpper(data["xc-potential"].getStr());
    {
        const char nLastChar = this->m_sXCFunctional[this->m_sXCFunctional.length() -1];
        this->m_bIsXCFitting = (nLastChar == '~') ? true : false;
    }
    this->m_bIsUpdateXC = (TlUtils::toUpper(data["xc-update"].getStr()) == "NO") ? false : true;

    // Grimme empirical dispersion check
    {
        this->enableGrimmeDispersion_ = false;
        int len = this->m_sXCFunctional.size();
        if (len > 2) {
            const std::string suffix = this->m_sXCFunctional.substr(len -2, 2);
            if (suffix == "-D") {
                this->enableGrimmeDispersion_ = true;
            }
        }
    }

    // RI_J
    this->isRI_J_ = data["RI_J"].getBoolean();
    
    // RI-K
    this->isRI_K_ = data["RI_K"].getBoolean();
    
    // matrix operation
    this->m_bUsingSCALAPACK = false;
#ifdef HAVE_SCALAPACK
    this->m_bUsingSCALAPACK = (TlUtils::toUpper(data["linear_algebra_package"].getStr())
                               == "SCALAPACK") ? true : false;
#endif // HAVE_SCALAPACK

    this->isSaveDistributedMatrixToLocalDisk_ = false;
#ifdef HAVE_SCALAPACK
    if (data.hasKey("save_distributed_matrix_to_local_disk") == true) {
        this->isSaveDistributedMatrixToLocalDisk_ = data["save_distributed_matrix_to_local_disk"].getBoolean();
    }
#endif // HAVE_SCALAPACK

    this->localDiskPath_ = "/tmp";
    if (data.hasKey("local_disk_path") == true) {
        this->localDiskPath_ = data["local_disk_path"].getStr();
    }

    // for HPC ============================================================
    this->grainSize_ = 100 * numOfThreads;
    if (data.hasKey("omp_grain_size") == true) {
        this->grainSize_ = data["omp_grain_size"].getInt() * numOfThreads;
    }
    
    this->isMasterSlave_ = false;
    {
        const std::string parallelProcessingType = TlUtils::toUpper(data["parallel_processing_type"].getStr());
        if ((parallelProcessingType == "MASTER_SLAVE") ||
            (parallelProcessingType == "MASTER-SLAVE") ||
            (parallelProcessingType == "MS")) {
            this->isMasterSlave_ = true;
        }
    }

    
    // for DEBUG ===============================================================
    this->isFileWarning = data["debug/file_warning"].getBoolean();
    this->isSaveJMatrix_ = data["debug/save_J"].getBoolean();
    this->enableExperimentalCode_ = data["experimental_code"].getBoolean();
    
    // for memory ==============================================================
    this->isUseCache_ = (*(this->pPdfParam_))["use_matrix_cache"].getBoolean();
    this->matrixCache_.setMaxMemSize(this->procMaxMemSize_);
    const bool isForceLoadingFromDisk = (*(this->pPdfParam_))["force_loading_from_disk"].getBoolean();
    this->matrixCache_.forceLoadingFromDisk(isForceLoadingFromDisk);
    
    // setup
    TlSerializeData& paramFileBaseName = (*(this->pPdfParam_))["control"]["file_base_name"];
    paramFileBaseName["Hpq_matrix"]    = "Hpq.mat";
    paramFileBaseName["Hpq2_matrix"]   = "Hpq2.mat";
    paramFileBaseName["Spq_matrix"]    = "Spq.mat";
    paramFileBaseName["Sab_matrix"]    = "Sab.mat";
    paramFileBaseName["Sab2_matrix"]   = "Sab2.mat";
    paramFileBaseName["Sgd_matrix"]    = "Sgd.mat";
    if (paramFileBaseName["SabInv_matrix"].getStr().empty() == true) {
        paramFileBaseName["SabInv_matrix"] = "Sabinv.mat";
    }
    paramFileBaseName["SgdInv_matrix"] = "Sgdinv.mat";
    paramFileBaseName["X_matrix"]      = "X.mat";
    paramFileBaseName["Xinv_matrix"]   = "Xinv.mat";
    if (paramFileBaseName["diff_density_matrix"].getStr().empty() == true) {
        paramFileBaseName["diff_density_matrix"] = "dP.mat";
    }

    if (paramFileBaseName["Ppq_matrix"].getStr().empty() == true) {
        paramFileBaseName["Ppq_matrix"] = "Ppq%s.mat";
    }
    paramFileBaseName["P1pq_matrix"]    = "P1pq%s.mat";
    paramFileBaseName["P2pq_matrix"]    = "P2pq%s.mat";
    paramFileBaseName["HFx_matrix"]     = "HFx%s.mat";
    paramFileBaseName["Fpq_matrix"]     = "Fpq%s.mat";
    paramFileBaseName["Fprime_matrix"]  = "Fprime%s.mat";
    paramFileBaseName["Fxc_matrix"]     = "Fxc%s.mat";
    paramFileBaseName["FxcPure_matrix"] = "FxcPure%s.mat";
    paramFileBaseName["J_matrix"]       = "J%s.mat";
    paramFileBaseName["C_matrix"]       = "C%s.mat";
    paramFileBaseName["Cprime_matrix"]  = "Cprime%s.mat";
    paramFileBaseName["grid_matrix"]    = "grid.mat";
    paramFileBaseName["Talpha.vtr"]     = "Talpha.vtr";

    if (paramFileBaseName["rho_vector"].getStr().empty() == true) {
        paramFileBaseName["rho_vector"] = "rho%s.vtr";
    }
    if (paramFileBaseName["myu_vector"].getStr().empty() == true) {
        paramFileBaseName["myu_vector"] = "myu%s.vtr";
    }
    if (paramFileBaseName["nyu_vector"].getStr().empty() == true) {
        paramFileBaseName["nyu_vector"] = "nyu%s.vtr";
    }
    if (paramFileBaseName["eigenvalues"].getStr().empty() == true) {
        paramFileBaseName["eigenvalues"] = "eigenvalues%s.vtr";
    }
}


void DfObject::logger(const std::string& str) const
{
    //TlLogX& log = TlLogX::getInstance();
    //log << str;
    this->log_.info(str);
}

void DfObject::loggerTime(const std::string& str) const
{
    // std::string out = str;
    // int size = out.size();
    // if (size > 0) {
    //     if (out[size -1] == '\n') {
    //         out.erase(size -1, 1);
    //     }

    //     const std::string timeStr = "[" + TlTime::getNow() + "]";
    //     TlUtils::pad(out, (72 - timeStr.length()), ' ');
    //     out += (timeStr + "\n");
    //     this->logger(out);
    // }
    this->log_.info(str);
}

void DfObject::loggerStartTitle(const std::string& stepName, const char lineChar) const
{
    // std::string line = "";
    // TlUtils::pad(line, 72, lineChar);

    // const std::string timeString = TlUtils::format("[%s %s]",
    //                                                TlTime::getNowDate().c_str(),
    //                                                TlTime::getNowTime().c_str());

    // std::string title = ">>>> " + stepName + " ";
    // TlUtils::pad(title, (72 - timeString.length()), ' ');
    // title += timeString;

    // // 出力
    // const std::string str = line + "\n" + title + "\n";

    std::string line = "";
    TlUtils::pad(line, 50, lineChar);
    const std::string str = line + "\n" + ">>>> " + stepName + "\n";
    this->logger(str);
}

void DfObject::loggerEndTitle(const std::string& stepName, const char lineChar) const
{
    // const std::string timeString = TlUtils::format("[%s %s]",
    //                                                TlTime::getNowDate().c_str(),
    //                                                TlTime::getNowTime().c_str());

    // std::string title = "<<<< " + stepName + " ";
    // TlUtils::pad(title, (72 - timeString.length()), ' ');
    // title += timeString;

    // // 出力
    // const std::string str = "\n" + title + "\n";

    const std::string str = "<<<<\n";
    this->logger(str);
}

// =====================================================================
std::string DfObject::makeFilePath(const std::string& baseFileName,
                                   const std::string& suffix) const
{
    std::string base = (*(this->pPdfParam_))["control"]["file_base_name"][baseFileName].getStr();
    if (suffix.empty() != true) {
        base = TlUtils::format(base.c_str(), suffix.c_str());
    }

    std::string path;
    if (this->isSaveDistributedMatrixToLocalDisk_ == true) {
        path = this->localDiskPath_ + "/" + base + "." + TlUtils::xtos(this->rank_);
    } else {
        path = DfObject::m_sWorkDirPath + "/" + base;
    }

    return path;
}

std::string DfObject::getHpqMatrixPath()
{
    return this->makeFilePath("Hpq_matrix");
}


std::string DfObject::getHpq2MatrixPath()
{
    return this->makeFilePath("Hpq2_matrix");
}


std::string DfObject::getSpqMatrixPath()
{
    return this->makeFilePath("Spq_matrix");
}


std::string DfObject::getSabMatrixPath()
{
    return this->makeFilePath("Sab_matrix");
}


std::string DfObject::getSab2MatrixPath()
{
    return this->makeFilePath("Sab2_matrix");
}


std::string DfObject::getSgdMatrixPath()
{
    return this->makeFilePath("Sgd_matrix");
}


std::string DfObject::getSabInvMatrixPath()
{
    return this->makeFilePath("SabInv_matrix");
}


std::string DfObject::getSgdInvMatrixPath()
{
    return this->makeFilePath("SgdInv_matrix");
}


std::string DfObject::getXMatrixPath()
{
    return this->makeFilePath("X_matrix");
}


std::string DfObject::getInvXMatrixPath()
{
    return this->makeFilePath("Xinv_matrix");
}


std::string DfObject::getNalphaPath()
{
    return (DfObject::m_sWorkDirPath + "/fl_Vct_Nalpha");
}

std::string DfObject::getOccupationPath(const RUN_TYPE runType)
{
    std::string sFileName = DfObject::m_sWorkDirPath + "/fl_Occupation";
    if (runType == RUN_UKS_ALPHA) {
        sFileName += "_Alpha";
    } else if (runType == RUN_UKS_BETA) {
        sFileName += "_Beta";
    }

    return sFileName;
}

std::string DfObject::getGridDataFilePath() const
{
    return DfObject::m_sWorkDirPath + "/grids.dat";
}


std::string DfObject::getGridMatrixPath(const int iteration) const
{
    return this->makeFilePath("grid_matrix",
                              TlUtils::xtos(iteration));
}


std::string DfObject::getDiffDensityMatrixPath(const RUN_TYPE runType, const int iteration) const
{
    return this->makeFilePath("diff_density_matrix",
                              DfObject::m_sRunTypeSuffix[runType] + TlUtils::xtos(iteration));
}


std::string DfObject::getPpqMatrixPath(const RUN_TYPE nRunType, const int nIteration) const
{
    return this->makeFilePath("Ppq_matrix",
                              DfObject::m_sRunTypeSuffix[nRunType] + TlUtils::xtos(nIteration));
}

std::string DfObject::getP1pqMatrixPath(const int iteration)
{
    return this->makeFilePath("P1pq_matrix",
                              TlUtils::xtos(iteration));
}

std::string DfObject::getP2pqMatrixPath(const int iteration)
{
    return this->makeFilePath("P2pq_matrix",
                              TlUtils::xtos(iteration));
}

std::string DfObject::getHFxMatrixPath(const RUN_TYPE nRunType, const int nIteration)
{
    return this->makeFilePath("HFx_matrix",
                              DfObject::m_sRunTypeSuffix[nRunType] + TlUtils::xtos(nIteration));
}

std::string DfObject::getFpqMatrixPath(const RUN_TYPE nRunType, const int nIteration) const
{
    return this->makeFilePath("Fpq_matrix",
                              DfObject::m_sRunTypeSuffix[nRunType] + TlUtils::xtos(nIteration));
}


std::string DfObject::getFprimeMatrixPath(const RUN_TYPE runType, const int iteration,
                                          const std::string& fragment)
{
    std::string suffix = "";
    if (fragment.empty() != true) {
        suffix = fragment + ".";
    }
    suffix += DfObject::m_sRunTypeSuffix[runType] + TlUtils::xtos(iteration);
    return this->makeFilePath("Fprime_matrix",
                              suffix);
}


std::string DfObject::getFxcMatrixPath(const RUN_TYPE nRunType, const int nIteration)
{
    return this->makeFilePath("Fxc_matrix",
                              DfObject::m_sRunTypeSuffix[nRunType] + TlUtils::xtos(nIteration));
}


std::string DfObject::getFxcPureMatrixPath(const RUN_TYPE nRunType, const int nIteration)
{
    return this->makeFilePath("FxcPure_matrix",
                              DfObject::m_sRunTypeSuffix[nRunType] + TlUtils::xtos(nIteration));
}


std::string DfObject::getJMatrixPath(const int iteration)
{
    return this->makeFilePath("J_matrix",
                              TlUtils::xtos(iteration));
}


std::string DfObject::getCMatrixPath(const RUN_TYPE runType, int iteration, const std::string& fragment) const
{
    std::string suffix = "";
    if (fragment.empty() != true) {
        suffix = fragment + ".";
    }
    suffix += DfObject::m_sRunTypeSuffix[runType] + TlUtils::xtos(iteration);
    return this->makeFilePath("C_matrix",
                              suffix);
}


std::string DfObject::getCprimeMatrixPath(const RUN_TYPE runType, int iteration, const std::string& fragment)
{
    std::string suffix = "";
    if (fragment.empty() != true) {
        suffix = fragment + ".";
    }
    suffix += DfObject::m_sRunTypeSuffix[runType] + TlUtils::xtos(iteration);
    return this->makeFilePath("Cprime_matrix",
                              suffix);
}


std::string DfObject::getRhoPath(const RUN_TYPE nRunType, const int nIteration) const
{
    return this->makeFilePath("rho_vector",
                              DfObject::m_sRunTypeSuffix[nRunType] + TlUtils::xtos(nIteration));

}

std::string DfObject::getMyuPath(const RUN_TYPE nRunType, const int nIteration) const
{
    return this->makeFilePath("myu_vector",
                              DfObject::m_sRunTypeSuffix[nRunType] + TlUtils::xtos(nIteration));
}

std::string DfObject::getNyuPath(const RUN_TYPE nRunType, const int nIteration) const
{
    return this->makeFilePath("myu_vector",
                              DfObject::m_sRunTypeSuffix[nRunType] + TlUtils::xtos(nIteration));
}

std::string DfObject::getTalphaPath(const RUN_TYPE runType, const int iteration) const
{
    return this->makeFilePath("Talpha.vtr",
                              DfObject::m_sRunTypeSuffix[runType] + TlUtils::xtos(iteration));

}

std::string DfObject::getEigenvaluesPath(const RUN_TYPE runType,
                                         const int iteration) const
{
    return this->makeFilePath("eigenvalues",
                              DfObject::m_sRunTypeSuffix[runType] + TlUtils::xtos(iteration));
}

