#ifdef HAVE_CONFIG_H
#include "config.h"    // this file created by autotools
#endif // HAVE_CONFIG_H

#include <iostream>
#include <fstream>
#include <cassert>

#ifdef _OPENMP
#include "omp.h"
#endif // _OPENMP

#include "ProteinDF.h"

#include "DfInputdata.h"
#include "DfIntegrals.h"
#include "DfInitialGuess.h"
#include "DfScf.h"
#include "DfForce.h"

#include "TlTime.h"
#include "TlUtils.h"
#include "TlStringTokenizer.h"
#include "TlMemManager.h"
#include "TlMsgPack.h"
#include "TlMatrixCache.h"

ProteinDF::ProteinDF() 
    : log_(TlLogging::getInstance()) {
    this->showCacheReport_ = false;
}


ProteinDF::~ProteinDF()
{
}


void ProteinDF::run()
{
    this->startlogo();

    // 入力データの解析
    // リスタート時には解析しない
    this->inputData();
    this->saveParam();

    this->exec();

    this->endlogo();
}


void ProteinDF::restart(const std::string& restartParamFilePath)
{
    this->startlogo();

    this->log_.info("loading ProteinDF parameter for restart.");

    // リスタート時にはすでにあるパラメータファイルを読み取るのみ
    this->loadParam(restartParamFilePath);

    this->exec();

    this->endlogo();
}


void ProteinDF::exec()
{
    // ProteinDF class parameter
    this->showCacheReport_ = this->pdfParam_["show_cache_report"].getBoolean();

    std::string control = this->pdfParam_["step_control"].getStr();

    // setup condition
    this->setupGlobalCondition();

    TlStringTokenizer st(control);
    while (st.hasMoreTokens()) {
        std::string group = st.nextToken();

        if (group == "create") {
            // [create]
            this->stepCreate();
        } else if (group == "integral") {
            // [integral]
            this->loadParam();
            this->stepIntegral();
            this->saveParam();
        } else if (group == "guess") {
            // [guess]
            this->loadParam();
            this->stepGuess();
            this->saveParam();
        } else if (group == "scf" || group == "scfqclo") {
            // [scf] or [scfqclo]
            this->loadParam();
            this->stepScf();
            this->saveParam();
        } else if (group == "force") {
            // [force]
            this->loadParam();
            this->stepForce();
            this->saveParam();
        } else if (group != "") {
            //nothing
            this->stepStartTitle(group);
            this->log_.warn(TlUtils::format("unkown control keyword: %s, continue...", group.c_str()));
            this->stepEndTitle();
        }
    }
}

void ProteinDF::startlogo()
{
#ifdef _OPENMP
    const std::string ompInfo = TlUtils::format(" OpenMP threads: %d\n", omp_get_max_threads());
    this->startlogo("serial", ompInfo);
#else
    this->startlogo("serial");
#endif // _OPENMP
}

void ProteinDF::startlogo(const std::string& version,
                          const std::string& info)
{
    this->log_.info("************************************************************************");
    this->log_.info(TlUtils::format("ProteinDF version %s:%s (%s)",
                                    VERSION,
                                    GIT_DESCRIBE,
                                    version.c_str()));
    this->log_.info("\n");
    this->log_.info(info);
    this->log_.info("\n");
    this->log_.info("copyright(c) 1997-2012 ProteinDF development team.");
    this->log_.info("\n");
    this->log_.info("PLEASE CITE following:");
    this->log_.info(" F. Sato, Y. Shigemitsu, I. Okazaki, S. Yahiro, M. Fukue, S. Kozuru,");
    this->log_.info(" H. Kashiwagi, \"Development of a new density functional program for");
    this->log_.info(" all-electron calculation of proteins\",");
    this->log_.info(" Int. J. Quant. Chem., 63, 245-246 (1997).");
    this->log_.info("************************************************************************");
}

void ProteinDF::endlogo()
{
    std::string reports = "";

    // cache report
    if (this->showCacheReport_ == true) {
        reports += " matrix cache report:\n";
        reports += TlMatrixCache::reportStats();
        reports += "\n";
    }

    reports += TlUtils::format("CPU_TIME:    %9.0lf sec\n", g_GlobalTime.getCpuTime());
    reports += TlUtils::format("ELAPS_TIME:  %9.0lf sec", g_GlobalTime.getElapseTime());
    
    this->endlogo(reports);
}

void ProteinDF::endlogo(const std::string& reports)
{
    this->log_.info("************************************************************************");
    this->log_.info("ProteinDF successful completion");

    this->log_.info(reports);
    
    this->log_.info("************************************************************************");
}


void ProteinDF::stepStartTitle(const std::string& stepName)
{
    this->log_.info("========================================================================");
    this->log_.info(">>>> " + stepName);
    this->log_.info("========================================================================");
}

void ProteinDF::stepEndTitle()
{
    this->log_.info("========================================================================");
}

////////////////////////////////////////////////////////////////////////
// operation
//
void ProteinDF::inputData()
{
    this->stepStartTitle("INPUT DATA");
    
    DfInputdata dfInputData;
    this->pdfParam_ = dfInputData.main();

    this->stepEndTitle();
}

void ProteinDF::setupGlobalCondition()
{
    this->stepStartTitle("RESOURCE");

    this->manageMemory();
    this->setupGlobalCondition_extra();
    
    this->stepEndTitle();
}

void ProteinDF::manageMemory()
{
    std::string memSizeStr = TlUtils::toUpper(this->pdfParam_["memory_size"].getStr());
    if (memSizeStr.empty() == true) {
        memSizeStr = "1GB";
    }
    {
        std::size_t value = std::atol(memSizeStr.c_str());
        if (memSizeStr.rfind("MB") != std::string::npos) {
            value *= (1024UL * 1024UL);
        } else if (memSizeStr.rfind("GB") != std::string::npos) {
            value *= (1024UL * 1024UL * 1024UL);
        }
        this->pdfParam_["memory_size"] = value;
    }
    this->log_.info(TlUtils::format("allocatable memory: %ld byte",
                                    this->pdfParam_["memory_size"].getLong()));

    if (this->pdfParam_["use_mapfile"].getBoolean() == true) {
        std::string filePath = this->pdfParam_["mapfile_basename"].getStr();
        if (filePath == "") {
            filePath = "/tmp/pdfmmap";
        }

        std::size_t mapFileSize = std::size_t(1024UL * 1024UL * 1024UL); // 少なくとも 1 GBは欲しい。
        std::string mapFileSizeStr = TlUtils::toUpper(this->pdfParam_["mapfile_size"].getStr());
        if (mapFileSizeStr != "AUTO") {
            // mapfile_sizeはMB単位で指定のこと。
            mapFileSize = std::max<std::size_t>(mapFileSize, std::atoi(mapFileSizeStr.c_str()));
            this->log_.info("map file size is calculated automatically.");
        } else {
            const std::size_t numOfAOs = this->pdfParam_["num_of_AOs"].getInt();
            const std::size_t numOfAuxDen = this->pdfParam_["num_of_auxCDs"].getInt();
            const std::size_t numOfAuxXC = this->pdfParam_["num_of_auxXCs"].getInt();
            const std::size_t numOfAux = std::max(numOfAuxDen, numOfAuxXC);

            const std::size_t needMem_AO = numOfAOs * numOfAOs * 3; // full matrix
            const std::size_t needMem_Aux = (numOfAux * (numOfAux +1) / 2) * 3; // half matrix

            mapFileSize = std::max(mapFileSize, needMem_AO * sizeof(double));
            mapFileSize = std::max(mapFileSize, needMem_Aux * sizeof(double));
        }
        TlMemManager::setParam(mapFileSize, filePath);

        this->pdfParam_["mapfile_size"] = mapFileSize;
        this->pdfParam_["mapfile_basename"] = filePath;
        this->log_.info(TlUtils::format("map file size: %ld byte", mapFileSize));
        this->log_.info(TlUtils::format("map file basename: %s", filePath.c_str()));

        this->saveParam();
    }
}


void ProteinDF::setupGlobalCondition_extra()
{
    // do nothing 
}

void ProteinDF::stepCreate()
{
}


void ProteinDF::stepIntegral()
{
    this->loadParam();
    
    this->stepStartTitle("INTEGRAL");

    DfIntegrals* pDfIntegrals = this->getDfIntegralsObject();
    pDfIntegrals->main();
    
    delete pDfIntegrals;
    pDfIntegrals = NULL;
    
    this->stepEndTitle();
}


void ProteinDF::stepGuess()
{
    this->stepStartTitle("GUESS");

    DfInitialGuess* pDfInitialGuess = this->getDfInitialGuessObject();
    pDfInitialGuess->exec();
    this->pdfParam_["control"]["guess_finished"] = true;
    
    delete pDfInitialGuess;
    pDfInitialGuess = NULL;
    
    this->stepEndTitle();
}


void ProteinDF::stepScf()
{
    this->stepStartTitle("SCF");

    DfScf* pDfScf = this->createDfScfInstance();
    pDfScf->dfScfMain();
    
    delete pDfScf;
    pDfScf = NULL;
    
    this->stepEndTitle();
}

DfScf* ProteinDF::createDfScfInstance()
{
    DfScf* pObj = new DfScf(&(this->pdfParam_));
    return pObj;
}


void ProteinDF::stepForce()
{
    this->stepStartTitle("FORCE");

    DfForce* pDfForce = this->getDfForceObject();
    pDfForce->calcForce();
    
    this->stepEndTitle();
}


void ProteinDF::loadParam(const std::string& requestFilePath)
{
    std::string filePath = requestFilePath;
    if (requestFilePath.empty() == true) {
        filePath = this->pdfParam_["pdf_param_path"].getStr();
    }
    
    TlMsgPack mpac;
    mpac.load(filePath);
    this->pdfParam_ = mpac.getSerializeData();

    // check
    if (requestFilePath.empty() != true) {
        const std::string pdfParamPath = this->pdfParam_["pdf_param_path"].getStr();
        if (requestFilePath != pdfParamPath) {
            std::cerr << "the specified parameter file path is not consistent with pdf_param_path parameter."
                      << std::endl;
        }
    }
}


void ProteinDF::saveParam() const
{
    TlMsgPack mpac(this->pdfParam_);
    const std::string pdfParamPath = this->pdfParam_["pdf_param_path"].getStr();
    mpac.save(pdfParamPath);
}


DfIntegrals* ProteinDF::getDfIntegralsObject()
{
    DfIntegrals* pDfIntegrals = new DfIntegrals(&(this->pdfParam_));
    return pDfIntegrals;
}


DfInitialGuess* ProteinDF::getDfInitialGuessObject()
{
    DfInitialGuess* pDfInitialGuess = new DfInitialGuess(&this->pdfParam_);
    return pDfInitialGuess;
}


DfForce* ProteinDF::getDfForceObject()
{
    DfForce* pDfForce = new DfForce(&this->pdfParam_);
    return pDfForce;
}

