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
#include "TlLogX.h"
#include "Fl_GlobalinputX.h"
#include "TlStringTokenizer.h"
#include "TlMemManager.h"
#include "TlMsgPack.h"
#include "TlMatrixCache.h"

ProteinDF::ProteinDF()
{
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

    this->logger("loading ProteinDF parameter for restart.\n");

    // リスタート時にはすでにあるパラメータファイルを読み取るのみ
    this->loadParam(restartParamFilePath);

    this->exec();

    this->endlogo();
}


void ProteinDF::exec()
{
    TlLogX& log = TlLogX::getInstance();

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
            log.warn(TlUtils::format("unkown control keyword: %s, continue...\n", group.c_str()));
            this->stepEndTitle();
        }

        log.flush();
    }
}

void ProteinDF::logger(const std::string& str) const
{
    TlLogX& log = TlLogX::getInstance();

    log << str;
}


void ProteinDF::startlogo()
{
    TlLogX& log = TlLogX::getInstance();
    TlTime time;

    log << "************************************************************************\n";
    log << "ProteinDF version " << VERSION << "\n";
    log << g_GlobalTime.getNowDate() << " " << g_GlobalTime.getNowTime() << "\n";
    log << "\n";
#ifdef _OPENMP
    log << TlUtils::format(" OpenMP threads: %d\n", omp_get_max_threads());
#endif // _OPENMP
    log << "\n";
    log << "copyright(c) 1997-2011 ProteinDF development team.                      \n";
    log << "\n";
    log << "PLEASE CITE following:\n";
    log << " F. Sato, Y. Shigemitsu, I. Okazaki, S. Yahiro, M. Fukue, S. Kozuru,    \n";
    log << " H. Kashiwagi, \"Development of a new density functional program for    \n";
    log << " all-electron calculation of proteins\",                                \n";
    log << " Int. J. Quant. Chem., 63, 245-246 (1997).                              \n";
    log << "\n";
    log << "************************************************************************\n";
}


void ProteinDF::endlogo()
{
    TlLogX& log = TlLogX::getInstance();

    log << "************************************************************************\n";
    log << "ProteinDF Normal Termination\n";

    // cache report
    if (this->showCacheReport_ == true) {
        log << " matrix cache report:\n";
        log << TlMatrixCache::reportStats();
        log << "\n";
    }
    
    //log << "CPU_TIME:   " << TlUtils::format("== %9.0lf sec\n", g_GlobalTime.getCpuTime());
    log << "ELAPS_TIME: " << TlUtils::format("== %9.0lf sec\n", g_GlobalTime.getElapseTime());
    log << TlTime::getNowDate() << " " << TlTime::getNowTime() << "\n";
    log << "************************************************************************\n";
}


void ProteinDF::stepStartTitle(const std::string& stepName)
{
    TlLogX& log = TlLogX::getInstance();

    log << "========================================================================\n";
    log << ">>>>" << stepName << " ";
    log << "(" << TlTime::getNowDate() << " " << TlTime::getNowTime() << ")\n";
    log << "========================================================================\n";
}

void ProteinDF::stepEndTitle()
{
    TlLogX& log = TlLogX::getInstance();

    log << "(" << TlTime::getNowDate() << " " << TlTime::getNowTime() << ")\n";
    log << "========================================================================\n";
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
    this->manageMemory();
}

void ProteinDF::manageMemory()
{
    if (TlUtils::toUpper(this->pdfParam_["use_mapfile"].getStr()) == "YES") {
        std::string filePath = this->pdfParam_["mapfile_basename"].getStr();
        if (filePath == "") {
            filePath = "/tmp/pdfmmap";
        }

        std::size_t mapFileSize = std::size_t(1024UL * 1024UL * 1024UL); // 少なくとも 1 GBは欲しい。
        std::string mapFileSizeStr = TlUtils::toUpper(this->pdfParam_["mapfile_size"].getStr());
        if (mapFileSizeStr != "AUTO") {
            // mapfile_sizeはMB単位で指定のこと。
            mapFileSize = std::max<std::size_t>(mapFileSize, std::atoi(mapFileSizeStr.c_str()));
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

        this->saveParam();
    }
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
    
    DfScf dscf(&(this->pdfParam_));
    dscf.dfScfMain();
    
    this->pdfParam_["control"]["SCF_finished"] = true;
    
    this->stepEndTitle();
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

