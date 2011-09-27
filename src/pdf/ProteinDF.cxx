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

ProteinDF::ProteinDF() : pdfParamPath_("pdfparam.mpac")
{
    this->showCacheReport_ = false;
}



ProteinDF::~ProteinDF()
{
}


void ProteinDF::run()
{
    this->startlogo();

    // input data
    this->inputData();
    this->saveParam();

    this->exec();

    this->endlogo();
}


void ProteinDF::restart()
{
    this->startlogo();

    this->logger("loading ProteinDF parameter for restart.\n");
    this->loadParam();
    this->exec();

    this->endlogo();
}


void ProteinDF::exec()
{
    TlLogX& log = TlLogX::getInstance();

    // ProteinDF class parameter
    this->showCacheReport_ = this->pdfParam_["model"]["show_cache_report"].getBoolean();

    std::string control = this->pdfParam_["model"]["step_control"].getStr();

    // setup condition
    this->setupGlobalCondition();

    TlStringTokenizer st(control);
    while (st.hasMoreTokens()) {
        std::string group = st.nextToken();
        //std::cout << group << std::endl;

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
    if (TlUtils::toUpper(this->pdfParam_["model"]["use_mapfile"].getStr()) == "YES") {
        std::string filePath = this->pdfParam_["model"]["mapfile_basename"].getStr();
        if (filePath == "") {
            filePath = "/tmp/pdfmmap";
        }

        std::size_t mapFileSize = std::size_t(1024UL * 1024UL * 1024UL); // 少なくとも 1 GBは欲しい。
        std::string mapFileSizeStr = TlUtils::toUpper(this->pdfParam_["model"]["mapfile_size"].getStr());
        if (mapFileSizeStr != "AUTO") {
            // mapfile_sizeはMB単位で指定のこと。
            mapFileSize = std::max<std::size_t>(mapFileSize, std::atoi(mapFileSizeStr.c_str()));
        } else {
            const std::size_t numOfAOs = this->pdfParam_["model"]["AOs"].getInt();
            const std::size_t numOfAuxDen = this->pdfParam_["model"]["auxCDs"].getInt();
            const std::size_t numOfAuxXC = this->pdfParam_["model"]["auxXCs"].getInt();
            const std::size_t numOfAux = std::max(numOfAuxDen, numOfAuxXC);

            const std::size_t needMem_AO = numOfAOs * numOfAOs * 3; // full matrix
            const std::size_t needMem_Aux = (numOfAux * (numOfAux +1) / 2) * 3; // half matrix

            mapFileSize = std::max(mapFileSize, needMem_AO * sizeof(double));
            mapFileSize = std::max(mapFileSize, needMem_Aux * sizeof(double));
        }
        TlMemManager::setParam(mapFileSize, filePath);

        this->pdfParam_["model"]["mapfile_size"] = mapFileSize;
        this->pdfParam_["model"]["mapfile_basename"] = filePath;

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
    
    DfScf dscf(&(this->pdfParam_), this->pdfParamPath_);
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


void ProteinDF::loadParam()
{
    TlMsgPack mpac;
    mpac.load(this->pdfParamPath_);
    this->pdfParam_ = mpac.getSerializeData();
}


void ProteinDF::saveParam() const
{
    TlMsgPack mpac(this->pdfParam_);
    mpac.save(this->pdfParamPath_);

    // conventional
    this->save_Fl_Globalinput();
}


void ProteinDF::save_Fl_Globalinput() const
{
    // support conventional QCLO program
    Fl_GlobalinputX flGbi;
    flGbi["SCF"]["method"] = this->pdfParam_["model"]["method"].getStr();
    flGbi["SCF"]["control-norb"] = this->pdfParam_["model"]["AOs"].getStr();
    flGbi["SCF"]["control-norbcut"] = this->pdfParam_["model"]["MOs"].getStr();
    flGbi["SCF"]["control-iteration"] = this->pdfParam_["model"]["iterations"].getStr();
    flGbi["SCF"]["method/nsp/electron-number"]  = this->pdfParam_["model"]["RKS/electrons"].getStr();
    flGbi["SCF"]["method/sp/alpha-elec-number"] = this->pdfParam_["model"]["UKS/alphaElectrons"].getStr();
    flGbi["SCF"]["method/sp/beta-elec-number"]  = this->pdfParam_["model"]["UKS/betaElectrons"].getStr();
    flGbi.save();
}


DfIntegrals* ProteinDF::getDfIntegralsObject()
{
    DfIntegrals* pDfIntegrals = new DfIntegrals(&(this->pdfParam_), this->pdfParamPath_);
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

