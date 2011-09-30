#include <iostream>
#include <cassert>

#ifdef _OPENMP
#include "omp.h"
#endif // _OPENMP
#include "ProteinDF_Parallel.h"

#include "DfInputdata_Parallel.h"
#include "DfIntegrals_Parallel.h"
#include "DfInitialGuess_Parallel.h"
#include "DfScf_Parallel.h"
#include "DfForce_Parallel.h"

#include "TlTime.h"
#include "TlCommunicate.h"
#include "TlLogX.h"
#include "TlUtils.h"

#ifdef HAVE_SCALAPACK
#include "TlDistributeMatrix.h"
#include "TlDistributeVector.h"
#endif // HAVE_SCALAPACK

ProteinDF_Parallel::ProteinDF_Parallel()
{
}


ProteinDF_Parallel::~ProteinDF_Parallel()
{
}


void ProteinDF_Parallel::logger(const std::string& str) const
{
    TlCommunicate& rComm = TlCommunicate::getInstance();
    if (rComm.isMaster() == true) {
        ProteinDF::logger(str);
    }
}


void ProteinDF_Parallel::save_Fl_Globalinput() const
{
    TlCommunicate& rComm = TlCommunicate::getInstance();

    if (rComm.isMaster() == true) {
        ProteinDF::save_Fl_Globalinput();
    }
}


void ProteinDF_Parallel::loadParam(const std::string& requestFilePath)
{
    TlCommunicate& rComm = TlCommunicate::getInstance();
    if (rComm.isMaster() == true) {
        ProteinDF::loadParam(requestFilePath);
    }
    rComm.broadcast(this->pdfParam_);
}


void ProteinDF_Parallel::saveParam() const
{
    TlCommunicate& rComm = TlCommunicate::getInstance();
    if (rComm.isMaster() == true)
    {
        ProteinDF::saveParam();
    }
    rComm.barrier();
}


void ProteinDF_Parallel::inputData()
{
    this->stepStartTitle("INPUT DATA");

    DfInputdata_Parallel dfInputData;
    this->pdfParam_ = dfInputData.main();
    
    this->stepEndTitle();
}


void ProteinDF_Parallel::setupGlobalCondition()
{
#ifdef HAVE_SCALAPACK
    {
        bool useScaLAPACK = (TlUtils::toUpper(this->pdfParam_["model"]["linear_algebra_package"].getStr())
                             == "SCALAPACK") ? true : false;
        if (useScaLAPACK == true) {
            int scalapackBlockSize = 64;
            if (this->pdfParam_["model"].hasKey("scalapack_block_size") == true) {
                scalapackBlockSize = this->pdfParam_["model"]["scalapack_block_size"].getInt();
            }
            this->logger(TlUtils::format("ScaLAPACK block size(= %d) is set.\n", scalapackBlockSize));
            TlDistributeMatrix::setSystemBlockSize(scalapackBlockSize);
            TlDistributeVector::setSystemBlockSize(scalapackBlockSize);

            bool isUsingPartialIO = false;
            if (this->pdfParam_["model"].hasKey("save_distributed_matrix_to_local_disk") == true) {
                isUsingPartialIO = this->pdfParam_["model"]["save_distributed_matrix_to_local_disk"].getBoolean();
            }
            const std::string isUsingPartialIO_YN = (isUsingPartialIO == true) ? "YES" : "NO ";
            this->logger(TlUtils::format("partial I/O mode = %s\n", isUsingPartialIO_YN.c_str()));
            TlDistributeMatrix::setUsingPartialIO(isUsingPartialIO);
        }
    }
#endif // HAVE_SCALAPACK
    
    ProteinDF::setupGlobalCondition();
}


DfIntegrals* ProteinDF_Parallel::getDfIntegralsObject()
{
    DfIntegrals* pDfIntegrals = new DfIntegrals_Parallel(&(this->pdfParam_));
    return pDfIntegrals;
}


DfInitialGuess* ProteinDF_Parallel::getDfInitialGuessObject()
{
    DfInitialGuess* pDfInitialGuess = new DfInitialGuess_Parallel(&this->pdfParam_);
    return pDfInitialGuess;
}


void ProteinDF_Parallel::stepScf()
{
    TlCommunicate& rComm = TlCommunicate::getInstance();

    if (this->pdfParam_["control"]["SCF_finished"].getBoolean() != true) {
        this->stepStartTitle("SCF");

        DfScf_Parallel dscf(&(this->pdfParam_));
        dscf.dfScfMain();

        this->stepEndTitle();
    } else {
        this->logger(" SCF need not to execute.\n");
    }

    rComm.barrier();
}


DfForce* ProteinDF_Parallel::getDfForceObject()
{
    DfForce* pDfForce = new DfForce_Parallel(&this->pdfParam_);
    return pDfForce;
}


void ProteinDF_Parallel::startlogo()
{
    TlCommunicate& rComm = TlCommunicate::getInstance();

    if (rComm.isMaster() == true) {

        //ProteinDF::startlogo();
        TlLogX& log = TlLogX::getInstance();
        TlTime time;

        log << "************************************************************************\n";
        log << "ProteinDF version " << VERSION << "(parallel)\n";
        log << g_GlobalTime.getNowDate() << " " << g_GlobalTime.getNowTime() << "\n";
        log << "\n";
        log << TlUtils::format(" MPI process: %d\n", rComm.getNumOfProc());
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
}

void ProteinDF_Parallel::endlogo()
{
    TlCommunicate& rComm = TlCommunicate::getInstance();
    const int numOfProcs = rComm.getNumOfProcs();
    const int rank = rComm.getRank();
    
    this->logger("************************************************************************\n");
    this->logger("ProteinDF Normal Termination\n");
    this->logger(TlUtils::format(" MPI Process: %d\n", rComm.getNumOfProc()));
#ifdef _OPENMP
    this->logger(TlUtils::format(" OpenMP threads: %d\n", omp_get_max_threads()));
#endif // _OPENMP

    this->logger(TlUtils::format("\n"));
    this->logger(TlUtils::format(" performance information:\n"));
    rComm.barrier();
    for (int i = 0; i < numOfProcs; ++i) {
        if (i == rank) {
            ProteinDF::logger(TlUtils::format(" #%6d:\n", rank));
            ProteinDF::logger(TlUtils::format(" CPU TIME  : %9.0lf sec\n", g_GlobalTime.getCpuTime()));
            ProteinDF::logger(TlUtils::format(" ELAPS_TIME: %9.0lf sec\n", g_GlobalTime.getElapseTime()));
        }
        rComm.barrier();
    }

    // matrix cache report
    if (this->showCacheReport_ == true) {
        this->logger("\n");
        this->logger(" matrix cache report:\n");
        rComm.barrier();
        for (int i = 0; i < numOfProcs; ++i) {
            if (i == rank) {
                ProteinDF::logger(TlUtils::format(" #%6d:\n", i));
                ProteinDF::logger(TlMatrixCache::reportStats());
                ProteinDF::logger("\n");
            }
            rComm.barrier();
        }
    }

    // communication report
    {
        this->logger("\n");
        this->logger(" MPI communication report:\n");
        rComm.barrier();
        for (int i = 0; i < numOfProcs; ++i) {
            if (i == rank) {
                ProteinDF::logger(rComm.getReport());
            }
            rComm.barrier();
        }
    }
    
    rComm.barrier();
    this->logger(TlUtils::format(" %s %s\n", TlTime::getNowDate().c_str(), TlTime::getNowTime().c_str()));
    this->logger("************************************************************************\n");
}

void ProteinDF_Parallel::stepStartTitle(const std::string& stepName)
{
    TlCommunicate& rComm = TlCommunicate::getInstance();
    if (rComm.isMaster() == true) {
        ProteinDF::stepStartTitle(stepName);
    }
}

void ProteinDF_Parallel::stepEndTitle()
{
    TlCommunicate& rComm = TlCommunicate::getInstance();
    if (rComm.isMaster() == true) {
        ProteinDF::stepEndTitle();
    }
}
