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

#include "TlCommunicate.h"
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


void ProteinDF_Parallel::setupGlobalCondition_extra()
{
#ifdef HAVE_SCALAPACK
    {
        bool useScaLAPACK = (TlUtils::toUpper(this->pdfParam_["linear_algebra_package"].getStr())
                             == "SCALAPACK") ? true : false;
        if (useScaLAPACK == true) {
            int scalapackBlockSize = 64;
            if (this->pdfParam_.hasKey("scalapack_block_size") == true) {
                scalapackBlockSize = this->pdfParam_["scalapack_block_size"].getInt();
            }
            this->log_.info(TlUtils::format("ScaLAPACK block size: %d", scalapackBlockSize));
            TlDistributeMatrix::setSystemBlockSize(scalapackBlockSize);
            TlDistributeVector::setSystemBlockSize(scalapackBlockSize);

            bool isUsingPartialIO = false;
            if (this->pdfParam_.hasKey("save_distributed_matrix_to_local_disk") == true) {
                isUsingPartialIO = this->pdfParam_["save_distributed_matrix_to_local_disk"].getBoolean();
            }
            const std::string isUsingPartialIO_YN = (isUsingPartialIO == true) ? "YES" : "NO ";
            this->log_.info(TlUtils::format("partial I/O mode = %s\n", isUsingPartialIO_YN.c_str()));
            TlDistributeMatrix::setUsingPartialIO(isUsingPartialIO);

            // experimental 
            this->log_.info("[experimental parameters]");
            this->log_.info("use_matrix_cache parameter disabled.");
            this->pdfParam_["use_matrix_cache"] = false;

            this->pdfParam_["ERI_calcmode"] = 2;
            this->log_.info(TlUtils::format("ERI calc mode: %d", this->pdfParam_["ERI_calcmode"].getInt()));
        }
    }
#endif // HAVE_SCALAPACK
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


DfScf* ProteinDF_Parallel::createDfScfInstance()
{
    DfScf* pObj = new DfScf_Parallel(&(this->pdfParam_));
    return pObj;
}


DfForce* ProteinDF_Parallel::getDfForceObject()
{
    DfForce* pDfForce = new DfForce_Parallel(&this->pdfParam_);
    return pDfForce;
}


void ProteinDF_Parallel::startlogo()
{
    TlCommunicate& rComm = TlCommunicate::getInstance();

    const std::string version = "parallel";
    std::string info = "";
    info += TlUtils::format(" MPI process: %d\n", rComm.getNumOfProc());
#ifdef _OPENMP
    info += TlUtils::format(" OpenMP threads: %d\n", omp_get_max_threads());
#endif // _OPENMP
    
    ProteinDF::startlogo(version, info);
}

void ProteinDF_Parallel::endlogo()
{
    TlCommunicate& rComm = TlCommunicate::getInstance();
    const int rank = rComm.getRank();
    const int numOfProcs = rComm.getNumOfProcs();
    
    std::string performanceReports = "";
    performanceReports += TlUtils::format(" #%6d:\n", rank);
    performanceReports += TlUtils::format(" CPU TIME  : %9.0lf sec\n", g_GlobalTime.getCpuTime());
    performanceReports += TlUtils::format(" ELAPS_TIME: %9.0lf sec\n", g_GlobalTime.getElapseTime());

    std::string matrixCacheReports = "";
    if (this->showCacheReport_ == true) {
        matrixCacheReports += TlUtils::format(" #%6d:\n", rank);
        matrixCacheReports += TlMatrixCache::reportStats();
    }

    std::string commReports = rComm.getReport();
    
    // 集計
    std::string reports = "";
    if (rComm.isMaster() == true) {
        std::string recvMsg = "";
        reports += " performance information:\n";
        reports += performanceReports; // from '0'
        for (int i = 1; i < numOfProcs; ++i) {
            rComm.receiveData(recvMsg, i);
            reports += recvMsg;
        }
        
        if (this->showCacheReport_ == true) {
            reports += " matrix cache report:\n";
            reports += matrixCacheReports; // from '0'
            for (int i = 1; i < numOfProcs; ++i) {
                rComm.receiveData(recvMsg, i);
                reports += recvMsg;
            }
        }
        
        reports += " MPI communication report:\n";
        reports += commReports; // from '0'
        for (int i = 1; i < numOfProcs; ++i) {
            rComm.receiveData(recvMsg, i);
            reports += recvMsg;
        }
    } else {
        rComm.sendData(performanceReports, 0);
        if (this->showCacheReport_ == true) {
            rComm.sendData(matrixCacheReports, 0);
        }
        rComm.sendData(commReports, 0);
    }
    
    //
    ProteinDF::endlogo(reports);
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
