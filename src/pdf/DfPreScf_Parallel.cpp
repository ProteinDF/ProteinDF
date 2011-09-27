#ifdef HAVE_CONFIG_H
#include "config.h"    // this file created by autotools
#endif // HAVE_CONFIG_H

#include "DfPreScf_Parallel.h"
#include "DfDmatrix_Parallel.h"
#include "TlCommunicate.h"
#include "TlFileMatrix.h"
#include "TlDistributeMatrix.h"

DfPreScf_Parallel::DfPreScf_Parallel(TlSerializeData* pPdfParam)
    : DfPreScf(pPdfParam)
{
}


DfPreScf_Parallel::~DfPreScf_Parallel()
{
}


void DfPreScf_Parallel::logger(const std::string& str) const
{
    TlCommunicate& rComm = TlCommunicate::getInstance();
    if (rComm.isMaster() == true) {
        DfPreScf::logger(str);
    }
}


void DfPreScf_Parallel::createInitialGuessUsingLCAO(const RUN_TYPE runType)
{
#ifdef HAVE_SCALAPACK
    if (this->m_bUsingSCALAPACK == true) {
        this->createInitialGuessUsingLCAO_onScaLAPACK(runType);
    } else {
        TlCommunicate& rComm = TlCommunicate::getInstance();
        if (rComm.isMaster() == true) {
            this->createInitialGuessUsingLCAO_onDisk(runType);
        }
        rComm.barrier();
    }
    
#else
    { 
        TlCommunicate& rComm = TlCommunicate::getInstance();
        if (rComm.isMaster() == true) {
            this->createInitialGuessUsingLCAO_onDisk(runType);
        }
        rComm.barrier();
    }
#endif // HAVE_SCALAPACK

}


void DfPreScf_Parallel::createInitialGuessUsingLCAO_onScaLAPACK(const RUN_TYPE runType)
{
    // read guess lcao
    //const TlDistributeMatrix LCAO = this->getLCAO<TlDistributeMatrix>(runType);
    const TlDistributeMatrix LCAO = this->getLCAO_onScaLAPACK(runType);
    this->saveC0(runType, LCAO);

    // read guess occupation
    const TlVector aOccupation = this->getOccupation(runType);
    this->saveOccupation(runType, aOccupation);

    // output guess lcao in orthonormal basis to a files in fl_Work directory
    this->buildCprime0(runType, LCAO);

    {
        TlSerializeData tmpParam = *(this->pPdfParam_);
        tmpParam["model"]["orbital-overlap-correspondence-method"] = "keep";
        tmpParam["model"]["iterations"] = 0;
        DfDmatrix_Parallel dfDmatrix(&tmpParam);
        dfDmatrix.DfDmatrixMain(); // RKS only?
    }
}


void DfPreScf_Parallel::createInitialGuessUsingLCAO_onDisk(const RUN_TYPE runType)
{
    TlMatrix::useMemManager(true);
    // read guess lcao
    const TlMatrix LCAO = this->getLCAO<TlMatrix>(runType);
    this->saveC0(runType, LCAO);

    // read guess occupation
    const TlVector aOccupation = this->getOccupation(runType);
    this->saveOccupation(runType, aOccupation);

    // output guess lcao in orthonormal basis to a files in fl_Work directory
    this->buildCprime0(runType, LCAO);

    {
        TlSerializeData tmpParam = *(this->pPdfParam_);
        tmpParam["model"]["orbital-overlap-correspondence-method"] = "keep";
        tmpParam["model"]["iterations"] = 0;
        DfDmatrix dfDmatrix(&tmpParam);
        dfDmatrix.DfDmatrixMain(); // RKS only?
    }
}


TlDistributeMatrix DfPreScf_Parallel::getLCAO_onScaLAPACK(const RUN_TYPE runType)
{
    TlCommunicate& rComm = TlCommunicate::getInstance();
    TlDistributeMatrix lcaoMatrix(this->m_nNumOfAOs, this->m_nNumOfMOs);

    const int numOfRows = this->m_nNumOfAOs;
    const int numOfCols = this->m_nNumOfMOs;

    if (rComm.isMaster() == true) {
        // for MASTER
        std::ifstream fi;
        const std::string sFile = std::string("./guess.lcao.") + this->m_sRunTypeSuffix[runType];
        fi.open(sFile.c_str(), std::ios::in);
        if (fi.rdstate()) {
            std::cerr << "cannot open file " << ("./guess.lcao." + sFile) << std::endl;
            CnErr.abort();
        }

        std::string dummy_line;
        fi >> dummy_line;

        int row_dimension, col_dimension;
        fi >> row_dimension >> col_dimension;
        if (row_dimension != this->m_nNumOfAOs) {
            CnErr.abort("DfPreScf", "", "prepare_occupation_and_or_mo", "inputted guess lcao has illegal dimension");
        }

        if (this->m_nNumOfMOs < col_dimension) {
            this->logger("The number of column dimension in inputed LCAO is larger than independent basis.\n");
            this->logger("Excess elements are discarded.\n");
        }

        // read matrix and store it
        const int excessCols = col_dimension - numOfCols;
        double readBuf = 0.0;
        std::vector<double> tmpColVec(numOfCols);
        for (int i = 0; i < numOfRows; ++i) {
            for (int j = 0; j < numOfCols; ++j) {
                fi >> readBuf;
                tmpColVec[j] = readBuf;
            }

            rComm.broadcast(tmpColVec);
            for (int j = 0; j < numOfCols; ++j) {
                lcaoMatrix.set(i, j, tmpColVec[j]);
            }            

            for (int j = 0; j < excessCols; ++j) {
                fi >> readBuf; // spoil
            }
        }
        fi.close();
    } else {
        // for SLAVE
        std::vector<double> tmpColVec;
        for (index_type i = 0; i < numOfRows; ++i) {
            rComm.broadcast(tmpColVec);
            assert(tmpColVec.size() == (std::size_t)numOfCols);
            for (index_type j = 0; j < numOfCols; ++j) {
                lcaoMatrix.set(i, j, tmpColVec[j]);
            }
        }
    }

    return lcaoMatrix;
}


void DfPreScf_Parallel::createOccupation(const RUN_TYPE runType)
{
    TlCommunicate& rComm = TlCommunicate::getInstance();
    if (rComm.isMaster() == true) {
        DfPreScf::createOccupation(runType);
    }
    rComm.barrier();
}

