#ifdef HAVE_CONFIG_H
#include "config.h"    // this file created by autotools
#endif // HAVE_CONFIG_H

#include "DfInitialGuess_Parallel.h"
#include "DfInitialguess.h"
#include "DfDmatrix_Parallel.h"

DfInitialGuess_Parallel::DfInitialGuess_Parallel(TlSerializeData* pPdfParam)
    : DfInitialGuess(pPdfParam)
{
}


DfInitialGuess_Parallel::~DfInitialGuess_Parallel()
{
}


void DfInitialGuess_Parallel::logger(const std::string& str) const
{
    TlCommunicate& rComm = TlCommunicate::getInstance();
    if (rComm.isMaster() == true) {
        DfInitialGuess::logger(str);
    }
}


void DfInitialGuess_Parallel::createRho()
{
    TlCommunicate& rComm = TlCommunicate::getInstance();
    if (rComm.isMaster() == true) {
        DfInitialguess diguess(this->pPdfParam_);
        diguess.dfGusMain();
    }
    rComm.barrier();
}


void DfInitialGuess_Parallel::saveRho1(const RUN_TYPE runType)
{
    TlCommunicate& rComm = TlCommunicate::getInstance();
    if (rComm.isMaster() == true) {
        DfInitialGuess::saveRho1(runType);
    }
    rComm.barrier();
}


void DfInitialGuess_Parallel::createInitialGuessUsingLCAO(const RUN_TYPE runType)
{
#ifdef HAVE_SCALAPACK
    if (this->m_bUsingSCALAPACK == true) {
        this->createInitialGuessUsingLCAO_onScaLAPACK(runType);
    } else {
        this->createInitialGuessUsingLCAO_onLAPACK(runType);
    }
#else
    { 
        this->createInitialGuessUsingLCAO_onLAPACK(runType);
    }
#endif // HAVE_SCALAPACK
}


void DfInitialGuess_Parallel::createInitialGuessUsingLCAO_onScaLAPACK(const RUN_TYPE runType)
{
    // read guess lcao (all proc.)
    const TlDistributeMatrix LCAO = this->getLCAO_onScaLAPACK(runType);
    this->saveC0(runType, LCAO);

    // read guess occupation
    const TlVector aOccupation = this->getOccupation(runType);
    this->saveOccupation(runType, aOccupation);

    // output guess lcao in orthonormal basis to a files in fl_Work directory
    this->buildCprime0<TlDistributeMatrix>(runType, LCAO);

    {
        TlSerializeData tmpParam = *(this->pPdfParam_);
        tmpParam["orbital-overlap-correspondence-method"] = "keep";
        tmpParam["num_of_iterations"] = 0;
        DfDmatrix_Parallel dfDmatrix(&tmpParam);
        dfDmatrix.DfDmatrixMain(); // RKS only?
    }
}


void DfInitialGuess_Parallel::createInitialGuessUsingLCAO_onLAPACK(const RUN_TYPE runType)
{
     TlCommunicate& rComm = TlCommunicate::getInstance();

     if (rComm.isMaster() == true) {
         // read guess lcao
         const TlMatrix LCAO = DfInitialGuess::getLCAO<TlMatrix>(runType);
         this->saveC0(runType, LCAO);
         
         // read guess occupation
         const TlVector aOccupation = DfInitialGuess::getOccupation(runType);
         DfInitialGuess::saveOccupation(runType, aOccupation);
         
         {
             TlSerializeData tmpParam = *(this->pPdfParam_);
             tmpParam["orbital-overlap-correspondence-method"] = "keep";
             tmpParam["control-iteration"] = 0;
             
             // 密度行列の作成
             DfDmatrix dfDmatrix(&tmpParam);
             dfDmatrix.DfDmatrixMain(); // RKS only?
         }
     }
     
    rComm.barrier();
}


TlDistributeMatrix DfInitialGuess_Parallel::getLCAO_onScaLAPACK(const RUN_TYPE runType)
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
            CnErr.abort(TlUtils::format("cannot open file %s.\n", sFile.c_str()));
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


void DfInitialGuess_Parallel::createOccupation(const RUN_TYPE runType)
{
    TlCommunicate& rComm = TlCommunicate::getInstance();
    if (rComm.isMaster() == true) {
        DfInitialGuess::createOccupation(runType);
    }
    rComm.barrier();
}


TlVector DfInitialGuess_Parallel::getOccupation(const RUN_TYPE runType)
{
    TlVector v;
    TlCommunicate& rComm = TlCommunicate::getInstance();
    if (rComm.isMaster() == true) {
        v = DfInitialGuess::getOccupation(runType);
    }
    rComm.broadcast(v);
    return v;
}


void DfInitialGuess_Parallel::saveOccupation(const RUN_TYPE runType, const TlVector& rOccupation)
{
    TlCommunicate& rComm = TlCommunicate::getInstance();
    if (rComm.isMaster() == true) {
        const std::string sOccFileName = this->getOccupationPath(runType);
        rOccupation.save(sOccFileName);
    }
}



void DfInitialGuess_Parallel::createInitialGuessUsingHuckel()
{
    TlCommunicate& rComm = TlCommunicate::getInstance();
    if (rComm.isMaster() == true) {
        DfInitialGuess::createInitialGuessUsingHuckel();
    }
    rComm.barrier();
}


void DfInitialGuess_Parallel::createInitialGuessUsingCore()
{
    TlCommunicate& rComm = TlCommunicate::getInstance();
    if (rComm.isMaster() == true) {
        DfInitialGuess::createInitialGuessUsingCore();
    }
    rComm.barrier();
}


void DfInitialGuess_Parallel::createInitialGuessUsingHarris()
{
    TlCommunicate& rComm = TlCommunicate::getInstance();
    if (rComm.isMaster() == true) {
        DfInitialGuess::createInitialGuessUsingHarris();
    }
    rComm.barrier();
}
