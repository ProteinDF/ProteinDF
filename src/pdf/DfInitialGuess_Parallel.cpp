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

#ifdef HAVE_CONFIG_H
#include "config.h"    // this file created by autotools
#endif // HAVE_CONFIG_H

#include "DfInitialGuess_Parallel.h"
#include "DfInitialGuessHarris_Parallel.h"
#include "DfInitialGuessHuckel_Parallel.h"
#include "DfDmatrix_Parallel.h"

DfInitialGuess_Parallel::DfInitialGuess_Parallel(TlSerializeData* pPdfParam)
    : DfInitialGuess(pPdfParam)
{
}


DfInitialGuess_Parallel::~DfInitialGuess_Parallel()
{
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

    this->makeDensityMatrix();
}


void DfInitialGuess_Parallel::createInitialGuessUsingLCAO_onLAPACK(const RUN_TYPE runType)
{
     TlCommunicate& rComm = TlCommunicate::getInstance();
     if (rComm.isMaster() == true) {
         DfInitialGuess::createInitialGuessUsingLCAO(runType);
    }
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
            this->log_.info("The number of column dimension in inputed LCAO is larger than independent basis.");
            this->log_.info("Excess elements are discarded.");
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


TlVector DfInitialGuess_Parallel::createOccupation(const RUN_TYPE runType)
{
    TlVector answer;
    TlCommunicate& rComm = TlCommunicate::getInstance();
    if (rComm.isMaster() == true) {
        answer = DfInitialGuess::createOccupation(runType);
    }
    rComm.broadcast(answer);
    return answer;
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
#ifdef HAVE_SCALAPACK
    if (this->m_bUsingSCALAPACK) {
        DfInitialGuessHuckel_Parallel huckel(this->pPdfParam_);
        huckel.createGuess();
    } else {
        TlCommunicate& rComm = TlCommunicate::getInstance();
        if (rComm.isMaster() == true) {
            DfInitialGuess::createInitialGuessUsingHuckel();
        }
        rComm.barrier();
    }
#else
    {
        TlCommunicate& rComm = TlCommunicate::getInstance();
        if (rComm.isMaster() == true) {
            DfInitialGuess::createInitialGuessUsingHuckel();
        }
        rComm.barrier();
    }
#endif // HAVE_SCALAPACK
}


void DfInitialGuess_Parallel::createInitialGuessUsingCore()
{
#ifdef HAVE_SCALAPACK
    if (this->m_bUsingSCALAPACK) {
        DfInitialGuessHuckel_Parallel huckel(this->pPdfParam_);
        huckel.createGuess();
    } else {
        TlCommunicate& rComm = TlCommunicate::getInstance();
        if (rComm.isMaster() == true) {
            DfInitialGuess::createInitialGuessUsingCore();
        }
        rComm.barrier();
    }
#else
    {
        TlCommunicate& rComm = TlCommunicate::getInstance();
        if (rComm.isMaster() == true) {
            DfInitialGuess::createInitialGuessUsingCore();
        }
        rComm.barrier();
    }
#endif // HAVE_SCALAPACK
}


void DfInitialGuess_Parallel::createInitialGuessUsingHarris()
{
    switch (this->m_nMethodType) {
    case METHOD_RKS:
        {
            DfInitialGuessHarris_Parallel harris(this->pPdfParam_);
            harris.main();
        }
        break;

    case METHOD_UKS:
        CnErr.abort("Sorry. harris method is not supported except RKS. stop.\n");
        break;

    case METHOD_ROKS:
        CnErr.abort("Sorry. harris method is not supported except RKS. stop.\n");
        break;

    default:
        CnErr.abort();
        break;
    }
}

DfDmatrix* DfInitialGuess_Parallel::getDfDmatrixObject(TlSerializeData* param)
{
    DfDmatrix* obj = new DfDmatrix_Parallel(param);
    return obj;
}
