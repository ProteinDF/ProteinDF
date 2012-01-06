#include <cassert>
#include "DfGenerateGrid_Parallel.h"
#include "DfXCFunctional.h"
#include "TlCommunicate.h"
//#include "GridDataManager.h"
#include "TlUtils.h"
#include "TlFileMatrix.h"

DfGenerateGrid_Parallel::DfGenerateGrid_Parallel(TlSerializeData* pPdfParam)
    : DfGenerateGrid(pPdfParam)
{
}

DfGenerateGrid_Parallel::~DfGenerateGrid_Parallel()
{
}

void DfGenerateGrid_Parallel::logger(const std::string& str) const
{
    TlCommunicate& rComm = TlCommunicate::getInstance();

    if (rComm.isMaster() == true) {
        DfGenerateGrid::logger(str);
    }
}

void DfGenerateGrid_Parallel::makeTable()
{
    TlCommunicate& rComm = TlCommunicate::getInstance();
    if (rComm.isMaster() == true) {
        DfGenerateGrid::makeTable();
    }

    rComm.broadcast(this->maxRadii_);
}


void DfGenerateGrid_Parallel::generateGrid()
{
    TlCommunicate& rComm = TlCommunicate::getInstance();
    const int numOfProcs = rComm.getNumOfProc();
    
    // if (this->isMasterSlave_ == true) {
    //     this->generateGrid_MS();
    // } else {
        this->generateGrid_DC();
    // }

    // gather
    index_type numOfRowsOfGlobalGridMatrix = this->grdMat_.getNumOfRows();
    rComm.allReduce_SUM(numOfRowsOfGlobalGridMatrix);
    if (rComm.isMaster() == true) {
        const index_type numOfColsOfGlobalGridMatrix = this->grdMat_.getNumOfCols();
        TlFileMatrix grdMat(this->getGridMatrixPath(0),
                            numOfRowsOfGlobalGridMatrix,
                            numOfColsOfGlobalGridMatrix);

        index_type currentNumOfRows = 0;
        // from 0
        grdMat.setBlockMatrix(currentNumOfRows, 0,
                              this->grdMat_);
        currentNumOfRows += this->grdMat_.getNumOfRows();
        
        // from the others
        for (int i = 1; i < numOfProcs; ++i) {
            TlMatrix tmpGrdMat;
            rComm.receiveData(tmpGrdMat, i);
            grdMat.setBlockMatrix(currentNumOfRows, 0,
                                  tmpGrdMat);
            currentNumOfRows += tmpGrdMat.getNumOfRows();
        }
    } else {
        rComm.sendData(this->grdMat_, 0);
    }
}

void DfGenerateGrid_Parallel::generateGrid_DC()
{
    TlCommunicate& rComm = TlCommunicate::getInstance();

    const int nEndAtomNumber = this->m_nNumOfAtoms;

    const int nProc = rComm.getNumOfProc();
    const int nRank = rComm.getRank();
    const int nRange = nEndAtomNumber;
    const int nInterval = (nRange + (nProc -1)) / nProc; // +(nProc-1) は余り対策
    const int nLocalStart = nRank * nInterval; // nProc = 0, 1, 2, ...
    const int nLocalEnd   = std::min((nLocalStart + nInterval), nEndAtomNumber);

    std::size_t numOfGrids = 0;
    for (int atom = nLocalStart; atom < nLocalEnd; ++atom) {
        std::vector<double> coordX;
        std::vector<double> coordY;
        std::vector<double> coordZ;
        std::vector<double> weight;

        if (this->m_gridType == SG_1) {
            DfGenerateGrid::generateGrid_SG1(atom, &coordX, &coordY, &coordZ, &weight);
        } else {
            DfGenerateGrid::generateGrid(atom, &coordX, &coordY, &coordZ, &weight);
        }

        // store grid matrix
        const std::size_t numOfAtomGrids = weight.size();
        this->grdMat_.resize(numOfGrids + numOfAtomGrids, this->numOfColsOfGrdMat_);
        for (std::size_t i = 0; i < numOfAtomGrids; ++i) {
            this->grdMat_.set(numOfGrids, 0, coordX[i]);
            this->grdMat_.set(numOfGrids, 1, coordY[i]);
            this->grdMat_.set(numOfGrids, 2, coordZ[i]);
            this->grdMat_.set(numOfGrids, 3, weight[i]);
            this->grdMat_.set(numOfGrids, 4, atom);
            ++numOfGrids;
        }
    }

}

