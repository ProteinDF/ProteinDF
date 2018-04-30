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

#include "DfGenerateGrid_Parallel.h"
#include <cassert>
#include "DfXCFunctional.h"
#include "TlCommunicate.h"
#include "TlFileMatrix.h"
#include "TlUtils.h"

DfGenerateGrid_Parallel::DfGenerateGrid_Parallel(TlSerializeData* pPdfParam)
    : DfGenerateGrid(pPdfParam) {}

DfGenerateGrid_Parallel::~DfGenerateGrid_Parallel() {}

void DfGenerateGrid_Parallel::logger(const std::string& str) const {
  TlCommunicate& rComm = TlCommunicate::getInstance();

  if (rComm.isMaster() == true) {
    DfGenerateGrid::logger(str);
  }
}

void DfGenerateGrid_Parallel::makeTable() {
  TlCommunicate& rComm = TlCommunicate::getInstance();
  if (rComm.isMaster() == true) {
    DfGenerateGrid::makeTable();
  }

  rComm.broadcast(this->maxRadii_);
}

TlMatrix DfGenerateGrid_Parallel::getOMatrix() {
  TlCommunicate& rComm = TlCommunicate::getInstance();

  TlMatrix O;
  if (rComm.isMaster() == true) {
    O = DfGenerateGrid::getOMatrix();
  }
  rComm.broadcast(O);

  return O;
}

void DfGenerateGrid_Parallel::generateGrid(const TlMatrix& O) {
  // if (this->isMasterSlave_ == true) {
  //     this->generateGrid_MS();
  // } else {
  this->generateGrid_DC(O);
  // }

  this->gatherGridData();
}

void DfGenerateGrid_Parallel::generateGrid_DC(const TlMatrix& O) {
  this->logger("generate grid by DC");
  TlCommunicate& rComm = TlCommunicate::getInstance();

  const int nEndAtomNumber = this->m_nNumOfAtoms;

  const int nProc = rComm.getNumOfProc();
  const int nRank = rComm.getRank();
  const int nRange = nEndAtomNumber;
  const int nInterval =
      (nRange + (nProc - 1)) / nProc;         // +(nProc-1) は余り対策
  const int nLocalStart = nRank * nInterval;  // nProc = 0, 1, 2, ...
  const int nLocalEnd = std::min((nLocalStart + nInterval), nEndAtomNumber);

  std::size_t numOfGrids = 0;
  for (int atom = nLocalStart; atom < nLocalEnd; ++atom) {
    std::vector<double> coordX;
    std::vector<double> coordY;
    std::vector<double> coordZ;
    std::vector<double> weight;

    if ((this->m_gridType == SG_1) || (this->m_gridType == USER)) {
      DfGenerateGrid::generateGrid_SG1(O, atom, &coordX, &coordY, &coordZ,
                                       &weight);
    } else {
      DfGenerateGrid::generateGrid(O, atom, &coordX, &coordY, &coordZ, &weight);
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

void DfGenerateGrid_Parallel::gatherGridData() {
  this->logger("gather grid data");

  TlCommunicate& rComm = TlCommunicate::getInstance();
  const int numOfProcs = rComm.getNumOfProc();

  index_type numOfRowsOfGlobalGridMatrix = this->grdMat_.getNumOfRows();

  const int tag = TAG_GENGRID_GATHER_GRID_DATA;
  if (rComm.isMaster() == true) {
    const index_type numOfColsOfGlobalGridMatrix = this->grdMat_.getNumOfCols();
    // TODO: use TlFileMatrix instead of TlMatrix because of memory waste.
    this->log_.info(TlUtils::format("grid matrix size: %d, %d",
                                    numOfRowsOfGlobalGridMatrix,
                                    numOfRowsOfGlobalGridMatrix));
    TlMatrix grdMat(numOfRowsOfGlobalGridMatrix, numOfColsOfGlobalGridMatrix);

    index_type currentNumOfRows = 0;

    // set elements from 0
    grdMat.setBlockMatrix(currentNumOfRows, 0, this->grdMat_);
    currentNumOfRows += this->grdMat_.getNumOfRows();

    // set elements from the others
    std::vector<bool> recvCheck(numOfProcs, false);
    for (int i = 1; i < numOfProcs; ++i) {
      int proc = 0;
      TlMatrix tmpGrdMat;
      rComm.receiveDataFromAnySource(tmpGrdMat, &proc, tag);
      if (recvCheck[proc] != false) {
        this->log_.warn(
            TlUtils::format("already received grid data from %d", proc));
      }
      recvCheck[proc] = true;
      this->log_.debug(TlUtils::format("recv grid data from %d", proc));

      this->log_.info(TlUtils::format("recv block mat (%d, %d) from %d",
                                      tmpGrdMat.getNumOfRows(),
                                      tmpGrdMat.getNumOfCols(), proc));
      numOfRowsOfGlobalGridMatrix += tmpGrdMat.getNumOfRows();
      grdMat.resize(numOfRowsOfGlobalGridMatrix, numOfColsOfGlobalGridMatrix);

      grdMat.setBlockMatrix(currentNumOfRows, 0, tmpGrdMat);
      currentNumOfRows += tmpGrdMat.getNumOfRows();
    }

    // grdMat.save(this->getGridMatrixPath(0));
    this->saveGridMatrix(0, grdMat);
  } else {
    rComm.sendData(this->grdMat_, 0, tag);
    this->log_.debug("send grid data");
  }
}
