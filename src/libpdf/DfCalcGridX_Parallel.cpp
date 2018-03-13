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

#ifdef _OPENMP
#include <omp.h>
#endif  // _OPENMP

#include <unistd.h>
#include <cassert>
#include <map>
#include <queue>
#include <vector>

#include "DfCalcGridX_Parallel.h"
#include "TlCommunicate.h"
#include "TlTime.h"
#include "TlUtils.h"
#include "tl_dense_general_matrix_blacs.h"
#include "tl_dense_general_matrix_io.h"
#include "tl_dense_symmetric_matrix_blacs.h"
#include "tl_sparse_symmetric_matrix.h"

//#define USE_FILE_MATRIX
#define JOB_PROTOCOL_SIZE (4)

DfCalcGridX_Parallel::DfCalcGridX_Parallel(TlSerializeData* pPdfParam)
    : DfCalcGridX(pPdfParam) {
  const TlSerializeData& pdfParam = *pPdfParam;

  this->assignAtomRange_ = 1;
  if (pdfParam["xc-ms-atom-range"].getStr().empty() != true) {
    this->assignAtomRange_ = pdfParam["xc-ms-atom-range"].getInt();
  }

  this->assignAoRange_ = 64;
  if (pdfParam["xc_ms_ao_range"].getStr().empty() != true) {
    this->assignAoRange_ = pdfParam["xc_ms_ao_range"].getInt();
  }

  this->assignJobsPerProc_ = 2;
  if (pdfParam["xc_ms_job_per_proc"].getStr().empty() != true) {
    this->assignJobsPerProc_ = pdfParam["xc_ms_job_per_proc"].getInt();
  }

  this->densityMatrixCacheMemSize_ = 100 * 1024UL * 1024UL;  // 100MB

  this->calcMode_ = pdfParam["grid_calcmode"].getInt();
}

DfCalcGridX_Parallel::~DfCalcGridX_Parallel() {}

void DfCalcGridX_Parallel::defineCutOffValues(
    const TlDenseSymmetricMatrix_BLAS_Old& P) {
  // TODO: 並列化
  DfCalcGridX::defineCutOffValues(P);
}

void DfCalcGridX_Parallel::defineCutOffValues(
    const TlDenseSymmetricMatrix_BLAS_Old& PA,
    const TlDenseSymmetricMatrix_BLAS_Old& PB) {
  // TODO: 並列化
  DfCalcGridX::defineCutOffValues(PA, PB);
}

void DfCalcGridX_Parallel::defineCutOffValues(
    const TlDenseSymmetricMatrix_blacs& P) {
  const double maxValueOfP = std::max(P.getMaxAbsoluteElement(), 1.0E-16);
  if (maxValueOfP < 1.0) {
    this->m_densityCutOffValueA /= maxValueOfP;
  }
  this->log_.info(TlUtils::format("density cutoff value = %e",
                                  this->m_densityCutOffValueA));

  TlCommunicate& rComm = TlCommunicate::getInstance();
  rComm.broadcast(this->m_densityCutOffValueA);
}

void DfCalcGridX_Parallel::defineCutOffValues(
    const TlDenseSymmetricMatrix_blacs& PA,
    const TlDenseSymmetricMatrix_blacs& PB) {
  const double maxValueOfPA = std::max(PA.getMaxAbsoluteElement(), 1.0E-16);
  const double maxValueOfPB = std::max(PB.getMaxAbsoluteElement(), 1.0E-16);
  if (maxValueOfPA < 1.0) {
    this->m_densityCutOffValueA /= maxValueOfPA;
  }
  if (maxValueOfPB < 1.0) {
    this->m_densityCutOffValueB /= maxValueOfPB;
  }
  this->log_.info(TlUtils::format(" density cutoff value(alpha) = %e",
                                  this->m_densityCutOffValueA));
  this->log_.info(TlUtils::format(" density cutoff value(beta ) = %e",
                                  this->m_densityCutOffValueB));

  TlCommunicate& rComm = TlCommunicate::getInstance();
  rComm.broadcast(this->m_densityCutOffValueA);
  rComm.broadcast(this->m_densityCutOffValueB);
}

double DfCalcGridX_Parallel::calcXCIntegForFockAndEnergy(
    const TlDenseSymmetricMatrix_BLAS_Old& P_A, DfFunctional_LDA* pFunctional,
    TlDenseSymmetricMatrix_BLAS_Old* pF_A) {
  this->calcRho_LDA(P_A);
  double energy = this->buildVxc(pFunctional, pF_A);
  return energy;
}

double DfCalcGridX_Parallel::calcXCIntegForFockAndEnergy(
    const TlDenseSymmetricMatrix_BLAS_Old& P_A,
    const TlDenseSymmetricMatrix_BLAS_Old& P_B, DfFunctional_LDA* pFunctional,
    TlDenseSymmetricMatrix_BLAS_Old* pF_A, TlDenseSymmetricMatrix_BLAS_Old* pF_B) {
  this->calcRho_LDA(P_A, P_B);
  double energy = this->buildVxc(pFunctional, pF_A, pF_B);
  return energy;
}

double DfCalcGridX_Parallel::calcXCIntegForFockAndEnergy(
    const TlDenseSymmetricMatrix_BLAS_Old& P_A, DfFunctional_GGA* pFunctional,
    TlDenseSymmetricMatrix_BLAS_Old* pF_A) {
  this->calcRho_GGA(P_A);
  double energy = this->buildVxc(pFunctional, pF_A);
  return energy;
}

double DfCalcGridX_Parallel::calcXCIntegForFockAndEnergy(
    const TlDenseSymmetricMatrix_BLAS_Old& P_A,
    const TlDenseSymmetricMatrix_BLAS_Old& P_B, DfFunctional_GGA* pFunctional,
    TlDenseSymmetricMatrix_BLAS_Old* pF_A, TlDenseSymmetricMatrix_BLAS_Old* pF_B) {
  this->calcRho_GGA(P_A, P_B);
  double energy = this->buildVxc(pFunctional, pF_A, pF_B);
  return energy;
}

void DfCalcGridX_Parallel::calcRho_LDA(const TlDenseSymmetricMatrix_BLAS_Old& P_A) {
  TlDenseGeneralMatrix_BLAS_old gridMat =
      this->distributeGridMatrix(this->m_nIteration - 1);

  this->calcRho_LDA_part(P_A, &gridMat);

  this->gatherAndSaveGridMatrix(gridMat);
}

void DfCalcGridX_Parallel::calcRho_LDA(const TlDenseSymmetricMatrix_BLAS_Old& P_A,
                                       const TlDenseSymmetricMatrix_BLAS_Old& P_B) {
  TlDenseGeneralMatrix_BLAS_old gridMat =
      this->distributeGridMatrix(this->m_nIteration - 1);

  this->calcRho_LDA_part(P_A, P_B, &gridMat);

  this->gatherAndSaveGridMatrix(gridMat);
}

void DfCalcGridX_Parallel::calcRho_GGA(const TlDenseSymmetricMatrix_BLAS_Old& P_A) {
  TlDenseGeneralMatrix_BLAS_old gridMat =
      this->distributeGridMatrix(this->m_nIteration - 1);

  this->calcRho_GGA_part(P_A, &gridMat);

  this->gatherAndSaveGridMatrix(gridMat);
}

void DfCalcGridX_Parallel::calcRho_GGA(const TlDenseSymmetricMatrix_BLAS_Old& P_A,
                                       const TlDenseSymmetricMatrix_BLAS_Old& P_B) {
  TlDenseGeneralMatrix_BLAS_old gridMat =
      this->distributeGridMatrix(this->m_nIteration - 1);

  this->calcRho_GGA_part(P_A, P_B, &gridMat);

  this->gatherAndSaveGridMatrix(gridMat);
}

double DfCalcGridX_Parallel::buildVxc(DfFunctional_LDA* pFunctional,
                                      TlDenseSymmetricMatrix_BLAS_Old* pF_A) {
  TlDenseGeneralMatrix_BLAS_old gridMat =
      this->distributeGridMatrix(this->m_nIteration);
  double energy = DfCalcGridX::buildVxc(gridMat, pFunctional, pF_A);

  TlCommunicate& rComm = TlCommunicate::getInstance();
  rComm.allReduce_SUM(*pF_A);
  rComm.allReduce_SUM(energy);
  return energy;
}

double DfCalcGridX_Parallel::buildVxc(DfFunctional_LDA* pFunctional,
                                      TlDenseSymmetricMatrix_BLAS_Old* pF_A,
                                      TlDenseSymmetricMatrix_BLAS_Old* pF_B) {
  TlDenseGeneralMatrix_BLAS_old gridMat =
      this->distributeGridMatrix(this->m_nIteration);
  double energy = DfCalcGridX::buildVxc(gridMat, pFunctional, pF_A, pF_B);

  TlCommunicate& rComm = TlCommunicate::getInstance();
  rComm.allReduce_SUM(*pF_A);
  rComm.allReduce_SUM(energy);
  return energy;
}

double DfCalcGridX_Parallel::buildVxc(DfFunctional_GGA* pFunctional,
                                      TlDenseSymmetricMatrix_BLAS_Old* pF_A) {
  TlDenseGeneralMatrix_BLAS_old gridMat =
      this->distributeGridMatrix(this->m_nIteration);
  double energy = DfCalcGridX::buildVxc(gridMat, pFunctional, pF_A);

  TlCommunicate& rComm = TlCommunicate::getInstance();
  rComm.allReduce_SUM(*pF_A);
  rComm.allReduce_SUM(energy);
  return energy;
}

double DfCalcGridX_Parallel::buildVxc(DfFunctional_GGA* pFunctional,
                                      TlDenseSymmetricMatrix_BLAS_Old* pF_A,
                                      TlDenseSymmetricMatrix_BLAS_Old* pF_B) {
  TlDenseGeneralMatrix_BLAS_old gridMat =
      this->distributeGridMatrix(this->m_nIteration);
  double energy = DfCalcGridX::buildVxc(gridMat, pFunctional, pF_A, pF_B);

  TlCommunicate& rComm = TlCommunicate::getInstance();
  rComm.allReduce_SUM(*pF_A);
  rComm.allReduce_SUM(energy);
  return energy;
}

TlDenseGeneralMatrix_BLAS_old DfCalcGridX_Parallel::distributeGridMatrix(
    const int iteration) {
  this->log_.info("distribute grid matrix: start");

  TlCommunicate& rComm = TlCommunicate::getInstance();
  const int numOfProcs = rComm.getNumOfProc();
  // const int rank = rComm.getRank();

  const int tag = TAG_CALC_GRID_DISTRIBUTE;
  TlDenseGeneralMatrix_BLAS_old gridMat;
  if (rComm.isMaster() == true) {
    // #ifdef USE_FILE_MATRIX
    //     TlMatrixFile
    //     globalGridMat(DfObject::getGridMatrixPath(iteration));
    // #else
    TlDenseGeneralMatrix_BLAS_old globalGridMat;
    globalGridMat.load(DfObject::getGridMatrixPath(iteration));
    // #endif  // USE_FILE_MATRIX
    this->numOfRows_gridMatrix_ = globalGridMat.getNumOfRows();
    this->numOfCols_gridMatrix_ = globalGridMat.getNumOfCols();

    const index_type range =
        (this->numOfRows_gridMatrix_ + numOfProcs - 1) / numOfProcs;

    // for self(master)
    {
      const index_type startGrid = 0;
      const index_type endGrid =
          std::min<index_type>(startGrid + range, this->numOfRows_gridMatrix_);
      gridMat = globalGridMat.getBlockMatrix(startGrid, 0, endGrid - startGrid,
                                             globalGridMat.getNumOfCols());
      this->log_.debug(TlUtils::format("send grid data to 0 (%d, %d)",
                                       gridMat.getNumOfRows(),
                                       gridMat.getNumOfCols()));
    }

    // for slave
    for (int i = 1; i < numOfProcs; ++i) {
      const index_type startGrid = range * i;
      const index_type endGrid =
          std::min<index_type>(startGrid + range, this->numOfRows_gridMatrix_);
      TlDenseGeneralMatrix_BLAS_old tmpMat = globalGridMat.getBlockMatrix(
          startGrid, 0, endGrid - startGrid, globalGridMat.getNumOfCols());
      this->log_.debug(TlUtils::format("send grid data to %d (%d, %d)", i,
                                       tmpMat.getNumOfRows(),
                                       tmpMat.getNumOfCols()));
      rComm.sendData(tmpMat, i, tag);
    }
  } else {
    rComm.receiveData(gridMat, 0, tag);
    this->log_.debug("recv grid data");
  }

  this->log_.info("distribute grid matrix: end");
  return gridMat;
}

void DfCalcGridX_Parallel::gatherAndSaveGridMatrix(
    const TlDenseGeneralMatrix_BLAS_old& gridMat) {
  this->log_.info("gather grid matrix: start");

  TlCommunicate& rComm = TlCommunicate::getInstance();
  const int numOfProcs = rComm.getNumOfProc();
  const int tag = TAG_CALCGRID_GATHER;

  if (rComm.isMaster() == true) {
    // #ifdef USE_FILE_MATRIX
    //     TlMatrixFile globalGridMat(
    //         DfObject::getGridMatrixPath(this->m_nIteration),
    //         this->numOfRows_gridMatrix_, this->numOfCols_gridMatrix_);
    // #else
    TlDenseGeneralMatrix_BLAS_old globalGridMat(this->numOfRows_gridMatrix_,
                                            this->numOfCols_gridMatrix_);
    // #endif  // USE_FILE_MATRIX
    this->log_.debug(TlUtils::format("recv grid data from 0 (%d, %d)",
                                     gridMat.getNumOfRows(),
                                     gridMat.getNumOfCols()));
    globalGridMat.setBlockMatrix(0, 0, gridMat);
    index_type currentGridIndex = gridMat.getNumOfRows();
    this->log_.debug(TlUtils::format("currentGridIndex=%d", currentGridIndex));

    std::vector<bool> recvCheck(numOfProcs, false);
    for (int i = 1; i < numOfProcs; ++i) {
      int proc = 0;
      TlDenseGeneralMatrix_BLAS_old tmpMat;
      rComm.receiveDataFromAnySource(tmpMat, &proc, tag);
      if (recvCheck[proc] != false) {
        this->log_.warn(
            TlUtils::format("already receive grid data from %d", proc));
      }
      recvCheck[proc] = true;
      this->log_.debug(TlUtils::format("recv grid data from %d (%d, %d)", proc,
                                       tmpMat.getNumOfRows(),
                                       tmpMat.getNumOfCols()));

      assert(globalGridMat.getNumOfCols() == gridMat.getNumOfCols());
      globalGridMat.setBlockMatrix(currentGridIndex, 0, tmpMat);
      currentGridIndex += tmpMat.getNumOfRows();
      this->log_.debug(
          TlUtils::format("currentGridIndex=%d", currentGridIndex));
    }

    // #ifndef USE_FILE_MATRIX
    globalGridMat.save(DfObject::getGridMatrixPath(this->m_nIteration));
    // #endif  // USE_FILE_MATRIX
  } else {
    rComm.sendData(gridMat, 0, tag);
  }
  this->log_.info("gather grid matrix: end");
}

double DfCalcGridX_Parallel::calcXCIntegForFockAndEnergy(
    const TlDenseSymmetricMatrix_blacs& P_A, DfFunctional_LDA* pFunctional,
    TlDenseSymmetricMatrix_blacs* pF_A) {
  this->calcRho_LDA(P_A);
  double energy = this->buildVxc(pFunctional, pF_A);
  return energy;
}

double DfCalcGridX_Parallel::calcXCIntegForFockAndEnergy(
    const TlDenseSymmetricMatrix_blacs& P_A,
    const TlDenseSymmetricMatrix_blacs& P_B, DfFunctional_LDA* pFunctional,
    TlDenseSymmetricMatrix_blacs* pF_A, TlDenseSymmetricMatrix_blacs* pF_B) {
  this->calcRho_LDA(P_A, P_B);
  double energy = this->buildVxc(pFunctional, pF_A, pF_B);
  return energy;
}

double DfCalcGridX_Parallel::calcXCIntegForFockAndEnergy(
    const TlDenseSymmetricMatrix_blacs& P_A, DfFunctional_GGA* pFunctional,
    TlDenseSymmetricMatrix_blacs* pF_A) {
  this->calcRho_GGA(P_A);
  double energy = this->buildVxc(pFunctional, pF_A);
  return energy;
}

double DfCalcGridX_Parallel::calcXCIntegForFockAndEnergy(
    const TlDenseSymmetricMatrix_blacs& P_A,
    const TlDenseSymmetricMatrix_blacs& P_B, DfFunctional_GGA* pFunctional,
    TlDenseSymmetricMatrix_blacs* pF_A, TlDenseSymmetricMatrix_blacs* pF_B) {
  this->calcRho_GGA(P_A, P_B);
  double energy = this->buildVxc(pFunctional, pF_A, pF_B);
  return energy;
}

TlDenseGeneralMatrix_BLAS_old DfCalcGridX_Parallel::getGlobalGridMatrix(
    const int iteration) {
  TlCommunicate& rComm = TlCommunicate::getInstance();
  TlDenseGeneralMatrix_BLAS_old gridMat;
  TlDenseGeneralMatrix_BLAS_old rhoMat;
  index_type numOfGrids = 0;
  index_type numOfCols = 0;
  if (rComm.isMaster() == true) {
    gridMat = DfObject::getGridMatrix<TlDenseGeneralMatrix_BLAS_old>(iteration);
    numOfGrids = gridMat.getNumOfRows();
    numOfCols = gridMat.getNumOfCols();

    TlDenseGeneralMatrix_BLAS_old crdMat =
        gridMat.getBlockMatrix(0, 0, numOfGrids, GM_ATOM_INDEX + 1);
    rhoMat = gridMat.getBlockMatrix(0, GM_ATOM_INDEX + 1, numOfGrids,
                                    numOfCols - (GM_ATOM_INDEX + 1));
    gridMat = crdMat;
  }
  rComm.broadcast(numOfGrids);
  rComm.broadcast(numOfCols);
  rComm.broadcast(gridMat);
  gridMat.resize(numOfGrids, numOfCols);

  if (this->m_bIsUpdateXC == true) {
    if (rComm.isMaster() == true) {
      gridMat.setBlockMatrix(0, GM_ATOM_INDEX + 1, rhoMat);
    }
  }

  return gridMat;
}

void DfCalcGridX_Parallel::allReduceGridMatrix(
    const TlDenseGeneralMatrix_BLAS_old& gridMat) {
  TlCommunicate& rComm = TlCommunicate::getInstance();
  TlDenseGeneralMatrix_BLAS_old rhoMat =
      gridMat.getBlockMatrix(0, GM_LDA_RHO_ALPHA, gridMat.getNumOfRows(),
                             gridMat.getNumOfCols() - GM_LDA_RHO_ALPHA);
  rComm.allReduce_SUM(rhoMat);

  if (rComm.isMaster() == true) {
    TlFileGenericMatrix globalGridMat(
        DfObject::getGridMatrixPath(this->m_nIteration), gridMat.getNumOfRows(),
        gridMat.getNumOfCols());
    assert(globalGridMat.getNumOfRows() == gridMat.getNumOfRows());
    assert(globalGridMat.getNumOfCols() == gridMat.getNumOfCols());

    const TlDenseGeneralMatrix_BLAS_old crdMat =
        gridMat.getBlockMatrix(0, 0, gridMat.getNumOfRows(), GM_LDA_RHO_ALPHA);
    globalGridMat.setBlockMatrix(0, 0, crdMat);
    globalGridMat.setBlockMatrix(0, GM_LDA_RHO_ALPHA, rhoMat);
  }
}

void DfCalcGridX_Parallel::calcRho_LDA(
    const TlDenseSymmetricMatrix_blacs& P_A) {
  TlDenseGeneralMatrix_BLAS_old gridMat =
      this->getGlobalGridMatrix(this->m_nIteration - 1);
  this->calcRho_LDA(TlDenseGeneralMatrix_blacs(P_A), &gridMat);
  this->allReduceGridMatrix(gridMat);
}

void DfCalcGridX_Parallel::calcRho_LDA(
    const TlDenseSymmetricMatrix_blacs& P_A,
    const TlDenseSymmetricMatrix_blacs& P_B) {
  TlDenseGeneralMatrix_BLAS_old gridMat =
      this->getGlobalGridMatrix(this->m_nIteration - 1);
  this->calcRho_LDA(TlDenseGeneralMatrix_blacs(P_A),
                    TlDenseGeneralMatrix_blacs(P_B), &gridMat);
  this->allReduceGridMatrix(gridMat);
}

void DfCalcGridX_Parallel::calcRho_GGA(
    const TlDenseSymmetricMatrix_blacs& P_A) {
  TlDenseGeneralMatrix_BLAS_old gridMat =
      this->getGlobalGridMatrix(this->m_nIteration - 1);
  this->calcRho_GGA(TlDenseGeneralMatrix_blacs(P_A), &gridMat);
  this->allReduceGridMatrix(gridMat);
}

void DfCalcGridX_Parallel::calcRho_GGA(
    const TlDenseSymmetricMatrix_blacs& P_A,
    const TlDenseSymmetricMatrix_blacs& P_B) {
  TlDenseGeneralMatrix_BLAS_old gridMat =
      this->getGlobalGridMatrix(this->m_nIteration - 1);
  this->calcRho_GGA(TlDenseGeneralMatrix_blacs(P_A),
                    TlDenseGeneralMatrix_blacs(P_B), &gridMat);
  this->allReduceGridMatrix(gridMat);
}

void DfCalcGridX_Parallel::calcRho_LDA(const TlDenseGeneralMatrix_blacs& P_A,
                                       TlDenseGeneralMatrix_BLAS_old* pGridMat) {
  const TlDenseGeneralMatrix_BLAS_old localP_A = P_A.getLocalMatrix();
  const std::vector<index_type> rowIndeces = P_A.getRowIndexTable();
  const std::vector<index_type> colIndeces = P_A.getColIndexTable();

  const index_type numOfGrids = pGridMat->getNumOfRows();
#pragma omp parallel for schedule(runtime)
  for (index_type grid = 0; grid < numOfGrids; ++grid) {
    const TlPosition gridPosition(pGridMat->get(grid, GM_X),
                                  pGridMat->get(grid, GM_Y),
                                  pGridMat->get(grid, GM_Z));

    // calc phi table
    std::vector<double> row_AO_values;
    this->getAOs_core(gridPosition, rowIndeces, &row_AO_values);
    std::vector<double> col_AO_values;
    this->getAOs_core(gridPosition, colIndeces, &col_AO_values);

    // get rho at grid point
    double rhoA = 0.0;
    this->getRhoAtGridPoint(localP_A, row_AO_values, col_AO_values, &rhoA);

    pGridMat->add(grid, GM_LDA_RHO_ALPHA, rhoA);
  }
}

void DfCalcGridX_Parallel::calcRho_LDA(const TlDenseGeneralMatrix_blacs& P_A,
                                       const TlDenseGeneralMatrix_blacs& P_B,
                                       TlDenseGeneralMatrix_BLAS_old* pGridMat) {
  const TlDenseGeneralMatrix_BLAS_old localP_A = P_A.getLocalMatrix();
  const TlDenseGeneralMatrix_BLAS_old localP_B = P_B.getLocalMatrix();
  const std::vector<index_type> rowIndeces = P_A.getRowIndexTable();
  const std::vector<index_type> colIndeces = P_A.getColIndexTable();
  // rowIndexes, colIndexes はP_Bと共通

  const index_type numOfGrids = pGridMat->getNumOfRows();
#pragma omp parallel for schedule(runtime)
  for (index_type grid = 0; grid < numOfGrids; ++grid) {
    const TlPosition gridPosition(pGridMat->get(grid, GM_X),
                                  pGridMat->get(grid, GM_Y),
                                  pGridMat->get(grid, GM_Z));

    // calc phi table
    std::vector<double> row_AO_values;
    this->getAOs_core(gridPosition, rowIndeces, &row_AO_values);
    std::vector<double> col_AO_values;
    this->getAOs_core(gridPosition, colIndeces, &col_AO_values);

    // get rho at grid point
    double rhoA = 0.0;
    double rhoB = 0.0;
    this->getRhoAtGridPoint(localP_A, row_AO_values, col_AO_values, &rhoA);
    this->getRhoAtGridPoint(localP_B, row_AO_values, col_AO_values, &rhoB);

    pGridMat->add(grid, GM_LDA_RHO_ALPHA, rhoA);
    pGridMat->add(grid, GM_LDA_RHO_BETA, rhoB);
  }
}

void DfCalcGridX_Parallel::calcRho_GGA(const TlDenseGeneralMatrix_blacs& P_A,
                                       TlDenseGeneralMatrix_BLAS_old* pGridMat) {
  const TlDenseGeneralMatrix_BLAS_old localP_A = P_A.getLocalMatrix();
  const std::vector<index_type> rowIndeces = P_A.getRowIndexTable();
  const std::vector<index_type> colIndeces = P_A.getColIndexTable();

  const index_type numOfGrids = pGridMat->getNumOfRows();
#pragma omp parallel for schedule(runtime)
  for (index_type grid = 0; grid < numOfGrids; ++grid) {
    const TlPosition gridPosition(pGridMat->get(grid, GM_X),
                                  pGridMat->get(grid, GM_Y),
                                  pGridMat->get(grid, GM_Z));

    // calc phi table
    std::vector<double> row_AO_values;
    this->getAOs_core(gridPosition, rowIndeces, &row_AO_values);
    std::vector<double> col_AO_values;
    this->getAOs_core(gridPosition, colIndeces, &col_AO_values);

    std::vector<double> row_dAO_dx_values;
    std::vector<double> row_dAO_dy_values;
    std::vector<double> row_dAO_dz_values;
    this->getDAOs_core(gridPosition, rowIndeces, &row_dAO_dx_values,
                       &row_dAO_dy_values, &row_dAO_dz_values);
    std::vector<double> col_dAO_dx_values;
    std::vector<double> col_dAO_dy_values;
    std::vector<double> col_dAO_dz_values;
    this->getDAOs_core(gridPosition, colIndeces, &col_dAO_dx_values,
                       &col_dAO_dy_values, &col_dAO_dz_values);

    // get rho at grid point
    double rhoA = 0.0;
    double gradRhoAX = 0.0;
    double gradRhoAY = 0.0;
    double gradRhoAZ = 0.0;
    this->getRhoAtGridPoint(localP_A, row_AO_values, col_AO_values, &rhoA);
    this->getGradRhoAtGridPoint(
        localP_A, row_AO_values, row_dAO_dx_values, row_dAO_dy_values,
        row_dAO_dz_values, col_AO_values, col_dAO_dx_values, col_dAO_dy_values,
        col_dAO_dz_values, &gradRhoAX, &gradRhoAY, &gradRhoAZ);

    pGridMat->add(grid, GM_GGA_RHO_ALPHA, rhoA);
    pGridMat->add(grid, GM_GGA_GRAD_RHO_X_ALPHA, gradRhoAX);
    pGridMat->add(grid, GM_GGA_GRAD_RHO_Y_ALPHA, gradRhoAY);
    pGridMat->add(grid, GM_GGA_GRAD_RHO_Z_ALPHA, gradRhoAZ);
  }
}

void DfCalcGridX_Parallel::calcRho_GGA(const TlDenseGeneralMatrix_blacs& P_A,
                                       const TlDenseGeneralMatrix_blacs& P_B,
                                       TlDenseGeneralMatrix_BLAS_old* pGridMat) {
  const TlDenseGeneralMatrix_BLAS_old localP_A = P_A.getLocalMatrix();
  const TlDenseGeneralMatrix_BLAS_old localP_B = P_B.getLocalMatrix();
  const std::vector<index_type> rowIndeces = P_A.getRowIndexTable();
  const std::vector<index_type> colIndeces = P_A.getColIndexTable();
  // rowIndexes, colIndexes はP_Bと共通

  const index_type numOfGrids = pGridMat->getNumOfRows();
#pragma omp parallel for schedule(runtime)
  for (index_type grid = 0; grid < numOfGrids; ++grid) {
    const TlPosition gridPosition(pGridMat->get(grid, GM_X),
                                  pGridMat->get(grid, GM_Y),
                                  pGridMat->get(grid, GM_Z));

    // calc phi table
    std::vector<double> row_AO_values;
    this->getAOs_core(gridPosition, rowIndeces, &row_AO_values);
    std::vector<double> col_AO_values;
    this->getAOs_core(gridPosition, colIndeces, &col_AO_values);

    std::vector<double> row_dAO_dx_values;
    std::vector<double> row_dAO_dy_values;
    std::vector<double> row_dAO_dz_values;
    this->getDAOs_core(gridPosition, rowIndeces, &row_dAO_dx_values,
                       &row_dAO_dy_values, &row_dAO_dz_values);
    std::vector<double> col_dAO_dx_values;
    std::vector<double> col_dAO_dy_values;
    std::vector<double> col_dAO_dz_values;
    this->getDAOs_core(gridPosition, colIndeces, &col_dAO_dx_values,
                       &col_dAO_dy_values, &col_dAO_dz_values);

    // get rho at grid point
    double rhoA = 0.0;
    double gradRhoAX = 0.0;
    double gradRhoAY = 0.0;
    double gradRhoAZ = 0.0;
    double rhoB = 0.0;
    double gradRhoBX = 0.0;
    double gradRhoBY = 0.0;
    double gradRhoBZ = 0.0;
    this->getRhoAtGridPoint(localP_A, row_AO_values, col_AO_values, &rhoA);
    this->getGradRhoAtGridPoint(
        localP_A, row_AO_values, row_dAO_dx_values, row_dAO_dy_values,
        row_dAO_dz_values, col_AO_values, col_dAO_dx_values, col_dAO_dy_values,
        col_dAO_dz_values, &gradRhoAX, &gradRhoAY, &gradRhoAZ);
    this->getRhoAtGridPoint(localP_B, row_AO_values, col_AO_values, &rhoB);
    this->getGradRhoAtGridPoint(
        localP_B, row_AO_values, row_dAO_dx_values, row_dAO_dy_values,
        row_dAO_dz_values, col_AO_values, col_dAO_dx_values, col_dAO_dy_values,
        col_dAO_dz_values, &gradRhoBX, &gradRhoBY, &gradRhoBZ);

    pGridMat->add(grid, GM_GGA_RHO_ALPHA, rhoA);
    pGridMat->add(grid, GM_GGA_GRAD_RHO_X_ALPHA, gradRhoAX);
    pGridMat->add(grid, GM_GGA_GRAD_RHO_Y_ALPHA, gradRhoAY);
    pGridMat->add(grid, GM_GGA_GRAD_RHO_Z_ALPHA, gradRhoAZ);
    pGridMat->add(grid, GM_GGA_RHO_BETA, rhoB);
    pGridMat->add(grid, GM_GGA_GRAD_RHO_X_BETA, gradRhoBX);
    pGridMat->add(grid, GM_GGA_GRAD_RHO_Y_BETA, gradRhoBY);
    pGridMat->add(grid, GM_GGA_GRAD_RHO_Z_BETA, gradRhoBZ);
  }
}

double DfCalcGridX_Parallel::buildVxc(DfFunctional_LDA* pFunctional,
                                      TlDenseSymmetricMatrix_blacs* pF_A) {
  const TlDenseGeneralMatrix_BLAS_old gridMat =
      this->distributeGridMatrix(this->m_nIteration);
  TlSparseSymmetricMatrix tmpF_A(this->m_nNumOfAOs);
  double energy = DfCalcGridX::buildVxc(gridMat, pFunctional, &tmpF_A);

  TlCommunicate& rComm = TlCommunicate::getInstance();
  pF_A->mergeSparseMatrix(tmpF_A);
  rComm.allReduce_SUM(energy);
  return energy;
}

double DfCalcGridX_Parallel::buildVxc(DfFunctional_LDA* pFunctional,
                                      TlDenseSymmetricMatrix_blacs* pF_A,
                                      TlDenseSymmetricMatrix_blacs* pF_B) {
  const TlDenseGeneralMatrix_BLAS_old gridMat =
      this->distributeGridMatrix(this->m_nIteration);
  TlSparseSymmetricMatrix tmpF_A(this->m_nNumOfAOs);
  TlSparseSymmetricMatrix tmpF_B(this->m_nNumOfAOs);
  double energy = DfCalcGridX::buildVxc(gridMat, pFunctional, &tmpF_A, &tmpF_B);

  TlCommunicate& rComm = TlCommunicate::getInstance();
  pF_A->mergeSparseMatrix(tmpF_A);
  pF_B->mergeSparseMatrix(tmpF_B);
  rComm.allReduce_SUM(energy);
  return energy;
}

double DfCalcGridX_Parallel::buildVxc(DfFunctional_GGA* pFunctional,
                                      TlDenseSymmetricMatrix_blacs* pF_A) {
  const TlDenseGeneralMatrix_BLAS_old gridMat =
      this->distributeGridMatrix(this->m_nIteration);
  TlSparseSymmetricMatrix tmpF_A(this->m_nNumOfAOs);
  double energy = DfCalcGridX::buildVxc(gridMat, pFunctional, &tmpF_A);

  TlCommunicate& rComm = TlCommunicate::getInstance();
  pF_A->mergeSparseMatrix(tmpF_A);
  rComm.allReduce_SUM(energy);
  return energy;
}

double DfCalcGridX_Parallel::buildVxc(DfFunctional_GGA* pFunctional,
                                      TlDenseSymmetricMatrix_blacs* pF_A,
                                      TlDenseSymmetricMatrix_blacs* pF_B) {
  const TlDenseGeneralMatrix_BLAS_old gridMat =
      this->distributeGridMatrix(this->m_nIteration);
  TlSparseSymmetricMatrix tmpF_A(this->m_nNumOfAOs);
  TlSparseSymmetricMatrix tmpF_B(this->m_nNumOfAOs);
  double energy = DfCalcGridX::buildVxc(gridMat, pFunctional, &tmpF_A, &tmpF_B);

  TlCommunicate& rComm = TlCommunicate::getInstance();
  pF_A->mergeSparseMatrix(tmpF_A);
  pF_B->mergeSparseMatrix(tmpF_B);
  rComm.allReduce_SUM(energy);
  return energy;
}

void DfCalcGridX_Parallel::getWholeDensity(double* pRhoA, double* pRhoB) const {
  TlCommunicate& rComm = TlCommunicate::getInstance();
  if (rComm.isMaster() == true) {
    DfCalcGridX::getWholeDensity(pRhoA, pRhoB);
  }
}

TlDenseGeneralMatrix_BLAS_old DfCalcGridX_Parallel::energyGradient(
    const TlDenseSymmetricMatrix_BLAS_Old& P_A, DfFunctional_LDA* pFunctional) {
  TlCommunicate& rComm = TlCommunicate::getInstance();
  this->log_.info(
      "pure DFT XC energy (LDA) gradient by grid method (parallel)");
  const int numOfAOs = this->m_nNumOfAOs;
  const int numOfAtoms = this->m_nNumOfAtoms;

  TlDenseGeneralMatrix_BLAS_old gammaX(numOfAOs, numOfAOs);
  TlDenseGeneralMatrix_BLAS_old gammaY(numOfAOs, numOfAOs);
  TlDenseGeneralMatrix_BLAS_old gammaZ(numOfAOs, numOfAOs);

  TlDenseGeneralMatrix_BLAS_old Fxc_f(numOfAtoms, 3);  // func derivative term
  TlDenseGeneralMatrix_BLAS_old Fxc_w(numOfAtoms, 3);  // weight derivative term

  const int myRank = rComm.getRank();
  const int numOfProcs = rComm.getNumOfProcs();
  const int range = (numOfAtoms + numOfProcs) / numOfProcs;
  const int startIndex = range * myRank;
  const int endIndex = std::min(range * (myRank + 1), numOfAtoms);
  double ene_xc = this->energyGradient_part(P_A, pFunctional, startIndex,
                                            endIndex, &Fxc_f, &Fxc_w);
  rComm.allReduce_SUM(ene_xc);
  rComm.allReduce_SUM(Fxc_f);
  rComm.allReduce_SUM(Fxc_w);

  this->log_.info(TlUtils::format("XC ene = % 16.10f", ene_xc));

  Fxc_f *= 2.0;  // rks

  // debug
  if (this->isSaveGrad_) {
    if (rComm.isMaster()) {
      Fxc_w.save("Fxc_w.mat");
      Fxc_f.save("Fxc_f.mat");
    }
  }

  return Fxc_w + Fxc_f;
}

TlDenseGeneralMatrix_BLAS_old DfCalcGridX_Parallel::energyGradient(
    const TlDenseSymmetricMatrix_BLAS_Old& P_A, DfFunctional_GGA* pFunctional) {
  TlCommunicate& rComm = TlCommunicate::getInstance();
  this->log_.info(
      "pure DFT XC energy (LDA) gradient by grid method (parallel)");
  const int numOfAOs = this->m_nNumOfAOs;
  const int numOfAtoms = this->m_nNumOfAtoms;

  TlDenseGeneralMatrix_BLAS_old gammaX(numOfAOs, numOfAOs);
  TlDenseGeneralMatrix_BLAS_old gammaY(numOfAOs, numOfAOs);
  TlDenseGeneralMatrix_BLAS_old gammaZ(numOfAOs, numOfAOs);

  TlDenseGeneralMatrix_BLAS_old Fxc_f(numOfAtoms, 3);  // func derivative term
  TlDenseGeneralMatrix_BLAS_old Fxc_w(numOfAtoms, 3);  // weight derivative term

  const int myRank = rComm.getRank();
  const int numOfProcs = rComm.getNumOfProcs();
  const int range = (numOfAtoms + numOfProcs) / numOfProcs;
  const int startIndex = range * myRank;
  const int endIndex = std::min(range * (myRank + 1), numOfAtoms);
  double ene_xc = this->energyGradient_part(P_A, pFunctional, startIndex,
                                            endIndex, &Fxc_f, &Fxc_w);
  rComm.allReduce_SUM(ene_xc);
  rComm.allReduce_SUM(Fxc_f);
  rComm.allReduce_SUM(Fxc_w);

  this->log_.info(TlUtils::format("XC ene = % 16.10f", ene_xc));

  Fxc_f *= 2.0;  // rks

  // debug
  if (this->isSaveGrad_) {
    if (rComm.isMaster()) {
      Fxc_w.save("Fxc_w.mat");
      Fxc_f.save("Fxc_f.mat");
    }
  }

  return Fxc_w + Fxc_f;
}
