#ifdef _OPENMP
#include <omp.h>
#endif // _OPENMP

#include <unistd.h>
#include <cassert>
#include <vector>
#include <map>
#include <queue>
#include <cassert>

#include "DfCalcGridX_Parallel.h"
#include "TlCommunicate.h"
#include "TlTime.h"
#include "TlUtils.h"
#include "TlSparseSymmetricMatrix.h"
#include "TlFileMatrix.h"

//#define USE_FILE_MATRIX
#define JOB_PROTOCOL_SIZE (4)

DfCalcGridX_Parallel::DfCalcGridX_Parallel(TlSerializeData* pPdfParam)
    : DfCalcGridX(pPdfParam)
{
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
    
    this->densityMatrixCacheMemSize_ = 100 * 1024UL * 1024UL; // 100MB

    this->calcMode_ = pdfParam["grid_calcmode"].getInt();
}

DfCalcGridX_Parallel::~DfCalcGridX_Parallel()
{
}


void DfCalcGridX_Parallel::defineCutOffValues(const TlSymmetricMatrix& P)
{
    // TODO: 並列化
    DfCalcGridX::defineCutOffValues(P);
}


void DfCalcGridX_Parallel::defineCutOffValues(const TlSymmetricMatrix& PA,
                                              const TlSymmetricMatrix& PB)
{
    // TODO: 並列化
    DfCalcGridX::defineCutOffValues(PA, PB);
}


void DfCalcGridX_Parallel::defineCutOffValues(const TlDistributeSymmetricMatrix& P)
{
    const double maxValueOfP = std::max(P.getMaxAbsoluteElement(), 1.0E-16);
    if (maxValueOfP < 1.0) {
        this->m_densityCutOffValueA /= maxValueOfP;
    }
    this->log_.info(TlUtils::format("density cutoff value = %e",
                                    this->m_densityCutOffValueA));
    
    TlCommunicate& rComm = TlCommunicate::getInstance();
    rComm.broadcast(this->m_densityCutOffValueA);
}


void DfCalcGridX_Parallel::defineCutOffValues(const TlDistributeSymmetricMatrix& PA,
                                              const TlDistributeSymmetricMatrix& PB)
{
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


double DfCalcGridX_Parallel::calcXCIntegForFockAndEnergy(const TlSymmetricMatrix& P_A,
                                                         DfFunctional_LDA* pFunctional,
                                                         TlSymmetricMatrix* pF_A)
{
    this->calcRho_LDA(P_A);
    double energy = this->buildVxc(pFunctional, pF_A);
    return energy;
}

double DfCalcGridX_Parallel::calcXCIntegForFockAndEnergy(const TlSymmetricMatrix& P_A,
                                                         const TlSymmetricMatrix& P_B,
                                                         DfFunctional_LDA* pFunctional,
                                                         TlSymmetricMatrix* pF_A,
                                                         TlSymmetricMatrix* pF_B)
{
    this->calcRho_LDA(P_A, P_B);
    double energy = this->buildVxc(pFunctional, pF_A, pF_B);
    return energy;
}

double DfCalcGridX_Parallel::calcXCIntegForFockAndEnergy(const TlSymmetricMatrix& P_A,
                                                         DfFunctional_GGA* pFunctional,
                                                         TlSymmetricMatrix* pF_A)
{
    this->calcRho_GGA(P_A);
    double energy = this->buildVxc(pFunctional, pF_A);
    return energy;
}

double DfCalcGridX_Parallel::calcXCIntegForFockAndEnergy(const TlSymmetricMatrix& P_A,
                                                         const TlSymmetricMatrix& P_B,
                                                         DfFunctional_GGA* pFunctional,
                                                         TlSymmetricMatrix* pF_A,
                                                         TlSymmetricMatrix* pF_B)
{
    this->calcRho_GGA(P_A, P_B);
    double energy = this->buildVxc(pFunctional, pF_A, pF_B);
    return energy;
}

void DfCalcGridX_Parallel::calcRho_LDA(const TlSymmetricMatrix& P_A)
{
    TlMatrix gridMat = this->distributeGridMatrix(this->m_nIteration -1);
    DfCalcGridX::calcRho_LDA(P_A, &gridMat);
    this->gatherGridMatrix(gridMat);
}

void DfCalcGridX_Parallel::calcRho_LDA(const TlSymmetricMatrix& P_A,
                                       const TlSymmetricMatrix& P_B)
{
    TlMatrix gridMat = this->distributeGridMatrix(this->m_nIteration -1);
    DfCalcGridX::calcRho_LDA(P_A, P_B, &gridMat);
    this->gatherGridMatrix(gridMat);
}

void DfCalcGridX_Parallel::calcRho_GGA(const TlSymmetricMatrix& P_A)
{
    TlMatrix gridMat = this->distributeGridMatrix(this->m_nIteration -1);
    DfCalcGridX::calcRho_GGA(P_A, &gridMat);
    this->gatherGridMatrix(gridMat);
}

void DfCalcGridX_Parallel::calcRho_GGA(const TlSymmetricMatrix& P_A,
                                       const TlSymmetricMatrix& P_B)
{
    TlMatrix gridMat = this->distributeGridMatrix(this->m_nIteration -1);
    DfCalcGridX::calcRho_GGA(P_A, P_B, &gridMat);
    this->gatherGridMatrix(gridMat);
}

double DfCalcGridX_Parallel::buildVxc(DfFunctional_LDA* pFunctional,
                                      TlSymmetricMatrix* pF_A)
{
    TlMatrix gridMat = this->distributeGridMatrix(this->m_nIteration);
    double energy = DfCalcGridX::buildVxc(gridMat, pFunctional, pF_A);

    TlCommunicate& rComm = TlCommunicate::getInstance();
    rComm.allReduce_SUM(*pF_A);
    rComm.allReduce_SUM(energy);
    return energy;
}

double DfCalcGridX_Parallel::buildVxc(DfFunctional_LDA* pFunctional,
                                      TlSymmetricMatrix* pF_A,
                                      TlSymmetricMatrix* pF_B)
{
    TlMatrix gridMat = this->distributeGridMatrix(this->m_nIteration);
    double energy = DfCalcGridX::buildVxc(gridMat, pFunctional, pF_A, pF_B);

    TlCommunicate& rComm = TlCommunicate::getInstance();
    rComm.allReduce_SUM(*pF_A);
    rComm.allReduce_SUM(energy);
    return energy;
}

double DfCalcGridX_Parallel::buildVxc(DfFunctional_GGA* pFunctional,
                                      TlSymmetricMatrix* pF_A)
{
    TlMatrix gridMat = this->distributeGridMatrix(this->m_nIteration);
    double energy = DfCalcGridX::buildVxc(gridMat, pFunctional, pF_A);

    TlCommunicate& rComm = TlCommunicate::getInstance();
    rComm.allReduce_SUM(*pF_A);
    rComm.allReduce_SUM(energy);
    return energy;
}

double DfCalcGridX_Parallel::buildVxc(DfFunctional_GGA* pFunctional,
                                      TlSymmetricMatrix* pF_A,
                                      TlSymmetricMatrix* pF_B)
{
    TlMatrix gridMat = this->distributeGridMatrix(this->m_nIteration);
    double energy = DfCalcGridX::buildVxc(gridMat, pFunctional, pF_A, pF_B);

    TlCommunicate& rComm = TlCommunicate::getInstance();
    rComm.allReduce_SUM(*pF_A);
    rComm.allReduce_SUM(energy);
    return energy;
}

TlMatrix DfCalcGridX_Parallel::distributeGridMatrix(const int iteration)
{
    this->log_.info("distribute grid matrix: start");

    TlCommunicate& rComm = TlCommunicate::getInstance();
    const int numOfProcs = rComm.getNumOfProc();
    //const int rank = rComm.getRank();

    const int tag = TAG_CALC_GRID_DISTRIBUTE;
    TlMatrix gridMat;
    if (rComm.isMaster() == true) {
#ifdef USE_FILE_MATRIX
        TlFileMatrix globalGridMat(DfObject::getGridMatrixPath(iteration));
#else
        TlMatrix globalGridMat;
        globalGridMat.load(DfObject::getGridMatrixPath(iteration));
#endif // USE_FILE_MATRIX
        this->numOfRows_gridMatrix_ = globalGridMat.getNumOfRows();
        this->numOfCols_gridMatrix_ = globalGridMat.getNumOfCols();
        
        const index_type range = (this->numOfRows_gridMatrix_ + numOfProcs -1) / numOfProcs;

        // for self(master)
        {
            const index_type startGrid = 0;
            const index_type endGrid = std::min<index_type>(startGrid + range,
                                                            this->numOfRows_gridMatrix_);
            gridMat = globalGridMat.getBlockMatrix(startGrid, 0,
                                                   endGrid - startGrid,
                                                   globalGridMat.getNumOfCols());
            this->log_.debug(TlUtils::format("send grid data to 0 (%d, %d)",
                                             gridMat.getNumOfRows(),
                                             gridMat.getNumOfCols()));
        }
        
        // for slave
        for (int i = 1; i < numOfProcs; ++i) {
            const index_type startGrid = range * i;
            const index_type endGrid = std::min<index_type>(startGrid + range,
                                                            this->numOfRows_gridMatrix_);
            TlMatrix tmpMat = globalGridMat.getBlockMatrix(startGrid, 0,
                                                           endGrid - startGrid,
                                                           globalGridMat.getNumOfCols());
            this->log_.debug(TlUtils::format("send grid data to %d (%d, %d)",
                                             i,
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

void DfCalcGridX_Parallel::gatherGridMatrix(const TlMatrix& gridMat)
{
    this->log_.info("gather grid matrix: start");
    
    TlCommunicate& rComm = TlCommunicate::getInstance();
    const int numOfProcs = rComm.getNumOfProc();
    const int tag = TAG_CALCGRID_GATHER;
    
    if (rComm.isMaster() == true) {
#ifdef USE_FILE_MATRIX
        TlFileMatrix globalGridMat(DfObject::getGridMatrixPath(this->m_nIteration),
                                   this->numOfRows_gridMatrix_,
                                   this->numOfCols_gridMatrix_);
#else
        TlMatrix globalGridMat(this->numOfRows_gridMatrix_,
                               this->numOfCols_gridMatrix_);
#endif // USE_FILE_MATRIX
        this->log_.debug(TlUtils::format("recv grid data from 0 (%d, %d)",
                                         gridMat.getNumOfRows(),
                                         gridMat.getNumOfCols()));
        globalGridMat.setBlockMatrix(0, 0,
                                     gridMat);
        index_type currentGridIndex = gridMat.getNumOfRows();
        this->log_.debug(TlUtils::format("currentGridIndex=%d",
                                         currentGridIndex));
        
        std::vector<bool> recvCheck(numOfProcs, false);
        for (int i = 1; i < numOfProcs; ++i) {
            int proc = 0;
            TlMatrix tmpMat;
            rComm.receiveDataFromAnySource(tmpMat, &proc, tag);
            if (recvCheck[proc] != false) {
                this->log_.warn(TlUtils::format("already receive grid data from %d",
                                                proc));
            }
            recvCheck[proc] = true;
            this->log_.debug(TlUtils::format("recv grid data from %d (%d, %d)",
                                             proc,
                                             tmpMat.getNumOfRows(),
                                             tmpMat.getNumOfCols()));
            
            assert(globalGridMat.getNumOfCols() == gridMat.getNumOfCols());
            globalGridMat.setBlockMatrix(currentGridIndex, 0,
                                         tmpMat);
            currentGridIndex += tmpMat.getNumOfRows();
            this->log_.debug(TlUtils::format("currentGridIndex=%d",
                                             currentGridIndex));
        }
        
#ifndef USE_FILE_MATRIX
        globalGridMat.save(DfObject::getGridMatrixPath(this->m_nIteration));
#endif // USE_FILE_MATRIX
    } else {
        rComm.sendData(gridMat, 0, tag);
    }
    this->log_.info("gather grid matrix: end");
}


double DfCalcGridX_Parallel::calcXCIntegForFockAndEnergy(const TlDistributeSymmetricMatrix& P_A,
                                                         DfFunctional_LDA* pFunctional,
                                                         TlDistributeSymmetricMatrix* pF_A)
{
    this->calcRho_LDA(P_A);
    double energy = this->buildVxc(pFunctional, pF_A);
    return energy;
}

double DfCalcGridX_Parallel::calcXCIntegForFockAndEnergy(const TlDistributeSymmetricMatrix& P_A,
                                                         const TlDistributeSymmetricMatrix& P_B,
                                                         DfFunctional_LDA* pFunctional,
                                                         TlDistributeSymmetricMatrix* pF_A,
                                                         TlDistributeSymmetricMatrix* pF_B)
{
    this->calcRho_LDA(P_A, P_B);
    double energy = this->buildVxc(pFunctional, pF_A, pF_B);
    return energy;
}

double DfCalcGridX_Parallel::calcXCIntegForFockAndEnergy(const TlDistributeSymmetricMatrix& P_A,
                                                         DfFunctional_GGA* pFunctional,
                                                         TlDistributeSymmetricMatrix* pF_A)
{
    this->calcRho_GGA(P_A);
    double energy = this->buildVxc(pFunctional, pF_A);
    return energy;
}

double DfCalcGridX_Parallel::calcXCIntegForFockAndEnergy(const TlDistributeSymmetricMatrix& P_A,
                                                         const TlDistributeSymmetricMatrix& P_B,
                                                         DfFunctional_GGA* pFunctional,
                                                         TlDistributeSymmetricMatrix* pF_A,
                                                         TlDistributeSymmetricMatrix* pF_B)
{
    this->calcRho_GGA(P_A, P_B);
    double energy = this->buildVxc(pFunctional, pF_A, pF_B);
    return energy;
}

TlMatrix DfCalcGridX_Parallel::getGlobalGridMatrix(const int iteration)
{
    TlCommunicate& rComm = TlCommunicate::getInstance();
    TlMatrix gridMat;
    TlMatrix rhoMat;
    index_type numOfGrids = 0;
    index_type numOfCols = 0;
    if (rComm.isMaster() == true) {
        gridMat = DfObject::getGridMatrix<TlMatrix>(iteration);
        numOfGrids = gridMat.getNumOfRows();
        numOfCols = gridMat.getNumOfCols();

        TlMatrix crdMat = gridMat.getBlockMatrix(0, 0,
                                                 numOfGrids,
                                                 GM_ATOM_INDEX +1);
        rhoMat = gridMat.getBlockMatrix(0, GM_ATOM_INDEX +1,
                                        numOfGrids,
                                        numOfCols - (GM_ATOM_INDEX +1));
        gridMat = crdMat;
    }
    rComm.broadcast(numOfGrids);
    rComm.broadcast(numOfCols);
    rComm.broadcast(gridMat);
    gridMat.resize(numOfGrids, numOfCols);

    if (this->m_bIsUpdateXC == true) {
        if (rComm.isMaster() == true) {
            gridMat.setBlockMatrix(0, GM_ATOM_INDEX +1,
                                   rhoMat);
        }
    }
    
    return gridMat;
}

void DfCalcGridX_Parallel::allReduceGridMatrix(const TlMatrix& gridMat)
{
    TlCommunicate& rComm = TlCommunicate::getInstance();
    //const index_type numOfAllGrids = gridMat.getNumOfRows();
    //const index_type numOfCols = gridMat.getNumOfCols();
    
    TlMatrix rhoMat = gridMat.getBlockMatrix(0, GM_LDA_RHO_ALPHA,
                                             gridMat.getNumOfRows(),
                                             gridMat.getNumOfCols() - GM_LDA_RHO_ALPHA);
    rComm.allReduce_SUM(rhoMat);

    if (rComm.isMaster() == true) {
        TlFileMatrix globalGridMat(DfObject::getGridMatrixPath(this->m_nIteration),
                                   gridMat.getNumOfRows(),
                                   gridMat.getNumOfCols());
        assert(globalGridMat.getNumOfRows() == gridMat.getNumOfRows());
        assert(globalGridMat.getNumOfCols() == gridMat.getNumOfCols());
        
        const TlMatrix crdMat = gridMat.getBlockMatrix(0, 0,
                                                       gridMat.getNumOfRows(),
                                                       GM_LDA_RHO_ALPHA);
        globalGridMat.setBlockMatrix(0, 0,
                                     crdMat);
        globalGridMat.setBlockMatrix(0, GM_LDA_RHO_ALPHA,
                                     rhoMat);
    }
}

void DfCalcGridX_Parallel::calcRho_LDA(const TlDistributeSymmetricMatrix& P_A)
{
    TlMatrix gridMat = this->getGlobalGridMatrix(this->m_nIteration -1);
    this->calcRho_LDA(TlDistributeMatrix(P_A), &gridMat);
    this->allReduceGridMatrix(gridMat);
}

void DfCalcGridX_Parallel::calcRho_LDA(const TlDistributeSymmetricMatrix& P_A,
                                       const TlDistributeSymmetricMatrix& P_B)
{
    TlMatrix gridMat = this->getGlobalGridMatrix(this->m_nIteration -1);
    this->calcRho_LDA(TlDistributeMatrix(P_A),
                      TlDistributeMatrix(P_B), &gridMat);
    this->allReduceGridMatrix(gridMat);
}

void DfCalcGridX_Parallel::calcRho_GGA(const TlDistributeSymmetricMatrix& P_A)
{
    TlMatrix gridMat = this->getGlobalGridMatrix(this->m_nIteration -1);
    this->calcRho_GGA(TlDistributeMatrix(P_A), &gridMat);
    this->allReduceGridMatrix(gridMat);
}

void DfCalcGridX_Parallel::calcRho_GGA(const TlDistributeSymmetricMatrix& P_A,
                                       const TlDistributeSymmetricMatrix& P_B)
{
    TlMatrix gridMat = this->getGlobalGridMatrix(this->m_nIteration -1);
    this->calcRho_GGA(TlDistributeMatrix(P_A),
                      TlDistributeMatrix(P_B), &gridMat);
    this->allReduceGridMatrix(gridMat);
}

void DfCalcGridX_Parallel::calcRho_LDA(const TlDistributeMatrix& P_A,
                                       TlMatrix* pGridMat)
{
    const TlMatrix localP_A = P_A.getLocalMatrix();
    const std::vector<index_type> rowIndexes = P_A.getRowIndexTable();
    const std::vector<index_type> colIndexes = P_A.getColIndexTable();

    const index_type numOfGrids = pGridMat->getNumOfRows();
#pragma omp parallel for schedule(runtime)
    for (index_type grid = 0; grid < numOfGrids; ++grid) {
        const TlPosition gridPosition(pGridMat->get(grid, GM_X),
                                      pGridMat->get(grid, GM_Y),
                                      pGridMat->get(grid, GM_Z));

        // calc phi table
        std::vector<WFGrid> rowPhi;
        std::vector<WFGrid> colPhi;
        this->getPhiTable(gridPosition, rowIndexes,
                          &rowPhi);
        this->getPhiTable(gridPosition, colIndexes,
                          &colPhi);

        // get rho at grid point
        double rhoA = 0.0;
        this->getRhoAtGridPoint(localP_A,
                                rowPhi, colPhi,
                                &rhoA);

        pGridMat->add(grid, GM_LDA_RHO_ALPHA, rhoA);
    }
}

void DfCalcGridX_Parallel::calcRho_LDA(const TlDistributeMatrix& P_A,
                                       const TlDistributeMatrix& P_B,
                                       TlMatrix* pGridMat)
{
    const TlMatrix localP_A = P_A.getLocalMatrix();
    const TlMatrix localP_B = P_B.getLocalMatrix();
    const std::vector<index_type> rowIndexes = P_A.getRowIndexTable();
    const std::vector<index_type> colIndexes = P_A.getColIndexTable();
    // rowIndexes, colIndexes はP_Bと共通

    const index_type numOfGrids = pGridMat->getNumOfRows();
#pragma omp parallel for schedule(runtime)
    for (index_type grid = 0; grid < numOfGrids; ++grid) {
        const TlPosition gridPosition(pGridMat->get(grid, GM_X),
                                      pGridMat->get(grid, GM_Y),
                                      pGridMat->get(grid, GM_Z));

        // calc phi table
        std::vector<WFGrid> rowPhi;
        std::vector<WFGrid> colPhi;
        this->getPhiTable(gridPosition, rowIndexes,
                          &rowPhi);
        this->getPhiTable(gridPosition, colIndexes,
                          &colPhi);

        // get rho at grid point
        double rhoA = 0.0;
        double rhoB = 0.0;
        this->getRhoAtGridPoint(localP_A,
                                rowPhi, colPhi,
                                &rhoA);
        this->getRhoAtGridPoint(localP_B,
                                rowPhi, colPhi,
                                &rhoB);

        pGridMat->add(grid, GM_LDA_RHO_ALPHA, rhoA);
        pGridMat->add(grid, GM_LDA_RHO_BETA, rhoB);
    }
}

void DfCalcGridX_Parallel::calcRho_GGA(const TlDistributeMatrix& P_A,
                                       TlMatrix* pGridMat)
{
    const TlMatrix localP_A = P_A.getLocalMatrix();
    const std::vector<index_type> rowIndexes = P_A.getRowIndexTable();
    const std::vector<index_type> colIndexes = P_A.getColIndexTable();

    const index_type numOfGrids = pGridMat->getNumOfRows();
#pragma omp parallel for schedule(runtime)
    for (index_type grid = 0; grid < numOfGrids; ++grid) {
        const TlPosition gridPosition(pGridMat->get(grid, GM_X),
                                      pGridMat->get(grid, GM_Y),
                                      pGridMat->get(grid, GM_Z));

        // calc phi table
        std::vector<WFGrid> rowPhi;
        std::vector<WFGrid> colPhi;
        std::vector<WFGrid> colGradPhiX;
        std::vector<WFGrid> colGradPhiY;
        std::vector<WFGrid> colGradPhiZ;
        this->getPhiTable(gridPosition, rowIndexes,
                          &rowPhi);
        this->getPhiTable(gridPosition, colIndexes,
                          &colPhi, &colGradPhiX, &colGradPhiY, &colGradPhiZ);

        // get rho at grid point
        double rhoA = 0.0;
        double gradRhoXA = 0.0;
        double gradRhoYA = 0.0;
        double gradRhoZA = 0.0;
        this->getRhoAtGridPoint(localP_A,
                                rowPhi,
                                colPhi, colGradPhiX, colGradPhiY, colGradPhiZ,
                                &rhoA, &gradRhoXA, &gradRhoYA, &gradRhoZA);

        pGridMat->add(grid, GM_GGA_RHO_ALPHA, rhoA);
        pGridMat->add(grid, GM_GGA_GRAD_RHO_X_ALPHA, gradRhoXA);
        pGridMat->add(grid, GM_GGA_GRAD_RHO_Y_ALPHA, gradRhoYA);
        pGridMat->add(grid, GM_GGA_GRAD_RHO_Z_ALPHA, gradRhoZA);
    }
}

void DfCalcGridX_Parallel::calcRho_GGA(const TlDistributeMatrix& P_A,
                                       const TlDistributeMatrix& P_B,
                                       TlMatrix* pGridMat)
{
    const TlMatrix localP_A = P_A.getLocalMatrix();
    const TlMatrix localP_B = P_B.getLocalMatrix();
    const std::vector<index_type> rowIndexes = P_A.getRowIndexTable();
    const std::vector<index_type> colIndexes = P_A.getColIndexTable();
    // rowIndexes, colIndexes はP_Bと共通

    const index_type numOfGrids = pGridMat->getNumOfRows();
#pragma omp parallel for schedule(runtime)
    for (index_type grid = 0; grid < numOfGrids; ++grid) {
        const TlPosition gridPosition(pGridMat->get(grid, GM_X),
                                      pGridMat->get(grid, GM_Y),
                                      pGridMat->get(grid, GM_Z));

        // calc phi table
        std::vector<WFGrid> rowPhi;
        std::vector<WFGrid> colPhi;
        std::vector<WFGrid> colGradPhiX;
        std::vector<WFGrid> colGradPhiY;
        std::vector<WFGrid> colGradPhiZ;
        this->getPhiTable(gridPosition, rowIndexes,
                          &rowPhi);
        this->getPhiTable(gridPosition, colIndexes,
                          &colPhi, &colGradPhiX, &colGradPhiY, &colGradPhiZ);

        // get rho at grid point
        double rhoA = 0.0;
        double gradRhoXA = 0.0;
        double gradRhoYA = 0.0;
        double gradRhoZA = 0.0;
        double rhoB = 0.0;
        double gradRhoXB = 0.0;
        double gradRhoYB = 0.0;
        double gradRhoZB = 0.0;
        this->getRhoAtGridPoint(localP_A,
                                rowPhi,
                                colPhi, colGradPhiX, colGradPhiY, colGradPhiZ,
                                &rhoA, &gradRhoXA, &gradRhoYA, &gradRhoZA);
        this->getRhoAtGridPoint(localP_B,
                                rowPhi,
                                colPhi, colGradPhiY, colGradPhiY, colGradPhiZ,
                                &rhoB, &gradRhoXB, &gradRhoYB, &gradRhoZB);

        pGridMat->add(grid, GM_GGA_RHO_ALPHA, rhoA);
        pGridMat->add(grid, GM_GGA_GRAD_RHO_X_ALPHA, gradRhoXA);
        pGridMat->add(grid, GM_GGA_GRAD_RHO_Y_ALPHA, gradRhoYA);
        pGridMat->add(grid, GM_GGA_GRAD_RHO_Z_ALPHA, gradRhoZA);
        pGridMat->add(grid, GM_GGA_RHO_BETA, rhoB);
        pGridMat->add(grid, GM_GGA_GRAD_RHO_X_BETA, gradRhoXB);
        pGridMat->add(grid, GM_GGA_GRAD_RHO_Y_BETA, gradRhoYB);
        pGridMat->add(grid, GM_GGA_GRAD_RHO_Z_BETA, gradRhoZB);
    }
}

void DfCalcGridX_Parallel::getPhiTable(const TlPosition& gridPosition,
                                       const std::vector<index_type>& AO_List,
                                       std::vector<WFGrid>* pPhis)
{
    assert(pPhis != NULL);
    double densityCutOffValue = 1.0E-16;
    
    const std::size_t AO_ListSize = AO_List.size();
    pPhis->clear();
    pPhis->reserve(AO_ListSize);

    // orbital loop
    for (std::size_t AO_index = 0; AO_index < AO_ListSize; ++AO_index) {
        const index_type AO = AO_List[AO_index];

        double phi = 0.0;

        const TlPosition pos = gridPosition - this->m_tlOrbInfo.getPosition(AO);
        const double distance2 = pos.squareDistanceFrom();
        const int basisType = this->m_tlOrbInfo.getBasisType(AO);
        const double prefactor = DfCalcGridX::getPrefactor(basisType, pos);
        const int contraction = this->m_tlOrbInfo.getCgtoContraction(AO);

        for (int PGTO = 0; PGTO < contraction; ++PGTO) {
            // double prefactorX = 0.0;
            // double prefactorY = 0.0;
            // double prefactorZ = 0.0;
            const double alpha = this->m_tlOrbInfo.getExponent(AO, PGTO);
            const double shoulder = alpha * distance2;

            if (shoulder <= TOOBIG) {
                const double g = prefactor
                    * this->m_tlOrbInfo.getCoefficient(AO, PGTO)
                    * std::exp(-shoulder);
                phi += g;
            }
        }
        
        if (std::fabs(phi) > densityCutOffValue) {
            //const WFGrid wfGrid(AO, phi);
            const WFGrid wfGrid(AO_index, phi);
            pPhis->push_back(wfGrid);
        }
    }

    std::sort(pPhis->begin(), pPhis->end(), WFGrid_sort_functional());
}

void DfCalcGridX_Parallel::getPhiTable(const TlPosition& gridPosition,
                                       const std::vector<index_type>& AO_List,
                                       std::vector<WFGrid>* pPhis,
                                       std::vector<WFGrid>* pGradPhiXs,
                                       std::vector<WFGrid>* pGradPhiYs,
                                       std::vector<WFGrid>* pGradPhiZs)
{
    assert(pPhis != NULL);
    assert(pGradPhiXs != NULL);
    assert(pGradPhiYs != NULL);
    assert(pGradPhiZs != NULL);
    const double densityCutOffValue = 1.0E-16;
    
    const std::size_t AO_ListSize = AO_List.size();
    pPhis->clear();
    pPhis->reserve(AO_ListSize);
    pGradPhiXs->clear();
    pGradPhiXs->reserve(AO_ListSize);
    pGradPhiYs->clear();
    pGradPhiYs->reserve(AO_ListSize);
    pGradPhiZs->clear();
    pGradPhiZs->reserve(AO_ListSize);

    // orbital loop
    for (std::size_t AO_index = 0; AO_index < AO_ListSize; ++AO_index) {
        const index_type AO = AO_List[AO_index];

        double phi = 0.0;
        double gradPhiX = 0.0;
        double gradPhiY = 0.0;
        double gradPhiZ = 0.0;

        const TlPosition pos = gridPosition - this->m_tlOrbInfo.getPosition(AO);
        const double distance2 = pos.squareDistanceFrom();
        const int basisType = this->m_tlOrbInfo.getBasisType(AO);
        const double prefactor = DfCalcGridX::getPrefactor(basisType, pos);
        const int contraction = this->m_tlOrbInfo.getCgtoContraction(AO);

        for (int PGTO = 0; PGTO < contraction; ++PGTO) {
            double prefactorX = 0.0;
            double prefactorY = 0.0;
            double prefactorZ = 0.0;
            const double alpha = this->m_tlOrbInfo.getExponent(AO, PGTO);
            const double shoulder = alpha * distance2;

            if (shoulder <= TOOBIG) {
                const double g = this->m_tlOrbInfo.getCoefficient(AO, PGTO) * std::exp(-shoulder);
                phi += prefactor * g;
                DfCalcGridX::getPrefactorForDerivative(basisType, alpha, pos,
                                                       &prefactorX, &prefactorY, &prefactorZ);
                gradPhiX += prefactorX * g;
                gradPhiY += prefactorY * g;
                gradPhiZ += prefactorZ * g;
            }
        }
        
        if (std::fabs(phi) > densityCutOffValue) {
            //const WFGrid wfGrid(AO, phi);
            const WFGrid wfGrid(AO_index, phi);
            pPhis->push_back(wfGrid);
        }

        if (std::fabs(gradPhiX) > densityCutOffValue) {
            //const WFGrid wfGrid(AO, gradPhiX);
            const WFGrid wfGrid(AO_index, gradPhiX);
            pGradPhiXs->push_back(wfGrid);
        }
        if (std::fabs(gradPhiY) > densityCutOffValue) {
            //const WFGrid wfGrid(AO, gradPhiY);
            const WFGrid wfGrid(AO_index, gradPhiY);
            pGradPhiYs->push_back(wfGrid);
        }
        if (std::fabs(gradPhiZ) > densityCutOffValue) {
            //const WFGrid wfGrid(AO, gradPhiZ);
            const WFGrid wfGrid(AO_index, gradPhiZ);
            pGradPhiZs->push_back(wfGrid);
        }
    }

    std::sort(pPhis->begin(), pPhis->end(), WFGrid_sort_functional());
    std::sort(pGradPhiXs->begin(), pGradPhiXs->end(), WFGrid_sort_functional());
    std::sort(pGradPhiYs->begin(), pGradPhiYs->end(), WFGrid_sort_functional());
    std::sort(pGradPhiZs->begin(), pGradPhiZs->end(), WFGrid_sort_functional());
}

// 分散行列専用
void DfCalcGridX_Parallel::getRhoAtGridPoint(const TlMatrix& P,
                                             const std::vector<WFGrid>& rowPhi,
                                             const std::vector<WFGrid>& colPhi,
                                             double* pRhoA)
{
    const double densityCutOffValue = this->m_densityCutOffValueA;
    double dRho = 0.0;
    // double dGradRhoX = 0.0;
    // double dGradRhoY = 0.0;
    // double dGradRhoZ = 0.0;

    // 密度行列のループ(index = p)の最大を求める
    double value = densityCutOffValue;

    std::vector<WFGrid>::const_iterator pEnd = std::upper_bound(rowPhi.begin(),
                                                                rowPhi.end(),
                                                                WFGrid(0, value),
                                                                WFGrid_sort_functional());
    const std::size_t max_p = std::distance(rowPhi.begin(), pEnd);
    for (std::size_t p = 0; p < max_p; ++p) {
        const std::size_t nOrb_p = rowPhi[p].index;
        const double phi_p = rowPhi[p].value;
        const double cutValue = std::fabs(densityCutOffValue / phi_p);

        // 高速化
        const TlVector P_row = P.getRowVector(nOrb_p);

        std::vector<WFGrid>::const_iterator qEnd = std::upper_bound(colPhi.begin(),
                                                                    colPhi.end(),
                                                                    WFGrid(0, cutValue),
                                                                    WFGrid_sort_functional());
        const std::size_t max_q = std::distance(colPhi.begin(), qEnd);
        for (std::size_t q = 0; q < max_q; ++q) {
            const std::size_t nOrb_q = colPhi[q].index;
            const double phi_q = colPhi[q].value;
            dRho += P_row[nOrb_q] * phi_p * phi_q;
        }
    }

    *pRhoA = dRho;
}


void DfCalcGridX_Parallel::getRhoAtGridPoint(const TlMatrix& P,
                                             const std::vector<WFGrid>& rowPhi,
                                             const std::vector<WFGrid>& colPhi,
                                             const std::vector<WFGrid>& colGradPhiX,
                                             const std::vector<WFGrid>& colGradPhiY,
                                             const std::vector<WFGrid>& colGradPhiZ,
                                             double* pRhoA,
                                             double* pGradRhoAX,
                                             double* pGradRhoAY,
                                             double* pGradRhoAZ)
{
    const double densityCutOffValue = 1.0E-16;
    double dRho = 0.0;
    double dGradRhoX = 0.0;
    double dGradRhoY = 0.0;
    double dGradRhoZ = 0.0;

    // 密度行列のループ(index = p)の最大を求める
    //double value = 1.0E-16;
    // if (aGradPhiX.size() > 0) {
    //     value = std::min(value, std::fabs(aGradPhiX[0].value));
    // }
    // if (aGradPhiY.size() > 0) {
    //     value = std::min(value, std::fabs(aGradPhiY[0].value));
    // }
    // if (aGradPhiZ.size() > 0) {
    //     value = std::min(value, std::fabs(aGradPhiZ[0].value));
    // }

    // std::size_t max_q = aPhi.size();
    // std::size_t max_qx = aGradPhiX.size();
    // std::size_t max_qy = aGradPhiY.size();
    // std::size_t max_qz = aGradPhiZ.size();

    std::vector<WFGrid>::const_iterator pEnd = std::upper_bound(rowPhi.begin(),
                                                                rowPhi.end(),
                                                                WFGrid(0, densityCutOffValue),
                                                                WFGrid_sort_functional());
    const std::size_t max_p = std::distance(rowPhi.begin(), pEnd);
    for (std::size_t p = 0; p < max_p; ++p) {
        const index_type nOrb_p = rowPhi[p].index;
        const double phi_p = rowPhi[p].value;
        const double cutValue = std::fabs(densityCutOffValue / phi_p);

        // 高速化
        const TlVector P_row = P.getRowVector(nOrb_p);

        std::vector<WFGrid>::const_iterator qEnd = std::upper_bound(colPhi.begin(),
                                                                    colPhi.end(),
                                                                    WFGrid(0, cutValue),
                                                                    WFGrid_sort_functional());
        const std::size_t max_q = std::distance(colPhi.begin(), qEnd);
        for (std::size_t q = 0; q < max_q; ++q) {
            const std::size_t nOrb_q = colPhi[q].index;
            const double phi_q = colPhi[q].value;
            dRho += P_row[nOrb_q] * phi_p * phi_q;
        }
        
        std::vector<WFGrid>::const_iterator qxEnd = std::upper_bound(colGradPhiX.begin(),
                                                                     colGradPhiX.end(),
                                                                     WFGrid(0, cutValue),
                                                                     WFGrid_sort_functional());
        const std::size_t max_qx = std::distance(colGradPhiX.begin(), qxEnd);
        for (std::size_t qx = 0; qx < max_qx; ++qx) {
            const index_type nOrb_q = colGradPhiX[qx].index;
            const double phi_q = colGradPhiX[qx].value;
            dGradRhoX += P_row[nOrb_q] * phi_p * phi_q;
        }

        std::vector<WFGrid>::const_iterator qyEnd = std::upper_bound(colGradPhiY.begin(),
                                                                     colGradPhiY.end(),
                                                                     WFGrid(0, cutValue),
                                                                     WFGrid_sort_functional());
        const std::size_t max_qy = std::distance(colGradPhiY.begin(), qyEnd);
        for (std::size_t qy = 0; qy < max_qy; ++qy) {
            const std::size_t nOrb_q = colGradPhiY[qy].index;
            const double phi_q = colGradPhiY[qy].value;
            dGradRhoY += P_row[nOrb_q] * phi_p * phi_q;
        }

        std::vector<WFGrid>::const_iterator qzEnd = std::upper_bound(colGradPhiZ.begin(),
                                                                     colGradPhiZ.end(),
                                                                     WFGrid(0, cutValue),
                                                                     WFGrid_sort_functional());
        const std::size_t max_qz = std::distance(colGradPhiZ.begin(), qzEnd);
        for (std::size_t qz = 0; qz < max_qz; ++qz) {
            const std::size_t nOrb_q = colGradPhiZ[qz].index;
            const double phi_q = colGradPhiZ[qz].value;
            dGradRhoZ += P_row[nOrb_q] * phi_p * phi_q;
        }
    }

    *pRhoA = dRho;
    *pGradRhoAX = dGradRhoX;
    *pGradRhoAY = dGradRhoY;
    *pGradRhoAZ = dGradRhoZ;
}

double DfCalcGridX_Parallel::buildVxc(DfFunctional_LDA* pFunctional,
                                      TlDistributeSymmetricMatrix* pF_A)
{
    const TlMatrix gridMat = this->distributeGridMatrix(this->m_nIteration);
    TlSparseSymmetricMatrix tmpF_A(this->m_nNumOfAOs);
    double energy = DfCalcGridX::buildVxc(gridMat, pFunctional, &tmpF_A);

    TlCommunicate& rComm = TlCommunicate::getInstance();
    pF_A->mergeSparseMatrix(tmpF_A);
    rComm.allReduce_SUM(energy);
    return energy;
}

double DfCalcGridX_Parallel::buildVxc(DfFunctional_LDA* pFunctional,
                                      TlDistributeSymmetricMatrix* pF_A,
                                      TlDistributeSymmetricMatrix* pF_B)
{
    const TlMatrix gridMat = this->distributeGridMatrix(this->m_nIteration);
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
                                      TlDistributeSymmetricMatrix* pF_A)
{
    const TlMatrix gridMat = this->distributeGridMatrix(this->m_nIteration);
    TlSparseSymmetricMatrix tmpF_A(this->m_nNumOfAOs);
    double energy = DfCalcGridX::buildVxc(gridMat, pFunctional, &tmpF_A);

    TlCommunicate& rComm = TlCommunicate::getInstance();
    pF_A->mergeSparseMatrix(tmpF_A);
    rComm.allReduce_SUM(energy);
    return energy;
}

double DfCalcGridX_Parallel::buildVxc(DfFunctional_GGA* pFunctional,
                                      TlDistributeSymmetricMatrix* pF_A,
                                      TlDistributeSymmetricMatrix* pF_B)
{
    const TlMatrix gridMat = this->distributeGridMatrix(this->m_nIteration);
    TlSparseSymmetricMatrix tmpF_A(this->m_nNumOfAOs);
    TlSparseSymmetricMatrix tmpF_B(this->m_nNumOfAOs);
    double energy = DfCalcGridX::buildVxc(gridMat, pFunctional, &tmpF_A, &tmpF_B);

    TlCommunicate& rComm = TlCommunicate::getInstance();
    pF_A->mergeSparseMatrix(tmpF_A);
    pF_B->mergeSparseMatrix(tmpF_B);
    rComm.allReduce_SUM(energy);
    return energy;
}

void DfCalcGridX_Parallel::getWholeDensity(double* pRhoA, double* pRhoB) const
{
    TlCommunicate& rComm = TlCommunicate::getInstance();
    if (rComm.isMaster() == true) {
        DfCalcGridX::getWholeDensity(pRhoA, pRhoB);
    }
}

