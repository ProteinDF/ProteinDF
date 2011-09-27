#include <cmath>
#include <iostream>
#include <algorithm>
#include <limits>
#include <cassert>

#include "DfCalcGridX.h"
#include "CnError.h"

#include "Fl_Geometry.h"
#include "TlFile.h"
#include "TlUtils.h"

////////////////////////////////////////////////////////////////////////
// DfCalcGridX
//

const double DfCalcGridX::TOOBIG = 30.0;
const double DfCalcGridX::EPS = std::numeric_limits<double>::epsilon();
const double DfCalcGridX::INV_SQRT3 = 1.0 / std::sqrt(3.0);
const double DfCalcGridX::INV_SQRT12 = 1.0 / std::sqrt(12.0);

DfCalcGridX::DfCalcGridX(TlSerializeData* pPdfParam)
    : DfObject(pPdfParam), m_tlOrbInfo((*pPdfParam)["model"]["coordinates"], (*pPdfParam)["model"]["basis_set"])
{
    const TlSerializeData& pdfParam = *pPdfParam;
    
    this->inputtedDensityCutoffValue_ = 1.0E-16;
    if (!(pdfParam["model"]["xc-density-threshold"].getStr().empty())) {
        this->inputtedDensityCutoffValue_ = pdfParam["model"]["xc-density-threshold"].getDouble();
    }
    this->m_densityCutOffValueA = this->inputtedDensityCutoffValue_;
    this->m_densityCutOffValueB = this->inputtedDensityCutoffValue_;

    this->m_inputedCutoffThreshold = pdfParam["model"]["cut-value"].getDouble();

    this->physicalValues_.clear();

    this->readTable();

    // for debug
    this->isDebugOutPhiTable_ = (TlUtils::toUpper(pdfParam["model"]["debug_out_phi_table"].getStr()) == "YES") ? true : false;
}

DfCalcGridX::~DfCalcGridX()
{
}


void DfCalcGridX::defineCutOffValues(const TlSymmetricMatrix& P)
{
    const double maxValueOfP = std::max(P.getMaxAbsoluteElement(), 1.0E-16);
    if (maxValueOfP < 1.0) {
        this->m_densityCutOffValueA /= maxValueOfP;
    }
    this->logger(TlUtils::format(" density cutoff value = %e\n",
                                 this->m_densityCutOffValueA));
}


void DfCalcGridX::defineCutOffValues(const TlSymmetricMatrix& PA,
                                     const TlSymmetricMatrix& PB)
{
    const double maxValueOfPA = std::max(PA.getMaxAbsoluteElement(), 1.0E-16);
    const double maxValueOfPB = std::max(PB.getMaxAbsoluteElement(), 1.0E-16);
    if (maxValueOfPA < 1.0) {
        this->m_densityCutOffValueA /= maxValueOfPA;
    }
    if (maxValueOfPB < 1.0) {
        this->m_densityCutOffValueB /= maxValueOfPB;
    }
    this->logger(TlUtils::format(" density cutoff value(alpha) = %e\n",
                                 this->m_densityCutOffValueA));
    this->logger(TlUtils::format(" density cutoff value(beta ) = %e\n",
                                 this->m_densityCutOffValueB));
}


void DfCalcGridX::backupGridData()
{
    const std::string prevGridDataFilePath = TlUtils::format("%s.itr%d", this->getGridDataFilePath().c_str(), this->m_nIteration -1);
    TlFile::copy(this->getGridDataFilePath(), prevGridDataFilePath);
}


// Fockの交換相関項と、エネルギーを同時に求める
// for LDA and RKS
double DfCalcGridX::calcXCIntegForFockAndEnergy(const TlSymmetricMatrix& P,
                                                DfFunctional_LDA* pFunctional,
                                                TlSymmetricMatrix* pF)
{
    if (this->enableExperimentalCode_ != true) {
        return this->calcXCIntegForFockAndEnergy1(P, pFunctional, pF);
    } else {
        this->logger(" !!! experimental code !!!\n");
        return this->calcXCIntegForFockAndEnergy2(P, pFunctional, pF);
    }
}


double DfCalcGridX::calcXCIntegForFockAndEnergy1(const TlSymmetricMatrix& P,
                                                 DfFunctional_LDA* pFunctional,
                                                 TlSymmetricMatrix* pF)
{
    assert(pFunctional != NULL);
    assert(pF != NULL);

    this->backupGridData();

    // setup
    const int nNumOfAtoms = this->numOfRealAtoms_;
    this->physicalValues_.clear();
    this->defineCutOffValues(P);

    // calc
    const double dEnergy = this->calcXCIntegForFockAndEnergy(0, nNumOfAtoms,
                                                             P, pFunctional, pF);

    if (this->m_bIsUpdateXC == true) {
        this->flushGridData();
    }

    return dEnergy;
}


// Fockの交換相関項と、エネルギーを同時に求める
// for LDA and UKS
double DfCalcGridX::calcXCIntegForFockAndEnergy(const TlSymmetricMatrix& PA,
                                                const TlSymmetricMatrix& PB,
                                                DfFunctional_LDA* pFunctional,
                                                TlSymmetricMatrix* pFA,
                                                TlSymmetricMatrix* pFB)
{
    assert(pFunctional != NULL);
    assert(pFA != NULL);
    assert(pFB != NULL);

    this->backupGridData();

    // setup
    const int nNumOfAtoms = this->numOfRealAtoms_;
    this->physicalValues_.clear();
    this->defineCutOffValues(PA, PB);

    // calc
    const double dEnergy = this->calcXCIntegForFockAndEnergy(0, nNumOfAtoms,
                                                             PA, PB, pFunctional,
                                                             pFA, pFB);
    if (this->m_bIsUpdateXC == true) {
        this->flushGridData();
    }

    return dEnergy;
}


// Fockの交換相関項と、エネルギーを同時に求める
// for LDA and RKS (共通部分)
double DfCalcGridX::calcXCIntegForFockAndEnergy(const int nStartAtom, const int nEndAtom,
                                                const TlSymmetricMatrix& P,
                                                DfFunctional_LDA* pFunctional,
                                                TlSymmetricMatrix* pF)
{
    assert(pFunctional != NULL);
    assert(pF != NULL);

    double dEnergy = 0.0;
    const double densityCutOffValue = this->m_densityCutOffValueA;
    const GridDataManager gdm(this->getGridDataFilePath());

    for (int nAtom = nStartAtom; nAtom < nEndAtom; ++nAtom) {
        std::vector<GridDataManager::GridInfo> grids;
        {
            const std::vector<double> coordX = gdm.getData(nAtom, GridDataManager::COORD_X);
            const std::vector<double> coordY = gdm.getData(nAtom, GridDataManager::COORD_Y);
            const std::vector<double> coordZ = gdm.getData(nAtom, GridDataManager::COORD_Z);
            const std::vector<double> weight = gdm.getData(nAtom, GridDataManager::GRID_WEIGHT);
            grids = GridDataManager::composeGridInfo(coordX, coordY, coordZ, weight);
        }
        const int nNumOfGrids = grids.size();

        // calc rho
        std::vector<double> rhoA(nNumOfGrids, 0.0);
        if ((this->m_bIsUpdateXC == true) && (this->m_nIteration >= 2)) {
            rhoA = gdm.getData(nAtom, GridDataManager::DENSITY);
        }

#pragma omp parallel for schedule(runtime)
        for (int nGrid = 0; nGrid < nNumOfGrids; ++nGrid) {
            const TlPosition gridPosition = grids[nGrid].position;
            const double dWeight = grids[nGrid].weight;

            // calc phi table
            std::vector<WFGrid> aPhi;
            this->getPhiTable(gridPosition, aPhi);
            std::sort(aPhi.begin(), aPhi.end(), WFGrid_sort_functional());
            // debug
            // for (std::vector<WFGrid>::iterator it = aPhi.begin(); it != aPhi.end(); ++it) {
            //     std::cerr << it->value << std::endl;
            // }
            
            // get rho at grid point
            double dRhoA = 0.0;
            this->getRhoAtGridPoint(P, aPhi, &dRhoA);
            // debug
            // std::cerr << dRhoA << std::endl;

            dRhoA *= 0.5;
            dRhoA += rhoA[nGrid];

            // build fock matrix
            if (dRhoA > densityCutOffValue) {
#pragma omp critical (DfCalcGridX_calcXCIntegForFockAndEnergy_LDA_R)
                {
                    this->buildFock(dRhoA, aPhi, pFunctional, dWeight, pF); // RKS code
                    dEnergy += dWeight * pFunctional->getFunctional(dRhoA); // RKS code
                }
            }

            // save rho data
            rhoA[nGrid] = dRhoA;
        }

        if (this->m_bIsUpdateXC == true) {
            this->physicalValues_[GridDataManager::DENSITY][nAtom] = rhoA;
        }
    }

    return dEnergy;
}


// Fockの交換相関項と、エネルギーを同時に求める
// for LDA and UKS (共通部分)
double DfCalcGridX::calcXCIntegForFockAndEnergy(const int nStartAtom, const int nEndAtom,
                                                const TlSymmetricMatrix& PA,
                                                const TlSymmetricMatrix& PB,
                                                DfFunctional_LDA* pFunctional,
                                                TlSymmetricMatrix* pFA,
                                                TlSymmetricMatrix* pFB)
{
    assert(pFunctional != NULL);
    assert(pFA != NULL);
    assert(pFB != NULL);

    double dEnergy = 0.0;
    const double densityCutOffValueA = this->m_densityCutOffValueA;
    const double densityCutOffValueB = this->m_densityCutOffValueB;
    const GridDataManager gdm(this->getGridDataFilePath());
    std::map<int, std::vector<double> > storedRhoA;
    std::map<int, std::vector<double> > storedRhoB;

    for (int nAtom = nStartAtom; nAtom < nEndAtom; ++nAtom) {
        // read data of grid and weight
        std::vector<GridDataManager::GridInfo> grids;
        {
            const std::vector<double> coordX = gdm.getData(nAtom, GridDataManager::COORD_X);
            const std::vector<double> coordY = gdm.getData(nAtom, GridDataManager::COORD_Y);
            const std::vector<double> coordZ = gdm.getData(nAtom, GridDataManager::COORD_Z);
            const std::vector<double> weight = gdm.getData(nAtom, GridDataManager::GRID_WEIGHT);
            grids = GridDataManager::composeGridInfo(coordX, coordY, coordZ, weight);
        }
        const int nNumOfGrids = grids.size();
        
        // calc rho
        std::vector<double> rhoA(nNumOfGrids, 0.0);
        std::vector<double> rhoB(nNumOfGrids, 0.0);
        if ((this->m_bIsUpdateXC == true) && (this->m_nIteration >= 2)) {
            rhoA = gdm.getData(nAtom, GridDataManager::DENSITY_ALPHA);
            rhoB = gdm.getData(nAtom, GridDataManager::DENSITY_BETA);
        }

#pragma omp parallel for schedule(runtime)
        for (int nGrid = 0; nGrid < nNumOfGrids; ++nGrid) {
            const TlPosition gridPosition = grids[nGrid].position;
            const double dWeight = grids[nGrid].weight;

            // calc phi table
            std::vector<WFGrid> aPhi;
            this->getPhiTable(gridPosition, aPhi);
            std::sort(aPhi.begin(), aPhi.end(), WFGrid_sort_functional());

            // get rho at grid point
            double dRhoA, dRhoB;
            this->getRhoAtGridPoint(PA, aPhi, &dRhoA);
            this->getRhoAtGridPoint(PB, aPhi, &dRhoB);
            dRhoA += rhoA[nGrid];
            dRhoB += rhoB[nGrid];

            // build fock matrix
            if ((dRhoA > densityCutOffValueA) || (dRhoB > densityCutOffValueB)) {
                const double rA = std::max(dRhoA, 0.0);
                const double rB = std::max(dRhoB, 0.0);

#pragma omp critical (DfCalcGridX_calcXCIntegForFockAndEnergy_LDA_U)
                {
                    this->buildFock(rA, rB, aPhi, pFunctional, dWeight, pFA, pFB);
                    dEnergy += dWeight * pFunctional->getFunctional(rA, rB);
                }
            }

            // save rho data
            rhoA[nGrid] = dRhoA;
            rhoB[nGrid] = dRhoB;
        }

        if (this->m_bIsUpdateXC == true) {
            this->physicalValues_[GridDataManager::DENSITY_ALPHA][nAtom] = rhoA;
            this->physicalValues_[GridDataManager::DENSITY_BETA][nAtom] = rhoB;
        }
    }

    return dEnergy;
}


void DfCalcGridX::flushGridData()
{
    typedef std::map<GridDataManager::ChunkType, std::map<int, std::vector<double> > > physicalValuesType;
    GridDataManager gdm(this->getGridDataFilePath());

    physicalValuesType::const_iterator pEnd = this->physicalValues_.end();
    for (physicalValuesType::const_iterator p = this->physicalValues_.begin(); p != pEnd; ++p) {
        const GridDataManager::ChunkType chunkType = p->first;
        std::map<int, std::vector<double> >::const_iterator qEnd = p->second.end();
        for (std::map<int, std::vector<double> >::const_iterator q = p->second.begin(); q != qEnd; ++q) {
            gdm.setData(q->first, chunkType, q->second);
        }
    }

    this->physicalValues_.clear();
}


// GGA, RKS 版
double DfCalcGridX::calcXCIntegForFockAndEnergy(const TlSymmetricMatrix& P,
                                                DfFunctional_GGA* pFunctional,
                                                TlSymmetricMatrix* pF)
{
    if (this->enableExperimentalCode_ != true) {
        return this->calcXCIntegForFockAndEnergy1(P, pFunctional, pF);
    } else {
        this->logger(" !!! experimental code !!!\n");
        return this->calcXCIntegForFockAndEnergy2(P, pFunctional, pF);
    }
}


double DfCalcGridX::calcXCIntegForFockAndEnergy1(const TlSymmetricMatrix& P,
                                                 DfFunctional_GGA* pFunctional,
                                                 TlSymmetricMatrix* pF)
{
    assert(pFunctional != NULL);
    assert(pF != NULL);

    this->backupGridData();

    // setup
    const int nNumOfAtoms = this->numOfRealAtoms_; //this->m_nNumOfAtoms - this->m_nNumOfDummyAtoms;
    this->physicalValues_.clear();
    this->defineCutOffValues(P);

    // calc
    const double dEnergy = this->calcXCIntegForFockAndEnergy(0, nNumOfAtoms,
                                                             P, pFunctional,
                                                             pF);
    if (this->m_bIsUpdateXC == true) {
        this->flushGridData();
    }

    return dEnergy;
}


// GGA, UKS
double DfCalcGridX::calcXCIntegForFockAndEnergy(const TlSymmetricMatrix& PA,
                                                const TlSymmetricMatrix& PB,
                                                DfFunctional_GGA* pFunctional,
                                                TlSymmetricMatrix* pFA,
                                                TlSymmetricMatrix* pFB)
{
    assert(pFunctional != NULL);
    assert(pFA != NULL);
    assert(pFB != NULL);

    // backup grid data
    this->backupGridData();

    // setup
    const int nNumOfAtoms = this->numOfRealAtoms_; //this->m_nNumOfAtoms - this->m_nNumOfDummyAtoms;
    this->physicalValues_.clear();
    this->defineCutOffValues(PA, PB);

    // calc
    const double dEnergy = this->calcXCIntegForFockAndEnergy(0, nNumOfAtoms,
                                                             PA, PB, pFunctional,
                                                             pFA, pFB);
    if (this->m_bIsUpdateXC == true) {
        this->flushGridData();
    }

    return dEnergy;
}


// for GGA and RKS
double DfCalcGridX::calcXCIntegForFockAndEnergy(const int nStartAtom, const int nEndAtom,
                                                const TlSymmetricMatrix& P,
                                                DfFunctional_GGA* pFunctional,
                                                TlSymmetricMatrix* pF)
{
    assert(pFunctional != NULL);
    assert(pF != NULL);

    double dEnergy = 0.0;
    const double densityCutOffValue = this->m_densityCutOffValueA;
    const GridDataManager gdm(this->getGridDataFilePath());

    for (int nAtom = nStartAtom; nAtom < nEndAtom; ++nAtom) {
        std::vector<GridDataManager::GridInfo> grids;
        {
            const std::vector<double> coordX = gdm.getData(nAtom, GridDataManager::COORD_X);
            const std::vector<double> coordY = gdm.getData(nAtom, GridDataManager::COORD_Y);
            const std::vector<double> coordZ = gdm.getData(nAtom, GridDataManager::COORD_Z);
            const std::vector<double> weight = gdm.getData(nAtom, GridDataManager::GRID_WEIGHT);
            grids = GridDataManager::composeGridInfo(coordX, coordY, coordZ, weight);
        }
        const int nNumOfGrids = grids.size();

        // calc rho
        std::vector<double> rhoA(nNumOfGrids, 0.0);
        std::vector<double> gradRhoAX(nNumOfGrids, 0.0);
        std::vector<double> gradRhoAY(nNumOfGrids, 0.0);
        std::vector<double> gradRhoAZ(nNumOfGrids, 0.0);
        if ((this->m_bIsUpdateXC == true) && (this->m_nIteration >= 2)) {
            rhoA = gdm.getData(nAtom, GridDataManager::DENSITY);
            gradRhoAX = gdm.getData(nAtom, GridDataManager::GRAD_DENSITY_X);
            gradRhoAY = gdm.getData(nAtom, GridDataManager::GRAD_DENSITY_Y);
            gradRhoAZ = gdm.getData(nAtom, GridDataManager::GRAD_DENSITY_Z);
        }

#pragma omp parallel for schedule(runtime)
        for (int nGrid = 0; nGrid < nNumOfGrids; ++nGrid) {
            const TlPosition gridPosition = grids[nGrid].position;
            const double dWeight = grids[nGrid].weight;

            // calc phi table
            std::vector<WFGrid> aPhi;
            std::vector<WFGrid> aGradPhiX;
            std::vector<WFGrid> aGradPhiY;
            std::vector<WFGrid> aGradPhiZ;
            this->getPhiTable(gridPosition, aPhi, aGradPhiX, aGradPhiY, aGradPhiZ);
            std::sort(aPhi.begin(), aPhi.end(), WFGrid_sort_functional());
            std::sort(aGradPhiX.begin(), aGradPhiX.end(), WFGrid_sort_functional());
            std::sort(aGradPhiY.begin(), aGradPhiY.end(), WFGrid_sort_functional());
            std::sort(aGradPhiZ.begin(), aGradPhiZ.end(), WFGrid_sort_functional());
            // // debug
            // for (std::vector<WFGrid>::iterator it = aPhi.begin(); it != aPhi.end(); ++it) {
            //     std::cerr << it->value << std::endl;
            // }

            // get rho at grid point
            double dRhoA = 0.0;
            double dGradRhoAX = 0.0;
            double dGradRhoAY = 0.0;
            double dGradRhoAZ = 0.0;
            this->getRhoAtGridPoint(P, aPhi, aGradPhiX, aGradPhiY, aGradPhiZ,
                                    &dRhoA, &dGradRhoAX, &dGradRhoAY, &dGradRhoAZ);
            // debug
            // std::cerr << dRhoA << std::endl;
            
            dRhoA *= 0.5;
            dGradRhoAX *= 0.5;
            dGradRhoAY *= 0.5;
            dGradRhoAZ *= 0.5;
            dRhoA += rhoA[nGrid];
            dGradRhoAX += gradRhoAX[nGrid];
            dGradRhoAY += gradRhoAY[nGrid];
            dGradRhoAZ += gradRhoAZ[nGrid];

            // calc
            if (dRhoA > densityCutOffValue) {
                const double gammaAA =  dGradRhoAX*dGradRhoAX + dGradRhoAY*dGradRhoAY + dGradRhoAZ*dGradRhoAZ;

#pragma omp critical (DfCalcGridX_calcXCIntegForFockAndEnergy_GGA_R)
                {
                    this->buildFock(dRhoA, dGradRhoAX, dGradRhoAY, dGradRhoAZ,
                                    aPhi, aGradPhiX, aGradPhiY, aGradPhiZ,
                                    pFunctional, dWeight, pF); // RKS code
                    dEnergy += dWeight * pFunctional->getFunctional(dRhoA, gammaAA); // RKS code
                }
            }

            rhoA[nGrid] = dRhoA;
            gradRhoAX[nGrid] = dGradRhoAX;
            gradRhoAY[nGrid] = dGradRhoAY;
            gradRhoAZ[nGrid] = dGradRhoAZ;
        }

        if (this->m_bIsUpdateXC == true) {
            this->physicalValues_[GridDataManager::DENSITY][nAtom] = rhoA;
            this->physicalValues_[GridDataManager::GRAD_DENSITY_X][nAtom] = gradRhoAX;
            this->physicalValues_[GridDataManager::GRAD_DENSITY_Y][nAtom] = gradRhoAY;
            this->physicalValues_[GridDataManager::GRAD_DENSITY_Z][nAtom] = gradRhoAZ;
        }
    }

    return dEnergy;
}


// for GGA and UKS
double DfCalcGridX::calcXCIntegForFockAndEnergy(const int nStartAtom, const int nEndAtom,
                                                const TlSymmetricMatrix& PA,
                                                const TlSymmetricMatrix& PB,
                                                DfFunctional_GGA* pFunctional,
                                                TlSymmetricMatrix* pFA,
                                                TlSymmetricMatrix* pFB)
{
    assert(pFunctional != NULL);
    assert(pFA != NULL);
    assert(pFB != NULL);

    double dEnergy = 0.0;
    const double densityCutOffValueA = this->m_densityCutOffValueA;
    const double densityCutOffValueB = this->m_densityCutOffValueB;
    const GridDataManager gdm(this->getGridDataFilePath());
    std::map<int, std::vector<double> > storedRhoA;
    std::map<int, std::vector<double> > storedGradRhoAX;
    std::map<int, std::vector<double> > storedGradRhoAY;
    std::map<int, std::vector<double> > storedGradRhoAZ;
    std::map<int, std::vector<double> > storedRhoB;
    std::map<int, std::vector<double> > storedGradRhoBX;
    std::map<int, std::vector<double> > storedGradRhoBY;
    std::map<int, std::vector<double> > storedGradRhoBZ;

    for (int nAtom = nStartAtom; nAtom < nEndAtom; ++nAtom) {
        std::vector<GridDataManager::GridInfo> grids;
        {
            const std::vector<double> coordX = gdm.getData(nAtom, GridDataManager::COORD_X);
            const std::vector<double> coordY = gdm.getData(nAtom, GridDataManager::COORD_Y);
            const std::vector<double> coordZ = gdm.getData(nAtom, GridDataManager::COORD_Z);
            const std::vector<double> weight = gdm.getData(nAtom, GridDataManager::GRID_WEIGHT);
            grids = GridDataManager::composeGridInfo(coordX, coordY, coordZ, weight);
        }
        const int nNumOfGrids = grids.size();

        // calc rho
        std::vector<double> rhoA(nNumOfGrids, 0.0);
        std::vector<double> rhoB(nNumOfGrids, 0.0);
        std::vector<double> gradRhoAX(nNumOfGrids, 0.0);
        std::vector<double> gradRhoAY(nNumOfGrids, 0.0);
        std::vector<double> gradRhoAZ(nNumOfGrids, 0.0);
        std::vector<double> gradRhoBX(nNumOfGrids, 0.0);
        std::vector<double> gradRhoBY(nNumOfGrids, 0.0);
        std::vector<double> gradRhoBZ(nNumOfGrids, 0.0);
        if ((this->m_bIsUpdateXC == true) && (this->m_nIteration >= 2)) {
            rhoA = gdm.getData(nAtom, GridDataManager::DENSITY_ALPHA);
            rhoB = gdm.getData(nAtom, GridDataManager::DENSITY_BETA);
            gradRhoAX = gdm.getData(nAtom, GridDataManager::GRAD_DENSITY_X_ALPHA);
            gradRhoAY = gdm.getData(nAtom, GridDataManager::GRAD_DENSITY_Y_ALPHA);
            gradRhoAZ = gdm.getData(nAtom, GridDataManager::GRAD_DENSITY_Z_ALPHA);
            gradRhoBX = gdm.getData(nAtom, GridDataManager::GRAD_DENSITY_X_BETA);
            gradRhoBY = gdm.getData(nAtom, GridDataManager::GRAD_DENSITY_Y_BETA);
            gradRhoBZ = gdm.getData(nAtom, GridDataManager::GRAD_DENSITY_Z_BETA);
        }

#pragma omp parallel for schedule(runtime)
        for (int nGrid = 0; nGrid < nNumOfGrids; ++nGrid) {
            const TlPosition gridPosition = grids[nGrid].position;
            const double dWeight = grids[nGrid].weight;

            // calc phi table
            std::vector<WFGrid> aPhi;
            std::vector<WFGrid> aGradPhiX;
            std::vector<WFGrid> aGradPhiY;
            std::vector<WFGrid> aGradPhiZ;
            this->getPhiTable(gridPosition, aPhi, aGradPhiX, aGradPhiY, aGradPhiZ);

            std::sort(aPhi.begin(), aPhi.end(), WFGrid_sort_functional());
            std::sort(aGradPhiX.begin(), aGradPhiX.end(), WFGrid_sort_functional());
            std::sort(aGradPhiY.begin(), aGradPhiY.end(), WFGrid_sort_functional());
            std::sort(aGradPhiZ.begin(), aGradPhiZ.end(), WFGrid_sort_functional());

            // get rho at grid point
            double dRhoA = 0.0;
            double dRhoB = 0.0;
            double dGradRhoAX = 0.0;
            double dGradRhoAY = 0.0;
            double dGradRhoAZ = 0.0;
            double dGradRhoBX = 0.0;
            double dGradRhoBY = 0.0;
            double dGradRhoBZ = 0.0;
            this->getRhoAtGridPoint(PA,
                                    aPhi, aGradPhiX, aGradPhiY, aGradPhiZ,
                                    &dRhoA, &dGradRhoAX, &dGradRhoAY, &dGradRhoAZ);
            this->getRhoAtGridPoint(PB,
                                    aPhi, aGradPhiX, aGradPhiY, aGradPhiZ,
                                    &dRhoB, &dGradRhoBX, &dGradRhoBY, &dGradRhoBZ);
            dRhoA += rhoA[nGrid];
            dGradRhoAX += gradRhoAX[nGrid];
            dGradRhoAY += gradRhoAY[nGrid];
            dGradRhoAZ += gradRhoAZ[nGrid];
            dRhoB += rhoB[nGrid];
            dGradRhoBX += gradRhoBX[nGrid];
            dGradRhoBY += gradRhoBY[nGrid];
            dGradRhoBZ += gradRhoBZ[nGrid];

            // calc
            if ((dRhoA > densityCutOffValueA) || (dRhoB > densityCutOffValueB)) {
                const double rA = std::max(dRhoA, 0.0);
                const double rB = std::max(dRhoB, 0.0);
                assert(rA >= 0.0);
                assert(rB >= 0.0);

                const double gammaAA =  dGradRhoAX*dGradRhoAX + dGradRhoAY*dGradRhoAY + dGradRhoAZ*dGradRhoAZ;
                const double gammaAB =  dGradRhoAX*dGradRhoBX + dGradRhoAY*dGradRhoBY + dGradRhoAZ*dGradRhoBZ;
                const double gammaBB =  dGradRhoBX*dGradRhoBX + dGradRhoBY*dGradRhoBY + dGradRhoBZ*dGradRhoBZ;

#pragma omp critical (DfCalcGridX_calcXCIntegForFockAndEnergy_GGA_U)
                {
                    this->buildFock(rA, rB,
                                    dGradRhoAX, dGradRhoAY, dGradRhoAZ,
                                    dGradRhoBX, dGradRhoBY, dGradRhoBZ,
                                    aPhi, aGradPhiX, aGradPhiY, aGradPhiZ,
                                    pFunctional, dWeight, pFA, pFB); // UKS code
                    dEnergy += dWeight * pFunctional->getFunctional(rA, rB, gammaAA, gammaAB, gammaBB);
                }
            }

            rhoA[nGrid] = dRhoA;
            gradRhoAX[nGrid] = dGradRhoAX;
            gradRhoAY[nGrid] = dGradRhoAY;
            gradRhoAZ[nGrid] = dGradRhoAZ;
            rhoB[nGrid] = dRhoB;
            gradRhoBX[nGrid] = dGradRhoBX;
            gradRhoBY[nGrid] = dGradRhoBY;
            gradRhoBZ[nGrid] = dGradRhoBZ;
        }

        if (this->m_bIsUpdateXC == true) {
            this->physicalValues_[GridDataManager::DENSITY_ALPHA][nAtom] = rhoA;
            this->physicalValues_[GridDataManager::GRAD_DENSITY_X_ALPHA][nAtom] = gradRhoAX;
            this->physicalValues_[GridDataManager::GRAD_DENSITY_Y_ALPHA][nAtom] = gradRhoAY;
            this->physicalValues_[GridDataManager::GRAD_DENSITY_Z_ALPHA][nAtom] = gradRhoAZ;
            this->physicalValues_[GridDataManager::DENSITY_BETA][nAtom] = rhoB;
            this->physicalValues_[GridDataManager::GRAD_DENSITY_X_BETA][nAtom] = gradRhoBX;
            this->physicalValues_[GridDataManager::GRAD_DENSITY_Y_BETA][nAtom] = gradRhoBY;
            this->physicalValues_[GridDataManager::GRAD_DENSITY_Z_BETA][nAtom] = gradRhoBZ;
        }
    }

    return dEnergy;
}


double DfCalcGridX::getPrefactor(const int nType, const TlPosition& pos)
{
    double prefactor = 1.0;
    switch (nType) {
    case 0:
        //prefactor = 1.0;
        break;
    case 1:
        prefactor = pos.x();
        break;
    case 2:
        prefactor = pos.y();
        break;
    case 3:
        prefactor = pos.z();
        break;
    case 4:
        prefactor = pos.x() * pos.y();
        break;
    case 5:
        prefactor = pos.z() * pos.x();
        break;
    case 6:
        prefactor = pos.y() * pos.z();
        break;
    case 7:
        //prefactor = pos.x() * pos.x() - pos.y() * pos.y();
        prefactor = 0.5 * (pos.x() * pos.x() - pos.y() * pos.y());
        break;
    case 8:
        //prefactor = 3.0 * pos.z() * pos.z() - pos.squareDistanceFrom();
        //prefactor = 2.0*pos.z()*pos.z() - (pos.x()*pos.x() + pos.y()*pos.y());
        //prefactor = INV_SQRT12 * (2.0 * pos.z()*pos.z() - (pos.x()*pos.x() + pos.y()*pos.y()));
        prefactor = INV_SQRT3 * (pos.z()*pos.z() - 0.5*(pos.x()*pos.x() + pos.y()*pos.y()));
        break;
    default:
        std::cout << "Basis Type is Wrong." << std::endl;
        break;
    }

    return prefactor;
}


void DfCalcGridX::getPrefactorForDerivative(const int nType, const double alpha, const TlPosition& pos,
                                            double* pPrefactorX, double* pPrefactorY, double* pPrefactorZ)
{
    assert(pPrefactorX != NULL);
    assert(pPrefactorY != NULL);
    assert(pPrefactorZ != NULL);

    const double alpha2 = 2.0 * alpha;
    
    switch (nType) {
    case 0:
        *pPrefactorX = alpha2 * pos.x();
        *pPrefactorY = alpha2 * pos.y();
        *pPrefactorZ = alpha2 * pos.z();
        break;
    case 1:
        *pPrefactorX = alpha2 * pos.x() * pos.x() -1.0;
        *pPrefactorY = alpha2 * pos.x() * pos.y();
        *pPrefactorZ = alpha2 * pos.x() * pos.z();
        break;
    case 2:
        *pPrefactorX = alpha2 * pos.y() * pos.x();
        *pPrefactorY = alpha2 * pos.y() * pos.y() -1.0;
        *pPrefactorZ = alpha2 * pos.y() * pos.z();
        break;
    case 3:
        *pPrefactorX = alpha2 * pos.z() * pos.x();
        *pPrefactorY = alpha2 * pos.z() * pos.y();
        *pPrefactorZ = alpha2 * pos.z() * pos.z() -1.0;
        break;
    case 4: {
        const double xy = pos.x() * pos.y();
        *pPrefactorX = alpha2 * xy * pos.x() - pos.y(); 
        *pPrefactorY = alpha2 * xy * pos.y() - pos.x();
        *pPrefactorZ = alpha2 * xy * pos.z();
    }
    break;
    case 5: {
        const double xz = pos.x() * pos.z();
        *pPrefactorX = alpha2 * xz * pos.x() - pos.z(); 
        *pPrefactorY = alpha2 * xz * pos.y();
        *pPrefactorZ = alpha2 * xz * pos.z() - pos.x();
    }
    break;
    case 6: {
        const double yz = pos.y() * pos.z();
        *pPrefactorX = alpha2 * yz * pos.x(); 
        *pPrefactorY = alpha2 * yz * pos.y() - pos.z();
        *pPrefactorZ = alpha2 * yz * pos.z() - pos.y();
    }
    break;
    case 7: {
        const double xx = pos.x() * pos.x();
        const double xx_X = alpha2 * xx * pos.x() -2.0 * pos.x();
        const double xx_Y = alpha2 * xx * pos.y();
        const double xx_Z = alpha2 * xx * pos.z();

        const double yy = pos.y() * pos.y();
        const double yy_X = alpha2 * yy * pos.x();
        const double yy_Y = alpha2 * yy * pos.y() -2.0 * pos.y();
        const double yy_Z = alpha2 * yy * pos.z();

        *pPrefactorX = 0.5 * (xx_X - yy_X);
        *pPrefactorY = 0.5 * (xx_Y - yy_Y);
        *pPrefactorZ = 0.5 * (xx_Z - yy_Z);
    }
    break;
    case 8: {
        const double xx = pos.x() * pos.x();
        const double xx_X = alpha2 * xx * pos.x() -2.0 * pos.x();
        const double xx_Y = alpha2 * xx * pos.y();
        const double xx_Z = alpha2 * xx * pos.z();

        const double yy = pos.y() * pos.y();
        const double yy_X = alpha2 * yy * pos.x();
        const double yy_Y = alpha2 * yy * pos.y() -2.0 * pos.y();
        const double yy_Z = alpha2 * yy * pos.z();

        const double zz = pos.z() * pos.z();
        const double zz_X = alpha2 * zz * pos.x();
        const double zz_Y = alpha2 * zz * pos.y();
        const double zz_Z = alpha2 * zz * pos.z() -2.0 * pos.z();

        *pPrefactorX = INV_SQRT3 * (zz_X - 0.5 * (xx_X + yy_X));
        *pPrefactorY = INV_SQRT3 * (zz_Y - 0.5 * (xx_Y + yy_Y));
        *pPrefactorZ = INV_SQRT3 * (zz_Z - 0.5 * (xx_Z + yy_Z));
    }
    break;
    default:
        std::cout << "Basis Type is Wrong." << std::endl;
        break;
    }
}

////////////////////////////////////////////////////////////////////////

void DfCalcGridX::readTable()
{
}

// φの値を求める
// for NSD
void DfCalcGridX::getPhiTable(const TlPosition& gridPosition, std::vector<WFGrid>& aPhi)
{
    this->getPhiTable(gridPosition, 0, this->m_tlOrbInfo.getNumOfOrbitals(), aPhi);
}

void DfCalcGridX::getPhiTable(const TlPosition& gridPosition,
                              const int startOrbIndex, const int endOrbIndex,
                              std::vector<WFGrid>& aPhi)
{
    const double densityCutOffValue = std::min(this->m_densityCutOffValueA, this->m_densityCutOffValueB);

    aPhi.clear();
    aPhi.reserve(endOrbIndex - startOrbIndex);

    // orbital loop
    for (int nOrb = startOrbIndex; nOrb < endOrbIndex; ++nOrb) {
        double dPhi = 0.0;

        const TlPosition pos = gridPosition - this->m_tlOrbInfo.getPosition(nOrb);
        const double distance2 = pos.squareDistanceFrom();
        const int nBasisType = this->m_tlOrbInfo.getBasisType(nOrb);
        const double prefactor = this->getPrefactor(nBasisType, pos);

        const int nContract = this->m_tlOrbInfo.getCgtoContraction(nOrb);
        for (int nPGTO = 0; nPGTO < nContract; ++nPGTO) {
            const double alpha = this->m_tlOrbInfo.getExponent(nOrb, nPGTO);
            const double shoulder = alpha * distance2;

            if (shoulder <= TOOBIG) {
                const double gtmp = this->m_tlOrbInfo.getCoefficient(nOrb, nPGTO) * std::exp(-shoulder);
                dPhi += prefactor * gtmp;
            }
        }

        if (std::fabs(dPhi) > densityCutOffValue) {
            const WFGrid wfGrid(nOrb, dPhi);
            aPhi.push_back(wfGrid);
        }
    }
}

void DfCalcGridX::getPhiTable(const TlPosition& gridPosition, std::vector<WFGrid>& aPhi,
                              std::vector<WFGrid>& aGradPhiX, std::vector<WFGrid>& aGradPhiY,
                              std::vector<WFGrid>& aGradPhiZ)
{
    this->getPhiTable(gridPosition,
                      0, this->m_tlOrbInfo.getNumOfOrbitals(),
                      aPhi, aGradPhiX, aGradPhiY, aGradPhiZ);
}


void DfCalcGridX::getPhiTable(const TlPosition& gridPosition,
                              const int startOrbIndex, const int endOrbIndex,
                              std::vector<WFGrid>& aPhi,
                              std::vector<WFGrid>& aGradPhiX, std::vector<WFGrid>& aGradPhiY,
                              std::vector<WFGrid>& aGradPhiZ)
{
    // initialize
    const double densityCutOffValue = std::min(this->m_densityCutOffValueA, this->m_densityCutOffValueB);

    const int range = endOrbIndex - startOrbIndex;
    aPhi.clear();
    aPhi.reserve(range);
    aGradPhiX.clear();
    aGradPhiX.reserve(range);
    aGradPhiY.clear();
    aGradPhiY.reserve(range);
    aGradPhiZ.clear();
    aGradPhiZ.reserve(range);

    // orbital loop
    for (int nOrb = startOrbIndex; nOrb < endOrbIndex; ++nOrb) {
        double dPhi = 0.0;
        double dGradPhiX = 0.0;
        double dGradPhiY = 0.0;
        double dGradPhiZ = 0.0;

        const TlPosition pos = gridPosition - this->m_tlOrbInfo.getPosition(nOrb);
        const double distance2 = pos.squareDistanceFrom();
        const int nBasisType = this->m_tlOrbInfo.getBasisType(nOrb);
        const double prefactor = this->getPrefactor(nBasisType, pos);

        const int nContract = this->m_tlOrbInfo.getCgtoContraction(nOrb);
        for (int nPGTO = 0; nPGTO < nContract; ++nPGTO) {
            double prefactorX = 0.0;
            double prefactorY = 0.0;
            double prefactorZ = 0.0;
            const double alpha = this->m_tlOrbInfo.getExponent(nOrb, nPGTO);
            const double shoulder = alpha * distance2;

            if (shoulder <= TOOBIG) {
                const double gtmp = this->m_tlOrbInfo.getCoefficient(nOrb, nPGTO) * std::exp(-shoulder);
                dPhi += prefactor * gtmp;
                this->getPrefactorForDerivative(nBasisType, alpha, pos, &prefactorX, &prefactorY, &prefactorZ);
                dGradPhiX += prefactorX * gtmp;
                dGradPhiY += prefactorY * gtmp;
                dGradPhiZ += prefactorZ * gtmp;
            }
        }
        
        if (std::fabs(dPhi) > densityCutOffValue) {
            const WFGrid wfGrid(nOrb, dPhi);
            aPhi.push_back(wfGrid);
        }

        if (std::fabs(dGradPhiX) > densityCutOffValue) {
            const WFGrid wfGrid(nOrb, dGradPhiX);
            aGradPhiX.push_back(wfGrid);
        }
        if (std::fabs(dGradPhiY) > densityCutOffValue) {
            const WFGrid wfGrid(nOrb, dGradPhiY);
            aGradPhiY.push_back(wfGrid);
        }
        if (std::fabs(dGradPhiZ) > densityCutOffValue) {
            const WFGrid wfGrid(nOrb, dGradPhiZ);
            aGradPhiZ.push_back(wfGrid);
        }
    }
}


void DfCalcGridX::getPhiTable(const TlPosition& gridPosition,
                              const std::vector<int>& AO_list,
                              std::vector<WFGrid>* pPhis)
{
    assert(pPhis != NULL);

    // initialize
    const double densityCutOffValue = std::min(this->m_densityCutOffValueA,
                                               this->m_densityCutOffValueB);

    const int max_AO_index = AO_list.size();
    pPhis->clear();
    pPhis->reserve(max_AO_index);

    // orbital loop
    for (int AO_index = 0; AO_index < max_AO_index; ++AO_index) {
        const int AO = AO_list[AO_index];

        double phi = 0.0;
//         double gradPhiX = 0.0;
//         double gradPhiY = 0.0;
//         double gradPhiZ = 0.0;

        const TlPosition pos = gridPosition - this->m_tlOrbInfo.getPosition(AO);
        const double distance2 = pos.squareDistanceFrom();
        const int basisType = this->m_tlOrbInfo.getBasisType(AO);
        const double prefactor = this->getPrefactor(basisType, pos);
        const int contraction = this->m_tlOrbInfo.getCgtoContraction(AO);

        for (int PGTO = 0; PGTO < contraction; ++PGTO) {
            const double alpha = this->m_tlOrbInfo.getExponent(AO, PGTO);
            const double shoulder = alpha * distance2;

            if (shoulder <= TOOBIG) {
                const double g = prefactor
                    * this->m_tlOrbInfo.getCoefficient(AO, PGTO)
                    * std::exp(-shoulder);
                phi += g;
            }
        }
        
        if (fabs(phi) > densityCutOffValue) {
            const WFGrid wfGrid(AO, phi);
            pPhis->push_back(wfGrid);
        }
    }
}


void DfCalcGridX::getPhiTable(const TlPosition& gridPosition,
                              const index_type* pAO_List,
                              const std::size_t AO_ListSize,
                              std::vector<WFGrid>* pPhis,
                              std::vector<WFGrid>* pGradPhiXs,
                              std::vector<WFGrid>* pGradPhiYs,
                              std::vector<WFGrid>* pGradPhiZs)
{
    assert(pPhis != NULL);
    assert(pGradPhiXs != NULL);
    assert(pGradPhiYs != NULL);
    assert(pGradPhiZs != NULL);

    // initialize
    const double densityCutOffValue = std::min(this->m_densityCutOffValueA,
                                               this->m_densityCutOffValueB);

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
        const index_type AO = pAO_List[AO_index];

        double phi = 0.0;
        double gradPhiX = 0.0;
        double gradPhiY = 0.0;
        double gradPhiZ = 0.0;

        const TlPosition pos = gridPosition - this->m_tlOrbInfo.getPosition(AO);
        const double distance2 = pos.squareDistanceFrom();
        const int basisType = this->m_tlOrbInfo.getBasisType(AO);
        const double prefactor = this->getPrefactor(basisType, pos);
        const int contraction = this->m_tlOrbInfo.getCgtoContraction(AO);

        for (int PGTO = 0; PGTO < contraction; ++PGTO) {
            double prefactorX = 0.0;
            double prefactorY = 0.0;
            double prefactorZ = 0.0;
            const double alpha = this->m_tlOrbInfo.getExponent(AO, PGTO);
            const double shoulder = alpha * distance2;

            if (shoulder <= TOOBIG) {
                const double g = prefactor
                    * this->m_tlOrbInfo.getCoefficient(AO, PGTO)
                    * std::exp(-shoulder);
                phi += g;
                this->getPrefactorForDerivative(basisType, alpha, pos,
                                                &prefactorX, &prefactorY, &prefactorZ);
                gradPhiX += prefactorX * g;
                gradPhiY += prefactorY * g;
                gradPhiZ += prefactorZ * g;
            }
        }
        
        if (fabs(phi) > densityCutOffValue) {
            const WFGrid wfGrid(AO, phi);
            pPhis->push_back(wfGrid);
        }

        if (fabs(gradPhiX) > densityCutOffValue) {
            const WFGrid wfGrid(AO, gradPhiX);
            pGradPhiXs->push_back(wfGrid);
        }
        if (fabs(gradPhiY) > densityCutOffValue) {
            const WFGrid wfGrid(AO, gradPhiY);
            pGradPhiYs->push_back(wfGrid);
        }
        if (fabs(gradPhiZ) > densityCutOffValue) {
            const WFGrid wfGrid(AO, gradPhiZ);
            pGradPhiZs->push_back(wfGrid);
        }
    }
}


// get rho at grid point
void DfCalcGridX::getRhoAtGridPoint(const TlMatrixObject& P,
                                    const std::vector<WFGrid>& aPhi,
                                    double* pRhoA)
{
    const double densityCutOffValue = this->m_densityCutOffValueA;
    double dRho = 0.0;

    // 密度行列のループ(index = p)の最大を求める
    double value = std::sqrt(densityCutOffValue);
    std::vector<WFGrid>::const_iterator pEnd = std::upper_bound(aPhi.begin(), aPhi.end(),
                                                                WFGrid(0, value),
                                                                WFGrid_sort_functional());
    const std::size_t max_p = std::distance(aPhi.begin(), pEnd);
    std::size_t max_q = aPhi.size();
    for (std::size_t p = 0; p < max_p; ++p) {
        const std::size_t nOrb_p = aPhi[p].index;
        const double phi_p = aPhi[p].value;
        const double cutValue = std::fabs(densityCutOffValue / phi_p);

        // 高速化
        const TlVector P_row = P.getRowVector(nOrb_p);

        dRho += P_row[nOrb_p] * phi_p * phi_p;

//     std::vector<WFGrid>::const_iterator qEnd = std::upper_bound(aPhi.begin() + p +1, pEnd,
//                                  WFGrid(0, cutValue),
//                                  WFGrid_sort_functional());
        std::vector<WFGrid>::const_iterator qEnd = std::upper_bound(aPhi.begin() + p +1,
                                                                    aPhi.begin() + max_q,
                                                                    WFGrid(0, cutValue),
                                                                    WFGrid_sort_functional());
        max_q = std::distance(aPhi.begin(), qEnd);
        for (std::size_t q = p + 1; q < max_q; ++q) {
            const std::size_t nOrb_q = aPhi[q].index;
            const double phi_q = aPhi[q].value;
            dRho += 2.0 * P_row[nOrb_q] * phi_p * phi_q;
        }
    }

    *pRhoA = dRho;
}


void DfCalcGridX::getRhoAtGridPoint(const TlMatrixObject& P,
                                    const std::vector<WFGrid>& aPhi,
                                    const std::vector<WFGrid>& aGradPhiX,
                                    const std::vector<WFGrid>& aGradPhiY,
                                    const std::vector<WFGrid>& aGradPhiZ,
                                    double* pRhoA,
                                    double* pGradRhoAX,
                                    double* pGradRhoAY,
                                    double* pGradRhoAZ)
{
    const double densityCutOffValue = this->m_densityCutOffValueA;
    double dRho = 0.0;
    double dGradRhoX = 0.0;
    double dGradRhoY = 0.0;
    double dGradRhoZ = 0.0;

    // 密度行列のループ(index = p)の最大を求める
    double value = std::sqrt(densityCutOffValue);
    if (aGradPhiX.size() > 0) {
        value = std::min(value, std::fabs(aGradPhiX[0].value));
    }
    if (aGradPhiY.size() > 0) {
        value = std::min(value, std::fabs(aGradPhiY[0].value));
    }
    if (aGradPhiZ.size() > 0) {
        value = std::min(value, std::fabs(aGradPhiZ[0].value));
    }

    std::size_t max_q = aPhi.size();
    std::size_t max_qx = aGradPhiX.size();
    std::size_t max_qy = aGradPhiY.size();
    std::size_t max_qz = aGradPhiZ.size();

    std::vector<WFGrid>::const_iterator pEnd = std::upper_bound(aPhi.begin(), aPhi.end(),
                                                                WFGrid(0, value),
                                                                WFGrid_sort_functional());
    const std::size_t max_p = std::distance(aPhi.begin(), pEnd);
    for (std::size_t p = 0; p < max_p; ++p) {
        const std::size_t nOrb_p = aPhi[p].index;
        const double phi_p = aPhi[p].value;
        const double cutValue = std::fabs(densityCutOffValue / phi_p);

        // 高速化
        const TlVector P_row = P.getRowVector(nOrb_p);

        dRho += P_row[nOrb_p] * phi_p * phi_p;

//     std::vector<WFGrid>::const_iterator qEnd = std::upper_bound(aPhi.begin() + p +1, pEnd,
//                                  WFGrid(0, cutValue),
//                                  WFGrid_sort_functional());
        std::vector<WFGrid>::const_iterator qEnd = std::upper_bound(aPhi.begin() + p +1,
                                                                    aPhi.begin() + max_q,
                                                                    WFGrid(0, cutValue),
                                                                    WFGrid_sort_functional());
        max_q = std::distance(aPhi.begin(), qEnd);
        for (std::size_t q = p + 1; q < max_q; ++q) {
            const std::size_t nOrb_q = aPhi[q].index;
            const double phi_q = aPhi[q].value;
            dRho += 2.0 * P_row[nOrb_q] * phi_p * phi_q;
        }

//     std::vector<WFGrid>::const_iterator qxEnd = std::upper_bound(aGradPhiX.begin(), aGradPhiX.end(),
//                               WFGrid(0, cutValue),
//                               WFGrid_sort_functional());
        std::vector<WFGrid>::const_iterator qxEnd = std::upper_bound(aGradPhiX.begin(),
                                                                     aGradPhiX.begin() + max_qx,
                                                                     WFGrid(0, cutValue),
                                                                     WFGrid_sort_functional());
        max_qx = std::distance(aGradPhiX.begin(), qxEnd);
        for (std::size_t qx = 0; qx < max_qx; ++qx) {
            const std::size_t nOrb_q = aGradPhiX[qx].index;
            const double phi_q = aGradPhiX[qx].value;
            dGradRhoX += 2.0 * P_row[nOrb_q] * phi_p * phi_q;
        }

//     std::vector<WFGrid>::const_iterator qyEnd = std::upper_bound(aGradPhiY.begin(), aGradPhiY.end(),
//                                   WFGrid(0, cutValue),
//                                   WFGrid_sort_functional());
        std::vector<WFGrid>::const_iterator qyEnd = std::upper_bound(aGradPhiY.begin(),
                                                                     aGradPhiY.begin() + max_qy,
                                                                     WFGrid(0, cutValue),
                                                                     WFGrid_sort_functional());
        max_qy = std::distance(aGradPhiY.begin(), qyEnd);
        for (std::size_t qy = 0; qy < max_qy; ++qy) {
            const std::size_t nOrb_q = aGradPhiY[qy].index;
            const double phi_q = aGradPhiY[qy].value;
            dGradRhoY += 2.0 * P_row[nOrb_q] * phi_p * phi_q;
        }

//     std::vector<WFGrid>::const_iterator qzEnd = std::upper_bound(aGradPhiZ.begin(), aGradPhiZ.end(),
//                               WFGrid(0, cutValue),
//                               WFGrid_sort_functional());
        std::vector<WFGrid>::const_iterator qzEnd = std::upper_bound(aGradPhiZ.begin(),
                                                                     aGradPhiZ.begin() + max_qz,
                                                                     WFGrid(0, cutValue),
                                                                     WFGrid_sort_functional());
        max_qz = std::distance(aGradPhiZ.begin(), qzEnd);
        for (std::size_t qz = 0; qz < max_qz; ++qz) {
            const std::size_t nOrb_q = aGradPhiZ[qz].index;
            const double phi_q = aGradPhiZ[qz].value;
            dGradRhoZ += 2.0 * P_row[nOrb_q] * phi_p * phi_q;
        }
    }

    *pRhoA = dRho;
    *pGradRhoAX = dGradRhoX;
    *pGradRhoAY = dGradRhoY;
    *pGradRhoAZ = dGradRhoZ;
}


////////////////////////////////////////////////////////////////////////////////
// build Fock matrix

// for LDA and RKS
void DfCalcGridX::buildFock(const double dRhoA, const std::vector<WFGrid>& aPhi,
                            DfFunctional_LDA* pFunctional, const double dWeight,
                            TlMatrixObject* pF)
{
    double dRoundF_roundRhoA;
    pFunctional->getDerivativeFunctional(dRhoA, &dRoundF_roundRhoA);

    const double dCoeff1_A = dWeight * dRoundF_roundRhoA;
    const double dCutOffValue = this->m_inputedCutoffThreshold;

    // phi v.s. phi
    if ((aPhi.size() > 0) &&
        (std::fabs(dCoeff1_A * aPhi[0].value * aPhi[0].value) > dCutOffValue)) {
        std::vector<WFGrid>::const_iterator pEnd =
            std::upper_bound(aPhi.begin(), aPhi.end(),
                             WFGrid(0, std::sqrt(std::fabs(dCutOffValue / dCoeff1_A))),
                             WFGrid_sort_functional());
        this->buildFock(aPhi.begin(), pEnd, dCoeff1_A, dCutOffValue, pF);
    }
}

// LDA and UKS
void DfCalcGridX::buildFock(const double dRhoA, const double dRhoB,
                            const std::vector<WFGrid>& aPhi,
                            DfFunctional_LDA* pFunctional, const double dWeight,
                            TlMatrixObject* pFA, TlMatrixObject* pFB)
{
    double dRoundF_roundRhoA;
    double dRoundF_roundRhoB;
    pFunctional->getDerivativeFunctional(dRhoA, dRhoB, &dRoundF_roundRhoA, &dRoundF_roundRhoB);

    const double dCoeff1_A = dWeight * dRoundF_roundRhoA;
    const double dCoeff1_B = dWeight * dRoundF_roundRhoB;
    const double dCutOffValue = this->m_inputedCutoffThreshold;

    if (aPhi.size() > 0) {
        const double aPhi_0 = aPhi[0].value;
        if (std::fabs(dCoeff1_A * aPhi_0 * aPhi_0) > dCutOffValue) {
            std::vector<WFGrid>::const_iterator pEnd =
                std::upper_bound(aPhi.begin(), aPhi.end(),
                                 WFGrid(0, std::sqrt(std::fabs(dCutOffValue/dCoeff1_A))),
                                 WFGrid_sort_functional());
            this->buildFock(aPhi.begin(), pEnd, dCoeff1_A, dCutOffValue, pFA);
        }
        if (std::fabs(dCoeff1_B * aPhi_0 * aPhi_0) > dCutOffValue) {
            std::vector<WFGrid>::const_iterator pEnd =
                std::upper_bound(aPhi.begin(), aPhi.end(),
                                 WFGrid(0, std::sqrt(std::fabs(dCutOffValue/dCoeff1_B))),
                                 WFGrid_sort_functional());
            this->buildFock(aPhi.begin(), pEnd, dCoeff1_B, dCutOffValue, pFB);
        }
    }
}


// for GGA and RKS
void DfCalcGridX::buildFock(const double dRhoA,
                            const double dGradRhoAX, const double dGradRhoAY, const double dGradRhoAZ,
                            const std::vector<WFGrid>& aPhi,
                            const std::vector<WFGrid>& aGradPhiX,
                            const std::vector<WFGrid>& aGradPhiY,
                            const std::vector<WFGrid>& aGradPhiZ,
                            DfFunctional_GGA* pFunctional, const double dWeight,
                            TlMatrixObject *pF)
{

    const double dGammaAA = dGradRhoAX * dGradRhoAX + dGradRhoAY * dGradRhoAY + dGradRhoAZ * dGradRhoAZ;

    double dRoundF_roundRhoA;
    double dRoundF_roundGammaAA, dRoundF_roundGammaAB;
    pFunctional->getDerivativeFunctional(dRhoA, dGammaAA,
                                         &dRoundF_roundRhoA, &dRoundF_roundGammaAA, &dRoundF_roundGammaAB);

    const double dCoeff1_A = dWeight * dRoundF_roundRhoA;
    const double dRoundF_roundGammaAA2 = 2.0 * dRoundF_roundGammaAA;

    // remark! dGradRhoBX = dGradRhoAX
    const double dCoeff2_AX =
        dWeight * (dRoundF_roundGammaAA2 + dRoundF_roundGammaAB) * dGradRhoAX;
    const double dCoeff2_AY =
        dWeight * (dRoundF_roundGammaAA2 + dRoundF_roundGammaAB) * dGradRhoAY;
    const double dCoeff2_AZ = 
        dWeight * (dRoundF_roundGammaAA2 + dRoundF_roundGammaAB) * dGradRhoAZ;
    
    const double dCutOffValue = this->m_inputedCutoffThreshold;

    if (aPhi.size() > 0) {
        const double aPhi_0 = aPhi[0].value;
        // phi v.s. phi
        if (std::fabs(dCoeff1_A * aPhi_0 * aPhi_0) > dCutOffValue) {
            std::vector<WFGrid>::const_iterator pEnd =
                std::upper_bound(aPhi.begin(), aPhi.end(),
                                 WFGrid(0, std::sqrt(std::fabs(dCutOffValue / dCoeff1_A))),
                                 WFGrid_sort_functional());
            this->buildFock(aPhi.begin(), pEnd, dCoeff1_A, dCutOffValue, pF);
        }

        // phi v.s. grad-phi(x)
        if ((aGradPhiX.size() > 0) &&
            (std::fabs(dCoeff2_AX * aPhi_0 * aGradPhiX[0].value) > dCutOffValue)) {
            std::vector<WFGrid>::const_iterator pEnd =
                std::upper_bound(aPhi.begin(), aPhi.end(),
                                 WFGrid(0, dCutOffValue / (dCoeff2_AX * aGradPhiX[0].value)),
                                 WFGrid_sort_functional());

        this->buildFock(aPhi.begin(), pEnd, aGradPhiX.begin(), aGradPhiX.end(),
                        dCoeff2_AX, dCutOffValue, pF);
        this->buildFock(aGradPhiX.begin(), aGradPhiX.end(), aPhi.begin(), pEnd,
                        dCoeff2_AX, dCutOffValue, pF);
        }

        // phi v.s. grad-phi(y)
        if ((aGradPhiY.size() > 0) &&
            (std::fabs(dCoeff2_AY * aPhi_0 * aGradPhiY[0].value) > dCutOffValue)) {
            std::vector<WFGrid>::const_iterator pEnd =
                std::upper_bound(aPhi.begin(), aPhi.end(),
                                 WFGrid(0, dCutOffValue / (dCoeff2_AY * aGradPhiY[0].value)),
                                 WFGrid_sort_functional());
        
            this->buildFock(aPhi.begin(), pEnd, aGradPhiY.begin(), aGradPhiY.end(),
                            dCoeff2_AY, dCutOffValue, pF);
            this->buildFock(aGradPhiY.begin(), aGradPhiY.end(), aPhi.begin(), pEnd,
                            dCoeff2_AY, dCutOffValue, pF);
        }

        // phi v.s. grad-phi(z)
        if ((aGradPhiZ.size() > 0) &&
            (std::fabs(dCoeff2_AZ * aPhi_0 * aGradPhiZ[0].value) > dCutOffValue)) {
            std::vector<WFGrid>::const_iterator pEnd =
                std::upper_bound(aPhi.begin(), aPhi.end(),
                                 WFGrid(0, dCutOffValue / (dCoeff2_AZ * aGradPhiZ[0].value)),
                                 WFGrid_sort_functional());
            
            this->buildFock(aPhi.begin(), pEnd, aGradPhiZ.begin(), aGradPhiZ.end(),
                            dCoeff2_AZ, dCutOffValue, pF);
            this->buildFock(aGradPhiZ.begin(), aGradPhiZ.end(), aPhi.begin(), pEnd,
                            dCoeff2_AZ, dCutOffValue, pF);
        }
    }
}


// for GGA and UKS
void DfCalcGridX::buildFock(const double dRhoA, const double dRhoB,
                            const double dGradRhoAX, const double dGradRhoAY, const double dGradRhoAZ,
                            const double dGradRhoBX, const double dGradRhoBY, const double dGradRhoBZ,
                            const std::vector<WFGrid>& aPhi,
                            const std::vector<WFGrid>& aGradPhiX,
                            const std::vector<WFGrid>& aGradPhiY,
                            const std::vector<WFGrid>& aGradPhiZ,
                            DfFunctional_GGA* pFunctional, const double dWeight,
                            TlMatrixObject* pFA, TlMatrixObject* pFB)
{
    assert(dRhoA >= 0.0);
    assert(dRhoB >= 0.0);
    const double dGammaAA = dGradRhoAX*dGradRhoAX + dGradRhoAY*dGradRhoAY + dGradRhoAZ*dGradRhoAZ;
    const double dGammaAB = dGradRhoAX*dGradRhoBX + dGradRhoAY*dGradRhoBY + dGradRhoAZ*dGradRhoBZ;
    const double dGammaBB = dGradRhoBX*dGradRhoBX + dGradRhoBY*dGradRhoBY + dGradRhoBZ*dGradRhoBZ;

    double dRoundF_roundRhoA, dRoundF_roundRhoB;
    double dRoundF_roundGammaAA, dRoundF_roundGammaAB, dRoundF_roundGammaBB;
    pFunctional->getDerivativeFunctional(dRhoA, dRhoB, dGammaAA, dGammaAB, dGammaBB,
                                         &dRoundF_roundRhoA, &dRoundF_roundRhoB,
                                         &dRoundF_roundGammaAA, &dRoundF_roundGammaAB, &dRoundF_roundGammaBB);

    const double dCoeff1_A = dWeight * dRoundF_roundRhoA;
    const double dCoeff1_B = dWeight * dRoundF_roundRhoB;
    const double dRoundF_roundGammaAA2 = 2.0 * dRoundF_roundGammaAA;
    const double dRoundF_roundGammaBB2 = 2.0 * dRoundF_roundGammaBB;
    const double dCoeff2_AX = dWeight
                              * (dRoundF_roundGammaAA2 * dGradRhoAX + dRoundF_roundGammaAB * dGradRhoAX);
    const double dCoeff2_AY = dWeight
                              * (dRoundF_roundGammaAA2 * dGradRhoAY + dRoundF_roundGammaAB * dGradRhoAY);
    const double dCoeff2_AZ = dWeight
                              * (dRoundF_roundGammaAA2 * dGradRhoAZ + dRoundF_roundGammaAB * dGradRhoAZ);
    const double dCoeff2_BX = dWeight
                              * (dRoundF_roundGammaBB2 * dGradRhoBX + dRoundF_roundGammaAB * dGradRhoBX);
    const double dCoeff2_BY = dWeight
                              * (dRoundF_roundGammaBB2 * dGradRhoBY + dRoundF_roundGammaAB * dGradRhoBY);
    const double dCoeff2_BZ = dWeight
                              * (dRoundF_roundGammaBB2 * dGradRhoBZ + dRoundF_roundGammaAB * dGradRhoBZ);

    const double dCutOffValue = this->m_inputedCutoffThreshold;

    if (aPhi.size() > 0) {
        const double aPhi_0 = aPhi[0].value;
        // phi v.s. phi
        if (std::fabs(dCoeff1_A * aPhi_0 * aPhi_0) > dCutOffValue) {
            std::vector<WFGrid>::const_iterator pEnd = std::upper_bound(aPhi.begin(), aPhi.end(),
                                                                        WFGrid(0, std::sqrt(std::fabs(dCutOffValue / dCoeff1_A))),
                                                                        WFGrid_sort_functional());
            this->buildFock(aPhi.begin(), pEnd, dCoeff1_A, dCutOffValue, pFA);
        }
        if (std::fabs(dCoeff1_B * aPhi_0 * aPhi_0) > dCutOffValue) {
            std::vector<WFGrid>::const_iterator pEnd = std::upper_bound(aPhi.begin(), aPhi.end(),
                                                                        WFGrid(0, std::sqrt(std::fabs(dCutOffValue / dCoeff1_B))),
                                                                        WFGrid_sort_functional());
            this->buildFock(aPhi.begin(), pEnd, dCoeff1_B, dCutOffValue, pFB);
        }

        // phi v.s. grad-phi(x)
        if ((aGradPhiX.size() > 0) &&
                (std::fabs(dCoeff2_AX * aPhi_0 * aGradPhiX[0].value) > dCutOffValue)) {
            std::vector<WFGrid>::const_iterator pEnd = std::upper_bound(aPhi.begin(), aPhi.end(),
                                                                        WFGrid(0, dCutOffValue / (dCoeff2_AX * aGradPhiX[0].value)),
                                                                        WFGrid_sort_functional());
            this->buildFock(aPhi.begin(), pEnd, aGradPhiX.begin(), aGradPhiX.end(),
                            dCoeff2_AX, dCutOffValue, pFA);
        }
        if ((aGradPhiX.size() > 0) &&
                (std::fabs(dCoeff2_BX * aPhi_0 * aGradPhiX[0].value) > dCutOffValue)) {
            std::vector<WFGrid>::const_iterator pEnd = std::upper_bound(aPhi.begin(), aPhi.end(),
                                                                        WFGrid(0, dCutOffValue / (dCoeff2_BX * aGradPhiX[0].value)),
                                                                        WFGrid_sort_functional());
            this->buildFock(aPhi.begin(), pEnd, aGradPhiX.begin(), aGradPhiX.end(),
                            dCoeff2_BX, dCutOffValue, pFB);
        }

        // phi v.s. grad-phi(y)
        if ((aGradPhiY.size() > 0) &&
                (std::fabs(dCoeff2_AY * aPhi_0 * aGradPhiY[0].value) > dCutOffValue)) {
            std::vector<WFGrid>::const_iterator pEnd = std::upper_bound(aPhi.begin(), aPhi.end(),
                                                                        WFGrid(0, dCutOffValue / (dCoeff2_AY * aGradPhiY[0].value)),
                                                                        WFGrid_sort_functional());
            this->buildFock(aPhi.begin(), pEnd, aGradPhiY.begin(), aGradPhiY.end(),
                            dCoeff2_AY, dCutOffValue, pFA);
        }
        if ((aGradPhiY.size() > 0) &&
                (std::fabs(dCoeff2_BY * aPhi_0 * aGradPhiY[0].value) > dCutOffValue)) {
            std::vector<WFGrid>::const_iterator pEnd = std::upper_bound(aPhi.begin(), aPhi.end(),
                                                                        WFGrid(0, dCutOffValue / (dCoeff2_BY * aGradPhiY[0].value)),
                                                                        WFGrid_sort_functional());
            this->buildFock(aPhi.begin(), pEnd, aGradPhiY.begin(), aGradPhiY.end(),
                            dCoeff2_BY, dCutOffValue, pFB);
        }

        // phi v.s. grad-phi(z)
        if ((aGradPhiZ.size() > 0) &&
                (std::fabs(dCoeff2_AZ * aPhi_0 * aGradPhiZ[0].value) > dCutOffValue)) {
            std::vector<WFGrid>::const_iterator pEnd = std::upper_bound(aPhi.begin(), aPhi.end(),
                                                                        WFGrid(0, dCutOffValue / (dCoeff2_AZ * aGradPhiZ[0].value)),
                                                                        WFGrid_sort_functional());
            this->buildFock(aPhi.begin(), pEnd, aGradPhiZ.begin(), aGradPhiZ.end(),
                            dCoeff2_AZ, dCutOffValue, pFA);
        }
        if ((aGradPhiZ.size() > 0) &&
                (std::fabs(dCoeff2_BZ * aPhi_0 * aGradPhiZ[0].value) > dCutOffValue)) {
            std::vector<WFGrid>::const_iterator pEnd = std::upper_bound(aPhi.begin(), aPhi.end(),
                                                                        WFGrid(0, dCutOffValue / (dCoeff2_BZ * aGradPhiZ[0].value)),
                                                                        WFGrid_sort_functional());
            this->buildFock(aPhi.begin(), pEnd, aGradPhiZ.begin(), aGradPhiZ.end(),
                            dCoeff2_BZ, dCutOffValue, pFB);
        }
    }
}


//
// void DfCalcGridX::buildFock(const double rhoA,
//                             const double gradRhoAX,
//                             const double gradRhoAY,
//                             const double gradRhoAZ,
//                             const std::vector<WFGrid>& phis,
//                             const std::vector<WFGrid>& gradPhiAXs,
//                             const std::vector<WFGrid>& gradPhiAYs,
//                             const std::vector<WFGrid>& gradPhiAZs,
//                             DfFunctional_GGA* pFunctional,
//                             const double dWeight,
//                             TlMatrixObject* pF)
// {

//     const double gammaAA =
//           gradRhoAX * gradRhoAX
//         + gradRhoAY * gradRhoAY
//         + gradRhoAZ * gradRhoAZ;
//     double roundF_roundRhoA;
//     double roundF_roundGammaAA, roundF_roundGammaAB;
//     pFunctional->getDerivativeFunctional(rhoA, gammaAA,
//                                          &roundF_roundRhoA,
//                                          &roundF_roundGammaAA,
//                                          &roundF_roundGammaAB);

//     const double coeff1_A = weight * roundF_roundRhoA;
//     const double roundF_roundGammaAA2 = 2.0 * roundF_roundGammaAA;
//     const double coeff2_AX = weight
//         * (  roundF_roundGammaAA2 * gradRhoAX
//            + roundF_roundGammaAB  * gradRhoAX); // dGradRhoBX = dGradRhoAX
//     const double coeff2_AY = weight
//         * (  roundF_roundGammaAA2 * gradRhoAY
//            + roundF_roundGammaAB  * gradRhoAY); // dGradRhoBY = dGradRhoAY
//     const double coeff2_AZ = weight
//         * (  roundF_roundGammaAA2 * gradRhoAZ
//            + roundF_roundGammaAB    gradRhoAZ); // dGradRhoBZ = dGradRhoAZ
//     const double cutOffValue = this->m_inputedCutoffThreshold;

//     if (phis.size() > 0) {
//         const double phi_0 = phis[0].value;
//         // phi v.s. phi
//         if (std::fabs(coeff1_A * phi_0 * phi_0) > cutOffValue) {
//             std::vector<WFGrid>::const_iterator pEnd =
//                 std::upper_bound(phis.begin(), phis.end(),
//                                  WFGrid(0, std::sqrt(std::fabs(cutOffValue / coeff1_A))),
//                                  WFGrid_sort_functional());
//             this->buildFock(phis.begin(), pEnd, coeff1_A, cutOffValue, pF);
//         }

//         // phi v.s. grad-phi(x)
//         if ((gradPhiAXs.size() > 0) &&
//             (std::fabs(coeff2_AX * phi_0 * gradPhiAXs[0].value) > cutOffValue)) {
//             std::vector<WFGrid>::const_iterator pEnd =
//                 std::upper_bound(phis.begin(), phis.end(),
//                                  WFGrid(0, cutOffValue / (coeff2_AX * gradPhiAXs[0].value)),
//                                  WFGrid_sort_functional());
//             this->buildFock(phis.begin(), pEnd, gradPhiAXs.begin(), gradPhiAXs.end(),
//                             coeff2_AX, cutOffValue, pF);
//         }

//         // phi v.s. grad-phi(y)
//         if ((gradPhiAYs.size() > 0) &&
//             (std::fabs(coeff2_AY * phi_0 * gradPhiAYs[0].value) > cutOffValue)) {
//             std::vector<WFGrid>::const_iterator pEnd =
//                 std::upper_bound(phis.begin(), phis.end(),
//                                  WFGrid(0, cutOffValue / (coeff2_AY * gradPhiAYs[0].value)),
//                                  WFGrid_sort_functional());
//             this->buildFock(phis.begin(), pEnd, gradPhiAYs.begin(), gradPhiAYs.end(),
//                             coeff2_AY, cutOffValue, pF);
//         }

//         // phi v.s. grad-phi(z)
//         if ((gradPhiAZs.size() > 0) &&
//             (std::fabs(coeff2_AZ * phi_0 * gradPhiAZs[0].value) > cutOffValue)) {
//             std::vector<WFGrid>::const_iterator pEnd =
//                 std::upper_bound(phis.begin(), phis.end(),
//                                  WFGrid(0, cutOffValue / (coeff2_AZ * gradPhiAZs[0].value)),
//                                  WFGrid_sort_functional());
//             this->buildFock(phis.begin(), pEnd, gradPhiAZs.begin(), gradPhiAZs.end(),
//                             coeff2_AZ, cutOffValue, pF);
//         }
//     }
// }


void DfCalcGridX::buildFock(std::vector<WFGrid>::const_iterator pBegin,
                            std::vector<WFGrid>::const_iterator pEnd,
                            const double coef, const double cutoffValue,
                            TlMatrixObject* pF)
{
    assert(pF != NULL);

    for (std::vector<WFGrid>::const_iterator p = pBegin; p != pEnd; ++p) {
        const int u_index = p->index;
        const double u_value = p->value;
        const double tmp1A = coef * u_value;

        // 対角要素
        pF->add(u_index, u_index, tmp1A * u_value);
//         if (u_index == 4) {
//             std::cerr << TlUtils::format("C4 % f", tmp1A * u_value)
//                       << std::endl;
//         }
//         if (u_index == 9) {
//             std::cerr << TlUtils::format("C9 % f", tmp1A * u_value)
//                       << std::endl;
//         }

        // 対角要素以外
        std::vector<WFGrid>::const_iterator qEnd = std::upper_bound(p +1, pEnd,
                                                                    WFGrid(0, cutoffValue / tmp1A),
                                                                    WFGrid_sort_functional());
        for (std::vector<WFGrid>::const_iterator q = p +1; q != qEnd; ++q) {
            const int v_index = q->index;
            const double v_value = q->value;
            pF->add(u_index, v_index, tmp1A * v_value);

//             if (u_index == 4 && v_index == 4) {
//                 std::cerr << TlUtils::format("A4 % f", tmp1A * v_value)
//                           << std::endl;
//             }
//             if (u_index == 9 && v_index == 9) {
//                 std::cerr << TlUtils::format("A9 % f",  tmp1A * v_value)
//                           << std::endl;
//             }
        }
    }
}


void DfCalcGridX::buildFock(std::vector<WFGrid>::const_iterator pBegin,
                            std::vector<WFGrid>::const_iterator pEnd,
                            std::vector<WFGrid>::const_iterator qBegin,
                            std::vector<WFGrid>::const_iterator qEnd,
                            const double coef, const double cutoffValue,
                            TlMatrixObject* pF)
{
    assert(pF != NULL);

    for (std::vector<WFGrid>::const_iterator p = pBegin; p != pEnd; ++p) {
        const int u_index = p->index;
        const double u_value = p->value;
        const double tmp2AX = coef * u_value;

        for (std::vector<WFGrid>::const_iterator q = qBegin; q != qEnd; ++q) {
            const int v_index = q->index;
            const double v_value = q->value;

             if (u_index >= v_index) {
                 pF->add(u_index, v_index, tmp2AX * v_value);
             }
             // else {
             //     //pF->add(u_index, v_index, tmp2AX * v_value);
             //    pF->add(v_index, u_index, tmp2AX * v_value);
             // }

//             if (u_index == 4 && v_index == 4) {
//                 std::cerr << TlUtils::format("B4 % f, % f, % f, % f, % f",
//                                              tmp2AX * v_value,
//                                              coef,
//                                              u_value,
//                                              v_value,
//                                              pF->get(4, 4))
//                           << std::endl;
//             } else if (u_index == 9 && v_index == 9) {
//                 std::cerr << TlUtils::format("B9 % f, % f, % f, % f, % f",
//                                              tmp2AX * v_value,
//                                              coef,
//                                              u_value,
//                                              v_value,
//                                              pF->get(9, 9))
//                           << std::endl;
//             }
        }
    }
}


void DfCalcGridX::gridDensity(const TlSymmetricMatrix& P,
                              const TlPosition& gridPosition,
                              double* pRho)
{
    // calc phi table
    std::vector<WFGrid> aPhi;
    this->getPhiTable(gridPosition, aPhi);
    std::sort(aPhi.begin(), aPhi.end(), WFGrid_sort_functional());

    // get rho at grid point
    //double dRhoA;
    this->getRhoAtGridPoint(P, aPhi, pRho);
}


void DfCalcGridX::calcForceFromXC(const TlSymmetricMatrix& P,
                                  DfFunctional_GGA* pFunctional,
                                  TlVector* pFx, TlVector* pFy, TlVector* pFz)
{
    const int numOfAtoms = this->numOfRealAtoms_;
    const int numOfAOs = this->m_nNumOfAOs;

    // grid derivative
    // to implement!
    
    // functional derivative
    TlMatrix GX(numOfAOs, numOfAtoms);
    TlMatrix GY(numOfAOs, numOfAtoms);
    TlMatrix GZ(numOfAOs, numOfAtoms);
    this->makeGammaMatrix(P, pFunctional,
                          &GX, &GY, &GZ);

    for (int atomIndex = 0; atomIndex < numOfAtoms; ++atomIndex) {
        double term1X = 0.0;
        double term1Y = 0.0;
        double term1Z = 0.0;
        double term2X = 0.0;
        double term2Y = 0.0;
        double term2Z = 0.0;

        for (int aoIndex = 0; aoIndex < numOfAOs; ++aoIndex) {
            if (this->m_tlOrbInfo.getAtomIndex(aoIndex) != atomIndex) {
                term1X += GX.get(atomIndex, aoIndex);
                term1Y += GY.get(atomIndex, aoIndex);
                term1Z += GZ.get(atomIndex, aoIndex);
            } else {
                for (int atomIndex2 = 0; atomIndex2 < numOfAtoms; ++atomIndex2) {
                    if (atomIndex != atomIndex2) {
                        term2X += GX.get(atomIndex2, aoIndex);
                        term2Y += GY.get(atomIndex2, aoIndex);
                        term2Z += GZ.get(atomIndex2, aoIndex);
                    }
                }
            }
        }

        pFx->set(atomIndex, (term1X - term2X));
        pFy->set(atomIndex, (term1Y - term2Y));
        pFz->set(atomIndex, (term1Z - term2Z));
    }
}


/// Brown et al., J. Comput. Chem., 31, 2008 (2010).
/// EQ.13 参照
/// Gamma行列の次元: (AO, atoms)
void DfCalcGridX::makeGammaMatrix(const TlSymmetricMatrix& P,
                                  DfFunctional_LDA* pFunctional,
                                  TlMatrix* pGX, TlMatrix* pGY, TlMatrix* pGZ)
{
    assert(pFunctional != NULL);
    assert(pGX != NULL);
    assert(pGY != NULL);
    assert(pGZ != NULL);

    const GridDataManager gdm(this->getGridDataFilePath());

    const int numOfAtoms = this->numOfRealAtoms_;
    const int numOfAOs = this->m_nNumOfAOs;
    
    pGX->resize(numOfAOs, numOfAtoms);
    pGY->resize(numOfAOs, numOfAtoms);
    pGZ->resize(numOfAOs, numOfAtoms);

    for (int atomIndex = 0; atomIndex < numOfAtoms; ++atomIndex) {
        // get grid information
        std::vector<GridDataManager::GridInfo> grids;
        {
            const std::vector<double> coordX = gdm.getData(atomIndex, GridDataManager::COORD_X);
            const std::vector<double> coordY = gdm.getData(atomIndex, GridDataManager::COORD_Y);
            const std::vector<double> coordZ = gdm.getData(atomIndex, GridDataManager::COORD_Z);
            const std::vector<double> weight = gdm.getData(atomIndex, GridDataManager::GRID_WEIGHT);
            grids = GridDataManager::composeGridInfo(coordX, coordY, coordZ, weight);
        }
        const int numOfGrids = grids.size();
    
        std::vector<double> rhoA(numOfGrids, 0.0);
        std::vector<double> gradRhoAX(numOfGrids, 0.0);
        std::vector<double> gradRhoAY(numOfGrids, 0.0);
        std::vector<double> gradRhoAZ(numOfGrids, 0.0);
        // 電子密度を計算する
        // LDAではSCF中に微分を計算しないため
        for (int gridIndex = 0; gridIndex < numOfGrids; ++gridIndex) {
            const TlPosition gridPosition = grids[gridIndex].position;
            
            // calc phi table
            std::vector<WFGrid> phis;
            std::vector<WFGrid> gradPhisX;
            std::vector<WFGrid> gradPhisY;
            std::vector<WFGrid> gradPhisZ;
            this->getPhiTable(gridPosition, phis, gradPhisX, gradPhisY, gradPhisZ);
            std::sort(phis.begin(), phis.end(), WFGrid_sort_functional());
            std::sort(gradPhisX.begin(), gradPhisX.end(), WFGrid_sort_functional());
            std::sort(gradPhisY.begin(), gradPhisY.end(), WFGrid_sort_functional());
            std::sort(gradPhisZ.begin(), gradPhisZ.end(), WFGrid_sort_functional());
            
            // get rho at grid point
            double gridRhoA = 0.0;
            double gridGradRhoAX = 0.0;
            double gridGradRhoAY = 0.0;
            double gridGradRhoAZ = 0.0;
            this->getRhoAtGridPoint(P, phis, gradPhisX, gradPhisY, gradPhisZ,
                                    &gridRhoA, &gridGradRhoAX, &gridGradRhoAY, &gridGradRhoAZ);
            gridRhoA *= 0.5;
            gridGradRhoAX *= 0.5;
            gridGradRhoAY *= 0.5;
            gridGradRhoAZ *= 0.5;
            
            rhoA[gridIndex] = gridRhoA;
            gradRhoAX[gridIndex] = gridGradRhoAX;
            gradRhoAY[gridIndex] = gridGradRhoAY;
            gradRhoAZ[gridIndex] = gridGradRhoAZ;
        }

        TlMatrix Gx(numOfAOs, numOfAOs);
        TlMatrix Gy(numOfAOs, numOfAOs);
        TlMatrix Gz(numOfAOs, numOfAOs);
        
        std::vector<WFGrid> etas(numOfGrids);
        std::vector<WFGrid> gradEtasX(numOfGrids);
        std::vector<WFGrid> gradEtasY(numOfGrids);
        std::vector<WFGrid> gradEtasZ(numOfGrids);
        for (int gridIndex = 0; gridIndex < numOfGrids; ++gridIndex) {
            const TlPosition gridPosition = grids[gridIndex].position;
            this->getPhiTable(gridPosition, etas, gradEtasX, gradEtasY, gradEtasZ);
            std::sort(etas.begin(), etas.end(), WFGrid_sort_functional());
            std::sort(gradEtasX.begin(), gradEtasX.end(), WFGrid_sort_functional());
            std::sort(gradEtasY.begin(), gradEtasY.end(), WFGrid_sort_functional());
            std::sort(gradEtasZ.begin(), gradEtasZ.end(), WFGrid_sort_functional());
            
            const double weight = grids[gridIndex].weight;
            const double gridRhoA = rhoA[gridIndex];
            // const double gridGradRhoAX = gradRhoAX[gridIndex];
            // const double gridGradRhoAY = gradRhoAY[gridIndex];
            // const double gridGradRhoAZ = gradRhoAZ[gridIndex];
            // const double gridGammaAA = gridGradRhoAX*gridGradRhoAX
            //     + gridGradRhoAY*gridGradRhoAY
            //     + gridGradRhoAZ*gridGradRhoAZ;
            double roundF_roundRhoA;
            // double roundF_roundGammaAA, roundF_roundGammaAB;
            // pFunctional->getDerivativeFunctional(gridRhoA, gridGammaAA,
            //                                      &roundF_roundRhoA, &roundF_roundGammaAA, &roundF_roundGammaAB);
            pFunctional->getDerivativeFunctional(gridRhoA, &roundF_roundRhoA);
            
            const double coef = weight * roundF_roundRhoA;
            const int numOfEtas = etas.size();
            const int numOfGradEtasX = gradEtasX.size();
            const int numOfGradEtasY = gradEtasY.size();
            const int numOfGradEtasZ = gradEtasZ.size();
            for (int i = 0; i < numOfEtas; ++i) {
                const double eta = etas[i].value;
                const index_type q = etas[i].index;
                const double coef_eta = coef * eta;
                
                for (int j = 0; j < numOfGradEtasX; ++j) {
                    const double etaDash = gradEtasX[j].value;
                    const index_type p = gradEtasX[j].index;
                    
                    const double value = coef_eta * etaDash;
                    Gx.add(p, q, value);
                }
                for (int j = 0; j < numOfGradEtasY; ++j) {
                    const double etaDash = gradEtasY[j].value;
                    const index_type p = gradEtasY[j].index;
                    
                    const double value = coef_eta * etaDash;
                    Gy.add(p, q, value);
                }
                for (int j = 0; j < numOfGradEtasZ; ++j) {
                    const double etaDash = gradEtasZ[j].value;
                    const index_type p = gradEtasZ[j].index;
                    
                    const double value = coef_eta * etaDash;
                    Gz.add(p, q, value);
                }
            }
        }
        
        Gx.dot(P);
        Gy.dot(P);
        Gz.dot(P);
        
        for (index_type p = 0; p < numOfAOs; ++p) {
            const TlVector Gpx = Gx.getRowVector(p);
            const TlVector Gpy = Gy.getRowVector(p);
            const TlVector Gpz = Gz.getRowVector(p);
            
            pGX->set(p, atomIndex, 2.0 * Gpx.sum());
            pGY->set(p, atomIndex, 2.0 * Gpy.sum());
            pGZ->set(p, atomIndex, 2.0 * Gpz.sum());
        }
    }
}


void DfCalcGridX::makeGammaMatrix(const TlSymmetricMatrix& P,
                                  DfFunctional_GGA* pFunctional,
                                  TlMatrix* pGX, TlMatrix* pGY, TlMatrix* pGZ)
{
    assert(pFunctional != NULL);
    assert(pGX != NULL);
    assert(pGY != NULL);
    assert(pGZ != NULL);

    const double densityCutOffValue = this->m_densityCutOffValueA;
    const GridDataManager gdm(this->getGridDataFilePath());

    const int numOfAtoms = this->numOfRealAtoms_;
    const int numOfAOs = this->m_nNumOfAOs;
    
    pGX->resize(numOfAOs, numOfAtoms);
    pGY->resize(numOfAOs, numOfAtoms);
    pGZ->resize(numOfAOs, numOfAtoms);

    for (int atomIndex = 0; atomIndex < numOfAtoms; ++atomIndex) {
        // get grid information
        std::vector<GridDataManager::GridInfo> grids;
        {
            const std::vector<double> coordX = gdm.getData(atomIndex, GridDataManager::COORD_X);
            const std::vector<double> coordY = gdm.getData(atomIndex, GridDataManager::COORD_Y);
            const std::vector<double> coordZ = gdm.getData(atomIndex, GridDataManager::COORD_Z);
            const std::vector<double> weight = gdm.getData(atomIndex, GridDataManager::GRID_WEIGHT);
            grids = GridDataManager::composeGridInfo(coordX, coordY, coordZ, weight);
        }
        const int numOfGrids = grids.size();
    
        std::vector<double> rhoA(numOfGrids, 0.0);
        std::vector<double> gradRhoAX(numOfGrids, 0.0);
        std::vector<double> gradRhoAY(numOfGrids, 0.0);
        std::vector<double> gradRhoAZ(numOfGrids, 0.0);
        if (this->m_bIsUpdateXC == true) {
            // すでにSCF計算中に電子密度を計算済み
            rhoA = gdm.getData(atomIndex, GridDataManager::DENSITY);
            gradRhoAX = gdm.getData(atomIndex, GridDataManager::GRAD_DENSITY_X);
            gradRhoAY = gdm.getData(atomIndex, GridDataManager::GRAD_DENSITY_Y);
            gradRhoAZ = gdm.getData(atomIndex, GridDataManager::GRAD_DENSITY_Z);
        } else {
            // 電子密度を計算し直す
            for (int gridIndex = 0; gridIndex < numOfGrids; ++gridIndex) {
                const TlPosition gridPosition = grids[gridIndex].position;
                
                // calc phi table
                std::vector<WFGrid> phis;
                std::vector<WFGrid> gradPhisX;
                std::vector<WFGrid> gradPhisY;
                std::vector<WFGrid> gradPhisZ;
                this->getPhiTable(gridPosition, phis, gradPhisX, gradPhisY, gradPhisZ);
                std::sort(phis.begin(), phis.end(), WFGrid_sort_functional());
                std::sort(gradPhisX.begin(), gradPhisX.end(), WFGrid_sort_functional());
                std::sort(gradPhisY.begin(), gradPhisY.end(), WFGrid_sort_functional());
                std::sort(gradPhisZ.begin(), gradPhisZ.end(), WFGrid_sort_functional());
                
                // get rho at grid point
                double gridRhoA = 0.0;
                double gridGradRhoAX = 0.0;
                double gridGradRhoAY = 0.0;
                double gridGradRhoAZ = 0.0;
                this->getRhoAtGridPoint(P, phis, gradPhisX, gradPhisY, gradPhisZ,
                                        &gridRhoA, &gridGradRhoAX, &gridGradRhoAY, &gridGradRhoAZ);
                gridRhoA *= 0.5;
                gridGradRhoAX *= 0.5;
                gridGradRhoAY *= 0.5;
                gridGradRhoAZ *= 0.5;
                
                rhoA[gridIndex] = gridRhoA;
                gradRhoAX[gridIndex] = gridGradRhoAX;
                gradRhoAY[gridIndex] = gridGradRhoAY;
                gradRhoAZ[gridIndex] = gridGradRhoAZ;
            }
        }

        TlMatrix Gx(numOfAOs, numOfAOs);
        TlMatrix Gy(numOfAOs, numOfAOs);
        TlMatrix Gz(numOfAOs, numOfAOs);
        
        std::vector<WFGrid> etas(numOfGrids);
        std::vector<WFGrid> gradEtasX(numOfGrids);
        std::vector<WFGrid> gradEtasY(numOfGrids);
        std::vector<WFGrid> gradEtasZ(numOfGrids);
        for (int gridIndex = 0; gridIndex < numOfGrids; ++gridIndex) {
            const TlPosition gridPosition = grids[gridIndex].position;
            this->getPhiTable(gridPosition, etas, gradEtasX, gradEtasY, gradEtasZ);
            std::sort(etas.begin(), etas.end(), WFGrid_sort_functional());
            std::sort(gradEtasX.begin(), gradEtasX.end(), WFGrid_sort_functional());
            std::sort(gradEtasY.begin(), gradEtasY.end(), WFGrid_sort_functional());
            std::sort(gradEtasZ.begin(), gradEtasZ.end(), WFGrid_sort_functional());
            
            const double weight = grids[gridIndex].weight;
            const double gridRhoA = rhoA[gridIndex];
            const double gridGradRhoAX = gradRhoAX[gridIndex];
            const double gridGradRhoAY = gradRhoAY[gridIndex];
            const double gridGradRhoAZ = gradRhoAZ[gridIndex];
            const double gridGammaAA = gridGradRhoAX*gridGradRhoAX
                + gridGradRhoAY*gridGradRhoAY
                + gridGradRhoAZ*gridGradRhoAZ;
            double roundF_roundRhoA;
            double roundF_roundGammaAA, roundF_roundGammaAB;

            
            if (gridRhoA > densityCutOffValue) {
                pFunctional->getDerivativeFunctional(gridRhoA, gridGammaAA,
                                                     &roundF_roundRhoA, &roundF_roundGammaAA, &roundF_roundGammaAB);
                
                const double coef = weight * roundF_roundRhoA;
                const int numOfEtas = etas.size();
                const int numOfGradEtasX = gradEtasX.size();
                const int numOfGradEtasY = gradEtasY.size();
                const int numOfGradEtasZ = gradEtasZ.size();
                for (int i = 0; i < numOfEtas; ++i) {
                    const double eta = etas[i].value;
                    const index_type q = etas[i].index;
                    const double coef_eta = coef * eta;
                    
                    for (int j = 0; j < numOfGradEtasX; ++j) {
                        const double etaDash = gradEtasX[j].value;
                        const index_type p = gradEtasX[j].index;
                        
                        const double value = coef_eta * etaDash;
                        Gx.add(p, q, value);
                    }
                    for (int j = 0; j < numOfGradEtasY; ++j) {
                        const double etaDash = gradEtasY[j].value;
                        const index_type p = gradEtasY[j].index;
                        
                        const double value = coef_eta * etaDash;
                        Gy.add(p, q, value);
                    }
                    for (int j = 0; j < numOfGradEtasZ; ++j) {
                        const double etaDash = gradEtasZ[j].value;
                        const index_type p = gradEtasZ[j].index;
                        
                        const double value = coef_eta * etaDash;
                        Gz.add(p, q, value);
                    }
                }
            }
        }
        
        Gx.dot(P);
        Gy.dot(P);
        Gz.dot(P);
        
        for (index_type p = 0; p < numOfAOs; ++p) {
            const TlVector Gpx = Gx.getRowVector(p);
            const TlVector Gpy = Gy.getRowVector(p);
            const TlVector Gpz = Gz.getRowVector(p);
            
            pGX->set(p, atomIndex, 2.0 * Gpx.sum());
            pGY->set(p, atomIndex, 2.0 * Gpy.sum());
            pGZ->set(p, atomIndex, 2.0 * Gpz.sum());
        }
    }
}


// experimental code
double DfCalcGridX::calcXCIntegForFockAndEnergy2(const TlSymmetricMatrix& P,
                                                 DfFunctional_LDA* pFunctional,
                                                 TlSymmetricMatrix* pF)
{
    assert(pFunctional != NULL);
    assert(pF != NULL);

    this->backupGridData();

    // setup
    //const int numOfAtoms = this->numOfRealAtoms_;
    this->physicalValues_.clear();
    this->defineCutOffValues(P);

    const index_type numOfAOs = this->m_nNumOfAOs;
    std::vector<index_type> AO_indexes(numOfAOs);
    for (index_type i = 0; i < numOfAOs; ++i) {
        AO_indexes[i] = i;
    }

    TlMatrix gridMatrix = this->getGridMatrix<TlMatrix>();
    gridMatrix.resize(gridMatrix.getNumOfRows(), 5);
    
    this->loggerTime(" calculate population of each grid");
    this->calcRhoVals_LDA(AO_indexes, AO_indexes, 0.5 * TlMatrix(P),
                          &gridMatrix);
    

    this->loggerTime(" build the KS matrix of XC term");
    // const double energy = this->buildK(gridMatrix,
    //                                    pFunctional, pF);
    const double energy = this->buildK_2(gridMatrix,
                                         AO_indexes, AO_indexes,
                                         pFunctional, pF);

    if (this->m_bIsUpdateXC == true) {
        this->loggerTime(" save grid data");
        this->saveGridMatrix(gridMatrix);
    }
    this->logger("\n");

    return energy;
}


double DfCalcGridX::calcXCIntegForFockAndEnergy2(const TlSymmetricMatrix& P,
                                                 DfFunctional_GGA* pFunctional,
                                                 TlSymmetricMatrix* pF)
{
    assert(pFunctional != NULL);
    assert(pF != NULL);

    this->backupGridData();

    // setup
    this->physicalValues_.clear();
    this->defineCutOffValues(P);

    const index_type numOfAOs = this->m_nNumOfAOs;
    std::vector<index_type> AO_indexes(numOfAOs);
    for (index_type i = 0; i < numOfAOs; ++i) {
        AO_indexes[i] = i;
    }

    TlMatrix gridMatrix = this->getGridMatrix<TlMatrix>();
    gridMatrix.resize(gridMatrix.getNumOfRows(), 8);
    
    this->loggerTime(" calculate population of each grid");
    this->calcRhoVals_GGA(AO_indexes, AO_indexes, 0.5 * TlMatrix(P),
                          &gridMatrix);
    

    this->loggerTime(" build the KS matrix of XC term");
    TlMatrix tmpF(numOfAOs, numOfAOs);
    const double energy = this->buildK_2(gridMatrix,
                                         AO_indexes, AO_indexes,
                                         pFunctional, &tmpF);
    *pF = tmpF;

    if (this->m_bIsUpdateXC == true) {
        this->loggerTime(" save grid data");
        this->saveGridMatrix(gridMatrix);
    }
    this->loggerTime("\n");

    return energy;
}


void DfCalcGridX::calcRhoVals_LDA(const std::vector<index_type>& P_rowIndexes,
                                  const std::vector<index_type>& P_colIndexes,
                                  const TlMatrix& P,
                                  TlMatrix* pGridMatrix)
{
    assert(pGridMatrix != NULL);
    
    const std::size_t numOfRows = P.getNumOfRows();
    const std::size_t numOfCols = P.getNumOfCols();
    assert(P_rowIndexes.size() == numOfRows);
    assert(P_colIndexes.size() == numOfCols);
    const int calcMode = pGridMatrix->getNumOfCols();
    
    const std::size_t numOfGrids = pGridMatrix->getNumOfRows();
#pragma omp parallel for schedule(runtime)
    for (std::size_t grid = 0; grid < numOfGrids; ++grid) {
        const double x = pGridMatrix->get(grid, 0);
        const double y = pGridMatrix->get(grid, 1);
        const double z = pGridMatrix->get(grid, 2);
        const TlPosition r(x, y, z);
        
        TlMatrix wf_row;
        this->getWaveFunctionValues(P_rowIndexes, r, &wf_row);
        {
            TlMatrix wf_col;
            this->getWaveFunctionValues(P_colIndexes, r, &wf_col);
            wf_col.transpose();
            TlMatrix coefMatrix = wf_row * wf_col;
            assert(coefMatrix.getNumOfRows() == numOfRows);
            assert(coefMatrix.getNumOfCols() == numOfCols);
            coefMatrix.dot(P);
            pGridMatrix->add(grid, 4, coefMatrix.sum());
        }
    }
}


void DfCalcGridX::calcRhoVals_GGA(const std::vector<index_type>& P_rowIndexes,
                                  const std::vector<index_type>& P_colIndexes,
                                  const TlMatrix& P,
                                  TlMatrix* pGridMatrix)
{
    assert(pGridMatrix != NULL);
    
    const std::size_t numOfRows = P.getNumOfRows();
    const std::size_t numOfCols = P.getNumOfCols();
    assert(P_rowIndexes.size() == numOfRows);
    assert(P_colIndexes.size() == numOfCols);
    const int calcMode = pGridMatrix->getNumOfCols();
    
    const std::size_t numOfGrids = pGridMatrix->getNumOfRows();
#pragma omp parallel for schedule(runtime)
    for (std::size_t grid = 0; grid < numOfGrids; ++grid) {
        const double x = pGridMatrix->get(grid, 0);
        const double y = pGridMatrix->get(grid, 1);
        const double z = pGridMatrix->get(grid, 2);
        const TlPosition r(x, y, z);
        
        TlMatrix wf_row;
        this->getWaveFunctionValues(P_rowIndexes, r, &wf_row);

        TlMatrix wf_col, wf_dX, wf_dY, wf_dZ;
        this->getWaveFunctionValues(P_colIndexes, r,
                                        &wf_col, &wf_dX, &wf_dY, &wf_dZ);
        wf_col.transpose();
        wf_dX.transpose();
        wf_dY.transpose();
        wf_dZ.transpose();
        {
            TlMatrix wf_rc = wf_row * wf_col;
            assert(wf_rc.getNumOfRows() == numOfRows);
            assert(wf_rc.getNumOfCols() == numOfCols);
            wf_rc.dot(P);
            pGridMatrix->add(grid, 4, wf_rc.sum());
        }
        {
            TlMatrix wf_rc = wf_row * wf_dX;
            assert(wf_rc.getNumOfRows() == numOfRows);
            assert(wf_rc.getNumOfCols() == numOfCols);
            wf_rc.dot(P);
            pGridMatrix->add(grid, 5, 2.0 * wf_rc.sum());
        }
        {
            TlMatrix wf_rc = wf_row * wf_dY;
            assert(wf_rc.getNumOfRows() == numOfRows);
            assert(wf_rc.getNumOfCols() == numOfCols);
            wf_rc.dot(P);
            pGridMatrix->add(grid, 6, 2.0 * wf_rc.sum());
        }
        {
            TlMatrix wf_rc = wf_row * wf_dZ;
            assert(wf_rc.getNumOfRows() == numOfRows);
            assert(wf_rc.getNumOfCols() == numOfCols);
            wf_rc.dot(P);
            pGridMatrix->add(grid, 7, 2.0 * wf_rc.sum());
        }
    }
}


void DfCalcGridX::getWaveFunctionValues(const std::vector<index_type>& AO_indexes,
                                        const TlPosition& gridPosition,
                                        TlMatrix* pWF)
{
    assert(pWF != NULL);
    const index_type size = AO_indexes.size();
    pWF->resize(size, 1);
    for (index_type aoIndex = 0; aoIndex < size; ++aoIndex) {
        const index_type orb = AO_indexes[aoIndex];

        const TlPosition r = gridPosition - this->m_tlOrbInfo.getPosition(orb);
        const double r2 = r.squareDistanceFrom();
        const int basisType = this->m_tlOrbInfo.getBasisType(orb);
        const double prefactor = this->getPrefactor(basisType, r);

        double value = 0.0;
        const int numOfPGTOs = this->m_tlOrbInfo.getCgtoContraction(orb);
        for (int pgtoIndex = 0; pgtoIndex < numOfPGTOs; ++pgtoIndex) {
            const double alpha = this->m_tlOrbInfo.getExponent(orb, pgtoIndex);
            const double exponent = alpha * r2;
            const double coef = this->m_tlOrbInfo.getCoefficient(orb, pgtoIndex);
            value += coef * std::exp(- exponent);
        }

        pWF->set(aoIndex, 0, prefactor * value);
    }
}


void DfCalcGridX::getWaveFunctionValues(const std::vector<index_type>& AO_indexes,
                                        const TlPosition& gridPosition,
                                        TlMatrix* pWF,
                                        TlMatrix* pGradWF_X,
                                        TlMatrix* pGradWF_Y,
                                        TlMatrix* pGradWF_Z)
{
    const index_type size = AO_indexes.size();
    pWF->resize(size, 1);
    pGradWF_X->resize(size, 1);
    pGradWF_Y->resize(size, 1);
    pGradWF_Z->resize(size, 1);

    double prefactorX = 0.0;
    double prefactorY = 0.0;
    double prefactorZ = 0.0;
    for (index_type aoIndex = 0; aoIndex < size; ++aoIndex) {
        const index_type orb = AO_indexes[aoIndex];

        const TlPosition r = gridPosition - this->m_tlOrbInfo.getPosition(orb);
        const double r2 = r.squareDistanceFrom();
        const int basisType = this->m_tlOrbInfo.getBasisType(orb);
        const double prefactor = this->getPrefactor(basisType, r);

        double value = 0.0;
        double gradX = 0.0;
        double gradY = 0.0;
        double gradZ = 0.0;
        const int numOfPGTOs = this->m_tlOrbInfo.getCgtoContraction(orb);
        for (int pgtoIndex = 0; pgtoIndex < numOfPGTOs; ++pgtoIndex) {
            const double alpha = this->m_tlOrbInfo.getExponent(orb, pgtoIndex);
            const double exponent = alpha * r2;
            const double e = std::exp(- exponent);
            const double coef = this->m_tlOrbInfo.getCoefficient(orb, pgtoIndex);
            this->getPrefactorForDerivative(basisType, alpha, r,
                                            &prefactorX, &prefactorY, &prefactorZ);

            value += coef * e;
            gradX += coef * prefactorX * e;
            gradY += coef * prefactorY * e;
            gradZ += coef * prefactorZ * e;
        }
        pWF->set(aoIndex, 0, prefactor* value);
        pGradWF_X->set(aoIndex, 0, gradX);
        pGradWF_Y->set(aoIndex, 0, gradY);
        pGradWF_Z->set(aoIndex, 0, gradZ);
    }
}


double DfCalcGridX::buildK(const TlMatrix& gridMatrix,
                           DfFunctional_LDA* pFunctional,
                           TlSymmetricMatrix* pF)
{
    const index_type numOfAOs = this->m_nNumOfAOs;
    const double densityCutOffValue = this->m_densityCutOffValueA;

    double energy = 0.0;
    const std::size_t numOfAllGrids = gridMatrix.getNumOfRows();
    for (std::size_t grid = 0; grid < numOfAllGrids; ++grid) {
        const double x = gridMatrix.get(grid, 0);
        const double y = gridMatrix.get(grid, 1);
        const double z = gridMatrix.get(grid, 2);
        const double w = gridMatrix.get(grid, 3);
        const TlPosition position(x, y, z);
        const double rhoA = gridMatrix.get(grid, 4);

        // calc phi table
        std::vector<WFGrid> phis;
        this->getPhiTable(position, 0, numOfAOs, phis);
        std::sort(phis.begin(), phis.end(), WFGrid_sort_functional());

        if (rhoA > densityCutOffValue) {
            this->buildFock(rhoA,
                            phis,
                            pFunctional, w, pF);
            energy += w * pFunctional->getFunctional(rhoA);
        }
    }

    return energy;
}


double DfCalcGridX::buildK_2(const TlMatrix& gridMatrix,
                             const std::vector<index_type>& rowIndexes,
                             const std::vector<index_type>& colIndexes,
                             DfFunctional_LDA* pFunctional,
                             TlSymmetricMatrix* pF)
{
    const double densityCutOffValue = this->m_densityCutOffValueA;

    TlMatrix tmpF(*pF);
    double energy = 0.0;
    const std::size_t numOfAllGrids = gridMatrix.getNumOfRows();
    for (std::size_t grid = 0; grid < numOfAllGrids; ++grid) {
        const double x = gridMatrix.get(grid, 0);
        const double y = gridMatrix.get(grid, 1);
        const double z = gridMatrix.get(grid, 2);
        const double w = gridMatrix.get(grid, 3);
        const TlPosition r(x, y, z);
        const double rhoA = gridMatrix.get(grid, 4);

        if (rhoA > densityCutOffValue) {
            TlMatrix wf_row;
            this->getWaveFunctionValues(rowIndexes, r,
                                        &wf_row);
            TlMatrix wf_col;
            this->getWaveFunctionValues(colIndexes, r,
                                        &wf_col);
            wf_col.transpose();
            
            // build K
            double roundF_roundRhoA;
            pFunctional->getDerivativeFunctional(rhoA, &roundF_roundRhoA);
            
            const double coef1_A = w * roundF_roundRhoA;
            {
                TlMatrix wf_rc = wf_row * wf_col;
                wf_rc *= coef1_A;
                tmpF += wf_rc;
            }
            
            // energy
            energy += w * pFunctional->getFunctional(rhoA); // RKS code
        }
    }

    *pF = tmpF;

    return energy;
}


double DfCalcGridX::buildK(const TlMatrix& gridMatrix,
                           DfFunctional_GGA* pFunctional,
                           TlSymmetricMatrix* pF)
{
    //const index_type numOfAOs = this->m_nNumOfAOs;
    const double densityCutOffValue = this->m_densityCutOffValueA;

    double energy = 0.0;
    const std::size_t numOfAllGrids = gridMatrix.getNumOfRows();
    for (std::size_t grid = 0; grid < numOfAllGrids; ++grid) {
        const double x = gridMatrix.get(grid, 0);
        const double y = gridMatrix.get(grid, 1);
        const double z = gridMatrix.get(grid, 2);
        const double w = gridMatrix.get(grid, 3);
        const TlPosition position(x, y, z);
        const double rhoA = gridMatrix.get(grid, 4);
        const double gradRhoX_A = gridMatrix.get(grid, 5);
        const double gradRhoY_A = gridMatrix.get(grid, 6);
        const double gradRhoZ_A = gridMatrix.get(grid, 7);

        // calc phi table
        std::vector<WFGrid> aPhi;
        std::vector<WFGrid> aGradPhiX;
        std::vector<WFGrid> aGradPhiY;
        std::vector<WFGrid> aGradPhiZ;
        this->getPhiTable(position, aPhi, aGradPhiX, aGradPhiY, aGradPhiZ);
        std::sort(aPhi.begin(), aPhi.end(), WFGrid_sort_functional());
        std::sort(aGradPhiX.begin(), aGradPhiX.end(), WFGrid_sort_functional());
        std::sort(aGradPhiY.begin(), aGradPhiY.end(), WFGrid_sort_functional());
        std::sort(aGradPhiZ.begin(), aGradPhiZ.end(), WFGrid_sort_functional());

        if (rhoA > densityCutOffValue) {
            this->buildFock(rhoA, gradRhoX_A, gradRhoY_A, gradRhoZ_A,
                            aPhi, aGradPhiX, aGradPhiY, aGradPhiZ,
                            pFunctional, w, pF); // RKS code
            const double gammaAA =  gradRhoX_A*gradRhoX_A + gradRhoY_A*gradRhoY_A + gradRhoZ_A*gradRhoZ_A;
            energy += w * pFunctional->getFunctional(rhoA, gammaAA); // RKS code
        }
    }

    return energy;
}


double DfCalcGridX::buildK_2(const TlMatrix& gridMatrix,
                             const std::vector<index_type>& rowIndexes,
                             const std::vector<index_type>& colIndexes,
                             DfFunctional_GGA* pFunctional,
                             TlMatrix* pF)
{
    const double densityCutOffValue = this->m_densityCutOffValueA;

    TlMatrix tmpF(*pF);
    double energy = 0.0;
    const std::size_t numOfAllGrids = gridMatrix.getNumOfRows();
#pragma omp parallel for schedule(runtime) reduction(+:energy)
    for (std::size_t grid = 0; grid < numOfAllGrids; ++grid) {
        const double x = gridMatrix.get(grid, 0);
        const double y = gridMatrix.get(grid, 1);
        const double z = gridMatrix.get(grid, 2);
        const double w = gridMatrix.get(grid, 3);
        const TlPosition r(x, y, z);
        const double rhoA = gridMatrix.get(grid, 4);

        if (rhoA > densityCutOffValue) {
            const double gradRhoX_A = gridMatrix.get(grid, 5);
            const double gradRhoY_A = gridMatrix.get(grid, 6);
            const double gradRhoZ_A = gridMatrix.get(grid, 7);

            TlMatrix wf_r, dfdx_r, dfdy_r, dfdz_r;
            this->getWaveFunctionValues(rowIndexes, r,
                                        &wf_r, &dfdx_r, &dfdy_r, &dfdz_r);
            TlMatrix wf_c, dfdx_c, dfdy_c, dfdz_c;
            this->getWaveFunctionValues(colIndexes, r,
                                        &wf_c, &dfdx_c, &dfdy_c, &dfdz_c);
            wf_c.transpose();
            dfdx_c.transpose();
            dfdy_c.transpose();
            dfdz_c.transpose();
            
            const double gammaAA =  gradRhoX_A*gradRhoX_A + gradRhoY_A*gradRhoY_A + gradRhoZ_A*gradRhoZ_A;
            // build K
            double rF_rR_A;
            double rF_rG_AA, rF_rG_AB;
            pFunctional->getDerivativeFunctional(rhoA, gammaAA,
                                                 &rF_rR_A,
                                                 &rF_rG_AA,
                                                 &rF_rG_AB);
            
            const double coef1_A = w * rF_rR_A;
            const double rF_rG_AA2 = 2.0 * rF_rG_AA;
            const double coef2_AX = w * (rF_rG_AA2 * gradRhoX_A + rF_rG_AB * gradRhoX_A);
            const double coef2_AY = w * (rF_rG_AA2 * gradRhoY_A + rF_rG_AB * gradRhoY_A);
            const double coef2_AZ = w * (rF_rG_AA2 * gradRhoZ_A + rF_rG_AB * gradRhoZ_A);

            {
                TlMatrix wf_rc = wf_r * wf_c;
                wf_rc *= coef1_A;
                *pF += wf_rc;
            }
            {
                TlMatrix wf_rc = wf_r * dfdx_c;
                TlMatrix wf_cr = dfdx_r * wf_c;
                wf_rc *= coef2_AX;
                wf_cr *= coef2_AX;
                *pF += (wf_rc + wf_cr);
            }
            {
                TlMatrix wf_rc = wf_r * dfdy_c;
                TlMatrix wf_cr = dfdy_r * wf_c;
                wf_rc *= coef2_AY;
                wf_cr *= coef2_AY;
                *pF += (wf_rc + wf_cr);
            }
            {
                TlMatrix wf_rc = wf_r * dfdz_c;
                TlMatrix wf_cr = dfdz_r * wf_c;
                wf_rc *= coef2_AZ;
                wf_cr *= coef2_AZ;
                *pF += (wf_rc + wf_cr);
            }
            
            // energy
            energy += w * pFunctional->getFunctional(rhoA, gammaAA); // RKS code
        }
    }

    return energy;
}


