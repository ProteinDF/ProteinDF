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

#include <cmath>
#include <iostream>
#include <algorithm>
#include <limits>
#include <cassert>

#include "DfCalcGridX.h"
#include "CnError.h"

#include "DfGenerateGrid.h"

#include "Fl_Geometry.h"
#include "TlFile.h"
#include "TlUtils.h"
#include "TlPrdctbl.h"

////////////////////////////////////////////////////////////////////////
// DfCalcGridX
//

const double DfCalcGridX::TOOBIG = 30.0;
const double DfCalcGridX::EPS = std::numeric_limits<double>::epsilon();
const double DfCalcGridX::INV_SQRT3 = 1.0 / std::sqrt(3.0);
const double DfCalcGridX::INV_SQRT12 = 1.0 / std::sqrt(12.0);

DfCalcGridX::DfCalcGridX(TlSerializeData* pPdfParam)
    : DfObject(pPdfParam), m_tlOrbInfo((*pPdfParam)["coordinates"], (*pPdfParam)["basis_set"])
{
    const TlSerializeData& pdfParam = *pPdfParam;

    DfXCFunctional dfXcFunc(pPdfParam);
    this->functionalType_ = dfXcFunc.getFunctionalType();
    
    this->inputtedDensityCutoffValue_ = 1.0E-16;
    if (!(pdfParam["xc-density-threshold"].getStr().empty())) {
        this->inputtedDensityCutoffValue_ = pdfParam["xc-density-threshold"].getDouble();
    }
    this->m_densityCutOffValueA = this->inputtedDensityCutoffValue_;
    this->m_densityCutOffValueB = this->inputtedDensityCutoffValue_;

    this->m_inputedCutoffThreshold = pdfParam["cut_value"].getDouble();

    //this->physicalValues_.clear();

    // for debug
    this->isDebugOutPhiTable_ = (TlUtils::toUpper(pdfParam["debug_out_phi_table"].getStr()) == "YES") ? true : false;
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
    this->log_.info(TlUtils::format(" density cutoff value = %e",
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
    this->log_.info(TlUtils::format("density cutoff value(alpha) = %e",
                                    this->m_densityCutOffValueA));
    this->log_.info(TlUtils::format("density cutoff value(beta ) = %e",
                                    this->m_densityCutOffValueB));
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
    const double x = pos.x();
    const double y = pos.y();
    const double z = pos.z();
    
    switch (nType) {
    case 0:
        {
            *pPrefactorX = - alpha2 * x;
            *pPrefactorY = - alpha2 * y;
            *pPrefactorZ = - alpha2 * z;
        }
        break;
    case 1:
        {
            *pPrefactorX = - alpha2 * x * x + 1.0;
            *pPrefactorY = - alpha2 * x * y;
            *pPrefactorZ = - alpha2 * x * z;
        }
        break;
    case 2:
        {
            *pPrefactorX = - alpha2 * y * x;
            *pPrefactorY = - alpha2 * y * y + 1.0;
            *pPrefactorZ = - alpha2 * y * z;
        }
        break;
    case 3:
        {
            *pPrefactorX = - alpha2 * z * x;
            *pPrefactorY = - alpha2 * z * y;
            *pPrefactorZ = - alpha2 * z * z + 1.0;
        }
        break;
    case 4:
        {
            const double xy = x * y;
            *pPrefactorX = - alpha2 * xy * x + y; 
            *pPrefactorY = - alpha2 * xy * y + x; 
            *pPrefactorZ = - alpha2 * xy * z; 
        }
        break;
    case 5:
        {
            const double xz = x * z;
            *pPrefactorX = - alpha2 * xz * x + z; 
            *pPrefactorY = - alpha2 * xz * y;
            *pPrefactorZ = - alpha2 * xz * z + x;
        }
        break;
    case 6:
        {
            const double yz = y * z;
            *pPrefactorX = - alpha2 * yz * x; 
            *pPrefactorY = - alpha2 * yz * y + z;
            *pPrefactorZ = - alpha2 * yz * z + y;
        }
        break;
    case 7:
        {
            const double xx = x * x;
            const double xx_X = - alpha2 * xx * x + 2.0 * x;
            const double xx_Y = - alpha2 * xx * y;
            const double xx_Z = - alpha2 * xx * z;

            const double yy = y * y;
            const double yy_X = - alpha2 * yy * x;
            const double yy_Y = - alpha2 * yy * y + 2.0 * y;
            const double yy_Z = - alpha2 * yy * z;
            
            *pPrefactorX = 0.5 * (xx_X - yy_X);
            *pPrefactorY = 0.5 * (xx_Y - yy_Y);
            *pPrefactorZ = 0.5 * (xx_Z - yy_Z);
        }
        break;
    case 8:
        {
            const double xx = x * x;
            const double xx_X = - alpha2 * xx * x + 2.0 * x;
            const double xx_Y = - alpha2 * xx * y;
            const double xx_Z = - alpha2 * xx * z;

            const double yy = y * y;
            const double yy_X = - alpha2 * yy * x;
            const double yy_Y = - alpha2 * yy * y + 2.0 * y;
            const double yy_Z = - alpha2 * yy * z;

            const double zz = z * z;
            const double zz_X = - alpha2 * zz * x;
            const double zz_Y = - alpha2 * zz * y;
            const double zz_Z = - alpha2 * zz * z + 2.0 * z;

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


void DfCalcGridX::getPrefactorForSecondDerivative(const int nType, const double a, const TlPosition& pos,
                                                  double* pXX, double* pXY, double* pXZ,
                                                  double* pYY, double* pYZ, double* pZZ)
{
    assert(pXX != NULL);
    assert(pXY != NULL);
    assert(pXZ != NULL);
    assert(pYY != NULL);
    assert(pYZ != NULL);
    assert(pZZ != NULL);

    const double aa = a * a;
    const double x = pos.x();
    const double y = pos.y();
    const double z = pos.z();
    
    switch (nType) {
    case 0: // s
        {
            *pXX = 4.0 * aa * x * x - 2.0 * a;
            *pXY = 4.0 * aa * x * y;
            *pXZ = 4.0 * aa * x * z;
            *pYY = 4.0 * aa * y * y - 2.0 * a;
            *pYZ = 4.0 * aa * y * z;
            *pZZ = 4.0 * aa * z * z - 2.0 * a;
        }
        break;
    case 1: // px 
        {
            *pXX = 4.0 * aa * x * x * x - 6.0 * a * x;
            *pXY = 4.0 * aa * x * x * y - 2.0 * a * y;
            *pXZ = 4.0 * aa * x * x * z - 2.0 * a * z;
            *pYY = 4.0 * aa * x * y * y - 2.0 * a * x;
            *pYZ = 4.0 * aa * x * y * z;
            *pZZ = 4.0 * aa * x * z * z - 2.0 * a * x;
        }
        break;
    case 2: // py
        {
            *pXX = 4.0 * aa * x * x * y - 2.0 * a * y;
            *pXY = 4.0 * aa * x * y * y - 2.0 * a * x;
            *pXZ = 4.0 * aa * x * y * z;
            *pYY = 4.0 * aa * y * y * y - 6.0 * a * y;
            *pYZ = 4.0 * aa * y * y * z - 2.0 * a * z;
            *pZZ = 4.0 * aa * y * z * z - 2.0 * a * y;
        }
        break;
    case 3: // pz
        {
            *pXX = 4.0 * aa * x * x * z - 2.0 * a * z;
            *pXY = 4.0 * aa * x * y * z;
            *pXZ = 4.0 * aa * x * z * z - 2.0 * a * x;
            *pYY = 4.0 * aa * y * y * z - 2.0 * a * z;
            *pYZ = 4.0 * aa * y * z * z - 2.0 * a * y;
            *pZZ = 4.0 * aa * z * z * z - 6.0 * a * z;
        }
        break;
    case 4:
        { // dxy
            *pXX = 4.0 * aa * x * x * x * y - 6.0 * a * x * y;
            *pXY = 4.0 * aa * x * x * y * y - 2.0 * a * x * x - 2.0 * a * y * y + 1.0;
            *pXZ = 4.0 * aa * x * x * y * z - 2.0 * a * y * z;
            *pYY = 4.0 * aa * x * y * y * y - 6.0 * a * x * y;
            *pYZ = 4.0 * aa * x * y * y * z - 2.0 * a * x * z;
            *pZZ = 4.0 * aa * x * y * z * z - 2.0 * a * x * y;
        }
        break;
    case 5: // dxz
        {
            *pXX = 4.0 * aa * x * x * x * z - 6.0 * a * x * z;
            *pXY = 4.0 * aa * x * x * y * z - 2.0 * a * y * z;
            *pXZ = 4.0 * aa * x * x * z * z - 2.0 * a * x * x - 2.0 * a * z * z + 1.0;
            *pYY = 4.0 * aa * x * y * y * z - 2.0 * a * x * z;
            *pYZ = 4.0 * aa * x * y * z * z - 2.0 * a * x * y;
            *pZZ = 4.0 * aa * x * z * z * z - 6.0 * a * x * z;
        }
        break;
    case 6: // yz
        {
            *pXX = 4.0 * aa * x * x * y * z - 2.0 * a * y * z;
            *pXY = 4.0 * aa * x * y * y * z - 2.0 * a * x * z;
            *pXZ = 4.0 * aa * x * y * z * z - 2.0 * a * x * y;
            *pYY = 4.0 * aa * y * y * y * z - 6.0 * a * y * z;
            *pYZ = 4.0 * aa * y * y * z * z - 2.0 * a * y * y - 2.0 * a * z * z + 1.0;
            *pZZ = 4.0 * aa * y * z * z * z - 6.0 * a * y * z;
        }
        break;
    case 7: // xx-yy
        {
            // xx
            const double xx_xx = 4.0 * aa * x * x * x * x - 10.0 * a * x * x + 2.0;
            const double xx_xy = 4.0 * aa * x * x * x * y -  4.0 * a * x * y;
            const double xx_xz = 4.0 * aa * x * x * x * z -  4.0 * a * x * z;
            const double xx_yy = 4.0 * aa * x * x * y * y -  2.0 * a * x * x;
            const double xx_yz = 4.0 * aa * x * x * y * z;
            const double xx_zz = 4.0 * aa * x * x * z * z -  2.0 * a * x * x;

            // yy
            const double yy_xx = 4.0 * aa * x * x * y * y -  2.0 * a * y * y;
            const double yy_xy = 4.0 * aa * x * y * y * y -  4.0 * a * x * y;
            const double yy_xz = 4.0 * aa * x * y * y * z;
            const double yy_yy = 4.0 * aa * y * y * y * y - 10.0 * a * y * y + 2.0;
            const double yy_yz = 4.0 * aa * y * y * y * z -  4.0 * a * y * z;
            const double yy_zz = 4.0 * aa * y * y * z * z -  2.0 * a * y * y;

            *pXX = 0.5 * (xx_xx - yy_xx);
            *pXY = 0.5 * (xx_xy - yy_xy);
            *pXZ = 0.5 * (xx_xz - yy_xz);
            *pYY = 0.5 * (xx_yy - yy_yy);
            *pYZ = 0.5 * (xx_yz - yy_yz);
            *pZZ = 0.5 * (xx_zz - yy_zz);
        }
        break;
    case 8: // r2
        {
            // xx
            const double xx_xx = 4.0 * aa * x * x * x * x - 10.0 * a * x * x + 2.0;
            const double xx_xy = 4.0 * aa * x * x * x * y -  4.0 * a * x * y;
            const double xx_xz = 4.0 * aa * x * x * x * z -  4.0 * a * x * z;
            const double xx_yy = 4.0 * aa * x * x * y * y -  2.0 * a * x * x;
            const double xx_yz = 4.0 * aa * x * x * y * z;
            const double xx_zz = 4.0 * aa * x * x * z * z -  2.0 * a * x * x;

            // yy
            const double yy_xx = 4.0 * aa * x * x * y * y -  2.0 * a * y * y;
            const double yy_xy = 4.0 * aa * x * y * y * y -  4.0 * a * x * y;
            const double yy_xz = 4.0 * aa * x * y * y * z;
            const double yy_yy = 4.0 * aa * y * y * y * y - 10.0 * a * y * y + 2.0;
            const double yy_yz = 4.0 * aa * y * y * y * z -  4.0 * a * y * z;
            const double yy_zz = 4.0 * aa * y * y * z * z -  2.0 * a * y * y;

            // zz
            const double zz_xx = 4.0 * aa * x * x * z * z -  2.0 * a * z * z;
            const double zz_xy = 4.0 * aa * x * y * z * z;
            const double zz_xz = 4.0 * aa * x * z * z * z -  4.0 * a * x * z;
            const double zz_yy = 4.0 * aa * y * y * z * z -  2.0 * a * z * z;
            const double zz_yz = 4.0 * aa * y * z * z * z -  4.0 * a * y * z;
            const double zz_zz = 4.0 * aa * z * z * z * z - 10.0 * a * z * z + 2.0;

            *pXX = INV_SQRT3 * (zz_xx - 0.5 * (xx_xx + yy_xx));
            *pXY = INV_SQRT3 * (zz_xy - 0.5 * (xx_xy + yy_xy));
            *pXZ = INV_SQRT3 * (zz_xz - 0.5 * (xx_xz + yy_xz));
            *pYY = INV_SQRT3 * (zz_yy - 0.5 * (xx_yy + yy_yy));
            *pYZ = INV_SQRT3 * (zz_yz - 0.5 * (xx_yz + yy_yz));
            *pZZ = INV_SQRT3 * (zz_zz - 0.5 * (xx_zz + yy_zz));
        }
        break;
    default:
        std::cout << "Basis Type is Wrong." << std::endl;
        break;
    }
}



////////////////////////////////////////////////////////////////////////

// φの値を求める
// for NSD
void DfCalcGridX::getPhiTable(const TlPosition& gridPosition, std::vector<WFGrid>& aPhi)
{
    this->getPhiTable(gridPosition, 0, this->m_tlOrbInfo.getNumOfOrbitals(), aPhi);

    std::sort(aPhi.begin(), aPhi.end(), WFGrid_sort_functional());
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

    std::sort(aPhi.begin(), aPhi.end(), WFGrid_sort_functional());
    std::sort(aGradPhiX.begin(), aGradPhiX.end(), WFGrid_sort_functional());
    std::sort(aGradPhiY.begin(), aGradPhiY.end(), WFGrid_sort_functional());
    std::sort(aGradPhiZ.begin(), aGradPhiZ.end(), WFGrid_sort_functional());
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


void DfCalcGridX::getDerivative2AOs(const TlPosition& gridPosition,
                                    std::vector<WFGrid>* p_d2AO_dxdx_values,
                                    std::vector<WFGrid>* p_d2AO_dxdy_values,
                                    std::vector<WFGrid>* p_d2AO_dxdz_values,
                                    std::vector<WFGrid>* p_d2AO_dydy_values,
                                    std::vector<WFGrid>* p_d2AO_dydz_values,
                                    std::vector<WFGrid>* p_d2AO_dzdz_values)
{
    this->getDerivative2AOs(gridPosition, 0, this->m_tlOrbInfo.getNumOfOrbitals(),
                            p_d2AO_dxdx_values, p_d2AO_dxdy_values, p_d2AO_dxdz_values,
                            p_d2AO_dydy_values, p_d2AO_dydz_values, p_d2AO_dzdz_values);
}

void DfCalcGridX::getDerivative2AOs(const TlPosition& gridPosition,
                                    const int start_AO_index,
                                    const int end_AO_index,
                                    std::vector<WFGrid>* p_d2AO_dxdx_values,
                                    std::vector<WFGrid>* p_d2AO_dxdy_values,
                                    std::vector<WFGrid>* p_d2AO_dxdz_values,
                                    std::vector<WFGrid>* p_d2AO_dydy_values,
                                    std::vector<WFGrid>* p_d2AO_dydz_values,
                                    std::vector<WFGrid>* p_d2AO_dzdz_values)
{
    assert(p_d2AO_dxdx_values != NULL);
    assert(p_d2AO_dxdy_values != NULL);
    assert(p_d2AO_dxdz_values != NULL);
    assert(p_d2AO_dydy_values != NULL);
    assert(p_d2AO_dydz_values != NULL);
    assert(p_d2AO_dzdz_values != NULL);

    // initialize
    const double cutOffValue = std::sqrt(std::min(this->m_densityCutOffValueA,
                                                  this->m_densityCutOffValueB));

    const int AO_ListSize = end_AO_index - start_AO_index +1;
    p_d2AO_dxdx_values->clear();
    p_d2AO_dxdx_values->reserve(AO_ListSize);
    p_d2AO_dxdy_values->clear();
    p_d2AO_dxdy_values->reserve(AO_ListSize);
    p_d2AO_dxdz_values->clear();
    p_d2AO_dxdz_values->reserve(AO_ListSize);
    p_d2AO_dydy_values->clear();
    p_d2AO_dydy_values->reserve(AO_ListSize);
    p_d2AO_dydz_values->clear();
    p_d2AO_dydz_values->reserve(AO_ListSize);

    // orbital loop
    for (int AO_index = start_AO_index; AO_index < end_AO_index; ++AO_index) {
        double d2AO_dxdx = 0.0;
        double d2AO_dxdy = 0.0;
        double d2AO_dxdz = 0.0;
        double d2AO_dydy = 0.0;
        double d2AO_dydz = 0.0;
        double d2AO_dzdz = 0.0;

        const TlPosition pos = gridPosition - this->m_tlOrbInfo.getPosition(AO_index);
        const double distance2 = pos.squareDistanceFrom();
        const int basisType = this->m_tlOrbInfo.getBasisType(AO_index);

        const int contraction = this->m_tlOrbInfo.getCgtoContraction(AO_index);
        for (int PGTO = 0; PGTO < contraction; ++PGTO) {
            double pXX = 0.0;
            double pXY = 0.0;
            double pXZ = 0.0;
            double pYY = 0.0;
            double pYZ = 0.0;
            double pZZ = 0.0;
            const double alpha = this->m_tlOrbInfo.getExponent(AO_index, PGTO);
            const double shoulder = alpha * distance2;

            if (shoulder <= TOOBIG) {
                const double g = this->m_tlOrbInfo.getCoefficient(AO_index, PGTO) * std::exp(- shoulder);
                this->getPrefactorForSecondDerivative(basisType, alpha, pos,
                                                      &pXX, &pXY, &pXZ, &pYY, &pYZ, &pZZ);
                d2AO_dxdx += pXX * g;
                d2AO_dxdy += pXY * g;
                d2AO_dxdz += pXZ * g;
                d2AO_dydy += pYY * g;
                d2AO_dydz += pYZ * g;
                d2AO_dzdz += pZZ * g;
            }
        }
        
        if (std::fabs(d2AO_dxdx) > cutOffValue) {
            const WFGrid wfGrid(AO_index, d2AO_dxdx);
            p_d2AO_dxdx_values->push_back(wfGrid);
        }
        if (std::fabs(d2AO_dxdy) > cutOffValue) {
            const WFGrid wfGrid(AO_index, d2AO_dxdy);
            p_d2AO_dxdy_values->push_back(wfGrid);
        }
        if (std::fabs(d2AO_dxdz) > cutOffValue) {
            const WFGrid wfGrid(AO_index, d2AO_dxdz);
            p_d2AO_dxdz_values->push_back(wfGrid);
        }
        if (std::fabs(d2AO_dydy) > cutOffValue) {
            const WFGrid wfGrid(AO_index, d2AO_dydy);
            p_d2AO_dydy_values->push_back(wfGrid);
        }
        if (std::fabs(d2AO_dydz) > cutOffValue) {
            const WFGrid wfGrid(AO_index, d2AO_dydz);
            p_d2AO_dydz_values->push_back(wfGrid);
        }
        if (std::fabs(d2AO_dzdz) > cutOffValue) {
            const WFGrid wfGrid(AO_index, d2AO_dzdz);
            p_d2AO_dzdz_values->push_back(wfGrid);
        }
    }
}


void DfCalcGridX::getRhoAtGridPoint(const TlMatrixObject& P,
                                    const std::vector<WFGrid>& aPhi,
                                    double* pRhoA)
{
    const double densityCutOffValue = this->m_densityCutOffValueA;
    double dRho = 0.0;

    // 密度行列のループ(index = p)の最大を求める
    std::vector<WFGrid>::const_iterator pEnd = std::upper_bound(aPhi.begin(), aPhi.end(),
                                                                WFGrid(0, densityCutOffValue),
                                                                WFGrid_sort_functional());
    const std::size_t max_p = std::distance(aPhi.begin(), pEnd);
    for (std::size_t p = 0; p < max_p; ++p) {
        const std::size_t nOrb_p = aPhi[p].index;
        const double phi_p = aPhi[p].value;
        //const double cutValue = std::fabs(densityCutOffValue / phi_p);

        // 高速化
        const TlVector P_row = P.getRowVector(nOrb_p);

        for (std::size_t q = 0; q < max_p; ++q) {
            const std::size_t nOrb_q = aPhi[q].index;
            const double phi_q = aPhi[q].value;
            dRho += P_row[nOrb_q] * phi_p * phi_q;
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
    const double densityCutOffValue = std::min(this->m_densityCutOffValueA,
                                               this->m_densityCutOffValueB);
    double dRho = 0.0;
    double dGradRhoX = 0.0;
    double dGradRhoY = 0.0;
    double dGradRhoZ = 0.0;

    std::vector<WFGrid>::const_iterator pEnd = std::upper_bound(aPhi.begin(), aPhi.end(),
                                                                WFGrid(0, densityCutOffValue),
                                                                WFGrid_sort_functional());
    const std::size_t max_p = std::distance(aPhi.begin(), pEnd);
    for (std::size_t p = 0; p < max_p; ++p) {
        const std::size_t nOrb_p = aPhi[p].index;
        const double phi_p = aPhi[p].value;
        const double cutValue = std::fabs(densityCutOffValue / phi_p);

        // 高速化
        const TlVector P_row = P.getRowVector(nOrb_p);

        for (std::size_t q = 0; q < max_p; ++q) {
            const std::size_t nOrb_q = aPhi[q].index;
            const double phi_q = aPhi[q].value;
            dRho += P_row[nOrb_q] * phi_p * phi_q;
        }

        std::vector<WFGrid>::const_iterator qxEnd = std::upper_bound(aGradPhiX.begin(),
                                                                     aGradPhiX.end(),
                                                                     WFGrid(0, cutValue),
                                                                     WFGrid_sort_functional());
        std::size_t max_qx = std::distance(aGradPhiX.begin(), qxEnd);
        for (std::size_t qx = 0; qx < max_qx; ++qx) {
            const std::size_t nOrb_q = aGradPhiX[qx].index;
            const double phi_q = aGradPhiX[qx].value;
            dGradRhoX += P_row[nOrb_q] * phi_p * phi_q;
        }

        std::vector<WFGrid>::const_iterator qyEnd = std::upper_bound(aGradPhiY.begin(),
                                                                     aGradPhiY.end(),
                                                                     WFGrid(0, cutValue),
                                                                     WFGrid_sort_functional());
        std::size_t max_qy = std::distance(aGradPhiY.begin(), qyEnd);
        for (std::size_t qy = 0; qy < max_qy; ++qy) {
            const std::size_t nOrb_q = aGradPhiY[qy].index;
            const double phi_q = aGradPhiY[qy].value;
            dGradRhoY += P_row[nOrb_q] * phi_p * phi_q;
        }

        std::vector<WFGrid>::const_iterator qzEnd = std::upper_bound(aGradPhiZ.begin(),
                                                                     aGradPhiZ.end(),
                                                                     WFGrid(0, cutValue),
                                                                     WFGrid_sort_functional());
        std::size_t max_qz = std::distance(aGradPhiZ.begin(), qzEnd);
        for (std::size_t qz = 0; qz < max_qz; ++qz) {
            const std::size_t nOrb_q = aGradPhiZ[qz].index;
            const double phi_q = aGradPhiZ[qz].value;
            dGradRhoZ += P_row[nOrb_q] * phi_p * phi_q;
        }
    }

    *pRhoA = dRho;
    *pGradRhoAX = 2.0 * dGradRhoX;
    *pGradRhoAY = 2.0 * dGradRhoY;
    *pGradRhoAZ = 2.0 * dGradRhoZ;
}


////////////////////////////////////////////////////////////////////////////////
// build Fock matrix

// for LDA and RKS
// void DfCalcGridX::buildFock(const double dRhoA, const std::vector<WFGrid>& aPhi,
//                             DfFunctional_LDA* pFunctional, const double dWeight,
//                             TlMatrixObject* pF)
// {
//     double dRoundF_roundRhoA;
//     pFunctional->getDerivativeFunctional(dRhoA, &dRoundF_roundRhoA);

//     const double dCoeff1_A = dWeight * dRoundF_roundRhoA;
//     const double dCutOffValue = this->m_inputedCutoffThreshold;

//     // phi v.s. phi
//     if ((aPhi.size() > 0) &&
//         (std::fabs(dCoeff1_A * aPhi[0].value * aPhi[0].value) > dCutOffValue)) {
//         std::vector<WFGrid>::const_iterator pEnd =
//             std::upper_bound(aPhi.begin(), aPhi.end(),
//                              WFGrid(0, std::sqrt(std::fabs(dCutOffValue / dCoeff1_A))),
//                              WFGrid_sort_functional());
//         this->buildFock(aPhi.begin(), pEnd, dCoeff1_A, dCutOffValue, pF);
//     }
// }

// LDA and UKS
// void DfCalcGridX::buildFock(const double dRhoA, const double dRhoB,
//                             const std::vector<WFGrid>& aPhi,
//                             DfFunctional_LDA* pFunctional, const double dWeight,
//                             TlMatrixObject* pFA, TlMatrixObject* pFB)
// {
//     double dRoundF_roundRhoA;
//     double dRoundF_roundRhoB;
//     pFunctional->getDerivativeFunctional(dRhoA, dRhoB, &dRoundF_roundRhoA, &dRoundF_roundRhoB);

//     const double dCoeff1_A = dWeight * dRoundF_roundRhoA;
//     const double dCoeff1_B = dWeight * dRoundF_roundRhoB;
//     const double dCutOffValue = this->m_inputedCutoffThreshold;

//     if (aPhi.size() > 0) {
//         const double aPhi_0 = aPhi[0].value;
//         if (std::fabs(dCoeff1_A * aPhi_0 * aPhi_0) > dCutOffValue) {
//             std::vector<WFGrid>::const_iterator pEnd =
//                 std::upper_bound(aPhi.begin(), aPhi.end(),
//                                  WFGrid(0, std::sqrt(std::fabs(dCutOffValue/dCoeff1_A))),
//                                  WFGrid_sort_functional());
//             this->buildFock(aPhi.begin(), pEnd, dCoeff1_A, dCutOffValue, pFA);
//         }
//         if (std::fabs(dCoeff1_B * aPhi_0 * aPhi_0) > dCutOffValue) {
//             std::vector<WFGrid>::const_iterator pEnd =
//                 std::upper_bound(aPhi.begin(), aPhi.end(),
//                                  WFGrid(0, std::sqrt(std::fabs(dCutOffValue/dCoeff1_B))),
//                                  WFGrid_sort_functional());
//             this->buildFock(aPhi.begin(), pEnd, dCoeff1_B, dCutOffValue, pFB);
//         }
//     }
// }


// for GGA and RKS
// void DfCalcGridX::buildFock(const double dRhoA,
//                             const double dGradRhoAX, const double dGradRhoAY, const double dGradRhoAZ,
//                             const std::vector<WFGrid>& aPhi,
//                             const std::vector<WFGrid>& aGradPhiX,
//                             const std::vector<WFGrid>& aGradPhiY,
//                             const std::vector<WFGrid>& aGradPhiZ,
//                             DfFunctional_GGA* pFunctional, const double dWeight,
//                             TlMatrixObject *pF)
// {

//     const double dGammaAA = dGradRhoAX * dGradRhoAX + dGradRhoAY * dGradRhoAY + dGradRhoAZ * dGradRhoAZ;

//     double dRoundF_roundRhoA;
//     double dRoundF_roundGammaAA, dRoundF_roundGammaAB;
//     pFunctional->getDerivativeFunctional(dRhoA, dGammaAA,
//                                          &dRoundF_roundRhoA, &dRoundF_roundGammaAA, &dRoundF_roundGammaAB);

//     const double dCoeff1_A = dWeight * dRoundF_roundRhoA;
//     const double dRoundF_roundGammaAA2 = 2.0 * dRoundF_roundGammaAA;

//     // remark! dGradRhoBX = dGradRhoAX
//     const double dCoeff2_AX =
//         dWeight * (dRoundF_roundGammaAA2 + dRoundF_roundGammaAB) * dGradRhoAX;
//     const double dCoeff2_AY =
//         dWeight * (dRoundF_roundGammaAA2 + dRoundF_roundGammaAB) * dGradRhoAY;
//     const double dCoeff2_AZ = 
//         dWeight * (dRoundF_roundGammaAA2 + dRoundF_roundGammaAB) * dGradRhoAZ;
    
//     const double dCutOffValue = this->m_inputedCutoffThreshold;

//     if (aPhi.size() > 0) {
//         const double aPhi_0 = aPhi[0].value;
//         // phi v.s. phi
//         if (std::fabs(dCoeff1_A * aPhi_0 * aPhi_0) > dCutOffValue) {
//             std::vector<WFGrid>::const_iterator pEnd =
//                 std::upper_bound(aPhi.begin(), aPhi.end(),
//                                  WFGrid(0, std::sqrt(std::fabs(dCutOffValue / dCoeff1_A))),
//                                  WFGrid_sort_functional());
//             this->buildFock(aPhi.begin(), pEnd, dCoeff1_A, dCutOffValue, pF);
//         }

//         // phi v.s. grad-phi(x)
//         if ((aGradPhiX.size() > 0) &&
//             (std::fabs(dCoeff2_AX * aPhi_0 * aGradPhiX[0].value) > dCutOffValue)) {
//             std::vector<WFGrid>::const_iterator pEnd =
//                 std::upper_bound(aPhi.begin(), aPhi.end(),
//                                  WFGrid(0, dCutOffValue / (dCoeff2_AX * aGradPhiX[0].value)),
//                                  WFGrid_sort_functional());
            
//             this->buildFock(aPhi.begin(), pEnd, aGradPhiX.begin(), aGradPhiX.end(),
//                             dCoeff2_AX, dCutOffValue, pF);
//             this->buildFock(aGradPhiX.begin(), aGradPhiX.end(), aPhi.begin(), pEnd,
//                             dCoeff2_AX, dCutOffValue, pF);
//         }
        
//         // phi v.s. grad-phi(y)
//         if ((aGradPhiY.size() > 0) &&
//             (std::fabs(dCoeff2_AY * aPhi_0 * aGradPhiY[0].value) > dCutOffValue)) {
//             std::vector<WFGrid>::const_iterator pEnd =
//                 std::upper_bound(aPhi.begin(), aPhi.end(),
//                                  WFGrid(0, dCutOffValue / (dCoeff2_AY * aGradPhiY[0].value)),
//                                  WFGrid_sort_functional());
            
//             this->buildFock(aPhi.begin(), pEnd, aGradPhiY.begin(), aGradPhiY.end(),
//                             dCoeff2_AY, dCutOffValue, pF);
//             this->buildFock(aGradPhiY.begin(), aGradPhiY.end(), aPhi.begin(), pEnd,
//                             dCoeff2_AY, dCutOffValue, pF);
//         }

//         // phi v.s. grad-phi(z)
//         if ((aGradPhiZ.size() > 0) &&
//             (std::fabs(dCoeff2_AZ * aPhi_0 * aGradPhiZ[0].value) > dCutOffValue)) {
//             std::vector<WFGrid>::const_iterator pEnd =
//                 std::upper_bound(aPhi.begin(), aPhi.end(),
//                                  WFGrid(0, dCutOffValue / (dCoeff2_AZ * aGradPhiZ[0].value)),
//                                  WFGrid_sort_functional());
            
//             this->buildFock(aPhi.begin(), pEnd, aGradPhiZ.begin(), aGradPhiZ.end(),
//                             dCoeff2_AZ, dCutOffValue, pF);
//             this->buildFock(aGradPhiZ.begin(), aGradPhiZ.end(), aPhi.begin(), pEnd,
//                             dCoeff2_AZ, dCutOffValue, pF);
//         }
//     }
// }


// for GGA and UKS
// void DfCalcGridX::buildFock(const double dRhoA, const double dRhoB,
//                             const double dGradRhoAX, const double dGradRhoAY, const double dGradRhoAZ,
//                             const double dGradRhoBX, const double dGradRhoBY, const double dGradRhoBZ,
//                             const std::vector<WFGrid>& aPhi,
//                             const std::vector<WFGrid>& aGradPhiX,
//                             const std::vector<WFGrid>& aGradPhiY,
//                             const std::vector<WFGrid>& aGradPhiZ,
//                             DfFunctional_GGA* pFunctional, const double dWeight,
//                             TlMatrixObject* pFA, TlMatrixObject* pFB)
// {
//     assert(dRhoA >= 0.0);
//     assert(dRhoB >= 0.0);
//     const double dGammaAA = dGradRhoAX*dGradRhoAX + dGradRhoAY*dGradRhoAY + dGradRhoAZ*dGradRhoAZ;
//     const double dGammaAB = dGradRhoAX*dGradRhoBX + dGradRhoAY*dGradRhoBY + dGradRhoAZ*dGradRhoBZ;
//     const double dGammaBB = dGradRhoBX*dGradRhoBX + dGradRhoBY*dGradRhoBY + dGradRhoBZ*dGradRhoBZ;

//     double dRoundF_roundRhoA, dRoundF_roundRhoB;
//     double dRoundF_roundGammaAA, dRoundF_roundGammaAB, dRoundF_roundGammaBB;
//     pFunctional->getDerivativeFunctional(dRhoA, dRhoB, dGammaAA, dGammaAB, dGammaBB,
//                                          &dRoundF_roundRhoA, &dRoundF_roundRhoB,
//                                          &dRoundF_roundGammaAA, &dRoundF_roundGammaAB, &dRoundF_roundGammaBB);

//     const double dCoeff1_A = dWeight * dRoundF_roundRhoA;
//     const double dCoeff1_B = dWeight * dRoundF_roundRhoB;
//     const double dRoundF_roundGammaAA2 = 2.0 * dRoundF_roundGammaAA;
//     const double dRoundF_roundGammaBB2 = 2.0 * dRoundF_roundGammaBB;
//     const double dCoeff2_AX = dWeight
//                               * (dRoundF_roundGammaAA2 * dGradRhoAX + dRoundF_roundGammaAB * dGradRhoAX);
//     const double dCoeff2_AY = dWeight
//                               * (dRoundF_roundGammaAA2 * dGradRhoAY + dRoundF_roundGammaAB * dGradRhoAY);
//     const double dCoeff2_AZ = dWeight
//                               * (dRoundF_roundGammaAA2 * dGradRhoAZ + dRoundF_roundGammaAB * dGradRhoAZ);
//     const double dCoeff2_BX = dWeight
//                               * (dRoundF_roundGammaBB2 * dGradRhoBX + dRoundF_roundGammaAB * dGradRhoBX);
//     const double dCoeff2_BY = dWeight
//                               * (dRoundF_roundGammaBB2 * dGradRhoBY + dRoundF_roundGammaAB * dGradRhoBY);
//     const double dCoeff2_BZ = dWeight
//                               * (dRoundF_roundGammaBB2 * dGradRhoBZ + dRoundF_roundGammaAB * dGradRhoBZ);

//     const double dCutOffValue = this->m_inputedCutoffThreshold;

//     if (aPhi.size() > 0) {
//         const double aPhi_0 = aPhi[0].value;
//         // phi v.s. phi
//         if (std::fabs(dCoeff1_A * aPhi_0 * aPhi_0) > dCutOffValue) {
//             std::vector<WFGrid>::const_iterator pEnd = std::upper_bound(aPhi.begin(), aPhi.end(),
//                                                                         WFGrid(0, std::sqrt(std::fabs(dCutOffValue / dCoeff1_A))),
//                                                                         WFGrid_sort_functional());
//             this->buildFock(aPhi.begin(), pEnd, dCoeff1_A, dCutOffValue, pFA);
//         }
//         if (std::fabs(dCoeff1_B * aPhi_0 * aPhi_0) > dCutOffValue) {
//             std::vector<WFGrid>::const_iterator pEnd = std::upper_bound(aPhi.begin(), aPhi.end(),
//                                                                         WFGrid(0, std::sqrt(std::fabs(dCutOffValue / dCoeff1_B))),
//                                                                         WFGrid_sort_functional());
//             this->buildFock(aPhi.begin(), pEnd, dCoeff1_B, dCutOffValue, pFB);
//         }

//         // phi v.s. grad-phi(x)
//         if ((aGradPhiX.size() > 0) &&
//                 (std::fabs(dCoeff2_AX * aPhi_0 * aGradPhiX[0].value) > dCutOffValue)) {
//             std::vector<WFGrid>::const_iterator pEnd = std::upper_bound(aPhi.begin(), aPhi.end(),
//                                                                         WFGrid(0, dCutOffValue / (dCoeff2_AX * aGradPhiX[0].value)),
//                                                                         WFGrid_sort_functional());
//             this->buildFock(aPhi.begin(), pEnd, aGradPhiX.begin(), aGradPhiX.end(),
//                             dCoeff2_AX, dCutOffValue, pFA);
//         }
//         if ((aGradPhiX.size() > 0) &&
//                 (std::fabs(dCoeff2_BX * aPhi_0 * aGradPhiX[0].value) > dCutOffValue)) {
//             std::vector<WFGrid>::const_iterator pEnd = std::upper_bound(aPhi.begin(), aPhi.end(),
//                                                                         WFGrid(0, dCutOffValue / (dCoeff2_BX * aGradPhiX[0].value)),
//                                                                         WFGrid_sort_functional());
//             this->buildFock(aPhi.begin(), pEnd, aGradPhiX.begin(), aGradPhiX.end(),
//                             dCoeff2_BX, dCutOffValue, pFB);
//         }

//         // phi v.s. grad-phi(y)
//         if ((aGradPhiY.size() > 0) &&
//                 (std::fabs(dCoeff2_AY * aPhi_0 * aGradPhiY[0].value) > dCutOffValue)) {
//             std::vector<WFGrid>::const_iterator pEnd = std::upper_bound(aPhi.begin(), aPhi.end(),
//                                                                         WFGrid(0, dCutOffValue / (dCoeff2_AY * aGradPhiY[0].value)),
//                                                                         WFGrid_sort_functional());
//             this->buildFock(aPhi.begin(), pEnd, aGradPhiY.begin(), aGradPhiY.end(),
//                             dCoeff2_AY, dCutOffValue, pFA);
//         }
//         if ((aGradPhiY.size() > 0) &&
//                 (std::fabs(dCoeff2_BY * aPhi_0 * aGradPhiY[0].value) > dCutOffValue)) {
//             std::vector<WFGrid>::const_iterator pEnd = std::upper_bound(aPhi.begin(), aPhi.end(),
//                                                                         WFGrid(0, dCutOffValue / (dCoeff2_BY * aGradPhiY[0].value)),
//                                                                         WFGrid_sort_functional());
//             this->buildFock(aPhi.begin(), pEnd, aGradPhiY.begin(), aGradPhiY.end(),
//                             dCoeff2_BY, dCutOffValue, pFB);
//         }

//         // phi v.s. grad-phi(z)
//         if ((aGradPhiZ.size() > 0) &&
//                 (std::fabs(dCoeff2_AZ * aPhi_0 * aGradPhiZ[0].value) > dCutOffValue)) {
//             std::vector<WFGrid>::const_iterator pEnd = std::upper_bound(aPhi.begin(), aPhi.end(),
//                                                                         WFGrid(0, dCutOffValue / (dCoeff2_AZ * aGradPhiZ[0].value)),
//                                                                         WFGrid_sort_functional());
//             this->buildFock(aPhi.begin(), pEnd, aGradPhiZ.begin(), aGradPhiZ.end(),
//                             dCoeff2_AZ, dCutOffValue, pFA);
//         }
//         if ((aGradPhiZ.size() > 0) &&
//                 (std::fabs(dCoeff2_BZ * aPhi_0 * aGradPhiZ[0].value) > dCutOffValue)) {
//             std::vector<WFGrid>::const_iterator pEnd = std::upper_bound(aPhi.begin(), aPhi.end(),
//                                                                         WFGrid(0, dCutOffValue / (dCoeff2_BZ * aGradPhiZ[0].value)),
//                                                                         WFGrid_sort_functional());
//             this->buildFock(aPhi.begin(), pEnd, aGradPhiZ.begin(), aGradPhiZ.end(),
//                             dCoeff2_BZ, dCutOffValue, pFB);
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

        // 対角要素以外
        std::vector<WFGrid>::const_iterator qEnd = std::upper_bound(p +1, pEnd,
                                                                    WFGrid(0, cutoffValue / tmp1A),
                                                                    WFGrid_sort_functional());
        for (std::vector<WFGrid>::const_iterator q = p +1; q != qEnd; ++q) {
            const int v_index = q->index;
            const double v_value = q->value;

            pF->add(u_index, v_index, tmp1A * v_value);
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

    // get rho at grid point
    //double dRhoA;
    this->getRhoAtGridPoint(P, aPhi, pRho);
}


// void DfCalcGridX::calcForceFromXC(const TlSymmetricMatrix& P,
//                                   DfFunctional_GGA* pFunctional,
//                                   TlVector* pFx, TlVector* pFy, TlVector* pFz)
// {
//     const int numOfAtoms = this->numOfRealAtoms_;
//     const int numOfAOs = this->m_nNumOfAOs;

//     // grid derivative
//     // to implement!
    
//     // functional derivative
//     TlMatrix GX(numOfAOs, numOfAtoms);
//     TlMatrix GY(numOfAOs, numOfAtoms);
//     TlMatrix GZ(numOfAOs, numOfAtoms);
//     this->makeGammaMatrix(P, pFunctional,
//                           &GX, &GY, &GZ);

//     for (int atomIndex = 0; atomIndex < numOfAtoms; ++atomIndex) {
//         double term1X = 0.0;
//         double term1Y = 0.0;
//         double term1Z = 0.0;
//         double term2X = 0.0;
//         double term2Y = 0.0;
//         double term2Z = 0.0;

//         for (int aoIndex = 0; aoIndex < numOfAOs; ++aoIndex) {
//             if (this->m_tlOrbInfo.getAtomIndex(aoIndex) != atomIndex) {
//                 term1X += GX.get(atomIndex, aoIndex);
//                 term1Y += GY.get(atomIndex, aoIndex);
//                 term1Z += GZ.get(atomIndex, aoIndex);
//             } else {
//                 for (int atomIndex2 = 0; atomIndex2 < numOfAtoms; ++atomIndex2) {
//                     if (atomIndex != atomIndex2) {
//                         term2X += GX.get(atomIndex2, aoIndex);
//                         term2Y += GY.get(atomIndex2, aoIndex);
//                         term2Z += GZ.get(atomIndex2, aoIndex);
//                     }
//                 }
//             }
//         }

//         pFx->set(atomIndex, (term1X - term2X));
//         pFy->set(atomIndex, (term1Y - term2Y));
//         pFz->set(atomIndex, (term1Z - term2Z));
//     }
// }


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

    // const double densityCutOffValue = this->m_densityCutOffValueA;
    const Fl_Geometry flGeom((*this->pPdfParam_)["coordinates"]);
    TlMatrix gridMat = this->getGridMatrix<TlMatrix>(this->m_nIteration);

    const int numOfAtoms = this->m_nNumOfAtoms;
    const int numOfAOs = this->m_nNumOfAOs;
    
    pGX->resize(numOfAOs, numOfAtoms);
    pGY->resize(numOfAOs, numOfAtoms);
    pGZ->resize(numOfAOs, numOfAtoms);
    
    for (int atomIndex = 0; atomIndex < numOfAtoms; ++atomIndex) {
        const int atomicNumber = TlPrdctbl::getAtomicNumber(flGeom.getAtomSymbol(atomIndex));
        if (atomicNumber == 0) {
            // pass for dummy charge(X)
            continue;
        }
        
        // get grid information
        const TlMatrix atomGridMat = this->selectGridMatrixByAtom(gridMat, atomIndex);
        const index_type numOfGrids = atomGridMat.getNumOfRows();
        // std::vector<double> rhoA(numOfGrids, 0.0);
        // std::vector<double> gradRhoAX(numOfGrids, 0.0);
        // std::vector<double> gradRhoAY(numOfGrids, 0.0);
        // std::vector<double> gradRhoAZ(numOfGrids, 0.0);
        // 電子密度を計算する
        // LDAではSCF中に微分を計算しないため
        // for (int gridIndex = 0; gridIndex < numOfGrids; ++gridIndex) {
        //     const TlPosition gridPosition(atomGridMat.get(gridIndex, GM_X),
        //                                   atomGridMat.get(gridIndex, GM_Y),
        //                                   atomGridMat.get(gridIndex, GM_Z));
        //     // calc phi table
        //     std::vector<WFGrid> phis;
        //     std::vector<WFGrid> gradPhisX;
        //     std::vector<WFGrid> gradPhisY;
        //     std::vector<WFGrid> gradPhisZ;
        //     this->getPhiTable(gridPosition, phis, gradPhisX, gradPhisY, gradPhisZ);
            
        //     // get rho at grid point
        //     double gridRhoA = 0.0;
        //     double gridGradRhoAX = 0.0;
        //     double gridGradRhoAY = 0.0;
        //     double gridGradRhoAZ = 0.0;
        //     this->getRhoAtGridPoint(P, phis, gradPhisX, gradPhisY, gradPhisZ,
        //                             &gridRhoA, &gridGradRhoAX, &gridGradRhoAY, &gridGradRhoAZ);
            
        //     rhoA[gridIndex] = gridRhoA;
        //     gradRhoAX[gridIndex] = gridGradRhoAX;
        //     gradRhoAY[gridIndex] = gridGradRhoAY;
        //     gradRhoAZ[gridIndex] = gridGradRhoAZ;
        // }


        TlMatrix Gx(numOfAOs, numOfAOs);
        TlMatrix Gy(numOfAOs, numOfAOs);
        TlMatrix Gz(numOfAOs, numOfAOs);
        
        std::vector<WFGrid> etas(numOfGrids);
        std::vector<WFGrid> gradEtasX(numOfGrids);
        std::vector<WFGrid> gradEtasY(numOfGrids);
        std::vector<WFGrid> gradEtasZ(numOfGrids);
        for (int gridIndex = 0; gridIndex < numOfGrids; ++gridIndex) {

            const TlPosition gridPosition(atomGridMat.get(gridIndex, GM_X),
                                          atomGridMat.get(gridIndex, GM_Y),
                                          atomGridMat.get(gridIndex, GM_Z));
            this->getPhiTable(gridPosition, etas, gradEtasX, gradEtasY, gradEtasZ);
            
            const double weight = atomGridMat.get(gridIndex, GM_WEIGHT);
            const double gridRhoA = atomGridMat.get(gridIndex, GM_LDA_RHO_ALPHA);
            // const double gridRhoA = rhoA[gridIndex];

            double roundF_roundRhoA;
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

        for (index_type p = 0; p < numOfAOs; ++p) {
            double x = 0.0;
            double y = 0.0;
            double z = 0.0;
            for (index_type q = 0; q < numOfAOs; ++q) {
                const double Ppq = P.get(p, q);
                x += Ppq * Gx.get(p, q);
                y += Ppq * Gy.get(p, q);
                z += Ppq * Gz.get(p, q);
            }

            pGX->set(p, atomIndex, 2.0 * x);
            pGY->set(p, atomIndex, 2.0 * y);
            pGZ->set(p, atomIndex, 2.0 * z);
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
    const Fl_Geometry flGeom((*this->pPdfParam_)["coordinates"]);
    TlMatrix gridMat = this->getGridMatrix<TlMatrix>(this->m_nIteration);

    const int numOfAtoms = this->m_nNumOfAtoms;
    const int numOfAOs = this->m_nNumOfAOs;
    
    pGX->resize(numOfAOs, numOfAtoms);
    pGY->resize(numOfAOs, numOfAtoms);
    pGZ->resize(numOfAOs, numOfAtoms);

    for (int atomIndex = 0; atomIndex < numOfAtoms; ++atomIndex) {
        const int atomicNumber = TlPrdctbl::getAtomicNumber(flGeom.getAtomSymbol(atomIndex));
        if (atomicNumber == 0) {
            // pass for dummy charge(X)
            continue;
        }

        // get grid information
        TlMatrix atomGridMat = this->selectGridMatrixByAtom(gridMat, atomIndex);
        const int numOfGrids = atomGridMat.getNumOfRows();
        std::vector<double> rhoA(numOfGrids, 0.0);
        std::vector<double> gradRhoAX(numOfGrids, 0.0);
        std::vector<double> gradRhoAY(numOfGrids, 0.0);
        std::vector<double> gradRhoAZ(numOfGrids, 0.0);
        if (this->m_bIsUpdateXC == true) {
            // すでにSCF計算中に電子密度を計算済み
            for (index_type i = 0; i < numOfGrids; ++i) {
                rhoA[i] = atomGridMat.get(i, GM_GGA_RHO_ALPHA);
                gradRhoAX[i] = atomGridMat.get(i, GM_GGA_GRAD_RHO_X_ALPHA);
                gradRhoAY[i] = atomGridMat.get(i, GM_GGA_GRAD_RHO_Y_ALPHA);
                gradRhoAZ[i] = atomGridMat.get(i, GM_GGA_GRAD_RHO_Z_ALPHA);
            }
        } else {
            // 電子密度を計算し直す
            for (int gridIndex = 0; gridIndex < numOfGrids; ++gridIndex) {
                //const TlPosition gridPosition = grids[gridIndex].position;
                TlPosition gridPosition(atomGridMat.get(gridIndex, GM_X),
                                        atomGridMat.get(gridIndex, GM_Y),
                                        atomGridMat.get(gridIndex, GM_Z));
                
                // calc phi table
                std::vector<WFGrid> phis;
                std::vector<WFGrid> gradPhisX;
                std::vector<WFGrid> gradPhisY;
                std::vector<WFGrid> gradPhisZ;
                this->getPhiTable(gridPosition, phis, gradPhisX, gradPhisY, gradPhisZ);
                
                // get rho at grid point
                double gridRhoA = 0.0;
                double gridGradRhoAX = 0.0;
                double gridGradRhoAY = 0.0;
                double gridGradRhoAZ = 0.0;
                this->getRhoAtGridPoint(P, phis, gradPhisX, gradPhisY, gradPhisZ,
                                        &gridRhoA, &gridGradRhoAX, &gridGradRhoAY, &gridGradRhoAZ);

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
            const TlPosition gridPosition(atomGridMat.get(gridIndex, GM_X),
                                          atomGridMat.get(gridIndex, GM_Y),
                                          atomGridMat.get(gridIndex, GM_Z));
            this->getPhiTable(gridPosition, etas, gradEtasX, gradEtasY, gradEtasZ);
            
            const double weight = atomGridMat.get(gridIndex, GM_WEIGHT);
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

TlMatrix DfCalcGridX::selectGridMatrixByAtom(const TlMatrix& globalGridMat,
                                             const int atomIndex)
{
    const index_type numOfGrids = globalGridMat.getNumOfRows();
    const double atomIndex_real = double(atomIndex);
    std::vector<index_type> finder;
    for (index_type i = 0; i < numOfGrids; ++i) {
        const double validation = std::fabs(globalGridMat.get(i, GM_ATOM_INDEX) - atomIndex_real);
        if (validation < 1.0E-5) {
            finder.push_back(i);
        }
    }

    const index_type numOfFinds = finder.size();
    const index_type numOfCols = globalGridMat.getNumOfCols();
    TlMatrix answer(numOfFinds, numOfCols);
    for (index_type i = 0; i < numOfFinds; ++i) {
        const index_type globalRow = finder[i];
        for (index_type col = 0; col < numOfCols; ++col) {
            const double value = globalGridMat.get(globalRow, col);
            answer.set(i, col, value);
        }
    }

    return answer;
}


TlMatrix DfCalcGridX::energyGradient(const TlSymmetricMatrix& P_A,
                                     DfFunctional_LDA* pFunctional)
{
    this->log_.info("DfCalcGridX::energyGradient() start");
    const int numOfAtoms = this->m_nNumOfAtoms;
    const int numOfAOs = this->m_nNumOfAOs;
    const TlOrbitalInfo orbitalInfo((*this->pPdfParam_)["coordinates"],
                                    (*this->pPdfParam_)["basis_set"]);
    const double densityCutOffValue = this->m_densityCutOffValueA;

    DfGenerateGrid dfGenGrid(this->pPdfParam_);

    TlMatrix gammaX(numOfAOs, numOfAOs);
    TlMatrix gammaY(numOfAOs, numOfAOs);
    TlMatrix gammaZ(numOfAOs, numOfAOs);

    double ene_xc = 0.0;
    this->log_.info("DfCalcGridX::energyGradient() Fxc1");
    TlMatrix Fxc1(numOfAtoms, 3); // weight derivative term
    for (int atomIndexA = 0; atomIndexA < numOfAtoms; ++atomIndexA) {
        std::vector<TlPosition> grids;
        std::vector<double> singleCenterWeights;
        std::vector<double> partitioningWeights;
        dfGenGrid.getGrids(atomIndexA, &grids, &singleCenterWeights, &partitioningWeights);

        const int numOfGrids = grids.size();
        assert(static_cast<int>(singleCenterWeights.size()) == numOfGrids);
        assert(static_cast<int>(partitioningWeights.size()) == numOfGrids);

        TlMatrix tmpGammaX(numOfAOs, numOfAOs);
        TlMatrix tmpGammaY(numOfAOs, numOfAOs);
        TlMatrix tmpGammaZ(numOfAOs, numOfAOs);
        this->log_.info(TlUtils::format("DfCalcGridX::energyGradient() atomindex %d/%d", atomIndexA, numOfAtoms));
        for (int gridIndex = 0; gridIndex < numOfGrids; ++gridIndex) {
            const TlPosition gridPosition = grids[gridIndex];
            const double singleCenterWeight = singleCenterWeights[gridIndex];
            const double partitioningWeight = partitioningWeights[gridIndex];
            
            // calc AO
            std::vector<WFGrid> AO_values;
            std::vector<WFGrid> dAO_dx_values;
            std::vector<WFGrid> dAO_dy_values;
            std::vector<WFGrid> dAO_dz_values;
            this->getPhiTable(gridPosition, AO_values, dAO_dx_values, dAO_dy_values, dAO_dz_values);

            // calc rho
            double rhoA = 0.0;
            this->getRhoAtGridPoint(P_A, AO_values, &rhoA);
            
            if (rhoA > densityCutOffValue) {
                // calc energy
                const double ene = pFunctional->getFunctional(rhoA, rhoA);
                ene_xc += singleCenterWeight * partitioningWeight * ene;
            
                // calc func derivative
                double dF_dRhoA = 0.0;
                pFunctional->getDerivativeFunctional(rhoA, &dF_dRhoA);
                
                // functional derivative term
                const double wAi_dF_dRhoA = partitioningWeight * singleCenterWeight * dF_dRhoA;
                const int numOf_AO_values = AO_values.size();
                const int numOf_dAO_dx_values = dAO_dx_values.size();
                const int numOf_dAO_dy_values = dAO_dy_values.size();
                const int numOf_dAO_dz_values = dAO_dz_values.size();
                for (int i = 0; i < numOf_AO_values; ++i) {
                    const double AO_value = AO_values[i].value;
                    const index_type q = AO_values[i].index;
                    const double coef = wAi_dF_dRhoA * AO_value;
                    
                    for (int j = 0; j < numOf_dAO_dx_values; ++j) {
                        const index_type p = dAO_dx_values[j].index;
                        tmpGammaX(p, q) += coef * dAO_dx_values[j].value;
                    }
                    for (int j = 0; j < numOf_dAO_dy_values; ++j) {
                        const index_type p = dAO_dy_values[j].index;
                        tmpGammaY(p, q) += coef * dAO_dy_values[j].value;
                    }
                    for (int j = 0; j < numOf_dAO_dz_values; ++j) {
                        const index_type p = dAO_dz_values[j].index;
                        tmpGammaZ(p, q) += coef * dAO_dz_values[j].value;
                    }
                }
                
                // weight derivative term
                for (int atomIndexB = 0; atomIndexB < numOfAtoms; ++atomIndexB) {
                    if (atomIndexA != atomIndexB) {
                        const TlVector nablaB_omegaA = dfGenGrid.JGP_nablaB_omegaA(atomIndexA, atomIndexB, gridPosition);
                        const TlVector E_nablaB_omegaAi = nablaB_omegaA * singleCenterWeight * ene;
                        
                        for (int i = 0; i < 3; ++i) {
                            Fxc1(atomIndexB, i) += E_nablaB_omegaAi[i];
                        }
                    }
                }
            }
        }

        // Gamma matrix
        {
            tmpGammaX.dot(P_A);
            tmpGammaY.dot(P_A);
            tmpGammaZ.dot(P_A);
            for (index_type p = 0; p < numOfAOs; ++p) {
                const double x = tmpGammaX.getRowVector(p).sum();
                const double y = tmpGammaY.getRowVector(p).sum();
                const double z = tmpGammaZ.getRowVector(p).sum();
        
                gammaX(p, atomIndexA) = 2.0 * x;
                gammaY(p, atomIndexA) = 2.0 * y;
                gammaZ(p, atomIndexA) = 2.0 * z;
            }
        }
    }
    this->log_.info(TlUtils::format("XC ene = % 16.10f", ene_xc));

    this->log_.info("DfCalcGridX::energyGradient() Fxc2");
    TlMatrix Fxc2(numOfAtoms, 3); // functional derivative term formed by GammaMatrix
    for (int mu = 0; mu < numOfAtoms; ++mu) {
        for (index_type alpha = 0; alpha < numOfAOs; ++alpha) {
            const index_type orbAtomId = orbitalInfo.getAtomIndex(alpha);
            for (int nu = 0; nu < numOfAtoms; ++nu) {
                if (orbAtomId != nu) {
                    if (mu != nu) {
                        const double fx = gammaX.get(alpha, nu);
                        const double fy = gammaY.get(alpha, nu);
                        const double fz = gammaZ.get(alpha, nu);
                        
                        Fxc2.add(mu, 0, -fx); // 0 means x
                        Fxc2.add(mu, 1, -fy); // 1 means y
                        Fxc2.add(mu, 2, -fz); // 2 means z
                    } else {
                        const double fx = gammaX.get(alpha, nu);
                        const double fy = gammaY.get(alpha, nu);
                        const double fz = gammaZ.get(alpha, nu);
                        
                        Fxc2.add(mu, 0,  fx); // 0 means x
                        Fxc2.add(mu, 1,  fy); // 1 means y
                        Fxc2.add(mu, 2,  fz); // 2 means z
                    }
                }
            }
        }
    }

    this->log_.info("DfCalcGridX::energyGradient() end");
    // Fxc1.save("Fxc1.mat");
    // Fxc2.save("Fxc2.mat");
    return Fxc1 + Fxc2;
}

TlMatrix DfCalcGridX::energyGradient(const TlSymmetricMatrix& P_A,
                                     DfFunctional_GGA* pFunctional)
{
    this->log_.info("DfCalcGridX::energyGradient() start");
    const int numOfAtoms = this->m_nNumOfAtoms;
    const int numOfAOs = this->m_nNumOfAOs;
    const TlOrbitalInfo orbitalInfo((*this->pPdfParam_)["coordinates"],
                                    (*this->pPdfParam_)["basis_set"]);
    const double densityCutOffValue = this->m_densityCutOffValueA;
    const double cutOffValue = this->m_inputedCutoffThreshold;

    DfGenerateGrid dfGenGrid(this->pPdfParam_);

    TlMatrix gammaX(numOfAOs, numOfAOs);
    TlMatrix gammaY(numOfAOs, numOfAOs);
    TlMatrix gammaZ(numOfAOs, numOfAOs);

    double ene_xc = 0.0;
    this->log_.info("DfCalcGridX::energyGradient() Fxc1");
    TlMatrix Fxc1(numOfAtoms, 3); // weight derivative term
    for (int atomIndexA = 0; atomIndexA < numOfAtoms; ++atomIndexA) {
        std::vector<TlPosition> grids;
        std::vector<double> singleCenterWeights;
        std::vector<double> partitioningWeights;
        dfGenGrid.getGrids(atomIndexA, &grids, &singleCenterWeights, &partitioningWeights);

        const int numOfGrids = grids.size();
        assert(static_cast<int>(singleCenterWeights.size()) == numOfGrids);
        assert(static_cast<int>(partitioningWeights.size()) == numOfGrids);

        TlMatrix tmpGammaX(numOfAOs, numOfAOs);
        TlMatrix tmpGammaY(numOfAOs, numOfAOs);
        TlMatrix tmpGammaZ(numOfAOs, numOfAOs);
        this->log_.info(TlUtils::format("DfCalcGridX::energyGradient() atomindex %d/%d", atomIndexA, numOfAtoms));
        for (int gridIndex = 0; gridIndex < numOfGrids; ++gridIndex) {
            const TlPosition gridPosition = grids[gridIndex];
            const double singleCenterWeight = singleCenterWeights[gridIndex];
            const double partitioningWeight = partitioningWeights[gridIndex];
            
            // calc AO
            std::vector<WFGrid> AO_values;
            std::vector<WFGrid> dAO_dx_values;
            std::vector<WFGrid> dAO_dy_values;
            std::vector<WFGrid> dAO_dz_values;
            this->getPhiTable(gridPosition, AO_values, dAO_dx_values, dAO_dy_values, dAO_dz_values);

            std::vector<WFGrid> d2AO_dxdx_values;
            std::vector<WFGrid> d2AO_dxdy_values;
            std::vector<WFGrid> d2AO_dxdz_values;
            std::vector<WFGrid> d2AO_dydy_values;
            std::vector<WFGrid> d2AO_dydz_values;
            std::vector<WFGrid> d2AO_dzdz_values;
            this->getDerivative2AOs(gridPosition,
                                    &d2AO_dxdx_values, &d2AO_dxdy_values, &d2AO_dxdz_values,
                                    &d2AO_dydy_values, &d2AO_dydz_values, &d2AO_dzdz_values);
            
            // calc rho
            double rhoA = 0.0;
            double gradRhoAX = 0.0;
            double gradRhoAY = 0.0;
            double gradRhoAZ = 0.0;
            this->getRhoAtGridPoint(P_A, AO_values, dAO_dx_values, dAO_dy_values, dAO_dz_values,
                                    &rhoA, &gradRhoAX, &gradRhoAY, &gradRhoAZ);
            const double& rhoB = rhoA;
            const double& gradRhoBX = gradRhoAX;
            const double& gradRhoBY = gradRhoAY;
            const double& gradRhoBZ = gradRhoAZ;
            
            if (rhoA > densityCutOffValue) {
                const double gammaAA = gradRhoAX*gradRhoAX + gradRhoAY*gradRhoAY + gradRhoAZ*gradRhoAZ;
                const double& gammaAB = gammaAA;
                const double& gammaBB = gammaAA;

                // calc energy
                const double ene = pFunctional->getFunctional(rhoA, rhoB, gammaAA, gammaAB, gammaBB);
                ene_xc += singleCenterWeight * partitioningWeight * ene;
            
                // calc func derivative
                double dF_dRhoA = 0.0;
                double dF_dRhoB = 0.0;
                double dF_dGammaAA = 0.0;
                double dF_dGammaAB = 0.0;
                double dF_dGammaBB = 0.0;
                pFunctional->getDerivativeFunctional(rhoA, rhoB,
                                                     gammaAA, gammaAB, gammaBB,
                                                     &dF_dRhoA, &dF_dRhoB,
                                                     &dF_dGammaAA, &dF_dGammaAB, &dF_dGammaBB);
                
                // functional derivative term(LDA)
                const double w = partitioningWeight * singleCenterWeight;

                const double wAi_dF_dRhoA = w * dF_dRhoA;
                const int numOf_AO_values = AO_values.size();
                const int numOf_dAO_dx_values = dAO_dx_values.size();
                const int numOf_dAO_dy_values = dAO_dy_values.size();
                const int numOf_dAO_dz_values = dAO_dz_values.size();
                for (int i = 0; i < numOf_AO_values; ++i) {
                    const double AO_value = AO_values[i].value;
                    const index_type q = AO_values[i].index;
                    const double coef = wAi_dF_dRhoA * AO_value;
                    
                    for (int j = 0; j < numOf_dAO_dx_values; ++j) {
                        const index_type p = dAO_dx_values[j].index;
                        tmpGammaX(p, q) += coef * dAO_dx_values[j].value;
                    }
                    for (int j = 0; j < numOf_dAO_dy_values; ++j) {
                        const index_type p = dAO_dy_values[j].index;
                        tmpGammaY(p, q) += coef * dAO_dy_values[j].value;
                    }
                    for (int j = 0; j < numOf_dAO_dz_values; ++j) {
                        const index_type p = dAO_dz_values[j].index;
                        tmpGammaZ(p, q) += coef * dAO_dz_values[j].value;
                    }
                }

                // functional derivative term(GGA) 
                {
                    const double coef_ax = 2.0 * dF_dGammaAA * gradRhoAX + dF_dGammaAB * gradRhoBX;
                    const double coef_ay = 2.0 * dF_dGammaAA * gradRhoAY + dF_dGammaAB * gradRhoBY;
                    const double coef_az = 2.0 * dF_dGammaAA * gradRhoAZ + dF_dGammaAB * gradRhoBZ;

                    const int numOf_d2AO_dxdx_values = d2AO_dxdx_values.size();
                    const int numOf_d2AO_dxdy_values = d2AO_dxdy_values.size();
                    const int numOf_d2AO_dxdz_values = d2AO_dxdz_values.size();
                    const int numOf_d2AO_dydy_values = d2AO_dydy_values.size();
                    const int numOf_d2AO_dydz_values = d2AO_dydz_values.size();
                    const int numOf_d2AO_dzdz_values = d2AO_dzdz_values.size();

                    TlMatrix Xx(numOfAOs, numOfAOs);
                    TlMatrix Xy(numOfAOs, numOfAOs);
                    TlMatrix Xz(numOfAOs, numOfAOs);
                    for (int j = 0; j < numOf_AO_values; ++j) {
                        const int v = AO_values[j].index;
                        const double v_v = AO_values[j].value;

                        for (int i = 0; i < numOf_d2AO_dxdx_values; ++i) {
                            const int u = d2AO_dxdx_values[i].index;
                            const double v_u = d2AO_dxdx_values[i].value;
                            Xx(u, v) += coef_ax * v_u * v_v;
                        }
                        for (int i = 0; i < numOf_d2AO_dxdy_values; ++i) {
                            const int u = d2AO_dxdy_values[i].index;
                            const double v_u = d2AO_dxdy_values[i].value;
                            Xx(u, v) += coef_ay * v_u * v_v;
                        }
                        for (int i = 0; i < numOf_d2AO_dxdz_values; ++i) {
                            const int u = d2AO_dxdz_values[i].index;
                            const double v_u = d2AO_dxdz_values[i].value;
                            Xx(u, v) += coef_az * v_u * v_v;
                        }

                        for (int i = 0; i < numOf_d2AO_dxdy_values; ++i) {
                            const int u = d2AO_dxdy_values[i].index;
                            const double v_u = d2AO_dxdy_values[i].value;
                            Xy(u, v) += coef_ax * v_u * v_v;
                        }
                        for (int i = 0; i < numOf_d2AO_dydy_values; ++i) {
                            const int u = d2AO_dydy_values[i].index;
                            const double v_u = d2AO_dydy_values[i].value;
                            Xy(u, v) += coef_ay * v_u * v_v;
                        }
                        for (int i = 0; i < numOf_d2AO_dydz_values; ++i) {
                            const int u = d2AO_dydz_values[i].index;
                            const double v_u = d2AO_dydz_values[i].value;
                            Xy(u, v) += coef_az * v_u * v_v;
                        }

                        for (int i = 0; i < numOf_d2AO_dxdz_values; ++i) {
                            const int u = d2AO_dxdz_values[i].index;
                            const double v_u = d2AO_dxdz_values[i].value;
                            Xz(u, v) += coef_ax * v_u * v_v;
                        }
                        for (int i = 0; i < numOf_d2AO_dydz_values; ++i) {
                            const int u = d2AO_dydz_values[i].index;
                            const double v_u = d2AO_dydz_values[i].value;
                            Xz(u, v) += coef_ay * v_u * v_v;
                        }
                        for (int i = 0; i < numOf_d2AO_dzdz_values; ++i) {
                            const int u = d2AO_dzdz_values[i].index;
                            const double v_u = d2AO_dzdz_values[i].value;
                            Xz(u, v) += coef_az * v_u * v_v;
                        }
                    }

                    for (int i = 0; i < numOf_dAO_dx_values; ++i) {
                        const int u = dAO_dx_values[i].index;
                        const double v_u = dAO_dx_values[i].value;

                        for (int j = 0; j < numOf_dAO_dx_values; ++j) {
                            const int v = dAO_dx_values[j].index;
                            const double v_v = dAO_dx_values[j].value;
                            Xx(u, v) += coef_ax * v_u * v_v;
                        }
                        for (int j = 0; j < numOf_dAO_dy_values; ++j) {
                            const int v = dAO_dy_values[j].index;
                            const double v_v = dAO_dy_values[j].value;
                            Xx(u, v) += coef_ay * v_u * v_v;
                        }
                        for (int j = 0; j < numOf_dAO_dz_values; ++j) {
                            const int v = dAO_dz_values[j].index;
                            const double v_v = dAO_dz_values[j].value;
                            Xx(u, v) += coef_az * v_u * v_v;
                        }
                    }

                    for (int i = 0; i < numOf_dAO_dy_values; ++i) {
                        const int u = dAO_dy_values[i].index;
                        const double v_u = dAO_dy_values[i].value;

                        for (int j = 0; j < numOf_dAO_dx_values; ++j) {
                            const int v = dAO_dx_values[j].index;
                            const double v_v = dAO_dx_values[j].value;
                            Xy(u, v) += coef_ax * v_u * v_v;
                        }
                        for (int j = 0; j < numOf_dAO_dy_values; ++j) {
                            const int v = dAO_dy_values[j].index;
                            const double v_v = dAO_dy_values[j].value;
                            Xy(u, v) += coef_ay * v_u * v_v;
                        }
                        for (int j = 0; j < numOf_dAO_dz_values; ++j) {
                            const int v = dAO_dz_values[j].index;
                            const double v_v = dAO_dz_values[j].value;
                            Xy(u, v) += coef_az * v_u * v_v;
                        }
                    }

                    for (int i = 0; i < numOf_dAO_dz_values; ++i) {
                        const int u = dAO_dz_values[i].index;
                        const double v_u = dAO_dz_values[i].value;

                        for (int j = 0; j < numOf_dAO_dx_values; ++j) {
                            const int v = dAO_dx_values[j].index;
                            const double v_v = dAO_dx_values[j].value;
                            Xz(u, v) += coef_ax * v_u * v_v;
                        }
                        for (int j = 0; j < numOf_dAO_dy_values; ++j) {
                            const int v = dAO_dy_values[j].index;
                            const double v_v = dAO_dy_values[j].value;
                            Xz(u, v) += coef_ay * v_u * v_v;
                        }
                        for (int j = 0; j < numOf_dAO_dz_values; ++j) {
                            const int v = dAO_dz_values[j].index;
                            const double v_v = dAO_dz_values[j].value;
                            Xz(u, v) += coef_az * v_u * v_v;
                        }
                    }

                    {
                        tmpGammaX += w * Xx;
                        tmpGammaY += w * Xy;
                        tmpGammaZ += w * Xz;
                    }
                }
                
                // weight derivative term
                for (int atomIndexB = 0; atomIndexB < numOfAtoms; ++atomIndexB) {
                    if (atomIndexA != atomIndexB) {
                        const TlVector nablaB_omegaA = dfGenGrid.JGP_nablaB_omegaA(atomIndexA, atomIndexB, gridPosition);
                        const TlVector E_nablaB_omegaAi = nablaB_omegaA * singleCenterWeight * ene;
                        
                        for (int i = 0; i < 3; ++i) {
                            Fxc1(atomIndexB, i) += E_nablaB_omegaAi[i];
                        }
                    }
                }
            }
        }

        // Gamma matrix
        {
            tmpGammaX.dot(P_A);
            tmpGammaY.dot(P_A);
            tmpGammaZ.dot(P_A);
            for (index_type p = 0; p < numOfAOs; ++p) {
                const double x = tmpGammaX.getRowVector(p).sum();
                const double y = tmpGammaY.getRowVector(p).sum();
                const double z = tmpGammaZ.getRowVector(p).sum();
        
                gammaX(p, atomIndexA) = 2.0 * x;
                gammaY(p, atomIndexA) = 2.0 * y;
                gammaZ(p, atomIndexA) = 2.0 * z;
            }
        }
    }
    this->log_.info(TlUtils::format("XC ene = % 16.10f", ene_xc));

    this->log_.info("DfCalcGridX::energyGradient() Fxc2");
    TlMatrix Fxc2(numOfAtoms, 3); // functional derivative term formed by GammaMatrix
    for (int mu = 0; mu < numOfAtoms; ++mu) {
        for (index_type alpha = 0; alpha < numOfAOs; ++alpha) {
            const index_type orbAtomId = orbitalInfo.getAtomIndex(alpha);
            for (int nu = 0; nu < numOfAtoms; ++nu) {
                if (orbAtomId != nu) {
                    if (mu != nu) {
                        const double fx = gammaX.get(alpha, nu);
                        const double fy = gammaY.get(alpha, nu);
                        const double fz = gammaZ.get(alpha, nu);
                        
                        Fxc2.add(mu, 0, -fx); // 0 means x
                        Fxc2.add(mu, 1, -fy); // 1 means y
                        Fxc2.add(mu, 2, -fz); // 2 means z
                    } else {
                        const double fx = gammaX.get(alpha, nu);
                        const double fy = gammaY.get(alpha, nu);
                        const double fz = gammaZ.get(alpha, nu);
                        
                        Fxc2.add(mu, 0,  fx); // 0 means x
                        Fxc2.add(mu, 1,  fy); // 1 means y
                        Fxc2.add(mu, 2,  fz); // 2 means z
                    }
                }
            }
        }
    }

    this->log_.info("DfCalcGridX::energyGradient() end");
    Fxc1.save("Fxc1.mat");
    Fxc2.save("Fxc2.mat");
    return Fxc1 + Fxc2;
}

// experimental code -----------------------------------------------------------
void DfCalcGridX::calcRhoVals_LDA(const std::vector<index_type>& P_rowIndexes,
                                  const std::vector<index_type>& P_colIndexes,
                                  const TlMatrix& P,
                                  TlMatrix* pGridMatrix)
{
    assert(pGridMatrix != NULL);
    
    // const std::size_t numOfRows = P.getNumOfRows();
    // const std::size_t numOfCols = P.getNumOfCols();
    assert((index_type)P_rowIndexes.size() == P.getNumOfRows());
    assert((index_type)P_colIndexes.size() == P.getNumOfCols());
    // const int calcMode = pGridMatrix->getNumOfCols();
    
    const index_type numOfGrids = pGridMatrix->getNumOfRows();
#pragma omp parallel for schedule(runtime)
    for (int grid = 0; grid < numOfGrids; ++grid) {
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
            assert(coefMatrix.getNumOfRows() == P.getNumOfRows());
            assert(coefMatrix.getNumOfCols() == P.getNumOfCols());
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
    
    // const std::size_t numOfRows = P.getNumOfRows();
    // const std::size_t numOfCols = P.getNumOfCols();
    assert((index_type)P_rowIndexes.size() == P.getNumOfRows());
    assert((index_type)P_colIndexes.size() == P.getNumOfCols());
    // const int calcMode = pGridMatrix->getNumOfCols();
    
    const index_type numOfGrids = pGridMatrix->getNumOfRows();
#pragma omp parallel for schedule(runtime)
    for (int grid = 0; grid < numOfGrids; ++grid) {
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
            assert(wf_rc.getNumOfRows() == P.getNumOfRows());
            assert(wf_rc.getNumOfCols() == P.getNumOfCols());
            wf_rc.dot(P);
            pGridMatrix->add(grid, 4, wf_rc.sum());
        }
        {
            TlMatrix wf_rc = wf_row * wf_dX;
            assert(wf_rc.getNumOfRows() == P.getNumOfRows());
            assert(wf_rc.getNumOfCols() == P.getNumOfCols());
            wf_rc.dot(P);
            pGridMatrix->add(grid, 5, 2.0 * wf_rc.sum());
        }
        {
            TlMatrix wf_rc = wf_row * wf_dY;
            assert(wf_rc.getNumOfRows() == P.getNumOfRows());
            assert(wf_rc.getNumOfCols() == P.getNumOfCols());
            wf_rc.dot(P);
            pGridMatrix->add(grid, 6, 2.0 * wf_rc.sum());
        }
        {
            TlMatrix wf_rc = wf_row * wf_dZ;
            assert(wf_rc.getNumOfRows() == P.getNumOfRows());
            assert(wf_rc.getNumOfCols() == P.getNumOfCols());
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


// double DfCalcGridX::buildK(const TlMatrix& gridMatrix,
//                            DfFunctional_LDA* pFunctional,
//                            TlSymmetricMatrix* pF)
// {
//     const index_type numOfAOs = this->m_nNumOfAOs;
//     const double densityCutOffValue = this->m_densityCutOffValueA;

//     double energy = 0.0;
//     const std::size_t numOfAllGrids = gridMatrix.getNumOfRows();
//     for (std::size_t grid = 0; grid < numOfAllGrids; ++grid) {
//         const double x = gridMatrix.get(grid, 0);
//         const double y = gridMatrix.get(grid, 1);
//         const double z = gridMatrix.get(grid, 2);
//         const double w = gridMatrix.get(grid, 3);
//         const TlPosition position(x, y, z);
//         const double rhoA = gridMatrix.get(grid, 4);

//         // calc phi table
//         std::vector<WFGrid> phis;
//         this->getPhiTable(position, 0, numOfAOs, phis);

//         if (rhoA > densityCutOffValue) {
//             this->buildFock(rhoA,
//                             phis,
//                             pFunctional, w, pF);
//             energy += w * pFunctional->getFunctional(rhoA);
//         }
//     }

//     return energy;
// }


// double DfCalcGridX::buildK_2(const TlMatrix& gridMatrix,
//                              const std::vector<index_type>& rowIndexes,
//                              const std::vector<index_type>& colIndexes,
//                              DfFunctional_LDA* pFunctional,
//                              TlSymmetricMatrix* pF)
// {
//     const double densityCutOffValue = this->m_densityCutOffValueA;

//     TlMatrix tmpF(*pF);
//     double energy = 0.0;
//     const std::size_t numOfAllGrids = gridMatrix.getNumOfRows();
//     for (std::size_t grid = 0; grid < numOfAllGrids; ++grid) {
//         const double x = gridMatrix.get(grid, 0);
//         const double y = gridMatrix.get(grid, 1);
//         const double z = gridMatrix.get(grid, 2);
//         const double w = gridMatrix.get(grid, 3);
//         const TlPosition r(x, y, z);
//         const double rhoA = gridMatrix.get(grid, 4);

//         if (rhoA > densityCutOffValue) {
//             TlMatrix wf_row;
//             this->getWaveFunctionValues(rowIndexes, r,
//                                         &wf_row);
//             TlMatrix wf_col;
//             this->getWaveFunctionValues(colIndexes, r,
//                                         &wf_col);
//             wf_col.transpose();
            
//             // build K
//             double roundF_roundRhoA;
//             pFunctional->getDerivativeFunctional(rhoA, &roundF_roundRhoA);
            
//             const double coef1_A = w * roundF_roundRhoA;
//             {
//                 TlMatrix wf_rc = wf_row * wf_col;
//                 wf_rc *= coef1_A;
//                 tmpF += wf_rc;
//             }
            
//             // energy
//             energy += w * pFunctional->getFunctional(rhoA); // RKS code
//         }
//     }

//     *pF = tmpF;

//     return energy;
// }


// double DfCalcGridX::buildK(const TlMatrix& gridMatrix,
//                            DfFunctional_GGA* pFunctional,
//                            TlSymmetricMatrix* pF)
// {
//     //const index_type numOfAOs = this->m_nNumOfAOs;
//     const double densityCutOffValue = this->m_densityCutOffValueA;

//     double energy = 0.0;
//     const std::size_t numOfAllGrids = gridMatrix.getNumOfRows();
//     for (std::size_t grid = 0; grid < numOfAllGrids; ++grid) {
//         const double x = gridMatrix.get(grid, 0);
//         const double y = gridMatrix.get(grid, 1);
//         const double z = gridMatrix.get(grid, 2);
//         const double w = gridMatrix.get(grid, 3);
//         const TlPosition position(x, y, z);
//         const double rhoA = gridMatrix.get(grid, 4);
//         const double gradRhoX_A = gridMatrix.get(grid, 5);
//         const double gradRhoY_A = gridMatrix.get(grid, 6);
//         const double gradRhoZ_A = gridMatrix.get(grid, 7);

//         // calc phi table
//         std::vector<WFGrid> aPhi;
//         std::vector<WFGrid> aGradPhiX;
//         std::vector<WFGrid> aGradPhiY;
//         std::vector<WFGrid> aGradPhiZ;
//         this->getPhiTable(position, aPhi, aGradPhiX, aGradPhiY, aGradPhiZ);

//         if (rhoA > densityCutOffValue) {
//             this->buildFock(rhoA, gradRhoX_A, gradRhoY_A, gradRhoZ_A,
//                             aPhi, aGradPhiX, aGradPhiY, aGradPhiZ,
//                             pFunctional, w, pF); // RKS code
//             const double gammaAA =  gradRhoX_A*gradRhoX_A + gradRhoY_A*gradRhoY_A + gradRhoZ_A*gradRhoZ_A;
//             energy += w * pFunctional->getFunctional(rhoA, gammaAA); // RKS code
//         }
//     }

//     return energy;
// }


// double DfCalcGridX::buildK_2(const TlMatrix& gridMatrix,
//                              const std::vector<index_type>& rowIndexes,
//                              const std::vector<index_type>& colIndexes,
//                              DfFunctional_GGA* pFunctional,
//                              TlMatrix* pF)
// {
//     const double densityCutOffValue = this->m_densityCutOffValueA;

//     TlMatrix tmpF(*pF);
//     double energy = 0.0;
//     const index_type numOfAllGrids = gridMatrix.getNumOfRows();
// #pragma omp parallel for schedule(runtime) reduction(+:energy)
//     for (int grid = 0; grid < numOfAllGrids; ++grid) {
//         const double x = gridMatrix.get(grid, 0);
//         const double y = gridMatrix.get(grid, 1);
//         const double z = gridMatrix.get(grid, 2);
//         const double w = gridMatrix.get(grid, 3);
//         const TlPosition r(x, y, z);
//         const double rhoA = gridMatrix.get(grid, 4);

//         if (rhoA > densityCutOffValue) {
//             const double gradRhoX_A = gridMatrix.get(grid, 5);
//             const double gradRhoY_A = gridMatrix.get(grid, 6);
//             const double gradRhoZ_A = gridMatrix.get(grid, 7);

//             TlMatrix wf_r, dfdx_r, dfdy_r, dfdz_r;
//             this->getWaveFunctionValues(rowIndexes, r,
//                                         &wf_r, &dfdx_r, &dfdy_r, &dfdz_r);
//             TlMatrix wf_c, dfdx_c, dfdy_c, dfdz_c;
//             this->getWaveFunctionValues(colIndexes, r,
//                                         &wf_c, &dfdx_c, &dfdy_c, &dfdz_c);
//             wf_c.transpose();
//             dfdx_c.transpose();
//             dfdy_c.transpose();
//             dfdz_c.transpose();
            
//             const double gammaAA =  gradRhoX_A*gradRhoX_A + gradRhoY_A*gradRhoY_A + gradRhoZ_A*gradRhoZ_A;
//             // build K
//             double rF_rR_A;
//             double rF_rG_AA, rF_rG_AB;
//             pFunctional->getDerivativeFunctional(rhoA, gammaAA,
//                                                  &rF_rR_A,
//                                                  &rF_rG_AA,
//                                                  &rF_rG_AB);
            
//             const double coef1_A = w * rF_rR_A;
//             const double rF_rG_AA2 = 2.0 * rF_rG_AA;
//             const double coef2_AX = w * (rF_rG_AA2 * gradRhoX_A + rF_rG_AB * gradRhoX_A);
//             const double coef2_AY = w * (rF_rG_AA2 * gradRhoY_A + rF_rG_AB * gradRhoY_A);
//             const double coef2_AZ = w * (rF_rG_AA2 * gradRhoZ_A + rF_rG_AB * gradRhoZ_A);

//             {
//                 TlMatrix wf_rc = wf_r * wf_c;
//                 wf_rc *= coef1_A;
//                 *pF += wf_rc;
//             }
//             {
//                 TlMatrix wf_rc = wf_r * dfdx_c;
//                 TlMatrix wf_cr = dfdx_r * wf_c;
//                 wf_rc *= coef2_AX;
//                 wf_cr *= coef2_AX;
//                 *pF += (wf_rc + wf_cr);
//             }
//             {
//                 TlMatrix wf_rc = wf_r * dfdy_c;
//                 TlMatrix wf_cr = dfdy_r * wf_c;
//                 wf_rc *= coef2_AY;
//                 wf_cr *= coef2_AY;
//                 *pF += (wf_rc + wf_cr);
//             }
//             {
//                 TlMatrix wf_rc = wf_r * dfdz_c;
//                 TlMatrix wf_cr = dfdz_r * wf_c;
//                 wf_rc *= coef2_AZ;
//                 wf_cr *= coef2_AZ;
//                 *pF += (wf_rc + wf_cr);
//             }
            
//             // energy
//             energy += w * pFunctional->getFunctional(rhoA, gammaAA); // RKS code
//         }
//     }

//     return energy;
// }

// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
void DfCalcGridX::calcRho_LDA(const TlSymmetricMatrix& P_A)
{
    TlMatrix gridMat = DfObject::getGridMatrix<TlMatrix>(this->m_nIteration -1);
    this->calcRho_LDA(P_A, &gridMat);
    DfObject::saveGridMatrix(this->m_nIteration, gridMat);
}

void DfCalcGridX::calcRho_LDA(const TlSymmetricMatrix& P_A,
                              const TlSymmetricMatrix& P_B)
{
    TlMatrix gridMat = DfObject::getGridMatrix<TlMatrix>(this->m_nIteration -1);
    this->calcRho_LDA(P_A, P_B, &gridMat);
    DfObject::saveGridMatrix(this->m_nIteration, gridMat);
}

void DfCalcGridX::calcRho_GGA(const TlSymmetricMatrix& P_A)
{
    TlMatrix gridMat = DfObject::getGridMatrix<TlMatrix>(this->m_nIteration -1);
    this->calcRho_GGA(P_A, &gridMat);
    DfObject::saveGridMatrix(this->m_nIteration, gridMat);
}

void DfCalcGridX::calcRho_GGA(const TlSymmetricMatrix& P_A,
                              const TlSymmetricMatrix& P_B)
{
    TlMatrix gridMat = DfObject::getGridMatrix<TlMatrix>(this->m_nIteration -1);
    this->calcRho_GGA(P_A, P_B, &gridMat);
    DfObject::saveGridMatrix(this->m_nIteration, gridMat);
}

void DfCalcGridX::calcRho_LDA(const TlSymmetricMatrix& P_A,
                              TlMatrix* pGridMat)
{
    assert(pGridMat != NULL);
    const index_type numOfGrids = pGridMat->getNumOfRows();
    assert(pGridMat->getNumOfCols() == GM_LDA_RHO_ALPHA +1);
    
    if (this->m_bIsUpdateXC != true) {
        // initialize
        TlMatrix zero(pGridMat->getNumOfRows(), 1);
        pGridMat->setBlockMatrix(0, GM_LDA_RHO_ALPHA, zero);
    }
    
#pragma omp parallel for schedule(runtime)
    for (index_type grid = 0; grid < numOfGrids; ++grid) {
        const TlPosition gridPosition(pGridMat->get(grid, GM_X),
                                      pGridMat->get(grid, GM_Y),
                                      pGridMat->get(grid, GM_Z));

        // calc phi table
        std::vector<WFGrid> phi;
        this->getPhiTable(gridPosition, phi);

        // get rho at grid point
        double rhoA = 0.0;
        this->getRhoAtGridPoint(P_A, phi, &rhoA);

        assert((0 <= grid) && (grid < pGridMat->getNumOfRows()));
        pGridMat->add(grid, GM_LDA_RHO_ALPHA, rhoA);
    }
}

void DfCalcGridX::calcRho_LDA(const TlSymmetricMatrix& P_A,
                              const TlSymmetricMatrix& P_B,
                              TlMatrix* pGridMat)
{
    const index_type numOfGrids = pGridMat->getNumOfRows();
    assert(pGridMat->getNumOfCols() == GM_LDA_RHO_BETA +1);

    if (this->m_bIsUpdateXC != true) {
        // initialize
        TlMatrix zero(pGridMat->getNumOfRows(), 2);
        pGridMat->setBlockMatrix(0, GM_LDA_RHO_ALPHA, zero);
    }
    
#pragma omp parallel for schedule(runtime)
    for (index_type grid = 0; grid < numOfGrids; ++grid) {
        const TlPosition gridPosition(pGridMat->get(grid, GM_X),
                                      pGridMat->get(grid, GM_Y),
                                      pGridMat->get(grid, GM_Z));

        // calc phi table
        std::vector<WFGrid> phi;
        this->getPhiTable(gridPosition, phi);

        // get rho at grid point
        double rhoA = 0.0;
        double rhoB = 0.0;
        this->getRhoAtGridPoint(P_A, phi, &rhoA);
        this->getRhoAtGridPoint(P_B, phi, &rhoB);

        pGridMat->add(grid, GM_LDA_RHO_ALPHA, rhoA);
        pGridMat->add(grid, GM_LDA_RHO_BETA,  rhoB);
    }
}

void DfCalcGridX::calcRho_GGA(const TlSymmetricMatrix& P_A,
                              TlMatrix* pGridMat)
{
    const index_type numOfGrids = pGridMat->getNumOfRows();
    assert(pGridMat->getNumOfCols() == GM_GGA_GRAD_RHO_Z_ALPHA +1);
    
    if (this->m_bIsUpdateXC != true) {
        // initialize
        TlMatrix zero(pGridMat->getNumOfRows(), 4);
        pGridMat->setBlockMatrix(0, GM_GGA_RHO_ALPHA, zero);
    }
    
#pragma omp parallel for schedule(runtime)
    for (index_type grid = 0; grid < numOfGrids; ++grid) {
        const TlPosition gridPosition(pGridMat->get(grid, GM_X),
                                      pGridMat->get(grid, GM_Y),
                                      pGridMat->get(grid, GM_Z));

        // calc phi table
        std::vector<WFGrid> phi;
        std::vector<WFGrid> gradPhiX;
        std::vector<WFGrid> gradPhiY;
        std::vector<WFGrid> gradPhiZ;
        this->getPhiTable(gridPosition, phi, gradPhiX, gradPhiY, gradPhiZ);
        
        // get rho at grid point
        double rhoA = 0.0;
        double gradRhoXA = 0.0;
        double gradRhoYA = 0.0;
        double gradRhoZA = 0.0;
        this->getRhoAtGridPoint(P_A, phi, gradPhiX, gradPhiY, gradPhiZ,
                                &rhoA, &gradRhoXA, &gradRhoYA, &gradRhoZA);

        pGridMat->add(grid, GM_GGA_RHO_ALPHA, rhoA);
        pGridMat->add(grid, GM_GGA_GRAD_RHO_X_ALPHA, gradRhoXA);
        pGridMat->add(grid, GM_GGA_GRAD_RHO_Y_ALPHA, gradRhoYA);
        pGridMat->add(grid, GM_GGA_GRAD_RHO_Z_ALPHA, gradRhoZA);
    }
}

void DfCalcGridX::calcRho_GGA(const TlSymmetricMatrix& P_A,
                              const TlSymmetricMatrix& P_B,
                              TlMatrix* pGridMat)
{
    const index_type numOfGrids = pGridMat->getNumOfRows();
    assert(pGridMat->getNumOfCols() == GM_GGA_GRAD_RHO_Z_BETA +1);

    if (this->m_bIsUpdateXC != true) {
        // initialize
        TlMatrix zero(pGridMat->getNumOfRows(), 8);
        pGridMat->setBlockMatrix(0, GM_GGA_RHO_ALPHA, zero);
    }
    
#pragma omp parallel for schedule(runtime)
    for (index_type grid = 0; grid < numOfGrids; ++grid) {
        const TlPosition gridPosition(pGridMat->get(grid, GM_X),
                                      pGridMat->get(grid, GM_Y),
                                      pGridMat->get(grid, GM_Z));

        // calc phi table
        std::vector<WFGrid> phi;
        std::vector<WFGrid> gradPhiX;
        std::vector<WFGrid> gradPhiY;
        std::vector<WFGrid> gradPhiZ;
        this->getPhiTable(gridPosition, phi, gradPhiX, gradPhiY, gradPhiZ);

        // get rho at grid point
        double rhoA = 0.0;
        double gradRhoXA = 0.0;
        double gradRhoYA = 0.0;
        double gradRhoZA = 0.0;
        double rhoB = 0.0;
        double gradRhoXB = 0.0;
        double gradRhoYB = 0.0;
        double gradRhoZB = 0.0;
        this->getRhoAtGridPoint(P_A, phi, gradPhiX, gradPhiY, gradPhiZ,
                                &rhoA, &gradRhoXA, &gradRhoYA, &gradRhoZA);
        this->getRhoAtGridPoint(P_B, phi, gradPhiX, gradPhiY, gradPhiZ,
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

double DfCalcGridX::buildVxc(DfFunctional_LDA* pFunctional,
                             TlSymmetricMatrix* pF_A)
{
    TlMatrix gridMat = DfObject::getGridMatrix<TlMatrix>(this->m_nIteration);
    double energy = this->buildVxc(gridMat, pFunctional, pF_A);
    return energy;
}

double DfCalcGridX::buildVxc(DfFunctional_LDA* pFunctional,
                             TlSymmetricMatrix* pF_A,
                             TlSymmetricMatrix* pF_B)
{
    TlMatrix gridMat = DfObject::getGridMatrix<TlMatrix>(this->m_nIteration);
    double energy = this->buildVxc(gridMat, pFunctional,
                                   pF_A, pF_B);
    return energy;
}

double DfCalcGridX::buildVxc(DfFunctional_GGA* pFunctional,
                             TlSymmetricMatrix* pF_A)
{
    TlMatrix gridMat = DfObject::getGridMatrix<TlMatrix>(this->m_nIteration);
    double energy = this->buildVxc(gridMat, pFunctional, pF_A);
    return energy;
}

double DfCalcGridX::buildVxc(DfFunctional_GGA* pFunctional,
                             TlSymmetricMatrix* pF_A,
                             TlSymmetricMatrix* pF_B)
{
    TlMatrix gridMat = DfObject::getGridMatrix<TlMatrix>(this->m_nIteration);
    double energy = this->buildVxc(gridMat, pFunctional,
                                   pF_A, pF_B);
    return energy;
}

double DfCalcGridX::buildVxc(const TlMatrix& gridMatrix,
                             DfFunctional_LDA* pFunctional,
                             TlMatrixObject* pF_A)
{
    double energy = 0.0;
    const double densityCutOffValue = this->m_densityCutOffValueA;
    const index_type numOfGrids = gridMatrix.getNumOfRows();
    
#pragma omp parallel for schedule(runtime) reduction(+:energy)
    for (index_type grid = 0; grid < numOfGrids; ++grid) {
        const TlPosition gridPosition(gridMatrix.get(grid, GM_X),
                                      gridMatrix.get(grid, GM_Y),
                                      gridMatrix.get(grid, GM_Z));
        const double weight = gridMatrix.get(grid, GM_WEIGHT);

        // calc phi table
        std::vector<WFGrid> phi;
        this->getPhiTable(gridPosition, phi);

        const double rhoA = gridMatrix.get(grid, GM_LDA_RHO_ALPHA);
        double roundF_roundRhoA = 0.0;

        if (rhoA > densityCutOffValue) {
            pFunctional->getDerivativeFunctional(rhoA,
                                                 &roundF_roundRhoA);
            
            this->build_XC_Matrix(roundF_roundRhoA, phi, pFunctional, weight, pF_A);
            energy += weight * pFunctional->getFunctional(rhoA, rhoA);
        }
    }

    return energy;
}

double DfCalcGridX::buildVxc(const TlMatrix& gridMatrix,
                             DfFunctional_LDA* pFunctional,
                             TlMatrixObject* pF_A,
                             TlMatrixObject* pF_B)
{
    double energy = 0.0;
    const double densityCutOffValue = std::min(this->m_densityCutOffValueA,
                                               this->m_densityCutOffValueB);
    const index_type numOfGrids = gridMatrix.getNumOfRows();
    
#pragma omp parallel for schedule(runtime) reduction(+:energy)
    for (index_type grid = 0; grid < numOfGrids; ++grid) {
        const TlPosition gridPosition(gridMatrix.get(grid, GM_X),
                                      gridMatrix.get(grid, GM_Y),
                                      gridMatrix.get(grid, GM_Z));
        const double weight = gridMatrix.get(grid, GM_WEIGHT);

        // calc phi table
        std::vector<WFGrid> phi;
        this->getPhiTable(gridPosition, phi);

        const double rhoA = gridMatrix.get(grid, GM_LDA_RHO_ALPHA);
        const double rhoB = gridMatrix.get(grid, GM_LDA_RHO_BETA);
        double roundF_roundRhoA = 0.0;
        double roundF_roundRhoB = 0.0;

        if ((rhoA > densityCutOffValue) || (rhoB > densityCutOffValue)) {
            pFunctional->getDerivativeFunctional(rhoA, rhoB,
                                                 &roundF_roundRhoA, &roundF_roundRhoB);
            
            this->build_XC_Matrix(roundF_roundRhoA, phi, pFunctional, weight, pF_A);
            this->build_XC_Matrix(roundF_roundRhoB, phi, pFunctional, weight, pF_B);
            energy += weight * pFunctional->getFunctional(rhoA, rhoB);
        }
    }

    return energy;
}

double DfCalcGridX::buildVxc(const TlMatrix& gridMatrix,
                             DfFunctional_GGA* pFunctional,
                             TlMatrixObject* pF_A)
{
    double energy = 0.0;
    const double densityCutOffValue = this->m_densityCutOffValueA;
    const index_type numOfGrids = gridMatrix.getNumOfRows();

#pragma omp parallel for schedule(runtime) reduction(+:energy)
    for (index_type grid = 0; grid < numOfGrids; ++grid) {
        const TlPosition gridPosition(gridMatrix.get(grid, GM_X),
                                      gridMatrix.get(grid, GM_Y),
                                      gridMatrix.get(grid, GM_Z));
        const double weight = gridMatrix.get(grid, GM_WEIGHT);
        
        // calc phi table
        std::vector<WFGrid> phi;
        std::vector<WFGrid> gradPhiX;
        std::vector<WFGrid> gradPhiY;
        std::vector<WFGrid> gradPhiZ;
        this->getPhiTable(gridPosition, phi, gradPhiX, gradPhiY, gradPhiZ);

        // get rho at grid point
        double rhoA = gridMatrix.get(grid, GM_GGA_RHO_ALPHA);
        double gradRhoXA = gridMatrix.get(grid, GM_GGA_GRAD_RHO_X_ALPHA);
        double gradRhoYA = gridMatrix.get(grid, GM_GGA_GRAD_RHO_Y_ALPHA);
        double gradRhoZA = gridMatrix.get(grid, GM_GGA_GRAD_RHO_Z_ALPHA);

        if (rhoA > densityCutOffValue) {
            const double gammaAA = gradRhoXA*gradRhoXA + gradRhoYA*gradRhoYA + gradRhoZA*gradRhoZA;
            double roundF_roundRhoA;
            double roundF_roundGammaAA;
            double roundF_roundGammaAB;
            pFunctional->getDerivativeFunctional(rhoA, gammaAA,
                                                 &roundF_roundRhoA,
                                                 &roundF_roundGammaAA,
                                                 &roundF_roundGammaAB);
            
            this->build_XC_Matrix(roundF_roundRhoA,
                                  roundF_roundGammaAA,
                                  roundF_roundGammaAB,
                                  gradRhoXA, gradRhoYA, gradRhoZA,
                                  phi, gradPhiX, gradPhiY, gradPhiZ,
                                  pFunctional, weight, pF_A);
            energy += weight * pFunctional->getFunctional(rhoA,
                                                          gammaAA);
        }
    }

    return energy;
}

double DfCalcGridX::buildVxc(const TlMatrix& gridMatrix,
                             DfFunctional_GGA* pFunctional,
                             TlMatrixObject* pF_A,
                             TlMatrixObject* pF_B)
{
    double energy = 0.0;
    const double densityCutOffValue = std::min(this->m_densityCutOffValueA,
                                               this->m_densityCutOffValueB);
    const index_type numOfGrids = gridMatrix.getNumOfRows();

#pragma omp parallel for schedule(runtime) reduction(+:energy)
    for (index_type grid = 0; grid < numOfGrids; ++grid) {
        const TlPosition gridPosition(gridMatrix.get(grid, GM_X),
                                      gridMatrix.get(grid, GM_Y),
                                      gridMatrix.get(grid, GM_Z));
        const double weight = gridMatrix.get(grid, GM_WEIGHT);
        
        // calc phi table
        std::vector<WFGrid> phi;
        std::vector<WFGrid> gradPhiX;
        std::vector<WFGrid> gradPhiY;
        std::vector<WFGrid> gradPhiZ;
        this->getPhiTable(gridPosition, phi, gradPhiX, gradPhiY, gradPhiZ);

        // get rho at grid point
        double rhoA = gridMatrix.get(grid, GM_GGA_RHO_ALPHA);
        double gradRhoXA = gridMatrix.get(grid, GM_GGA_GRAD_RHO_X_ALPHA);
        double gradRhoYA = gridMatrix.get(grid, GM_GGA_GRAD_RHO_Y_ALPHA);
        double gradRhoZA = gridMatrix.get(grid, GM_GGA_GRAD_RHO_Z_ALPHA);
        double rhoB = gridMatrix.get(grid, GM_GGA_RHO_BETA);
        double gradRhoXB = gridMatrix.get(grid, GM_GGA_GRAD_RHO_X_BETA);
        double gradRhoYB = gridMatrix.get(grid, GM_GGA_GRAD_RHO_X_BETA);
        double gradRhoZB = gridMatrix.get(grid, GM_GGA_GRAD_RHO_X_BETA);

        if ((rhoA > densityCutOffValue) || ((rhoB > densityCutOffValue))) {
            const double gammaAA = gradRhoXA*gradRhoXA + gradRhoYA*gradRhoYA + gradRhoZA*gradRhoZA;
            const double gammaAB = gradRhoXA*gradRhoXB + gradRhoYA*gradRhoYB + gradRhoZA*gradRhoZB;
            const double gammaBB = gradRhoXB*gradRhoXB + gradRhoYB*gradRhoYB + gradRhoZB*gradRhoZB;
            
            double roundF_roundRhoA, roundF_roundRhoB;
            double roundF_roundGammaAA, roundF_roundGammaAB, roundF_roundGammaBB;
            pFunctional->getDerivativeFunctional(rhoA, rhoB, gammaAA, gammaAB, gammaBB,
                                                 &roundF_roundRhoA, &roundF_roundRhoB,
                                                 &roundF_roundGammaAA,
                                                 &roundF_roundGammaAB,
                                                 &roundF_roundGammaBB);
            
            this->build_XC_Matrix(roundF_roundRhoA, roundF_roundGammaAA, roundF_roundGammaAB,
                                  gradRhoXA, gradRhoYA, gradRhoZA,
                                  phi, gradPhiX, gradPhiY, gradPhiZ,
                                  pFunctional, weight, pF_A);
            this->build_XC_Matrix(roundF_roundRhoB, roundF_roundGammaBB, roundF_roundGammaAB,
                                  gradRhoXB, gradRhoYB, gradRhoZB,
                                  phi, gradPhiX, gradPhiY, gradPhiZ,
                                  pFunctional, weight, pF_B);
            energy += weight * pFunctional->getFunctional(rhoA, rhoB,
                                                          gammaAA, gammaAB, gammaBB);
        }
    }

    return energy;
}

double DfCalcGridX::calcXCIntegForFockAndEnergy(const TlSymmetricMatrix& P_A,
                                                DfFunctional_LDA* pFunctional,
                                                TlSymmetricMatrix* pF_A)
{
    assert(pFunctional != NULL);
    assert(pF_A != NULL);

    this->calcRho_LDA(P_A);
    double energy = this->buildVxc(pFunctional, pF_A);

    return energy;
}

double DfCalcGridX::calcXCIntegForFockAndEnergy(const TlSymmetricMatrix& P_A,
                                                const TlSymmetricMatrix& P_B,
                                                DfFunctional_LDA* pFunctional,
                                                TlSymmetricMatrix* pF_A,
                                                TlSymmetricMatrix* pF_B)
{
    assert(pFunctional != NULL);
    assert(pF_A != NULL);
    assert(pF_B != NULL);

    this->calcRho_LDA(P_A, P_B);
    double energy = this->buildVxc(pFunctional, pF_A, pF_B);
    return energy;
}

double DfCalcGridX::calcXCIntegForFockAndEnergy(const TlSymmetricMatrix& P_A,
                                                DfFunctional_GGA* pFunctional,
                                                TlSymmetricMatrix* pF_A)
{
    assert(pFunctional != NULL);
    assert(pF_A != NULL);

    this->calcRho_GGA(P_A);
    double energy = this->buildVxc(pFunctional, pF_A);

    return energy;
}

double DfCalcGridX::calcXCIntegForFockAndEnergy(const TlSymmetricMatrix& P_A,
                                                const TlSymmetricMatrix& P_B,
                                                DfFunctional_GGA* pFunctional,
                                                TlSymmetricMatrix* pF_A,
                                                TlSymmetricMatrix* pF_B)
{
    assert(pFunctional != NULL);
    assert(pF_A != NULL);
    assert(pF_B != NULL);

    this->calcRho_GGA(P_A, P_B);
    double energy = this->buildVxc(pFunctional, pF_A, pF_B);
    return energy;
}

void DfCalcGridX::build_XC_Matrix(const double roundF_roundRhoA,
                                  const std::vector<WFGrid>& phi,
                                  DfFunctional_LDA* pFunctional, const double weight,
                                  TlMatrixObject* pF_A)
{
    const double coeff1_A = weight * roundF_roundRhoA;
    const double cutOffValue = this->m_inputedCutoffThreshold;

    if (phi.size() > 0) {
        const double phi_0 = phi[0].value;
        if (std::fabs(coeff1_A * phi_0 * phi_0) > cutOffValue) {
            std::vector<WFGrid>::const_iterator pEnd =
                std::upper_bound(phi.begin(), phi.end(),
                                 WFGrid(0, std::sqrt(std::fabs(cutOffValue/coeff1_A))),
                                 WFGrid_sort_functional());
            this->buildFock(phi.begin(), pEnd, coeff1_A, cutOffValue, pF_A);
        }
    }
}

void DfCalcGridX::build_XC_Matrix(const double roundF_roundRhoA,
                                  const double roundF_roundGammaAA,
                                  const double roundF_roundGammaAB,
                                  const double gradRhoAX, const double gradRhoAY, const double gradRhoAZ,
                                  const std::vector<WFGrid>& phi,
                                  const std::vector<WFGrid>& gradPhiX,
                                  const std::vector<WFGrid>& gradPhiY,
                                  const std::vector<WFGrid>& gradPhiZ,
                                  DfFunctional_GGA* pFunctional,
                                  const double weight,
                                  TlMatrixObject* pF_A)
{
    const double coeff1_A = weight * roundF_roundRhoA;
    const double roundF_roundGammaAA2 = 2.0 * roundF_roundGammaAA;
    const double coeff2_AX = weight
        * (roundF_roundGammaAA2 * gradRhoAX + roundF_roundGammaAB * gradRhoAX);
    const double coeff2_AY = weight
        * (roundF_roundGammaAA2 * gradRhoAY + roundF_roundGammaAB * gradRhoAY);
    const double coeff2_AZ = weight
        * (roundF_roundGammaAA2 * gradRhoAZ + roundF_roundGammaAB * gradRhoAZ);
    const double cutOffValue = this->m_inputedCutoffThreshold;

    if (phi.size() > 0) {
        const double phi_0 = phi[0].value;

        {
            // phi v.s. phi
            if (std::fabs(coeff1_A * phi_0 * phi_0) > cutOffValue) {
                std::vector<WFGrid>::const_iterator pEnd = std::upper_bound(phi.begin(), phi.end(),
                                                                            WFGrid(0, std::sqrt(std::fabs(cutOffValue / coeff1_A))),
                                                                            WFGrid_sort_functional());
                this->buildFock(phi.begin(), pEnd, coeff1_A, cutOffValue, pF_A);
            }
            // phi v.s. grad-phi(x)
            if ((gradPhiX.size() > 0) &&
                (std::fabs(coeff2_AX * phi_0 * gradPhiX[0].value) > cutOffValue)) {
                std::vector<WFGrid>::const_iterator pEnd = std::upper_bound(phi.begin(), phi.end(),
                                                                            WFGrid(0, cutOffValue / (coeff2_AX * gradPhiX[0].value)),
                                                                            WFGrid_sort_functional());
                this->buildFock(phi.begin(), pEnd, gradPhiX.begin(), gradPhiX.end(),
                                coeff2_AX, cutOffValue, pF_A);
                this->buildFock(gradPhiX.begin(), gradPhiX.end(), phi.begin(), pEnd,
                                coeff2_AX, cutOffValue, pF_A);
            }
            // phi v.s. grad-phi(y)
            if ((gradPhiY.size() > 0) &&
                (std::fabs(coeff2_AY * phi_0 * gradPhiY[0].value) > cutOffValue)) {
                std::vector<WFGrid>::const_iterator pEnd = std::upper_bound(phi.begin(), phi.end(),
                                                                            WFGrid(0, cutOffValue / (coeff2_AY * gradPhiY[0].value)),
                                                                            WFGrid_sort_functional());
                this->buildFock(phi.begin(), pEnd, gradPhiY.begin(), gradPhiY.end(),
                                coeff2_AY, cutOffValue, pF_A);
                this->buildFock(gradPhiY.begin(), gradPhiY.end(), phi.begin(), pEnd,
                                coeff2_AY, cutOffValue, pF_A);
            }
            // phi v.s. grad-phi(z)
            if ((gradPhiZ.size() > 0) &&
                (std::fabs(coeff2_AZ * phi_0 * gradPhiZ[0].value) > cutOffValue)) {
                std::vector<WFGrid>::const_iterator pEnd = std::upper_bound(phi.begin(), phi.end(),
                                                                            WFGrid(0, cutOffValue / (coeff2_AZ * gradPhiZ[0].value)),
                                                                            WFGrid_sort_functional());
                this->buildFock(phi.begin(), pEnd, gradPhiZ.begin(), gradPhiZ.end(),
                                coeff2_AZ, cutOffValue, pF_A);
                this->buildFock(gradPhiZ.begin(), gradPhiZ.end(), phi.begin(), pEnd,
                                coeff2_AZ, cutOffValue, pF_A);
            }
        }
    }
}

void DfCalcGridX::getWholeDensity(double* pRhoA, double* pRhoB) const
{
    TlMatrix gridMat;
    gridMat.load(this->getGridMatrixPath(this->m_nIteration));
    const index_type numOfGrids = gridMat.getNumOfRows();
    
    TlMatrix weightMat = gridMat.getBlockMatrix(0, GM_WEIGHT,
                                                numOfGrids, 1);

    const DfXCFunctional dfXcFunc(this->pPdfParam_);
    if (this->m_nMethodType == METHOD_RKS) {
        TlMatrix rhoMat_A;
        if (dfXcFunc.getFunctionalType() == DfXCFunctional::LDA) {
            rhoMat_A = gridMat.getBlockMatrix(0, GM_LDA_RHO_ALPHA,
                                              numOfGrids, 1);
        } else {
            rhoMat_A = gridMat.getBlockMatrix(0, GM_GGA_RHO_ALPHA,
                                              numOfGrids, 1);
        }

        assert(pRhoA != NULL);
        *pRhoA = weightMat.dot(rhoMat_A).sum();
    } else {
        TlMatrix rhoMat_A;
        TlMatrix rhoMat_B;
        if (dfXcFunc.getFunctionalType() == DfXCFunctional::LDA) {
            rhoMat_A = gridMat.getBlockMatrix(0, GM_LDA_RHO_ALPHA,
                                              numOfGrids, 1);
            rhoMat_B = gridMat.getBlockMatrix(0, GM_LDA_RHO_BETA,
                                              numOfGrids, 1);
        } else {
            rhoMat_A = gridMat.getBlockMatrix(0, GM_GGA_RHO_ALPHA,
                                              numOfGrids, 1);
            rhoMat_B = gridMat.getBlockMatrix(0, GM_GGA_RHO_BETA,
                                              numOfGrids, 1);
        }
        
        assert(pRhoA != NULL);
        assert(pRhoB != NULL);
        *pRhoA = weightMat.dot(rhoMat_A).sum();
        *pRhoB = weightMat.dot(rhoMat_B).sum();
    }
}

