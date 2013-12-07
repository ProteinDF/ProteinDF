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

#include "DfForce.h"
#include "Fl_Geometry.h"
#include "DfHpqX.h"
#include "DfOverlapX.h"
#include "DfEriX.h"
#include "DfCalcGridX.h"
#include "DfXCFunctional.h"

#include "DfFunctional_SVWN.h"
#include "DfFunctional_B88LYP.h"
#include "DfFunctional_B3LYP.h"
#include "TlMsgPack.h"
#include "TlTime.h"

DfForce::DfForce(TlSerializeData* pPdfParam)
    : DfObject(pPdfParam),
      orbitalInfo_((*pPdfParam)["coordinates"],
                   (*pPdfParam)["basis_sets"]),
      orbitalInfoDens_((*pPdfParam)["coordinates"],
                       (*pPdfParam)["basis_sets_j"]) {
    // initialize
    this->force_.resize(this->m_nNumOfAtoms, 3);

    this->storedCutoffValue_ = (*pPdfParam)["cut-value"].getDouble();
    if ((*pPdfParam)["force-cut-value"].getStr().empty() != true) {
        (*pPdfParam)["cut-value"] = (*pPdfParam)["force-cut-value"];
    }
    
    // debug
    this->isDebugOutMatrix_ = false;
    if ((*pPdfParam)["debug/save_forces"].getStr().empty() != true) {
        this->isDebugOutMatrix_ = (*pPdfParam)["debug/save_forces"].getBoolean();
    }
}


DfForce::~DfForce()
{
    (*(this->pPdfParam_))["cut-value"] = this->storedCutoffValue_;
}


void DfForce::calcForce()
{
    RUN_TYPE runType = RUN_RKS;
    const int iteration = this->m_nIteration;
    
    this->calcForceFromNuclei();
    this->calcForceFromWS(runType);
    
    TlSymmetricMatrix P = this->getPpqMatrix<TlSymmetricMatrix>(runType, iteration);
    this->loggerTime("core H");
    this->calcForceFromHpq(P);

    this->loggerTime("coulomb");
    this->calcForceFromCoulomb(runType);

    this->loggerTime("pureXC");
    this->calcForceFromPureXC(runType);
    
    this->loggerTime("Fock exchange");
    this->calcForceFromK(runType);

    this->force_ *= -1.0;
    this->saveForce();
}


void DfForce::saveForce()
{
    this->force_.save("force.mtx");
}


void DfForce::outputStartTitle(const std::string& stepName, const char lineChar)
{
    const std::string timeString = TlUtils::format("[%s %s]", TlTime::getNowDate().c_str(), TlTime::getNowTime().c_str());

    std::string title = ">>>> " + stepName + " ";
    TlUtils::pad(title, (72 - timeString.length()), ' ');
    title += timeString;

    // 出力
    this->logger(title + "\n");
}


void DfForce::outputEndTitle(const std::string& stepName, const char lineChar)
{
    const std::string timeString = TlUtils::format("[%s %s]", TlTime::getNowDate().c_str(), TlTime::getNowTime().c_str());

    std::string title = "<<<< " + stepName + " ";
    TlUtils::pad(title, (72 - timeString.length()), ' ');
    title += timeString;

    // 出力
    this->logger(title + "\n");
}


void DfForce::calcForceFromNuclei()
{
    const int numOfAtoms = this->m_nNumOfAtoms;

    const Fl_Geometry flGeom((*this->pPdfParam_)["coordinates"]);

    TlMatrix F_nuc(this->m_nNumOfAtoms, 3);
    for (int i = 0; i < numOfAtoms; ++i) {
        const TlPosition posI = flGeom.getCoordinate(i);
        const double chargeI = flGeom.getCharge(i);

        for (int j = i+1; j < numOfAtoms; ++j) {
            const TlPosition posJ = flGeom.getCoordinate(j);
            const double chargeJ = flGeom.getCharge(j);

            const TlPosition IJ = posJ - posI;
            const double distance = IJ.distanceFrom();
            const double invIJ3 = 1.0 / (distance * distance * distance);

            const double chargeIJ = chargeI * chargeJ;
            const double coef = - chargeIJ * invIJ3;
            const double fx = coef * IJ.x();
            const double fy = coef * IJ.y();
            const double fz = coef * IJ.z();

            F_nuc.add(i, X, fx);
            F_nuc.add(i, Y, fy);
            F_nuc.add(i, Z, fz);
            F_nuc.add(j, X, -fx);
            F_nuc.add(j, Y, -fy);
            F_nuc.add(j, Z, -fz);
        }
    }

    F_nuc *= -1.0;
    if (this->isDebugOutMatrix_ == true) {
        F_nuc.save("F_nuc.mtx");
    }

    this->force_ += F_nuc;
}


void DfForce::calcForceFromHpq(const TlSymmetricMatrix& P)
{
    // const index_type numOfAOs = this->m_nNumOfAOs;
    
    DfHpqX dfHpqX(this->pPdfParam_);
    TlMatrix force_Hpq(this->m_nNumOfAtoms, 3);
    dfHpqX.getForce(P, &force_Hpq);

    if (this->isDebugOutMatrix_ == true) {
        force_Hpq.save("F_h.mtx");
    }
    this->force_ += force_Hpq;
}


void DfForce::calcForceFromWS(RUN_TYPE runType)
{
    // const index_type numOfAOs = this->m_nNumOfAOs;
    
    DfOverlapX dfOvpX(this->pPdfParam_);

    const TlSymmetricMatrix W = this->getEnergyWeightedDensityMatrix(runType);
    if (this->isDebugOutMatrix_ == true) {
        W.save("W.mtx");
    }

    TlMatrix F_WS(this->m_nNumOfAtoms, 3);
    dfOvpX.getForce(W, &F_WS);

    F_WS *= -1.0;
    if (this->isDebugOutMatrix_ == true) {
        F_WS.save("F_WS.mtx");
    }

    this->force_ += F_WS;
}


TlMatrix DfForce::getEnergyWeightedDensityMatrix(RUN_TYPE runType)
{
    const int iteration = this->m_nIteration;
    const int numOfAOs = this->m_nNumOfAOs;
    const int numOfMOs = this->m_nNumOfMOs;

    TlVector eps;
    eps.load(this->getOccupationPath(runType));
    assert(eps.getSize() == numOfMOs);
    {
        TlVector eig;
        eig.load(this->getEigenvaluesPath(runType, iteration));
        assert(eig.getSize() == numOfMOs);

        // TODO: 高速化
        for (int i = 0; i < numOfMOs; ++i) {
            eps[i] *= eig[i];
        }
    }

    TlMatrix C = this->getCMatrix<TlMatrix>(runType, iteration);
    assert(C.getNumOfRows() == numOfAOs);
    assert(C.getNumOfCols() == numOfMOs);

    // TODO: 高速化
    TlMatrix W(numOfAOs, numOfAOs);
    for (int m = 0; m < numOfAOs; ++m) {
        for (int n = 0; n < numOfAOs; ++n) {
            double value = 0.0;
            for (int i = 0; i < numOfMOs; ++i) {
                value += eps[i] * C.get(m, i) * C.get(n, i);
            }
            W.set(m, n, value);
        }
    }

    return W;
}


void DfForce::calcForceFromCoulomb(RUN_TYPE runType)
{
    if (this->J_engine_ == J_ENGINE_RI_J) {
        this->calcForceFromCoulomb_RIJ(runType);
    } else {
        this->calcForceFromCoulomb_exact(runType);
    }
}


void DfForce::calcForceFromCoulomb_exact(RUN_TYPE runType)
{
    const int iteration = this->m_nIteration;
    // const int numOfAOs = this->m_nNumOfAOs;
    const int numOfAtoms = this->m_nNumOfAtoms;
    
    const TlSymmetricMatrix P = this->getPpqMatrix<TlSymmetricMatrix>(runType, iteration);

    DfEriX dfEri(this->pPdfParam_);

    // ((pq)'|(rs))
    TlMatrix F_J(numOfAtoms, 3);
    dfEri.getForceJ(P, &F_J);

    F_J *= 0.5;
    if (this->isDebugOutMatrix_ == true) {
        F_J.save("F_J.mtx");
    }

    this->force_ += F_J;
}


void DfForce::calcForceFromCoulomb_RIJ(const RUN_TYPE runType)
{
    const int iteration = this->m_nIteration;
    const int numOfAtoms = this->m_nNumOfAtoms;
    
    const TlVector rho = this->getRho<TlVector>(runType, iteration);
    const TlSymmetricMatrix P = this->getPpqMatrix<TlSymmetricMatrix>(runType, iteration);

    DfEriX dfEri(this->pPdfParam_);

    // ((pq)'|a)
    TlMatrix F_pqa(numOfAtoms, 3);
    dfEri.getForceJ(P, rho, &F_pqa);

    // (a'|b)
    TlMatrix F_ab(numOfAtoms, 3);
    dfEri.getForceJ(rho, &F_ab);

    if (this->isDebugOutMatrix_ == true) {
        F_pqa.save("F_pqa.mtx");
        F_ab.save("F_ab.mtx");
    }

    const TlMatrix F_J = (F_pqa - 0.5 * F_ab);
    this->force_ += F_J;
}


void DfForce::calcForceFromPureXC(RUN_TYPE runType)
{
    const int iteration = this->m_nIteration;
    const int numOfAOs = this->m_nNumOfAOs;
    const int numOfAtoms = this->m_nNumOfAtoms;
    const int numOfRealAtoms = this->numOfRealAtoms_;
    
    const TlSymmetricMatrix P = 0.5 * this->getPpqMatrix<TlSymmetricMatrix>(runType, iteration);

    DfCalcGridX calcGrid(this->pPdfParam_);

    TlMatrix Gx(numOfAOs, numOfAtoms);
    TlMatrix Gy(numOfAOs, numOfAtoms);
    TlMatrix Gz(numOfAOs, numOfAtoms);

    DfXCFunctional dfXCFunctional(this->pPdfParam_);
    DfXCFunctional::XC_TYPE xcType = dfXCFunctional.getXcType();
    switch (xcType) {
    case DfXCFunctional::SVWN:
        {
            DfFunctional_SVWN svwn;
            calcGrid.makeGammaMatrix(P, &svwn,
                                     &Gx, &Gy, &Gz);
        }
        break;

    case DfXCFunctional::BLYP:
        {
            DfFunctional_B88LYP blyp;
            calcGrid.makeGammaMatrix(P, &blyp,
                                     &Gx, &Gy, &Gz);
        }
        break;

    case DfXCFunctional::B3LYP:
        {
            DfFunctional_B3LYP b3lyp;
            calcGrid.makeGammaMatrix(P, &b3lyp,
                                     &Gx, &Gy, &Gz);
        }
        break;

    case DfXCFunctional::HF:
        {
            // do nothing
        }
        break;

    default:
        std::cerr << "unsupported functional. XC_TYPE=" << xcType << std::endl;
        abort();
        break;
    }
    if (this->isDebugOutMatrix_ == true) {
        Gx.save("Gx.mtx");
        Gy.save("Gy.mtx");
        Gz.save("Gz.mtx");
    }
    
    TlMatrix Fxc(numOfAtoms, 3);
    for (int mu = 0; mu < numOfRealAtoms; ++mu) {
        for (index_type p = 0; p < numOfAOs; ++p) {
            const index_type orbAtomId = this->orbitalInfo_.getAtomIndex(p);
            if (mu != orbAtomId) {
                const double fx = Gx.get(p, mu);
                const double fy = Gy.get(p, mu);
                const double fz = Gz.get(p, mu);

                Fxc.add(mu, X, fx);
                Fxc.add(mu, Y, fy);
                Fxc.add(mu, Z, fz);
            } else {
                for (int nu = 0; nu < numOfRealAtoms; ++nu) {
                    if (mu != nu) {
                        const double fx = Gx.get(p, nu);
                        const double fy = Gy.get(p, nu);
                        const double fz = Gz.get(p, nu);
                        Fxc.add(mu, X, -fx);
                        Fxc.add(mu, Y, -fy);
                        Fxc.add(mu, Z, -fz);
                    }
                }
            }
        }
    }

    Fxc *= -1.0;
    if (this->isDebugOutMatrix_ == true) {
        Fxc.save("F_xc.mtx");
    }
    this->force_ += Fxc;
}


void DfForce::calcForceFromK(RUN_TYPE runType)
{
    const DfXCFunctional dfXCFunctional(this->pPdfParam_);

    if (dfXCFunctional.isHybridFunctional() == true) {
        const int iteration = this->m_nIteration;
        // const int numOfAOs = this->m_nNumOfAOs;
        const int numOfAtoms = this->m_nNumOfAtoms;
        
        const TlSymmetricMatrix P = this->getPpqMatrix<TlSymmetricMatrix>(runType, iteration);
        
        DfEriX dfEri(this->pPdfParam_);
        
        TlMatrix F_K(numOfAtoms, 3);
        dfEri.getForceK(0.5 * P, &F_K);

        F_K *= -1.0;
        F_K *= dfXCFunctional.getFockExchangeCoefficient(); // for B3LYP

        if (this->isDebugOutMatrix_ == true) {
            F_K.save("F_K.mtx");
        }
        this->force_ += F_K;
    }
}
