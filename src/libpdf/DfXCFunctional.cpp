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

#include "DfXCFunctional.h"
#include <cassert>
#include "DfCalcGridX.h"
#include "DfEriX.h"

#include "TlTime.h"
#include "TlUtils.h"
#include "tl_matrix_object.h"

// #define CHECK_GRID_ACCURACY

// double DfXCFunctional::m_dXC_Energy = 0.0;
// double DfXCFunctional::m_dFockExchangeEnergyAlpha = 0.0; //
// Update法を使うため、直前のenergyを保存 double
// DfXCFunctional::m_dFockExchangeEnergyBeta = 0.0; //
// Update法を使うため、直前のenergyを保存

bool DfXCFunctional::isCalcd_E_disp_ = false;
double DfXCFunctional::E_disp_ = 0.0;

DfXCFunctional::DfXCFunctional(TlSerializeData* pPdfParam)
    : DfObject(pPdfParam), m_bIsHybrid(false), m_dFockExchangeCoef(1.0) {
    const TlSerializeData& pdfParam = *pPdfParam;

    // DfXCFunctional::m_dXC_Energy =
    // pdfParam["control"]["XC_Energy"].getDouble();
    this->XC_energy_ = pdfParam["control"]["XC_Energy"].getDouble();
    DfXCFunctional::E_disp_ = pdfParam["control"]["DFT_D_Energy"].getDouble();

    this->m_bRI_K =
        (TlUtils::toUpper(pdfParam["RI-K"].getStr()) == "YES") ? true : false;

    this->m_bUseRTmethod = pdfParam["RT_method"].getBoolean();
    this->isSaveFxcPure_ = pdfParam["debug/save_Fxc_pure"].getBoolean();

    this->m_bDebugTEI =
        pdfParam["debug_two_electron_integral_object"].getBoolean();
    this->m_bPrintTEI =
        pdfParam["print_two_electron_integral_value"].getBoolean();

    this->isUseNewEngine_ = pdfParam["new_engine"].getBoolean();

    // XC function
    std::string checkXC = this->m_sXCFunctional;
    if (this->m_bIsXCFitting == true) {
        checkXC = checkXC.substr(0, checkXC.size() - 1);  // remove "~"
    }
    if (this->enableGrimmeDispersion_ == true) {
        checkXC = checkXC.substr(0, checkXC.size() - 2);  // remove "-D"
    }
    if (checkXC == "HFS") {
        // SHF
        this->m_nXCFunctional = HFS;
        this->functionalType_ = LDA;
        this->m_bIsHybrid = false;
    } else if (checkXC == "SVWN") {
        // SVWN
        this->m_nXCFunctional = SVWN;
        this->functionalType_ = LDA;
        this->m_bIsHybrid = false;
    } else if (checkXC == "HFB") {
        // HFB
        this->m_nXCFunctional = HFB;
        this->functionalType_ = GGA;
        this->m_bIsHybrid = false;
    } else if ((checkXC == "BLYP") || (checkXC == "B88LYP")) {
        // BLYP
        this->m_nXCFunctional = BLYP;
        this->functionalType_ = GGA;
        this->m_bIsHybrid = false;
    } else if (checkXC == "B3LYP") {
        // B3LYP
        this->m_nXCFunctional = B3LYP;
        this->functionalType_ = GGA;
        this->m_bIsHybrid = true;
        this->m_dFockExchangeCoef = 0.2;
    } else if (checkXC == "HF") {
        // RI-HF
        this->functionalType_ = UNDEFINED_TYPE;
        this->m_nXCFunctional = HF;
        this->m_bIsHybrid = true;
        this->m_dFockExchangeCoef = 1.0;
    } else {
        this->m_nXCFunctional = UNDEFINED_XC;
    }
}

DfXCFunctional::~DfXCFunctional() {
    (*this->pPdfParam_)["control"]["XC_Energy"] = this->XC_energy_;
    (*this->pPdfParam_)["control"]["DFT_D_Energy"] = DfXCFunctional::E_disp_;
}

bool DfXCFunctional::isHybridFunctional() const { return this->m_bIsHybrid; }

double DfXCFunctional::getFockExchangeCoefficient() const {
    return this->m_dFockExchangeCoef;
}

DfXCFunctional::FUNCTIONAL_TYPE DfXCFunctional::getFunctionalType() const {
    return this->functionalType_;
}

DfXCFunctional::XC_TYPE DfXCFunctional::getXcType() const {
    return this->m_nXCFunctional;
}

void DfXCFunctional::buildXcMatrix() {
    DfCalcGridX dfCalcGrid(this->pPdfParam_);

    switch (this->m_nMethodType) {
        case METHOD_RKS: {
            TlDenseSymmetricMatrix_Lapack Ppq;
            if (this->m_bIsUpdateXC == true) {
                Ppq = 0.5 *
                      this->getDiffDensityMatrix<TlDenseSymmetricMatrix_Lapack>(
                          RUN_RKS, this->m_nIteration);
            } else {
                Ppq = 0.5 * this->getPpqMatrix<TlDenseSymmetricMatrix_Lapack>(
                                RUN_RKS, this->m_nIteration - 1);
            }

            TlDenseSymmetricMatrix_Lapack Fxc(this->m_nNumOfAOs);
            this->loggerTime(" start: pure XC term");
            this->getFxc(Ppq, &dfCalcGrid, &Fxc);
            this->loggerTime(" end: pure XC term");

            if (this->isSaveFxcPure_ == true) {
                this->saveFxcPureMatrix(RUN_RKS, this->m_nIteration, Fxc);
            }

            this->saveFxcMatrix<TlDenseSymmetricMatrix_Lapack>(
                RUN_RKS, this->m_nIteration, Fxc);
        } break;

        case METHOD_UKS: {
            TlDenseSymmetricMatrix_Lapack PApq;
            TlDenseSymmetricMatrix_Lapack PBpq;
            if (this->m_bIsUpdateXC == true) {
                PApq =
                    this->getDiffDensityMatrix<TlDenseSymmetricMatrix_Lapack>(
                        RUN_UKS_ALPHA, this->m_nIteration);
                PBpq =
                    this->getDiffDensityMatrix<TlDenseSymmetricMatrix_Lapack>(
                        RUN_UKS_BETA, this->m_nIteration);
            } else {
                PApq = this->getPpqMatrix<TlDenseSymmetricMatrix_Lapack>(
                    RUN_UKS_ALPHA, this->m_nIteration - 1);
                PBpq = this->getPpqMatrix<TlDenseSymmetricMatrix_Lapack>(
                    RUN_UKS_BETA, this->m_nIteration - 1);
            }

            TlDenseSymmetricMatrix_Lapack FxcA(this->m_nNumOfAOs);
            TlDenseSymmetricMatrix_Lapack FxcB(this->m_nNumOfAOs);
            this->loggerTime(" start: pure XC term");
            this->getFxc(PApq, PBpq, &dfCalcGrid, &FxcA, &FxcB);
            this->loggerTime(" end: pure XC term");
            if (this->isSaveFxcPure_ == true) {
                this->saveFxcPureMatrix(RUN_UKS_ALPHA, this->m_nIteration,
                                        FxcA);
                this->saveFxcPureMatrix(RUN_UKS_BETA, this->m_nIteration, FxcB);
            }

            this->saveFxcMatrix<TlDenseSymmetricMatrix_Lapack>(
                RUN_UKS_ALPHA, this->m_nIteration, FxcA);
            this->saveFxcMatrix<TlDenseSymmetricMatrix_Lapack>(
                RUN_UKS_BETA, this->m_nIteration, FxcB);
        } break;

        case METHOD_ROKS: {
            TlDenseSymmetricMatrix_Lapack PApq;
            TlDenseSymmetricMatrix_Lapack PBpq;
            if (this->m_bIsUpdateXC == true) {
                PApq =
                    0.5 *
                    this->getDiffDensityMatrix<TlDenseSymmetricMatrix_Lapack>(
                        RUN_ROKS_CLOSED, this->m_nIteration);
                PBpq = PApq;
                PApq +=
                    this->getDiffDensityMatrix<TlDenseSymmetricMatrix_Lapack>(
                        RUN_ROKS_OPEN, this->m_nIteration);
            } else {
                PApq = 0.5 * this->getPpqMatrix<TlDenseSymmetricMatrix_Lapack>(
                                 RUN_ROKS_CLOSED, this->m_nIteration - 1);
                PBpq = PApq;
                PApq += this->getPpqMatrix<TlDenseSymmetricMatrix_Lapack>(
                    RUN_ROKS_OPEN, this->m_nIteration - 1);
            }

            TlDenseSymmetricMatrix_Lapack FxcA(this->m_nNumOfAOs);
            TlDenseSymmetricMatrix_Lapack FxcB(this->m_nNumOfAOs);
            this->loggerTime(" start: pure XC term");
            this->getFxc(PApq, PBpq, &dfCalcGrid, &FxcA, &FxcB);
            this->loggerTime(" end: pure XC term");

            if (this->isSaveFxcPure_ == true) {
                this->saveFxcPureMatrix(RUN_ROKS_ALPHA, this->m_nIteration,
                                        FxcA);
                this->saveFxcPureMatrix(RUN_ROKS_BETA, this->m_nIteration,
                                        FxcB);
            }
            this->saveFxcMatrix<TlDenseSymmetricMatrix_Lapack>(
                RUN_ROKS_ALPHA, this->m_nIteration, FxcA);
            this->saveFxcMatrix<TlDenseSymmetricMatrix_Lapack>(
                RUN_ROKS_BETA, this->m_nIteration, FxcB);
        } break;

        default:
            break;
    }
}

void DfXCFunctional::checkGridAccuracy() {
    if (this->m_nXCFunctional != HF) {
        DfCalcGridX dfCalcGrid(this->pPdfParam_);

        this->log_.info(" grid accuracy check:");
        double rhoA = 0.0;
        double rhoB = 0.0;
        dfCalcGrid.getWholeDensity(&rhoA, &rhoB);
        if (this->m_nMethodType == METHOD_RKS) {
            this->log_.info(
                TlUtils::format(" number of electrons = % 16.10f (input: %d)",
                                rhoA * 2, this->m_nNumOfElectrons));
        } else {
            this->log_.info(TlUtils::format(
                " number of electrons(alpha) = % 16.10f (input: %d)", rhoA,
                this->m_nNumOfAlphaElectrons));
            this->log_.info(TlUtils::format(
                " number of electrons(beta ) = % 16.10f (input: %d)", rhoB,
                this->m_nNumOfBetaElectrons));
        }
    }
}

DfEriX* DfXCFunctional::getDfEriXObject() {
    DfEriX* pDfEriX = new DfEriX(this->pPdfParam_);

    return pDfEriX;
}

double DfXCFunctional::getGrimmeDispersionEnergy() {
    if (DfXCFunctional::isCalcd_E_disp_ != true) {
        // C6 parameters in Grimme2004.
        // the unit is 'J nm^6 mol^-1'
        // const double C6Coef = (1.0 / 4.35974417E-18) * std::pow((1.0E-9
        // / 5.291772108E-11), 6.0);
        const double C6Coef = 1.0;
        std::map<std::string, double> C6;
        C6["H"] = C6Coef * 0.14;
        C6["He"] = C6Coef * 0.08;
        C6["Li"] = C6Coef * 1.61;
        C6["Be"] = C6Coef * 1.61;
        C6["B"] = C6Coef * 3.13;
        C6["C"] = C6Coef * 1.75;
        C6["N"] = C6Coef * 1.23;
        C6["O"] = C6Coef * 0.70;
        C6["F"] = C6Coef * 0.75;
        C6["Ne"] = C6Coef * 0.63;
        C6["Na"] = C6Coef * 5.71;
        C6["Mg"] = C6Coef * 5.71;
        C6["Al"] = C6Coef * 10.79;
        C6["Si"] = C6Coef * 9.23;
        C6["P"] = C6Coef * 7.84;
        C6["S"] = C6Coef * 5.57;
        C6["Cl"] = C6Coef * 5.07;
        C6["Ar"] = C6Coef * 4.61;
        C6["K"] = C6Coef * 10.80;
        C6["Ca"] = C6Coef * 10.80;
        C6["Sc"] = C6Coef * 10.80;  // to Zn
        C6["Ga"] = C6Coef * 16.99;
        C6["Ge"] = C6Coef * 17.10;
        C6["As"] = C6Coef * 16.37;
        C6["Se"] = C6Coef * 12.64;
        C6["Br"] = C6Coef * 12.47;
        C6["Kr"] = C6Coef * 12.01;
        C6["Br"] = C6Coef * 12.47;
        C6["Kr"] = C6Coef * 12.01;
        C6["Rb"] = C6Coef * 24.67;
        C6["Sr"] = C6Coef * 24.67;
        C6["Y"] = C6Coef * 24.67;  // to Cd
        C6["In"] = C6Coef * 37.32;
        C6["Sn"] = C6Coef * 38.71;
        C6["Sb"] = C6Coef * 38.44;
        C6["Te"] = C6Coef * 31.74;
        C6["I"] = C6Coef * 31.50;
        C6["Xe"] = C6Coef * 29.99;

        // R0(vdw) parameters in Grimme2004.
        // the unit is 'pm'
        // convert coefficient = 1.0E-10 / 5.291772108E-11
        // const double R0Coef = 1.0E-10 / 5.291772108E-11;
        const double R0Coef = 1.0;
        std::map<std::string, double> R0;
        R0["H"] = R0Coef * 1.001;
        R0["He"] = R0Coef * 1.001;
        R0["Li"] = R0Coef * 0.825;
        R0["Li"] = R0Coef * 0.825;
        R0["Be"] = R0Coef * 1.408;
        R0["B"] = R0Coef * 1.485;
        R0["C"] = R0Coef * 1.452;
        R0["N"] = R0Coef * 1.397;
        R0["O"] = R0Coef * 1.342;
        R0["F"] = R0Coef * 1.287;
        R0["Ne"] = R0Coef * 1.273;
        R0["Na"] = R0Coef * 1.144;
        R0["Mg"] = R0Coef * 1.364;
        R0["Al"] = R0Coef * 1.639;
        R0["Si"] = R0Coef * 1.716;
        R0["P"] = R0Coef * 1.705;
        R0["S"] = R0Coef * 1.687;
        R0["Cl"] = R0Coef * 1.639;
        R0["Ar"] = R0Coef * 1.595;
        R0["K"] = R0Coef * 1.485;
        R0["Ca"] = R0Coef * 1.474;
        R0["Sc"] = R0Coef * 1.562;  // to Zn
        R0["Ga"] = R0Coef * 1.650;
        R0["Ge"] = R0Coef * 1.727;
        R0["As"] = R0Coef * 1.760;
        R0["Se"] = R0Coef * 1.771;
        R0["Br"] = R0Coef * 1.749;
        R0["Kr"] = R0Coef * 1.727;
        R0["Rb"] = R0Coef * 1.628;
        R0["Sr"] = R0Coef * 1.606;
        R0["Y"] = R0Coef * 1.639;  // to Cd
        R0["In"] = R0Coef * 1.672;
        R0["Sn"] = R0Coef * 1.804;
        R0["Sb"] = R0Coef * 1.881;
        R0["Te"] = R0Coef * 1.892;
        R0["I"] = R0Coef * 1.892;
        R0["Xe"] = R0Coef * 1.881;

        // alpha value
        const double alpha = 20.0;  // Grimme2006

        const Fl_Geometry geom((*(this->pPdfParam_))["coordinates"]);
        const int numOfAtoms = geom.getNumOfAtoms();
        double E_disp = 0.0;
        for (int i = 0; i < numOfAtoms; ++i) {
            const TlPosition posA = geom.getCoordinate(i);
            const std::string symbolA = geom.getAtomSymbol(i);
            if ((C6.find(symbolA) == C6.end()) ||
                (R0.find(symbolA) == R0.end())) {
                this->logger(
                    TlUtils::format(" Grimme dispersion parameters in atom(%s) "
                                    "is not regiestered.\n",
                                    symbolA.c_str()));
                continue;
            }
            const double C6i = C6[symbolA];
            const double R0i = R0[symbolA];

            for (int j = i + 1; j < numOfAtoms; ++j) {
                const TlPosition posB = geom.getCoordinate(j);
                const std::string symbolB = geom.getAtomSymbol(j);
                if ((C6.find(symbolA) == C6.end()) ||
                    (R0.find(symbolA) == R0.end())) {
                    continue;
                }
                const double C6j = C6[symbolB];
                const double R0j = R0[symbolB];

                const double Rij2 = posA.squareDistanceFrom(posB);
                const double Rij = std::sqrt(Rij2);
                const double Rij6 = Rij2 * Rij2 * Rij2;
                const double R0ij = R0i + R0j;
                const double f_dmp =
                    1.0 / (1.0 + std::exp(-alpha * (Rij / R0ij - 1.0)));
                const double C6ij = 2.0 * C6i * C6j / (C6i + C6j);
                E_disp += C6ij / Rij6 * f_dmp;
            }
        }

        // s6 depends on XC functional.
        double s6 = 1.0;
        switch (this->m_nXCFunctional) {
            case BLYP:
                s6 = 1.2;
                break;

            case B3LYP:
                s6 = 1.05;
                break;

                //     case BP86:
                //         s6 = 1.05;
                //         break;

                //     case PBE:
                //         s6 = 0.75;
                //         break;

            default:
                s6 = 1.0;
                this->logger(
                    " The Grimme s6 parameter is not registerd in the XC "
                    "functional.\n");
                this->logger(TlUtils::format(" use s6 = %f\n", s6));
                break;
        }
        E_disp *= (-s6);

        DfXCFunctional::E_disp_ = E_disp;
        DfXCFunctional::isCalcd_E_disp_ = true;

        // std::cerr << "E_disp = " << E_disp << std::endl;
    }

    return DfXCFunctional::E_disp_;
}
