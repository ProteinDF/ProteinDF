#include <cassert>
#include "DfXCFunctional.h"
#include "DfEri2.h"
#include "DfCalcGridX.h"
#include "DfTwoElectronIntegral.h"
#include "DfEriX.h"

#include "FileX.h"
#include "TlMatrixObject.h"
#include "TlUtils.h"
#include "TlTime.h"

double DfXCFunctional::m_dXC_Energy = 0.0;
double DfXCFunctional::m_dFockExchangeEnergyAlpha = 0.0; // Update法を使うため、直前のenergyを保存
double DfXCFunctional::m_dFockExchangeEnergyBeta = 0.0; // Update法を使うため、直前のenergyを保存

bool DfXCFunctional::isCalcd_E_disp_ = false;
double DfXCFunctional::E_disp_ = 0.0;

DfXCFunctional::DfXCFunctional(TlSerializeData* pPdfParam)
    : DfObject(pPdfParam), m_bIsHybrid(false), m_dFockExchangeCoef(1.0)
{
    const TlSerializeData& pdfParam = *pPdfParam;

    DfXCFunctional::m_dXC_Energy = pdfParam["control"]["XC_Energy"].getDouble();
    DfXCFunctional::E_disp_ = pdfParam["control"]["DFT_D_Energy"].getDouble();

    this->m_bRI_K = (TlUtils::toUpper(pdfParam["RI-K"].getStr()) == "YES") ? true : false;

    this->m_bUseRTmethod = pdfParam["RT_method"].getBoolean();
    this->isSaveFxcPure_ = pdfParam["save_Fxc_pure"].getBoolean();

    this->m_bDebugTEI = pdfParam["debug_two_electron_integral_object"].getBoolean();
    this->m_bPrintTEI = pdfParam["print_two_electron_integral_value"].getBoolean();
    this->m_bKMatrixDebugOut = pdfParam["debug-K-matrix"].getBoolean();
    this->isDebugOutFockExchangeMatrix_ = pdfParam["debugout_fock_exchange"].getBoolean();

    this->isUseNewEngine_ = pdfParam["new_engine"].getBoolean();
    
    // XC function
    std::string checkXC = this->m_sXCFunctional;
    if (this->enableGrimmeDispersion_ == true) {
        checkXC = checkXC.substr(0, checkXC.size() -2);
    }
    if (checkXC == "SHF") {
        // SHF
        this->m_nXCFunctional = SHF;
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
        this->functionalType_ = LDA;
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


DfXCFunctional::~DfXCFunctional()
{
    (*this->pPdfParam_)["control"]["XC_Energy"] = DfXCFunctional::m_dXC_Energy;
    (*this->pPdfParam_)["control"]["DFT_D_Energy"] = DfXCFunctional::E_disp_;
}


bool DfXCFunctional::isHybridFunctional() const
{
    return this->m_bIsHybrid;
}


double DfXCFunctional::getFockExchangeCoefficient() const
{
    return this->m_dFockExchangeCoef;
}


DfXCFunctional::FUNCTIONAL_TYPE DfXCFunctional::getFunctionalType() const
{
    return this->functionalType_;
}


DfXCFunctional::XC_TYPE DfXCFunctional::getXcType() const
{
    return this->m_nXCFunctional;
}


void DfXCFunctional::buildXcMatrix()
{
    DfCalcGridX dfCalcGrid(this->pPdfParam_);

    switch (this->m_nMethodType) {
    case METHOD_RKS: {
        // 密度行列の準備
        TlSymmetricMatrix Ppq;
        if (this->m_bIsUpdateXC == true) {
            Ppq = this->getDiffDensityMatrix<TlSymmetricMatrix>(RUN_RKS, this->m_nIteration);
        } else {
            Ppq = this->getPpqMatrix<TlSymmetricMatrix>(RUN_RKS, this->m_nIteration -1);
        }

        TlSymmetricMatrix Fxc(this->m_nNumOfAOs);
        this->loggerTime(" start: pure XC term");
        this->getFxc(Ppq, &dfCalcGrid, &Fxc);
        this->loggerTime(" end: pure XC term");

        if (this->isSaveFxcPure_ == true) {
            this->saveFxcPureMatrix(RUN_RKS, this->m_nIteration, Fxc);
        }
        
        if (this->m_bIsHybrid == true) {
            // Fockの交換項を求める
            this->loggerTime(" start: Fock exchange");
            Fxc += this->getFockExchange((0.5 * this->m_dFockExchangeCoef) * Ppq, RUN_RKS);
            this->m_dXC_Energy += DfXCFunctional::m_dFockExchangeEnergyAlpha;
            this->loggerTime(" end: Fock exchange");
        }

        this->saveFxcMatrix<TlSymmetricMatrix>(RUN_RKS, this->m_nIteration, Fxc);
    }
    break;

    case METHOD_UKS: {
        TlSymmetricMatrix PApq;
        TlSymmetricMatrix PBpq;
        if (this->m_bIsUpdateXC == true) {
            PApq = this->getDiffDensityMatrix<TlSymmetricMatrix>(RUN_UKS_ALPHA, this->m_nIteration);
            PBpq = this->getDiffDensityMatrix<TlSymmetricMatrix>(RUN_UKS_BETA, this->m_nIteration);
        } else {
            PApq = this->getPpqMatrix<TlSymmetricMatrix>(RUN_UKS_ALPHA, this->m_nIteration -1);
            PBpq = this->getPpqMatrix<TlSymmetricMatrix>(RUN_UKS_BETA, this->m_nIteration -1);
        }
        
        TlSymmetricMatrix FxcA(this->m_nNumOfAOs);
        TlSymmetricMatrix FxcB(this->m_nNumOfAOs);
        this->loggerTime(" start: pure XC term");
        this->getFxc(PApq, PBpq, &dfCalcGrid, &FxcA, &FxcB);
        this->loggerTime(" end: pure XC term");
        if (this->isSaveFxcPure_ == true) {
            this->saveFxcPureMatrix(RUN_UKS_ALPHA, this->m_nIteration, FxcA);
            this->saveFxcPureMatrix(RUN_UKS_BETA,  this->m_nIteration, FxcB);
        }

        if (this->m_bIsHybrid == true) {
            this->loggerTime(" start: Fock exchange");
            FxcA += this->getFockExchange(this->m_dFockExchangeCoef * PApq, RUN_UKS_ALPHA);
            FxcB += this->getFockExchange(this->m_dFockExchangeCoef * PBpq, RUN_UKS_BETA);
            this->m_dXC_Energy += (DfXCFunctional::m_dFockExchangeEnergyAlpha + DfXCFunctional::m_dFockExchangeEnergyBeta);
            this->loggerTime(" end: Fock exchange");
        }

        this->saveFxcMatrix<TlSymmetricMatrix>(RUN_UKS_ALPHA, this->m_nIteration, FxcA);
        this->saveFxcMatrix<TlSymmetricMatrix>(RUN_UKS_BETA, this->m_nIteration, FxcB);
    }
    break;

    case METHOD_ROKS:
        CnErr.abort("sorry. ROKS method is not implemented. @ DfXCFunctional::buildXcMatrix()");
        break;

    default:
        break;
    }
}


TlSymmetricMatrix DfXCFunctional::getFockExchange(const TlSymmetricMatrix& deltaP, const RUN_TYPE nRunType)
{
    this->loggerTime(" start to build fock exchange.");

    TlSymmetricMatrix FxHF = this->getFockExchange(deltaP);

    if (this->m_bIsUpdateXC == true) {
        if (this->m_nIteration >= 2) {
            FxHF += this->getHFxMatrix<TlSymmetricMatrix>(nRunType, this->m_nIteration -1);
        }
        this->saveHFxMatrix(nRunType, this->m_nIteration, FxHF);
    }

    // calc energy
    {
        TlSymmetricMatrix P = this->getPpqMatrix<TlSymmetricMatrix>(nRunType, this->m_nIteration -1);

        if (nRunType == RUN_RKS) {
            DfXCFunctional::m_dFockExchangeEnergyAlpha = 0.5 * this->getFockExchangeEnergy(P, FxHF);
        } else if (nRunType == RUN_UKS_ALPHA) {
            DfXCFunctional::m_dFockExchangeEnergyAlpha = 0.5 * this->getFockExchangeEnergy(P, FxHF);
        } else if (nRunType == RUN_UKS_BETA) {
            DfXCFunctional::m_dFockExchangeEnergyBeta = 0.5 * this->getFockExchangeEnergy(P, FxHF);
        }
    }

    this->loggerTime(" end to build fock exchange.");
    return FxHF;
}


TlSymmetricMatrix DfXCFunctional::getFockExchange(const TlSymmetricMatrix& P)
{
    TlSymmetricMatrix FxHF(this->m_nNumOfAOs);

    if (this->m_bRI_K == false) {
        // disable RI-K
        if (this->isUseNewEngine_ == true) {
            this->logger(" new ERI engine\n");
            DfEriX* pDfEri = this->getDfEriXObject();

            pDfEri->getK(P, &FxHF);

            delete pDfEri;
            pDfEri = NULL;
        } else {
            DfTwoElectronIntegral* pDfTei = this->getDfTwoElectronIntegral();
            
            if (this->m_bUseRTmethod == true) {
                this->logger(" using RT method\n");
                pDfTei->getContractKMatrixByRTmethod(P, &FxHF);
            } else {
                this->logger(" using integral-driven\n");
                pDfTei->getContractKMatrixByIntegralDriven(P, &FxHF);
            }
            
            delete pDfTei;
            pDfTei = NULL;
        }
    } else {
        // using RI-K
        DfEri2* pDfEri2 = this->getDfEri2();
        FxHF = pDfEri2->getKMatrix(P);

        delete pDfEri2;
        pDfEri2 = NULL;
    }

    return FxHF;
}


double DfXCFunctional::getFockExchangeEnergy(const TlSymmetricMatrix& P, const TlSymmetricMatrix& Ex)
{
    double dEnergy = 0.0;

    const int nNumOfAOs = this->m_nNumOfAOs;
    for (int i = 0; i < nNumOfAOs; ++i) {
        // case: i != j
        for (int j = 0; j < i; ++j) {
            dEnergy += 2.0 * P(i, j) * Ex(i, j);
        }
        // case: i == j
        dEnergy += P(i, i) * Ex(i, i);
    }

    //dEnergy *= 0.5;

    return dEnergy;
}


DfEri2* DfXCFunctional::getDfEri2()
{
    DfEri2* pDfEri2 = new DfEri2(this->pPdfParam_);

    return pDfEri2;
}


DfEriX* DfXCFunctional::getDfEriXObject()
{
    DfEriX* pDfEriX = new DfEriX(this->pPdfParam_);

    return pDfEriX;
}


DfTwoElectronIntegral* DfXCFunctional::getDfTwoElectronIntegral()
{
    DfTwoElectronIntegral* pDfTei = new DfTwoElectronIntegral(this->pPdfParam_);

    return pDfTei;
}


double DfXCFunctional::getGrimmeDispersionEnergy()
{
    if (DfXCFunctional::isCalcd_E_disp_ != true) {
        // C6 parameters in Grimme2004.
        // the unit is 'J nm^6 mol^-1'
        //const double C6Coef = (1.0 / 4.35974417E-18) * std::pow((1.0E-9 / 5.291772108E-11), 6.0);
        const double C6Coef = 1.0;
        std::map<std::string, double> C6;
        C6["H"]  = C6Coef *  0.14;
        C6["He"] = C6Coef *  0.08;
        C6["Li"] = C6Coef *  1.61;
        C6["Be"] = C6Coef *  1.61;
        C6["B"]  = C6Coef *  3.13;
        C6["C"]  = C6Coef *  1.75;
        C6["N"]  = C6Coef *  1.23;
        C6["O"]  = C6Coef *  0.70;
        C6["F"]  = C6Coef *  0.75;
        C6["Ne"] = C6Coef *  0.63;
        C6["Na"] = C6Coef *  5.71;
        C6["Mg"] = C6Coef *  5.71;
        C6["Al"] = C6Coef * 10.79;
        C6["Si"] = C6Coef *  9.23;
        C6["P"]  = C6Coef *  7.84;
        C6["S"]  = C6Coef *  5.57;
        C6["Cl"] = C6Coef *  5.07;
        C6["Ar"] = C6Coef *  4.61;
        C6["K"]  = C6Coef * 10.80;
        C6["Ca"] = C6Coef * 10.80;
        C6["Sc"] = C6Coef * 10.80; // to Zn
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
        C6["Y"]  = C6Coef * 24.67; // to Cd
        C6["In"] = C6Coef * 37.32;
        C6["Sn"] = C6Coef * 38.71;
        C6["Sb"] = C6Coef * 38.44;
        C6["Te"] = C6Coef * 31.74;
        C6["I"]  = C6Coef * 31.50;
        C6["Xe"] = C6Coef * 29.99;
        
        // R0(vdw) parameters in Grimme2004.
        // the unit is 'pm'
        // convert coefficient = 1.0E-10 / 5.291772108E-11
        //const double R0Coef = 1.0E-10 / 5.291772108E-11;
        const double R0Coef = 1.0;
        std::map<std::string, double> R0;
        R0["H"]  = R0Coef * 1.001;
        R0["He"] = R0Coef * 1.001;
        R0["Li"] = R0Coef * 0.825;
        R0["Li"] = R0Coef * 0.825;
        R0["Be"] = R0Coef * 1.408;
        R0["B"]  = R0Coef * 1.485;
        R0["C"]  = R0Coef * 1.452;
        R0["N"]  = R0Coef * 1.397;
        R0["O"]  = R0Coef * 1.342;
        R0["F"]  = R0Coef * 1.287;
        R0["Ne"] = R0Coef * 1.273;
        R0["Na"] = R0Coef * 1.144;
        R0["Mg"] = R0Coef * 1.364;
        R0["Al"] = R0Coef * 1.639;
        R0["Si"] = R0Coef * 1.716;
        R0["P"]  = R0Coef * 1.705;
        R0["S"]  = R0Coef * 1.687;
        R0["Cl"] = R0Coef * 1.639;
        R0["Ar"] = R0Coef * 1.595;
        R0["K"]  = R0Coef * 1.485;
        R0["Ca"] = R0Coef * 1.474;
        R0["Sc"] = R0Coef * 1.562; // to Zn
        R0["Ga"] = R0Coef * 1.650;
        R0["Ge"] = R0Coef * 1.727;
        R0["As"] = R0Coef * 1.760;
        R0["Se"] = R0Coef * 1.771;
        R0["Br"] = R0Coef * 1.749;
        R0["Kr"] = R0Coef * 1.727;
        R0["Rb"] = R0Coef * 1.628;
        R0["Sr"] = R0Coef * 1.606;
        R0["Y"]  = R0Coef * 1.639; // to Cd
        R0["In"] = R0Coef * 1.672;
        R0["Sn"] = R0Coef * 1.804;
        R0["Sb"] = R0Coef * 1.881;
        R0["Te"] = R0Coef * 1.892;
        R0["I"]  = R0Coef * 1.892;
        R0["Xe"] = R0Coef * 1.881;
        
        // alpha value
        const double alpha = 20.0; // Grimme2006
        
        const Fl_Geometry geom((*(this->pPdfParam_))["coordinates"]);
        const int numOfAtoms = geom.getNumOfAtoms();
        double E_disp = 0.0;
        for (int i = 0; i < numOfAtoms; ++i) {
            const TlPosition posA = geom.getCoordinate(i);
            const std::string symbolA = geom.getAtom(i);
            if ((C6.find(symbolA) == C6.end()) ||
                (R0.find(symbolA) == R0.end())) {
                this->logger(TlUtils::format(" Grimme dispersion parameters in atom(%s) is not regiestered.\n",
                                             symbolA.c_str()));
                continue;
            }
            const double C6i = C6[symbolA];
            const double R0i = R0[symbolA];
            
            for (int j = i + 1; j < numOfAtoms; ++j) {
                const TlPosition posB = geom.getCoordinate(j);
                const std::string symbolB = geom.getAtom(j);
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
                const double f_dmp = 1.0 / (1.0 + std::exp(- alpha * (Rij/R0ij - 1.0)));
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
            this->logger(" The Grimme s6 parameter is not registerd in the XC functional.\n");
            this->logger(TlUtils::format(" use s6 = %f\n", s6));
            break;
        }
        E_disp *= (- s6);
        
        DfXCFunctional::E_disp_ = E_disp;
        DfXCFunctional::isCalcd_E_disp_ = true;

        //std::cerr << "E_disp = " << E_disp << std::endl;
    }
    
    return DfXCFunctional::E_disp_;
}
