#include "DfXCFunctional_Parallel.h"
#include "DfCalcGridX_Parallel.h"
#include "DfEri2_Parallel.h"
#include "DfEriX_Parallel.h"
#include "DfTwoElectronIntegral_Parallel.h"
#include "TlCommunicate.h"
#include "TlTime.h"
#include "CnError.h"

#include "TlMemManager.h"

DfXCFunctional_Parallel::DfXCFunctional_Parallel(TlSerializeData* pPdfParam)
    : DfXCFunctional(pPdfParam)
{
}


DfXCFunctional_Parallel::~DfXCFunctional_Parallel()
{
//     std::cerr << "DfXCFunctional_Parallel::~DfXCFunctional_Parallel()" << std::endl;
}


void DfXCFunctional_Parallel::logger(const std::string& str) const
{
    TlCommunicate& rComm = TlCommunicate::getInstance();

    if (rComm.isMaster() == true) {
        DfXCFunctional::logger(str);
    }
}


void DfXCFunctional_Parallel::buildXcMatrix()
{
    const std::size_t needMem = this->m_nNumOfAOs * (this->m_nNumOfAOs + 1) * sizeof(double);
    if ((isWorkOnDisk_ == true) ||
        (this->procMaxMemSize_ < needMem)) {
        this->logger(" build XC matrix on disk.\n");
        TlMatrix::useMemManager(true);
    } else {
        this->logger(" build XC matrix on memory.\n");
        TlMatrix::useMemManager(false);
    }

    if (this->m_bUsingSCALAPACK == true) {
        this->buildXC_ScaLAPACK();
    } else {
        this->buildXC_LAPACK();
    }
}


void DfXCFunctional_Parallel::buildXC_LAPACK()
{
    TlCommunicate& rComm = TlCommunicate::getInstance();
    DfCalcGridX_Parallel dfCalcGrid(this->pPdfParam_);

    switch (this->m_nMethodType) {
    case METHOD_RKS: {
        TlSymmetricMatrix Ppq;
        if (rComm.isMaster() == true) {
            if (this->m_bIsUpdateXC == true) {
                Ppq = 0.5 * this->getDiffDensityMatrix<TlSymmetricMatrix>(RUN_RKS, this->m_nIteration);
            } else {
                Ppq = 0.5 * this->getPpqMatrix<TlSymmetricMatrix>(RUN_RKS, this->m_nIteration -1);
            }
        }
        rComm.broadcast(Ppq);

        TlSymmetricMatrix Fxc(this->m_nNumOfAOs);
        DfXCFunctional::getFxc(Ppq, &dfCalcGrid, &Fxc);
        if (this->isSaveFxcPure_ == true) {
            if (rComm.isMaster() == true) {
                this->saveFxcPureMatrix(RUN_RKS, this->m_nIteration, Fxc);
            }
        }

        // if (this->m_bIsHybrid == true) {
        //     Fxc += this->getFockExchange(this->m_dFockExchangeCoef * Ppq, RUN_RKS);
        //     this->XC_energy_ += DfXCFunctional::m_dFockExchangeEnergyAlpha;
        // }

        if (rComm.isMaster() == true) {
            this->saveFxcMatrix<TlSymmetricMatrix>(RUN_RKS, this->m_nIteration, Fxc);
        }
    }
    break;

    case METHOD_UKS: {
        TlSymmetricMatrix PApq, PBpq;
        if (rComm.isMaster() == true) {
            if (this->m_bIsUpdateXC == true) {
                PApq = this->getDiffDensityMatrix<TlSymmetricMatrix>(RUN_UKS_ALPHA, this->m_nIteration);
                PBpq = this->getDiffDensityMatrix<TlSymmetricMatrix>(RUN_UKS_BETA, this->m_nIteration);
            } else {
                PApq = this->getPpqMatrix<TlSymmetricMatrix>(RUN_UKS_ALPHA, this->m_nIteration -1);
                PBpq = this->getPpqMatrix<TlSymmetricMatrix>(RUN_UKS_BETA, this->m_nIteration -1);
            }
        }
        rComm.broadcast(PApq);
        rComm.broadcast(PBpq);

        TlSymmetricMatrix FxcA(this->m_nNumOfAOs);
        TlSymmetricMatrix FxcB(this->m_nNumOfAOs);
        DfXCFunctional::getFxc(PApq, PBpq, &dfCalcGrid, &FxcA, &FxcB);
        if (this->isSaveFxcPure_ == true) {
            if (rComm.isMaster() == true) {
                this->saveFxcPureMatrix(RUN_UKS_ALPHA, this->m_nIteration, FxcA);
                this->saveFxcPureMatrix(RUN_UKS_BETA,  this->m_nIteration, FxcB);
            }
        }

        // if (this->m_bIsHybrid == true) {
        //     FxcA += this->getFockExchange(this->m_dFockExchangeCoef * PApq, RUN_UKS_ALPHA);
        //     FxcB += this->getFockExchange(this->m_dFockExchangeCoef * PBpq, RUN_UKS_BETA);
        //     this->XC_energy_ += (DfXCFunctional::m_dFockExchangeEnergyAlpha + DfXCFunctional::m_dFockExchangeEnergyBeta);
        // }

        if (rComm.isMaster() == true) {
            this->saveFxcMatrix<TlSymmetricMatrix>(RUN_UKS_ALPHA, this->m_nIteration, FxcA);
            this->saveFxcMatrix<TlSymmetricMatrix>(RUN_UKS_BETA, this->m_nIteration, FxcB);
        }
    }
    break;

    case METHOD_ROKS:
        CnErr.abort("sorry. ROKS method is not implemented. @ DfXCFunctional_Parallel::buildXcMatrix()");
        break;
    default:
        break;
    }

}


void DfXCFunctional_Parallel::buildXC_ScaLAPACK()
{
    DfCalcGridX_Parallel dfCalcGrid(this->pPdfParam_);

    switch (this->m_nMethodType) {
    case METHOD_RKS: {
        TlDistributeSymmetricMatrix Ppq;
        if (this->m_bIsUpdateXC == true) {
            Ppq = 0.5 * this->getDiffDensityMatrix<TlDistributeSymmetricMatrix>(RUN_RKS, this->m_nIteration);
        } else {
            Ppq = 0.5 * this->getPpqMatrix<TlDistributeSymmetricMatrix>(RUN_RKS, this->m_nIteration -1);
        }

        TlDistributeSymmetricMatrix Fxc(this->m_nNumOfAOs);
        DfXCFunctional::getFxc(Ppq, &dfCalcGrid, &Fxc);
        if (this->isSaveFxcPure_ == true) {
            this->saveFxcPureMatrix(RUN_RKS, this->m_nIteration, Fxc);
        }
        
        // if (this->m_bIsHybrid == true) {
        //     Fxc += this->getFockExchange(this->m_dFockExchangeCoef * Ppq, RUN_RKS);
        //     this->XC_energy_ += DfXCFunctional::m_dFockExchangeEnergyAlpha;
        // }

        this->saveFxcMatrix<TlDistributeSymmetricMatrix>(RUN_RKS, this->m_nIteration, Fxc);
    }
    break;

    case METHOD_UKS: {
        TlDistributeSymmetricMatrix PApq, PBpq;
        if (this->m_bIsUpdateXC == true) {
            PApq = this->getDiffDensityMatrix<TlDistributeSymmetricMatrix>(RUN_UKS_ALPHA, this->m_nIteration);
            PBpq = this->getDiffDensityMatrix<TlDistributeSymmetricMatrix>(RUN_UKS_BETA, this->m_nIteration);
        } else {
            PApq = this->getPpqMatrix<TlDistributeSymmetricMatrix>(RUN_UKS_ALPHA, this->m_nIteration -1);
            PBpq = this->getPpqMatrix<TlDistributeSymmetricMatrix>(RUN_UKS_BETA, this->m_nIteration -1);
        }

        TlDistributeSymmetricMatrix FxcA(this->m_nNumOfAOs);
        TlDistributeSymmetricMatrix FxcB(this->m_nNumOfAOs);
        DfXCFunctional::getFxc(PApq, PBpq, &dfCalcGrid, &FxcA, &FxcB);
        if (this->isSaveFxcPure_ == true) {
            this->saveFxcPureMatrix(RUN_UKS_ALPHA, this->m_nIteration, FxcA);
            this->saveFxcPureMatrix(RUN_UKS_BETA,  this->m_nIteration, FxcB);
        }

        // if (this->m_bIsHybrid == true) {
        //     FxcA += this->getFockExchange(this->m_dFockExchangeCoef * PApq, RUN_UKS_ALPHA);
        //     FxcB += this->getFockExchange(this->m_dFockExchangeCoef * PBpq, RUN_UKS_BETA);
        //     this->XC_energy_ += (DfXCFunctional::m_dFockExchangeEnergyAlpha + DfXCFunctional::m_dFockExchangeEnergyBeta);
        // }

        this->saveFxcMatrix<TlDistributeSymmetricMatrix>(RUN_UKS_ALPHA, this->m_nIteration, FxcA);
        this->saveFxcMatrix<TlDistributeSymmetricMatrix>(RUN_UKS_BETA, this->m_nIteration, FxcB);
    }
    break;

    case METHOD_ROKS:
        CnErr.abort("sorry. ROKS method is not implemented. @ DfXCFunctional_Parallel::buildXcMatrix()");
        break;
    default:
        break;
    }
}


DfEri2* DfXCFunctional_Parallel::getDfEri2()
{
    DfEri2* pDfEri2 = new DfEri2_Parallel(this->pPdfParam_);

    return pDfEri2;
}


DfEriX* DfXCFunctional_Parallel::getDfEriXObject()
{
    this->logger(" use new engine\n");
    DfEriX* pDfEriX = new DfEriX_Parallel(this->pPdfParam_);

    return pDfEriX;
}


DfTwoElectronIntegral* DfXCFunctional_Parallel::getDfTwoElectronIntegral()
{
#ifdef USE_OLD_ERI_ENGINE
    DfTwoElectronIntegral* pDfTEI = new DfTwoElectronIntegral_Parallel(this->pPdfParam_);
    return pDfTEI;
#else
    this->log_.critical("cannot use old two-electron integral engine.");
    this->log_.critical("please check USE_OLD_TEI_ENGINE flag.");
    abort();
    return NULL;
#endif // USE_OLD_ERI_ENGINE
}


TlSymmetricMatrix DfXCFunctional_Parallel::getFockExchange(const TlSymmetricMatrix& deltaP,
                                                           const RUN_TYPE nRunType)
{
    this->loggerTime(" build fock exchange (on memory)");

    TlCommunicate& rComm = TlCommunicate::getInstance();
    TlSymmetricMatrix FxHF = DfXCFunctional::getFockExchange(deltaP);

    // update法を使用する場合は、前回のFock交換項を加算する
    if (this->m_bIsUpdateXC == true) {
        if (rComm.isMaster() == true) {
            if (this->m_nIteration >= 2) {
                FxHF += this->getHFxMatrix<TlSymmetricMatrix>(nRunType, this->m_nIteration -1);
            }
            this->saveHFxMatrix(nRunType, this->m_nIteration, FxHF);
        }
        rComm.broadcast(FxHF);
    }

    // calc energy
    {
        TlSymmetricMatrix P;
        if (rComm.isMaster() == true) {
            P = this->getPpqMatrix<TlSymmetricMatrix>(nRunType, this->m_nIteration -1);
        }
        rComm.broadcast(P);

        if (nRunType == RUN_RKS) {
            DfXCFunctional::m_dFockExchangeEnergyAlpha = 0.5 * DfXCFunctional::getFockExchangeEnergy(P, FxHF);
        } else if (nRunType == RUN_UKS_ALPHA) {
            DfXCFunctional::m_dFockExchangeEnergyAlpha = 0.5 * DfXCFunctional::getFockExchangeEnergy(P, FxHF);
        } else if (nRunType == RUN_UKS_BETA) {
            DfXCFunctional::m_dFockExchangeEnergyBeta = 0.5 * DfXCFunctional::getFockExchangeEnergy(P, FxHF);
        }
    }

    this->loggerTime(" end to build fock exchange.");
    return FxHF;
}


TlDistributeSymmetricMatrix
DfXCFunctional_Parallel::getFockExchange(const TlDistributeSymmetricMatrix& deltaP,
                                         const RUN_TYPE runType)
{
    this->loggerTime(" build fock exchange using distribute matrix");

    TlDistributeSymmetricMatrix FxHF = this->getFockExchange(deltaP);
    if (this->m_bIsUpdateXC == true) {
        if (this->m_nIteration >= 2) {
            FxHF += this->getHFxMatrix<TlDistributeSymmetricMatrix>(runType, this->m_nIteration -1);
        }
        this->saveHFxMatrix(runType, this->m_nIteration, FxHF);
    }

    // calc energy
    {
        TlDistributeSymmetricMatrix P = this->getPpqMatrix<TlDistributeSymmetricMatrix>(runType, this->m_nIteration -1);

        if (runType == RUN_RKS) {
            DfXCFunctional::m_dFockExchangeEnergyAlpha = 0.5 * this->getFockExchangeEnergy(P, FxHF);
        } else if (runType == RUN_UKS_ALPHA) {
            DfXCFunctional::m_dFockExchangeEnergyAlpha = 0.5 * this->getFockExchangeEnergy(P, FxHF);
        } else if (runType == RUN_UKS_BETA) {
            DfXCFunctional::m_dFockExchangeEnergyBeta = 0.5 * this->getFockExchangeEnergy(P, FxHF);
        }
    }

    this->loggerTime(" end to build fock exchange.");
    return FxHF;

  

//     deltaP.save("deltaP.tmp");

//     TlCommunicate& rComm = TlCommunicate::getInstance();
//     TlSymmetricMatrix tmpP;
//     if (rComm.isMaster() == true) {
//         tmpP.load("deltaP.tmp");
//     }
//     rComm.broadcast(tmpP);

//     TlSymmetricMatrix tmpF = this->getFockExchange(tmpP, runType);
//     tmpF.save("FxHF.tmp");

//     TlDistributeSymmetricMatrix FxHF;
//     FxHF.load("FxHF.tmp");
    
    // return FxHF;
}


double DfXCFunctional_Parallel::getFockExchangeEnergy(const TlDistributeSymmetricMatrix& P,
                                                      const TlDistributeSymmetricMatrix& HFx)
{
    TlDistributeSymmetricMatrix tmp(P);
    double energy = tmp.dot(HFx).sum();
    
    return energy;
}


TlDistributeSymmetricMatrix DfXCFunctional_Parallel::getFockExchange(const TlDistributeSymmetricMatrix& P)
{
    TlDistributeSymmetricMatrix FxHF(this->m_nNumOfAOs);

    if (this->isUseNewEngine_ == true) {
        this->logger(" use new engine\n");
        DfEriX_Parallel dfEriX(this->pPdfParam_);
        dfEriX.getK_D(P, &FxHF);
    } else {
#ifdef USE_OLD_ERI_ENGINE
        DfTwoElectronIntegral_Parallel dfTEI(this->pPdfParam_);
        if (this->m_bUseRTmethod != true) {
            this->logger(" force into using RT-method replaced by integral-driven.\n");
        }
        this->logger(" using RT method\n");
        dfTEI.getContractKMatrixByRTmethod(P, &FxHF);
#endif // USE_OLD_ERI_ENGINE
    }
    
    return FxHF;
}


double DfXCFunctional_Parallel::getGrimmeDispersionEnergy()
{
    if (this->isCalcd_E_disp_ != true) {
        TlCommunicate& rComm = TlCommunicate::getInstance();
        if (rComm.isMaster() == true) {
            DfXCFunctional::getGrimmeDispersionEnergy();
        }
        rComm.broadcast(this->isCalcd_E_disp_);
        rComm.broadcast(this->E_disp_);
    }
    
    return this->E_disp_;
}

