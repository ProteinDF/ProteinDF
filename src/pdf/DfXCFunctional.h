#ifndef DFXCFUNCTIONAL_H
#define DFXCFUNCTIONAL_H

#include "DfObject.h"
#include "TlSymmetricMatrix.h"

#include "DfFunctional_SHF.h"
#include "DfFunctional_SVWN.h"
#include "DfFunctional_Becke88.h"
#include "DfFunctional_B88LYP.h"
#include "DfFunctional_B3LYP.h"
#include "CnError.h"

class DfCalcGridX;
class DfEri2;
class DfEriX;
class DfTwoElectronIntegral;

// 交換・相関項の計算を行う
//
// 交換相関項はUpdate法で計算するが、
// 出力された行列Fxcは差分ではないことに注意すること。
class DfXCFunctional : public DfObject {
public:
    enum FUNCTIONAL_TYPE {
        UNDEFINED_TYPE,
        LDA,
        GGA
    };
    
    enum XC_TYPE {
        UNDEFINED_XC,
        HF,
        SHF,
        SVWN,
        HFB,
        BLYP,
        B3LYP
    };

public:
    DfXCFunctional(TlSerializeData* pPdfParam);
    virtual ~DfXCFunctional();

public:
    bool isHybridFunctional() const;
    double getFockExchangeCoefficient() const;

    FUNCTIONAL_TYPE getFunctionalType() const;
    XC_TYPE getXcType() const;
    
public:
    /// Kohn-Sham行列におけるXC項を作成する
    virtual void buildXcMatrix();

    /// XCエネルギーを返す
    ///
    /// XCエネルギー計算はKohn-Sham行列を作成するときに、
    /// その時点の密度行列を元に計算する
    double getEnergy() const {
        return this->XC_energy_;
    }

    virtual double getGrimmeDispersionEnergy();

protected:
    /// build F_xc matrix for RKS case
    template<typename SymmetricMatrixType, typename DfCalcGridClass>
    void getFxc(const SymmetricMatrixType& Ppq,
                DfCalcGridClass* pDfCalcGridObj,
                SymmetricMatrixType* pFxc);

    /// build F_xc matrix for UKS case
    template<typename SymmetricMatrixType, typename DfCalcGridClass>
    void getFxc(const SymmetricMatrixType& PpqA,
                const SymmetricMatrixType& PpqB,
                DfCalcGridClass* pDfCalcGridObj,
                SymmetricMatrixType* pFxcA,
                SymmetricMatrixType* pFxcB);

protected:
    /// Fockの交換項を求める
    virtual TlSymmetricMatrix getFockExchange(const TlSymmetricMatrix& prevP, RUN_TYPE nRunType);
    TlSymmetricMatrix getFockExchange(const TlSymmetricMatrix& P);

    double getFockExchangeEnergy(const TlSymmetricMatrix& P, const TlSymmetricMatrix& Ex);
    double getFockExchangeEnergy(const TlSymmetricMatrix& PA, const TlSymmetricMatrix& PB,
                                 const TlSymmetricMatrix& Ex);


    /// 電子数を数値積分により求めることにより、グリッドの精度をチェックする
    void checkGridAccuracy();
    
protected:
    virtual DfEri2* getDfEri2();
    virtual DfEriX* getDfEriXObject();
    virtual DfTwoElectronIntegral* getDfTwoElectronIntegral();

protected:
    FUNCTIONAL_TYPE functionalType_;
    XC_TYPE m_nXCFunctional;

    /// hybrid汎関数を使う(true)かどうか
    ///
    /// Fockの交換項を計算する必要がある場合はtrueになる
    bool m_bIsHybrid;

    /// Fockの交換項にかかる係数
    double m_dFockExchangeCoef;

    bool m_bRI_K;
    bool m_bUseRTmethod;
    bool isSaveFxcPure_;
    bool m_bDebugTEI;
    bool m_bPrintTEI;

    bool isUseNewEngine_;
    
    //static double m_dXC_Energy;
    double XC_energy_;
    static double m_dFockExchangeEnergyAlpha; // Update法を使うため、直前のenergyを保存
    static double m_dFockExchangeEnergyBeta;
    
    /// Grimmeの分散力補正
    static bool isCalcd_E_disp_;
    static double E_disp_;
};


// build F_xc matrix for RKS case
template<typename SymmetricMatrixType, typename DfCalcGridClass>
void DfXCFunctional::getFxc(const SymmetricMatrixType& P1,
                            DfCalcGridClass* pDfCalcGridObj,
                            SymmetricMatrixType* pFxc)
{
    assert(pDfCalcGridObj != NULL);
    assert(pFxc != NULL);

    switch (this->m_nXCFunctional) {
    case HF:
        this->XC_energy_ = 0.0;
        break;

    case SHF: {
        DfFunctional_SHF shf;
        this->XC_energy_ = pDfCalcGridObj->calcXCIntegForFockAndEnergy(P1,
                                                                       &shf,
                                                                       pFxc);
    }
        break;
        
    case SVWN: {
        DfFunctional_SVWN svwn;
        this->XC_energy_ = pDfCalcGridObj->calcXCIntegForFockAndEnergy(P1,
                                                                       &svwn,
                                                                       pFxc);
    }
        break;
        
    case HFB: {
        DfFunctional_Becke88 b88;
        this->XC_energy_ = pDfCalcGridObj->calcXCIntegForFockAndEnergy(P1,
                                                                       &b88,
                                                                       pFxc);
    }
        break;
        
    case BLYP: {
        DfFunctional_B88LYP blyp;
        this->XC_energy_ = pDfCalcGridObj->calcXCIntegForFockAndEnergy(P1,
                                                                       &blyp,
                                                                       pFxc);
    }
        break;
        
    case B3LYP: {
        // except HF exchange
        DfFunctional_B3LYP b3lyp;
        this->XC_energy_ = pDfCalcGridObj->calcXCIntegForFockAndEnergy(P1,
                                                                       &b3lyp,
                                                                       pFxc);
    }
        break;
        
    default:
        CnErr.abort("unknown XC functional. STOP. @DfXCFunctional::getFxc()");
    }
    
    this->log_.debug(TlUtils::format("DfXCFunctional::getFxc(): XC_energy=% 16.10f",
                                     this->XC_energy_));
}


template<typename SymmetricMatrixType, typename DfCalcGridClass>
void DfXCFunctional::getFxc(const SymmetricMatrixType& PpqA,
                            const SymmetricMatrixType& PpqB,
                            DfCalcGridClass* pDfCalcGridObj,
                            SymmetricMatrixType* pFxcA,
                            SymmetricMatrixType* pFxcB)
{
    assert(pDfCalcGridObj != NULL);
    assert(pFxcA != NULL);
    assert(pFxcB != NULL);
    
    switch (this->m_nXCFunctional) {
    case HF:
        this->XC_energy_ = 0.0;
        break;
        
    case SHF: {
        DfFunctional_SHF shf;
        this->XC_energy_ = pDfCalcGridObj->calcXCIntegForFockAndEnergy(PpqA, PpqB, &shf,
                                                                       pFxcA, pFxcB);
    }
        break;
        
    case SVWN: {
        DfFunctional_SVWN svwn;
        this->XC_energy_ = pDfCalcGridObj->calcXCIntegForFockAndEnergy(PpqA, PpqB, &svwn,
                                                                       pFxcA, pFxcB);
    }
        break;
        
    case HFB: {
        DfFunctional_Becke88 hfb;
        this->XC_energy_ = pDfCalcGridObj->calcXCIntegForFockAndEnergy(PpqA, PpqB, &hfb,
                                                                       pFxcA, pFxcB);
    }
        break;
        
    case BLYP: {
        DfFunctional_B88LYP blyp;
        this->XC_energy_ = pDfCalcGridObj->calcXCIntegForFockAndEnergy(PpqA, PpqB, &blyp,
                                                                       pFxcA, pFxcB);
    }
        break;
        
    case B3LYP: {
        // except HF exchange
        DfFunctional_B3LYP b3lyp;
        this->XC_energy_ = pDfCalcGridObj->calcXCIntegForFockAndEnergy(PpqA, PpqB, &b3lyp,
                                                                       pFxcA, pFxcB);
    }
        break;
        
    default:
        CnErr.abort("unknown XC functional. STOP. @DfXCFunctional::getFxc()");
    }
}


#endif // DFXCFUNCTIONAL_H
