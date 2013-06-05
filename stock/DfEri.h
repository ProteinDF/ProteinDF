#ifndef DFERI_H
#define DFERI_H

#include <iostream>
#include <list>
#include <vector>

#include "DfObject.h"
#include "TlVector.h"
#include "TlMatrix.h"
#include "TlSymmetricMatrix.h"
#include "TlSparseMatrix.h"
#include "TlSparseSymmetricMatrix.h"
#include "TlPosition.h"
#include "TlOrbitalInfo.h"
#include "TlOrbitalInfo_Density.h"

/// 3中心積分[pq|alpha]、2中心積分[alpha|beta]を行うクラス
class DfEri : public DfObject {
public:
    DfEri(TlSerializeData* pPdfParam);
    virtual ~DfEri();

public:
    /// [pq|alpha]とP からtalphaの計算を行う
    virtual void getDeltaT(const TlSymmetricMatrix& deltaPpq, TlVector* pDeltaT);

    /// [pq|alpha]の計算を行う
    virtual void getdeltaHpqA(const TlVector& deltaRho, TlSymmetricMatrix& deltaH);

    //void getDeltaHpqA2(const TlVector& deltaRho, TlSymmetricMatrix& deltaH);

    /// @param[in] dCutoffCoef カットオフにかける係数
    void getDeltaHpqAForEri2(const TlVector& deltaRho, TlSymmetricMatrix& deltaH, double dCutoffCoef = 1.0);

    /// [alpha|beta]の計算を行う
    virtual void getSab(TlSymmetricMatrix* pSab);

protected:
    void initialize();

    /// Schwartzの不等式で、cutoffできる場合はtrueを返す
    bool isCutoffUsingSchwartzInequality(const TlOrbitalInfo& orbitalInfo,
                                         const std::size_t p, const std::size_t q,
                                         const double dCutoffValue) const;

    void initializeCounter();
    void countupTotal(int ity, int jty);
    void countupCutoff(int ity, int jty);
    virtual void cutoffReport();
    void cutoffReport(const std::string& shell, std::size_t cutoffCount, std::size_t totalCount);

    void getSab_core(TlMatrixObject* pSab);

protected:
    void getDeltaT_core(const TlMatrixObject* deltaPpq, TlVectorObject* pDeltaT);

protected:
    /// [ss|s] type function
    void eriSSS(const int np,
                const double* px, const double* py, const double* pz,
                const double cx, const double cy, const double cz,
                const double cc, const double Gammainv, const double* Zpi, const double* HP, double** ERP);

    /// [ss|p] type function
    void eriSSP(const int np,
                const double* px, const double* py, const double* pz,
                const double cx, const double cy, const double cz,
                const double cc, const double Gammainv, const double* Zpi, const double* HP, double** ERP);

    /// [ss|d] type function
    void eriSSD(const int np,
                const double* px, const double* py, const double* pz,
                const double cx, const double cy, const double cz,
                const double cc, const double Gammainv, const double* Zpi, const double* HP, double** ERP);

    /// [ps|s] type function
    void eriPSS(const int np,
                const double* px, const double* py, const double* pz,
                const double* pax, const double* pay, const double* paz,
                const double cx, const double cy, const double cz,
                const double cc, const double Gammainv, const double* Zpi, const double* HP, double** ERP);

    /// [ps|p] type function
    void eriPSP(int np,
                const double* px, const double* py, const double* pz,
                const double* pax, const double* pay, const double* paz,
                const double cx, const double cy, const double cz,
                const double cc, const double Gammainv, const double* Zpi, const double* HP, double** ERP);

    /// [ps|d] type function
    void eriPSD(const int np,
                const double* px, const double* py, const double* pz,
                const double* pax, const double* pay, const double* paz,
                const double cx, const double cy, const double cz,
                const double cc, const double Gammainv, const double* Zpi, const double* HP, double** ERP);

    /// [pp|s] type function
    void eriPPS(const int np,
                const double* px, const double* py, const double* pz,
                const double* pax, const double* pay, const double* paz,
                const double* pbx, const double* pby, const double* pbz,
                const double cx, const double cy, const double cz,
                const double cc, const double Gammainv, const double* Zpi, const double* HP, double** ERP);

    /// [pp|p] type function
    void eriPPP(const int np,
                const double* px, const double* py, const double* pz,
                const double* pax, const double* pay, const double* paz,
                const double* pbx, const double* pby, const double* pbz,
                const double cx, const double cy, const double cz,
                const double cc, const double Gammainv, const double* Zpi, const double* HP, double** ERP);

    /// [pp|d] type function
    void eriPPD(const int np,
                const double* px, const double* py, const double* pz,
                const double* pax, const double* pay, const double* paz,
                const double* pbx, const double* pby, const double* pbz,
                const double cx, const double cy, const double cz,
                const double cc, const double Gammainv, const double* Zpi, const double* HP, double** ERP);

    /// [ds|s] type function
    void eriDSS(const int np,
                const double* px, const double* py, const double* pz,
                const double* pax, const double* pay, const double* paz,
                const double cx, const double cy, const double cz,
                const double cc, const double Gammainv, const double* Zpi, const double* HP, double** ERP);

    /// [ds|p] type function
    void eriDSP(const int np,
                const double* px, const double* py, const double* pz,
                const double* pax, const double* pay, const double* paz,
                const double cx, const double cy, const double cz,
                const double cc, const double Gammainv, const double* Zpi, const double* HP, double** ERP);

    /// [ds|d] type function
    void eriDSD(const int np,
                const double* px, const double* py, const double* pz,
                const double* pax, const double* pay, const double* paz,
                const double cx, const double cy, const double cz,
                const double cc, const double Gammainv, const double* Zpi, const double* HP, double** ERP);

    /// [dp|s] type function
    void eriDPS(const int np,
                const double* px, const double* py, const double* pz,
                const double* pax, const double* pay, const double* paz,
                const double* pbx, const double* pby, const double* pbz,
                const double cx, const double cy, const double cz,
                const double cc, const double Gammainv, const double* Zpi, const double* HP, double** ERP);

    /// [dp|p] type function
    void eriDPP(const int np,
                const double* px, const double* py, const double* pz,
                const double* pax, const double* pay, const double* paz,
                const double* pbx, const double* pby, const double* pbz,
                const double cx, const double cy, const double cz,
                const double cc, const double Gammainv, const double* Zpi, const double* HP, double** ERP);

    /// [dp|d] type function
    void eriDPD(const int np,
                const double* px, const double* py, const double* pz,
                const double* pax, const double* pay, const double* paz,
                const double* pbx, const double* pby, const double* pbz,
                const double cx, const double cy, const double cz,
                const double cc, const double Gammainv, const double* Zpi, const double* HP, double** ERP);

    /// [dd|s] type function
    void eriDDS(const int np,
                const double* px, const double* py, const double* pz,
                const double* pax, const double* pay, const double* paz,
                const double* pbx, const double* pby, const double* pbz,
                const double cx, const double cy, const double cz,
                const double cc, const double Gammainv, const double* Zpi, const double* HP, double** ERP);

    /// [dd|p] type function
    void eriDDP(const int np,
                const double* px, const double* py, const double* pz,
                const double* pax, const double* pay, const double* paz,
                const double* pbx, const double* pby, const double* pbz,
                const double cx, const double cy, const double cz,
                const double cc, const double Gammainv, const double* Zpi, const double* HP, double** ERP);

    /// [dd|d] type function
    void eriDDD(const int np,
                const double* px, const double* py, const double* pz,
                const double* pax, const double* pay, const double* paz,
                const double* pbx, const double* pby, const double* pbz,
                const double cx, const double cy, const double cz,
                const double cc, const double Gammainv, const double* Zpi, const double* HP, double** ERP);

    void getMemory();

    /// 分子積分を行うためのテーブルを作る
    void makeTable();

    /// 補助積分メインルーチン
    int auxSet();

    /// 補助積分の計算
    void fmtSet();

    /// 補助積分の計算
    void fmtRecursive(int, int, int);

    /// 3 中心積分のメインルーチン
    void ericalc(const int nqA, const int nqB, const int nqC,
                 const int npA, const double* CA, const double* ZA, const double* A,
                 const int npB, const double* CB, const double* ZB, const double* B,
                 const double CC, const double ZC, const double* C,
                 double* ERI);

    /// [alpha|beta]の分子積分を行う
    void sabcalc(const std::vector<IJShellPair>& IJShellList, TlMatrixObject* Sab);

    /// [pq|alpha]と⊿P から⊿talphaの計算を行う
    void ericalcDT(const std::vector<IJShellPair>& IJShellList, const TlMatrixObject*, TlVectorObject*);

    /// delta rho[pq|alpha]の計算を行う．
    void ericalcDH(const std::vector<IJShellPair>& IJShellList, const TlVector& deltaRho, TlMatrixObject* deltaH);

protected:
    static const double FPAI; /// 2*PAI**(5/2)
    static const double D33;  /// 1/0.03
    static const double SQR3I;    /// 1/sqrt(3)
    static const double SQR3I2;   /// 2/sqrt(3)
    static const int MAX_TYPE;
    
    TlOrbitalInfo* pOrbitalInfo_;
    TlOrbitalInfo_Density* pOrbitalInfo_Density_;

    std::vector<std::vector<std::size_t> > shellList_;
    std::vector<std::vector<std::size_t> > shellList_Dens_;
    
    double m_dStoredCutValue; // 入力されたカットオフ値
    double m_dCutValue;       // 動的に変化するカットオフ値

    // for auxSet & fmtSet
    double* TF;
    double* RMI;
    double* GA;
    double* EDAT;
    double* SDAT;
    double** FDAT;
    double** ADAT;

    std::size_t TcountSS;    // Total number of SS shell
    std::size_t TcountPS;
    std::size_t TcountPP;
    std::size_t TcountDS;
    std::size_t TcountDP;
    std::size_t TcountDD;
    std::size_t CcountSS;    // Cut off number of SS shell
    std::size_t CcountPS;
    std::size_t CcountPP;
    std::size_t CcountDS;
    std::size_t CcountDP;
    std::size_t CcountDD;

protected:
    virtual std::vector<IJShellPair> getQueue(const TlOrbitalInfo* pOrbitalInfo,
                                              const std::vector<std::vector<std::size_t> >& shellList,
                                              const int maxNumOfPairs, const bool initialize = false);

    virtual void finalizeIntegral(TlSymmetricMatrix& rMatrix);
    virtual void finalizeIntegral(TlVector& rVector);

    bool isCutoff_IK(const std::size_t orbA, const std::size_t orbC) const;

protected:
    /// 計算順序を管理する
    ///
    /// @param[in] initialize 計算キューを初期化する(true)
    /// @retval true  計算キューに未計算項目がある
    /// @retval false 計算キューに未計算項目がない
    bool getQueueX(const std::vector<std::vector<std::size_t> >& shellList,
                   int* pShellType, index_type* pShell,
                   const bool initialize = false);
    
protected:
    int blockSize_;
};

#endif // DFERI_H

