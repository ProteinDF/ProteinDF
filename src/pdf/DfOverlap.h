#ifndef DFOVERLAP_H
#define DFOVERLAP_H

#include <iostream>
#include <vector>

#include "DfObject.h"
#include "Fl_Int_Pqg.h"
#include "TlVector.h"
#include "TlSymmetricMatrix.h"
#include "TlOrbitalInfo.h"
#include "TlOrbitalInfo_Density.h"
#include "TlOrbitalInfo_XC.h"

/// 3 中心積分[pq gamma], 2 中心積分[pq]、[alpha beta]、[gamma delta]、1 中心積分[alpha]の計算を行うクラス
class DfOverlap : public DfObject {
protected:
    typedef std::vector<std::vector<std::size_t> > ShellListType;

public:
    DfOverlap(TlSerializeData* pPdfParam);
    virtual ~DfOverlap();

public:
    /// myu[pq gamma], epsilon[pq gamma]の計算を行う
    virtual void getdeltaHpqG(const TlVector&, TlSymmetricMatrix&);

    /// myu[pq gamma], epsilon[pq gamma]の計算を行う
    void getdeltaHpqG(const TlVector&, const TlVector&, TlSymmetricMatrix&, TlSymmetricMatrix&);

    /// [pq]の計算を行う
    virtual void getSpq(TlSymmetricMatrix* pSpq);

    TlMatrix getSpq(const TlOrbitalInfo& orbInfo1,
                    const TlOrbitalInfo& orbInfo2);
    
    /// [gamma delta]の計算を行う
    virtual void getSgd(TlSymmetricMatrix* pSgd);

    /// [alpha beta]の計算を行う
    virtual void getSab2(TlSymmetricMatrix* pSab2);

    /// [alpha]の計算を行う
    virtual void getNa(TlVector* pNa);

protected:
    void getNa_core(TlVector* pNa);

    void getSpq_core(TlMatrixObject* pSpq);

    void getSpq_core(const TlOrbitalInfo& orbInfo1,
                     const TlOrbitalInfo& orbInfo2,
                     TlMatrixObject* pSpq);

    void getSab2_core(TlMatrixObject* pSab2);
    void getSgd_core(TlMatrixObject* pSgd);

protected:
    void countupTotal(int ity, int jty);
    void countupCutoff(int ity, int jty);

    virtual void cutoffReport();
    void cutoffReport(const std::string& shell, std::size_t cutoffCount, std::size_t totalCount);

    //virtual void clockUp(const std::string& str);

    // cut off by schwaltz inequality
    bool isCutoffUsingSchwartzInequality(const TlOrbitalInfo& orbInfo,
                                         const std::size_t orb1,
                                         const std::size_t orb2) const;

    bool isCutoffUsingSchwartzInequality(const TlOrbitalInfo& orbInfo1,
                                         const TlOrbitalInfo& orbInfo2,
                                         const std::size_t orb1,
                                         const std::size_t orb2) const;
    
    /// 分子積分を行うためのテーブルを作る
    void makeTable();

    ShellListType makeShellList(const TlOrbitalInfoObject& orbInfo);
    
    /// 3 中心積分のメインルーチン
    void ovpqgcalc(int nqA, int nqB, int nqC,
                   int npA, const double* CA, const double* ZA, const double* A,
                   int npB, const double* CB, const double* ZB, const double* B,
                   double CC, double ZC, const double* C, double* OVP);

    /// [pq] [gamma delta]の分子積分を行う
    void spqcalc(const TlOrbitalInfoObject& orbInfo,
                 const std::vector<IJShellPair>& IJShellList, TlMatrixObject* Spq);

    /// [pq] [gamma delta]の分子積分を行う
    /// 
    /// 2種類の基底関数を用いるバージョン。
    /// 受け取る行列(Spq)は非対称行列になる。 
    void spqcalc(const TlOrbitalInfoObject& orbInfo1,
                 const TlOrbitalInfoObject& orbInfo2,
                 const std::vector<IJShellPair>& IJShellList, TlMatrixObject* Spq);
    
    /// [alpha]の分子積分を行う
    void calcNa(const std::vector<std::size_t>&, TlVector*);

    /// myu[pq gamma], epsilon[pq gamma]の計算を行う
    void ovpqgDH(const std::vector<IJShellPair>&, const TlVector&, TlMatrixObject*);

    /// myu[pq gamma], epsilon[pq gamma]の計算を行う
    void ovpqgDH(const std::vector<IJShellPair>&, const TlVector&, const TlVector&,
                 TlMatrixObject*, TlMatrixObject*);

    void ovpSSS(const int np,
                const double* px, const double* py, const double* pz, const double Gamma,
                const double cx, const double cy, const double cz, const double cc,
                const double* Zp, const double* HP, double** ERP);
    void ovpSSP(const int np,
                const double* px, const double* py, const double* pz,
                const double Gamma,
                const double cx, const double cy, const double cz,
                const double cc,
                const double* Zp, const double* HP, double** ERP);
    void ovpSSD(const int np,
                const double* px, const double* py, const double* pz,
                const double Gamma,
                const double cx, const double cy, const double cz,
                const double cc,
                const double* Zp, const double* HP, double** ERP);
    void ovpPSS(const int np,
                const double* px, const double* py, const double* pz,
                const double* pax, const double* pay, const double* paz,
                const double Gamma,
                const double cx, const double cy, const double cz,
                const double cc,
                const double* Zp, const double* HP, double** ERP);
    void ovpPSP(const int np,
                const double* px, const double* py, const double* pz,
                const double* pax, const double* pay, const double* paz,
                const double Gamma,
                const double cx, const double cy, const double cz,
                const double cc,
                const double* Zp, const double* HP, double** ERP);
    void ovpPSD(const int np,
                const double* px, const double* py, const double* pz,
                const double* pax, const double* pay, const double* paz,
                const double Gamma,
                const double cx, const double cy, const double cz,
                const double cc,
                const double* Zp, const double* HP, double** ERP);
    void ovpPPS(const int np,
                const double* px, const double* py, const double* pz,
                const double* pax, const double* pay, const double* paz,
                const double* pbx, const double* pby, const double* pbz,
                const double Gamma,
                const double cx, const double cy, const double cz,
                const double cc,
                const double* Zp, const double* HP, double** ERP);
    void ovpPPP(const int np,
                const double* px, const double* py, const double* pz,
                const double* pax, const double* pay, const double* paz,
                const double* pbx, const double* pby, const double* pbz,
                const double Gamma,
                const double cx, const double cy, const double cz,
                const double cc,
                const double* Zp, const double* HP, double** ERP);
    void ovpPPD(const int np,
                const double* px, const double* py, const double* pz,
                const double* pax, const double* pay, const double* paz,
                const double* pbx, const double* pby, const double* pbz,
                const double Gamma,
                const double cx, const double cy, const double cz,
                const double cc,
                const double* Zp, const double* HP, double** ERP);
    void ovpDSS(const int np,
                const double* px, const double* py, const double* pz,
                const double* pax, const double* pay, const double* paz,
                const double Gamma,
                const double cx, const double cy, const double cz,
                const double cc,
                const double* Zp, const double* HP, double** ERP);
    void ovpDSP(const int np,
                const double* px, const double* py, const double* pz,
                const double* pax, const double* pay, const double* paz,
                const double Gamma,
                const double cx, const double cy, const double cz,
                const double cc,
                const double* Zp, const double* HP, double** ERP);
    void ovpDSD(const int np,
                const double* px, const double* py, const double* pz,
                const double* pax, const double* pay, const double* paz,
                const double Gamma,
                const double cx, const double cy, const double cz,
                const double cc,
                const double* Zp, const double* HP, double** ERP);
    void ovpDPS(const int np,
                const double* px, const double* py, const double* pz,
                const double* pax, const double* pay, const double* paz,
                const double* pbx, const double* pby, const double* pbz,
                const double Gamma,
                const double cx, const double cy, const double cz,
                const double cc,
                const double* Zp, const double* HP, double** ERP);
    void ovpDPP(const int np,
                const double* px, const double* py, const double* pz,
                const double* pax, const double* pay, const double* paz,
                const double* pbx, const double* pby, const double* pbz,
                const double Gamma,
                const double cx, const double cy, const double cz,
                const double cc,
                const double* Zp, const double* HP, double** ERP);
    void ovpDPD(const int np,
                const double* px, const double* py, const double* pz,
                const double* pax, const double* pay, const double* paz,
                const double* pbx, const double* pby, const double* pbz,
                const double Gamma,
                const double cx, const double cy, const double cz,
                const double cc,
                const double* Zp, const double* HP, double** ERP);
    void ovpDDS(const int np,
                const double* px, const double* py, const double* pz,
                const double* pax, const double* pay, const double* paz,
                const double* pbx, const double* pby, const double* pbz,
                const double Gamma,
                const double cx, const double cy, const double cz,
                const double cc,
                const double* Zp, const double* HP, double** ERP);
    void ovpDDP(const int np,
                const double* px, const double* py, const double* pz,
                const double* pax, const double* pay, const double* paz,
                const double* pbx, const double* pby, const double* pbz,
                const double Gamma,
                const double cx, const double cy, const double cz,
                const double cc,
                const double* Zp, const double* HP, double** ERP);
    void ovpDDD(const int np,
                const double* px, const double* py, const double* pz,
                const double* pax, const double* pay, const double* paz,
                const double* pbx, const double* pby, const double* pbz,
                const double Gamma,
                const double cx, const double cy, const double cz,
                const double cc,
                const double* Zp, const double* HP, double** ERP);

protected:
    static const double SQR3I;            // 1/sqrt(3)
    static const double SQR3I2;           // 2/sqrt(3)

    TlOrbitalInfo* pOrbitalInfo_;
    TlOrbitalInfo_Density* pOrbitalInfo_Density_;
    TlOrbitalInfo_XC* pOrbitalInfo_XC_;

    double cutvalue;      // cut value

    unsigned long TcountSS;     // Total number of SS shell
    unsigned long TcountPS;
    unsigned long TcountPP;
    unsigned long TcountDS;
    unsigned long TcountDP;
    unsigned long TcountDD;
    unsigned long CcountSS;     // Cut off number of SS shell
    unsigned long CcountPS;
    unsigned long CcountPP;
    unsigned long CcountDS;
    unsigned long CcountDP;
    unsigned long CcountDD;

protected:
    virtual std::vector<IJShellPair> getQueue(const TlOrbitalInfo* pOrbitalInfo,
                                              const ShellListType& shellList,
                                              const int maxScore,
                                              const bool initialize = false);

    /// jobの割り当て
    /// 非対称用
    virtual std::vector<IJShellPair> getQueue(const TlOrbitalInfo& orbitalInfo1,
                                              const TlOrbitalInfo& orbitalInfo2,
                                              const ShellListType& shellList1,
                                              const ShellListType& shellList2,
                                              const int maxScore,
                                              const bool initialize = false);

    virtual void finalizeIntegral(TlSymmetricMatrix& rMatrix);
    virtual void finalizeIntegral(TlVector& rVector);

protected:
    int blockSize_;

    ShellListType shellList_; // shellList_[shell type][index] = AO index
    ShellListType shellList_Dens_; // shellList_[shell type][index] = AO index
    ShellListType shellList_XC_; // shellList_[shell type][index] = AO index


private:
//     void getSpq_core2(TlMatrixObject* pSpq);
//     void getSpq_new();
    
};

#endif

