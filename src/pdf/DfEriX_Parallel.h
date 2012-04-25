#ifndef DFERIX_PARALLEL_H
#define DFERIX_PARALLEL_H

#include "DfEriX.h"
#include "TlDistributeSymmetricMatrix.h"

class DfEriX_Parallel : public DfEriX {
public:
    DfEriX_Parallel(TlSerializeData* pPdfParam);
    virtual ~DfEriX_Parallel();

public:
    void getJ(const TlDistributeSymmetricMatrix& P,
              TlDistributeVector* pRho) {
        this->getJ_D(P, pRho);
    }
    virtual void getJ(const TlSymmetricMatrix& P, TlVector* pRho) {
        DfEriX::getJ(P, pRho);
    }
                
    void getJ_D(const TlVector& rho, TlDistributeSymmetricMatrix* pJ);
    
    /// J([pq | rs])
    void getJpq_D(const TlDistributeSymmetricMatrix& P,
                  TlDistributeSymmetricMatrix* pJpq);

    /// J([alpha | beta])
    void getJab_D(TlDistributeSymmetricMatrix* pJab);

    /// K
    void getK_D(const TlDistributeSymmetricMatrix& P,
                TlDistributeSymmetricMatrix* pK);
    
public:
    void getdeltaHpqA(const TlVector& rho, TlSymmetricMatrix& P) {
        //this->getJ(rho, &P);
        DfEriX::getJ(rho, &P);
    }

    void getdeltaHpqA(const TlVector& rho, TlDistributeSymmetricMatrix& P) {
        this->getJ_D(rho, &P);
    }

protected:
    void getJ_D(const TlDistributeSymmetricMatrix& P,
                TlDistributeVector* pRho);

    virtual DfTaskCtrl* getDfTaskCtrlObject() const;

    virtual void finalize(TlMatrix* pMtx);
    virtual void finalize(TlSymmetricMatrix* pMtx);
    virtual void finalize(TlVector* pVct);

    virtual TlSparseSymmetricMatrix makeSchwarzTable(const TlOrbitalInfoObject& orbitalInfo);
    
    void waitAnotherProcs(const TlDistributeSymmetricMatrix& P);

protected:
    // 非同期通信により電子密度行列を送受信する
    void getJ_D_BG(const TlDistributeSymmetricMatrix& P,
                   TlDistributeVector* pRho);
    
    void getJ_D_local(const TlDistributeSymmetricMatrix& P,
                      TlDistributeVector* pRho);

    void getJ_part2(const TlOrbitalInfo& orbitalInfo,
                    const TlOrbitalInfo_Density& orbitalInfo_Density,
                    const ShellArrayTable& shellArrayTable_Density,
                    const std::vector<DfTaskCtrl::Task2>& taskList,
                    // const TlSparseSymmetricMatrix& schwarzTable,
                    const TlDistributeMatrix& P, TlVector* pRho);

    void getK_D_BG(const TlDistributeSymmetricMatrix& P,
                   TlDistributeSymmetricMatrix* pK);

    void getK_D_local(const TlDistributeSymmetricMatrix& P,
                      TlDistributeSymmetricMatrix* pK);
    
    // void getK_integralDriven_part2(const TlOrbitalInfoObject& orbitalInfo,
    //                                const std::vector<DfTaskCtrl::Task4>& taskList,
    //                                const TlDistributeMatrix& P, TlMatrixObject* pK);
    
    // void storeK_integralDriven2(const index_type shellIndexP, const int maxStepsP,
    //                             const index_type shellIndexQ, const int maxStepsQ,
    //                             const index_type shellIndexR, const int maxStepsR,
    //                             const index_type shellIndexS, const int maxStepsS,
    //                             const DfEriEngine& engine,
    //                             const TlDistributeMatrix& P,
    //                             TlMatrixObject* pK);
    
    void expandLocalDensityMatrix(const TlDistributeSymmetricMatrix& P,
                                  const TlOrbitalInfo& orbInfo,
                                  TlMatrix* pLocalP,
                                  std::vector<index_type>* pRowIndexes,
                                  std::vector<index_type>* pColIndexes);
    std::vector<index_type>
    getExpandIndexes(const std::vector<index_type>& refIndexes,
                     const TlOrbitalInfo& orbInfo);
    
protected:
    enum {
        TAG_ALL_PROC_FINISHED = 9999
    };

    enum CalcMode {
        CalcMode_Default = 0,
        CalcMode_BackGroundTransport = 1,
        CalcMode_UsingLocalMatrix = 2
    };
    
private:
    int calcMode_;
};

#endif // DFERIX_PARALLEL_H
