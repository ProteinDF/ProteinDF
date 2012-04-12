#ifndef DFCD_H
#define DFCD_H

#include "DfObject.h"
#include "DfTaskCtrl.h"
#include "TlSparseSymmetricMatrix.h"
#include "TlSymmetricMatrix.h"

class TlOrbitalInfo;
class DfEriEngine;

// #define CHECK_LOOP // 計算ループ構造のチェック

class DfCD : public DfObject
{
public:
    DfCD(TlSerializeData* pPdfParam);
    virtual ~DfCD();

public:
    void makeSuperMatrix();
    void makeSuperMatrix_exact();
    
    void getJ(TlSymmetricMatrix *pJ);
    void getK(const RUN_TYPE runType,
              TlSymmetricMatrix *pK);

protected:
    void createEngines();
    void destroyEngines();
    
    void makeSuperMatrix_kernel(const TlOrbitalInfo& orbitalInfo,
                                const std::vector<DfTaskCtrl::Task4>& taskList,
                                TlSymmetricMatrix* pG);
    void storeG(const index_type shellIndexP, const int maxStepsP,
                const index_type shellIndexQ, const int maxStepsQ,
                const index_type shellIndexR, const int maxStepsR,
                const index_type shellIndexS, const int maxStepsS,
                const DfEriEngine& engine,
                TlSymmetricMatrix* pG);

    TlSparseSymmetricMatrix makeSchwarzTable(const TlOrbitalInfoObject& orbitalInfo);

    std::size_t index(index_type p, index_type q) const;
    DfTaskCtrl* getDfTaskCtrlObject() const;
    void finalize(TlMatrix* pMtx);

protected:
    void makeSuperMatrix_screening();
    void makeSuperMatrix_noScreening();

protected:
    void calcPQPQ(const TlOrbitalInfoObject& orbitalInfo,
                  TlSparseSymmetricMatrix *pSchwarzTable,
                  std::vector<index_type> *pI2PQ);
    void calcPQPQ_kernel(const TlOrbitalInfoObject& orbitalInfo,
                         const std::vector<DfTaskCtrl::Task2>& taskList,
                         TlSparseSymmetricMatrix *pSchwarzTable,
                         std::vector<index_type> *pI2PQ);
    void makeSuperMatrix_kernel2(const TlOrbitalInfo& orbitalInfo,
                                 const std::vector<DfTaskCtrl::Task4>& taskList,
                                 const std::vector<index_type>& PQ2I,
                                 TlSymmetricMatrix* pG);
    void storeG2(const index_type shellIndexP, const int maxStepsP,
                 const index_type shellIndexQ, const int maxStepsQ,
                 const index_type shellIndexR, const int maxStepsR,
                 const index_type shellIndexS, const int maxStepsS,
                 const std::vector<index_type>& PQ2I,
                 const DfEriEngine& engine,
                 TlSymmetricMatrix* pG);
    void saveI2PQ(const std::vector<int>& I2PQ);
    std::vector<int> getI2PQ();

    TlSymmetricMatrix getCholeskyVector(const TlVector& L_col,
                                        const std::vector<int>& I2PQ);

protected: // for exact
    typedef std::vector<index_type> ShellArray;
    typedef std::vector<ShellArray> ShellArrayTable;
    struct ShellPair {
    public:
        ShellPair(index_type index1 =0, index_type index2 =0) : shellIndex1(index1), shellIndex2(index2) {
        }
        
    public:
        index_type shellIndex1;
        index_type shellIndex2;
    };
    typedef std::vector<ShellPair> ShellPairArray;
    typedef std::vector<ShellPairArray> ShellPairArrayTable;
    
    
    ShellArrayTable makeShellArrayTable(const TlOrbitalInfoObject& orbitalInfo);
    ShellPairArrayTable getShellPairArrayTable(const ShellArrayTable& shellArrayTable);
    static const int MAX_SHELL_TYPE;
    
protected:
    index_type numOfPQs_;

    DfEriEngine* pEriEngines_;

    double cutoffThreshold_;
    double cutoffEpsilon3_;
    double epsilon_;

private:
#ifdef CHECK_LOOP
    TlSymmetricMatrix check;
#endif // CHECK_LOOP
};

#endif // DFCD_H
