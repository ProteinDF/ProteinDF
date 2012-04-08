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
    struct CholeskyBasis {
        explicit CholeskyBasis(index_type in_p = 0, index_type in_q = 0)
            : p(in_p), q(in_q) {
        }
        index_type p;
        index_type q;
    };

    ///
    /// 0, 1, 2,... の順に対応するp, qのペアを格納する
    typedef std::vector<CholeskyBasis> FittingTbl;
    
    
public:
    DfCD(TlSerializeData* pPdfParam);
    virtual ~DfCD();

public:
    void makeSuperMatrix();
    void makeSuperMatrix_exact();

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

    void checkCholeskyFactorization(const TlSymmetricMatrix& G);
    
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
    DfEriEngine* pEriEngines_;
    double cutoffThreshold_;
    double cutoffEpsilon3_;

private:
#ifdef CHECK_LOOP
    TlSymmetricMatrix check;
#endif // CHECK_LOOP
};

#endif // DFCD_H
