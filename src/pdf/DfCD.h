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
    struct PQ_Pair {
    public:
        PQ_Pair(index_type index1 =0, index_type index2 =0) : shellIndex1(index1), shellIndex2(index2) {
            if (this->shellIndex1 > this->shellIndex2) {
                std::swap(this->shellIndex1, this->shellIndex2);
            }
        }

    public:
        index_type shellIndex1;
        index_type shellIndex2;
    };
    
    struct PQ_Pair_less {
    public:
        bool operator()(const PQ_Pair& a, const PQ_Pair& b) const {
            bool answer = false;
            if ((a.shellIndex1 < b.shellIndex1) ||
                ((a.shellIndex1 == b.shellIndex1) && (a.shellIndex2 < b.shellIndex2))) {
                answer = true;
            }
            return answer;
        }
    };
    typedef std::vector<PQ_Pair> PQ_PairArray;

    typedef std::vector<PQ_Pair> I2PQ_Type;
    //typedef std::map<PQ_Pair, std::size_t, PQ_Pair_less> PQ2I_Type;
    typedef std::vector<size_type> PQ2I_Type;

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

    virtual DfTaskCtrl* getDfTaskCtrlObject() const;

    virtual void finalize(TlSymmetricMatrix *pMat);
    virtual void finalize(TlSparseSymmetricMatrix *pMat);
    virtual void finalize_I2PQ(I2PQ_Type *pI2PQ);

protected:
    void makeSuperMatrix_screening();
    void makeSuperMatrix_noScreening();

protected:
    void calcPQPQ(const TlOrbitalInfoObject& orbitalInfo,
                  TlSparseSymmetricMatrix *pSchwarzTable,
                  PQ_PairArray *pI2PQ);
    void calcPQPQ_kernel(const TlOrbitalInfoObject& orbitalInfo,
                         const std::vector<DfTaskCtrl::Task2>& taskList,
                         TlSparseSymmetricMatrix *pSchwarzTable,
                         PQ_PairArray *pI2PQ);
    void makeSuperMatrix_kernel2(const TlOrbitalInfo& orbitalInfo,
                                 const std::vector<DfTaskCtrl::Task4>& taskList,
                                 const PQ2I_Type& PQ2I,
                                 TlSymmetricMatrix* pG);
    void storeG2(const index_type shellIndexP, const int maxStepsP,
                 const index_type shellIndexQ, const int maxStepsQ,
                 const index_type shellIndexR, const int maxStepsR,
                 const index_type shellIndexS, const int maxStepsS,
                 const PQ2I_Type& PQ2I,
                 const DfEriEngine& engine,
                 TlSymmetricMatrix* pG);
    void saveI2PQ(const PQ_PairArray& I2PQ);
    PQ_PairArray getI2PQ();
    
    virtual void saveL(const TlMatrix& L);
    virtual TlMatrix getL();

    TlSymmetricMatrix getCholeskyVector(const TlVector& L_col,
                                        const I2PQ_Type& I2PQ);

    virtual TlSymmetricMatrix getPMatrix();

    virtual void divideCholeskyBasis(const index_type numOfCBs,
                                     index_type *pStart, index_type *pEnd);

    size_type pqPairIndex(const PQ_Pair& pq) const {
        // 'U' format
        assert(pq.shellIndex1 <= pq.shellIndex2);
        size_type answer = pq.shellIndex1 +  pq.shellIndex2 * (pq.shellIndex2 +1) / 2;
        return answer;
    }

protected: // for exact
    struct ShellPair {
    public:
        ShellPair(index_type index1 =0, index_type index2 =0) : shellIndex1(index1), shellIndex2(index2) {
        }
        
    public:
        index_type shellIndex1;
        index_type shellIndex2;
    };
    typedef std::vector<index_type> ShellArray;
    typedef std::vector<ShellArray> ShellArrayTable;
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
