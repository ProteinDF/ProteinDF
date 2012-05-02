#ifndef DFCD_H
#define DFCD_H

#include "DfObject.h"
#include "DfTaskCtrl.h"
#include "TlOrbitalInfo.h"
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
    virtual void calcCholeskyVectors();

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

        bool operator<(const PQ_Pair& rhs) const {
            bool answer = false;
            if ((this->shellIndex1 < rhs.shellIndex1) ||
                ((this->shellIndex1 == rhs.shellIndex1) && (this->shellIndex2 < rhs.shellIndex2))) {
                answer = true;
            }
            return answer;
        }

    public:
        index_type shellIndex1;
        index_type shellIndex2;
    };
    
    typedef std::vector<PQ_Pair> PQ_PairArray;
    typedef std::vector<PQ_Pair> I2PQ_Type;
    typedef std::vector<size_type> PQ2I_Type;

protected:
    void createEngines();
    void destroyEngines();
    
protected:
    void makeSuperMatrix_screening();
    TlSparseSymmetricMatrix makeSchwarzTable(const TlOrbitalInfoObject& orbitalInfo);

    std::size_t index(index_type p, index_type q) const;

    virtual DfTaskCtrl* getDfTaskCtrlObject() const;

    virtual void finalize(TlSymmetricMatrix *pMat);
    virtual void finalize(TlSparseSymmetricMatrix *pMat);
    virtual void finalize_I2PQ(I2PQ_Type *pI2PQ);

protected:
    void calcPQPQ(const TlOrbitalInfoObject& orbitalInfo,
                  TlSparseSymmetricMatrix *pSchwarzTable,
                  PQ_PairArray *pI2PQ);
    void calcPQPQ_kernel(const TlOrbitalInfoObject& orbitalInfo,
                         const std::vector<DfTaskCtrl::Task2>& taskList,
                         TlSparseSymmetricMatrix *pSchwarzTable,
                         PQ_PairArray *pI2PQ);

    TlSymmetricMatrix getGMatrix(const TlOrbitalInfoObject& orbitalInfo, 
                                 const TlSparseSymmetricMatrix& schwarzTable,
                                 const index_type numOfItilde,
                                 const PQ2I_Type& PQ2I);

    void makeSuperMatrix_kernel2(const TlOrbitalInfoObject& orbitalInfo,
                                 const std::vector<DfTaskCtrl::Task4>& taskList,
                                 const PQ2I_Type& PQ2I,
                                 TlMatrixObject* pG);
    void storeG2(const index_type shellIndexP, const int maxStepsP,
                 const index_type shellIndexQ, const int maxStepsQ,
                 const index_type shellIndexR, const int maxStepsR,
                 const index_type shellIndexS, const int maxStepsS,
                 const PQ2I_Type& PQ2I,
                 const DfEriEngine& engine,
                 TlMatrixObject* pG);
    virtual void saveI2PQ(const I2PQ_Type& I2PQ);
    PQ_PairArray getI2PQ();

    void makeL(const TlSymmetricMatrix& G);
    
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

protected:
    // NEW ---------------------------------------------------------------------
    void calcCholeskyVectors_onTheFly();
    void calcDiagonals(TlSparseSymmetricMatrix *pSchwartzTable,
                       PQ_PairArray *pI2PQ,
                       TlVector *pDiagonals);
    void calcDiagonals_kernel(const std::vector<DfTaskCtrl::Task2>& taskList,
                              TlSparseSymmetricMatrix *pSchwartzTable,
                              TlSparseSymmetricMatrix *pDiagonalMat,
                              PQ_PairArray *pI2PQ);
    void calcERIs(const TlSparseSymmetricMatrix& schwartzTable,
                  const I2PQ_Type& I2PQ,
                  TlSparseSymmetricMatrix* pG);
    bool isAliveBySchwartzCutoff(const index_type shellIndexP,
                                 const index_type shellIndexQ,
                                 const index_type shellIndexR,
                                 const index_type shellIndexS,
                                 const int shellQuartetType,
                                 const TlSparseSymmetricMatrix& schwarzTable,
                                 const double threshold);
    void schwartzCutoffReport();
    mutable std::vector<unsigned long> cutoffAll_schwartz_;
    mutable std::vector<unsigned long> cutoffAlive_schwartz_;
    
protected:
    index_type numOfPQs_;

    DfEriEngine* pEriEngines_;
    TlOrbitalInfo orbitalInfo_;

    double cutoffThreshold_;
    double cutoffEpsilon3_;
    double CDAM_tau_;
    double epsilon_;

private:
#ifdef CHECK_LOOP
    TlSymmetricMatrix check;
#endif // CHECK_LOOP
};


#endif // DFCD_H
