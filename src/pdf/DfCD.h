#ifndef DFCD_H
#define DFCD_H

#include <deque>
#include "DfObject.h"
#include "DfTaskCtrl.h"
#include "TlOrbitalInfo.h"
#include "TlSparseSymmetricMatrix.h"
#include "TlSymmetricMatrix.h"
#include "TlStlUtils.h"

// #define CD_DEBUG

class TlOrbitalInfo;
class DfEriEngine;

// #define CHECK_LOOP // 計算ループ構造のチェック

class DfCD : public DfObject
{
public:
    DfCD(TlSerializeData* pPdfParam);
    virtual ~DfCD();

public:
    void calcCholeskyVectors();

    void getJ(TlSymmetricMatrix *pJ);
    void getK(const RUN_TYPE runType,
              TlSymmetricMatrix *pK);

protected:
    struct IndexPair2 {
    public:
        explicit IndexPair2(index_type i1 =0, index_type i2 =0) : index1_(i1), index2_(i2) {
            if (this->index1_ > this->index2_) {
                std::swap(this->index1_, this->index2_);
            }
        }

        std::size_t index() const {
            // 'U' format
            assert(this->index1_ <= this->index2_);
            return this->index1_ + this->index2_ * (this->index2_ +1) / 2;
        }

        bool operator<(const IndexPair2& rhs) const {
            return (this->index() < rhs.index());
        }

        index_type index1() const {
            return this->index1_;
        }

        index_type index2() const {
            return this->index2_;
        }

    private:
        index_type index1_;
        index_type index2_;
    };

    struct IndexPair4 {
    public:
        explicit IndexPair4(index_type i1 =0, index_type i2 =0,
                            index_type i3 =0, index_type i4 =0) 
            : indexPair1_(i1, i2), indexPair2_(i3, i4) {
            if (this->indexPair1_.index() > this->indexPair2_.index()) {
                std::swap(this->indexPair1_, this->indexPair2_);
            }
        }

        std::size_t index() const {
            std::size_t pair_index1 = this->indexPair1_.index();
            std::size_t pair_index2 = this->indexPair2_.index();
            assert(pair_index1 <= pair_index2);
            return pair_index1 + pair_index2 * (pair_index2 +1) / 2;
        }

        bool operator<(const IndexPair4& rhs) const {
            return (this->index() < rhs.index());
        }

        bool operator==(const IndexPair4& rhs) const {
            return (this->index() == rhs.index());
        }

        index_type index1() const {
            return this->indexPair1_.index1();
        }

        index_type index2() const {
            return this->indexPair1_.index2();
        }

        index_type index3() const {
            return this->indexPair2_.index1();
        }

        index_type index4() const {
            return this->indexPair2_.index2();
        }

    private:
        IndexPair2 indexPair1_;
        IndexPair2 indexPair2_;
    };


protected:
    typedef std::vector<IndexPair2> PQ_PairArray;
    typedef std::vector<IndexPair2> I2PQ_Type;
    typedef std::vector<size_type> PQ2I_Type;

protected:
    void createEngines();
    void destroyEngines();
    
protected:
    void makeSuperMatrix_screening();
    TlSparseSymmetricMatrix makeSchwarzTable(const TlOrbitalInfoObject& orbitalInfo);

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

protected:
    // NEW ---------------------------------------------------------------------
    virtual void calcCholeskyVectors_onTheFly();
    void calcDiagonals(TlSparseSymmetricMatrix *pSchwartzTable,
                       PQ_PairArray *pI2PQ,
                       TlVector *pDiagonals);
    void calcDiagonals_kernel(const std::vector<DfTaskCtrl::Task2>& taskList,
                              TlSparseSymmetricMatrix *pSchwartzTable,
                              TlSparseSymmetricMatrix *pDiagonalMat,
                              PQ_PairArray *pI2PQ);

    // void calcERIs(const TlSparseSymmetricMatrix& schwartzTable,
    //               const I2PQ_Type& I2PQ,
    //               TlSparseSymmetricMatrix* pG);
    bool isAliveBySchwartzCutoff(const index_type shellIndexP,
                                 const index_type shellIndexQ,
                                 const index_type shellIndexR,
                                 const index_type shellIndexS,
                                 const int shellQuartetType,
                                 const TlSparseSymmetricMatrix& schwarzTable,
                                 const double threshold);
    void initializeCutoffStats();
    void schwartzCutoffReport();
    mutable std::vector<unsigned long> cutoffAll_schwartz_;
    mutable std::vector<unsigned long> cutoffAlive_schwartz_;

protected:
    // 2電子積分キャッシュ -----------------------------------------------------
    // class ERI_CacheManager {
    // public:
    //     typedef std::map<IndexPair4, std::vector<double> > ERI_CacheType;

    // private:
    //     typedef std::deque<ERI_CacheType> ERI_CacheGenerations;

    // public:
    //     ERI_CacheManager() {
    //         this->eriCaches_.clear();
    //         this->createNewGeneration();
    //     }

    // public:
    //     void createNewGeneration() {
    //         this->eriCaches_.push_front(ERI_CacheType());
    //     }

    //     bool find(const IndexPair4& key) const {
    //         bool answer = false;

    //         ERI_CacheGenerations::iterator itEnd = this->eriCaches_.end();
    //         for (ERI_CacheGenerations::iterator it = this->eriCaches_.begin(); it != itEnd; ++it) {
    //             ERI_CacheType::iterator genIt = it->find(key);
    //             if (genIt != it->end()) {
    //                 if (it != this->eriCaches_.begin()) {
    //                     TlStlUtils::efficientAddOrUpdate(this->eriCaches_.front(),
    //                                                      genIt->first, genIt->second);
    //                     it->erase(genIt);
    //                 }
    //                 answer = true;
    //                 break;
    //             }
    //         }
    //         return answer;
    //     }

    //     void insert(const ERI_CacheType& caches) {
    //         this->eriCaches_.front().insert(caches.begin(), caches.end());
    //     }

    //     std::vector<double> get(const IndexPair4& key) const {
    //         std::vector<double> answer;
    //         ERI_CacheGenerations::iterator itEnd = this->eriCaches_.end();
    //         for (ERI_CacheGenerations::iterator it = this->eriCaches_.begin(); it != itEnd; ++it) {
    //             ERI_CacheType::iterator genIt = it->find(key);
    //             if (genIt != it->end()) {
    //                 answer = genIt->second;
    //                 break;
    //             }
    //         }
    //         return answer;
    //     }

    //     std::size_t memSize() const {
    //         std::size_t answer = 0;
    //         ERI_CacheGenerations::iterator itEnd = this->eriCaches_.end();
    //         for (ERI_CacheGenerations::iterator it = this->eriCaches_.begin(); it != itEnd; ++it) {
    //             answer += (sizeof(IndexPair4) + sizeof(double)) * it->size();
    //         }
    //         return answer;
    //     }

    //     void releaseMem(const std::size_t limit) {
    //         while ((this->eriCaches_.size() > 1) &&
    //                (this->memSize() > limit)) {
    //             this->eriCaches_.pop_back();
    //         }
    //     }

    // private:
    //     mutable ERI_CacheGenerations eriCaches_;
    // };

    /// 2電子積分キャッシュの型
    typedef std::map<IndexPair4, std::vector<double> > ERI_CacheType;

    /// 与えられたsuper matrix の要素に対し、2電子積分を計算して代入する。
    /// On-the-Fly時に使用する。
    virtual std::vector<double> getSuperMatrixElements(const index_type G_row,
                                                       const std::vector<index_type>& G_col_list,
                                                       const I2PQ_Type& I2PQ,
                                                       const TlSparseSymmetricMatrix& schwartzTable);

    /// 要求されたsuper matrixの行列要素のうち、必要なshell indexのリストを返す。
    ///
    /// @param G_row 必要なsuper matrixの行要素。
    /// @param G_col_list 必要なsuper matrixの列要素の配列。
    std::vector<DfCD::IndexPair4> getCalcList(const index_type G_row,
                                              const std::vector<index_type>& G_col_list,
                                              const I2PQ_Type& I2PQ);

    /// 計算リストの2電子積分を求め、キャッシュに代入して返す。
    void calcERIs(const std::vector<IndexPair4>& calcList,
                  const TlSparseSymmetricMatrix& schwartzTable);

    /// キャッシュから必要な行列要素を代入する。
    std::vector<double> setERIs(const index_type G_row,
                                const std::vector<index_type> G_col_list,
                                const I2PQ_Type& I2PQ);

    /// 2電子積分をキャッシュするかどうか
    bool isStoreERIs_;

    /// 2電子積分キャッシュ
    // ERI_CacheManager ERI_cache_manager_;
    // TlCache<IndexPair4, std::vector<double> > ERI_cache_manager_;
    ERI_CacheType ERI_cache_;


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
