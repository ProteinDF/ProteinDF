#ifndef DFCD_H
#define DFCD_H

#include <deque>
#include "DfObject.h"
#include "DfTaskCtrl.h"
#include "TlOrbitalInfo.h"
#include "TlSparseSymmetricMatrix.h"
#include "TlSymmetricMatrix.h"
#include "TlStlUtils.h"

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
    void getJ_S(TlSymmetricMatrix *pJ);
    void getK_S(const RUN_TYPE runType,
                TlSymmetricMatrix *pK);

    bool useSymmetric_;

protected:
    class Index2 {
    public:
        explicit Index2(index_type i1 =0, index_type i2 =0) : index1_(i1), index2_(i2) {
        }

        bool operator<(const Index2& rhs) const {
            if (this->index1_ < rhs.index1_) {
                return true;
            } else if (this->index1_ == rhs.index1_) {
                return (this->index2_ < rhs.index2_);
            }

            return false;
        }

        bool operator==(const Index2& rhs) const {
            return ((this->index1_ == rhs.index1_) &&
                    (this->index2_ == rhs.index2_));
        }

        index_type index1() const {
            return this->index1_;
        }

        index_type index2() const {
            return this->index2_;
        }

    protected:
        index_type index1_;
        index_type index2_;
    };

    class Index4 {
    public:
        explicit Index4(index_type i1 =0, index_type i2 =0,
                        index_type i3 =0, index_type i4 =0) 
            : index2_1_(i1, i2), index2_2_(i3, i4) {
        }

        bool operator<(const Index4& rhs) const {
            if (this->index2_1_ < rhs.index2_1_) {
                return true;
            } else if (this->index2_1_ == rhs.index2_1_) {
                return (this->index2_2_ < rhs.index2_2_);
            }

            return false;
        }

        bool operator==(const Index4& rhs) const {
            return ((this->index2_1_ == rhs.index2_1_) &&
                    (this->index2_2_ == rhs.index2_2_));
        }

        index_type index1() const {
            return this->index2_1_.index1();
        }

        index_type index2() const {
            return this->index2_1_.index2();
        }

        index_type index3() const {
            return this->index2_2_.index1();
        }

        index_type index4() const {
            return this->index2_2_.index2();
        }

    private:
        Index2 index2_1_;
        Index2 index2_2_;
    };


protected:
    typedef std::vector<Index2> PQ_PairArray;
    typedef std::vector<size_type> PQ2I_Type;

protected:
    void createEngines();
    void destroyEngines();
    
protected:
    virtual DfTaskCtrl* getDfTaskCtrlObject() const;

    virtual void finalize(TlSymmetricMatrix *pMat);
    virtual void finalize(TlSparseSymmetricMatrix *pMat);
    virtual void finalize(TlMatrix *pMat);
    virtual void finalize(TlSparseMatrix *pMat);
    virtual void finalize_I2PQ(PQ_PairArray *pI2PQ);

protected:
    virtual void saveI2PQ(const PQ_PairArray& I2PQ);
    PQ_PairArray getI2PQ();

    virtual void saveL(const TlMatrix& L);
    virtual TlMatrix getL();

    TlSymmetricMatrix getCholeskyVector(const TlVector& L_col,
                                        const PQ_PairArray& I2PQ);

    virtual TlSymmetricMatrix getPMatrix();

    virtual void divideCholeskyBasis(const index_type numOfCBs,
                                     index_type *pStart, index_type *pEnd);

protected:
    virtual void calcCholeskyVectors_onTheFly();

    void calcDiagonals(TlSparseSymmetricMatrix *pSchwartzTable,
                       PQ_PairArray *pI2PQ,
                       TlVector *pDiagonals);
    void calcDiagonals_kernel(const std::vector<DfTaskCtrl::Task2>& taskList,
                              TlSparseSymmetricMatrix *pSchwartzTable,
                              TlSparseSymmetricMatrix *pDiagonalMat,
                              PQ_PairArray *pI2PQ);

    bool isAliveBySchwartzCutoff(const index_type shellIndexP,
                                 const index_type shellIndexQ,
                                 const index_type shellIndexR,
                                 const index_type shellIndexS,
                                 const int shellQuartetType,
                                 const TlSparseMatrix& schwarzTable,
                                 const double threshold);
    void initializeCutoffStats();
    void schwartzCutoffReport();
    mutable std::vector<unsigned long> cutoffAll_schwartz_;
    mutable std::vector<unsigned long> cutoffAlive_schwartz_;

protected:
    // 2電子積分キャッシュ -----------------------------------------------------
    class IndexPair2S : public Index2 {
     public:
         IndexPair2S(index_type i1, index_type i2) 
             : Index2(i1, i2) {
             if (this->index1_ < this->index2_) {
                 std::swap(this->index1_, this->index2_);
             }
         }
    };

    class IndexPair4S {
    public:
        IndexPair4S(index_type i1, index_type i2, index_type i3, index_type i4) 
            : ip2_1_(i1, i2), ip2_2_(i3, i4) {
            if (this->ip2_1_ < this->ip2_2_) {
                std::swap(this->ip2_1_, this->ip2_2_);
            }
        }

        IndexPair4S(const Index4& i4) 
            : ip2_1_(i4.index1(), i4.index2()), ip2_2_(i4.index3(), i4.index4()) {
            if (this->ip2_1_ < this->ip2_2_) {
                std::swap(this->ip2_1_, this->ip2_2_);
            }
        }

        bool operator<(const IndexPair4S& rhs) const {
            if (this->ip2_1_ < rhs.ip2_1_) {
                return true;
            } else if (this->ip2_1_ == rhs.ip2_1_) {
                return (this->ip2_2_ < rhs.ip2_2_);
            }

            return false;
        }

        bool operator==(const IndexPair4S& rhs) const {
            return ((this->ip2_1_ == rhs.ip2_1_) &&
                    (this->ip2_2_ == rhs.ip2_2_));
        }

        index_type index1() const {
            return this->ip2_1_.index1();
        }

        index_type index2() const {
            return this->ip2_1_.index2();
        }

        index_type index3() const {
            return this->ip2_2_.index1();
        }

        index_type index4() const {
            return this->ip2_2_.index2();
        }

    private:
        IndexPair2S ip2_1_;
        IndexPair2S ip2_2_;
    };

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
    typedef std::map<IndexPair4S, std::vector<double> > ERI_CacheType;

    /// 与えられたsuper matrix の要素に対し、2電子積分を計算して代入する。
    /// On-the-Fly時に使用する。
    virtual std::vector<double> getSuperMatrixElements(const index_type G_row,
                                                       const std::vector<index_type>& G_col_list,
                                                       const PQ_PairArray& I2PQ,
                                                       const TlSparseSymmetricMatrix& schwartzTable);

    /// 要求されたsuper matrixの行列要素のうち、必要なshell indexのリストを返す。
    ///
    /// @param G_row 必要なsuper matrixの行要素。
    /// @param G_col_list 必要なsuper matrixの列要素の配列。
    std::vector<DfCD::Index4> getCalcList(const index_type G_row,
                                          const std::vector<index_type>& G_col_list,
                                          const PQ_PairArray& I2PQ);

    /// 計算リストの2電子積分を求め、キャッシュに代入して返す。
    void calcERIs(const std::vector<Index4>& calcList,
                  const TlSparseSymmetricMatrix& schwartzTable);

    /// キャッシュから必要な行列要素を代入する。
    std::vector<double> setERIs(const index_type G_row,
                                const std::vector<index_type> G_col_list,
                                const PQ_PairArray& I2PQ);

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

// =====================================================================
protected:
    /// 非対称形2index管理クラス
    ///
    /// (pq)の組を管理するクラス
    /// pの基底関数とqの基底関数が異なる場合に使用する
    class IndexPair2A {
    public:
        explicit IndexPair2A(index_type i1 =0, index_type i2 =0) : index1_(i1), index2_(i2) {
        }

        bool operator<(const IndexPair2A& rhs) const {
            if (this->index1_ < rhs.index1_) {
                return true;
            } else if (this->index1_ == rhs.index1_) {
                return (this->index2_ < rhs.index2_);
            }

            return false;
        }

        bool operator==(const IndexPair2A& rhs) const {
            return ((this->index1_ == rhs.index1_) &&
                    (this->index2_ == rhs.index2_));
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
    typedef std::vector<IndexPair2A> PQ_PairArray_A;

    /// 非対称形4index管理クラス
    ///
    /// (pq, rs)の組を管理するクラス
    /// p,rの基底関数とq, sの基底関数が異なる場合に使用する
    struct IndexPair4A {
    public:
        explicit IndexPair4A(index_type i1 =0, index_type i2 =0,
                             index_type i3 =0, index_type i4 =0) 
            : indexPair1_(i1, i2), indexPair2_(i3, i4) {
        }

        bool operator<(const IndexPair4A& rhs) const {
            if (this->indexPair1_ < rhs.indexPair1_) {
                return true;
            } else if (this->indexPair1_ == rhs.indexPair1_) {
                return (this->indexPair2_ < rhs.indexPair2_);
            }

            return false;
        }

        bool operator==(const IndexPair4A& rhs) const {
            return ((this->indexPair1_ == rhs.indexPair1_) &&
                    (this->indexPair2_ == rhs.indexPair2_));
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
        IndexPair2A indexPair1_;
        IndexPair2A indexPair2_;
    };

    /// 2電子積分キャッシュの型
    typedef std::map<IndexPair4A, std::vector<double> > ERI_CacheType_A;
    ERI_CacheType_A ERI_cache_A_;

    virtual void calcCholeskyVectorsA_onTheFly(const TlOrbitalInfoObject& orbInfo_p,
                                               const TlOrbitalInfoObject& orbInfo_q);
    void calcDiagonalsA(const TlOrbitalInfoObject& orbInfo_p,
                        const TlOrbitalInfoObject& orbInfo_q,
                        PQ_PairArray_A *pI2PQ,
                        TlSparseMatrix *pSchwartzTable,
                        TlVector *pDiagonals);
    void calcDiagonalsA_kernel(const TlOrbitalInfoObject& orbInfo_p,
                               const TlOrbitalInfoObject& orbInfo_q,
                               const std::vector<DfTaskCtrl::Task2>& taskList,
                               PQ_PairArray_A *pI2PQ,
                               TlSparseMatrix *pSchwartzTable,
                               TlSparseMatrix *pDiagonalMat);

    void finalizeI2PQ_A(PQ_PairArray_A *pI2PQ);
    virtual void saveI2PQ_A(const PQ_PairArray_A& I2PQ);

    virtual std::vector<double>
    getSuperMatrixElementsA(const TlOrbitalInfoObject& orbInfo_p,
                            const TlOrbitalInfoObject& orbInfo_q,
                            const index_type G_row,
                            const std::vector<index_type>& G_col_list,
                            const PQ_PairArray_A& I2PQ,
                            const TlSparseMatrix& schwartzTable);
    std::vector<DfCD::IndexPair4A> 
    getCalcListA(const TlOrbitalInfoObject& orbInfo_p,
                 const TlOrbitalInfoObject& orbInfo_q,
                 const index_type G_row,
                 const std::vector<index_type>& G_col_list,
                 const PQ_PairArray_A& I2PQ);
    void calcERIsA(const TlOrbitalInfoObject& orbInfo_p,
                   const TlOrbitalInfoObject& orbInfo_q,
                   const std::vector<IndexPair4A>& calcList,
                   const TlSparseMatrix& schwartzTable);
    std::vector<double>
    setERIsA(const TlOrbitalInfoObject& orbInfo_p,
             const TlOrbitalInfoObject& orbInfo_q,
             const index_type G_row,
             const std::vector<index_type> G_col_list,
             const PQ_PairArray_A& I2PQ);
    PQ_PairArray_A getI2PQ_A();
    TlMatrix getCholeskyVectorA(const TlOrbitalInfoObject& orbInfo_p,
                                const TlOrbitalInfoObject& orbInfo_q,
                                const TlVector& L_col,
                                const PQ_PairArray_A& I2PQ);

    /// デバッグ用にSuperMatrixを作成します
    ///
    /// V_pq,rs = (pq|rs) or (pqrs) 
    /// @param[in] orbInfo_p pまたはrで示される軌道の軌道情報オブジェクト
    /// @param[in] orbInfo_q qまたはsで示される軌道の軌道情報オブジェクト
    /// @retval supermatrix
    TlSymmetricMatrix getSuperMatrix(const TlOrbitalInfoObject& orbInfo_p,
                                     const TlOrbitalInfoObject& orbInfo_q);

    /// コレスキー分解を行います(デバッグ用)
    ///
    /// V = L * L
    /// @param[in] V supermatrix
    /// @retval コレスキーベクトル(L)
    TlMatrix calcCholeskyVectors(const TlSymmetricMatrix& V);

public:
    void getJ_A(TlSymmetricMatrix* pJ);
    void getK_A(const RUN_TYPE runType,
                TlSymmetricMatrix *pK);
};


#endif // DFCD_H
