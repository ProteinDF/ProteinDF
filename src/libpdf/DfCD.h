// Copyright (C) 2002-2014 The ProteinDF project
// see also AUTHORS and README.
// 
// This file is part of ProteinDF.
// 
// ProteinDF is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
// ProteinDF is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with ProteinDF.  If not, see <http://www.gnu.org/licenses/>.

#ifndef DFCD_H
#define DFCD_H

#include <deque>
#include "DfObject.h"
#include "DfTaskCtrl.h"
#include "TlOrbitalInfo.h"
#include "TlSparseSymmetricMatrix.h"
#include "TlSymmetricMatrix.h"
#include "TlStlUtils.h"
#include "TlRowVectorMatrix.h"
#include "TlColVectorMatrix.h"

class TlOrbitalInfo;
class DfEngineObject;

// #define CHECK_LOOP // 計算ループ構造のチェック

class DfCD : public DfObject
{
public:
    DfCD(TlSerializeData* pPdfParam);
    virtual ~DfCD();

public:
    virtual void calcCholeskyVectorsForJK();
    virtual void calcCholeskyVectorsForK();
    virtual void calcCholeskyVectorsForGridFree();

    void getJ(TlSymmetricMatrix *pJ);
    virtual void getK(const RUN_TYPE runType,
                      TlSymmetricMatrix *pK);
    virtual void getM(const TlSymmetricMatrix& P,
                      TlSymmetricMatrix* pM);

protected:
    enum FastCDK_MODE {
        FASTCDK_NONE,
        FASTCDK_PRODUCTIVE,
        FASTCDK_PRODUCTIVE_FULL,
        FASTCDK_DEBUG_FULL_SUPERMATRIX, // V_pr,qs = <pq|rs>
        FASTCDK_DEBUG_SUPERMATRIX       // V_pr,qs = <pq|rs> + <ps|rq>
    };
    
protected:
    void getJ_S(TlSymmetricMatrix *pJ);
    void getJ_S_mmap(TlSymmetricMatrix* pJ);

    virtual TlSymmetricMatrix getPMatrix(const RUN_TYPE runType, const int iteration);
    
    virtual void getK_S_woCD(const RUN_TYPE runType,
                             TlSymmetricMatrix *pK);
    virtual void getK_S_woCD_mmap(const RUN_TYPE runType,
                                  TlSymmetricMatrix *pK);

    virtual void getK_S_fast(const RUN_TYPE runType,
                             TlSymmetricMatrix *pK);

    virtual void getM_S(const TlSymmetricMatrix& P,
                        TlSymmetricMatrix* pM);
    virtual void getM_A(const TlSymmetricMatrix& P,
                        TlSymmetricMatrix* pM);

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
    template<class EngineClass>
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
    virtual void saveI2PQ(const PQ_PairArray& I2PQ, const std::string& filepath);
    virtual PQ_PairArray getI2PQ(const std::string& filepath);

    virtual void saveLjk(const TlRowVectorMatrix& Ljk);
    virtual void saveLk(const TlRowVectorMatrix& Lk);
    virtual void debugOutLjk(const TlMatrix& Ljk);
    virtual void debugOutLk(const TlMatrix& Lk);
    virtual TlColVectorMatrix getLjk();
    virtual TlColVectorMatrix getLk();

    virtual void saveLxc(const TlRowVectorMatrix& Ljk);
    virtual void debugOutLxc(const TlMatrix& Lxc);
    virtual TlColVectorMatrix getLxc();

    TlSymmetricMatrix getCholeskyVector(const TlVector& L_col,
                                        const PQ_PairArray& I2PQ);

    virtual TlSymmetricMatrix getPMatrix();

    virtual void divideCholeskyBasis(const index_type numOfCBs,
                                     index_type *pStart, index_type *pEnd);

protected:
    template<class EngineClass>
    TlRowVectorMatrix calcCholeskyVectorsOnTheFly(const TlOrbitalInfoObject& orbInfo_p,
                                                   const TlOrbitalInfoObject& orbInfo_q,
                                                   const std::string& I2PQ_path);

    virtual TlRowVectorMatrix calcCholeskyVectorsOnTheFlyA(const TlOrbitalInfoObject& orbInfo_p,
                                                            const TlOrbitalInfoObject& orbInfo_q,
                                                            const std::string& I2PQ_path);

public:
    void calcDiagonals(const TlOrbitalInfoObject& orbInfo,
                       PQ_PairArray *pI2PQ,
                       TlVector *pDiagonals);
    void calcDiagonals_K_full(const TlOrbitalInfoObject& orbInfo,
                              PQ_PairArray *pI2PR,
                              TlVector *pDiagonals);
    void calcDiagonals_K_half(const TlOrbitalInfoObject& orbInfo,
                              PQ_PairArray *pI2PR,
                              TlVector *pDiagonals);
protected:
    void calcDiagonals_kernel(const TlOrbitalInfoObject& orbInfo,
                              const std::vector<DfTaskCtrl::Task2>& taskList,
                              TlSparseSymmetricMatrix *pDiagonalMat,
                              PQ_PairArray *pI2PQ);
    void calcDiagonals_K_full_kernel(const TlOrbitalInfoObject& orbInfo,
                                     const std::vector<DfTaskCtrl::Task2>& taskList,
                                     TlSparseMatrix *pDiagonalMat,
                                     PQ_PairArray *pI2PR);
    void calcDiagonals_K_half_kernel(const TlOrbitalInfoObject& orbInfo,
                                     const std::vector<DfTaskCtrl::Task2>& taskList,
                                     TlSparseSymmetricMatrix *pDiagonalMat,
                                     PQ_PairArray *pI2PR);

    bool isAliveBySchwartzCutoff(const index_type shellIndexP,
                                 const index_type shellIndexQ,
                                 const index_type shellIndexR,
                                 const index_type shellIndexS,
                                 const int shellQuartetType,
                                 const TlSparseMatrix& schwarzTable,
                                 const double threshold);
    void initializeCutoffStats(const int maxShellType);
    void schwartzCutoffReport(const int maxShellType);
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
        explicit IndexPair4S(index_type i1 =0, index_type i2 =0, index_type i3 =0, index_type i4 =0) 
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

    class IndexPair4A {
    public:
        explicit IndexPair4A(index_type i1 =0, index_type i2 =0, index_type i3 =0, index_type i4 =0) 
            : index2_1_(i1, i2), index2_2_(i3, i4) {
            // if (this->index2_1_ < this->index2_2_) {
            //     std::swap(this->index2_1_, this->index2_2_);
            // }
        }

        IndexPair4A(const Index4& i4) 
            : index2_1_(i4.index1(), i4.index2()), index2_2_(i4.index3(), i4.index4()) {
            // if (this->index2_1_ < this->index2_2_) {
            //     std::swap(this->index2_1_, this->index2_2_);
            // }
        }

        bool operator<(const IndexPair4A& rhs) const {
            if (this->index2_1_ < rhs.index2_1_) {
                return true;
            } else if (this->index2_1_ == rhs.index2_1_) {
                return (this->index2_2_ < rhs.index2_2_);
            }

            return false;
        }

        bool operator==(const IndexPair4A& rhs) const {
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

public:
    /// 与えられたsuper matrix の要素に対し、2電子積分を計算して代入する。
    /// On-the-Fly時に使用する。
    virtual std::vector<double> getSuperMatrixElements(const TlOrbitalInfoObject& orbInfo,
                                                       const index_type G_row,
                                                       const std::vector<index_type>& G_col_list,
                                                       const PQ_PairArray& I2PQ);

protected:
    /// 要求されたsuper matrixの行列要素のうち、必要なshell indexのリストを返す。
    ///
    /// @param G_row 必要なsuper matrixの行要素。
    /// @param G_col_list 必要なsuper matrixの列要素の配列。
    virtual std::vector<IndexPair4S> 
    getCalcList(const TlOrbitalInfoObject& orbInfo,
                const index_type G_row,
                const std::vector<index_type>& G_col_list,
                const index_type start,
                const index_type end,
                const PQ_PairArray& I2PQ);

    /// 計算リストの2電子積分を求め、キャッシュに代入して返す。
    void calcERIs(const TlOrbitalInfoObject& orbInfo,
                  const std::vector<IndexPair4S>& calcList);

    /// キャッシュから必要な行列要素を代入する。
    std::vector<double> setERIs(const TlOrbitalInfoObject& orbInfo,
                                const index_type G_row,
                                const std::vector<index_type> G_col_list,
                                const index_type start,
                                const index_type end,
                                const PQ_PairArray& I2PQ);

    /// 2電子積分をキャッシュするかどうか
    // bool isStoreERIs_;

    /// 2電子積分キャッシュ
    // ERI_CacheManager ERI_cache_manager_;
    // TlCache<IndexPair4, std::vector<double> > ERI_cache_manager_;
    ERI_CacheType ERI_cache_;

protected:
    void getJ_S_v2(TlSymmetricMatrix* pJ);
    TlVector getScreenedDensityMatrix(const PQ_PairArray& I2PQ);
    TlVector getScreenedDensityMatrix(const RUN_TYPE runTyoe, const PQ_PairArray& I2PQ);
    // TlVector getScreenedDensityMatrix2(const RUN_TYPE runTyoe, const PQ_PairArray& I2PQ);
    void expandJMatrix(const TlVector& vJ, const PQ_PairArray& I2PQ, TlSymmetricMatrix* pJ);
    void expandKMatrix(const TlVector& vK, const PQ_PairArray& I2PR, TlSymmetricMatrix* pK);

protected:
    //index_type numOfPQs_;

    DfEngineObject** pEngines_;

    double cutoffThreshold_;
    double cutoffThreshold_primitive_;
    double CDAM_tau_;
    double epsilon_;
    double CDAM_tau_K_;
    double epsilon_K_;

    FastCDK_MODE fastCDK_mode_;

    bool useMmapMatrix_;

// =====================================================================
protected:
    /// 非対称形2index管理クラス
    ///
    /// (pq)の組を管理するクラス
    /// pの基底関数とqの基底関数が異なる場合に使用する
    // class IndexPair2A {
    // public:
    //     explicit IndexPair2A(index_type i1 =0, index_type i2 =0) : index1_(i1), index2_(i2) {
    //     }

    //     bool operator<(const IndexPair2A& rhs) const {
    //         if (this->index1_ < rhs.index1_) {
    //             return true;
    //         } else if (this->index1_ == rhs.index1_) {
    //             return (this->index2_ < rhs.index2_);
    //         }

    //         return false;
    //     }

    //     bool operator==(const IndexPair2A& rhs) const {
    //         return ((this->index1_ == rhs.index1_) &&
    //                 (this->index2_ == rhs.index2_));
    //     }

    //     index_type index1() const {
    //         return this->index1_;
    //     }

    //     index_type index2() const {
    //         return this->index2_;
    //     }

    // private:
    //     index_type index1_;
    //     index_type index2_;
    // };
    // typedef std::vector<Index2> PQ_PairArray_A;

    /// 非対称形4index管理クラス
    ///
    /// (pq, rs)の組を管理するクラス
    /// p,rの基底関数とq, sの基底関数が異なる場合に使用する
    // struct IndexPair4A {
    // public:
    //     explicit IndexPair4A(index_type i1 =0, index_type i2 =0,
    //                          index_type i3 =0, index_type i4 =0) 
    //         : indexPair1_(i1, i2), indexPair2_(i3, i4) {
    //     }

    //     bool operator<(const IndexPair4A& rhs) const {
    //         if (this->indexPair1_ < rhs.indexPair1_) {
    //             return true;
    //         } else if (this->indexPair1_ == rhs.indexPair1_) {
    //             return (this->indexPair2_ < rhs.indexPair2_);
    //         }

    //         return false;
    //     }

    //     bool operator==(const IndexPair4A& rhs) const {
    //         return ((this->indexPair1_ == rhs.indexPair1_) &&
    //                 (this->indexPair2_ == rhs.indexPair2_));
    //     }

    //     index_type index1() const {
    //         return this->indexPair1_.index1();
    //     }

    //     index_type index2() const {
    //         return this->indexPair1_.index2();
    //     }

    //     index_type index3() const {
    //         return this->indexPair2_.index1();
    //     }

    //     index_type index4() const {
    //         return this->indexPair2_.index2();
    //     }

    // private:
    //     IndexPair2A indexPair1_;
    //     IndexPair2A indexPair2_;
    // };

    /// 2電子積分キャッシュの型
    typedef std::map<IndexPair4A, std::vector<double> > ERI_CacheType_A;
    ERI_CacheType_A ERI_cache_A_;

    void calcDiagonalsA(const TlOrbitalInfoObject& orbInfo_p,
                        const TlOrbitalInfoObject& orbInfo_q,
                        PQ_PairArray *pI2PQ,
                        TlSparseMatrix *pSchwartzTable,
                        TlVector *pDiagonals);
    void calcDiagonalsA_kernel(const TlOrbitalInfoObject& orbInfo_p,
                               const TlOrbitalInfoObject& orbInfo_q,
                               const std::vector<DfTaskCtrl::Task2>& taskList,
                               PQ_PairArray *pI2PQ,
                               TlSparseMatrix *pSchwartzTable,
                               TlSparseMatrix *pDiagonalMat);

    // void finalizeI2PQ_A(PQ_PairArray_A *pI2PQ);
    // virtual void saveI2PQ_A(const PQ_PairArray_A& I2PQ);

    virtual std::vector<double>
    getSuperMatrixElementsA(const TlOrbitalInfoObject& orbInfo_p,
                            const TlOrbitalInfoObject& orbInfo_q,
                            const index_type G_row,
                            const std::vector<index_type>& G_col_list,
                            const PQ_PairArray& I2PQ,
                            const TlSparseMatrix& schwartzTable);
    std::vector<DfCD::IndexPair4A> 
    getCalcListA(const TlOrbitalInfoObject& orbInfo_p,
                 const TlOrbitalInfoObject& orbInfo_q,
                 const index_type G_row,
                 const std::vector<index_type>& G_col_list,
                 const PQ_PairArray& I2PQ);
    void calcERIsA(const TlOrbitalInfoObject& orbInfo_p,
                   const TlOrbitalInfoObject& orbInfo_q,
                   const std::vector<IndexPair4A>& calcList,
                   const TlSparseMatrix& schwartzTable);
    std::vector<double>
    setERIsA(const TlOrbitalInfoObject& orbInfo_p,
             const TlOrbitalInfoObject& orbInfo_q,
             const index_type G_row,
             const std::vector<index_type> G_col_list,
             const PQ_PairArray& I2PQ);
    // PQ_PairArray_A getI2PQ_A();
    TlMatrix getCholeskyVectorA(const TlOrbitalInfoObject& orbInfo_p,
                                const TlOrbitalInfoObject& orbInfo_q,
                                const TlVector& L_col,
                                const PQ_PairArray& I2PQ);

    /// デバッグ用にSuperMatrixを作成します
    ///
    /// V_pq,rs = (pq|rs) or (pqrs) 
    /// @param[in] orbInfo 軌道情報オブジェクト
    /// @retval supermatrix
    TlSymmetricMatrix getSuperMatrix(const TlOrbitalInfoObject& orbInfo, PQ_PairArray* pI2PQ);

    /// デバッグ用にSuperMatrix(交換項用)を作成します
    ///
    /// V_pr,qs = (pq|rs) 
    /// @param[in] orbInfo 軌道情報オブジェクト
    /// @retval supermatrix
    TlSymmetricMatrix getSuperMatrix_K_full(const TlOrbitalInfoObject& orbInfo, PQ_PairArray* pI2PR);
    TlSymmetricMatrix getSuperMatrix_K_half(const TlOrbitalInfoObject& orbInfo, PQ_PairArray* pI2PR);

    /// デバッグ用にSuperMatrixを作成します
    ///
    /// V_pq,rs = (pq|rs) or (pqrs) 
    /// @param[in] orbInfo_p pまたはrで示される軌道の軌道情報オブジェクト
    /// @param[in] orbInfo_q qまたはsで示される軌道の軌道情報オブジェクト
    /// @retval supermatrix
    TlSymmetricMatrix getSuperMatrix(const TlOrbitalInfoObject& orbInfo_p,
                                     const TlOrbitalInfoObject& orbInfo_q,
                                     PQ_PairArray* pI2PQ);

    /// コレスキー分解を行います(デバッグ用)
    ///
    /// V = L * L
    /// @param[in] V supermatrix
    /// @retval コレスキーベクトル(L)
    TlMatrix calcCholeskyVectors(const TlSymmetricMatrix& V);

    /// デバッグ用
    void getJ_A(TlSymmetricMatrix* pJ);

    /// デバッグ用
    void getK_A(const RUN_TYPE runType,
                TlSymmetricMatrix *pK);

protected:
    bool debugBuildSuperMatrix_;
    bool debugCheckCD_;

protected:
    virtual TlRowVectorMatrix calcCholeskyVectorsOnTheFlyS_new(const TlOrbitalInfoObject& orbInfo,
                                                               const std::string& I2PQ_path,
                                                               const double epsilon,
                                                               void(DfCD::*calcDiagonalsFunc)(
                                                                   const TlOrbitalInfoObject&,
                                                                   PQ_PairArray*,
                                                                   TlVector*),
                                                               std::vector<double>(DfCD::*getSuperMatrixElements)(
                                                                   const TlOrbitalInfoObject&,
                                                                   const index_type,
                                                                   const std::vector<index_type>&,
                                                                   const PQ_PairArray&));

    // K full ----------------------------------------------------------
    std::vector<double>
    getSuperMatrixElements_K_full(const TlOrbitalInfoObject& orbInfo,
                                  const index_type G_row,
                                  const std::vector<index_type>& G_col_list,
                                  const PQ_PairArray& I2PQ);
    virtual std::vector<DfCD::IndexPair4S> 
    getCalcList_K_full(const TlOrbitalInfoObject& orbInfo,
                       const index_type G_row,
                       const std::vector<index_type>& G_col_list,
                       const index_type start,
                       const index_type end,
                       const PQ_PairArray& I2PQ);
    std::vector<double>
    setERIs_K_full(const TlOrbitalInfoObject& orbInfo,
                   const index_type G_row,
                   const std::vector<index_type> G_col_list,
                   const index_type start,
                   const index_type end,
                   const PQ_PairArray& I2PQ);

    // K half ----------------------------------------------------------
    std::vector<double>
    getSuperMatrixElements_K_half(const TlOrbitalInfoObject& orbInfo,
                                  const index_type G_row,
                                  const std::vector<index_type>& G_col_list,
                                  const PQ_PairArray& I2PQ);
    
    virtual std::vector<DfCD::IndexPair4S> 
    getCalcList_K_half(const TlOrbitalInfoObject& orbInfo,
                       const index_type G_row,
                       const std::vector<index_type>& G_col_list,
                       const index_type start,
                       const index_type end,
                       const PQ_PairArray& I2PQ);
    void calcERIs_K(const TlOrbitalInfoObject& orbInfo,
                    const std::vector<IndexPair4S>& calcList);
    std::vector<double>
    setERIs_K_half(const TlOrbitalInfoObject& orbInfo,
                   const index_type G_row,
                   const std::vector<index_type> G_col_list,
                   const index_type start,
                   const index_type end,
                   const PQ_PairArray& I2PQ);

    bool getCachedValue(const TlOrbitalInfoObject& orbInfo,
                        const index_type indexP,
                        const index_type indexQ,
                        const index_type indexR,
                        const index_type indexS,
                        const ERI_CacheType& cache,
                        double* pValue);

    // for debug
    bool get_I_index(const PQ_PairArray& I2PQ, 
                     const index_type p, const index_type q,
                     index_type* pI);

    // debug
    PQ_PairArray debug_I2PQ_;
    TlSymmetricMatrix debug_V_;
};


template<class EngineClass>
TlRowVectorMatrix DfCD::calcCholeskyVectorsOnTheFly(const TlOrbitalInfoObject& orbInfo_p,
                                                     const TlOrbitalInfoObject& orbInfo_q,
                                                     const std::string& I2PQ_path)
{
    this->createEngines<EngineClass>();
    const TlRowVectorMatrix L = this->calcCholeskyVectorsOnTheFlyA(orbInfo_p,
                                                                    orbInfo_q,
                                                                    I2PQ_path);
    this->destroyEngines();

    return L;
}


template<class EngineClass>
void DfCD::createEngines()
{
    assert(this->pEngines_ == NULL);
    const int numOfThreads = this->numOfThreads_;
    this->log_.info(TlUtils::format("create engine: %d", numOfThreads));

    this->pEngines_ = new DfEngineObject*[numOfThreads];
    for (int i = 0; i < numOfThreads; ++i) {
        this->pEngines_[i] = new EngineClass;
    }
}

#endif // DFCD_H
