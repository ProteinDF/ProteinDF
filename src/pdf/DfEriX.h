#ifndef DFERIX_H
#define DFERIX_H

#include <vector>
#include "DfObject.h"
#include "DfEriEngine.h"
#include "DfTaskCtrl.h"
#include "TlOrbitalInfoObject.h"
#include "TlOrbitalInfo.h"
#include "TlOrbitalInfo_Density.h"
#include "TlSymmetricMatrix.h"
#include "TlSparseSymmetricMatrix.h"

//#define DEBUG_J
//#deinfe DEBUG_K

class DfEriX : public DfObject {
public:
    DfEriX(TlSerializeData* pPdfParam);
    virtual ~DfEriX();

public:
    void getDeltaT(const TlSymmetricMatrix& P, TlVector* pRho) {
        this->getJ(P, pRho);
    }
    void getJ(const TlSymmetricMatrix& P, TlVector* pRho);

    void getdeltaHpqA(const TlVector& rho, TlSymmetricMatrix& P) {
        this->getJ(rho, &P);
    }
    void getJ(const TlVector& rho, TlSymmetricMatrix* pP);

    /// J([pq | rs])
    void getJpq(const TlSymmetricMatrix& P, TlSymmetricMatrix* pJpq);
    
    /// J([alpha | beta])
    void getJab(TlSymmetricMatrix* pJab);

    void getForceJ(const TlSymmetricMatrix& P, TlMatrix* pForce);
    void getForceJ(const TlSymmetricMatrix& P, const TlVector& rho,
                   TlMatrix* pForce);

    void getForceJ(const TlVector& rho,
                   TlMatrix* pForce);

    void getK(const TlSymmetricMatrix& P, TlSymmetricMatrix* pK);

    void getForceK(const TlSymmetricMatrix& P, TlMatrix* pForce);

   
protected:
    static const int MAX_SHELL_TYPE;
    static const int FORCE_K_BUFFER_SIZE;
    
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
    
protected:
    virtual DfTaskCtrl* getDfTaskCtrlObject() const;

    virtual void finalize(TlMatrix* pMtx);
    virtual void finalize(TlSymmetricMatrix* pMtx);
    virtual void finalize(TlVector* pVct);
    
protected:
    /// カットオフ用統計変数を初期化する
    //void clearCutoffStats();

    /// カットオフレポートを表示する
    //virtual void cutoffReport();

    /// shellの種類別に軌道番号リストを作成する
    ShellArrayTable makeShellArrayTable(const TlOrbitalInfoObject& orbitalInfo);

    /// クーロン項の計算を行う
    ///
    /// 高速化無しに、式に書かれた通りに実装されている。
    void getJpq_exact(const TlSymmetricMatrix& P, TlSymmetricMatrix* pJ);

    /// クーロン項の計算を行う
    ///
    /// integral-driven法を用いる。
    //void getJpq_integralDriven(const TlSymmetricMatrix& P, TlSymmetricMatrix* pJ);
    void getJpq_integralDriven2(const TlSymmetricMatrix& P, TlSymmetricMatrix* pJ);
    
    void getJ_integralDriven_part(const TlOrbitalInfoObject& orbitalInfo,
                                  const std::vector<DfTaskCtrl::Task4>& taskList,
                                  const TlMatrixObject& P, TlMatrixObject* pJ);

    void storeJ_integralDriven(const index_type shellIndexP, const int maxStepsP,
                               const index_type shellIndexQ, const int maxStepsQ,
                               const index_type shellIndexR, const int maxStepsR,
                               const index_type shellIndexS, const int maxStepsS,
                               const DfEriEngine& engine,
                               const TlMatrixObject& P,
                               TlMatrixObject* pJ);

    void getJab_part(const TlOrbitalInfoObject& orbitalInfo_Density,
                     const std::vector<DfTaskCtrl::Task2>& taskList,
                     TlMatrixObject* pJab);
    
    void getK_exact(const TlSymmetricMatrix& P, TlSymmetricMatrix* pK);
    void getK_integralDriven(const TlSymmetricMatrix& P, TlSymmetricMatrix* pK);
    void storeK_integralDriven(const index_type shellIndexP, const int maxStepsP,
                               const index_type shellIndexQ, const int maxStepsQ,
                               const index_type shellIndexR, const int maxStepsR,
                               const index_type shellIndexS, const int maxStepsS,
                               const DfEriEngine& engine,
                               const TlMatrixObject& P,
                               TlMatrixObject* pK);
    void debugoutK_integralDriven() const;
    
    ShellPairArrayTable getShellPairArrayTable(const ShellArrayTable& shellArrayTable);

    // 入力された原子軌道群から、対となる原子軌道と有効な大きさを持つ原子軌道群を抽出し、返す
    // ShellArray selectShellArrayByDistribution(const ShellArray& inShellArray,
    //                                           const index_type companionShellIndex,
    //                                           const TlOrbitalInfoObject& orbitalInfo);

    // ShellPairArrayTable selectShellPairArrayTableByDensity(
    //     const ShellPairArrayTable& inShellPairArrayTable,
    //     const TlOrbitalInfoObject& orbitalInfo);

    virtual TlSparseSymmetricMatrix makeSchwarzTable(const TlOrbitalInfoObject& orbitalInfo);

    // bool isAliveBySchwarzCutoff(const index_type shellIndexP,
    //                             const index_type shellIndexQ,
    //                             const index_type shellIndexR,
    //                             const index_type shellIndexS,
    //                             const int shellQuartetType,
    //                             const TlSparseSymmetricMatrix& schwarzTable,
    //                             const double threshold);

protected:
    void getJ_part(const TlOrbitalInfo& orbitalInfo,
                   const TlOrbitalInfo_Density& orbitalInfo_Density,
                   const ShellArrayTable& shellArrayTable_Density,
                   const std::vector<DfTaskCtrl::Task2>& taskList,
                   const TlSparseSymmetricMatrix& schwarzTable,
                   const TlMatrixObject& P, TlVector* pRho);

    void getJ_part(const TlOrbitalInfo& orbitalInfo,
                   const TlOrbitalInfo_Density& orbitalInfo_Density,
                   const ShellArrayTable& shellArrayTable_Density,
                   const std::vector<DfTaskCtrl::Task2>& taskList,
                   const TlSparseSymmetricMatrix& schwarzTable,
                   const TlVector& rho, TlMatrixObject* pP);

    void getK_integralDriven_part(const TlOrbitalInfoObject& orbitalInfo,
                                  const std::vector<DfTaskCtrl::Task4>& taskList,
                                  const TlMatrixObject& P, TlMatrixObject* pK);
    
    void getForceJ_part(const TlOrbitalInfoObject& orbitalInfo,
                        const TlOrbitalInfoObject& orbitalInfo_Density,
                        const ShellArrayTable& shellArrayTable_Density,
                        std::vector<DfTaskCtrl::Task2>& taskList,
                        const TlSymmetricMatrix& P, const TlVector& rho,
                        TlMatrix* pForce);

    void storeForceJ(const index_type atomIndexA,
                     const index_type atomIndexB,
                     const index_type atomIndexC,
                     const index_type shellIndexP, const int maxStepsP,
                     const index_type shellIndexQ, const int maxStepsQ,
                     const index_type shellIndexR, const int maxStepsR,
                     const double* p_dJdA, const double* p_dJdB,
                     const TlMatrixObject& P,
                     const TlVectorObject& rho,
                     TlMatrixObject* pForce,
                     const int target, int* pIndex);
    
    void getForceJ_part(const TlOrbitalInfoObject& orbitalInfo_Density,
                        std::vector<DfTaskCtrl::Task2>& taskList,
                        const TlVector& rho,
                        TlMatrix* pForce);
    
    void storeForceJ(const index_type atomIndexA,
                     const index_type atomIndexC,
                     const index_type shellIndexP, const int maxStepsP,
                     const index_type shellIndexR, const int maxStepsR,
                     const DfEriEngine& engine,
                     const TlVectorObject& rho,
                     TlMatrixObject* pForce,
                     const int target, int* pIndex);
    
    void getForceJ_part(const TlOrbitalInfoObject& orbitalInfo,
                        const std::vector<DfTaskCtrl::Task4>& taskList,
                        const TlMatrixObject& P, TlMatrixObject* pForce);

    void storeForceJ_integralDriven(const int atomIndexA, const int atomIndexB,
                                    const int atomIndexC, const int atomIndexD,
                                    const index_type shellIndexP, const int maxStepsP,
                                    const index_type shellIndexQ, const int maxStepsQ,
                                    const index_type shellIndexR, const int maxStepsR,
                                    const index_type shellIndexS, const int maxStepsS,
                                    const DfEriEngine& engine,
                                    const TlMatrixObject& P,
                                    TlMatrixObject* pForce);

    void storeForceJ_integralDriven(const int atomIndexA, const int atomIndexB,
                                    const int atomIndexC, const int atomIndexD,
                                    const index_type shellIndexP, const int maxStepsP,
                                    const index_type shellIndexQ, const int maxStepsQ,
                                    const index_type shellIndexR, const int maxStepsR,
                                    const index_type shellIndexS, const int maxStepsS,
                                    const DfEriEngine& engine,
                                    const TlMatrixObject& P,
                                    TlMatrixObject* pForce,
                                    const int target, int* pIndex);

    void getForceK_part(const TlOrbitalInfoObject& orbitalInfo,
                        const std::vector<DfTaskCtrl::Task4>& taskList,
                        const TlMatrixObject& P, TlMatrixObject* pForce);

    void storeForceK_integralDriven(const int atomIndexA, const int atomIndexB,
                                    const int atomIndexC, const int atomIndexD,
                                    const index_type shellIndexP, const int maxStepsP,
                                    const index_type shellIndexQ, const int maxStepsQ,
                                    const index_type shellIndexR, const int maxStepsR,
                                    const index_type shellIndexS, const int maxStepsS,
                                    const DfEriEngine& engine,
                                    const TlMatrixObject& P,
                                    TlMatrixObject* pForce);

    void storeForceK_integralDriven(const int atomIndexA, const int atomIndexB,
                                    const int atomIndexC, const int atomIndexD,
                                    const index_type shellIndexP, const int maxStepsP,
                                    const index_type shellIndexQ, const int maxStepsQ,
                                    const index_type shellIndexR, const int maxStepsR,
                                    const index_type shellIndexS, const int maxStepsS,
                                    const DfEriEngine& engine,
                                    const TlMatrixObject& P,
                                    TlMatrixObject* pForce,
                                    const int target, int* pIndex);
    
protected:
    enum {
        X = 0,
        Y = 1,
        Z = 2
    };


protected:
    double cutoffThreshold_;
    
    // std::vector<unsigned long> cutoffAll_schwarz_;
    // std::vector<unsigned long> cutoffAlive_schwarz_;


    double lengthScaleParameter_;
    /// カットオフ用閾値
    /// J. Chem. Phys.,105,2726 (1996) : eq.32
    double cutoffEpsilon1_;

    /// カットオフ用閾値
    /// J. Chem. Phys.,105,2726 (1996) : eq.32
    double cutoffEpsilon2_;

    /// カットオフ用閾値
    /// J. Chem. Phys.,105,2726 (1996) : eq.33
    double cutoffEpsilon3_;
    
    // mutable std::vector<unsigned long> cutoffAll_E1_;
    // mutable std::vector<unsigned long> cutoffAlive_E1_;
    // mutable std::vector<unsigned long> cutoffAll_E2_;
    // mutable std::vector<unsigned long> cutoffAlive_E2_;

    DfEriEngine* pEriEngines_;
    
protected:
    /// デバッグ用積分インデックス積算クラス
    ///
    /// integral-driven法等、ループに応じて必要な積分を格納できたか確認するためのクラス
    class IntegralAggregater {
    public:
        explicit IntegralAggregater(index_type N = 1) : maxIndex_(N) {
            this->resize(N);
        }

    public:
        void resize(index_type N) {
            this->maxIndex_ = N;
            const std::size_t size = N * (N +1) / 2;
            //const std::size_t size = N * N;
            this->value_.clear();
            this->value_.resize(size);
            for (std::size_t i = 0; i < size; ++i) {
                this->value_[i].resize(size, 0);
            }
        }
        
        int countUp(index_type i, index_type j,
                    index_type k, index_type l,
                    int value) {
            assert((0 <= i) && (i < this->maxIndex_));
            assert((0 <= j) && (j < this->maxIndex_));
            assert((0 <= k) && (k < this->maxIndex_));
            assert((0 <= l) && (l < this->maxIndex_));

            if (i < j) {
                std::swap(i, j);
            }
            if (k < l) {
                std::swap(k, l);
            }
            const std::size_t ijKey = i + (2 * this->maxIndex_ - (j + 1)) * j / 2;
            const std::size_t klKey = k + (2 * this->maxIndex_ - (l + 1)) * l / 2;
            //const std::size_t ijKey = this->maxIndex_ * i + j;
            //const std::size_t klKey = this->maxIndex_ * k + l;
            
            //assert(ijKey < (this->maxIndex_ * (this->maxIndex_ + 1) / 2));
            //assert(klKey < (this->maxIndex_ * (this->maxIndex_ + 1) / 2));
            // std::cerr << TlUtils::format("(%2d %2d %2d %2d)=(%2d %2d)",
            //                              i, j, k, l,
            //                              ijKey, klKey)
            //           << std::endl;
            this->value_[ijKey][klKey] += value;

            // std::cerr << TlUtils::format("countup (%2d,%2d,%2d,%2d)=%2d",
            //                              i, j, k, l, this->value_[ijKey][klKey])
            //           << std::endl;
            return this->value_[ijKey][klKey];
        }

        int getCount(index_type i, index_type j,
                     index_type k, index_type l) const {
            assert((0 <= i) && (i < this->maxIndex_));
            assert((0 <= j) && (j < this->maxIndex_));
            assert((0 <= k) && (k < this->maxIndex_));
            assert((0 <= l) && (l < this->maxIndex_));

            if (i < j) {
                std::swap(i, j);
            }
            if (k < l) {
                std::swap(k, l);
            }

            const std::size_t ijKey = i + (2 * this->maxIndex_ - (j + 1)) * j / 2;
            const std::size_t klKey = k + (2 * this->maxIndex_ - (l + 1)) * l / 2;
            //const std::size_t ijKey = this->maxIndex_ * i + j;
            //const std::size_t klKey = this->maxIndex_ * k + l;

            //assert(ijKey < (this->maxIndex_ * (this->maxIndex_ + 1) / 2));
            //assert(klKey < (this->maxIndex_ * (this->maxIndex_ + 1) / 2));
            return this->value_[ijKey][klKey];
        }

    private:
        index_type maxIndex_;
        std::vector<std::vector<int> > value_;
    };

    /// クーロン項におけるデバッグ出力用変数
    bool isDebugOutJ_;
#ifdef DEBUG_J
    IntegralAggregater IA_J_ID1_;
    IntegralAggregater IA_J_ID2_;
#endif // DEBUG_J

    /// 交換項におけるデバッグ出力用変数
    bool isDebugOutK_;
#ifdef DEBUG_K
    IntegralAggregater IA_K_ID1_;
    IntegralAggregater IA_K_ID2_;
    IntegralAggregater IA_K_ID3_;
    IntegralAggregater IA_K_ID4_;
#endif // DEBUG_K
    
    /// 4中心積分によるクーロン項(J)を求める(デバッグ用)
    bool isDebugExactJ_;

    /// 4中心積分による交換項(K)を求める(デバッグ用)
    bool isDebugExactK_;

};

#endif // DFERIX_H
