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

#ifndef DFERIX_H
#define DFERIX_H

#include <vector>
#include "DfEriEngine.h"
#include "DfObject.h"
#include "DfTaskCtrl.h"
#include "TlOrbitalInfo.h"
#include "TlOrbitalInfoObject.h"
#include "TlOrbitalInfo_Density.h"
#include "tl_dense_general_matrix_lapack.h"
#include "tl_dense_symmetric_matrix_lapack.h"
#include "tl_dense_vector_lapack.h"
#include "tl_sparse_symmetric_matrix.h"

//#define DEBUG_J
//#deinfe DEBUG_K

class DfEriX : public DfObject {
 public:
  DfEriX(TlSerializeData* pPdfParam);
  virtual ~DfEriX();

 public:
  virtual void getJ(const TlDenseSymmetricMatrix_Lapack& P,
                    TlDenseVector_Lapack* pRho);
  void getJ(const TlDenseVector_Lapack& rho, TlDenseSymmetricMatrix_Lapack* pP);

  /// J([pq | rs])
  void getJpq(const TlDenseSymmetricMatrix_Lapack& P,
              TlDenseSymmetricMatrix_Lapack* pJpq);

  /// J([alpha | beta])
  virtual void getJab(TlDenseSymmetricMatrix_Lapack* pJab);

  void getForceJ(const TlDenseSymmetricMatrix_Lapack& P,
                 TlDenseGeneralMatrix_Lapack* pForce);
  void getForceJ(const TlDenseSymmetricMatrix_Lapack& P,
                 const TlDenseVector_Lapack& rho,
                 TlDenseGeneralMatrix_Lapack* pForce);

  void getForceJ(const TlDenseVector_Lapack& rho,
                 TlDenseGeneralMatrix_Lapack* pForce);

  void getK(const TlDenseSymmetricMatrix_Lapack& P,
            TlDenseSymmetricMatrix_Lapack* pK);

  void getForceK(const TlDenseSymmetricMatrix_Lapack& P,
                 TlDenseGeneralMatrix_Lapack* pForce);

 protected:
  /// DfEriEngineオブジェクトを作成する
  ///
  /// OpenMPスレッド数のオブジェクトを作成する。
  /// 富士通コンパイラではコンストラクタ中で
  /// オブジェクトを作成できないため。
  void createEngines();

  /// DfEriEngineオブジェクトを破棄する
  void destroyEngines();

 protected:
  // static const int MAX_SHELL_TYPE;
  static const int FORCE_K_BUFFER_SIZE;

  typedef std::vector<index_type> ShellArray;
  typedef std::vector<ShellArray> ShellArrayTable;

  struct ShellPair {
   public:
    ShellPair(index_type index1 = 0, index_type index2 = 0)
        : shellIndex1(index1), shellIndex2(index2) {}

   public:
    index_type shellIndex1;
    index_type shellIndex2;
  };
  typedef std::vector<ShellPair> ShellPairArray;
  typedef std::vector<ShellPairArray> ShellPairArrayTable;

 protected:
  virtual DfTaskCtrl* getDfTaskCtrlObject() const;

  virtual void finalize(TlDenseGeneralMatrix_Lapack* pMtx);
  virtual void finalize(TlDenseSymmetricMatrix_Lapack* pMtx);
  virtual void finalize(TlDenseVector_Lapack* pVct);

 protected:
  /// カットオフ用統計変数を初期化する
  // void clearCutoffStats();

  /// カットオフレポートを表示する
  // virtual void cutoffReport();

  /// shellの種類別に軌道番号リストを作成する
  ShellArrayTable makeShellArrayTable(const TlOrbitalInfoObject& orbitalInfo);

  /// クーロン項の計算を行う
  ///
  /// 高速化無しに、式に書かれた通りに実装されている。
  void getJpq_exact(const TlDenseSymmetricMatrix_Lapack& P,
                    TlDenseSymmetricMatrix_Lapack* pJ);

  /// クーロン項の計算を行う
  ///
  /// integral-driven法を用いる。
  void getJpq_integralDriven(const TlDenseSymmetricMatrix_Lapack& P,
                             TlDenseSymmetricMatrix_Lapack* pJ);

  int getJ_integralDriven_part(const TlOrbitalInfoObject& orbitalInfo,
                               const std::vector<DfTaskCtrl::Task4>& taskList,
                               const TlMatrixObject& P, index_type* pIndexPairs,
                               double* pValues);

  int storeJ_integralDriven(const index_type shellIndexP, const int maxStepsP,
                            const index_type shellIndexQ, const int maxStepsQ,
                            const index_type shellIndexR, const int maxStepsR,
                            const index_type shellIndexS, const int maxStepsS,
                            const DfEriEngine& engine, const TlMatrixObject& P,
                            index_type* pIndexPairs, double* pValues);

  void getJab_part(const TlOrbitalInfoObject& orbitalInfo_Density,
                   const std::vector<DfTaskCtrl::Task2>& taskList,
                   TlMatrixObject* pJab);

  void getK_exact(const TlDenseSymmetricMatrix_Lapack& P,
                  TlDenseSymmetricMatrix_Lapack* pK);
  void getK_integralDriven(const TlDenseSymmetricMatrix_Lapack& P,
                           TlDenseSymmetricMatrix_Lapack* pK);
  int storeK_integralDriven(const index_type shellIndexP, const int maxStepsP,
                            const index_type shellIndexQ, const int maxStepsQ,
                            const index_type shellIndexR, const int maxStepsR,
                            const index_type shellIndexS, const int maxStepsS,
                            const DfEriEngine& engine, const TlMatrixObject& P,
                            index_type* pIndexPairs, double* pValues);

  void debugoutK_integralDriven() const;

  ShellPairArrayTable getShellPairArrayTable(
      const ShellArrayTable& shellArrayTable);

  // 入力された原子軌道群から、対となる原子軌道と有効な大きさを持つ原子軌道群を抽出し、返す
  // ShellArray selectShellArrayByDistribution(const ShellArray& inShellArray,
  //                                           const index_type
  //                                           companionShellIndex, const
  //                                           TlOrbitalInfoObject&
  //                                           orbitalInfo);

  // ShellPairArrayTable selectShellPairArrayTableByDensity(
  //     const ShellPairArrayTable& inShellPairArrayTable,
  //     const TlOrbitalInfoObject& orbitalInfo);

  virtual TlSparseSymmetricMatrix makeSchwarzTable(
      const TlOrbitalInfoObject& orbitalInfo);

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
                 const TlMatrixObject& P, TlDenseVector_Lapack* pRho);

  void getJ_part(const TlOrbitalInfo& orbitalInfo,
                 const TlOrbitalInfo_Density& orbitalInfo_Density,
                 const ShellArrayTable& shellArrayTable_Density,
                 const std::vector<DfTaskCtrl::Task2>& taskList,
                 const TlDenseVector_Lapack& rho, TlMatrixObject* pP);

  ///
  /// @param [out] pIndexPQ
  /// K行列の行、列インデックスを格納する。総数はtaskListのサイズの2倍
  /// @param [out] pValues  K行列の要素を格納する。総数はtaskListのサイズ
  int getK_integralDriven_part(const TlOrbitalInfoObject& orbitalInfo,
                               const std::vector<DfTaskCtrl::Task4>& taskList,
                               const TlMatrixObject& P, index_type* pIndexPairs,
                               double* pValues);

  void getForceJ_part(const TlOrbitalInfoObject& orbitalInfo,
                      const TlOrbitalInfoObject& orbitalInfo_Density,
                      const ShellArrayTable& shellArrayTable_Density,
                      std::vector<DfTaskCtrl::Task2>& taskList,
                      const TlDenseSymmetricMatrix_Lapack& P,
                      const TlDenseVector_Lapack& rho,
                      TlDenseGeneralMatrix_Lapack* pForce);

  void storeForceJ(const index_type atomIndexA, const index_type atomIndexB,
                   const index_type atomIndexC, const index_type shellIndexP,
                   const int maxStepsP, const index_type shellIndexQ,
                   const int maxStepsQ, const index_type shellIndexR,
                   const int maxStepsR, const double* p_dJdA,
                   const double* p_dJdB, const TlMatrixObject& P,
                   const TlDenseVectorObject& rho, TlMatrixObject* pForce,
                   const int target, int* pIndex);

  void getForceJ_part(const TlOrbitalInfoObject& orbitalInfo_Density,
                      std::vector<DfTaskCtrl::Task2>& taskList,
                      const TlDenseVector_Lapack& rho,
                      TlDenseGeneralMatrix_Lapack* pForce);

  void storeForceJ(const index_type atomIndexA, const index_type atomIndexC,
                   const index_type shellIndexP, const int maxStepsP,
                   const index_type shellIndexR, const int maxStepsR,
                   const DfEriEngine& engine, const TlDenseVectorObject& rho,
                   TlMatrixObject* pForce, const int target, int* pIndex);

  void getForceJ_part(const TlOrbitalInfoObject& orbitalInfo,
                      const std::vector<DfTaskCtrl::Task4>& taskList,
                      const TlMatrixObject& P, TlMatrixObject* pForce);

  void storeForceJ_integralDriven(
      const int atomIndexA, const int atomIndexB, const int atomIndexC,
      const int atomIndexD, const index_type shellIndexP, const int maxStepsP,
      const index_type shellIndexQ, const int maxStepsQ,
      const index_type shellIndexR, const int maxStepsR,
      const index_type shellIndexS, const int maxStepsS,
      const DfEriEngine& engine, const TlMatrixObject& P,
      TlMatrixObject* pForce);

  void storeForceJ_integralDriven(
      const int atomIndexA, const int atomIndexB, const int atomIndexC,
      const int atomIndexD, const index_type shellIndexP, const int maxStepsP,
      const index_type shellIndexQ, const int maxStepsQ,
      const index_type shellIndexR, const int maxStepsR,
      const index_type shellIndexS, const int maxStepsS,
      const DfEriEngine& engine, const TlMatrixObject& P,
      TlMatrixObject* pForce, const int target, int* pIndex);

  void getForceK_part(const TlOrbitalInfoObject& orbitalInfo,
                      const std::vector<DfTaskCtrl::Task4>& taskList,
                      const TlMatrixObject& P, TlMatrixObject* pForce);

  void storeForceK_integralDriven(
      const int atomIndexA, const int atomIndexB, const int atomIndexC,
      const int atomIndexD, const index_type shellIndexP, const int maxStepsP,
      const index_type shellIndexQ, const int maxStepsQ,
      const index_type shellIndexR, const int maxStepsR,
      const index_type shellIndexS, const int maxStepsS,
      const DfEriEngine& engine, const TlMatrixObject& P,
      TlMatrixObject* pForce);

  void storeForceK_integralDriven(
      const int atomIndexA, const int atomIndexB, const int atomIndexC,
      const int atomIndexD, const index_type shellIndexP, const int maxStepsP,
      const index_type shellIndexQ, const int maxStepsQ,
      const index_type shellIndexR, const int maxStepsR,
      const index_type shellIndexS, const int maxStepsS,
      const DfEriEngine& engine, const TlMatrixObject& P,
      TlMatrixObject* pForce, const int target, int* pIndex);

 protected:
  enum { X = 0, Y = 1, Z = 2 };

 protected:
  double cutoffThreshold_;
  double cutoffEpsilon_density_;
  double cutoffEpsilon_distribution_;

  // std::vector<unsigned long> cutoffAll_schwarz_;
  // std::vector<unsigned long> cutoffAlive_schwarz_;

  double lengthScaleParameter_;
  // /// カットオフ用閾値
  // /// J. Chem. Phys.,105,2726 (1996) : eq.32
  // double cutoffEpsilon1_;

  // /// カットオフ用閾値
  // /// J. Chem. Phys.,105,2726 (1996) : eq.32
  // double cutoffEpsilon2_;

  /// カットオフ用閾値
  /// J. Chem. Phys.,105,2726 (1996) : eq.33
  double cutoffThreshold_primitive_;

  // mutable std::vector<unsigned long> cutoffAll_E1_;
  // mutable std::vector<unsigned long> cutoffAlive_E1_;
  // mutable std::vector<unsigned long> cutoffAll_E2_;
  // mutable std::vector<unsigned long> cutoffAlive_E2_;

  /// cutoff threshold of the density matrix (P) element for gradient
  double cutoffThreshold_P_grad_;

  DfEriEngine* pEriEngines_;
  std::vector<index_type>* pThreadIndexPairs_;
  std::vector<double>* pThreadValues_;

  // statics
  // double elapsetime_calc_;
  // double elapsetime_makepair_;
  // double elapsetime_calc_eri_;
  // double elapsetime_store_;
  // double elapsetime_sumup_;

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
      const std::size_t size = N * (N + 1) / 2;
      // const std::size_t size = N * N;
      this->value_.clear();
      this->value_.resize(size);
      for (std::size_t i = 0; i < size; ++i) {
        this->value_[i].resize(size, 0);
      }
    }

    int countUp(index_type i, index_type j, index_type k, index_type l,
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
      // const std::size_t ijKey = this->maxIndex_ * i + j;
      // const std::size_t klKey = this->maxIndex_ * k + l;

      // assert(ijKey < (this->maxIndex_ * (this->maxIndex_ + 1) / 2));
      // assert(klKey < (this->maxIndex_ * (this->maxIndex_ + 1) / 2));
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

    int getCount(index_type i, index_type j, index_type k, index_type l) const {
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
      // const std::size_t ijKey = this->maxIndex_ * i + j;
      // const std::size_t klKey = this->maxIndex_ * k + l;

      // assert(ijKey < (this->maxIndex_ * (this->maxIndex_ + 1) / 2));
      // assert(klKey < (this->maxIndex_ * (this->maxIndex_ + 1) / 2));
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
#endif  // DEBUG_J

  /// 交換項におけるデバッグ出力用変数
  bool isDebugOutK_;
#ifdef DEBUG_K
  IntegralAggregater IA_K_ID1_;
  IntegralAggregater IA_K_ID2_;
  IntegralAggregater IA_K_ID3_;
  IntegralAggregater IA_K_ID4_;
#endif  // DEBUG_K

  /// 4中心積分によるクーロン項(J)を求める(デバッグ用)
  bool isDebugExactJ_;

  /// 4中心積分による交換項(K)を求める(デバッグ用)
  bool isDebugExactK_;
};

#endif  // DFERIX_H
