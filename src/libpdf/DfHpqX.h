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

#ifndef DFHPQX_H
#define DFHPQX_H

#ifdef HAVE_CONFIG_H
#include "config.h"  // this file created by autotools
#endif               // HAVE_CONFIG_H

#include <vector>
#include "DfHpqEngine.h"
#include "DfObject.h"
#include "DfTaskCtrl.h"
#include "TlOrbitalInfo.h"
#include "tl_dense_general_matrix_lapack.h"
#include "tl_dense_symmetric_matrix_lapack.h"

class DfHpqX : public DfObject {
 public:
  typedef std::vector<index_type> ShellArray;
  typedef std::vector<ShellArray> ShellArrayTable;

 public:
  DfHpqX(TlSerializeData* pPdfParam);
  virtual ~DfHpqX();

 public:
  void getHpq(TlDenseSymmetricMatrix_Lapack* pHpq,
              TlDenseSymmetricMatrix_Lapack* pHpq2);

  /// Hpq由来の力成分を求める
  ///
  /// @param [in] P 密度行列
  /// @param [out] pForce 原子数(dummy chargeを含む)×3(x, y, z成分)の行列
  /// @param [out] pForce_Xonly 原子数(dummy chargeを含む)×3(x, y,
  /// z成分)の行列(dummy charge由来)
  void getForce(const TlDenseSymmetricMatrix_Lapack& P,
                TlDenseGeneralMatrix_Lapack* pForce,
                TlDenseGeneralMatrix_Lapack* pForce_Xonly);

  /// ESP を求める
  ///
  /// @param [in] P 密度行列
  /// @param [in] grids 評価グリッド座標群
  /// @return 評価グリッド群順にESP値を代入した配列
  std::vector<double> getESP(const TlMatrixObject& P,
                             const std::vector<TlPosition>& grids);

 protected:
  /// DfHpqEngineオブジェクトを作成する
  ///
  /// OpenMPスレッド数のオブジェクトを作成する。
  /// 富士通コンパイラではコンストラクタ中で
  /// オブジェクトを作成できないため。
  void createEngines();

  /// DfHpqEngineオブジェクトを破棄する
  void destroyEngines();

  virtual DfTaskCtrl* getDfTaskCtrlObject() const;

  virtual void finalize(TlDenseSymmetricMatrix_Lapack* pHpq,
                        TlDenseSymmetricMatrix_Lapack* pHpq2);
  virtual void finalize(std::vector<double>* pValues);

 protected:
  void getHpq_part(const TlOrbitalInfoObject& orbitalInfo,
                   const std::vector<DfTaskCtrl::Task2>& taskList,
                   const std::vector<TlAtom>& Cs, const std::vector<TlAtom>& Xs,
                   TlMatrixObject* pHpq, TlMatrixObject* pHpq2);

  void getForce_partProc(const TlOrbitalInfoObject& orbitalInfo,
                         const int shellTypeP, const int shellTypeQ,
                         const index_type shellIndexP,
                         const ShellArray& shellArrayQ, const TlMatrixObject& P,
                         TlDenseGeneralMatrix_Lapack* pForce,
                         TlDenseGeneralMatrix_Lapack* pForceX);

  void getESP_part(const TlOrbitalInfoObject& orbitalInfo,
                   const std::vector<DfTaskCtrl::Task2>& taskList,
                   const TlMatrixObject& P,
                   const std::vector<TlPosition>& grids,
                   std::vector<double>* pValues);

 protected:
  void makeShellArrayTable();
  DfHpqEngine::PGTOs getPGTOs(const index_type shellIndex);

 protected:
  enum { X = 0, Y = 1, Z = 2 };

 protected:
  DfHpqEngine* pEngines_;
  TlOrbitalInfo orbitalInfo_;
  ShellArrayTable shellArrayTable_;
};

#endif  // DFHPQX_H
