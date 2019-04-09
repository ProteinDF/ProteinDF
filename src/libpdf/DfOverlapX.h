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

#ifndef DFOVERLAPX_H
#define DFOVERLAPX_H

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif  // HAVE_CONFIG_H

#include <vector>
#include "DfObject.h"
#include "DfOverlapEngine.h"
#include "DfTaskCtrl.h"

#include "tl_dense_general_matrix_lapack.h"
#include "tl_dense_symmetric_matrix_lapack.h"
#ifdef HAVE_EIGEN
#include "tl_dense_general_matrix_eigen.h"
#include "tl_dense_symmetric_matrix_eigen.h"
#endif  // HAVE_EIGEN

class TlOrbitalInfoObject;

class DfOverlapX : public DfObject {
   public:
    typedef std::vector<index_type> ShellArray;
    typedef std::vector<ShellArray> ShellArrayTable;

   public:
    DfOverlapX(TlSerializeData* pPdfParam);
    virtual ~DfOverlapX();

   public:
    void getSpq(TlDenseSymmetricMatrix_Lapack* pSpq);
    void getSab(TlDenseSymmetricMatrix_Lapack* pSab);
    virtual void getSgd(TlDenseSymmetricMatrix_Lapack* pSgd);
    void getNalpha(TlDenseVector_Lapack* pNalpha);
    void getForce(const TlDenseSymmetricMatrix_Lapack& W,
                  TlDenseGeneralMatrix_Lapack* pForce);

    /// 重なり行列を作成する
    virtual void getOvpMat(const TlOrbitalInfoObject& orbInfo,
                           TlDenseSymmetricMatrix_Lapack* pS);

    /// 変換行列を作成する
    virtual void getTransMat(const TlOrbitalInfoObject& orbInfo1,
                             const TlOrbitalInfoObject& orbInfo2,
                             TlDenseGeneralMatrix_Lapack* pS);

    /// <p|nabra|q> を求める
    void getGradient(const TlOrbitalInfoObject& orbitalInfo,
                     TlDenseGeneralMatrix_Lapack* pMatX,
                     TlDenseGeneralMatrix_Lapack* pMatY,
                     TlDenseGeneralMatrix_Lapack* pMatZ);

    /// M_pq = sum_rs{P_rs <pqrs>}
    void getM(const TlDenseSymmetricMatrix_Lapack& P,
              TlDenseSymmetricMatrix_Lapack* pM);
    void getM_A(const TlDenseSymmetricMatrix_Lapack& P,
                TlDenseSymmetricMatrix_Lapack* pM);

    void get_dM_exact(const TlDenseSymmetricMatrix_Eigen& P,
                      TlDenseGeneralMatrix_Eigen* pdMx,
                      TlDenseGeneralMatrix_Eigen* pdMy,
                      TlDenseGeneralMatrix_Eigen* pdMz, const int atomIndex);
    void get_dM_exact2(const TlOrbitalInfoObject& orbitalInfo_pq,
                       const TlOrbitalInfoObject& orbitalInfo_rs,
                       const TlDenseSymmetricMatrix_Eigen& P,
                       TlDenseGeneralMatrix_Eigen* pdMx,
                       TlDenseGeneralMatrix_Eigen* pdMy,
                       TlDenseGeneralMatrix_Eigen* pdMz, const int atomIndex);
    void get_dM(const TlOrbitalInfoObject& orbitalInfo_pq,
                const TlOrbitalInfoObject& orbitalInfo_rs,
                const TlDenseSymmetricMatrix_Eigen& P,
                TlDenseSymmetricMatrix_Eigen* pdMx,
                TlDenseSymmetricMatrix_Eigen* pdMy,
                TlDenseSymmetricMatrix_Eigen* pdMz, const int atomIndex);

   public:
    /// calc <pq gamma>
    virtual void get_pqg(const TlDenseVector_Lapack& myu,
                         TlDenseSymmetricMatrix_Lapack* pF);
    virtual void get_pqg(const TlDenseVector_Lapack& myu,
                         const TlDenseVector_Lapack& epsilon,
                         TlDenseSymmetricMatrix_Lapack* pF,
                         TlDenseSymmetricMatrix_Lapack* pE);

   protected:
    enum { X = 0, Y = 1, Z = 2 };

   protected:
    /// DfOverlapEngineオブジェクトを作成する
    ///
    /// OpenMPスレッド数のオブジェクトを作成する。
    /// 富士通コンパイラではコンストラクタ中で
    /// オブジェクトを作成できないため。
    void createEngines();

    /// DfOverlapEngineオブジェクトを破棄する
    void destroyEngines();

    virtual DfTaskCtrl* getDfTaskCtrlObject() const;

    virtual void finalize(TlDenseGeneralMatrix_Lapack* pMtx);
    virtual void finalize(TlDenseSymmetricMatrix_Lapack* pMtx);
    virtual void finalize(TlDenseVector_Lapack* pVct);

   protected:
    void calcOverlap(const TlOrbitalInfoObject& orbitalInfo,
                     TlMatrixObject* pMatrix);
    void calcOverlap(const TlOrbitalInfoObject& orbitalInfo1,
                     const TlOrbitalInfoObject& orbitalInfo2,
                     TlMatrixObject* pMatrix);

    void calcOverlap_part(const TlOrbitalInfoObject& orbitalInfo,
                          const std::vector<DfTaskCtrl::Task2>& taskList,
                          TlMatrixObject* pMatrix);
    void calcOverlap_part(const TlOrbitalInfoObject& orbitalInfo1,
                          const TlOrbitalInfoObject& orbitalInfo2,
                          const std::vector<DfTaskCtrl::Task2>& taskList,
                          TlMatrixObject* pMatrix);
    // <pq gamma>
    void calcOverlap(const TlOrbitalInfoObject& orbitalInfo_XC,
                     const TlDenseVector_Lapack& myu,
                     const TlOrbitalInfoObject& orbitalInfo,
                     TlMatrixObject* pF);
    void calcOverlap_part(const TlOrbitalInfoObject& orbitalInfo_XC,
                          const TlDenseVector_Lapack& myu,
                          const TlOrbitalInfoObject& orbitalInfo,
                          const std::vector<DfTaskCtrl::Task2>& taskList,
                          TlMatrixObject* pMatrix);

    void calcOverlap(const TlOrbitalInfoObject& orbitalInfo_XC,
                     const TlDenseVector_Lapack& myu,
                     const TlDenseVector_Lapack& epsilon,
                     const TlOrbitalInfoObject& orbitalInfo, TlMatrixObject* pF,
                     TlMatrixObject* pE);
    void calcOverlap_part(const TlOrbitalInfoObject& orbitalInfo_XC,
                          const TlDenseVector_Lapack& myu,
                          const TlDenseVector_Lapack& epsilon,
                          const TlOrbitalInfoObject& orbitalInfo,
                          const std::vector<DfTaskCtrl::Task2>& taskList,
                          TlMatrixObject* pF, TlMatrixObject* pE);

    // <alpha>
    void calcOverlap(const TlOrbitalInfoObject& orbitalInfo,
                     TlDenseVectorObject* pVector);
    void calcOverlap_part(const TlOrbitalInfoObject& orbitalInfo,
                          const std::vector<DfTaskCtrl::Task>& taskList,
                          TlDenseVectorObject* pVector);

    ShellArrayTable makeShellArrayTable(const TlOrbitalInfoObject& orbitalInfo);

    void getForce_partProc(const TlOrbitalInfoObject& orbitalInfo,
                           const int shellTypeP, const int shellTypeQ,
                           const index_type shellIndexP,
                           const ShellArray& shellArrayQ,
                           const TlDenseSymmetricMatrix_Lapack& W,
                           TlDenseGeneralMatrix_Lapack* pForce);

    void getGradient_partProc(const TlOrbitalInfoObject& orbitalInfo,
                              const std::vector<DfTaskCtrl::Task2>& taskList,
                              TlMatrixObject* pMatX, TlMatrixObject* pMatY,
                              TlMatrixObject* pMatZ);

    TlSparseSymmetricMatrix makeSchwarzTable(
        const TlOrbitalInfoObject& orbitalInfo);
    void getM_part(const TlOrbitalInfoObject& orbitalInfo,
                   const std::vector<DfTaskCtrl::Task4>& taskList,
                   const TlMatrixObject& P, TlMatrixObject* pM);
    void getM_part(const TlOrbitalInfoObject& orbitalInfo_PQ,
                   const TlOrbitalInfoObject& orbitalInfo_RS,
                   const std::vector<DfTaskCtrl::Task4>& taskList,
                   const TlMatrixObject& P, TlMatrixObject* pM);

    void storeM(const index_type shellIndexP, const int maxStepsP,
                const index_type shellIndexQ, const int maxStepsQ,
                const index_type shellIndexR, const int maxStepsR,
                const index_type shellIndexS, const int maxStepsS,
                const DfOverlapEngine& engine, const TlMatrixObject& P,
                TlMatrixObject* pM);
    void storeM_A(const index_type shellIndexP, const int maxStepsP,
                  const index_type shellIndexQ, const int maxStepsQ,
                  const index_type shellIndexR, const int maxStepsR,
                  const index_type shellIndexS, const int maxStepsS,
                  const DfOverlapEngine& engine, const TlMatrixObject& P,
                  TlMatrixObject* pM);

   protected:
    DfOverlapEngine* pEngines_;
};

#endif  // DFOVERLAPX_H
