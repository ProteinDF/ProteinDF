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

#ifndef DFERIX_PARALLEL_H
#define DFERIX_PARALLEL_H

#include "DfEriX.h"

class TlDenseSymmetricMatrix_Lapack;
class TlDenseGeneralMatrix_Scalapack;
class TlDenseSymmetricMatrix_Scalapack;
class TlDenseVector_Scalapack;
class TlDenseVector_Lapack;

class DfEriX_Parallel : public DfEriX {
   public:
    DfEriX_Parallel(TlSerializeData* pPdfParam);
    virtual ~DfEriX_Parallel();

   public:
    virtual void getJ(const TlDenseSymmetricMatrix_Lapack& P,
                      TlDenseVector_Lapack* pRho) {
        DfEriX::getJ(P, pRho);
    }
    void getJ(const TlDenseSymmetricMatrix_Scalapack& P,
              TlDenseVector_Scalapack* pRho);

    virtual void getJ(const TlDenseVector_Lapack& rho,
                      TlDenseSymmetricMatrix_Lapack* pJ) {
        DfEriX::getJ(rho, pJ);
    }
    void getJ(const TlDenseVector_Lapack& rho,
              TlDenseSymmetricMatrix_Scalapack* pJ);

    /// J([pq | rs])
    void getJpq_D(const TlDenseSymmetricMatrix_Scalapack& P,
                  TlDenseSymmetricMatrix_Scalapack* pJpq);

    /// J([alpha | beta])
    virtual void getJab(TlDenseSymmetricMatrix_Lapack* pJab) {
        DfEriX::getJab(pJab);
    }
    void getJab(TlDenseSymmetricMatrix_Scalapack* pJab);

    /// K
    void getK_D(const TlDenseSymmetricMatrix_Scalapack& P,
                TlDenseSymmetricMatrix_Scalapack* pK);

   protected:
    virtual DfTaskCtrl* getDfTaskCtrlObject() const;

    virtual void finalize(TlDenseGeneralMatrix_Lapack* pMtx);
    virtual void finalize(TlDenseSymmetricMatrix_Lapack* pMtx);
    virtual void finalize(TlDenseVector_Lapack* pVct);

    virtual TlSparseSymmetricMatrix makeSchwarzTable(
        const TlOrbitalInfoObject& orbitalInfo);

    void waitAnotherProcs(const TlDenseSymmetricMatrix_Scalapack& P);

   protected:
    // 非同期通信により電子密度行列を送受信する
    void getJ_D_BG(const TlDenseSymmetricMatrix_Scalapack& P,
                   TlDenseVector_Scalapack* pRho);

    void getJ_D_local(const TlDenseSymmetricMatrix_Scalapack& P,
                      TlDenseVector_Scalapack* pRho);

    void getJ_part2(const TlOrbitalInfo& orbitalInfo,
                    const TlOrbitalInfo_Density& orbitalInfo_Density,
                    const ShellArrayTable& shellArrayTable_Density,
                    const std::vector<DfTaskCtrl::Task2>& taskList,
                    // const TlSparseSymmetricMatrix& schwarzTable,
                    const TlDenseGeneralMatrix_Scalapack& P,
                    TlDenseVector_Lapack* pRho);

    void getK_D_BG(const TlDenseSymmetricMatrix_Scalapack& P,
                   TlDenseSymmetricMatrix_Scalapack* pK);

    void getK_D_local(const TlDenseSymmetricMatrix_Scalapack& P,
                      TlDenseSymmetricMatrix_Scalapack* pK);

    // void getK_integralDriven_part2(const TlOrbitalInfoObject& orbitalInfo,
    //                                const std::vector<DfTaskCtrl::Task4>&
    //                                taskList, const
    //                                TlDenseGeneralMatrix_Scalapack&
    //                                P,
    //                                TlMatrixObject* pK);

    // void storeK_integralDriven2(const index_type shellIndexP, const int
    // maxStepsP,
    //                             const index_type shellIndexQ, const int
    //                             maxStepsQ, const index_type shellIndexR,
    //                             const int maxStepsR, const index_type
    //                             shellIndexS, const int maxStepsS, const
    //                             DfEriEngine& engine, const
    //                             TlDenseGeneralMatrix_Scalapack& P,
    //                             TlMatrixObject* pK);

    void expandLocalDensityMatrix(const TlDenseSymmetricMatrix_Scalapack& P,
                                  const TlOrbitalInfo& orbInfo,
                                  TlDenseGeneralMatrix_Lapack* pLocalP,
                                  std::vector<index_type>* pRowIndexes,
                                  std::vector<index_type>* pColIndexes);
    std::vector<index_type> getExpandIndexes(
        const std::vector<index_type>& refIndexes,
        const TlOrbitalInfo& orbInfo);

   protected:
    enum { TAG_ALL_PROC_FINISHED = 9999 };

    enum CalcMode {
        CalcMode_Default = 0,
        CalcMode_BackGroundTransport = 1,
        CalcMode_UsingLocalMatrix = 2
    };

   private:
    int calcMode_;
};

#endif  // DFERIX_PARALLEL_H
