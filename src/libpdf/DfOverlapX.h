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

#include <vector>
#include "DfObject.h"
#include "DfOverlapEngine.h"
#include "DfTaskCtrl.h"

#include "TlMatrix.h"
#include "TlSymmetricMatrix.h"

class TlOrbitalInfoObject;

class DfOverlapX : public DfObject {
public:
    typedef std::vector<index_type> ShellArray;
    typedef std::vector<ShellArray> ShellArrayTable;
    
public:
    DfOverlapX(TlSerializeData* pPdfParam);
    virtual ~DfOverlapX();

public:
    void getSpq(TlSymmetricMatrix* pSpq);
    void getSab(TlSymmetricMatrix* pSab);
    virtual void getSgd(TlSymmetricMatrix* pSgd);
    void getNalpha(TlVector* pNalpha);
    void getForce(const TlSymmetricMatrix& W, TlMatrix* pForce);

    /// 重なり行列を作成する
    virtual void getOvpMat(const TlOrbitalInfoObject& orbInfo,
                           TlSymmetricMatrix* pS);

    /// 変換行列を作成する
    virtual void getTransMat(const TlOrbitalInfoObject& orbInfo1,
                             const TlOrbitalInfoObject& orbInfo2,
                             TlMatrix* pS);

    /// <p|nabra|q> を求める
    void getGradient(const TlOrbitalInfoObject& orbitalInfo,
                     TlMatrix* pMatX, TlMatrix* pMatY, TlMatrix* pMatZ);

    /// M_pq = sum_rs{P_rs <pqrs>}
    void getM(const TlSymmetricMatrix& P, TlSymmetricMatrix* pM);
    void getM_A(const TlSymmetricMatrix& P, TlSymmetricMatrix* pM);

public:
    /// calc <pq gamma>
    virtual void get_pqg(const TlVector& myu, TlSymmetricMatrix* pF);
    virtual void get_pqg(const TlVector& myu, const TlVector& epsilon,
                         TlSymmetricMatrix* pF,
                         TlSymmetricMatrix* pE);

protected:
    enum {
        X = 0,
        Y = 1,
        Z = 2
    };

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

    virtual void finalize(TlMatrix* pMtx);
    virtual void finalize(TlSymmetricMatrix* pMtx);
    virtual void finalize(TlVector* pVct);
    
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
                     const TlVector& myu,
                     const TlOrbitalInfoObject& orbitalInfo,
                     TlMatrixObject* pF);
    void calcOverlap_part(const TlOrbitalInfoObject& orbitalInfo_XC,
                          const TlVector& myu,
                          const TlOrbitalInfoObject& orbitalInfo,
                          const std::vector<DfTaskCtrl::Task2>& taskList,
                          TlMatrixObject* pMatrix);

    void calcOverlap(const TlOrbitalInfoObject& orbitalInfo_XC,
                     const TlVector& myu,
                     const TlVector& epsilon,
                     const TlOrbitalInfoObject& orbitalInfo,
                     TlMatrixObject* pF,
                     TlMatrixObject* pE);
    void calcOverlap_part(const TlOrbitalInfoObject& orbitalInfo_XC,
                          const TlVector& myu,
                          const TlVector& epsilon,
                          const TlOrbitalInfoObject& orbitalInfo,
                          const std::vector<DfTaskCtrl::Task2>& taskList,
                          TlMatrixObject* pF,
                          TlMatrixObject* pE);

    // <alpha>
    void calcOverlap(const TlOrbitalInfoObject& orbitalInfo,
                     TlVectorObject* pVector);
    void calcOverlap_part(const TlOrbitalInfoObject& orbitalInfo,
                          const std::vector<DfTaskCtrl::Task>& taskList,
                          TlVectorObject* pVector);
    
    ShellArrayTable makeShellArrayTable(const TlOrbitalInfoObject& orbitalInfo);

    void getForce_partProc(const TlOrbitalInfoObject& orbitalInfo,
                           const int shellTypeP, const int shellTypeQ,
                           const index_type shellIndexP,
                           const ShellArray& shellArrayQ,
                           const TlSymmetricMatrix& W,
                           TlMatrix* pForce);
    
    void getGradient_partProc(const TlOrbitalInfoObject& orbitalInfo,
                              const std::vector<DfTaskCtrl::Task2>& taskList,
                              TlMatrixObject* pMatX, TlMatrixObject* pMatY, TlMatrixObject* pMatZ);

    TlSparseSymmetricMatrix makeSchwarzTable(const TlOrbitalInfoObject& orbitalInfo);
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
                const DfOverlapEngine& engine,
                const TlMatrixObject& P,
                TlMatrixObject* pM);
    void storeM_A(const index_type shellIndexP, const int maxStepsP,
                  const index_type shellIndexQ, const int maxStepsQ,
                  const index_type shellIndexR, const int maxStepsR,
                  const index_type shellIndexS, const int maxStepsS,
                  const DfOverlapEngine& engine,
                  const TlMatrixObject& P,
                  TlMatrixObject* pM);

protected:
    DfOverlapEngine* pEngines_;
};

#endif // DFOVERLAPX_H
