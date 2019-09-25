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

#ifndef DFOVERLAPX_PARALLEL_H
#define DFOVERLAPX_PARALLEL_H

#include "DfOverlapX.h"

class TlDenseGeneralMatrix_Scalapack;
class TlDenseSymmetricMatrix_Scalapack;
class TlDenseVector_Scalapack;

class DfOverlapX_Parallel : public DfOverlapX {
   public:
    DfOverlapX_Parallel(TlSerializeData* pPdfParam);
    virtual ~DfOverlapX_Parallel();

   public:
    void getSpqD(TlDenseSymmetricMatrix_Scalapack* pSpq);
    void getSabD(TlDenseSymmetricMatrix_Scalapack* pSab);

    virtual void getSgd(TlDenseSymmetricMatrix_Lapack* pSgd) {
        DfOverlapX::getSgd(pSgd);
    }
    void getSgd(TlDenseSymmetricMatrix_Scalapack* pSab);

    /// 変換行列を作成する
    virtual void getTransMat(const TlOrbitalInfoObject& orbInfo1,
                             const TlOrbitalInfoObject& orbInfo2,
                             TlDenseGeneralMatrix_Lapack* pS) {
        DfOverlapX::getTransMat(orbInfo1, orbInfo2, pS);
    }
    void getTransMat(const TlOrbitalInfoObject& orbInfo1,
                     const TlOrbitalInfoObject& orbInfo2,
                     TlDenseGeneralMatrix_Scalapack* pS);

    /// calc <pq gamma>
    virtual void get_pqg(const TlDenseVector_Lapack& myu,
                         TlDenseSymmetricMatrix_Lapack* pF) {
        DfOverlapX::get_pqg(myu, pF);
    }
    virtual void get_pqg(const TlDenseVector_Lapack& myu,
                         const TlDenseVector_Lapack& epsilon,
                         TlDenseSymmetricMatrix_Lapack* pF,
                         TlDenseSymmetricMatrix_Lapack* pE) {
        DfOverlapX::get_pqg(myu, epsilon, pF, pE);
    }
    virtual void get_pqg(const TlDenseVector_Scalapack& myu,
                         TlDenseSymmetricMatrix_Scalapack* pF);

    /// 重なり行列を作成する
    virtual void getOvpMat(const TlOrbitalInfoObject& orbInfo,
                           TlDenseSymmetricMatrix_Lapack* pS) {
        DfOverlapX::getOvpMat(orbInfo, pS);
    }
    void getOvpMat(const TlOrbitalInfoObject& orbitalInfo,
                   TlDenseSymmetricMatrix_Scalapack* pS);

    /// <p|nabra|q> を求める
    virtual void getGradient(const TlOrbitalInfoObject& orbitalInfo,
                             TlDenseGeneralMatrix_Lapack* pMatX,
                             TlDenseGeneralMatrix_Lapack* pMatY,
                             TlDenseGeneralMatrix_Lapack* pMatZ) {
        DfOverlapX::getGradient(orbitalInfo, pMatX, pMatY, pMatZ);
    }
    void getGradient(const TlOrbitalInfoObject& orbitalInfo,
                     TlDenseGeneralMatrix_Scalapack* pMatX,
                     TlDenseGeneralMatrix_Scalapack* pMatY,
                     TlDenseGeneralMatrix_Scalapack* pMatZ);

    virtual void getM(const TlDenseSymmetricMatrix_Lapack& P,
                      TlDenseSymmetricMatrix_Lapack* pM);
    virtual void getM_A(const TlDenseSymmetricMatrix_Lapack& P,
                        TlDenseSymmetricMatrix_Lapack* pM);

    void getM(const TlDenseSymmetricMatrix_Scalapack& P,
              TlDenseSymmetricMatrix_Scalapack* pM);
    void getM_A(const TlDenseSymmetricMatrix_Scalapack& P,
                TlDenseSymmetricMatrix_Scalapack* pM);

   protected:
    virtual void logger(const std::string& str) const;

    virtual DfTaskCtrl* getDfTaskCtrlObject() const;

    virtual void finalize(TlDenseGeneralMatrix_Lapack* pMtx);
    virtual void finalize(TlDenseSymmetricMatrix_Lapack* pMtx);
    virtual void finalize(TlDenseVector_Lapack* pVct);
};

#endif  // DFOVERLAPX_PARALLEL_H
