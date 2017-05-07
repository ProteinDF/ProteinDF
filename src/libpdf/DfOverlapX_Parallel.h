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
#include "TlDistributeSymmetricMatrix.h"

class DfOverlapX_Parallel : public DfOverlapX {
public:
    DfOverlapX_Parallel(TlSerializeData* pPdfParam);
    virtual ~DfOverlapX_Parallel();

public:
    void getSpqD(TlDistributeSymmetricMatrix* pSpq);
    void getSabD(TlDistributeSymmetricMatrix* pSab);

    virtual void getSgd(TlSymmetricMatrix* pSgd) {
        DfOverlapX::getSgd(pSgd);
    }
    void getSgd(TlDistributeSymmetricMatrix* pSab);

    /// 変換行列を作成する
    virtual void getTransMat(const TlOrbitalInfoObject& orbInfo1,
                             const TlOrbitalInfoObject& orbInfo2,
                             TlMatrix* pS) {
        DfOverlapX::getTransMat(orbInfo1, orbInfo2, pS);
    }
    void getTransMat(const TlOrbitalInfoObject& orbInfo1,
                     const TlOrbitalInfoObject& orbInfo2,
                     TlDistributeMatrix* pS);

    /// calc <pq gamma>
    virtual void get_pqg(const TlVector& myu, TlSymmetricMatrix* pF) {
        DfOverlapX::get_pqg(myu, pF);
    }
    virtual void get_pqg(const TlVector& myu, const TlVector& epsilon,
                         TlSymmetricMatrix* pF,
                         TlSymmetricMatrix* pE) {
        DfOverlapX::get_pqg(myu, epsilon, pF, pE);
    }
    virtual void get_pqg(const TlDistributeVector& myu,
                         TlDistributeSymmetricMatrix* pF);


    /// 重なり行列を作成する
    virtual void getOvpMat(const TlOrbitalInfoObject& orbInfo,
                           TlSymmetricMatrix* pS) {
        DfOverlapX::getOvpMat(orbInfo, pS);
    }
    void getOvpMat(const TlOrbitalInfoObject& orbitalInfo,
                   TlDistributeSymmetricMatrix* pS);

    /// <p|nabra|q> を求める
    virtual void getGradient(const TlOrbitalInfoObject& orbitalInfo,
                             TlMatrix* pMatX,
                             TlMatrix* pMatY,
                             TlMatrix* pMatZ) {
        DfOverlapX::getGradient(orbitalInfo, pMatX, pMatY, pMatZ);
    }
    void getGradient(const TlOrbitalInfoObject& orbitalInfo,
                     TlDistributeMatrix* pMatX,
                     TlDistributeMatrix* pMatY,
                     TlDistributeMatrix* pMatZ);

    virtual void getM(const TlSymmetricMatrix& P,
                      TlSymmetricMatrix* pM);
    virtual void getM_A(const TlSymmetricMatrix& P,
                        TlSymmetricMatrix* pM);

    void getM(const TlDistributeSymmetricMatrix& P,
              TlDistributeSymmetricMatrix* pM);
    void getM_A(const TlDistributeSymmetricMatrix& P,
                TlDistributeSymmetricMatrix* pM);
    
protected:
    virtual void logger(const std::string& str) const;

    virtual DfTaskCtrl* getDfTaskCtrlObject() const;

    virtual void finalize(TlMatrix* pMtx);
    virtual void finalize(TlSymmetricMatrix* pMtx);
    virtual void finalize(TlVector* pVct);
};

#endif // DFOVERLAPX_PARALLEL_H
