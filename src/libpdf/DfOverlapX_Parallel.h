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

class TlDenseGeneralMatrix_blacs;
class TlDenseSymmetricMatrix_blacs;
class TlDistributedVector;

class DfOverlapX_Parallel : public DfOverlapX {
 public:
  DfOverlapX_Parallel(TlSerializeData* pPdfParam);
  virtual ~DfOverlapX_Parallel();

 public:
  void getSpqD(TlDenseSymmetricMatrix_blacs* pSpq);
  void getSabD(TlDenseSymmetricMatrix_blacs* pSab);

  virtual void getSgd(TlDenseSymmetricMatrix_BLAS_Old* pSgd) {
    DfOverlapX::getSgd(pSgd);
  }
  void getSgd(TlDenseSymmetricMatrix_blacs* pSab);

  /// 変換行列を作成する
  virtual void getTransMat(const TlOrbitalInfoObject& orbInfo1,
                           const TlOrbitalInfoObject& orbInfo2,
                           TlDenseGeneralMatrix_BLAS_old* pS) {
    DfOverlapX::getTransMat(orbInfo1, orbInfo2, pS);
  }
  void getTransMat(const TlOrbitalInfoObject& orbInfo1,
                   const TlOrbitalInfoObject& orbInfo2,
                   TlDenseGeneralMatrix_blacs* pS);

  /// calc <pq gamma>
  virtual void get_pqg(const TlVector_BLAS& myu,
                       TlDenseSymmetricMatrix_BLAS_Old* pF) {
    DfOverlapX::get_pqg(myu, pF);
  }
  virtual void get_pqg(const TlVector_BLAS& myu, const TlVector_BLAS& epsilon,
                       TlDenseSymmetricMatrix_BLAS_Old* pF,
                       TlDenseSymmetricMatrix_BLAS_Old* pE) {
    DfOverlapX::get_pqg(myu, epsilon, pF, pE);
  }
  virtual void get_pqg(const TlDistributedVector& myu,
                       TlDenseSymmetricMatrix_blacs* pF);

  /// 重なり行列を作成する
  virtual void getOvpMat(const TlOrbitalInfoObject& orbInfo,
                         TlDenseSymmetricMatrix_BLAS_Old* pS) {
    DfOverlapX::getOvpMat(orbInfo, pS);
  }
  void getOvpMat(const TlOrbitalInfoObject& orbitalInfo,
                 TlDenseSymmetricMatrix_blacs* pS);

  /// <p|nabra|q> を求める
  virtual void getGradient(const TlOrbitalInfoObject& orbitalInfo,
                           TlDenseGeneralMatrix_BLAS_old* pMatX,
                           TlDenseGeneralMatrix_BLAS_old* pMatY,
                           TlDenseGeneralMatrix_BLAS_old* pMatZ) {
    DfOverlapX::getGradient(orbitalInfo, pMatX, pMatY, pMatZ);
  }
  void getGradient(const TlOrbitalInfoObject& orbitalInfo,
                   TlDenseGeneralMatrix_blacs* pMatX,
                   TlDenseGeneralMatrix_blacs* pMatY,
                   TlDenseGeneralMatrix_blacs* pMatZ);

  virtual void getM(const TlDenseSymmetricMatrix_BLAS_Old& P,
                    TlDenseSymmetricMatrix_BLAS_Old* pM);
  virtual void getM_A(const TlDenseSymmetricMatrix_BLAS_Old& P,
                      TlDenseSymmetricMatrix_BLAS_Old* pM);

  void getM(const TlDenseSymmetricMatrix_blacs& P,
            TlDenseSymmetricMatrix_blacs* pM);
  void getM_A(const TlDenseSymmetricMatrix_blacs& P,
              TlDenseSymmetricMatrix_blacs* pM);

 protected:
  virtual void logger(const std::string& str) const;

  virtual DfTaskCtrl* getDfTaskCtrlObject() const;

  virtual void finalize(TlDenseGeneralMatrix_BLAS_old* pMtx);
  virtual void finalize(TlDenseSymmetricMatrix_BLAS_Old* pMtx);
  virtual void finalize(TlVector_BLAS* pVct);
};

#endif  // DFOVERLAPX_PARALLEL_H
