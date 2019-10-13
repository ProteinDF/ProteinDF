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

#ifndef TLCOMBINEDENSITYMATRIX_H
#define TLCOMBINEDENSITYMATRIX_H

#include <vector>
#include "TlOrbitalInfoObject.h"

class TlCombineDensityMatrix {
   public:
    TlCombineDensityMatrix(double nearThreshold = 1.0E-3, bool verbose = false);
    virtual ~TlCombineDensityMatrix();

   public:
    void make(const TlOrbitalInfoObject& orbInfo_ref,
              const TlMatrixObject& P_ref,
              const TlOrbitalInfoObject& orbInfo_target,
              TlMatrixObject* pP_target);

    bool check(const TlMatrixObject& P);

    void setCheckPosition(bool yn) { this->isCheckAtomPosition_ = yn; }

   protected:
    std::vector<int> getMatchTable(const TlOrbitalInfoObject& refOrbInfo,
                                   const TlOrbitalInfoObject& targetOrbInfo);

   protected:
    /// 計算機epsilon
    static const double EPSILON;

    /// 同一原子と判断する原子間距離の閾値
    double nearThreshold_;

    /// 原子座標をチェックする(true)かどうか
    bool isCheckAtomPosition_;

    /// 冗長出力の有無
    bool isVerbose_;
};

#endif  // TLCOMBINEDENSITYMATRIX_H
