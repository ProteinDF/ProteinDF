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

#ifndef DFFORCE_H
#define DFFORCE_H

#include "DfObject.h"
#include "TlMatrix.h"
#include "TlVector.h"
#include "TlOrbitalInfo.h"
#include "TlOrbitalInfo_Density.h"

class DfForce : public DfObject {
protected:
    enum {
        X = 0,
        Y = 1,
        Z = 2
    };

public:
    DfForce(TlSerializeData* pPdfParam);
    virtual ~DfForce();

public:
    void calcForce();
    
protected:
    /// 力を保存する
    void saveForce();

    void outputStartTitle(const std::string& stepName, const char lineChar = '-');

    void outputEndTitle(const std::string& stepName, const char lineChar = '-');
    
    /// 核-核反発の力を計算する
    void calcForceFromNuclei();

    /// Weighted ...
    void calcForceFromWS(RUN_TYPE runType);
    
    TlMatrix getEnergyWeightedDensityMatrix(RUN_TYPE runType);

    void calcForceFromHpq(const TlSymmetricMatrix& P);
    
    ///
    void calcForceFromCoulomb(const RUN_TYPE runType);
    void calcForceFromCoulomb_exact(const RUN_TYPE runType);
    virtual void calcForceFromCoulomb_RIJ(const RUN_TYPE runType);

    void calcForceFromPureXC(RUN_TYPE runType);
    void calcForceFromPureXC_gridfree(RUN_TYPE runType);

    virtual void calcForceFromK(RUN_TYPE runType);
    
    TlMatrix getTransformMatrix(const TlMatrix& force);
   
protected:
    TlMatrix force_;
    TlOrbitalInfo orbitalInfo_;
    TlOrbitalInfo_Density orbitalInfoDens_;

    /// 全ての行列を出力するかどうか(debug用)
    bool isDebugOutMatrix_;

    double storedCutoffValue_;
};

#endif // DFFORCE_H
