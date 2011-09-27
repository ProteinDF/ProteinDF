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

    virtual void calcForceFromK(RUN_TYPE runType);
    
   
protected:
    TlMatrix force_;
    TlOrbitalInfo orbitalInfo_;
    TlOrbitalInfo_Density orbitalInfoDens_;

    /// 全ての行列を出力するかどうか(debug用)
    bool isDebugOutMatrix_;

    double storedCutoffValue_;
};

#endif // DFFORCE_H
