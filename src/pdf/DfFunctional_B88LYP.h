#ifndef DFFUNCTIONAL_B88LYP_H
#define DFFUNCTIONAL_B88LYP_H

#include "DfFunctional.h"
#include "DfFunctional_Becke88.h"
#include "DfFunctional_LYP.h"

class DfFunctional_B88LYP : public DfFunctional_GGA {
public:
    DfFunctional_B88LYP();
    virtual ~DfFunctional_B88LYP();

public:
    virtual double getFunctional(double dRhoA, double dRhoB, double dGammaAA, double dGammaAB, double dGammaBB);
    virtual void getDerivativeFunctional(double dRhoA, double dRhoB, double dGammaAA, double dGammaAB, double dGammaBB,
                                         double* pRoundF_roundRhoA, double* pRoundF_roundRhoB,
                                         double* pRoundF_roundGammaAA, double* pRoundF_roundGammaAB, double* pRoundF_roundGammaBB);
// public:
//   virtual double getEnergy(double dRhoA, double dRhoB, double dGammaAA, double dGammaAB, double dGammaBB);

//   virtual void buildFock(double dRhoA, double dRhoB,
//           double dGradRhoAx, double dGradRhoAy, double dGradRhoAz,
//           double dGradRhoBx, double dGradRhoBy, double dGradRhoBz,
//           const std::vector<DfCalcGridX::WFGrid>& aPhi,
//           const std::vector<DfCalcGridX::WFGrid>& aGradPhiX,
//           const std::vector<DfCalcGridX::WFGrid>& aGradPhiY,
//           const std::vector<DfCalcGridX::WFGrid>& aGradPhiZ,
//           double dWeight, TlMatrix_Symmetric& F);

// for Grid-Free =======================================================
protected:
    virtual TlVector getFunctionalTermCoef_GF();
    virtual TlVector getDerivativeFunctionalTermCoef_GF();

    virtual TlMatrix getFunctionalCore(const double rhoA, 
                                       const double rhoB,
                                       const double xA,
                                       const double xB);

    virtual TlMatrix getDerivativeFunctionalCore(const double rhoA,
                                                 const double rhoB,
                                                 const double xA,
                                                 const double xB);

protected:
    DfFunctional_Becke88 m_Becke88;
    DfFunctional_LYP m_LYP;
};

#endif // DFFUNCTIONAL_B88LYP_H
