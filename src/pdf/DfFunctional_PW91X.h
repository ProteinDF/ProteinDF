#ifndef DFFUNCTIONAL_PW91X_H
#define DFFUNCTIONAL_PW91X_H

#include "DfFunctional.h"

class DfFunctional_PW91X : public DfFunctional_GGA {
public:
    DfFunctional_PW91X();
    virtual ~DfFunctional_PW91X();

public:
    double getFunctional(double dRhoA, double dGammaAA);
    virtual double getFunctional(double dRhoA, double dRhoB, double dGammaAA, double dGammaAB, double dGammaBB);

    void getDerivativeFunctional(double dRhoA, double dGammaAA,
                                 double* pRoundF_roundRhoA, double* pRoundF_roundGammaAA, double* pRoundF_roundGammaAB);
    virtual void getDerivativeFunctional(double dRhoA, double dRhoB,
                                         double dGammaAA, double dGammaAB, double dGammaBB,
                                         double* pRoundF_roundRhoA, double* pRoundF_roundRhoB,
                                         double* pRoundF_roundGammaAA, double* pRoundF_roundGammaAB, double* pRoundF_roundGammaBB);

protected:
    double pw91x(double r, double s);
    double pw91x_roundRho(double r, double s);
    double pw91x_roundGamma(double r, double s);
};

#endif // DFFUNCTIONAL_PW91X_H
