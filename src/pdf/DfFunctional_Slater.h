#ifndef DFFUNCTIONAL_SLATER_H
#define DFFUNCTIONAL_SLATER_H

#include "DfFunctional.h"

class DfFunctional_Slater : public DfFunctional_LDA {
public:
    DfFunctional_Slater();
    virtual ~DfFunctional_Slater();

public:
    virtual double getFunctional(double dRhoA, double dRhoB);
    virtual double getFunctional(double dRhoA);

    virtual void getDerivativeFunctional(double dRhoA, double dRhoB,
                                         double* pRoundF_roundRhoA, double* pRoundF_roundRhoB);
    virtual void getDerivativeFunctional(double dRhoA, double* pRoundF_roundRhoA);

private:
    static const double INV_3;
    static const double M_4_3;

    //static const double SLATER_ALPHA;
    static const double SLATER_COEF;
    static const double ROUND_SLATER_ROUND_RHO_COEF;
};

#endif // DFFUNCTIONAL_SLATER_H
