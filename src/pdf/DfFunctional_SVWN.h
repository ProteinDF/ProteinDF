#ifndef DFFUNCTIONAL_SVWN_H
#define DFFUNCTIONAL_SVWN_H

#include "DfFunctional.h"
#include "DfFunctional_Slater.h"
#include "DfFunctional_VWN.h"

class DfFunctional_SVWN : public DfFunctional_LDA {
public:
    DfFunctional_SVWN();
    virtual ~DfFunctional_SVWN();

public:
    virtual double getFunctional(double dRhoA, double dRhoB);
    virtual void getDerivativeFunctional(double dRhoA, double dRhoB,
                                         double* pRoundFunctional_roundRhoA, double* pRoundFunctional_roundRhoB);

protected:
    DfFunctional_Slater m_Slater;
    DfFunctional_VWN m_VWN;
};

#endif // DFFUNCTIONAL_SVWN_H
