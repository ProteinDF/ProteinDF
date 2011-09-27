#ifndef DFFUNCTIONAL_SHF_H
#define DFFUNCTIONAL_SHF_H

#include "DfFunctional.h"
#include "DfFunctional_Slater.h"

class DfFunctional_SHF : public DfFunctional_LDA {
public:
    DfFunctional_SHF();
    virtual ~DfFunctional_SHF();

public:
    virtual double getFunctional(double dRhoA, double dRhoB);
    virtual void getDerivativeFunctional(double dRhoA, double dRhoB,
                                         double* pRoundFunctional_roundRhoA, double* pRoundFunctional_roundRhoB);

protected:
    DfFunctional_Slater m_Slater;
};

#endif // DFFUNCTIONAL_SHF_H
