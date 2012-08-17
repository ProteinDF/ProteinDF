#ifndef DFFUNCTIONAL_HFS_H
#define DFFUNCTIONAL_HFS_H

#include "DfFunctional.h"
#include "DfFunctional_Slater.h"

class DfFunctional_HFS : public DfFunctional_LDA {
public:
    DfFunctional_HFS();
    virtual ~DfFunctional_HFS();

public:
    virtual double getFunctional(double dRhoA, double dRhoB);
    virtual void getDerivativeFunctional(double dRhoA, double dRhoB,
                                         double* pRoundFunctional_roundRhoA,
                                         double* pRoundFunctional_roundRhoB);

protected:
    DfFunctional_Slater m_Slater;
};

#endif // DFFUNCTIONAL_HFS_H
