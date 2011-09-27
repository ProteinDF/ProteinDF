#include "DfFunctional_SHF.h"

DfFunctional_SHF::DfFunctional_SHF()
{
}

DfFunctional_SHF::~DfFunctional_SHF()
{
}

double DfFunctional_SHF::getFunctional(double dRhoA, double dRhoB)
{
    const double ex = this->m_Slater.getFunctional(dRhoA, dRhoB);

    return ex;
}

void DfFunctional_SHF::getDerivativeFunctional(double dRhoA, double dRhoB,
                                               double* pRoundFunctional_roundRhoA, double* pRoundFunctional_roundRhoB)
{
    double dTmpXA, dTmpXB;
    this->m_Slater.getDerivativeFunctional(dRhoA, dRhoB, &dTmpXA, &dTmpXB);

    *pRoundFunctional_roundRhoA = dTmpXA;
    *pRoundFunctional_roundRhoB = dTmpXB;
}


