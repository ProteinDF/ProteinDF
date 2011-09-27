#include "DfFunctional_SVWN.h"

DfFunctional_SVWN::DfFunctional_SVWN()
{
}

DfFunctional_SVWN::~DfFunctional_SVWN()
{
}

double DfFunctional_SVWN::getFunctional(double dRhoA, double dRhoB)
{
    const double ex = this->m_Slater.getFunctional(dRhoA, dRhoB);
    const double ec = this->m_VWN.getFunctional(dRhoA, dRhoB);

    return ex + ec;
}

void DfFunctional_SVWN::getDerivativeFunctional(double dRhoA, double dRhoB,
                                                double* pRoundFunctional_roundRhoA, double* pRoundFunctional_roundRhoB)
{
    double dTmpXA, dTmpXB, dTmpCA, dTmpCB;
    this->m_Slater.getDerivativeFunctional(dRhoA, dRhoB, &dTmpXA, &dTmpXB);
    this->m_VWN.getDerivativeFunctional(dRhoA, dRhoB, &dTmpCA, &dTmpCB);

    *pRoundFunctional_roundRhoA = dTmpXA + dTmpCA;
    *pRoundFunctional_roundRhoB = dTmpXB + dTmpCB;
}


