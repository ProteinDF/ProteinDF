#include "DfFunctional_HFS.h"

DfFunctional_HFS::DfFunctional_HFS()
{
}

DfFunctional_HFS::~DfFunctional_HFS()
{
}

double DfFunctional_HFS::getFunctional(double dRhoA, double dRhoB)
{
    const double ex = this->m_Slater.getFunctional(dRhoA, dRhoB);

    return ex;
}

void DfFunctional_HFS::getDerivativeFunctional(double dRhoA, double dRhoB,
                                               double* pRoundFunctional_roundRhoA,
                                               double* pRoundFunctional_roundRhoB)
{
    double dTmpXA, dTmpXB;
    this->m_Slater.getDerivativeFunctional(dRhoA, dRhoB, &dTmpXA, &dTmpXB);

    *pRoundFunctional_roundRhoA = dTmpXA;
    *pRoundFunctional_roundRhoB = dTmpXB;
}


