#include <cmath>
#include "DfFunctional_Slater.h"

#include "TlUtils.h"

const double DfFunctional_Slater::INV_3 = 1.0 / 3.0;
const double DfFunctional_Slater::M_4_3 = 4.0 / 3.0;

//const double DfFunctional_Slater::SLATER_ALPHA = 2.0 / 3.0;
const double DfFunctional_Slater::SLATER_COEF = - (9.0 / 4.0) * std::pow((3.0 / (4.0 * M_PI)), (1.0 / 3.0))
                                                * (2.0 / 3.0); // (2.0 / 3.0) means Slater's alpha
const double DfFunctional_Slater::ROUND_SLATER_ROUND_RHO_COEF = - 3.0 * std::pow((3.0 / (4.0 * M_PI)), (1.0 / 3.0))
                                                                * (2.0 / 3.0); // (2.0 / 3.0) means Slater's alpha

DfFunctional_Slater::DfFunctional_Slater()
{
}

DfFunctional_Slater::~DfFunctional_Slater()
{
}

// for UKS
double DfFunctional_Slater::getFunctional(const double dRhoA, const double dRhoB)
{
    const double dRhoATo4_3 = std::pow(dRhoA, M_4_3);
    const double dRhoBTo4_3 = std::pow(dRhoB, M_4_3);

    const double dAnswer = SLATER_COEF * (dRhoATo4_3 + dRhoBTo4_3);

    return dAnswer;
}

// specialized for RKS
double DfFunctional_Slater::getFunctional(const double dRhoA)
{
    const double dRhoATo4_3 = std::pow(dRhoA, M_4_3);

    const double dAnswer = SLATER_COEF * (2.0 * dRhoATo4_3);

    return dAnswer;
}

// for UKS
void DfFunctional_Slater::getDerivativeFunctional(const double dRhoA, const double dRhoB,
                                                  double* pRoundF_roundRhoA, double* pRoundF_roundRhoB)
{
    assert(pRoundF_roundRhoA != NULL);
    assert(pRoundF_roundRhoB != NULL);

    const double dRhoA_for_1_3 = std::pow(dRhoA, INV_3);
    const double dRhoB_for_1_3 = std::pow(dRhoB, INV_3);

    *pRoundF_roundRhoA = ROUND_SLATER_ROUND_RHO_COEF * dRhoA_for_1_3;
    *pRoundF_roundRhoB = ROUND_SLATER_ROUND_RHO_COEF * dRhoB_for_1_3;
}

// specialized for RKS
void DfFunctional_Slater::getDerivativeFunctional(const double dRhoA, double* pRoundF_roundRhoA)
{
    assert(pRoundF_roundRhoA != NULL);

    const double dRhoA_for_1_3 = std::pow(dRhoA, INV_3);

    *pRoundF_roundRhoA = ROUND_SLATER_ROUND_RHO_COEF * dRhoA_for_1_3;
}
