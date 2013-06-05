#include "DfFunctional_LYP.h"
#include "TlUtils.h"

const double DfFunctional_LYP::TOLERANCE = 1.0E-20;
const double DfFunctional_LYP::INV_3 = 1.0 / 3.0;
const double DfFunctional_LYP::INV_9 = 1.0 / 9.0;
const double DfFunctional_LYP::M_5_3 = 5.0 / 3.0;
const double DfFunctional_LYP::M_8_3 = 8.0 / 3.0;
const double DfFunctional_LYP::M_4_3 = 4.0 / 3.0;
const double DfFunctional_LYP::M_7_9 = 7.0 / 9.0;
const double DfFunctional_LYP::M_11_3 = 11.0 / 3.0;
const double DfFunctional_LYP::LYP_COEF = std::pow(2.0, 11.0 / 3.0) * 0.3 * std::pow((3.0 * M_PI * M_PI), 2.0 / 3.0);

const double DfFunctional_LYP::LYP_PARAM_A = 0.04918;
const double DfFunctional_LYP::LYP_PARAM_B = 0.132;
const double DfFunctional_LYP::LYP_PARAM_C = 0.2533;
const double DfFunctional_LYP::LYP_PARAM_D = 0.349;
const double DfFunctional_LYP::LYP_PARAM_AB = DfFunctional_LYP::LYP_PARAM_A * DfFunctional_LYP::LYP_PARAM_B;
const double DfFunctional_LYP::LYP_PARAM_DD = DfFunctional_LYP::LYP_PARAM_D * DfFunctional_LYP::LYP_PARAM_D;

DfFunctional_LYP::DfFunctional_LYP()
{
}

DfFunctional_LYP::~DfFunctional_LYP()
{
}

double DfFunctional_LYP::getFunctional(const double dRhoA, const double dRhoB,
                                       const double dGammaAA, const double dGammaAB, const double dGammaBB)
{
    return this->LYP(dRhoA, dRhoB, dGammaAA, dGammaAB, dGammaBB);
}

double DfFunctional_LYP::getFunctional(const double dRhoA, const double dGammaAA)
{
    return this->LYP(dRhoA, dGammaAA);
}

void DfFunctional_LYP::getDerivativeFunctional(const double dRhoA,
                                               const double dRhoB,
                                               const double dGammaAA,
                                               const double dGammaAB,
                                               const double dGammaBB,
                                               double* pRoundF_roundRhoA,
                                               double* pRoundF_roundRhoB,
                                               double* pRoundF_roundGammaAA,
                                               double* pRoundF_roundGammaAB,
                                               double* pRoundF_roundGammaBB)
{
    assert(pRoundF_roundRhoA != NULL);
    assert(pRoundF_roundRhoB != NULL);
    assert(pRoundF_roundGammaAA != NULL);
    assert(pRoundF_roundGammaAB != NULL);
    assert(pRoundF_roundGammaBB != NULL);

    // initialize
    *pRoundF_roundRhoA = 0.0;
    *pRoundF_roundRhoB = 0.0;
    *pRoundF_roundGammaAA = 0.0;
    *pRoundF_roundGammaAB = 0.0;
    *pRoundF_roundGammaBB = 0.0;

    const double dRho = dRhoA + dRhoB;
    if ((std::fabs(dRhoA) > TOLERANCE) && (std::fabs(dRhoB) > TOLERANCE)) {
        const double dRhoAB = dRhoA * dRhoB;
        const double dInvRho = 1.0 / dRho;
        const double dInvRhoA = 1.0 / dRhoA;
        const double dInvRhoB = 1.0 / dRhoB;

        const double dRhoTo1_3 = std::pow(dRho, INV_3); // rho^(1/3)
        const double dRhoToM1_3 = 1.0 / dRhoTo1_3; // rho^(-1/3)
        const double dRhoToM4_3 = std::pow(dRho, - 4.0 / 3.0); // rho^(-4/3)

        const double dRhoATo8_3 = std::pow(dRhoA, 8.0 / 3.0); // rhoA^(8/3)
        const double dRhoBTo8_3 = std::pow(dRhoB, 8.0 / 3.0); // rhoB^(8/3)

        const double dOmega = this->omega(dRho);
        const double ABOmega = LYP_PARAM_AB * dOmega;
        const double dOmegaPrime = this->omega_prime(dRho, dOmega);
        const double dDelta = this->delta(dRho);
        const double dDeltaPrime = this->delta_prime(dRho, dDelta);

        // roundF_roundRho
        double dRoundRoundLYP_roundRhoARoundGammaAA, dRoundRoundLYP_roundRhoBRoundGammaAA;
        double dRoundRoundLYP_roundRhoARoundGammaAB, dRoundRoundLYP_roundRhoBRoundGammaAB;
        double dRoundRoundLYP_roundRhoARoundGammaBB, dRoundRoundLYP_roundRhoBRoundGammaBB;
        this->roundRoundLYP_roundRhoRoundGamma(dRho, dRhoA, dRhoB, dInvRho, dRhoAB, dOmega, dOmegaPrime, dDelta, dDeltaPrime,
                                               &dRoundRoundLYP_roundRhoARoundGammaAA, &dRoundRoundLYP_roundRhoBRoundGammaAA,
                                               &dRoundRoundLYP_roundRhoARoundGammaAB, &dRoundRoundLYP_roundRhoBRoundGammaAB,
                                               &dRoundRoundLYP_roundRhoARoundGammaBB, &dRoundRoundLYP_roundRhoBRoundGammaBB);

        const double dArg1_coef = - 4.0 * LYP_PARAM_A / (1.0 + LYP_PARAM_D * dRhoToM1_3) * dRhoAB / dRho;
        const double dArg1A = dArg1_coef * (INV_3 * LYP_PARAM_D * dRhoToM4_3 / (1.0 + LYP_PARAM_D * dRhoToM1_3) + dInvRhoA - dInvRho);
        const double dArg1B = dArg1_coef * (INV_3 * LYP_PARAM_D * dRhoToM4_3 / (1.0 + LYP_PARAM_D * dRhoToM1_3) + dInvRhoB - dInvRho);

        const double dArg2A = - LYP_COEF * LYP_PARAM_AB * (dOmegaPrime * dRhoAB * (dRhoATo8_3 + dRhoBTo8_3)
                                                           + dOmega * dRhoB * (M_11_3 * dRhoATo8_3 + dRhoBTo8_3));
        const double dArg2B = - LYP_COEF * LYP_PARAM_AB * (dOmegaPrime * dRhoAB * (dRhoBTo8_3 + dRhoATo8_3)
                                                           + dOmega * dRhoA * (M_11_3 * dRhoBTo8_3 + dRhoATo8_3));

        // for debug
        //   std::cout << "\n"
        //      << TlUtils::format("dRoundRoundLYP_roundRhoARoundGammaAA = %+e\n", dRoundRoundLYP_roundRhoARoundGammaAA)
        //      << TlUtils::format("dRoundRoundLYP_roundRhoARoundGammaAB = %+e\n", dRoundRoundLYP_roundRhoARoundGammaAB)
        //      << TlUtils::format("dRoundRoundLYP_roundRhoARoundGammaBB = %+e\n", dRoundRoundLYP_roundRhoARoundGammaBB)
        //      << TlUtils::format("dRoundRoundLYP_roundRhoBRoundGammaAA = %+e\n", dRoundRoundLYP_roundRhoBRoundGammaAA)
        //      << TlUtils::format("dRoundRoundLYP_roundRhoBRoundGammaAB = %+e\n", dRoundRoundLYP_roundRhoBRoundGammaAB)
        //      << TlUtils::format("dRoundRoundLYP_roundRhoBRoundGammaBB = %+e\n", dRoundRoundLYP_roundRhoBRoundGammaBB)
        //      << std::endl;

        // alpha spin ======================================================
        if (dRhoA > TOLERANCE) {
            // roundF_roundRho
            *pRoundF_roundRhoA = dArg1A + dArg2A
                                 + dRoundRoundLYP_roundRhoARoundGammaAA * dGammaAA
                                 + dRoundRoundLYP_roundRhoARoundGammaAB * dGammaAB
                                 + dRoundRoundLYP_roundRhoARoundGammaBB * dGammaBB;

            // roundF_roundGamma
            *pRoundF_roundGammaAA = - ABOmega * (INV_9 * dRhoAB * (1.0 - 3.0 * dDelta - (dDelta - 11.0) * dRhoA / dRho) - dRhoB * dRhoB);
        }

        // beta spin =======================================================
        if (dRhoB > TOLERANCE) {
            // roundF_roundRho
            *pRoundF_roundRhoB = dArg1B + dArg2B
                                 + dRoundRoundLYP_roundRhoBRoundGammaAA * dGammaAA
                                 + dRoundRoundLYP_roundRhoBRoundGammaAB * dGammaAB
                                 + dRoundRoundLYP_roundRhoBRoundGammaBB * dGammaBB;

            // roundF_roundGamma
            *pRoundF_roundGammaBB = - ABOmega * (INV_9 * dRhoAB * (1.0 - 3.0 * dDelta - (dDelta - 11.0) * dRhoB / dRho) - dRhoA * dRhoA);
        }

        // both ============================================================
        if ((dRhoA > TOLERANCE) && (dRhoB > TOLERANCE)) {
            // roundF_roundGamma
            *pRoundF_roundGammaAB = - ABOmega * (INV_9 * dRhoAB * (47.0 - 7.0 * dDelta) - M_4_3 * dRho * dRho);
        }
    }
}

void DfFunctional_LYP::getDerivativeFunctional(const double dRhoA, const double dGammaAA,
                                               double* pRoundF_roundRhoA, double* pRoundF_roundGammaAA, double* pRoundF_roundGammaAB)
{
    double dRoundF_roundRhoB;
    double dRoundF_roundGammaBB;

    this->getDerivativeFunctional(dRhoA, dRhoA, dGammaAA, dGammaAA, dGammaAA,
                                  pRoundF_roundRhoA, &dRoundF_roundRhoB,
                                  pRoundF_roundGammaAA, pRoundF_roundGammaAB, &dRoundF_roundGammaBB);
}

// void DfFunctional_LYP::getDerivativeFunctional(const double dRhoA, const double dGammaAA,
//                         double* pRoundF_roundRhoA, double* pRoundF_roundGammaAA, double* pRoundF_roundGammaAB){
//   assert(pRoundF_roundRhoA != NULL);
//   assert(pRoundF_roundGammaAA != NULL);
//   assert(pRoundF_roundGammaAB != NULL);

//   // initialize
//   *pRoundF_roundRhoA = 0.0;
//   *pRoundF_roundGammaAA = 0.0;
//   *pRoundF_roundGammaAB = 0.0;

//   const double dRho = 2.0 * dRhoA;
//   if (dRho > TOLERANCE){
//     const double dRhoAB = dRhoA * dRhoA;
//     const double dInvRho = 1.0 / dRho;
//     const double dInvRhoA = 1.0 / dRhoA;
//     //     const double dInvRhoB = 1.0 / dRhoB;

//     const double dRhoTo1_3 = std::pow(dRho, INV_3); // rho^(1/3)
//     const double dRhoToM1_3 = 1.0 / dRhoTo1_3; // rho^(-1/3)
//     const double dRhoToM4_3 = std::pow(dRho, - 4.0 / 3.0); // rho^(-4/3)

//     const double dRhoATo8_3 = std::pow(dRhoA, 8.0 / 3.0); // rhoA^(8/3)
//     //     const double dRhoBTo8_3 = std::pow(dRhoB, 8.0 / 3.0); // rhoB^(8/3)

//     const double dOmega = this->omega(dRho);
//     const double ABOmega = LYP_PARAM_AB * dOmega;
//     const double dOmegaPrime = this->omega_prime(dRho, dOmega);
//     const double dDelta = this->delta(dRho);
//     const double dDeltaPrime = this->delta_prime(dRho, dDelta);

//     // roundF_roundRho
//     double dRoundRoundLYP_roundRhoARoundGammaAA;//, dRoundRoundLYP_roundRhoBRoundGammaAA;
//     double dRoundRoundLYP_roundRhoARoundGammaAB;//, dRoundRoundLYP_roundRhoBRoundGammaAB;
//     //double dRoundRoundLYP_roundRhoARoundGammaBB, dRoundRoundLYP_roundRhoBRoundGammaBB;
//     this->roundRoundLYP_roundRhoRoundGamma(dRho, dRhoA, dInvRho, dRhoAB, dOmega, dOmegaPrime, dDelta, dDeltaPrime,
//                     &dRoundRoundLYP_roundRhoARoundGammaAA, &dRoundRoundLYP_roundRhoARoundGammaAB);


//     const double dArg1_coef = - 4.0 * LYP_PARAM_A / (1.0 + LYP_PARAM_D * dRhoToM1_3) * dRhoAB / dRho;
//     const double dArg1A = dArg1_coef * (INV_3 * LYP_PARAM_D * dRhoToM4_3 / (1.0 + LYP_PARAM_D * dRhoToM1_3) + dInvRhoA - dInvRho);
// //     const double dArg1B = dArg1_coef * (INV_3 * LYP_PARAM_D * dRhoToM4_3 / (1.0 + LYP_PARAM_D * dRhoToM1_3) + dInvRhoB - dInvRho);

//     const double dArg2A = - LYP_COEF * LYP_PARAM_AB * (dOmegaPrime * dRhoAB * (2.0 * dRhoATo8_3)
//                             + dOmega * dRhoA * (M_11_3 * (2.0 * dRhoATo8_3)));
// //     const double dArg2B = - LYP_COEF * LYP_PARAM_AB * (dOmegaPrime * dRhoAB * (dRhoBTo8_3 + dRhoATo8_3)
// //                              + dOmega * dRhoA * (M_11_3 * dRhoBTo8_3 + dRhoATo8_3));

//     // roundF_roundRho
//     *pRoundF_roundRhoA = dArg1A + dArg2A
//       + 2.0 * dRoundRoundLYP_roundRhoARoundGammaAA * dGammaAA
//       + dRoundRoundLYP_roundRhoARoundGammaAB * dGammaAA;

//     // roundF_roundGamma
//     *pRoundF_roundGammaAA = - ABOmega * (INV_9 * dRhoAB * (1.0 - 3.0 * dDelta - (dDelta - 11.0) * dRhoA / dRho) - dRhoA * dRhoA);
//     *pRoundF_roundGammaAB = - ABOmega * (INV_9 * dRhoAB * (47.0 - 7.0 * dDelta) - M_4_3 * dRho * dRho);
//   }
// }

double DfFunctional_LYP::LYP(const double dRhoA, const double dRhoB,
                             const double dGammaAA, const double dGammaAB, const double dGammaBB)
{
    const double dRho = dRhoA + dRhoB;
    const double dRhoAB = dRhoA * dRhoB;

    const double dPow_rho_m1_3 = std::pow(dRho, - INV_3);
    const double dPow_rhoA_8_3 = std::pow(dRhoA, M_8_3);
    const double dPow_rhoB_8_3 = std::pow(dRhoB, M_8_3);

    const double dOmega = this->omega(dRho);
    const double dDelta = this->delta(dRho);

    double dRoundLYP_roundGammaAA, dRoundLYP_roundGammaAB, dRoundLYP_roundGammaBB;
    this->roundLYP_roundGamma(dRho, dRhoA, dRhoB, dOmega, dDelta,
                              dRoundLYP_roundGammaAA, dRoundLYP_roundGammaAB, dRoundLYP_roundGammaBB);

    const double term1 = - 4.0 * LYP_PARAM_A / (1.0 + LYP_PARAM_D * dPow_rho_m1_3) * dRhoAB / dRho;
    const double term2 = - LYP_COEF * LYP_PARAM_AB * dOmega * dRhoAB * (dPow_rhoA_8_3 + dPow_rhoB_8_3);

    const double dAnswer = term1 + term2
                           + dRoundLYP_roundGammaAA * dGammaAA + dRoundLYP_roundGammaAB * dGammaAB + dRoundLYP_roundGammaBB * dGammaBB;

    return dAnswer;
}

// for RKS
double DfFunctional_LYP::LYP(const double dRhoA, const double dGammaAA)
{
    const double dRho = 2.0 * dRhoA;
    const double dRhoAB = dRhoA * dRhoA;

    const double dPow_rho_m1_3 = std::pow(dRho, - INV_3);
    const double dPow_rhoA_8_3 = std::pow(dRhoA, M_8_3);

    const double dOmega = this->omega(dRho);
    const double dDelta = this->delta(dRho);

    double dRoundLYP_roundGammaAA, dRoundLYP_roundGammaAB;
    this->roundLYP_roundGamma(dRho, dRhoA, dOmega, dDelta,
                              dRoundLYP_roundGammaAA, dRoundLYP_roundGammaAB);

    const double term1 = - 4.0 * LYP_PARAM_A / (1.0 + LYP_PARAM_D * dPow_rho_m1_3) * dRhoAB / dRho;
    const double term2 = - LYP_COEF * LYP_PARAM_AB * dOmega * dRhoAB * (2.0 * dPow_rhoA_8_3);

    const double dAnswer = term1 + term2
                           + 2.0 * dRoundLYP_roundGammaAA * dGammaAA + dRoundLYP_roundGammaAB * dGammaAA;

    return dAnswer;
}

double DfFunctional_LYP::omega(const double dRho)
{
    const double dRho_for_m1_3 = std::pow(dRho, - INV_3);
    const double dRho_for_m11_3 = std::pow(dRho, - M_11_3);

    const double term1 = std::exp(- LYP_PARAM_C * dRho_for_m1_3);
    const double term2 = 1.0 + LYP_PARAM_D * dRho_for_m1_3;

    const double dAnswer =  term1 * dRho_for_m11_3 / term2;

    return dAnswer;
}

double DfFunctional_LYP::delta(const double dRho)
{
    const double dRho_for_m1_3 = std::pow(dRho, - INV_3);

    const double dAnswer = LYP_PARAM_C * dRho_for_m1_3 + LYP_PARAM_D * dRho_for_m1_3 / (1.0 + LYP_PARAM_D * dRho_for_m1_3);

    return dAnswer;
}

void DfFunctional_LYP::roundLYP_roundGamma(const double dRho, const double dRhoA, const double dRhoB,
                                           const double dOmega, const double dDelta,
                                           double& dRoundLYP_roundGammaAA, double& dRoundLYP_roundGammaAB, double& dRoundLYP_roundGammaBB)
{
    const double dRhoAB = dRhoA * dRhoB;
    const double ABOmega = LYP_PARAM_AB * dOmega;

    dRoundLYP_roundGammaAA = - ABOmega * (INV_9 * dRhoAB * (1.0 - 3.0 * dDelta - (dDelta - 11.0) * dRhoA / dRho) - dRhoB * dRhoB);
    dRoundLYP_roundGammaBB = - ABOmega * (INV_9 * dRhoAB * (1.0 - 3.0 * dDelta - (dDelta - 11.0) * dRhoB / dRho) - dRhoA * dRhoA);
    dRoundLYP_roundGammaAB = - ABOmega * (INV_9 * dRhoAB * (47.0 - 7.0 * dDelta) - M_4_3 * dRho * dRho);
}

// for RKS
void DfFunctional_LYP::roundLYP_roundGamma(const double dRho, const double dRhoA,
                                           const double dOmega, const double dDelta,
                                           double& dRoundLYP_roundGammaAA, double& dRoundLYP_roundGammaAB)
{
    const double dRhoAB = dRhoA * dRhoA;
    const double ABOmega = LYP_PARAM_AB * dOmega;

    dRoundLYP_roundGammaAA = - ABOmega * (INV_9 * dRhoAB * (1.0 - 3.0 * dDelta - (dDelta - 11.0) * dRhoA / dRho) - dRhoA * dRhoA);
    //dRoundLYP_roundGammaBB = - ABOmega * (INV_9 * dRhoAB * (1.0 - 3.0 * dDelta - (dDelta - 11.0) * dRhoB / dRho) - dRhoA * dRhoA);
    dRoundLYP_roundGammaAB = - ABOmega * (INV_9 * dRhoAB * (47.0 - 7.0 * dDelta) - M_4_3 * dRho * dRho);
}

void DfFunctional_LYP::roundLYP_roundRhoAB(const double dRhoA, const double dRhoB,
                                           double& roundLYP_roundRhoA, double& roundLYP_roundRhoB)
{
    const double dRho = dRhoA + dRhoB;
    const double dInv_Rho = 1.0 / dRho;
    const double dInv_RhoA = 1.0 / dRhoA;
    const double dInv_RhoB = 1.0 / dRhoB;

    const double dPow_rho_inv3 = std::pow(dRho, INV_3); // rho^(1/3)
    const double dPow_rho_mInv3 = 1.0 / dPow_rho_inv3; // rho^(-1/3)
    const double dPow_rho_m4_3 = dInv_Rho * dPow_rho_mInv3; // rho^(-4/3)

    const double dPow_rhoA_8_3 = std::pow(dRhoA, 8.0 / 3.0); // rhoA^(8/3)
    const double dPow_rhoB_8_3 = std::pow(dRhoB, 8.0 / 3.0); // rhoB^(8/3)

    const double dRhoAB = dRhoA * dRhoB;
    const double dOmega = this->omega(dRho);
    const double dOmegaPrime = this->omega_prime(dRho, dOmega);
    const double dDelta = this->delta(dRho);
    const double dDeltaPrime = this->delta_prime(dRho, dDelta);

    // for debug
//   std::cout << std::endl;
//   std::cout << ">>>> DfFunctional_LYP::roundLYP_roundRhoAB()" << std::endl;
//   std::cout << TlUtils::format("dRhoA = %+e, dRhoB = %+e", dRhoA, dRhoB)
//      << std::endl;
//   std::cout << TlUtils::format("omega = %+e, omega' = %+e, delta = %+e, delta' = %+e",
//                 dOmega, dOmegaPrime, dDelta, dDeltaPrime)
//      << std::endl;

    double dRoundRoundLYP_roundRhoARoundGammaAA, dRoundRoundLYP_roundRhoBRoundGammaAA;
    double dRoundRoundLYP_roundRhoARoundGammaAB, dRoundRoundLYP_roundRhoBRoundGammaAB;
    double dRoundRoundLYP_roundRhoARoundGammaBB, dRoundRoundLYP_roundRhoBRoundGammaBB;
    this->roundRoundLYP_roundRhoRoundGamma(dRho, dRhoA, dRhoB, dInv_Rho, dRhoAB, dOmega, dOmegaPrime, dDelta, dDeltaPrime,
                                           &dRoundRoundLYP_roundRhoARoundGammaAA, &dRoundRoundLYP_roundRhoBRoundGammaAA,
                                           &dRoundRoundLYP_roundRhoARoundGammaAB, &dRoundRoundLYP_roundRhoBRoundGammaAB,
                                           &dRoundRoundLYP_roundRhoARoundGammaBB, &dRoundRoundLYP_roundRhoBRoundGammaBB);

    const double dArg1_coef = - 4.0 * LYP_PARAM_A / (1.0 + LYP_PARAM_D * dPow_rho_mInv3) * dRhoAB / dRho;
    const double dArg1A = dArg1_coef * (INV_3 * (LYP_PARAM_D * dPow_rho_m4_3) / (1.0 + LYP_PARAM_D * dPow_rho_mInv3) + dInv_RhoA - dInv_Rho);
    const double dArg1B = dArg1_coef * (INV_3 * (LYP_PARAM_D * dPow_rho_m4_3) / (1.0 + LYP_PARAM_D * dPow_rho_mInv3) + dInv_RhoB - dInv_Rho);

    const double dArg2A = - LYP_COEF * LYP_PARAM_AB * (dOmegaPrime * dRhoAB * (dPow_rhoA_8_3 + dPow_rhoB_8_3)
                                                       + dOmega * dRhoB * (M_11_3 * dPow_rhoA_8_3 + dPow_rhoB_8_3));
    const double dArg2B = - LYP_COEF * LYP_PARAM_AB * (dOmegaPrime * dRhoAB * (dPow_rhoB_8_3 + dPow_rhoA_8_3)
                                                       + dOmega * dRhoA * (M_11_3 * dPow_rhoB_8_3 + dPow_rhoA_8_3));

    roundLYP_roundRhoA = dArg1A + dArg2A
                         + dRoundRoundLYP_roundRhoARoundGammaAA + dRoundRoundLYP_roundRhoARoundGammaAB + dRoundRoundLYP_roundRhoARoundGammaBB;
    roundLYP_roundRhoB = dArg1B + dArg2B
                         + dRoundRoundLYP_roundRhoBRoundGammaAA + dRoundRoundLYP_roundRhoBRoundGammaAB + dRoundRoundLYP_roundRhoBRoundGammaBB;

    // for debug
//   std::cout << TlUtils::format("ddLYP/dRAdGAA = %+e, ddLYP/dRAdGAB = %+e, ddLYP/dRAdGBB = %+e",
//                 dRoundRoundLYP_roundRhoARoundGammaAA, dRoundRoundLYP_roundRhoARoundGammaAB, dRoundRoundLYP_roundRhoARoundGammaBB)
//      << std::endl;
//   std::cout << TlUtils::format("ddLYP/dRBdGAA = %+e, ddLYP/dRBdGAB = %+e, ddLYP/dRBdGBB = %+e",
//                 dRoundRoundLYP_roundRhoBRoundGammaAA, dRoundRoundLYP_roundRhoBRoundGammaAB, dRoundRoundLYP_roundRhoBRoundGammaBB)
//      << std::endl;
//   std::cout << TlUtils::format("dLYP/dRhoA = %+e, dLYP/dRhoB = %+e", roundLYP_roundRhoA, roundLYP_roundRhoB)
//      << std::endl;
}

// gamma derivatives
void DfFunctional_LYP::roundRoundLYP_roundRhoRoundGamma(const double dRho, const double dRhoA, const double dRhoB,
                                                        const double dInvRho, const double dRhoAB,
                                                        const double dOmega, const double dOmegaPrime,
                                                        const double dDelta, const double dDeltaPrime,
                                                        double* pRoundRoundLYP_roundRhoARoundGammaAA, double* pRoundRoundLYP_roundRhoBRoundGammaAA,
                                                        double* pRoundRoundLYP_roundRhoARoundGammaAB, double* pRoundRoundLYP_roundRhoBRoundGammaAB,
                                                        double* pRoundRoundLYP_roundRhoARoundGammaBB, double* pRoundRoundLYP_roundRhoBRoundGammaBB)
{
    const double ABOmega = LYP_PARAM_AB * dOmega;

    const double dOmegaPrime_Omega = - INV_3 * std::pow(dRho, - M_4_3) * (11.0 * std::pow(dRho, INV_3)
                                                                          - LYP_PARAM_C - LYP_PARAM_D / (1.0 + LYP_PARAM_D * std::pow(dRho, -INV_3)));

    double dRoundLYP_roundGammaAA, dRoundLYP_roundGammaAB, dRoundLYP_roundGammaBB;
    this->roundLYP_roundGamma(dRho, dRhoA, dRhoB, dOmega, dDelta,
                              dRoundLYP_roundGammaAA, dRoundLYP_roundGammaAB, dRoundLYP_roundGammaBB);

    // dRoundRoundLYP_roundRhoXRoundGammaAA
    {
        // for alpha spin
        const double term1A = dOmegaPrime_Omega * dRoundLYP_roundGammaAA;
        const double term2A_1 =   INV_9 * dRhoB * (1.0 - 3.0 * dDelta - (dDelta - 11.0) * dRhoA * dInvRho);
        const double term2A_2 = - INV_9 * dRhoAB * ((3.0 + dRhoA * dInvRho) * dDeltaPrime
                                                    + (dDelta - 11.0) * dRhoB * dInvRho * dInvRho);
        const double term2A = - ABOmega * (term2A_1 + term2A_2);

        // for beta spin
        const double term1B = dOmegaPrime_Omega * dRoundLYP_roundGammaBB;
        const double term2B_1 =   INV_9 * dRhoA * (1.0 - 3.0 * dDelta - (dDelta - 11.0) * dRhoB * dInvRho);
        const double term2B_2 = - INV_9 * dRhoAB * ((3.0 + dRhoB * dInvRho) * dDeltaPrime
                                                    + (dDelta - 11.0) * dRhoA * dInvRho * dInvRho);
        const double term2B = - ABOmega * (term2B_1 + term2B_2);

        *pRoundRoundLYP_roundRhoARoundGammaAA = term1A + term2A;
        *pRoundRoundLYP_roundRhoBRoundGammaBB = term1B + term2B;
    }

    // roundRoundLYP_roundRhoARoundGammaAB
    {
        const double term1 = dOmegaPrime_Omega * dRoundLYP_roundGammaAB;

        const double term2_1 = INV_9 * (47.0 - 7.0 * dDelta);
        const double term2_2 = - M_7_9 * dRhoAB * dDeltaPrime - M_8_3 * dRho;
        const double term2A = - ABOmega * (term2_1 * dRhoB + term2_2);
        const double term2B = - ABOmega * (term2_1 * dRhoA + term2_2);

        *pRoundRoundLYP_roundRhoARoundGammaAB = term1 + term2A;
        *pRoundRoundLYP_roundRhoBRoundGammaAB = term1 + term2B;
    }

    // dRoundRoundLYP_roundRhoXRoundGammaBB
    {
        const double term1A = dOmegaPrime_Omega * dRoundLYP_roundGammaBB;
        const double term2A_1 =   INV_9 * dRhoB * (1.0 - 3.0 * dDelta - (dDelta - 11.0) * dRhoB * dInvRho);
        const double term2A_2 = - INV_9 * dRhoAB * ((3.0 + dRhoB * dInvRho) * dDeltaPrime
                                                    - (dDelta - 11.0) * dRhoB * dInvRho * dInvRho);
        const double term2A_3 = - 2.0 * dRhoA;
        const double term2A = - ABOmega * (term2A_1 + term2A_2 + term2A_3);

        const double term1B = dOmegaPrime_Omega * dRoundLYP_roundGammaAA;
        const double term2B_1 = INV_9 * dRhoA * (1.0 - 3.0 * dDelta - (dDelta - 11.0) * dRhoA * dInvRho);
        const double term2B_2 = - INV_9 * dRhoAB * ((3.0 + dRhoA * dInvRho) * dDeltaPrime
                                                    - (dDelta - 11.0) * dRhoA * dInvRho * dInvRho);
        const double term2B_3 = - 2.0 * dRhoB;
        const double term2B = - ABOmega * (term2B_1 + term2B_2 + term2B_3);

        *pRoundRoundLYP_roundRhoARoundGammaBB = term1A + term2A;
        *pRoundRoundLYP_roundRhoBRoundGammaAA = term1B + term2B;
    }
}

// for RKS
void DfFunctional_LYP::roundRoundLYP_roundRhoRoundGamma(const double dRho, const double dRhoA,
                                                        const double dInvRho, const double dRhoAB,
                                                        const double dOmega, const double dOmegaPrime,
                                                        const double dDelta, const double dDeltaPrime,
                                                        double* pRoundRoundLYP_roundRhoARoundGammaAA,
                                                        double* pRoundRoundLYP_roundRhoARoundGammaAB)
{
    const double ABOmega = LYP_PARAM_AB * dOmega;

    const double dOmegaPrime_Omega = - INV_3 * std::pow(dRho, - M_4_3) * (11.0 * std::pow(dRho, INV_3)
                                                                          - LYP_PARAM_C - LYP_PARAM_D / (1.0 + LYP_PARAM_D * std::pow(dRho, -INV_3)));

    double dRoundLYP_roundGammaAA, dRoundLYP_roundGammaAB;
    this->roundLYP_roundGamma(dRho, dRhoA, dOmega, dDelta,
                              dRoundLYP_roundGammaAA, dRoundLYP_roundGammaAB);

    // dRoundRoundLYP_roundRhoXRoundGammaAA
    {
        const double term1 = dOmegaPrime_Omega * dRoundLYP_roundGammaAA;

        // for alpha spin
        const double term2A_1 =   INV_9 * dRhoA * (1.0 - 3.0 * dDelta - (dDelta - 11.0) * dRhoA * dInvRho);
        const double term2A_2 = - INV_9 * dRhoAB * ((3.0 + dRhoA * dInvRho) * dDeltaPrime
                                                    + (dDelta - 11.0) * dRhoA * dInvRho * dInvRho);
        const double term2A = - ABOmega * (term2A_1 + term2A_2);

        // for beta spin
//     const double term2B_1 =   INV_9 * dRhoA * (1.0 - 3.0 * dDelta - (dDelta - 11.0) * dRhoB * dInvRho);
//     const double term2B_2 = - INV_9 * dRhoAB * ((3.0 + dRhoB * dInvRho) * dDeltaPrime
//                       + (dDelta - 11.0) * dRhoA * dInvRho * dInvRho);
//     const double term2B = - ABOmega * (term2B_1 + term2B_2);

        *pRoundRoundLYP_roundRhoARoundGammaAA = term1 + term2A;
//     *pRoundRoundLYP_roundRhoBRoundGammaBB = term1 + term2B;
    }

    // roundRoundLYP_roundRhoARoundGammaAB
    {
        const double term1 = dOmegaPrime_Omega * dRoundLYP_roundGammaAB;

        const double term2_1 = INV_9 * (47.0 - 7.0 * dDelta);
        const double term2_2 = - M_7_9 * dRhoAB * dDeltaPrime - M_8_3 * dRho;
        const double term2A = - ABOmega * (term2_1 * dRhoA + term2_2);
//     const double term2B = - ABOmega * (term2_1 * dRhoA + term2_2);

        *pRoundRoundLYP_roundRhoARoundGammaAB = term1 + term2A;
//     *pRoundRoundLYP_roundRhoBRoundGammaAB = term1 + term2B;
    }
}

double DfFunctional_LYP::omega_prime(const double dRho, const double dOmega)
{
    const double dRho_for_1_3 = std::pow(dRho, INV_3); // rho~(1/3)
    const double dRho_for_m_1_3 = 1.0 / dRho_for_1_3;  // rho~(-1/3)
    const double dRho_for_m_4_3 = (1.0 / dRho) * dRho_for_m_1_3; // rho~(-4/3)

    const double dAnswer = - INV_3 * dRho_for_m_4_3 * dOmega * (11.0 * dRho_for_1_3
                                                                - LYP_PARAM_C
                                                                - (LYP_PARAM_D / (1.0 + LYP_PARAM_D * dRho_for_m_1_3)));

    return dAnswer;
}

double DfFunctional_LYP::delta_prime(const double dRho, const double dDelta)
{
    const double dArg1_1 = 1.0 + LYP_PARAM_D * std::pow(dRho, - INV_3);
    const double dArg1 = (LYP_PARAM_DD * std::pow(dRho, - M_5_3)) / (dArg1_1 * dArg1_1);
    const double dArg2 = - 1.0 / dRho * dDelta;

    const double dAnswer = INV_3 * (dArg1 + dArg2);

    return dAnswer;
}

