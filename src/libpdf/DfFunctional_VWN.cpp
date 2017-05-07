// Copyright (C) 2002-2014 The ProteinDF project
// see also AUTHORS and README.
// 
// This file is part of ProteinDF.
// 
// ProteinDF is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
// ProteinDF is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with ProteinDF.  If not, see <http://www.gnu.org/licenses/>.

#include <cmath>
#include <cassert>
#include "DfFunctional_VWN.h"

const double DfFunctional_VWN::TOLERANCE = 1.0E-20;

const double DfFunctional_VWN::M_3_4PI = 3.0 / (4.0 * M_PI);
const double DfFunctional_VWN::M_3_2 = 3.0 / 2.0;
const double DfFunctional_VWN::INV_3 = 1.0 / 3.0;
const double DfFunctional_VWN::INV_6 = 1.0 / 6.0;
const double DfFunctional_VWN::EC_COEF = 4.0 / (9.0 * (pow(2.0, 1.0/3.0) -1.0)); // epsilon_c の係数
const double DfFunctional_VWN::M_4_3 = 4.0 / 3.0;
const double DfFunctional_VWN::M_9_8 = 9.0 / 8.0;

const double DfFunctional_VWN::VWN5_A_PARA  =   0.0310907;
const double DfFunctional_VWN::VWN5_B_PARA  =   3.72744;
const double DfFunctional_VWN::VWN5_C_PARA  =  12.9352;
const double DfFunctional_VWN::VWN5_X0_PARA = - 0.10498;
const double DfFunctional_VWN::VWN5_A_FERR  =   0.01554535;
const double DfFunctional_VWN::VWN5_B_FERR  =   7.06042;
const double DfFunctional_VWN::VWN5_C_FERR  =  18.0578;
const double DfFunctional_VWN::VWN5_X0_FERR = - 0.325;

const double DfFunctional_VWN::VWN_A_ANTI  =   -1.0 / (6.0 * M_PI * M_PI);
const double DfFunctional_VWN::VWN_B_ANTI  =   1.13107;
const double DfFunctional_VWN::VWN_C_ANTI  =  13.0045;
const double DfFunctional_VWN::VWN_X0_ANTI = - 0.0047584;

DfFunctional_VWN::DfFunctional_VWN()
{
    this->numOfFunctionalTerms_ = 1;
    this->numOfDerivativeFunctionalTerms_ = 1;   
}

DfFunctional_VWN::~DfFunctional_VWN()
{
}

double DfFunctional_VWN::getFunctional(const double dRhoA, const double dRhoB)
{
    return this->VWN(dRhoA, dRhoB);
}

double DfFunctional_VWN::getFunctional(const double dRhoA)
{
    return this->VWN(dRhoA);
}

void DfFunctional_VWN::getDerivativeFunctional(const double dRhoA, const double dRhoB,
                                               double* pRoundF_roundRhoA, double* pRoundF_roundRhoB)
{
    this->roundVWN_roundRho(dRhoA, dRhoB, pRoundF_roundRhoA, pRoundF_roundRhoB);
}

void DfFunctional_VWN::getDerivativeFunctional(const double dRhoA, double* pRoundF_roundRhoA)
{
    this->roundVWN_roundRho(dRhoA, pRoundF_roundRhoA);
}

// for UKS
double DfFunctional_VWN::VWN(const double dRhoA, const double dRhoB)
{
    // (A10)
    const double dRho = dRhoA + dRhoB;
    const double dInvRho = 1.0 / dRho;

    const double x = pow(M_3_4PI * dInvRho, INV_6);
    const double zeta = (dRhoA - dRhoB) * dInvRho;

    // (A9)
    const double dAnswer = dRho * this->epsilonC(x, zeta);

    return dAnswer;
}

// specialized for RKS
double DfFunctional_VWN::VWN(const double dRhoA)
{
    // (A10)
    const double dRho = 2.0 * dRhoA;
    const double dInvRho = 1.0 / dRho;

    const double x = pow(M_3_4PI * dInvRho, INV_6);
    //const double zeta = 0.0; // (dRhoA - dRhoB) * dInvRho = 0.0 for RKS;

    // (A9)
    const double dAnswer = dRho * this->epsilonC(x); // call RKS version

    return dAnswer;
}

// VWN eq.A11 for UKS
double DfFunctional_VWN::epsilonC(const double x, const double zeta)
{
    const double dEc_p = this->epsilonC_PARA(x);
    const double dEc_f = this->epsilonC_FERR(x);
    const double dEc_a = this->epsilonC_ANTI(x);

    const double g = this->g(zeta);

    const double zeta4 = zeta * zeta * zeta * zeta; // zeta^4

    const double dArg2Coef = 1.0 + (EC_COEF * (dEc_f - dEc_p) / dEc_a - 1.0) * zeta4;

    const double dAnswer = dEc_p + dEc_a * g * dArg2Coef;

    //std::cout << "zeta=" << zeta << ", ECp=" << dEc_p << ", ECa=" << dEc_a << ", g=" << g << ", coef=" << dArg2Coef << std::endl;

    return dAnswer;
}

// VWN eq.A11 for RKS
double DfFunctional_VWN::epsilonC(const double x)
{
    const double dEc_p = this->epsilonC_PARA(x);
    //const double dEc_f = this->epsilonC_FERR(x);
    //const double dEc_a = this->epsilonC_ANTI(x);

    // g = 0, when rks calculation.
    //const double g = this->g(zeta);

    //const double zeta4 = zeta * zeta * zeta * zeta; // zeta^4

    //const double dArg2Coef = 1.0 + (EC_COEF * (dEc_f - dEc_p) / dEc_a - 1.0) * zeta4;

    const double dAnswer = dEc_p;

    return dAnswer;
}

double DfFunctional_VWN::epsilonC_PARA(const double x)
{
    return this->epsilonC(VWN5_A_PARA, VWN5_B_PARA, VWN5_C_PARA, VWN5_X0_PARA, x);
}

double DfFunctional_VWN::epsilonC_FERR(const double x)
{
    return this->epsilonC(VWN5_A_FERR, VWN5_B_FERR, VWN5_C_FERR, VWN5_X0_FERR, x);
}

double DfFunctional_VWN::epsilonC_ANTI(const double x)
{
    return this->epsilonC(VWN_A_ANTI, VWN_B_ANTI, VWN_C_ANTI, VWN_X0_ANTI, x);
}

// g(zeta)
// (A12)
double DfFunctional_VWN::g(const double zeta)
{
    const double dArg1 = pow((1.0 + zeta), M_4_3);
    const double dArg2 = pow((1.0 - zeta), M_4_3);

    const double dAnswer = M_9_8 * (dArg1 + dArg2 - 2.0);

    return dAnswer;
}

// (A13)
double DfFunctional_VWN::epsilonC(const double A, const double b, const double c, const double x0, const double x)
{
    // X(x) = x^2 + bx + c (A14)
    const double X = x * (x + b) + c;
    const double dInvX = 1.0 / X;
    const double X0 = x0 * (x0 + b) + c;
    const double dInvX0 = 1.0 / X0;

    // Q (A14)
    const double Q = sqrt(4 * c - b * b);
    const double dInvQ = 1.0 / Q;

    // (tan(Q / (2x+b)))^-1
    //const double dInvTan = (1.0 / tan(Q / (2.0 * x + b)));
    const double dInvTan = std::atan(Q / (2.0 * x + b));

    // epsilon_c(x)
    const double arg1 = std::log(x * x * dInvX);
    const double arg2 = 2.0 * b * dInvQ * dInvTan;

    const double x_x0 = x - x0;
    const double arg3_1 = std::log(x_x0 * x_x0 * dInvX);
    const double arg3_2 = 2.0 * (2.0 * x0 + b) * dInvQ * dInvTan;
    const double arg3 = - b * x0 * dInvX0 * (arg3_1 + arg3_2);

    const double dAnswer = A * (arg1 + arg2 + arg3);

    return dAnswer;
}

// h(x)
// (A15)
double DfFunctional_VWN::h(const double dEc_p, const double dEc_f, const double dEc_a)
{
    return EC_COEF *((dEc_f - dEc_p) / dEc_a) - 1.0;
}


// (A18)
double DfFunctional_VWN::epsilonCPrime(const double A, const double b, const double c, const double x0, const double x)
{
    // X(x) = x^2 + bx + c (A14)
    const double X = x * (x + b) + c;
    const double dInvX = 1.0 / X;
    const double X0 = x0 * (x0 + b) + c;
    const double dInvX0 = 1.0 / X0;

    // Q (A14)
    const double Q = sqrt(4 * c - b * b);
    //const double dInvQ = 1.0 / Q;

    // 2x+b
    const double x2_b = 2.0 * x + b;

    // 1.0 / ((2x +b)^2 + Q^2)
    const double inv_x2bQ = 1.0 / (x2_b * x2_b + Q * Q);

    const double dArg1 = 2.0 / x;
    const double dArg2 = - x2_b * dInvX;
    const double dArg3 = - 4.0 * b * inv_x2bQ;
    const double dArg4 = - b * x0 * dInvX0 * (2.0 / (x - x0) - x2_b * dInvX - 4.0 * (2.0 * x0 + b) * inv_x2bQ);

    const double dAnswer = A * (dArg1 + dArg2 + dArg3 + dArg4);
    return dAnswer;
}

// (A19)
double DfFunctional_VWN::h_prime(const double dEc_p, const double dEc_f, const double dEc_a,
                                 const double dEc_p_prime, const double dEc_f_prime, const double dEc_a_prime)
{
    const double dInvEc_a = 1.0 / dEc_a;

    const double dAnswer = EC_COEF * dInvEc_a * (dEc_f_prime - dEc_p_prime - (dEc_f - dEc_p) * dInvEc_a * dEc_a_prime);

    return dAnswer;
}

// (A20)
double DfFunctional_VWN::g_prime(const double zeta)
{
    const double term1 = pow((1.0 + zeta), INV_3);
    const double term2 = pow((1.0 - zeta), INV_3);

    return (M_3_2 *(term1 - term2));
}

void DfFunctional_VWN::roundVWN_roundRho(const double dRhoA, const double dRhoB,
                                         double* pRoundF_roundRhoA, double* pRoundF_roundRhoB)
{
    //assert(dRhoA >= dRhoB);
    assert(pRoundF_roundRhoA != NULL);
    assert(pRoundF_roundRhoB != NULL);

    // initialize
    *pRoundF_roundRhoA = 0.0;
    *pRoundF_roundRhoB = 0.0;

    // (A10)
    const double dRho = dRhoA + dRhoB;
    const double dInvRho = 1.0 / dRho;

    const double x = pow(M_3_4PI * dInvRho, INV_6);
    const double zeta = (dRhoA - dRhoB) * dInvRho;

    const double EC = this->epsilonC(x, zeta);

    // roundEC_roundRho
    const double dEc_p = this->epsilonC_PARA(x);
    const double dEc_f = this->epsilonC_FERR(x);
    const double dEc_a = this->epsilonC_ANTI(x);
    const double dEc_p_prime = this->epsilonCPrime_PARA(x);
    const double dEc_f_prime = this->epsilonCPrime_FERR(x);
    const double dEc_a_prime = this->epsilonCPrime_ANTI(x);

    const double g = this->g(zeta);
    const double g_prime = this->g_prime(zeta);
    const double h = this->h(dEc_p, dEc_f, dEc_a);
    const double h_prime = this->h_prime(dEc_p, dEc_f, dEc_a, dEc_p_prime, dEc_f_prime, dEc_a_prime);
    const double zeta3 = zeta * zeta * zeta; // zeta^3
    const double zeta4 = zeta3 * zeta; // zeta^4

    const double term1 = - x / (6.0 * dRho) * (dEc_p_prime
                                               + dEc_a_prime * g * (1.0 + h * zeta4)
                                               + dEc_a * g * h_prime * zeta4);

    const double term2coef = dEc_a * (g_prime * (1.0 + h * zeta4) + 4.0 * g * h * zeta3);
    //const double term2coef = dEc_a * (g_prime * (1.0 + h * zeta) + 4.0 * g * h * zeta);

    const double dRoundZeta_roundRhoA =   dInvRho * (1.0 - zeta);
    const double dRoundZeta_roundRhoB = - dInvRho * (1.0 + zeta);
    const double term2A = term2coef * dRoundZeta_roundRhoA;
    const double term2B = term2coef * dRoundZeta_roundRhoB;
    const double dRoundEC_roundRhoA = term1 + term2A;
    const double dRoundEC_roundRhoB = term1 + term2B;

    // roundVWN_roundRho
    if (dRhoA > TOLERANCE) {
        *pRoundF_roundRhoA = EC + dRho * dRoundEC_roundRhoA;
    }
    if (dRhoB > TOLERANCE) {
        *pRoundF_roundRhoB = EC + dRho * dRoundEC_roundRhoB;
    }

    // for debug
//   std::cout << "\n>>>>DfFunctional_VWN::roundVWN_roundRho()\n"
//      << TlUtils::format("dInvRho = %+e, zeta=%+e\n", dInvRho, zeta)
//      << TlUtils::format("(dEc_p,  dEc_f,  dEc_a ) = (%+e, %+e, %+e)\n", dEc_p, dEc_f, dEc_a)
//      << TlUtils::format("(dEc_p', dEc_f', dEc_a') = (%+e, %+e, %+e)\n", dEc_p_prime, dEc_f_prime, dEc_a_prime)
//      << TlUtils::format("(g, g') = %+e, %+e\n", g, g_prime)
//      << TlUtils::format("(h, h') = %+e, %+e\n", h, h_prime)
//      << TlUtils::format("term2coef = %+e\n", term2coef)
//      << TlUtils::format(" RoundZeta_roundRhoA = %+e\n", dRoundZeta_roundRhoA)
//      << TlUtils::format(" RoundEC_roundRhoA= %+e + %+e * %+e = %e\n", term1, term2coef, dRoundZeta_roundRhoA, dRoundEC_roundRhoA)
//      << TlUtils::format(" RoundF_roundRhoA = %+e + %+e * %+e = %+e", EC, dRho, dRoundEC_roundRhoA, *pRoundF_roundRhoA)
//      << "\n"
//      << TlUtils::format(" RoundZeta_roundRhoB = %+e\n", dRoundZeta_roundRhoB)
//      << TlUtils::format(" RoundEC_roundRhoB= %+e + %+e * %+e = %+e\n", term1, term2coef, dRoundZeta_roundRhoB, dRoundEC_roundRhoB)
//      << TlUtils::format(" RoundF_roundRhoB = %+e + %+e * %+e = %+e", EC, dRho, dRoundEC_roundRhoB, *pRoundF_roundRhoB) << std::endl;
}

// specialized for RKS
void DfFunctional_VWN::roundVWN_roundRho(const double dRhoA, double* pRoundF_roundRhoA)
{
    assert(pRoundF_roundRhoA != NULL);

    // initialize
    *pRoundF_roundRhoA = 0.0;

    // (A10)
    const double dRho = 2.0 * dRhoA;
    const double dInvRho = 1.0 / dRho;

    const double x = pow(M_3_4PI * dInvRho, INV_6);
    //const double zeta = (dRhoA - dRhoB) * dInvRho; // =0 when RKS.

    const double EC = this->epsilonC(x); // call RKS version

    // roundEC_roundRho
    //const double dEc_p = this->epsilonC_PARA(x);
    //const double dEc_f = this->epsilonC_FERR(x);
    //const double dEc_a = this->epsilonC_ANTI(x);
    const double dEc_p_prime = this->epsilonCPrime_PARA(x);
    //const double dEc_f_prime = this->epsilonCPrime_FERR(x);
    //const double dEc_a_prime = this->epsilonCPrime_ANTI(x);

    //const double g = this->g(zeta); // =0 when RKS.
    //const double g_prime = this->g_prime(zeta); // =0 when RKS.
    //const double h = this->h(dEc_p, dEc_f, dEc_a);
    //const double h_prime = this->h_prime(dEc_p, dEc_f, dEc_a, dEc_p_prime, dEc_f_prime, dEc_a_prime);
    //const double zeta3 = zeta * zeta * zeta; // zeta^3
    //const double zeta4 = zeta3 * zeta; // zeta^4

//   const double term1 = - x / (6.0 * dRho) * (dEc_p_prime
//                       + dEc_a_prime * g * (1.0 + h * zeta4)
//                       + dEc_a * g * h_prime * zeta4);
    //const double term1 = - x / (6.0 * dRho) * dEc_p_prime;
    const double term1 = - x * INV_6 * dInvRho * dEc_p_prime;

    //const double term2coef = dEc_a * (g_prime * (1.0 + h * zeta4) + 4.0 * g * h * zeta3); // =0 when RKS.

    //const double dRoundZeta_roundRhoA =   dInvRho * (1.0 - zeta);
    //const double dRoundZeta_roundRhoB = - dInvRho * (1.0 + zeta);
    //const double term2A = term2coef * dRoundZeta_roundRhoA;
    //const double term2B = term2coef * dRoundZeta_roundRhoB;
    //const double dRoundEC_roundRhoA = term1 + term2A;
    //const double dRoundEC_roundRhoB = term1 + term2B;
    const double dRoundEC_roundRhoA = term1;

    // roundVWN_roundRho
    if (dRhoA > TOLERANCE) {
        *pRoundF_roundRhoA = EC + dRho * dRoundEC_roundRhoA;
    }
}

double  DfFunctional_VWN::epsilonCPrime_PARA(const double x)
{
    return this->epsilonCPrime(VWN5_A_PARA, VWN5_B_PARA, VWN5_C_PARA, VWN5_X0_PARA, x);
}

double  DfFunctional_VWN::epsilonCPrime_FERR(const double x)
{
    return this->epsilonCPrime(VWN5_A_FERR, VWN5_B_FERR, VWN5_C_FERR, VWN5_X0_FERR, x);
}

double  DfFunctional_VWN::epsilonCPrime_ANTI(const double x)
{
    return this->epsilonCPrime(VWN_A_ANTI, VWN_B_ANTI, VWN_C_ANTI, VWN_X0_ANTI, x);
}

// ----------------
TlMatrix DfFunctional_VWN::getFunctionalCore(const double rhoA, 
                                             const double rhoB)
{
    TlMatrix answer(F_DIM, this->getNumOfFunctionalTerms());
    assert(this->getNumOfFunctionalTerms() == 1);

    // (A10)
    const double rho = rhoA + rhoB;
    const double invRho = 1.0 / rho;

    const double x = pow(M_3_4PI * invRho, INV_6);
    const double zeta = (rhoA - rhoB) * invRho;

    // (A9)
    const double epsilonC = this->epsilonC(x, zeta);
    // const double vwn = rho * epsilonC;

    const double FA_termR = rhoA * epsilonC;
    const double FA_termX = 1.0;
    const double FB_termR = rhoB * epsilonC;
    const double FB_termX = 1.0;
    answer.set(FA_R, 0, FA_termR);
    answer.set(FA_X, 0, FA_termX);
    answer.set(FB_R, 0, FB_termR);
    answer.set(FB_X, 0, FB_termX);

    return answer;
}

TlMatrix DfFunctional_VWN::getDerivativeFunctionalCore(const double rhoA,
                                                       const double rhoB)
{
    TlMatrix answer(D_DIM, this->getNumOfDerivativeFunctionalTerms());
    assert(this->getNumOfFunctionalTerms() == 1);

    // roundF_roundRho =========================================================
    double roundF_roundRhoA = 0.0;
    double roundF_roundRhoB = 0.0;
    this->roundVWN_roundRho(rhoA, rhoB, &roundF_roundRhoA, &roundF_roundRhoB);

    const double roundF_roundRhoA_termR = roundF_roundRhoA;
    const double roundF_roundRhoA_termX = 1.0;
    const double roundF_roundRhoB_termR = roundF_roundRhoB;
    const double roundF_roundRhoB_termX = 1.0;
    answer.set(RA_R, 0, roundF_roundRhoA_termR);
    answer.set(RA_X, 0, roundF_roundRhoA_termX);
    answer.set(RB_R, 0, roundF_roundRhoB_termR);
    answer.set(RB_X, 0, roundF_roundRhoB_termX);
    
    return answer;
}


