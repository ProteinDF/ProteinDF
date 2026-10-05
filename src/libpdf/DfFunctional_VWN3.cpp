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

#include "DfFunctional_VWN3.h"
#include <cassert>
#include <cmath>

const double DfFunctional_VWN3::VWN3_A_PARA = 0.0310907;
const double DfFunctional_VWN3::VWN3_B_PARA = 13.0720;
const double DfFunctional_VWN3::VWN3_C_PARA = 42.7198;
const double DfFunctional_VWN3::VWN3_X0_PARA = -0.409286;
const double DfFunctional_VWN3::VWN3_A_FERR = 0.01554535;
const double DfFunctional_VWN3::VWN3_B_FERR = 20.1231;
const double DfFunctional_VWN3::VWN3_C_FERR = 101.578;
const double DfFunctional_VWN3::VWN3_X0_FERR = -0.743294;

// 1 / (2 * (2^(1/3) - 1)), the normalization of the von Barth-Hedin f(zeta)
const double DfFunctional_VWN3::F_ZETA_COEF =
    1.0 / (2.0 * (pow(2.0, 1.0 / 3.0) - 1.0));

DfFunctional_VWN3::DfFunctional_VWN3() {}

DfFunctional_VWN3::~DfFunctional_VWN3() {}

double DfFunctional_VWN3::epsilonC_PARA(const double x) {
    return this->epsilonC(VWN3_A_PARA, VWN3_B_PARA, VWN3_C_PARA, VWN3_X0_PARA,
                          x);
}

double DfFunctional_VWN3::epsilonC_FERR(const double x) {
    return this->epsilonC(VWN3_A_FERR, VWN3_B_FERR, VWN3_C_FERR, VWN3_X0_FERR,
                          x);
}

double DfFunctional_VWN3::epsilonCPrime_PARA(const double x) {
    return this->epsilonCPrime(VWN3_A_PARA, VWN3_B_PARA, VWN3_C_PARA,
                               VWN3_X0_PARA, x);
}

double DfFunctional_VWN3::epsilonCPrime_FERR(const double x) {
    return this->epsilonCPrime(VWN3_A_FERR, VWN3_B_FERR, VWN3_C_FERR,
                               VWN3_X0_FERR, x);
}

// von Barth-Hedin interpolation function
// f(zeta) = [(1+zeta)^(4/3) + (1-zeta)^(4/3) - 2] / (2 * (2^(1/3) - 1))
double DfFunctional_VWN3::f_zeta(const double zeta) {
    const double term1 = pow(1.0 + zeta, 4.0 / 3.0);
    const double term2 = pow(1.0 - zeta, 4.0 / 3.0);

    return F_ZETA_COEF * (term1 + term2 - 2.0);
}

// f'(zeta) = (4/3) * [(1+zeta)^(1/3) - (1-zeta)^(1/3)] / (2 * (2^(1/3) - 1))
double DfFunctional_VWN3::f_zeta_prime(const double zeta) {
    const double term1 = pow(1.0 + zeta, 1.0 / 3.0);
    const double term2 = pow(1.0 - zeta, 1.0 / 3.0);

    return F_ZETA_COEF * (4.0 / 3.0) * (term1 - term2);
}

// RPA paramagnetic/ferromagnetic energies with the von Barth-Hedin linear
// interpolation (no spin-stiffness term), used instead of VWN::epsilonC's
// eq.A11.
double DfFunctional_VWN3::epsilonC(const double x, const double zeta) {
    const double dEc_p = this->epsilonC_PARA(x);
    const double dEc_f = this->epsilonC_FERR(x);
    const double f = this->f_zeta(zeta);

    return dEc_p + f * (dEc_f - dEc_p);
}

void DfFunctional_VWN3::roundVWN_roundRho(const double dRhoA,
                                          const double dRhoB,
                                          double* pRoundF_roundRhoA,
                                          double* pRoundF_roundRhoB) {
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

    const double dEc_p = this->epsilonC_PARA(x);
    const double dEc_f = this->epsilonC_FERR(x);
    const double dEc_p_prime = this->epsilonCPrime_PARA(x);
    const double dEc_f_prime = this->epsilonCPrime_FERR(x);

    const double f = this->f_zeta(zeta);
    const double f_prime = this->f_zeta_prime(zeta);

    // dEC/dx (zeta fixed)
    const double term1 =
        -x / (6.0 * dRho) * (dEc_p_prime + f * (dEc_f_prime - dEc_p_prime));

    // dEC/dzeta
    const double term2coef = f_prime * (dEc_f - dEc_p);

    const double dRoundZeta_roundRhoA = dInvRho * (1.0 - zeta);
    const double dRoundZeta_roundRhoB = -dInvRho * (1.0 + zeta);
    const double dRoundEC_roundRhoA = term1 + term2coef * dRoundZeta_roundRhoA;
    const double dRoundEC_roundRhoB = term1 + term2coef * dRoundZeta_roundRhoB;

    if (dRhoA > TOLERANCE) {
        *pRoundF_roundRhoA = EC + dRho * dRoundEC_roundRhoA;
    }
    if (dRhoB > TOLERANCE) {
        *pRoundF_roundRhoB = EC + dRho * dRoundEC_roundRhoB;
    }
}
