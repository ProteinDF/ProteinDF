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

#include <cassert>
#include <cmath>

#include "DfFunctional_Slater.h"
#include "TlUtils.h"

const double DfFunctional_Slater::INV_3 = 1.0 / 3.0;
const double DfFunctional_Slater::M_4_3 = 4.0 / 3.0;

// const double DfFunctional_Slater::SLATER_ALPHA = 2.0 / 3.0;
const double DfFunctional_Slater::SLATER_COEF =
    -(9.0 / 4.0) * std::pow((3.0 / (4.0 * M_PI)), (1.0 / 3.0)) *
    (2.0 / 3.0);  // (2.0 / 3.0) means Slater's alpha
const double DfFunctional_Slater::ROUND_SLATER_ROUND_RHO_COEF =
    -3.0 * std::pow((3.0 / (4.0 * M_PI)), (1.0 / 3.0)) *
    (2.0 / 3.0);  // (2.0 / 3.0) means Slater's alpha

DfFunctional_Slater::DfFunctional_Slater() {
    this->numOfFunctionalTerms_ = 1;
    this->numOfDerivativeFunctionalTerms_ = 1;
}

DfFunctional_Slater::~DfFunctional_Slater() {}

// for UKS
double DfFunctional_Slater::getFunctional(const double dRhoA,
                                          const double dRhoB) {
    const double dRhoATo4_3 = std::pow(dRhoA, M_4_3);
    const double dRhoBTo4_3 = std::pow(dRhoB, M_4_3);

    const double dAnswer = SLATER_COEF * (dRhoATo4_3 + dRhoBTo4_3);

    return dAnswer;
}

// specialized for RKS
double DfFunctional_Slater::getFunctional(const double dRhoA) {
    const double dRhoATo4_3 = std::pow(dRhoA, M_4_3);

    const double dAnswer = SLATER_COEF * (2.0 * dRhoATo4_3);

    return dAnswer;
}

// for UKS
void DfFunctional_Slater::getDerivativeFunctional(const double dRhoA,
                                                  const double dRhoB,
                                                  double* pRoundF_roundRhoA,
                                                  double* pRoundF_roundRhoB) {
    assert(pRoundF_roundRhoA != NULL);
    assert(pRoundF_roundRhoB != NULL);

    const double dRhoA_for_1_3 = std::pow(dRhoA, INV_3);
    const double dRhoB_for_1_3 = std::pow(dRhoB, INV_3);

    *pRoundF_roundRhoA = ROUND_SLATER_ROUND_RHO_COEF * dRhoA_for_1_3;
    *pRoundF_roundRhoB = ROUND_SLATER_ROUND_RHO_COEF * dRhoB_for_1_3;
}

// specialized for RKS
void DfFunctional_Slater::getDerivativeFunctional(const double dRhoA,
                                                  double* pRoundF_roundRhoA) {
    assert(pRoundF_roundRhoA != NULL);

    const double dRhoA_for_1_3 = std::pow(dRhoA, INV_3);

    *pRoundF_roundRhoA = ROUND_SLATER_ROUND_RHO_COEF * dRhoA_for_1_3;
}

// -------
TlDenseGeneralMatrix_Lapack DfFunctional_Slater::getFunctionalCore(
    const double rhoA, const double rhoB) {
    TlDenseGeneralMatrix_Lapack answer(F_DIM, this->getNumOfFunctionalTerms());
    assert(this->getNumOfFunctionalTerms() == 1);

    const double rhoA43 = std::pow(rhoA, M_4_3);
    const double rhoB43 = std::pow(rhoB, M_4_3);

    const double FA_termR = SLATER_COEF * rhoA43;
    const double FA_termX = 1.0;

    const double FB_termR = SLATER_COEF * rhoB43;
    const double FB_termX = 1.0;

    answer.set(FA_R, 0, FA_termR);
    answer.set(FA_X, 0, FA_termX);
    answer.set(FB_R, 0, FB_termR);
    answer.set(FB_X, 0, FB_termX);

    return answer;
}

TlDenseGeneralMatrix_Lapack DfFunctional_Slater::getDerivativeFunctionalCore(
    const double rhoA, const double rhoB) {
    TlDenseGeneralMatrix_Lapack answer(
        D_DIM, this->getNumOfDerivativeFunctionalTerms());
    assert(this->getNumOfFunctionalTerms() == 1);

    // roundF_roundRho =========================================================
    const double roundF_roundRhoA_termR =
        ROUND_SLATER_ROUND_RHO_COEF * std::pow(rhoA, INV_3);
    const double roundF_roundRhoA_termX = 1.0;
    const double roundF_roundRhoB_termR =
        ROUND_SLATER_ROUND_RHO_COEF * std::pow(rhoB, INV_3);
    const double roundF_roundRhoB_termX = 1.0;
    answer.set(RA_R, 0, roundF_roundRhoA_termR);
    answer.set(RA_X, 0, roundF_roundRhoA_termX);
    answer.set(RB_R, 0, roundF_roundRhoB_termR);
    answer.set(RB_X, 0, roundF_roundRhoB_termX);

    return answer;
}
