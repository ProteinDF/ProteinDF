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

#include "DfFunctional_B3LYP.h"
#include "TlUtils.h"

// HF 交換項の係数は0.2
const double DfFunctional_B3LYP::LDA_COEF = 0.80;
const double DfFunctional_B3LYP::B88_COEF = 0.72;
const double DfFunctional_B3LYP::VWN_COEF = 0.19;
const double DfFunctional_B3LYP::LYP_COEF = 0.81;

DfFunctional_B3LYP::DfFunctional_B3LYP() {
  const index_type numOfFunctionalTerms_LDA =
      this->m_LDA.getNumOfFunctionalTerms();
  const index_type numOfFunctionalTerms_B88 =
      this->m_B88.getNumOfFunctionalTerms();
  const index_type numOfFunctionalTerms_VWN =
      this->m_VWN.getNumOfFunctionalTerms();
  const index_type numOfFunctionalTerms_LYP =
      this->m_LYP.getNumOfFunctionalTerms();
  this->numOfFunctionalTerms_ =
      numOfFunctionalTerms_LDA + numOfFunctionalTerms_B88 +
      numOfFunctionalTerms_VWN + numOfFunctionalTerms_LYP;

  const index_type numOfDerivativeFunctionalTerms_LDA =
      this->m_LDA.getNumOfDerivativeFunctionalTerms();
  const index_type numOfDerivativeFunctionalTerms_B88 =
      this->m_B88.getNumOfDerivativeFunctionalTerms();
  const index_type numOfDerivativeFunctionalTerms_VWN =
      this->m_VWN.getNumOfDerivativeFunctionalTerms();
  const index_type numOfDerivativeFunctionalTerms_LYP =
      this->m_LYP.getNumOfDerivativeFunctionalTerms();
  this->numOfDerivativeFunctionalTerms_ =
      numOfDerivativeFunctionalTerms_LDA + numOfDerivativeFunctionalTerms_B88 +
      numOfDerivativeFunctionalTerms_VWN + numOfDerivativeFunctionalTerms_LYP;
}

DfFunctional_B3LYP::~DfFunctional_B3LYP() {}

double DfFunctional_B3LYP::getFunctional(const double dRhoA, const double dRhoB,
                                         const double dGammaAA,
                                         const double dGammaAB,
                                         const double dGammaBB) {
  const double ex_lda = this->m_LDA.getFunctional(dRhoA, dRhoB);
  const double ex_b88 =
      this->m_B88.getFunctional(dRhoA, dRhoB, dGammaAA, dGammaAB, dGammaBB);

  const double ec_vwn = this->m_VWN.getFunctional(dRhoA, dRhoB);
  const double ec_lyp =
      this->m_LYP.getFunctional(dRhoA, dRhoB, dGammaAA, dGammaAB, dGammaBB);

  const double answer = LDA_COEF * ex_lda + B88_COEF * ex_b88 +
                        VWN_COEF * ec_vwn + LYP_COEF * ec_lyp;

  return answer;
}

double DfFunctional_B3LYP::getFunctional(const double dRhoA,
                                         const double dGammaAA) {
  const double ex_lda = this->m_LDA.getFunctional(dRhoA);
  const double ex_b88 = this->m_B88.getFunctional(dRhoA, dGammaAA);

  const double ec_vwn = this->m_VWN.getFunctional(dRhoA);
  const double ec_lyp = this->m_LYP.getFunctional(dRhoA, dGammaAA);

  const double answer = LDA_COEF * ex_lda + B88_COEF * ex_b88 +
                        VWN_COEF * ec_vwn + LYP_COEF * ec_lyp;

  return answer;
}

void DfFunctional_B3LYP::getDerivativeFunctional(
    const double dRhoA, const double dRhoB, const double dGammaAA,
    const double dGammaAB, const double dGammaBB, double* pRoundF_roundRhoA,
    double* pRoundF_roundRhoB, double* pRoundF_roundGammaAA,
    double* pRoundF_roundGammaAB, double* pRoundF_roundGammaBB) {
  assert(pRoundF_roundRhoA != NULL);
  assert(pRoundF_roundRhoB != NULL);
  assert(pRoundF_roundGammaAA != NULL);
  assert(pRoundF_roundGammaAB != NULL);
  assert(pRoundF_roundGammaBB != NULL);

  // LDA
  double dRoundF_roundRhoA_LDA, dRoundF_roundRhoB_LDA;
  double dRoundF_roundRhoA_VWN, dRoundF_roundRhoB_VWN;
  // GGA
  double dRoundF_roundRhoA_B88, dRoundF_roundRhoB_B88;
  double dRoundF_roundRhoA_LYP, dRoundF_roundRhoB_LYP;
  double dRoundF_roundGammaAA_B88, dRoundF_roundGammaAB_B88,
      dRoundF_roundGammaBB_B88;
  double dRoundF_roundGammaAA_LYP, dRoundF_roundGammaAB_LYP,
      dRoundF_roundGammaBB_LYP;

  this->m_LDA.getDerivativeFunctional(dRhoA, dRhoB, &dRoundF_roundRhoA_LDA,
                                      &dRoundF_roundRhoB_LDA);
  this->m_VWN.getDerivativeFunctional(dRhoA, dRhoB, &dRoundF_roundRhoA_VWN,
                                      &dRoundF_roundRhoB_VWN);

  this->m_B88.getDerivativeFunctional(
      dRhoA, dRhoB, dGammaAA, dGammaAB, dGammaBB, &dRoundF_roundRhoA_B88,
      &dRoundF_roundRhoB_B88, &dRoundF_roundGammaAA_B88,
      &dRoundF_roundGammaAB_B88, &dRoundF_roundGammaBB_B88);
  this->m_LYP.getDerivativeFunctional(
      dRhoA, dRhoB, dGammaAA, dGammaAB, dGammaBB, &dRoundF_roundRhoA_LYP,
      &dRoundF_roundRhoB_LYP, &dRoundF_roundGammaAA_LYP,
      &dRoundF_roundGammaAB_LYP, &dRoundF_roundGammaBB_LYP);

  *pRoundF_roundRhoA =
      LDA_COEF * dRoundF_roundRhoA_LDA + B88_COEF * dRoundF_roundRhoA_B88 +
      VWN_COEF * dRoundF_roundRhoA_VWN + LYP_COEF * dRoundF_roundRhoA_LYP;

  *pRoundF_roundRhoB =
      LDA_COEF * dRoundF_roundRhoB_LDA + B88_COEF * dRoundF_roundRhoB_B88 +
      VWN_COEF * dRoundF_roundRhoB_VWN + LYP_COEF * dRoundF_roundRhoB_LYP;

  *pRoundF_roundGammaAA =
      B88_COEF * dRoundF_roundGammaAA_B88 + LYP_COEF * dRoundF_roundGammaAA_LYP;
  *pRoundF_roundGammaAB =
      B88_COEF * dRoundF_roundGammaAB_B88 + LYP_COEF * dRoundF_roundGammaAB_LYP;
  *pRoundF_roundGammaBB =
      B88_COEF * dRoundF_roundGammaBB_B88 + LYP_COEF * dRoundF_roundGammaBB_LYP;
}

void DfFunctional_B3LYP::getDerivativeFunctional(const double dRhoA,
                                                 const double dGammaAA,
                                                 double* pRoundF_roundRhoA,
                                                 double* pRoundF_roundGammaAA,
                                                 double* pRoundF_roundGammaAB) {
  assert(pRoundF_roundRhoA != NULL);
  assert(pRoundF_roundGammaAA != NULL);
  assert(pRoundF_roundGammaAB != NULL);

  // LDA
  double dRoundF_roundRhoA_LDA;
  double dRoundF_roundRhoA_VWN;
  // GGA
  double dRoundF_roundRhoA_B88;
  double dRoundF_roundRhoA_LYP;
  double dRoundF_roundGammaAA_B88, dRoundF_roundGammaAB_B88;
  double dRoundF_roundGammaAA_LYP, dRoundF_roundGammaAB_LYP;

  this->m_LDA.getDerivativeFunctional(dRhoA, &dRoundF_roundRhoA_LDA);
  this->m_VWN.getDerivativeFunctional(dRhoA, &dRoundF_roundRhoA_VWN);

  this->m_B88.getDerivativeFunctional(dRhoA, dGammaAA, &dRoundF_roundRhoA_B88,
                                      &dRoundF_roundGammaAA_B88,
                                      &dRoundF_roundGammaAB_B88);
  this->m_LYP.getDerivativeFunctional(dRhoA, dGammaAA, &dRoundF_roundRhoA_LYP,
                                      &dRoundF_roundGammaAA_LYP,
                                      &dRoundF_roundGammaAB_LYP);

  *pRoundF_roundRhoA =
      LDA_COEF * dRoundF_roundRhoA_LDA + B88_COEF * dRoundF_roundRhoA_B88 +
      VWN_COEF * dRoundF_roundRhoA_VWN + LYP_COEF * dRoundF_roundRhoA_LYP;

  *pRoundF_roundGammaAA =
      B88_COEF * dRoundF_roundGammaAA_B88 + LYP_COEF * dRoundF_roundGammaAA_LYP;
  *pRoundF_roundGammaAB =
      B88_COEF * dRoundF_roundGammaAB_B88 + LYP_COEF * dRoundF_roundGammaAB_LYP;
}

TlMatrix DfFunctional_B3LYP::getFunctionalCore(const double rhoA,
                                               const double rhoB,
                                               const double xA,
                                               const double xB) {
  const TlMatrix x_lda = this->m_LDA.getFunctionalCore(rhoA, rhoB);
  const TlMatrix x_b88 = this->m_B88.getFunctionalCore(rhoA, rhoB, xA, xB);
  const TlMatrix c_vwn = this->m_VWN.getFunctionalCore(rhoA, rhoB);
  const TlMatrix c_lyp = this->m_LYP.getFunctionalCore(rhoA, rhoB, xA, xB);

  const index_type numOfFunctionalTerms_LDA =
      this->m_LDA.getNumOfFunctionalTerms();
  const index_type numOfFunctionalTerms_B88 =
      this->m_B88.getNumOfFunctionalTerms();
  const index_type numOfFunctionalTerms_VWN =
      this->m_VWN.getNumOfFunctionalTerms();
  const index_type numOfFunctionalTerms_LYP =
      this->m_LYP.getNumOfFunctionalTerms();
  assert((numOfFunctionalTerms_LDA + numOfFunctionalTerms_B88 +
          numOfFunctionalTerms_VWN + numOfFunctionalTerms_LYP) ==
         this->getNumOfFunctionalTerms());
  TlMatrix xc(F_DIM, this->getNumOfFunctionalTerms());
  for (index_type i = 0; i < F_DIM; ++i) {
    index_type base = 0;
    for (index_type j = 0; j < numOfFunctionalTerms_LDA; ++j) {
      xc(i, j) = DfFunctional_B3LYP::LDA_COEF * x_lda(i, j);
    }
    base += numOfFunctionalTerms_LDA;
    for (index_type j = 0; j < numOfFunctionalTerms_B88; ++j) {
      xc(i, base + j) = DfFunctional_B3LYP::B88_COEF * x_b88(i, j);
    }
    base += numOfFunctionalTerms_B88;
    for (index_type j = 0; j < numOfFunctionalTerms_VWN; ++j) {
      xc(i, base + j) = DfFunctional_B3LYP::VWN_COEF * c_vwn(i, j);
    }
    base += numOfFunctionalTerms_VWN;
    for (index_type j = 0; j < numOfFunctionalTerms_LYP; ++j) {
      xc(i, base + j) = DfFunctional_B3LYP::LYP_COEF * c_lyp(i, j);
    }
  }

  return xc;
}

TlMatrix DfFunctional_B3LYP::getDerivativeFunctionalCore(const double rhoA,
                                                         const double rhoB,
                                                         const double xA,
                                                         const double xB) {
  const TlMatrix x_lda = this->m_LDA.getDerivativeFunctionalCore(rhoA, rhoB);
  const TlMatrix x_b88 =
      this->m_B88.getDerivativeFunctionalCore(rhoA, rhoB, xA, xB);
  const TlMatrix c_vwn = this->m_VWN.getDerivativeFunctionalCore(rhoA, rhoB);
  const TlMatrix c_lyp =
      this->m_LYP.getDerivativeFunctionalCore(rhoA, rhoB, xA, xB);

  const index_type numOfDerivativeFunctionalTerms_LDA =
      this->m_LDA.getNumOfDerivativeFunctionalTerms();
  const index_type numOfDerivativeFunctionalTerms_B88 =
      this->m_B88.getNumOfDerivativeFunctionalTerms();
  const index_type numOfDerivativeFunctionalTerms_VWN =
      this->m_VWN.getNumOfDerivativeFunctionalTerms();
  const index_type numOfDerivativeFunctionalTerms_LYP =
      this->m_LYP.getNumOfDerivativeFunctionalTerms();
  assert((numOfDerivativeFunctionalTerms_LDA +
          numOfDerivativeFunctionalTerms_B88 +
          numOfDerivativeFunctionalTerms_VWN +
          numOfDerivativeFunctionalTerms_LYP) ==
         this->getNumOfDerivativeFunctionalTerms());
  TlMatrix xc(D_DIM, this->getNumOfDerivativeFunctionalTerms());
  for (index_type i = 0; i < D_DIM; ++i) {
    index_type base = 0;
    for (index_type j = 0; j < numOfDerivativeFunctionalTerms_LDA; ++j) {
      xc(i, j) = DfFunctional_B3LYP::LDA_COEF * x_lda(i, j);
    }
    base += numOfDerivativeFunctionalTerms_LDA;
    for (index_type j = 0; j < numOfDerivativeFunctionalTerms_B88; ++j) {
      xc(i, base + j) = DfFunctional_B3LYP::B88_COEF * x_b88(i, j);
    }
    base += numOfDerivativeFunctionalTerms_B88;
    for (index_type j = 0; j < numOfDerivativeFunctionalTerms_VWN; ++j) {
      xc(i, base + j) = DfFunctional_B3LYP::VWN_COEF * c_vwn(i, j);
    }
    base += numOfDerivativeFunctionalTerms_VWN;
    for (index_type j = 0; j < numOfDerivativeFunctionalTerms_LYP; ++j) {
      xc(i, base + j) = DfFunctional_B3LYP::LYP_COEF * c_lyp(i, j);
    }
  }

  return xc;
}
