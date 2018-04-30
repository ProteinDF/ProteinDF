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

#include "DfFunctional_SVWN.h"

DfFunctional_SVWN::DfFunctional_SVWN() {}

DfFunctional_SVWN::~DfFunctional_SVWN() {}

double DfFunctional_SVWN::getFunctional(double dRhoA, double dRhoB) {
  const double ex = this->m_Slater.getFunctional(dRhoA, dRhoB);
  const double ec = this->m_VWN.getFunctional(dRhoA, dRhoB);

  return ex + ec;
}

void DfFunctional_SVWN::getDerivativeFunctional(
    double dRhoA, double dRhoB, double* pRoundFunctional_roundRhoA,
    double* pRoundFunctional_roundRhoB) {
  double dTmpXA, dTmpXB, dTmpCA, dTmpCB;
  this->m_Slater.getDerivativeFunctional(dRhoA, dRhoB, &dTmpXA, &dTmpXB);
  this->m_VWN.getDerivativeFunctional(dRhoA, dRhoB, &dTmpCA, &dTmpCB);

  *pRoundFunctional_roundRhoA = dTmpXA + dTmpCA;
  *pRoundFunctional_roundRhoB = dTmpXB + dTmpCB;
}
