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
