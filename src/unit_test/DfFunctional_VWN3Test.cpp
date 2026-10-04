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
#include "gtest/gtest.h"

static const double EPS = 1.0E-10;
static const double NUM_DERIV_EPS = 1.0E-5;

// Reference values for the spin-polarized correlation energy/potential were
// obtained from libxc's LDA_C_VWN_RPA (the correlation part PySCF's B3LYP
// uses), via:
//   from pyscf.dft import libxc
//   exc, vxc, fxc, kxc = libxc.eval_xc(
//       'LDA_C_VWN_RPA', (rhoa, rhob), spin=1, deriv=1)
//   zk = exc * (rhoa + rhob); vrhoa, vrhob = vxc[0]
TEST(DfFunctional_VWN3, pointwise_closed_shell) {
    // input
    const double dRhoA = 1.7;
    const double dRhoB = 1.7;

    // expected value (zeta = 0)
    const double zk = -3.499718043733E-01;
    const double vRhoA = -1.121820560891E-01;
    const double vRhoB = -1.121820560891E-01;

    DfFunctional_VWN3 f;

    const double dFunctionalValue = f.getFunctional(dRhoA, dRhoB);
    EXPECT_NEAR(zk, dFunctionalValue, EPS);

    double dRoundF_roundRhoA, dRoundF_roundRhoB;
    f.getDerivativeFunctional(dRhoA, dRhoB, &dRoundF_roundRhoA,
                              &dRoundF_roundRhoB);
    EXPECT_NEAR(vRhoA, dRoundF_roundRhoA, EPS);
    EXPECT_NEAR(vRhoB, dRoundF_roundRhoB, EPS);
}

TEST(DfFunctional_VWN3, pointwise_partially_polarized) {
    // input
    const double dRhoA = 1.0;
    const double dRhoB = 0.5;

    // expected value (zeta = 1/3)
    const double zk = -1.381523240446E-01;
    const double vRhoA = -8.721853335312E-02;
    const double vRhoB = -1.277899671720E-01;

    DfFunctional_VWN3 f;

    const double dFunctionalValue = f.getFunctional(dRhoA, dRhoB);
    EXPECT_NEAR(zk, dFunctionalValue, EPS);

    double dRoundF_roundRhoA, dRoundF_roundRhoB;
    f.getDerivativeFunctional(dRhoA, dRhoB, &dRoundF_roundRhoA,
                              &dRoundF_roundRhoB);
    EXPECT_NEAR(vRhoA, dRoundF_roundRhoA, EPS);
    EXPECT_NEAR(vRhoB, dRoundF_roundRhoB, EPS);
}

TEST(DfFunctional_VWN3, pointwise_nearly_fully_polarized) {
    // input
    const double dRhoA = 1.9;
    const double dRhoB = 0.1;

    // expected value (zeta close to 1)
    const double zk = -1.406800947336E-01;
    const double vRhoA = -6.901931372758E-02;
    const double vRhoB = -2.122356397547E-01;

    DfFunctional_VWN3 f;

    const double dFunctionalValue = f.getFunctional(dRhoA, dRhoB);
    EXPECT_NEAR(zk, dFunctionalValue, EPS);

    double dRoundF_roundRhoA, dRoundF_roundRhoB;
    f.getDerivativeFunctional(dRhoA, dRhoB, &dRoundF_roundRhoA,
                              &dRoundF_roundRhoB);
    EXPECT_NEAR(vRhoA, dRoundF_roundRhoA, EPS);
    EXPECT_NEAR(vRhoB, dRoundF_roundRhoB, EPS);
}

// the analytical derivative must agree with the central-difference numerical
// derivative for an asymmetric, partially polarized point.
TEST(DfFunctional_VWN3, derivative_matches_numerical) {
    const double dRhoA = 0.8;
    const double dRhoB = 0.2;
    const double dDelta = 1.0E-6;

    DfFunctional_VWN3 f;

    double dRoundF_roundRhoA, dRoundF_roundRhoB;
    f.getDerivativeFunctional(dRhoA, dRhoB, &dRoundF_roundRhoA,
                              &dRoundF_roundRhoB);

    const double dNumRhoA = (f.getFunctional(dRhoA + dDelta, dRhoB) -
                             f.getFunctional(dRhoA - dDelta, dRhoB)) /
                            (2.0 * dDelta);
    const double dNumRhoB = (f.getFunctional(dRhoA, dRhoB + dDelta) -
                             f.getFunctional(dRhoA, dRhoB - dDelta)) /
                            (2.0 * dDelta);

    EXPECT_NEAR(dNumRhoA, dRoundF_roundRhoA, NUM_DERIV_EPS);
    EXPECT_NEAR(dNumRhoB, dRoundF_roundRhoB, NUM_DERIV_EPS);
}

// the RKS (closed-shell) path must still agree with the UKS path at zeta = 0,
// and VWN3's zeta-independent VWN5 base class must be unaffected (checked in
// DfFunctional_VWNTest.cpp).
TEST(DfFunctional_VWN3, rks_matches_uks_at_zeta_zero) {
    const double dRho = 1.3;

    DfFunctional_VWN3 f;

    const double dFunctionalValueRKS = f.getFunctional(dRho);
    const double dFunctionalValueUKS = f.getFunctional(dRho, dRho);
    EXPECT_NEAR(dFunctionalValueUKS, dFunctionalValueRKS, EPS);

    double dRoundF_roundRhoRKS;
    f.getDerivativeFunctional(dRho, &dRoundF_roundRhoRKS);

    double dRoundF_roundRhoA, dRoundF_roundRhoB;
    f.getDerivativeFunctional(dRho, dRho, &dRoundF_roundRhoA,
                              &dRoundF_roundRhoB);
    EXPECT_NEAR(dRoundF_roundRhoA, dRoundF_roundRhoRKS, EPS);
    EXPECT_NEAR(dRoundF_roundRhoB, dRoundF_roundRhoRKS, EPS);
}
