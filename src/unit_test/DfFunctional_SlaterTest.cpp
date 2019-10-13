#include <string>
#include <vector>
#include "DfFunctional_Slater.h"
#include "gtest/gtest.h"

const double EPS = 1.0E-8;

// =====================================================================
// data used from: http://www.cse.scitech.ac.uk/ccg/dft/data_pt_x_lda.html
//

// test1
//  rhoa= 0.17E+01 rhob= 0.17E+01 sigmaaa= 0.81E-11 sigmaab= 0.81E-11 sigmabb=
//  0.81E-11 zk            = -0.377592720836E+01 vrhoa         =
//  -0.148075576798E+01 vrhob         = -0.148075576798E+01 vsigmaaa      =
//  0.000000000000E+00 vsigmaab      =  0.000000000000E+00 vsigmabb      =
//  0.000000000000E+00

//  v2rhoa2       = -0.290344268232E+00
//  v2rhoab       =  0.000000000000E+00
//  v2rhob2       = -0.290344268232E+00

//  v2rhoasigmaaa =  0.000000000000E+00
//  v2rhoasigmaab =  0.000000000000E+00
//  v2rhoasigmabb =  0.000000000000E+00
//  v2rhobsigmaaa =  0.000000000000E+00
//  v2rhobsigmaab =  0.000000000000E+00
//  v2rhobsigmabb =  0.000000000000E+00

//  v2sigmaaa2    =  0.000000000000E+00
//  v2sigmaaaab   =  0.000000000000E+00
//  v2sigmaaabb   =  0.000000000000E+00
//  v2sigmaab2    =  0.000000000000E+00
//  v2sigmaabbb   =  0.000000000000E+00
//  v2sigmabb2    =  0.000000000000E+00
TEST(DfFunctional_Slater, pointwise1) {
    // input
    const double dRhoA = 0.17E+01;
    const double dRhoB = 0.17E+01;
    //   const double dGammaAA = 0.81E-11;
    //   const double dGammaAB = 0.81E-11;
    //   const double dGammaBB = 0.81E-11;

    // expected value
    const double zk = -0.377592720836E+01;
    const double vRhoA = -0.148075576798E+01;
    const double vRhoB = -0.148075576798E+01;
    //   const double vGammaAA = 0.000000000000E+00;
    //   const double vGammaAB = 0.000000000000E+00;
    //   const double vGammaBB = 0.000000000000E+00;

    // execute test
    DfFunctional_Slater f;

    double dFunctionalValue = f.getFunctional(dRhoA, dRhoB);
    EXPECT_NEAR(zk, dFunctionalValue, EPS);

    double dRoundF_roundRhoA, dRoundF_roundRhoB;
    f.getDerivativeFunctional(dRhoA, dRhoB, &dRoundF_roundRhoA,
                              &dRoundF_roundRhoB);
    EXPECT_NEAR(vRhoA, dRoundF_roundRhoA, EPS);
    EXPECT_NEAR(vRhoB, dRoundF_roundRhoB, EPS);
}

TEST(DfFunctional_Slater, pointwise1_RKS) {
    // input
    const double dRhoA = 0.17E+01;
    //   const double dRhoB = 0.17E+01;
    //   const double dGammaAA = 0.81E-11;
    //   const double dGammaAB = 0.81E-11;
    //   const double dGammaBB = 0.81E-11;

    // expected value
    const double zk = -0.377592720836E+01;
    const double vRhoA = -0.148075576798E+01;
    //   const double vRhoB = -0.148075576798E+01;
    //   const double vGammaAA = 0.000000000000E+00;
    //   const double vGammaAB = 0.000000000000E+00;
    //   const double vGammaBB = 0.000000000000E+00;

    // execute test
    DfFunctional_Slater f;

    double dFunctionalValue = f.getFunctional(dRhoA);
    EXPECT_NEAR(zk, dFunctionalValue, EPS);

    double dRoundF_roundRhoA;
    f.getDerivativeFunctional(dRhoA, &dRoundF_roundRhoA);
    EXPECT_NEAR(vRhoA, dRoundF_roundRhoA, EPS);
    //   EXPECT_NEAR(vRhoB, dRoundF_roundRhoB, EPS);
}

// test2
//  rhoa= 0.17E+01 rhob= 0.17E+01 sigmaaa= 0.17E+01 sigmaab= 0.17E+01 sigmabb=
//  0.17E+01 zk            = -0.377592720836E+01 vrhoa         =
//  -0.148075576798E+01 vrhob         = -0.148075576798E+01 vsigmaaa      =
//  0.000000000000E+00 vsigmaab      =  0.000000000000E+00 vsigmabb      =
//  0.000000000000E+00

//  v2rhoa2       = -0.290344268232E+00
//  v2rhoab       =  0.000000000000E+00
//  v2rhob2       = -0.290344268232E+00

//  v2rhoasigmaaa =  0.000000000000E+00
//  v2rhoasigmaab =  0.000000000000E+00
//  v2rhoasigmabb =  0.000000000000E+00
//  v2rhobsigmaaa =  0.000000000000E+00
//  v2rhobsigmaab =  0.000000000000E+00
//  v2rhobsigmabb =  0.000000000000E+00

//  v2sigmaaa2    =  0.000000000000E+00
//  v2sigmaaaab   =  0.000000000000E+00
//  v2sigmaaabb   =  0.000000000000E+00
//  v2sigmaab2    =  0.000000000000E+00
//  v2sigmaabbb   =  0.000000000000E+00
//  v2sigmabb2    =  0.000000000000E+00
TEST(DfFunctional_Slater, pointwise2) {
    // input
    const double dRhoA = 0.17E+01;
    const double dRhoB = 0.17E+01;

    // expected value
    const double zk = -0.377592720836E+01;
    const double vRhoA = -0.148075576798E+01;
    const double vRhoB = -0.148075576798E+01;

    // execute test
    DfFunctional_Slater f;

    double dFunctionalValue = f.getFunctional(dRhoA, dRhoB);
    EXPECT_NEAR(zk, dFunctionalValue, EPS);

    double dRoundF_roundRhoA, dRoundF_roundRhoB;
    f.getDerivativeFunctional(dRhoA, dRhoB, &dRoundF_roundRhoA,
                              &dRoundF_roundRhoB);
    EXPECT_NEAR(vRhoA, dRoundF_roundRhoA, EPS);
    EXPECT_NEAR(vRhoB, dRoundF_roundRhoB, EPS);
}

TEST(DfFunctional_Slater, pointwise2_RKS) {
    // input
    const double dRhoA = 0.17E+01;
    //   const double dRhoB = 0.17E+01;

    // expected value
    const double zk = -0.377592720836E+01;
    const double vRhoA = -0.148075576798E+01;
    //   const double vRhoB = -0.148075576798E+01;

    // execute test
    DfFunctional_Slater f;

    double dFunctionalValue = f.getFunctional(dRhoA);
    EXPECT_NEAR(zk, dFunctionalValue, EPS);

    double dRoundF_roundRhoA;
    f.getDerivativeFunctional(dRhoA, &dRoundF_roundRhoA);
    EXPECT_NEAR(vRhoA, dRoundF_roundRhoA, EPS);
    //   EXPECT_NEAR(vRhoB, dRoundF_roundRhoB, EPS);
}

// test3
//  rhoa= 0.15E+01 rhob= 0.15E+01 sigmaaa= 0.36E+02 sigmaab= 0.36E+02 sigmabb=
//  0.36E+02 zk            = -0.319555819038E+01 vrhoa         =
//  -0.142024808461E+01 vrhob         = -0.142024808461E+01 vsigmaaa      =
//  0.000000000000E+00 vsigmaab      =  0.000000000000E+00 vsigmabb      =
//  0.000000000000E+00

//  v2rhoa2       = -0.315610685470E+00
//  v2rhoab       =  0.000000000000E+00
//  v2rhob2       = -0.315610685470E+00

//  v2rhoasigmaaa =  0.000000000000E+00
//  v2rhoasigmaab =  0.000000000000E+00
//  v2rhoasigmabb =  0.000000000000E+00
//  v2rhobsigmaaa =  0.000000000000E+00
//  v2rhobsigmaab =  0.000000000000E+00
//  v2rhobsigmabb =  0.000000000000E+00

//  v2sigmaaa2    =  0.000000000000E+00
//  v2sigmaaaab   =  0.000000000000E+00
//  v2sigmaaabb   =  0.000000000000E+00
//  v2sigmaab2    =  0.000000000000E+00
//  v2sigmaabbb   =  0.000000000000E+00
//  v2sigmabb2    =  0.000000000000E+00
TEST(DfFunctional_Slater, pointwise3) {
    // input
    const double dRhoA = 0.15E+01;
    const double dRhoB = 0.15E+01;

    // expected value
    const double zk = -0.319555819038E+01;
    const double vRhoA = -0.142024808461E+01;
    const double vRhoB = -0.142024808461E+01;

    // execute test
    DfFunctional_Slater f;

    double dFunctionalValue = f.getFunctional(dRhoA, dRhoB);
    EXPECT_NEAR(zk, dFunctionalValue, EPS);

    double dRoundF_roundRhoA, dRoundF_roundRhoB;
    f.getDerivativeFunctional(dRhoA, dRhoB, &dRoundF_roundRhoA,
                              &dRoundF_roundRhoB);
    EXPECT_NEAR(vRhoA, dRoundF_roundRhoA, EPS);
    EXPECT_NEAR(vRhoB, dRoundF_roundRhoB, EPS);
}

TEST(DfFunctional_Slater, pointwise3_RKS) {
    // input
    const double dRhoA = 0.15E+01;
    //   const double dRhoB = 0.15E+01;

    // expected value
    const double zk = -0.319555819038E+01;
    const double vRhoA = -0.142024808461E+01;
    //   const double vRhoB = -0.142024808461E+01;

    // execute test
    DfFunctional_Slater f;

    double dFunctionalValue = f.getFunctional(dRhoA);
    EXPECT_NEAR(zk, dFunctionalValue, EPS);

    double dRoundF_roundRhoA;
    f.getDerivativeFunctional(dRhoA, &dRoundF_roundRhoA);
    EXPECT_NEAR(vRhoA, dRoundF_roundRhoA, EPS);
    //   EXPECT_NEAR(vRhoB, dRoundF_roundRhoB, EPS);
}

// test4
//  rhoa= 0.88E-01 rhob= 0.88E-01 sigmaaa= 0.87E-01 sigmaab= 0.87E-01 sigmabb=
//  0.87E-01 zk            = -0.728453690414E-01 vrhoa         =
//  -0.551858856374E+00 vrhob         = -0.551858856374E+00 vsigmaaa      =
//  0.000000000000E+00 vsigmaab      =  0.000000000000E+00 vsigmabb      =
//  0.000000000000E+00

//  v2rhoa2       = -0.209037445596E+01
//  v2rhoab       =  0.000000000000E+00
//  v2rhob2       = -0.209037445596E+01

//  v2rhoasigmaaa =  0.000000000000E+00
//  v2rhoasigmaab =  0.000000000000E+00
//  v2rhoasigmabb =  0.000000000000E+00
//  v2rhobsigmaaa =  0.000000000000E+00
//  v2rhobsigmaab =  0.000000000000E+00
//  v2rhobsigmabb =  0.000000000000E+00

//  v2sigmaaa2    =  0.000000000000E+00
//  v2sigmaaaab   =  0.000000000000E+00
//  v2sigmaaabb   =  0.000000000000E+00
//  v2sigmaab2    =  0.000000000000E+00
//  v2sigmaabbb   =  0.000000000000E+00
//  v2sigmabb2    =  0.000000000000E+00
TEST(DfFunctional_Slater, pointwise4) {
    // input
    const double dRhoA = 0.88E-01;
    const double dRhoB = 0.88E-01;

    // expected value
    const double zk = -0.728453690414E-01;
    const double vRhoA = -0.551858856374E+00;
    const double vRhoB = -0.551858856374E+00;

    // execute test
    DfFunctional_Slater f;

    double dFunctionalValue = f.getFunctional(dRhoA, dRhoB);
    EXPECT_NEAR(zk, dFunctionalValue, EPS);

    double dRoundF_roundRhoA, dRoundF_roundRhoB;
    f.getDerivativeFunctional(dRhoA, dRhoB, &dRoundF_roundRhoA,
                              &dRoundF_roundRhoB);
    EXPECT_NEAR(vRhoA, dRoundF_roundRhoA, EPS);
    EXPECT_NEAR(vRhoB, dRoundF_roundRhoB, EPS);
}

// test5
//  rhoa= 0.18E+04 rhob= 0.18E+04 sigmaaa= 0.55E+00 sigmaab= 0.55E+00 sigmabb=
//  0.55E+00 zk            = -0.407494475320E+05 vrhoa         =
//  -0.150923879748E+02 vrhob         = -0.150923879748E+02 vsigmaaa      =
//  0.000000000000E+00 vsigmaab      =  0.000000000000E+00 vsigmabb      =
//  0.000000000000E+00

//  v2rhoa2       = -0.279488666200E-02
//  v2rhoab       =  0.000000000000E+00
//  v2rhob2       = -0.279488666200E-02

//  v2rhoasigmaaa =  0.000000000000E+00
//  v2rhoasigmaab =  0.000000000000E+00
//  v2rhoasigmabb =  0.000000000000E+00
//  v2rhobsigmaaa =  0.000000000000E+00
//  v2rhobsigmaab =  0.000000000000E+00
//  v2rhobsigmabb =  0.000000000000E+00

//  v2sigmaaa2    =  0.000000000000E+00
//  v2sigmaaaab   =  0.000000000000E+00
//  v2sigmaaabb   =  0.000000000000E+00
//  v2sigmaab2    =  0.000000000000E+00
//  v2sigmaabbb   =  0.000000000000E+00
//  v2sigmabb2    =  0.000000000000E+00

// test6
//  rhoa= 0.18E+04 rhob= 0.18E+04 sigmaaa= 0.86E+04 sigmaab= 0.86E+04 sigmabb=
//  0.86E+04 zk            = -0.407494475320E+05 vrhoa         =
//  -0.150923879748E+02 vrhob         = -0.150923879748E+02 vsigmaaa      =
//  0.000000000000E+00 vsigmaab      =  0.000000000000E+00 vsigmabb      =
//  0.000000000000E+00

//  v2rhoa2       = -0.279488666200E-02
//  v2rhoab       =  0.000000000000E+00
//  v2rhob2       = -0.279488666200E-02

//  v2rhoasigmaaa =  0.000000000000E+00
//  v2rhoasigmaab =  0.000000000000E+00
//  v2rhoasigmabb =  0.000000000000E+00
//  v2rhobsigmaaa =  0.000000000000E+00
//  v2rhobsigmaab =  0.000000000000E+00
//  v2rhobsigmabb =  0.000000000000E+00

//  v2sigmaaa2    =  0.000000000000E+00
//  v2sigmaaaab   =  0.000000000000E+00
//  v2sigmaaabb   =  0.000000000000E+00
//  v2sigmaab2    =  0.000000000000E+00
//  v2sigmaabbb   =  0.000000000000E+00
//  v2sigmabb2    =  0.000000000000E+00

// test7
//  rhoa= 0.16E+04 rhob= 0.16E+04 sigmaaa= 0.37E+10 sigmaab= 0.37E+10 sigmabb=
//  0.37E+10 zk            = -0.348271841145E+05 vrhoa         =
//  -0.145113267144E+02 vrhob         = -0.145113267144E+02 vsigmaaa      =
//  0.000000000000E+00 vsigmaab      =  0.000000000000E+00 vsigmabb      =
//  0.000000000000E+00

//  v2rhoa2       = -0.302319306550E-02
//  v2rhoab       =  0.000000000000E+00
//  v2rhob2       = -0.302319306550E-02

//  v2rhoasigmaaa =  0.000000000000E+00
//  v2rhoasigmaab =  0.000000000000E+00
//  v2rhoasigmabb =  0.000000000000E+00
//  v2rhobsigmaaa =  0.000000000000E+00
//  v2rhobsigmaab =  0.000000000000E+00
//  v2rhobsigmabb =  0.000000000000E+00

//  v2sigmaaa2    =  0.000000000000E+00
//  v2sigmaaaab   =  0.000000000000E+00
//  v2sigmaaabb   =  0.000000000000E+00
//  v2sigmaab2    =  0.000000000000E+00
//  v2sigmaabbb   =  0.000000000000E+00
//  v2sigmabb2    =  0.000000000000E+00

// test8
//  rhoa= 0.26E+00 rhob= 0.26E+00 sigmaaa= 0.28E+00 sigmaab= 0.28E+00 sigmabb=
//  0.28E+00 zk            = -0.308832394647E+00 vrhoa         =
//  -0.791877934993E+00 vrhob         = -0.791877934993E+00 vsigmaaa      =
//  0.000000000000E+00 vsigmaab      =  0.000000000000E+00 vsigmabb      =
//  0.000000000000E+00

//  v2rhoa2       = -0.101522812179E+01
//  v2rhoab       =  0.000000000000E+00
//  v2rhob2       = -0.101522812179E+01

//  v2rhoasigmaaa =  0.000000000000E+00
//  v2rhoasigmaab =  0.000000000000E+00
//  v2rhoasigmabb =  0.000000000000E+00
//  v2rhobsigmaaa =  0.000000000000E+00
//  v2rhobsigmaab =  0.000000000000E+00
//  v2rhobsigmabb =  0.000000000000E+00

//  v2sigmaaa2    =  0.000000000000E+00
//  v2sigmaaaab   =  0.000000000000E+00
//  v2sigmaaabb   =  0.000000000000E+00
//  v2sigmaab2    =  0.000000000000E+00
//  v2sigmaabbb   =  0.000000000000E+00
//  v2sigmabb2    =  0.000000000000E+00

// test9
//  rhoa= 0.53E+05 rhob= 0.53E+05 sigmaaa= 0.96E+05 sigmaab= 0.96E+05 sigmabb=
//  0.96E+05 zk            = -0.370503980143E+07 vrhoa         =
//  -0.466042742318E+02 vrhob         = -0.466042742318E+02 vsigmaaa      =
//  0.000000000000E+00 vsigmaab      =  0.000000000000E+00 vsigmabb      =
//  0.000000000000E+00

//  v2rhoa2       = -0.293108642967E-03
//  v2rhoab       =  0.000000000000E+00
//  v2rhob2       = -0.293108642967E-03

//  v2rhoasigmaaa =  0.000000000000E+00
//  v2rhoasigmaab =  0.000000000000E+00
//  v2rhoasigmabb =  0.000000000000E+00
//  v2rhobsigmaaa =  0.000000000000E+00
//  v2rhobsigmaab =  0.000000000000E+00
//  v2rhobsigmabb =  0.000000000000E+00

//  v2sigmaaa2    =  0.000000000000E+00
//  v2sigmaaaab   =  0.000000000000E+00
//  v2sigmaaabb   =  0.000000000000E+00
//  v2sigmaab2    =  0.000000000000E+00
//  v2sigmaabbb   =  0.000000000000E+00
//  v2sigmabb2    =  0.000000000000E+00

// test10
//  rhoa= 0.47E+05 rhob= 0.47E+05 sigmaaa= 0.29E+14 sigmaab= 0.29E+14 sigmabb=
//  0.29E+14 zk            = -0.315661921284E+07 vrhoa         =
//  -0.447747406077E+02 vrhob         = -0.447747406077E+02 vsigmaaa      =
//  0.000000000000E+00 vsigmaab      =  0.000000000000E+00 vsigmabb      =
//  0.000000000000E+00

//  v2rhoa2       = -0.317551351828E-03
//  v2rhoab       =  0.000000000000E+00
//  v2rhob2       = -0.317551351828E-03

//  v2rhoasigmaaa =  0.000000000000E+00
//  v2rhoasigmaab =  0.000000000000E+00
//  v2rhoasigmabb =  0.000000000000E+00
//  v2rhobsigmaaa =  0.000000000000E+00
//  v2rhobsigmaab =  0.000000000000E+00
//  v2rhobsigmabb =  0.000000000000E+00

//  v2sigmaaa2    =  0.000000000000E+00
//  v2sigmaaaab   =  0.000000000000E+00
//  v2sigmaaabb   =  0.000000000000E+00
//  v2sigmaab2    =  0.000000000000E+00
//  v2sigmaabbb   =  0.000000000000E+00
//  v2sigmabb2    =  0.000000000000E+00

// test11
//  rhoa= 0.15E+00 rhob= 0.15E+00 sigmaaa= 0.16E+00 sigmaab= 0.16E+00 sigmabb=
//  0.16E+00 zk            = -0.148324672136E+00 vrhoa         =
//  -0.659220765051E+00 vrhob         = -0.659220765051E+00 vsigmaaa      =
//  0.000000000000E+00 vsigmaab      =  0.000000000000E+00 vsigmabb      =
//  0.000000000000E+00

//  v2rhoa2       = -0.146493503345E+01
//  v2rhoab       =  0.000000000000E+00
//  v2rhob2       = -0.146493503345E+01

//  v2rhoasigmaaa =  0.000000000000E+00
//  v2rhoasigmaab =  0.000000000000E+00
//  v2rhoasigmabb =  0.000000000000E+00
//  v2rhobsigmaaa =  0.000000000000E+00
//  v2rhobsigmaab =  0.000000000000E+00
//  v2rhobsigmabb =  0.000000000000E+00

//  v2sigmaaa2    =  0.000000000000E+00
//  v2sigmaaaab   =  0.000000000000E+00
//  v2sigmaaabb   =  0.000000000000E+00
//  v2sigmaab2    =  0.000000000000E+00
//  v2sigmaabbb   =  0.000000000000E+00
//  v2sigmabb2    =  0.000000000000E+00

// test12
//  rhoa= 0.35E+01 rhob= 0.00E+00 sigmaaa= 0.46E-10 sigmaab= 0.00E+00 sigmabb=
//  0.00E+00 zk            = -0.494484233083E+01 vrhoa         =
//  -0.188374945936E+01 vrhob         =  0.000000000000E+00 vsigmaaa      =
//  0.000000000000E+00 vsigmaab      =  0.000000000000E+00 vsigmabb      =
//  0.000000000000E+00

//  v2rhoa2       = -0.179404710416E+00
//  v2rhoab       =  0.000000000000E+00
//  v2rhob2       =  0.000000000000E+00

//  v2rhoasigmaaa =  0.000000000000E+00
//  v2rhoasigmaab =  0.000000000000E+00
//  v2rhoasigmabb =  0.000000000000E+00
//  v2rhobsigmaaa =  0.000000000000E+00
//  v2rhobsigmaab =  0.000000000000E+00
//  v2rhobsigmabb =  0.000000000000E+00

//  v2sigmaaa2    =  0.000000000000E+00
//  v2sigmaaaab   =  0.000000000000E+00
//  v2sigmaaabb   =  0.000000000000E+00
//  v2sigmaab2    =  0.000000000000E+00
//  v2sigmaabbb   =  0.000000000000E+00
//  v2sigmabb2    =  0.000000000000E+00
TEST(DfFunctional_Slater, pointwise12) {
    // input
    const double dRhoA = 0.35E+01;
    const double dRhoB = 0.00E+00;

    // expected value
    const double zk = -0.494484233083E+01;
    const double vRhoA = -0.188374945936E+01;
    const double vRhoB = 0.000000000000E+00;

    // execute test
    DfFunctional_Slater f;

    double dFunctionalValue = f.getFunctional(dRhoA, dRhoB);
    EXPECT_NEAR(zk, dFunctionalValue, EPS);

    double dRoundF_roundRhoA, dRoundF_roundRhoB;
    f.getDerivativeFunctional(dRhoA, dRhoB, &dRoundF_roundRhoA,
                              &dRoundF_roundRhoB);
    EXPECT_NEAR(vRhoA, dRoundF_roundRhoA, EPS);
    EXPECT_NEAR(vRhoB, dRoundF_roundRhoB, EPS);
}

// test13
//  rhoa= 0.35E+01 rhob= 0.00E+00 sigmaaa= 0.34E+01 sigmaab= 0.00E+00 sigmabb=
//  0.00E+00 zk            = -0.494484233083E+01 vrhoa         =
//  -0.188374945936E+01 vrhob         =  0.000000000000E+00 vsigmaaa      =
//  0.000000000000E+00 vsigmaab      =  0.000000000000E+00 vsigmabb      =
//  0.000000000000E+00

//  v2rhoa2       = -0.179404710416E+00
//  v2rhoab       =  0.000000000000E+00
//  v2rhob2       =  0.000000000000E+00

//  v2rhoasigmaaa =  0.000000000000E+00
//  v2rhoasigmaab =  0.000000000000E+00
//  v2rhoasigmabb =  0.000000000000E+00
//  v2rhobsigmaaa =  0.000000000000E+00
//  v2rhobsigmaab =  0.000000000000E+00
//  v2rhobsigmabb =  0.000000000000E+00

//  v2sigmaaa2    =  0.000000000000E+00
//  v2sigmaaaab   =  0.000000000000E+00
//  v2sigmaaabb   =  0.000000000000E+00
//  v2sigmaab2    =  0.000000000000E+00
//  v2sigmaabbb   =  0.000000000000E+00
//  v2sigmabb2    =  0.000000000000E+00
TEST(DfFunctional_Slater, pointwise13) {
    // input
    const double dRhoA = 0.35E+01;
    const double dRhoB = 0.00E+00;

    // expected value
    const double zk = -0.494484233083E+01;
    const double vRhoA = -0.188374945936E+01;
    const double vRhoB = 0.000000000000E+00;

    // execute test
    DfFunctional_Slater f;

    double dFunctionalValue = f.getFunctional(dRhoA, dRhoB);
    EXPECT_NEAR(zk, dFunctionalValue, EPS);

    double dRoundF_roundRhoA, dRoundF_roundRhoB;
    f.getDerivativeFunctional(dRhoA, dRhoB, &dRoundF_roundRhoA,
                              &dRoundF_roundRhoB);
    EXPECT_NEAR(vRhoA, dRoundF_roundRhoA, EPS);
    EXPECT_NEAR(vRhoB, dRoundF_roundRhoB, EPS);
}

// test14
//  rhoa= 0.30E+01 rhob= 0.00E+00 sigmaaa= 0.20E+03 sigmaab= 0.00E+00 sigmabb=
//  0.00E+00 zk            = -0.402615103023E+01 vrhoa         =
//  -0.178940045788E+01 vrhob         =  0.000000000000E+00 vsigmaaa      =
//  0.000000000000E+00 vsigmaab      =  0.000000000000E+00 vsigmabb      =
//  0.000000000000E+00

//  v2rhoa2       = -0.198822273098E+00
//  v2rhoab       =  0.000000000000E+00
//  v2rhob2       =  0.000000000000E+00

//  v2rhoasigmaaa =  0.000000000000E+00
//  v2rhoasigmaab =  0.000000000000E+00
//  v2rhoasigmabb =  0.000000000000E+00
//  v2rhobsigmaaa =  0.000000000000E+00
//  v2rhobsigmaab =  0.000000000000E+00
//  v2rhobsigmabb =  0.000000000000E+00

//  v2sigmaaa2    =  0.000000000000E+00
//  v2sigmaaaab   =  0.000000000000E+00
//  v2sigmaaabb   =  0.000000000000E+00
//  v2sigmaab2    =  0.000000000000E+00
//  v2sigmaabbb   =  0.000000000000E+00
//  v2sigmabb2    =  0.000000000000E+00
TEST(DfFunctional_Slater, pointwise14) {
    // input
    const double dRhoA = 0.30E+01;
    const double dRhoB = 0.00E+00;

    // expected value
    const double zk = -0.402615103023E+01;
    const double vRhoA = -0.178940045788E+01;
    const double vRhoB = 0.000000000000E+00;

    // execute test
    DfFunctional_Slater f;

    double dFunctionalValue = f.getFunctional(dRhoA, dRhoB);
    EXPECT_NEAR(zk, dFunctionalValue, EPS);

    double dRoundF_roundRhoA, dRoundF_roundRhoB;
    f.getDerivativeFunctional(dRhoA, dRhoB, &dRoundF_roundRhoA,
                              &dRoundF_roundRhoB);
    EXPECT_NEAR(vRhoA, dRoundF_roundRhoA, EPS);
    EXPECT_NEAR(vRhoB, dRoundF_roundRhoB, EPS);
}

// test15
//  rhoa= 0.58E-01 rhob= 0.00E+00 sigmaaa= 0.47E-01 sigmaab= 0.00E+00 sigmabb=
//  0.00E+00 zk            = -0.208913119508E-01 vrhoa         =
//  -0.480260044845E+00 vrhob         =  0.000000000000E+00 vsigmaaa      =
//  0.000000000000E+00 vsigmaab      =  0.000000000000E+00 vsigmabb      =
//  0.000000000000E+00

//  v2rhoa2       = -0.276011520026E+01
//  v2rhoab       =  0.000000000000E+00
//  v2rhob2       =  0.000000000000E+00

//  v2rhoasigmaaa =  0.000000000000E+00
//  v2rhoasigmaab =  0.000000000000E+00
//  v2rhoasigmabb =  0.000000000000E+00
//  v2rhobsigmaaa =  0.000000000000E+00
//  v2rhobsigmaab =  0.000000000000E+00
//  v2rhobsigmabb =  0.000000000000E+00

//  v2sigmaaa2    =  0.000000000000E+00
//  v2sigmaaaab   =  0.000000000000E+00
//  v2sigmaaabb   =  0.000000000000E+00
//  v2sigmaab2    =  0.000000000000E+00
//  v2sigmaabbb   =  0.000000000000E+00
//  v2sigmabb2    =  0.000000000000E+00
TEST(DfFunctional_Slater, pointwise15) {
    // input
    const double dRhoA = 0.58E-01;
    const double dRhoB = 0.00E+00;

    // expected value
    const double zk = -0.208913119508E-01;
    const double vRhoA = -0.480260044845E+00;
    const double vRhoB = 0.000000000000E+00;

    // execute test
    DfFunctional_Slater f;

    double dFunctionalValue = f.getFunctional(dRhoA, dRhoB);
    EXPECT_NEAR(zk, dFunctionalValue, EPS);

    double dRoundF_roundRhoA, dRoundF_roundRhoB;
    f.getDerivativeFunctional(dRhoA, dRhoB, &dRoundF_roundRhoA,
                              &dRoundF_roundRhoB);
    EXPECT_NEAR(vRhoA, dRoundF_roundRhoA, EPS);
    EXPECT_NEAR(vRhoB, dRoundF_roundRhoB, EPS);
}

// test16
//  rhoa= 0.82E+02 rhob= 0.81E+02 sigmaaa= 0.49E+07 sigmaab= 0.49E+07 sigmabb=
//  0.49E+07 zk            = -0.657615683804E+03 vrhoa         =
//  -0.539020244480E+01 vrhob         = -0.536820137364E+01 vsigmaaa      =
//  0.000000000000E+00 vsigmaab      =  0.000000000000E+00 vsigmabb      =
//  0.000000000000E+00

//  v2rhoa2       = -0.219113920520E-01
//  v2rhoab       =  0.000000000000E+00
//  v2rhob2       = -0.220913636775E-01

//  v2rhoasigmaaa =  0.000000000000E+00
//  v2rhoasigmaab =  0.000000000000E+00
//  v2rhoasigmabb =  0.000000000000E+00
//  v2rhobsigmaaa =  0.000000000000E+00
//  v2rhobsigmaab =  0.000000000000E+00
//  v2rhobsigmabb =  0.000000000000E+00

//  v2sigmaaa2    =  0.000000000000E+00
//  v2sigmaaaab   =  0.000000000000E+00
//  v2sigmaaabb   =  0.000000000000E+00
//  v2sigmaab2    =  0.000000000000E+00
//  v2sigmaabbb   =  0.000000000000E+00
//  v2sigmabb2    =  0.000000000000E+00
TEST(DfFunctional_Slater, pointwise16) {
    // input
    const double dRhoA = 0.82E+02;
    const double dRhoB = 0.81E+02;

    // expected value
    const double zk = -0.657615683804E+03;
    const double vRhoA = -0.539020244480E+01;
    const double vRhoB = -0.536820137364E+01;

    // execute test
    DfFunctional_Slater f;

    double dFunctionalValue = f.getFunctional(dRhoA, dRhoB);
    EXPECT_NEAR(zk, dFunctionalValue, EPS);

    double dRoundF_roundRhoA, dRoundF_roundRhoB;
    f.getDerivativeFunctional(dRhoA, dRhoB, &dRoundF_roundRhoA,
                              &dRoundF_roundRhoB);
    EXPECT_NEAR(vRhoA, dRoundF_roundRhoA, EPS);
    EXPECT_NEAR(vRhoB, dRoundF_roundRhoB, EPS);
}

// test17
//  rhoa= 0.39E+02 rhob= 0.38E+02 sigmaaa= 0.81E+06 sigmaab= 0.82E+06 sigmabb=
//  0.82E+06 zk            = -0.241948147838E+03 vrhoa         =
//  -0.420747936684E+01 vrhob         = -0.417120618800E+01 vsigmaaa      =
//  0.000000000000E+00 vsigmaab      =  0.000000000000E+00 vsigmabb      =
//  0.000000000000E+00

//  v2rhoa2       = -0.359613621097E-01
//  v2rhoab       =  0.000000000000E+00
//  v2rhob2       = -0.365895279649E-01

//  v2rhoasigmaaa =  0.000000000000E+00
//  v2rhoasigmaab =  0.000000000000E+00
//  v2rhoasigmabb =  0.000000000000E+00
//  v2rhobsigmaaa =  0.000000000000E+00
//  v2rhobsigmaab =  0.000000000000E+00
//  v2rhobsigmabb =  0.000000000000E+00

//  v2sigmaaa2    =  0.000000000000E+00
//  v2sigmaaaab   =  0.000000000000E+00
//  v2sigmaaabb   =  0.000000000000E+00
//  v2sigmaab2    =  0.000000000000E+00
//  v2sigmaabbb   =  0.000000000000E+00
//  v2sigmabb2    =  0.000000000000E+00
TEST(DfFunctional_Slater, pointwise17) {
    // input
    const double dRhoA = 0.39E+02;
    const double dRhoB = 0.38E+02;

    // expected value
    const double zk = -0.241948147838E+03;
    const double vRhoA = -0.420747936684E+01;
    const double vRhoB = -0.417120618800E+01;

    // execute test
    DfFunctional_Slater f;

    double dFunctionalValue = f.getFunctional(dRhoA, dRhoB);
    EXPECT_NEAR(zk, dFunctionalValue, EPS);

    double dRoundF_roundRhoA, dRoundF_roundRhoB;
    f.getDerivativeFunctional(dRhoA, dRhoB, &dRoundF_roundRhoA,
                              &dRoundF_roundRhoB);
    EXPECT_NEAR(vRhoA, dRoundF_roundRhoA, EPS);
    EXPECT_NEAR(vRhoB, dRoundF_roundRhoB, EPS);
}

// test18
//  rhoa= 0.13E+00 rhob= 0.95E-01 sigmaaa= 0.15E+00 sigmaab= 0.18E+00 sigmabb=
//  0.22E+00 zk            = -0.101616142698E+00 vrhoa         =
//  -0.628513933519E+00 vrhob         = -0.566119777958E+00 vsigmaaa      =
//  0.000000000000E+00 vsigmaab      =  0.000000000000E+00 vsigmabb      =
//  0.000000000000E+00

//  v2rhoa2       = -0.161157418851E+01
//  v2rhoab       =  0.000000000000E+00
//  v2rhob2       = -0.198638518582E+01

//  v2rhoasigmaaa =  0.000000000000E+00
//  v2rhoasigmaab =  0.000000000000E+00
//  v2rhoasigmabb =  0.000000000000E+00
//  v2rhobsigmaaa =  0.000000000000E+00
//  v2rhobsigmaab =  0.000000000000E+00
//  v2rhobsigmabb =  0.000000000000E+00

//  v2sigmaaa2    =  0.000000000000E+00
//  v2sigmaaaab   =  0.000000000000E+00
//  v2sigmaaabb   =  0.000000000000E+00
//  v2sigmaab2    =  0.000000000000E+00
//  v2sigmaabbb   =  0.000000000000E+00
//  v2sigmabb2    =  0.000000000000E+00
TEST(DfFunctional_Slater, pointwise18) {
    // input
    const double dRhoA = 0.13E+00;
    const double dRhoB = 0.95E-01;

    // expected value
    const double zk = -0.101616142698E+00;
    const double vRhoA = -0.628513933519E+00;
    const double vRhoB = -0.566119777958E+00;

    // execute test
    DfFunctional_Slater f;

    double dFunctionalValue = f.getFunctional(dRhoA, dRhoB);
    EXPECT_NEAR(zk, dFunctionalValue, EPS);

    double dRoundF_roundRhoA, dRoundF_roundRhoB;
    f.getDerivativeFunctional(dRhoA, dRhoB, &dRoundF_roundRhoA,
                              &dRoundF_roundRhoB);
    EXPECT_NEAR(vRhoA, dRoundF_roundRhoA, EPS);
    EXPECT_NEAR(vRhoB, dRoundF_roundRhoB, EPS);
}

// test19
//  rhoa= 0.78E-01 rhob= 0.31E-01 sigmaaa= 0.41E-02 sigmaab= 0.38E-02 sigmabb=
//  0.36E-02 zk            = -0.400731073431E-01 vrhoa         =
//  -0.530109182127E+00 vrhob         = -0.389751405963E+00 vsigmaaa      =
//  0.000000000000E+00 vsigmaab      =  0.000000000000E+00 vsigmabb      =
//  0.000000000000E+00

//  v2rhoa2       = -0.226542385525E+01
//  v2rhoab       =  0.000000000000E+00
//  v2rhob2       = -0.419087533293E+01

//  v2rhoasigmaaa =  0.000000000000E+00
//  v2rhoasigmaab =  0.000000000000E+00
//  v2rhoasigmabb =  0.000000000000E+00
//  v2rhobsigmaaa =  0.000000000000E+00
//  v2rhobsigmaab =  0.000000000000E+00
//  v2rhobsigmabb =  0.000000000000E+00

//  v2sigmaaa2    =  0.000000000000E+00
//  v2sigmaaaab   =  0.000000000000E+00
//  v2sigmaaabb   =  0.000000000000E+00
//  v2sigmaab2    =  0.000000000000E+00
//  v2sigmaabbb   =  0.000000000000E+00
//  v2sigmabb2    =  0.000000000000E+00
TEST(DfFunctional_Slater, pointwise19) {
    // input
    const double dRhoA = 0.78E-01;
    const double dRhoB = 0.31E-01;

    // expected value
    const double zk = -0.400731073431E-01;
    const double vRhoA = -0.530109182127E+00;
    const double vRhoB = -0.389751405963E+00;

    // execute test
    DfFunctional_Slater f;

    double dFunctionalValue = f.getFunctional(dRhoA, dRhoB);
    EXPECT_NEAR(zk, dFunctionalValue, EPS);

    double dRoundF_roundRhoA, dRoundF_roundRhoB;
    f.getDerivativeFunctional(dRhoA, dRhoB, &dRoundF_roundRhoA,
                              &dRoundF_roundRhoB);
    EXPECT_NEAR(vRhoA, dRoundF_roundRhoA, EPS);
    EXPECT_NEAR(vRhoB, dRoundF_roundRhoB, EPS);
}

// test20
//  rhoa= 0.50E+02 rhob= 0.49E+02 sigmaaa= 0.11E+06 sigmaab= 0.11E+06 sigmabb=
//  0.11E+06 zk            = -0.338253135027E+03 vrhoa         =
//  -0.457078149734E+01 vrhob         = -0.454010418713E+01 vsigmaaa      =
//  0.000000000000E+00 vsigmaab      =  0.000000000000E+00 vsigmabb      =
//  0.000000000000E+00

//  v2rhoa2       = -0.304718766489E-01
//  v2rhoab       =  0.000000000000E+00
//  v2rhob2       = -0.308850624975E-01

//  v2rhoasigmaaa =  0.000000000000E+00
//  v2rhoasigmaab =  0.000000000000E+00
//  v2rhoasigmabb =  0.000000000000E+00
//  v2rhobsigmaaa =  0.000000000000E+00
//  v2rhobsigmaab =  0.000000000000E+00
//  v2rhobsigmabb =  0.000000000000E+00

//  v2sigmaaa2    =  0.000000000000E+00
//  v2sigmaaaab   =  0.000000000000E+00
//  v2sigmaaabb   =  0.000000000000E+00
//  v2sigmaab2    =  0.000000000000E+00
//  v2sigmaabbb   =  0.000000000000E+00
//  v2sigmabb2    =  0.000000000000E+00
TEST(DfFunctional_Slater, pointwise20) {
    // input
    const double dRhoA = 0.50E+02;
    const double dRhoB = 0.49E+02;

    // expected value
    const double zk = -0.338253135027E+03;
    const double vRhoA = -0.457078149734E+01;
    const double vRhoB = -0.454010418713E+01;

    // execute test
    DfFunctional_Slater f;

    double dFunctionalValue = f.getFunctional(dRhoA, dRhoB);
    EXPECT_NEAR(zk, dFunctionalValue, EPS);

    double dRoundF_roundRhoA, dRoundF_roundRhoB;
    f.getDerivativeFunctional(dRhoA, dRhoB, &dRoundF_roundRhoA,
                              &dRoundF_roundRhoB);
    EXPECT_NEAR(vRhoA, dRoundF_roundRhoA, EPS);
    EXPECT_NEAR(vRhoB, dRoundF_roundRhoB, EPS);
}

// test21
//  rhoa= 0.40E+02 rhob= 0.40E+02 sigmaaa= 0.99E+05 sigmaab= 0.98E+05 sigmabb=
//  0.98E+05 zk            = -0.254588260307E+03 vrhoa         =
//  -0.424313767179E+01 vrhob         = -0.424313767179E+01 vsigmaaa      =
//  0.000000000000E+00 vsigmaab      =  0.000000000000E+00 vsigmabb      =
//  0.000000000000E+00

//  v2rhoa2       = -0.353594805982E-01
//  v2rhoab       =  0.000000000000E+00
//  v2rhob2       = -0.353594805982E-01

//  v2rhoasigmaaa =  0.000000000000E+00
//  v2rhoasigmaab =  0.000000000000E+00
//  v2rhoasigmabb =  0.000000000000E+00
//  v2rhobsigmaaa =  0.000000000000E+00
//  v2rhobsigmaab =  0.000000000000E+00
//  v2rhobsigmabb =  0.000000000000E+00

//  v2sigmaaa2    =  0.000000000000E+00
//  v2sigmaaaab   =  0.000000000000E+00
//  v2sigmaaabb   =  0.000000000000E+00
//  v2sigmaab2    =  0.000000000000E+00
//  v2sigmaabbb   =  0.000000000000E+00
//  v2sigmabb2    =  0.000000000000E+00

// test22
//  rhoa= 0.12E+00 rhob= 0.10E+00 sigmaaa= 0.12E+00 sigmaab= 0.13E+00 sigmabb=
//  0.14E+00 zk            = -0.982681500273E-01 vrhoa         =
//  -0.611966348389E+00 vrhob         = -0.575882382297E+00 vsigmaaa      =
//  0.000000000000E+00 vsigmaab      =  0.000000000000E+00 vsigmabb      =
//  0.000000000000E+00

//  v2rhoa2       = -0.169990652330E+01
//  v2rhoab       =  0.000000000000E+00
//  v2rhob2       = -0.191960794099E+01

//  v2rhoasigmaaa =  0.000000000000E+00
//  v2rhoasigmaab =  0.000000000000E+00
//  v2rhoasigmabb =  0.000000000000E+00
//  v2rhobsigmaaa =  0.000000000000E+00
//  v2rhobsigmaab =  0.000000000000E+00
//  v2rhobsigmabb =  0.000000000000E+00

//  v2sigmaaa2    =  0.000000000000E+00
//  v2sigmaaaab   =  0.000000000000E+00
//  v2sigmaaabb   =  0.000000000000E+00
//  v2sigmaab2    =  0.000000000000E+00
//  v2sigmaabbb   =  0.000000000000E+00
//  v2sigmabb2    =  0.000000000000E+00

// test23
//  rhoa= 0.48E-01 rhob= 0.25E-01 sigmaaa= 0.46E-02 sigmaab= 0.44E-02 sigmabb=
//  0.41E-02 zk            = -0.230346081831E-01 vrhoa         =
//  -0.450900660715E+00 vrhob         = -0.362783167860E+00 vsigmaaa      =
//  0.000000000000E+00 vsigmaab      =  0.000000000000E+00 vsigmabb      =
//  0.000000000000E+00

//  v2rhoa2       = -0.313125458830E+01
//  v2rhoab       =  0.000000000000E+00
//  v2rhob2       = -0.483710890480E+01

//  v2rhoasigmaaa =  0.000000000000E+00
//  v2rhoasigmaab =  0.000000000000E+00
//  v2rhoasigmabb =  0.000000000000E+00
//  v2rhobsigmaaa =  0.000000000000E+00
//  v2rhobsigmaab =  0.000000000000E+00
//  v2rhobsigmabb =  0.000000000000E+00

//  v2sigmaaa2    =  0.000000000000E+00
//  v2sigmaaaab   =  0.000000000000E+00
//  v2sigmaaabb   =  0.000000000000E+00
//  v2sigmaab2    =  0.000000000000E+00
//  v2sigmaabbb   =  0.000000000000E+00
//  v2sigmabb2    =  0.000000000000E+00
