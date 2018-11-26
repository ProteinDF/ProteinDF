#include "DfFunctional_Becke88.h"
#include <string>
#include <vector>
#include "gtest/gtest.h"

const double EPS = 1.0E-5;

// =====================================================================
// data used from: http://www.cse.scitech.ac.uk/ccg/dft/data_pt_x_lda.html
//

// test1
//  rhoa = 0.17E+01
//  rhob = 0.17E+01
//  sigmaaa = 0.81E-11
//  sigmaab = 0.81E-11
//  sigmabb =  0.81E-11
//  zk            = -0.377592720836E+01
//  vrhoa         =  -0.148075576798E+01
//  vrhob         = -0.148075576798E+01
//  vsigmaaa      =  -0.207006537839E-02
//  vsigmaab      =  0.000000000000E+00
//  vsigmabb      =  -0.207006537839E-02

//  v2rhoa2       = -0.290344268232E+00
//  v2rhoab       =  0.000000000000E+00
//  v2rhob2       = -0.290344268232E+00

//  v2rhoasigmaaa =  0.162358068893E-02
//  v2rhoasigmaab =  0.000000000000E+00
//  v2rhoasigmabb =  0.000000000000E+00
//  v2rhobsigmaaa =  0.000000000000E+00
//  v2rhobsigmaab =  0.000000000000E+00
//  v2rhobsigmabb =  0.162358068893E-02

//  v2sigmaaa2    =  0.253445241319E-04
//  v2sigmaaaab   =  0.000000000000E+00
//  v2sigmaaabb   =  0.000000000000E+00
//  v2sigmaab2    =  0.000000000000E+00
//  v2sigmaabbb   =  0.000000000000E+00
//  v2sigmabb2    =  0.253445241319E-04
TEST(DfFunctional_Becke88, pointwise1) {
  // input
  const double dRhoA = 0.17E+01;
  const double dRhoB = 0.17E+01;
  const double dGammaAA = 0.81E-11;
  const double dGammaAB = 0.81E-11;
  const double dGammaBB = 0.81E-11;

  // expected value
  const double zk = -0.377592720836E+01;
  const double vRhoA = -0.148075576798E+01;
  const double vRhoB = -0.148075576798E+01;
  const double vGammaAA = -0.207006537839E-02;
  const double vGammaAB = 0.000000000000E+00;
  const double vGammaBB = -0.207006537839E-02;

  // execute test
  DfFunctional_Becke88 f;

  double dFunctionalValue =
      f.getFunctional(dRhoA, dRhoB, dGammaAA, dGammaAB, dGammaBB);
  EXPECT_NEAR(zk, dFunctionalValue, EPS);

  double dRoundF_roundRhoA, dRoundF_roundRhoB;
  double dRoundF_roundGammaAA, dRoundF_roundGammaAB, dRoundF_roundGammaBB;
  f.getDerivativeFunctional(dRhoA, dRhoB, dGammaAA, dGammaAB, dGammaBB,
                            &dRoundF_roundRhoA, &dRoundF_roundRhoB,
                            &dRoundF_roundGammaAA, &dRoundF_roundGammaAB,
                            &dRoundF_roundGammaBB);

  EXPECT_NEAR(vRhoA, dRoundF_roundRhoA, EPS);
  EXPECT_NEAR(vRhoB, dRoundF_roundRhoB, EPS);
  EXPECT_NEAR(vGammaAA, dRoundF_roundGammaAA, EPS);
  EXPECT_NEAR(vGammaAB, dRoundF_roundGammaAB, EPS);
  EXPECT_NEAR(vGammaBB, dRoundF_roundGammaBB, EPS);
}

TEST(DfFunctional_Becke88, pointwise1_RKS) {
  // input
  const double dRhoA = 0.17E+01;
  //   const double dRhoB = 0.17E+01;
  const double dGammaAA = 0.81E-11;
  //   const double dGammaAB = 0.81E-11;
  //   const double dGammaBB = 0.81E-11;

  // expected value
  const double zk = -0.377592720836E+01;
  const double vRhoA = -0.148075576798E+01;
  //   const double vRhoB = -0.148075576798E+01;
  const double vGammaAA = -0.207006537839E-02;
  const double vGammaAB = 0.000000000000E+00;
  //   const double vGammaBB = -0.207006537839E-02;

  // execute test
  DfFunctional_Becke88 f;

  double dFunctionalValue = f.getFunctional(dRhoA, dGammaAA);
  EXPECT_NEAR(zk, dFunctionalValue, EPS);

  double dRoundF_roundRhoA;
  double dRoundF_roundGammaAA, dRoundF_roundGammaAB;
  f.getDerivativeFunctional(dRhoA, dGammaAA, &dRoundF_roundRhoA,
                            &dRoundF_roundGammaAA, &dRoundF_roundGammaAB);

  EXPECT_NEAR(vRhoA, dRoundF_roundRhoA, EPS);
  //   EXPECT_NEAR(vRhoB, dRoundF_roundRhoB, EPS);
  EXPECT_NEAR(vGammaAA, dRoundF_roundGammaAA, EPS);
  EXPECT_NEAR(vGammaAB, dRoundF_roundGammaAB, EPS);
  //   EXPECT_NEAR(vGammaBB, dRoundF_roundGammaBB, EPS);
}

// test2
//  rhoa= 0.17E+01
//  rhob= 0.17E+01
//  sigmaaa= 0.17E+01
//  sigmaab= 0.17E+01
//  sigmabb=  0.17E+01
//  zk            = -0.378289713911E+01
//  vrhoa         = -0.147807268065E+01
//  vrhob         = -0.147807268065E+01
//  vsigmaaa      = -0.203114756676E-02
//  vsigmaab      =  0.000000000000E+00
//  vsigmabb      = -0.203114756676E-02

//  v2rhoa2       = -0.293917869186E+00
//  v2rhoab       =  0.000000000000E+00
//  v2rhob2       = -0.293917869186E+00

//  v2rhoasigmaaa =  0.153738619102E-02
//  v2rhoasigmaab =  0.000000000000E+00
//  v2rhoasigmabb =  0.000000000000E+00
//  v2rhobsigmaaa =  0.000000000000E+00
//  v2rhobsigmaab =  0.000000000000E+00
//  v2rhobsigmabb =  0.153738619102E-02

//  v2sigmaaa2    =  0.208765215311E-04
//  v2sigmaaaab   =  0.000000000000E+00
//  v2sigmaaabb   =  0.000000000000E+00
//  v2sigmaab2    =  0.000000000000E+00
//  v2sigmaabbb   =  0.000000000000E+00
//  v2sigmabb2    =  0.208765215311E-04
TEST(DfFunctional_Becke88, pointwise2) {
  // input
  const double dRhoA = 0.17E+01;
  const double dRhoB = 0.17E+01;
  const double dGammaAA = 0.17E+01;
  const double dGammaAB = 0.17E+01;
  const double dGammaBB = 0.17E+01;

  // expected value
  const double zk = -0.378289713911E+01;
  const double vRhoA = -0.147807268065E+01;
  const double vRhoB = -0.147807268065E+01;
  const double vGammaAA = -0.203114756676E-02;
  const double vGammaAB = 0.000000000000E+00;
  const double vGammaBB = -0.203114756676E-02;

  // execute test
  DfFunctional_Becke88 f;

  double dFunctionalValue =
      f.getFunctional(dRhoA, dRhoB, dGammaAA, dGammaAB, dGammaBB);
  EXPECT_NEAR(zk, dFunctionalValue, EPS);

  double dRoundF_roundRhoA, dRoundF_roundRhoB;
  double dRoundF_roundGammaAA, dRoundF_roundGammaAB, dRoundF_roundGammaBB;
  f.getDerivativeFunctional(dRhoA, dRhoB, dGammaAA, dGammaAB, dGammaBB,
                            &dRoundF_roundRhoA, &dRoundF_roundRhoB,
                            &dRoundF_roundGammaAA, &dRoundF_roundGammaAB,
                            &dRoundF_roundGammaBB);

  EXPECT_NEAR(vRhoA, dRoundF_roundRhoA, EPS);
  EXPECT_NEAR(vRhoB, dRoundF_roundRhoB, EPS);
  EXPECT_NEAR(vGammaAA, dRoundF_roundGammaAA, EPS);
  EXPECT_NEAR(vGammaAB, dRoundF_roundGammaAB, EPS);
  EXPECT_NEAR(vGammaBB, dRoundF_roundGammaBB, EPS);
}

TEST(DfFunctional_Becke88, pointwise2_RKS) {
  // input
  const double dRhoA = 0.17E+01;
  //   const double dRhoB = 0.17E+01;
  const double dGammaAA = 0.17E+01;
  //   const double dGammaAB = 0.17E+01;
  //   const double dGammaBB = 0.17E+01;

  // expected value
  const double zk = -0.378289713911E+01;
  const double vRhoA = -0.147807268065E+01;
  //   const double vRhoB = -0.147807268065E+01;
  const double vGammaAA = -0.203114756676E-02;
  const double vGammaAB = 0.000000000000E+00;
  //   const double vGammaBB = -0.203114756676E-02;

  // execute test
  DfFunctional_Becke88 f;

  double dFunctionalValue = f.getFunctional(dRhoA, dGammaAA);
  EXPECT_NEAR(zk, dFunctionalValue, EPS);

  double dRoundF_roundRhoA;
  double dRoundF_roundGammaAA, dRoundF_roundGammaAB;
  f.getDerivativeFunctional(dRhoA, dGammaAA, &dRoundF_roundRhoA,
                            &dRoundF_roundGammaAA, &dRoundF_roundGammaAB);

  EXPECT_NEAR(vRhoA, dRoundF_roundRhoA, EPS);
  //   EXPECT_NEAR(vRhoB, dRoundF_roundRhoB, EPS);
  EXPECT_NEAR(vGammaAA, dRoundF_roundGammaAA, EPS);
  EXPECT_NEAR(vGammaAB, dRoundF_roundGammaAB, EPS);
  //   EXPECT_NEAR(vGammaBB, dRoundF_roundGammaBB, EPS);
}

// test3
//  rhoa= 0.88E-01 rhob= 0.88E-01 sigmaaa= 0.87E-01 sigmaab= 0.87E-01 sigmabb=
//  0.87E-01 zk            = -0.851611545044E-01 vrhoa         =
//  -0.501899165865E+00 vrhob         = -0.501899165865E+00 vsigmaaa      =
//  -0.543404155466E-01 vsigmaab      =  0.000000000000E+00 vsigmabb      =
//  -0.543404155466E-01

//  v2rhoa2       = -0.253321231850E+01
//  v2rhoab       =  0.000000000000E+00
//  v2rhob2       = -0.253321231850E+01

//  v2rhoasigmaaa =  0.239754146866E+00
//  v2rhoasigmaab =  0.000000000000E+00
//  v2rhoasigmabb =  0.000000000000E+00
//  v2rhobsigmaaa =  0.000000000000E+00
//  v2rhobsigmaab =  0.000000000000E+00
//  v2rhobsigmabb =  0.239754146866E+00

//  v2sigmaaa2    =  0.221360010652E+00
//  v2sigmaaaab   =  0.000000000000E+00
//  v2sigmaaabb   =  0.000000000000E+00
//  v2sigmaab2    =  0.000000000000E+00
//  v2sigmaabbb   =  0.000000000000E+00
//  v2sigmabb2    =  0.221360010652E+00
TEST(DfFunctional_Becke88, pointwise3) {
  // input
  const double dRhoA = 0.88E-01;
  const double dRhoB = 0.88E-01;
  const double dGammaAA = 0.87E-01;
  const double dGammaAB = 0.87E-01;
  const double dGammaBB = 0.87E-01;

  // expected value
  const double zk = -0.851611545044E-01;
  const double vRhoA = -0.501899165865E+00;
  const double vRhoB = -0.501899165865E+00;
  const double vGammaAA = -0.543404155466E-01;
  const double vGammaAB = 0.000000000000E+00;
  const double vGammaBB = -0.543404155466E-01;

  // execute test
  DfFunctional_Becke88 f;

  double dFunctionalValue =
      f.getFunctional(dRhoA, dRhoB, dGammaAA, dGammaAB, dGammaBB);
  EXPECT_NEAR(zk, dFunctionalValue, EPS);

  double dRoundF_roundRhoA, dRoundF_roundRhoB;
  double dRoundF_roundGammaAA, dRoundF_roundGammaAB, dRoundF_roundGammaBB;
  f.getDerivativeFunctional(dRhoA, dRhoB, dGammaAA, dGammaAB, dGammaBB,
                            &dRoundF_roundRhoA, &dRoundF_roundRhoB,
                            &dRoundF_roundGammaAA, &dRoundF_roundGammaAB,
                            &dRoundF_roundGammaBB);

  EXPECT_NEAR(vRhoA, dRoundF_roundRhoA, EPS);
  EXPECT_NEAR(vRhoB, dRoundF_roundRhoB, EPS);
  EXPECT_NEAR(vGammaAA, dRoundF_roundGammaAA, EPS);
  EXPECT_NEAR(vGammaAB, dRoundF_roundGammaAB, EPS);
  EXPECT_NEAR(vGammaBB, dRoundF_roundGammaBB, EPS);
}

TEST(DfFunctional_Becke88, pointwise3_RKS) {
  // input
  const double dRhoA = 0.88E-01;
  //   const double dRhoB = 0.88E-01;
  const double dGammaAA = 0.87E-01;
  //   const double dGammaAB = 0.87E-01;
  //   const double dGammaBB = 0.87E-01;

  // expected value
  const double zk = -0.851611545044E-01;
  const double vRhoA = -0.501899165865E+00;
  //   const double vRhoB = -0.501899165865E+00;
  const double vGammaAA = -0.543404155466E-01;
  const double vGammaAB = 0.000000000000E+00;
  //   const double vGammaBB = -0.543404155466E-01;

  // execute test
  DfFunctional_Becke88 f;

  double dFunctionalValue = f.getFunctional(dRhoA, dGammaAA);
  EXPECT_NEAR(zk, dFunctionalValue, EPS);

  double dRoundF_roundRhoA;
  double dRoundF_roundGammaAA, dRoundF_roundGammaAB;
  f.getDerivativeFunctional(dRhoA, dGammaAA, &dRoundF_roundRhoA,
                            &dRoundF_roundGammaAA, &dRoundF_roundGammaAB);

  EXPECT_NEAR(vRhoA, dRoundF_roundRhoA, EPS);
  //   EXPECT_NEAR(vRhoB, dRoundF_roundRhoB, EPS);
  EXPECT_NEAR(vGammaAA, dRoundF_roundGammaAA, EPS);
  EXPECT_NEAR(vGammaAB, dRoundF_roundGammaAB, EPS);
  //   EXPECT_NEAR(vGammaBB, dRoundF_roundGammaBB, EPS);
}

// test4
//  rhoa= 0.18E+04 rhob= 0.18E+04 sigmaaa= 0.55E+00 sigmaab= 0.55E+00 sigmabb=
//  0.55E+00 zk            = -0.407494475322E+05 vrhoa         =
//  -0.150923879747E+02 vrhob         = -0.150923879747E+02 vsigmaaa      =
//  -0.191816494659E-06 vsigmaab      =  0.000000000000E+00 vsigmabb      =
//  -0.191816494659E-06

//  v2rhoa2       = -0.279488666210E-02
//  v2rhoab       =  0.000000000000E+00
//  v2rhob2       = -0.279488666210E-02

//  v2rhoasigmaaa =  0.142086292324E-09
//  v2rhoasigmaab =  0.000000000000E+00
//  v2rhoasigmabb =  0.000000000000E+00
//  v2rhobsigmaaa =  0.000000000000E+00
//  v2rhobsigmaab =  0.000000000000E+00
//  v2rhobsigmabb =  0.142086292324E-09

//  v2sigmaaa2    =  0.201646090402E-16
//  v2sigmaaaab   =  0.000000000000E+00
//  v2sigmaaabb   =  0.000000000000E+00
//  v2sigmaab2    =  0.000000000000E+00
//  v2sigmaabbb   =  0.000000000000E+00
//  v2sigmabb2    =  0.201646090402E-16
TEST(DfFunctional_Becke88, pointwise4) {
  // input
  const double dRhoA = 0.18E+04;
  const double dRhoB = 0.18E+04;
  const double dGammaAA = 0.55E+00;
  const double dGammaAB = 0.55E+00;
  const double dGammaBB = 0.55E+00;

  // expected value
  const double zk = -0.407494475322E+05;
  const double vRhoA = -0.150923879747E+02;
  const double vRhoB = -0.150923879747E+02;
  const double vGammaAA = -0.191816494659E-06;
  const double vGammaAB = 0.000000000000E+00;
  const double vGammaBB = -0.191816494659E-06;

  // execute test
  DfFunctional_Becke88 f;

  double dFunctionalValue =
      f.getFunctional(dRhoA, dRhoB, dGammaAA, dGammaAB, dGammaBB);
  EXPECT_NEAR(zk, dFunctionalValue, EPS);

  double dRoundF_roundRhoA, dRoundF_roundRhoB;
  double dRoundF_roundGammaAA, dRoundF_roundGammaAB, dRoundF_roundGammaBB;
  f.getDerivativeFunctional(dRhoA, dRhoB, dGammaAA, dGammaAB, dGammaBB,
                            &dRoundF_roundRhoA, &dRoundF_roundRhoB,
                            &dRoundF_roundGammaAA, &dRoundF_roundGammaAB,
                            &dRoundF_roundGammaBB);

  EXPECT_NEAR(vRhoA, dRoundF_roundRhoA, EPS);
  EXPECT_NEAR(vRhoB, dRoundF_roundRhoB, EPS);
  EXPECT_NEAR(vGammaAA, dRoundF_roundGammaAA, EPS);
  EXPECT_NEAR(vGammaAB, dRoundF_roundGammaAB, EPS);
  EXPECT_NEAR(vGammaBB, dRoundF_roundGammaBB, EPS);
}

// test5
//  rhoa= 0.18E+04 rhob= 0.18E+04 sigmaaa= 0.86E+04 sigmaab= 0.86E+04 sigmabb=
//  0.86E+04 zk            = -0.407494508312E+05 vrhoa         =
//  -0.150923867529E+02 vrhob         = -0.150923867529E+02 vsigmaaa      =
//  -0.191816321255E-06 vsigmaab      =  0.000000000000E+00 vsigmabb      =
//  -0.191816321255E-06

//  v2rhoa2       = -0.279488824600E-02
//  v2rhoab       =  0.000000000000E+00
//  v2rhob2       = -0.279488824600E-02

//  v2rhoasigmaaa =  0.142085906983E-09
//  v2rhoasigmaab =  0.000000000000E+00
//  v2rhoasigmabb =  0.000000000000E+00
//  v2rhobsigmaaa =  0.000000000000E+00
//  v2rhobsigmaab =  0.000000000000E+00
//  v2rhobsigmabb =  0.142085906983E-09

//  v2sigmaaa2    =  0.201644008560E-16
//  v2sigmaaaab   =  0.000000000000E+00
//  v2sigmaaabb   =  0.000000000000E+00
//  v2sigmaab2    =  0.000000000000E+00
//  v2sigmaabbb   =  0.000000000000E+00
//  v2sigmabb2    =  0.201644008560E-16
TEST(DfFunctional_Becke88, pointwise5) {
  // input
  const double dRhoA = 0.18E+04;
  const double dRhoB = 0.18E+04;
  const double dGammaAA = 0.86E+04;
  const double dGammaAB = 0.86E+04;
  const double dGammaBB = 0.86E+04;

  // expected value
  const double zk = -0.407494508312E+05;
  const double vRhoA = -0.150923867529E+02;
  const double vRhoB = -0.150923867529E+02;
  const double vGammaAA = -0.191816321255E-06;
  const double vGammaAB = 0.000000000000E+00;
  const double vGammaBB = -0.191816321255E-06;

  // execute test
  DfFunctional_Becke88 f;

  double dFunctionalValue =
      f.getFunctional(dRhoA, dRhoB, dGammaAA, dGammaAB, dGammaBB);
  EXPECT_NEAR(zk, dFunctionalValue, EPS);

  double dRoundF_roundRhoA, dRoundF_roundRhoB;
  double dRoundF_roundGammaAA, dRoundF_roundGammaAB, dRoundF_roundGammaBB;
  f.getDerivativeFunctional(dRhoA, dRhoB, dGammaAA, dGammaAB, dGammaBB,
                            &dRoundF_roundRhoA, &dRoundF_roundRhoB,
                            &dRoundF_roundGammaAA, &dRoundF_roundGammaAB,
                            &dRoundF_roundGammaBB);

  EXPECT_NEAR(vRhoA, dRoundF_roundRhoA, EPS);
  EXPECT_NEAR(vRhoB, dRoundF_roundRhoB, EPS);
  EXPECT_NEAR(vGammaAA, dRoundF_roundGammaAA, EPS);
  EXPECT_NEAR(vGammaAB, dRoundF_roundGammaAB, EPS);
  EXPECT_NEAR(vGammaBB, dRoundF_roundGammaBB, EPS);
}

// test6
//  rhoa= 0.16E+04 rhob= 0.16E+04 sigmaaa= 0.37E+10 sigmaab= 0.37E+10 sigmabb=
//  0.37E+10 zk            = -0.362648637930E+05 vrhoa         =
//  -0.140333722784E+02 vrhob         = -0.140333722784E+02 vsigmaaa      =
//  -0.174646643568E-06 vsigmaab      =  0.000000000000E+00 vsigmabb      =
//  -0.174646643568E-06

//  v2rhoa2       = -0.352244394438E-02
//  v2rhoab       =  0.000000000000E+00
//  v2rhob2       = -0.352244394438E-02

//  v2rhoasigmaaa =  0.971067113036E-10
//  v2rhoasigmaab =  0.000000000000E+00
//  v2rhoasigmabb =  0.000000000000E+00
//  v2rhobsigmaaa =  0.000000000000E+00
//  v2rhobsigmaab =  0.000000000000E+00
//  v2rhobsigmabb =  0.971067113036E-10

//  v2sigmaaa2    =  0.785386351399E-17
//  v2sigmaaaab   =  0.000000000000E+00
//  v2sigmaaabb   =  0.000000000000E+00
//  v2sigmaab2    =  0.000000000000E+00
//  v2sigmaabbb   =  0.000000000000E+00
//  v2sigmabb2    =  0.785386351399E-17
TEST(DfFunctional_Becke88, pointwise6) {
  // input
  const double dRhoA = 0.16E+04;
  const double dRhoB = 0.16E+04;
  const double dGammaAA = 0.37E+10;
  const double dGammaAB = 0.37E+10;
  const double dGammaBB = 0.37E+10;

  // expected value
  const double zk = -0.362648637930E+05;
  const double vRhoA = -0.140333722784E+02;
  const double vRhoB = -0.140333722784E+02;
  const double vGammaAA = -0.174646643568E-06;
  const double vGammaAB = 0.000000000000E+00;
  const double vGammaBB = -0.174646643568E-06;

  // execute test
  DfFunctional_Becke88 f;

  double dFunctionalValue =
      f.getFunctional(dRhoA, dRhoB, dGammaAA, dGammaAB, dGammaBB);
  EXPECT_NEAR(zk, dFunctionalValue, EPS);

  double dRoundF_roundRhoA, dRoundF_roundRhoB;
  double dRoundF_roundGammaAA, dRoundF_roundGammaAB, dRoundF_roundGammaBB;
  f.getDerivativeFunctional(dRhoA, dRhoB, dGammaAA, dGammaAB, dGammaBB,
                            &dRoundF_roundRhoA, &dRoundF_roundRhoB,
                            &dRoundF_roundGammaAA, &dRoundF_roundGammaAB,
                            &dRoundF_roundGammaBB);

  EXPECT_NEAR(vRhoA, dRoundF_roundRhoA, EPS);
  EXPECT_NEAR(vRhoB, dRoundF_roundRhoB, EPS);
  EXPECT_NEAR(vGammaAA, dRoundF_roundGammaAA, EPS);
  EXPECT_NEAR(vGammaAB, dRoundF_roundGammaAB, EPS);
  EXPECT_NEAR(vGammaBB, dRoundF_roundGammaBB, EPS);
}

// test7
//  rhoa= 0.26E+00 rhob= 0.26E+00 sigmaaa= 0.28E+00 sigmaab= 0.28E+00 sigmabb=
//  0.28E+00 zk            = -0.321148637763E+00 vrhoa         =
//  -0.766539815464E+00 vrhob         = -0.766539815464E+00 vsigmaaa      =
//  -0.198197408319E-01 vsigmaab      =  0.000000000000E+00 vsigmabb      =
//  -0.198197408319E-01

//  v2rhoa2       = -0.117949870779E+01
//  v2rhoab       =  0.000000000000E+00
//  v2rhob2       = -0.117949870779E+01

//  v2rhoasigmaaa =  0.685130252745E-01
//  v2rhoasigmaab =  0.000000000000E+00
//  v2rhoasigmabb =  0.000000000000E+00
//  v2rhobsigmaaa =  0.000000000000E+00
//  v2rhobsigmaab =  0.000000000000E+00
//  v2rhobsigmabb =  0.685130252745E-01

//  v2sigmaaa2    =  0.115351801846E-01
//  v2sigmaaaab   =  0.000000000000E+00
//  v2sigmaaabb   =  0.000000000000E+00
//  v2sigmaab2    =  0.000000000000E+00
//  v2sigmaabbb   =  0.000000000000E+00
//  v2sigmabb2    =  0.115351801846E-01
TEST(DfFunctional_Becke88, pointwise7) {
  // input
  const double dRhoA = 0.26E+00;
  const double dRhoB = 0.26E+00;
  const double dGammaAA = 0.28E+00;
  const double dGammaAB = 0.28E+00;
  const double dGammaBB = 0.28E+00;

  // expected value
  const double zk = -0.321148637763E+00;
  const double vRhoA = -0.766539815464E+00;
  const double vRhoB = -0.766539815464E+00;
  const double vGammaAA = -0.198197408319E-01;
  const double vGammaAB = 0.000000000000E+00;
  const double vGammaBB = -0.198197408319E-01;

  // execute test
  DfFunctional_Becke88 f;

  double dFunctionalValue =
      f.getFunctional(dRhoA, dRhoB, dGammaAA, dGammaAB, dGammaBB);
  EXPECT_NEAR(zk, dFunctionalValue, EPS);

  double dRoundF_roundRhoA, dRoundF_roundRhoB;
  double dRoundF_roundGammaAA, dRoundF_roundGammaAB, dRoundF_roundGammaBB;
  f.getDerivativeFunctional(dRhoA, dRhoB, dGammaAA, dGammaAB, dGammaBB,
                            &dRoundF_roundRhoA, &dRoundF_roundRhoB,
                            &dRoundF_roundGammaAA, &dRoundF_roundGammaAB,
                            &dRoundF_roundGammaBB);

  EXPECT_NEAR(vRhoA, dRoundF_roundRhoA, EPS);
  EXPECT_NEAR(vRhoB, dRoundF_roundRhoB, EPS);
  EXPECT_NEAR(vGammaAA, dRoundF_roundGammaAA, EPS);
  EXPECT_NEAR(vGammaAB, dRoundF_roundGammaAB, EPS);
  EXPECT_NEAR(vGammaBB, dRoundF_roundGammaBB, EPS);
}

// test8
//  rhoa= 0.53E+05 rhob= 0.53E+05 sigmaaa= 0.96E+05 sigmaab= 0.96E+05 sigmabb=
//  0.96E+05 zk            = -0.370503980183E+07 vrhoa         =
//  -0.466042742267E+02 vrhob         = -0.466042742267E+02 vsigmaaa      =
//  -0.210967131116E-08 vsigmaab      =  0.000000000000E+00 vsigmabb      =
//  -0.210967131116E-08

//  v2rhoa2       = -0.293108643192E-03
//  v2rhoab       =  0.000000000000E+00
//  v2rhob2       = -0.293108643192E-03

//  v2rhoasigmaaa =  0.530734919752E-13
//  v2rhoasigmaab =  0.000000000000E+00
//  v2rhoasigmabb =  0.000000000000E+00
//  v2rhobsigmaaa =  0.000000000000E+00
//  v2rhobsigmaab =  0.000000000000E+00
//  v2rhobsigmabb =  0.530734919752E-13

//  v2sigmaaa2    =  0.268272614874E-22
//  v2sigmaaaab   =  0.000000000000E+00
//  v2sigmaaabb   =  0.000000000000E+00
//  v2sigmaab2    =  0.000000000000E+00
//  v2sigmaabbb   =  0.000000000000E+00
//  v2sigmabb2    =  0.268272614874E-22
TEST(DfFunctional_Becke88, pointwise8) {
  // input
  const double dRhoA = 0.53E+05;
  const double dRhoB = 0.53E+05;
  const double dGammaAA = 0.96E+05;
  const double dGammaAB = 0.96E+05;
  const double dGammaBB = 0.96E+05;

  // expected value
  const double zk = -0.370503980183E+07;
  const double vRhoA = -0.466042742267E+02;
  const double vRhoB = -0.466042742267E+02;
  const double vGammaAA = -0.210967131116E-08;
  const double vGammaAB = 0.000000000000E+00;
  const double vGammaBB = -0.210967131116E-08;

  // execute test
  DfFunctional_Becke88 f;

  double dFunctionalValue =
      f.getFunctional(dRhoA, dRhoB, dGammaAA, dGammaAB, dGammaBB);
  EXPECT_NEAR(zk, dFunctionalValue, EPS);

  double dRoundF_roundRhoA, dRoundF_roundRhoB;
  double dRoundF_roundGammaAA, dRoundF_roundGammaAB, dRoundF_roundGammaBB;
  f.getDerivativeFunctional(dRhoA, dRhoB, dGammaAA, dGammaAB, dGammaBB,
                            &dRoundF_roundRhoA, &dRoundF_roundRhoB,
                            &dRoundF_roundGammaAA, &dRoundF_roundGammaAB,
                            &dRoundF_roundGammaBB);

  EXPECT_NEAR(vRhoA, dRoundF_roundRhoA, EPS);
  EXPECT_NEAR(vRhoB, dRoundF_roundRhoB, EPS);
  EXPECT_NEAR(vGammaAA, dRoundF_roundGammaAA, EPS);
  EXPECT_NEAR(vGammaAB, dRoundF_roundGammaAB, EPS);
  EXPECT_NEAR(vGammaBB, dRoundF_roundGammaBB, EPS);
}

// test9
//  rhoa= 0.47E+05 rhob= 0.47E+05 sigmaaa= 0.29E+14 sigmaab= 0.29E+14 sigmabb=
//  0.29E+14 zk            = -0.328152696735E+07 vrhoa         =
//  -0.433514250199E+02 vrhob         = -0.433514250199E+02 vsigmaaa      =
//  -0.194182330561E-08 vsigmaab      =  0.000000000000E+00 vsigmabb      =
//  -0.194182330561E-08

//  v2rhoa2       = -0.368694289315E-03
//  v2rhoab       =  0.000000000000E+00
//  v2rhob2       = -0.368694289315E-03

//  v2rhoasigmaaa =  0.372175421272E-13
//  v2rhoasigmaab =  0.000000000000E+00
//  v2rhoasigmabb =  0.000000000000E+00
//  v2rhobsigmaaa =  0.000000000000E+00
//  v2rhobsigmaab =  0.000000000000E+00
//  v2rhobsigmabb =  0.372175421272E-13

//  v2sigmaaa2    =  0.108604300970E-22
//  v2sigmaaaab   =  0.000000000000E+00
//  v2sigmaaabb   =  0.000000000000E+00
//  v2sigmaab2    =  0.000000000000E+00
//  v2sigmaabbb   =  0.000000000000E+00
//  v2sigmabb2    =  0.108604300970E-22
TEST(DfFunctional_Becke88, pointwise9) {
  // input
  const double dRhoA = 0.47E+05;
  const double dRhoB = 0.47E+05;
  const double dGammaAA = 0.29E+14;
  const double dGammaAB = 0.29E+14;
  const double dGammaBB = 0.29E+14;

  // expected value
  const double zk = -0.328152696735E+07;
  const double vRhoA = -0.433514250199E+02;
  const double vRhoB = -0.433514250199E+02;
  const double vGammaAA = -0.194182330561E-08;
  const double vGammaAB = 0.000000000000E+00;
  const double vGammaBB = -0.194182330561E-08;

  // execute test
  DfFunctional_Becke88 f;

  double dFunctionalValue =
      f.getFunctional(dRhoA, dRhoB, dGammaAA, dGammaAB, dGammaBB);
  EXPECT_NEAR(zk, dFunctionalValue, EPS);

  double dRoundF_roundRhoA, dRoundF_roundRhoB;
  double dRoundF_roundGammaAA, dRoundF_roundGammaAB, dRoundF_roundGammaBB;
  f.getDerivativeFunctional(dRhoA, dRhoB, dGammaAA, dGammaAB, dGammaBB,
                            &dRoundF_roundRhoA, &dRoundF_roundRhoB,
                            &dRoundF_roundGammaAA, &dRoundF_roundGammaAB,
                            &dRoundF_roundGammaBB);

  EXPECT_NEAR(vRhoA, dRoundF_roundRhoA, EPS);
  EXPECT_NEAR(vRhoB, dRoundF_roundRhoB, EPS);
  EXPECT_NEAR(vGammaAA, dRoundF_roundGammaAA, EPS);
  EXPECT_NEAR(vGammaAB, dRoundF_roundGammaAB, EPS);
  EXPECT_NEAR(vGammaBB, dRoundF_roundGammaBB, EPS);
}

// test10
//  rhoa= 0.15E+00 rhob= 0.15E+00 sigmaaa= 0.16E+00 sigmaab= 0.16E+00 sigmabb=
//  0.16E+00 zk            = -0.161367392847E+00 vrhoa         =
//  -0.619947650806E+00 vrhob         = -0.619947650806E+00 vsigmaaa      =
//  -0.341862053361E-01 vsigmaab      =  0.000000000000E+00 vsigmabb      =
//  -0.341862053361E-01

//  v2rhoa2       = -0.179936512300E+01
//  v2rhoab       =  0.000000000000E+00
//  v2rhob2       = -0.179936512300E+01

//  v2rhoasigmaaa =  0.148255198864E+00
//  v2rhoasigmaab =  0.000000000000E+00
//  v2rhoasigmabb =  0.000000000000E+00
//  v2rhobsigmaaa =  0.000000000000E+00
//  v2rhobsigmaab =  0.000000000000E+00
//  v2rhobsigmabb =  0.148255198864E+00

//  v2sigmaaa2    =  0.547109233247E-01
//  v2sigmaaaab   =  0.000000000000E+00
//  v2sigmaaabb   =  0.000000000000E+00
//  v2sigmaab2    =  0.000000000000E+00
//  v2sigmaabbb   =  0.000000000000E+00
//  v2sigmabb2    =  0.547109233247E-01
TEST(DfFunctional_Becke88, pointwise10) {
  // input
  const double dRhoA = 0.15E+00;
  const double dRhoB = 0.15E+00;
  const double dGammaAA = 0.16E+00;
  const double dGammaAB = 0.16E+00;
  const double dGammaBB = 0.16E+00;

  // expected value
  const double zk = -0.161367392847E+00;
  const double vRhoA = -0.619947650806E+00;
  const double vRhoB = -0.619947650806E+00;
  const double vGammaAA = -0.341862053361E-01;
  const double vGammaAB = 0.000000000000E+00;
  const double vGammaBB = -0.341862053361E-01;

  // execute test
  DfFunctional_Becke88 f;

  double dFunctionalValue =
      f.getFunctional(dRhoA, dRhoB, dGammaAA, dGammaAB, dGammaBB);
  EXPECT_NEAR(zk, dFunctionalValue, EPS);

  double dRoundF_roundRhoA, dRoundF_roundRhoB;
  double dRoundF_roundGammaAA, dRoundF_roundGammaAB, dRoundF_roundGammaBB;
  f.getDerivativeFunctional(dRhoA, dRhoB, dGammaAA, dGammaAB, dGammaBB,
                            &dRoundF_roundRhoA, &dRoundF_roundRhoB,
                            &dRoundF_roundGammaAA, &dRoundF_roundGammaAB,
                            &dRoundF_roundGammaBB);

  EXPECT_NEAR(vRhoA, dRoundF_roundRhoA, EPS);
  EXPECT_NEAR(vRhoB, dRoundF_roundRhoB, EPS);
  EXPECT_NEAR(vGammaAA, dRoundF_roundGammaAA, EPS);
  EXPECT_NEAR(vGammaAB, dRoundF_roundGammaAB, EPS);
  EXPECT_NEAR(vGammaBB, dRoundF_roundGammaBB, EPS);
}

// test11
//  rhoa= 0.35E+01 rhob= 0.00E+00 sigmaaa= 0.46E-10 sigmaab= 0.00E+00 sigmabb=
//  0.00E+00 zk            = -0.494484233083E+01 vrhoa         =
//  -0.188374945936E+01 vrhob         =  0.000000000000E+00 vsigmaaa      =
//  -0.790360507210E-03 vsigmaab      =  0.000000000000E+00 vsigmabb      =
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

//  v2sigmaaa2    =  0.141061224494E-05
//  v2sigmaaaab   =  0.000000000000E+00
//  v2sigmaaabb   =  0.000000000000E+00
//  v2sigmaab2    =  0.000000000000E+00
//  v2sigmaabbb   =  0.000000000000E+00
//  v2sigmabb2    =  0.000000000000E+00
TEST(DfFunctional_Becke88, pointwise11) {
  // input
  const double dRhoA = 0.35E+01;
  const double dRhoB = 0.00E+00;
  const double dGammaAA = 0.46E-10;
  const double dGammaAB = 0.00E+00;
  const double dGammaBB = 0.00E+00;

  // expected value
  const double zk = -0.494484233083E+01;
  const double vRhoA = -0.188374945936E+01;
  const double vRhoB = 0.000000000000E+00;
  const double vGammaAA = -0.790360507210E-03;
  const double vGammaAB = 0.000000000000E+00;
  const double vGammaBB = 0.000000000000E+00;

  // execute test
  DfFunctional_Becke88 f;

  double dFunctionalValue =
      f.getFunctional(dRhoA, dRhoB, dGammaAA, dGammaAB, dGammaBB);
  EXPECT_NEAR(zk, dFunctionalValue, EPS);

  double dRoundF_roundRhoA, dRoundF_roundRhoB;
  double dRoundF_roundGammaAA, dRoundF_roundGammaAB, dRoundF_roundGammaBB;
  f.getDerivativeFunctional(dRhoA, dRhoB, dGammaAA, dGammaAB, dGammaBB,
                            &dRoundF_roundRhoA, &dRoundF_roundRhoB,
                            &dRoundF_roundGammaAA, &dRoundF_roundGammaAB,
                            &dRoundF_roundGammaBB);

  EXPECT_NEAR(vRhoA, dRoundF_roundRhoA, EPS);
  EXPECT_NEAR(vRhoB, dRoundF_roundRhoB, EPS);
  EXPECT_NEAR(vGammaAA, dRoundF_roundGammaAA, EPS);
  EXPECT_NEAR(vGammaAB, dRoundF_roundGammaAB, EPS);
  EXPECT_NEAR(vGammaBB, dRoundF_roundGammaBB, EPS);
}

// test12
//  rhoa= 0.35E+01 rhob= 0.00E+00 sigmaaa= 0.34E+01 sigmaab= 0.00E+00 sigmabb=
//  0.00E+00 zk            = -0.494752158228E+01 vrhoa         =
//  -0.188273473796E+01 vrhob         =  0.000000000000E+00 vsigmaaa      =
//  -0.785719874797E-03 vsigmaab      =  0.000000000000E+00 vsigmabb      =
//  0.000000000000E+00

//  v2rhoa2       = -0.180074584533E+00
//  v2rhoab       =  0.000000000000E+00
//  v2rhob2       =  0.000000000000E+00

//  v2rhoasigmaaa =  0.000000000000E+00
//  v2rhoasigmaab =  0.000000000000E+00
//  v2rhoasigmabb =  0.000000000000E+00
//  v2rhobsigmaaa =  0.000000000000E+00
//  v2rhobsigmaab =  0.000000000000E+00
//  v2rhobsigmabb =  0.000000000000E+00

//  v2sigmaaa2    =  0.132207800833E-05
//  v2sigmaaaab   =  0.000000000000E+00
//  v2sigmaaabb   =  0.000000000000E+00
//  v2sigmaab2    =  0.000000000000E+00
//  v2sigmaabbb   =  0.000000000000E+00
//  v2sigmabb2    =  0.000000000000E+00
TEST(DfFunctional_Becke88, pointwise12) {
  // input
  const double dRhoA = 0.35E+01;
  const double dRhoB = 0.00E+00;
  const double dGammaAA = 0.34E+01;
  const double dGammaAB = 0.00E+00;
  const double dGammaBB = 0.00E+00;

  // expected value
  const double zk = -0.494752158228E+01;
  const double vRhoA = -0.188273473796E+01;
  const double vRhoB = 0.000000000000E+00;
  const double vGammaAA = -0.785719874797E-03;
  const double vGammaAB = 0.000000000000E+00;
  const double vGammaBB = 0.000000000000E+00;

  // execute test
  DfFunctional_Becke88 f;

  double dFunctionalValue =
      f.getFunctional(dRhoA, dRhoB, dGammaAA, dGammaAB, dGammaBB);
  EXPECT_NEAR(zk, dFunctionalValue, EPS);

  double dRoundF_roundRhoA, dRoundF_roundRhoB;
  double dRoundF_roundGammaAA, dRoundF_roundGammaAB, dRoundF_roundGammaBB;
  f.getDerivativeFunctional(dRhoA, dRhoB, dGammaAA, dGammaAB, dGammaBB,
                            &dRoundF_roundRhoA, &dRoundF_roundRhoB,
                            &dRoundF_roundGammaAA, &dRoundF_roundGammaAB,
                            &dRoundF_roundGammaBB);

  EXPECT_NEAR(vRhoA, dRoundF_roundRhoA, EPS);
  EXPECT_NEAR(vRhoB, dRoundF_roundRhoB, EPS);
  EXPECT_NEAR(vGammaAA, dRoundF_roundGammaAA, EPS);
  EXPECT_NEAR(vGammaAB, dRoundF_roundGammaAB, EPS);
  EXPECT_NEAR(vGammaBB, dRoundF_roundGammaBB, EPS);
}

//  rhoa= 0.30E+01 rhob= 0.00E+00 sigmaaa= 0.20E+03 sigmaab= 0.00E+00 sigmabb=
//  0.00E+00 zk            = -0.419401965265E+01 vrhoa         =
//  -0.172996985225E+01 vrhob         =  0.000000000000E+00 vsigmaaa      =
//  -0.753968712720E-03 vsigmaab      =  0.000000000000E+00 vsigmabb      =
//  0.000000000000E+00

//  v2rhoa2       = -0.231847776994E+00
//  v2rhoab       =  0.000000000000E+00
//  v2rhob2       =  0.000000000000E+00

//  v2rhoasigmaaa =  0.000000000000E+00
//  v2rhoasigmaab =  0.000000000000E+00
//  v2rhoasigmabb =  0.000000000000E+00
//  v2rhobsigmaaa =  0.000000000000E+00
//  v2rhobsigmaab =  0.000000000000E+00
//  v2rhobsigmabb =  0.000000000000E+00

//  v2sigmaaa2    =  0.631038474677E-06
//  v2sigmaaaab   =  0.000000000000E+00
//  v2sigmaaabb   =  0.000000000000E+00
//  v2sigmaab2    =  0.000000000000E+00
//  v2sigmaabbb   =  0.000000000000E+00
//  v2sigmabb2    =  0.000000000000E+00
TEST(DfFunctional_Becke88, pointwise13) {
  // input
  const double dRhoA = 0.30E+01;
  const double dRhoB = 0.00E+00;
  const double dGammaAA = 0.20E+03;
  const double dGammaAB = 0.00E+00;
  const double dGammaBB = 0.00E+00;

  // expected value
  const double zk = -0.419401965265E+01;
  const double vRhoA = -0.172996985225E+01;
  const double vRhoB = 0.000000000000E+00;
  const double vGammaAA = -0.753968712720E-03;
  const double vGammaAB = 0.000000000000E+00;
  const double vGammaBB = 0.000000000000E+00;

  // execute test
  DfFunctional_Becke88 f;

  double dFunctionalValue =
      f.getFunctional(dRhoA, dRhoB, dGammaAA, dGammaAB, dGammaBB);
  EXPECT_NEAR(zk, dFunctionalValue, EPS);

  double dRoundF_roundRhoA, dRoundF_roundRhoB;
  double dRoundF_roundGammaAA, dRoundF_roundGammaAB, dRoundF_roundGammaBB;
  f.getDerivativeFunctional(dRhoA, dRhoB, dGammaAA, dGammaAB, dGammaBB,
                            &dRoundF_roundRhoA, &dRoundF_roundRhoB,
                            &dRoundF_roundGammaAA, &dRoundF_roundGammaAB,
                            &dRoundF_roundGammaBB);

  EXPECT_NEAR(vRhoA, dRoundF_roundRhoA, EPS);
  EXPECT_NEAR(vRhoB, dRoundF_roundRhoB, EPS);
  EXPECT_NEAR(vGammaAA, dRoundF_roundGammaAA, EPS);
  EXPECT_NEAR(vGammaAB, dRoundF_roundGammaAB, EPS);
  EXPECT_NEAR(vGammaBB, dRoundF_roundGammaBB, EPS);
}

//  rhoa= 0.58E-01 rhob= 0.00E+00 sigmaaa= 0.47E-01 sigmaab= 0.00E+00 sigmabb=
//  0.00E+00 zk            = -0.259998774808E-01 vrhoa         =
//  -0.428541964579E+00 vrhob         =  0.000000000000E+00 vsigmaaa      =
//  -0.782798087404E-01 vsigmaab      =  0.000000000000E+00 vsigmabb      =
//  0.000000000000E+00

//  v2rhoa2       = -0.312635380295E+01
//  v2rhoab       =  0.000000000000E+00
//  v2rhob2       =  0.000000000000E+00

//  v2rhoasigmaaa =  0.000000000000E+00
//  v2rhoasigmaab =  0.000000000000E+00
//  v2rhoasigmabb =  0.000000000000E+00
//  v2rhobsigmaaa =  0.000000000000E+00
//  v2rhobsigmaab =  0.000000000000E+00
//  v2rhobsigmabb =  0.000000000000E+00

//  v2sigmaaa2    =  0.690680500539E+00
//  v2sigmaaaab   =  0.000000000000E+00
//  v2sigmaaabb   =  0.000000000000E+00
//  v2sigmaab2    =  0.000000000000E+00
//  v2sigmaabbb   =  0.000000000000E+00
//  v2sigmabb2    =  0.000000000000E+00
TEST(DfFunctional_Becke88, pointwise14) {
  // input
  const double dRhoA = 0.58E-01;
  const double dRhoB = 0.00E+00;
  const double dGammaAA = 0.47E-01;
  const double dGammaAB = 0.00E+00;
  const double dGammaBB = 0.00E+00;

  // expected value
  const double zk = -0.259998774808E-01;
  const double vRhoA = -0.428541964579E+00;
  const double vRhoB = 0.000000000000E+00;
  const double vGammaAA = -0.782798087404E-01;
  const double vGammaAB = 0.000000000000E+00;
  const double vGammaBB = 0.000000000000E+00;

  // execute test
  DfFunctional_Becke88 f;

  double dFunctionalValue =
      f.getFunctional(dRhoA, dRhoB, dGammaAA, dGammaAB, dGammaBB);
  EXPECT_NEAR(zk, dFunctionalValue, EPS);

  double dRoundF_roundRhoA, dRoundF_roundRhoB;
  double dRoundF_roundGammaAA, dRoundF_roundGammaAB, dRoundF_roundGammaBB;
  f.getDerivativeFunctional(dRhoA, dRhoB, dGammaAA, dGammaAB, dGammaBB,
                            &dRoundF_roundRhoA, &dRoundF_roundRhoB,
                            &dRoundF_roundGammaAA, &dRoundF_roundGammaAB,
                            &dRoundF_roundGammaBB);

  EXPECT_NEAR(vRhoA, dRoundF_roundRhoA, EPS);
  EXPECT_NEAR(vRhoB, dRoundF_roundRhoB, EPS);
  EXPECT_NEAR(vGammaAA, dRoundF_roundGammaAA, EPS);
  EXPECT_NEAR(vGammaAB, dRoundF_roundGammaAB, EPS);
  EXPECT_NEAR(vGammaBB, dRoundF_roundGammaBB, EPS);
}

//  rhoa= 0.82E+02 rhob= 0.81E+02 sigmaaa= 0.49E+07 sigmaab= 0.49E+07 sigmabb=
//  0.49E+07 zk            = -0.740814331850E+03 vrhoa         =
//  -0.498247046111E+01 vrhob         = -0.495532429311E+01 vsigmaaa      =
//  -0.678262141094E-05 vsigmaab      =  0.000000000000E+00 vsigmabb      =
//  -0.682517937334E-05

//  v2rhoa2       = -0.270432929879E-01
//  v2rhoab       =  0.000000000000E+00
//  v2rhob2       = -0.272490060294E-01

//  v2rhoasigmaaa =  0.426066228043E-07
//  v2rhoasigmaab =  0.000000000000E+00
//  v2rhoasigmabb =  0.000000000000E+00
//  v2rhobsigmaaa =  0.000000000000E+00
//  v2rhobsigmaab =  0.000000000000E+00
//  v2rhobsigmabb =  0.425046982664E-07

//  v2sigmaaa2    =  0.424725929436E-12
//  v2sigmaaaab   =  0.000000000000E+00
//  v2sigmaaabb   =  0.000000000000E+00
//  v2sigmaab2    =  0.000000000000E+00
//  v2sigmaabbb   =  0.000000000000E+00
//  v2sigmabb2    =  0.432961117719E-12
TEST(DfFunctional_Becke88, pointwise15) {
  // input
  const double dRhoA = 0.82E+02;
  const double dRhoB = 0.81E+02;
  const double dGammaAA = 0.49E+07;
  const double dGammaAB = 0.49E+07;
  const double dGammaBB = 0.49E+07;

  // expected value
  const double zk = -0.740814331850E+03;
  const double vRhoA = -0.498247046111E+01;
  const double vRhoB = -0.495532429311E+01;
  const double vGammaAA = -0.678262141094E-05;
  const double vGammaAB = 0.000000000000E+00;
  const double vGammaBB = -0.682517937334E-05;

  // execute test
  DfFunctional_Becke88 f;

  double dFunctionalValue =
      f.getFunctional(dRhoA, dRhoB, dGammaAA, dGammaAB, dGammaBB);
  EXPECT_NEAR(zk, dFunctionalValue, EPS);

  double dRoundF_roundRhoA, dRoundF_roundRhoB;
  double dRoundF_roundGammaAA, dRoundF_roundGammaAB, dRoundF_roundGammaBB;
  f.getDerivativeFunctional(dRhoA, dRhoB, dGammaAA, dGammaAB, dGammaBB,
                            &dRoundF_roundRhoA, &dRoundF_roundRhoB,
                            &dRoundF_roundGammaAA, &dRoundF_roundGammaAB,
                            &dRoundF_roundGammaBB);

  EXPECT_NEAR(vRhoA, dRoundF_roundRhoA, EPS);
  EXPECT_NEAR(vRhoB, dRoundF_roundRhoB, EPS);
  EXPECT_NEAR(vGammaAA, dRoundF_roundGammaAA, EPS);
  EXPECT_NEAR(vGammaAB, dRoundF_roundGammaAB, EPS);
  EXPECT_NEAR(vGammaBB, dRoundF_roundGammaBB, EPS);
}

//  rhoa= 0.39E+02 rhob= 0.38E+02 sigmaaa= 0.81E+06 sigmaab= 0.82E+06 sigmabb=
//  0.82E+06 zk            = -0.277987329958E+03 vrhoa         =
//  -0.385951846654E+01 vrhob         = -0.381309494319E+01 vsigmaaa      =
//  -0.172434478018E-04 vsigmaab      =  0.000000000000E+00 vsigmabb      =
//  -0.173712338362E-04

//  v2rhoa2       = -0.441426807406E-01
//  v2rhoab       =  0.000000000000E+00
//  v2rhob2       = -0.447245742260E-01

//  v2rhoasigmaaa =  0.201415922856E-06
//  v2rhoasigmaab =  0.000000000000E+00
//  v2rhoasigmabb =  0.000000000000E+00
//  v2rhobsigmaaa =  0.000000000000E+00
//  v2rhobsigmaab =  0.000000000000E+00
//  v2rhobsigmabb =  0.195961359539E-06

//  v2sigmaaa2    =  0.700742719647E-11
//  v2sigmaaaab   =  0.000000000000E+00
//  v2sigmaaabb   =  0.000000000000E+00
//  v2sigmaab2    =  0.000000000000E+00
//  v2sigmaabbb   =  0.000000000000E+00
//  v2sigmabb2    =  0.718678968862E-11
TEST(DfFunctional_Becke88, pointwise16) {
  // input
  const double dRhoA = 0.39E+02;
  const double dRhoB = 0.38E+02;
  const double dGammaAA = 0.81E+06;
  const double dGammaAB = 0.82E+06;
  const double dGammaBB = 0.82E+06;

  // expected value
  const double zk = -0.277987329958E+03;
  const double vRhoA = -0.385951846654E+01;
  const double vRhoB = -0.381309494319E+01;
  const double vGammaAA = -0.172434478018E-04;
  const double vGammaAB = 0.000000000000E+00;
  const double vGammaBB = -0.173712338362E-04;

  // execute test
  DfFunctional_Becke88 f;

  double dFunctionalValue =
      f.getFunctional(dRhoA, dRhoB, dGammaAA, dGammaAB, dGammaBB);
  EXPECT_NEAR(zk, dFunctionalValue, EPS);

  double dRoundF_roundRhoA, dRoundF_roundRhoB;
  double dRoundF_roundGammaAA, dRoundF_roundGammaAB, dRoundF_roundGammaBB;
  f.getDerivativeFunctional(dRhoA, dRhoB, dGammaAA, dGammaAB, dGammaBB,
                            &dRoundF_roundRhoA, &dRoundF_roundRhoB,
                            &dRoundF_roundGammaAA, &dRoundF_roundGammaAB,
                            &dRoundF_roundGammaBB);

  EXPECT_NEAR(vRhoA, dRoundF_roundRhoA, EPS);
  EXPECT_NEAR(vRhoB, dRoundF_roundRhoB, EPS);
  EXPECT_NEAR(vGammaAA, dRoundF_roundGammaAA, EPS);
  EXPECT_NEAR(vGammaAB, dRoundF_roundGammaAB, EPS);
  EXPECT_NEAR(vGammaBB, dRoundF_roundGammaBB, EPS);
}

//  rhoa= 0.13E+00 rhob= 0.95E-01 sigmaaa= 0.15E+00 sigmaab= 0.18E+00 sigmabb=
//  0.22E+00 zk            = -0.120208982576E+00 vrhoa         =
//  -0.583637510333E+00 vrhob         = -0.501672871724E+00 vsigmaaa      =
//  -0.379227871606E-01 vsigmaab      =  0.000000000000E+00 vsigmabb      =
//  -0.367802205908E-01

//  v2rhoa2       = -0.199082326261E+01
//  v2rhoab       =  0.000000000000E+00
//  v2rhob2       = -0.213690115573E+01

//  v2rhoasigmaaa =  0.160652968405E+00
//  v2rhoasigmaab =  0.000000000000E+00
//  v2rhoasigmabb =  0.000000000000E+00
//  v2rhobsigmaaa =  0.000000000000E+00
//  v2rhobsigmaab =  0.000000000000E+00
//  v2rhobsigmabb =  0.609908850345E-01

//  v2sigmaaa2    =  0.741970758035E-01
//  v2sigmaaaab   =  0.000000000000E+00
//  v2sigmaaabb   =  0.000000000000E+00
//  v2sigmaab2    =  0.000000000000E+00
//  v2sigmaabbb   =  0.000000000000E+00
//  v2sigmabb2    =  0.737150455275E-01
TEST(DfFunctional_Becke88, pointwise17) {
  // input
  const double dRhoA = 0.13E+00;
  const double dRhoB = 0.95E-01;
  const double dGammaAA = 0.15E+00;
  const double dGammaAB = 0.18E+00;
  const double dGammaBB = 0.22E+00;

  // expected value
  const double zk = -0.120208982576E+00;
  const double vRhoA = -0.583637510333E+00;
  const double vRhoB = -0.501672871724E+00;
  const double vGammaAA = -0.379227871606E-01;
  const double vGammaAB = 0.000000000000E+00;
  const double vGammaBB = -0.367802205908E-01;

  // execute test
  DfFunctional_Becke88 f;

  double dFunctionalValue =
      f.getFunctional(dRhoA, dRhoB, dGammaAA, dGammaAB, dGammaBB);
  EXPECT_NEAR(zk, dFunctionalValue, EPS);

  double dRoundF_roundRhoA, dRoundF_roundRhoB;
  double dRoundF_roundGammaAA, dRoundF_roundGammaAB, dRoundF_roundGammaBB;
  f.getDerivativeFunctional(dRhoA, dRhoB, dGammaAA, dGammaAB, dGammaBB,
                            &dRoundF_roundRhoA, &dRoundF_roundRhoB,
                            &dRoundF_roundGammaAA, &dRoundF_roundGammaAB,
                            &dRoundF_roundGammaBB);

  EXPECT_NEAR(vRhoA, dRoundF_roundRhoA, EPS);
  EXPECT_NEAR(vRhoB, dRoundF_roundRhoB, EPS);
  EXPECT_NEAR(vGammaAA, dRoundF_roundGammaAA, EPS);
  EXPECT_NEAR(vGammaAB, dRoundF_roundGammaAB, EPS);
  EXPECT_NEAR(vGammaBB, dRoundF_roundGammaBB, EPS);
}

//  rhoa= 0.78E-01 rhob= 0.31E-01 sigmaaa= 0.41E-02 sigmaab= 0.38E-02 sigmabb=
//  0.36E-02 zk            = -0.416730804241E-01 vrhoa         =
//  -0.522700471817E+00 vrhob         = -0.360524872778E+00 vsigmaaa      =
//  -0.111846105258E+00 vsigmaab      =  0.000000000000E+00 vsigmabb      =
//  -0.249411314940E+00

//  v2rhoa2       = -0.245380935846E+01
//  v2rhoab       =  0.000000000000E+00
//  v2rhob2       = -0.517371957290E+01

//  v2rhoasigmaaa =  0.156984506289E+01
//  v2rhoasigmaab =  0.000000000000E+00
//  v2rhoasigmabb =  0.000000000000E+00
//  v2rhobsigmaaa =  0.000000000000E+00
//  v2rhobsigmaab =  0.000000000000E+00
//  v2rhobsigmabb =  0.418857803825E+01

//  v2sigmaaa2    =  0.244026452179E+01
//  v2sigmaaaab   =  0.000000000000E+00
//  v2sigmaaabb   =  0.000000000000E+00
//  v2sigmaab2    =  0.000000000000E+00
//  v2sigmaabbb   =  0.000000000000E+00
//  v2sigmabb2    =  0.211148438264E+02
TEST(DfFunctional_Becke88, pointwise18) {
  // input
  const double dRhoA = 0.78E-01;
  const double dRhoB = 0.31E-01;
  const double dGammaAA = 0.41E-02;
  const double dGammaAB = 0.38E-02;
  const double dGammaBB = 0.36E-02;

  // expected value
  const double zk = -0.416730804241E-01;
  const double vRhoA = -0.522700471817E+00;
  const double vRhoB = -0.360524872778E+00;
  const double vGammaAA = -0.111846105258E+00;
  const double vGammaAB = 0.000000000000E+00;
  const double vGammaBB = -0.249411314940E+00;

  // execute test
  DfFunctional_Becke88 f;

  double dFunctionalValue =
      f.getFunctional(dRhoA, dRhoB, dGammaAA, dGammaAB, dGammaBB);
  EXPECT_NEAR(zk, dFunctionalValue, EPS);

  double dRoundF_roundRhoA, dRoundF_roundRhoB;
  double dRoundF_roundGammaAA, dRoundF_roundGammaAB, dRoundF_roundGammaBB;
  f.getDerivativeFunctional(dRhoA, dRhoB, dGammaAA, dGammaAB, dGammaBB,
                            &dRoundF_roundRhoA, &dRoundF_roundRhoB,
                            &dRoundF_roundGammaAA, &dRoundF_roundGammaAB,
                            &dRoundF_roundGammaBB);

  EXPECT_NEAR(vRhoA, dRoundF_roundRhoA, EPS);
  EXPECT_NEAR(vRhoB, dRoundF_roundRhoB, EPS);
  EXPECT_NEAR(vGammaAA, dRoundF_roundGammaAA, EPS);
  EXPECT_NEAR(vGammaAB, dRoundF_roundGammaAB, EPS);
  EXPECT_NEAR(vGammaBB, dRoundF_roundGammaBB, EPS);
}

//  rhoa= 0.50E+02 rhob= 0.49E+02 sigmaaa= 0.11E+06 sigmaab= 0.11E+06 sigmabb=
//  0.11E+06 zk            = -0.343037899309E+03 vrhoa         =
//  -0.451375531283E+01 vrhob         = -0.448071957650E+01 vsigmaaa      =
//  -0.204625762085E-04 vsigmaab      =  0.000000000000E+00 vsigmabb      =
//  -0.209266539812E-04

//  v2rhoa2       = -0.327661718483E-01
//  v2rhoab       =  0.000000000000E+00
//  v2rhob2       = -0.333090761894E-01

//  v2rhoasigmaaa =  0.455875527750E-06
//  v2rhoasigmaab =  0.000000000000E+00
//  v2rhoasigmabb =  0.000000000000E+00
//  v2rhobsigmaaa =  0.000000000000E+00
//  v2rhobsigmaab =  0.000000000000E+00
//  v2rhobsigmabb =  0.472402981076E-06

//  v2sigmaaa2    =  0.153056541719E-10
//  v2sigmaaaab   =  0.000000000000E+00
//  v2sigmaaabb   =  0.000000000000E+00
//  v2sigmaab2    =  0.000000000000E+00
//  v2sigmaabbb   =  0.000000000000E+00
//  v2sigmabb2    =  0.162083837575E-10
TEST(DfFunctional_Becke88, pointwise19) {
  // input
  const double dRhoA = 0.50E+02;
  const double dRhoB = 0.49E+02;
  const double dGammaAA = 0.11E+06;
  const double dGammaAB = 0.11E+06;
  const double dGammaBB = 0.11E+06;

  // expected value
  const double zk = -0.343037899309E+03;
  const double vRhoA = -0.451375531283E+01;
  const double vRhoB = -0.448071957650E+01;
  const double vGammaAA = -0.204625762085E-04;
  const double vGammaAB = 0.000000000000E+00;
  const double vGammaBB = -0.209266539812E-04;

  // execute test
  DfFunctional_Becke88 f;

  double dFunctionalValue =
      f.getFunctional(dRhoA, dRhoB, dGammaAA, dGammaAB, dGammaBB);
  EXPECT_NEAR(zk, dFunctionalValue, EPS);

  double dRoundF_roundRhoA, dRoundF_roundRhoB;
  double dRoundF_roundGammaAA, dRoundF_roundGammaAB, dRoundF_roundGammaBB;
  f.getDerivativeFunctional(dRhoA, dRhoB, dGammaAA, dGammaAB, dGammaBB,
                            &dRoundF_roundRhoA, &dRoundF_roundRhoB,
                            &dRoundF_roundGammaAA, &dRoundF_roundGammaAB,
                            &dRoundF_roundGammaBB);

  EXPECT_NEAR(vRhoA, dRoundF_roundRhoA, EPS);
  EXPECT_NEAR(vRhoB, dRoundF_roundRhoB, EPS);
  EXPECT_NEAR(vGammaAA, dRoundF_roundGammaAA, EPS);
  EXPECT_NEAR(vGammaAB, dRoundF_roundGammaAB, EPS);
  EXPECT_NEAR(vGammaBB, dRoundF_roundGammaBB, EPS);
}

//  rhoa= 0.40E+02 rhob= 0.40E+02 sigmaaa= 0.99E+05 sigmaab= 0.98E+05 sigmabb=
//  0.98E+05 zk            = -0.260133861611E+03 vrhoa         =
//  -0.416254426190E+01 vrhob         = -0.416322526434E+01 vsigmaaa      =
//  -0.262815715798E-04 vsigmaab      =  0.000000000000E+00 vsigmabb      =
//  -0.263113502629E-04

//  v2rhoa2       = -0.391759770307E-01
//  v2rhoab       =  0.000000000000E+00
//  v2rhob2       = -0.391492164486E-01

//  v2rhoasigmaaa =  0.680016391145E-06
//  v2rhoasigmaab =  0.000000000000E+00
//  v2rhoasigmabb =  0.000000000000E+00
//  v2rhobsigmaaa =  0.000000000000E+00
//  v2rhobsigmaab =  0.000000000000E+00
//  v2rhobsigmabb =  0.681990700893E-06

//  v2sigmaaa2    =  0.297024234619E-10
//  v2sigmaaaab   =  0.000000000000E+00
//  v2sigmaaabb   =  0.000000000000E+00
//  v2sigmaab2    =  0.000000000000E+00
//  v2sigmaabbb   =  0.000000000000E+00
//  v2sigmabb2    =  0.298552512045E-10
TEST(DfFunctional_Becke88, pointwise20) {
  // input
  const double dRhoA = 0.40E+02;
  const double dRhoB = 0.40E+02;
  const double dGammaAA = 0.99E+05;
  const double dGammaAB = 0.98E+05;
  const double dGammaBB = 0.98E+05;

  // expected value
  const double zk = -0.260133861611E+03;
  const double vRhoA = -0.416254426190E+01;
  const double vRhoB = -0.416322526434E+01;
  const double vGammaAA = -0.262815715798E-04;
  const double vGammaAB = 0.000000000000E+00;
  const double vGammaBB = -0.263113502629E-04;

  // execute test
  DfFunctional_Becke88 f;

  double dFunctionalValue =
      f.getFunctional(dRhoA, dRhoB, dGammaAA, dGammaAB, dGammaBB);
  EXPECT_NEAR(zk, dFunctionalValue, EPS);

  double dRoundF_roundRhoA, dRoundF_roundRhoB;
  double dRoundF_roundGammaAA, dRoundF_roundGammaAB, dRoundF_roundGammaBB;
  f.getDerivativeFunctional(dRhoA, dRhoB, dGammaAA, dGammaAB, dGammaBB,
                            &dRoundF_roundRhoA, &dRoundF_roundRhoB,
                            &dRoundF_roundGammaAA, &dRoundF_roundGammaAB,
                            &dRoundF_roundGammaBB);

  EXPECT_NEAR(vRhoA, dRoundF_roundRhoA, EPS);
  EXPECT_NEAR(vRhoB, dRoundF_roundRhoB, EPS);
  EXPECT_NEAR(vGammaAA, dRoundF_roundGammaAA, EPS);
  EXPECT_NEAR(vGammaAB, dRoundF_roundGammaAB, EPS);
  EXPECT_NEAR(vGammaBB, dRoundF_roundGammaBB, EPS);
}

//  rhoa= 0.12E+00 rhob= 0.10E+00 sigmaaa= 0.12E+00 sigmaab= 0.13E+00 sigmabb=
//  0.14E+00 zk            = -0.112603177863E+00 vrhoa         =
//  -0.568498809942E+00 vrhob         = -0.520860206545E+00 vsigmaaa      =
//  -0.423139064219E-01 vsigmaab      =  0.000000000000E+00 vsigmabb      =
//  -0.436372569148E-01

//  v2rhoa2       = -0.209994430549E+01
//  v2rhoab       =  0.000000000000E+00
//  v2rhob2       = -0.229642938280E+01

//  v2rhoasigmaaa =  0.195292854203E+00
//  v2rhoasigmaab =  0.000000000000E+00
//  v2rhoasigmabb =  0.000000000000E+00
//  v2rhobsigmaaa =  0.000000000000E+00
//  v2rhobsigmaab =  0.000000000000E+00
//  v2rhobsigmabb =  0.150061257405E+00

//  v2sigmaaa2    =  0.103073123098E+00
//  v2sigmaaaab   =  0.000000000000E+00
//  v2sigmaaabb   =  0.000000000000E+00
//  v2sigmaab2    =  0.000000000000E+00
//  v2sigmaabbb   =  0.000000000000E+00
//  v2sigmabb2    =  0.115652366462E+00
TEST(DfFunctional_Becke88, pointwise21) {
  // input
  const double dRhoA = 0.12E+00;
  const double dRhoB = 0.10E+00;
  const double dGammaAA = 0.12E+00;
  const double dGammaAB = 0.13E+00;
  const double dGammaBB = 0.14E+00;

  // expected value
  const double zk = -0.112603177863E+00;
  const double vRhoA = -0.568498809942E+00;
  const double vRhoB = -0.520860206545E+00;
  const double vGammaAA = -0.423139064219E-01;
  const double vGammaAB = 0.000000000000E+00;
  const double vGammaBB = -0.436372569148E-01;

  // execute test
  DfFunctional_Becke88 f;

  double dFunctionalValue =
      f.getFunctional(dRhoA, dRhoB, dGammaAA, dGammaAB, dGammaBB);
  EXPECT_NEAR(zk, dFunctionalValue, EPS);

  double dRoundF_roundRhoA, dRoundF_roundRhoB;
  double dRoundF_roundGammaAA, dRoundF_roundGammaAB, dRoundF_roundGammaBB;
  f.getDerivativeFunctional(dRhoA, dRhoB, dGammaAA, dGammaAB, dGammaBB,
                            &dRoundF_roundRhoA, &dRoundF_roundRhoB,
                            &dRoundF_roundGammaAA, &dRoundF_roundGammaAB,
                            &dRoundF_roundGammaBB);

  EXPECT_NEAR(vRhoA, dRoundF_roundRhoA, EPS);
  EXPECT_NEAR(vRhoB, dRoundF_roundRhoB, EPS);
  EXPECT_NEAR(vGammaAA, dRoundF_roundGammaAA, EPS);
  EXPECT_NEAR(vGammaAB, dRoundF_roundGammaAB, EPS);
  EXPECT_NEAR(vGammaBB, dRoundF_roundGammaBB, EPS);
}

//  rhoa= 0.48E-01 rhob= 0.25E-01 sigmaaa= 0.46E-02 sigmaab= 0.44E-02 sigmabb=
//  0.41E-02 zk            = -0.253983617946E-01 vrhoa         =
//  -0.431645138108E+00 vrhob         = -0.325991019475E+00 vsigmaaa      =
//  -0.175455092294E+00 vsigmaab      =  0.000000000000E+00 vsigmabb      =
//  -0.260075409570E+00

//  v2rhoa2       = -0.374315005133E+01
//  v2rhoab       =  0.000000000000E+00
//  v2rhob2       = -0.566470285703E+01

//  v2rhoasigmaaa =  0.291762144791E+01
//  v2rhoasigmaab =  0.000000000000E+00
//  v2rhoasigmabb =  0.000000000000E+00
//  v2rhobsigmaaa =  0.000000000000E+00
//  v2rhobsigmaab =  0.000000000000E+00
//  v2rhobsigmabb =  0.301407606104E+01

//  v2sigmaaa2    =  0.765442610536E+01
//  v2sigmaaaab   =  0.000000000000E+00
//  v2sigmaaabb   =  0.000000000000E+00
//  v2sigmaab2    =  0.000000000000E+00
//  v2sigmaabbb   =  0.000000000000E+00
//  v2sigmabb2    =  0.248245711494E+02
TEST(DfFunctional_Becke88, pointwise22) {
  // input
  const double dRhoA = 0.48E-01;
  const double dRhoB = 0.25E-01;
  const double dGammaAA = 0.46E-02;
  const double dGammaAB = 0.44E-02;
  const double dGammaBB = 0.41E-02;

  // expected value
  const double zk = -0.253983617946E-01;
  const double vRhoA = -0.431645138108E+00;
  const double vRhoB = -0.325991019475E+00;
  const double vGammaAA = -0.175455092294E+00;
  const double vGammaAB = 0.000000000000E+00;
  const double vGammaBB = -0.260075409570E+00;

  // execute test
  DfFunctional_Becke88 f;

  double dFunctionalValue =
      f.getFunctional(dRhoA, dRhoB, dGammaAA, dGammaAB, dGammaBB);
  EXPECT_NEAR(zk, dFunctionalValue, EPS);

  double dRoundF_roundRhoA, dRoundF_roundRhoB;
  double dRoundF_roundGammaAA, dRoundF_roundGammaAB, dRoundF_roundGammaBB;
  f.getDerivativeFunctional(dRhoA, dRhoB, dGammaAA, dGammaAB, dGammaBB,
                            &dRoundF_roundRhoA, &dRoundF_roundRhoB,
                            &dRoundF_roundGammaAA, &dRoundF_roundGammaAB,
                            &dRoundF_roundGammaBB);

  EXPECT_NEAR(vRhoA, dRoundF_roundRhoA, EPS);
  EXPECT_NEAR(vRhoB, dRoundF_roundRhoB, EPS);
  EXPECT_NEAR(vGammaAA, dRoundF_roundGammaAA, EPS);
  EXPECT_NEAR(vGammaAB, dRoundF_roundGammaAB, EPS);
  EXPECT_NEAR(vGammaBB, dRoundF_roundGammaBB, EPS);
}
