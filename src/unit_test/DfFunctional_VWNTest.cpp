#include <vector>
#include <string>
#include "digit.h"
#include "DfFunctional_VWNTest.h"

#define EFFCTIVE_DIGITS 12
#define ED_CHECK(x) (std::pow(10.0, -EFFCTIVE_DIGITS + digit(x)))

const double DfFunctional_VWNTest::EPS = 1.0E-10;

// CPPUNIT_ASSERT( condition );
// conditionが偽(false,0)であったとき、失敗します。
//
// CPPUNIT_ASSERT_MESSAGE( message, condition );
// conditionが偽であったとき、失敗します。このときmessageを出力します。
//
// CPPUNIT_FAIL( message );
// 必ず失敗します。messageを出力します。
//
// CPPUNIT_ASSERT_EQUAL( expected, actual );
// 得られた結果actualが期待する値expectedでなかったとき、すなわちexpected != actualのときに失敗します。
//
// CPPUNIT_ASSERT_EQUAL_MESSAGE( message, expected, actual );
// 得られた結果actualが期待する値expectedでなかったとき、すなわちexpected != actualのときに失敗します。このときmessageを出力します。
//
// CPPUNIT_ASSERT_DOUBLES_EQUAL( expected, actual, delta );
// 得られた結果actualと期待する値expectedとの差がdeltaより大きいとき、失敗します。

// =====================================================================
// data used from: http://www.cse.scitech.ac.uk/ccg/dft/data_pt_x_lda.html
//

void DfFunctional_VWNTest::testConstructer(){
  DfFunctional_VWN a;
}

// test1
//  rhoa= 0.17E+01 rhob= 0.17E+01 sigmaaa= 0.81E-11 sigmaab= 0.81E-11 sigmabb= 0.81E-11
//  zk            = -0.278978177367E+00
//  vrhoa         = -0.907896301530E-01
//  vrhob         = -0.907896301530E-01
//  vsigmaaa      =  0.000000000000E+00
//  vsigmaab      =  0.000000000000E+00
//  vsigmabb      =  0.000000000000E+00

//  v2rhoa2       =  0.129443214985E-01
//  v2rhoab       = -0.182559901422E-01
//  v2rhob2       =  0.129443214985E-01

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
void DfFunctional_VWNTest::testPointwise1(){
  // input
  const double dRhoA = 0.17E+01;
  const double dRhoB = 0.17E+01;

  // expected value
  const double zk = -0.278978177367E+00;
  const double vRhoA = -0.907896301530E-01;
  const double vRhoB = -0.907896301530E-01;

  // execute test
  DfFunctional_VWN f;
  
  double dFunctionalValue = f.getFunctional(dRhoA, dRhoB);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(zk, dFunctionalValue, ED_CHECK(zk));

  double dRoundF_roundRhoA, dRoundF_roundRhoB;
  f.getDerivativeFunctional(dRhoA, dRhoB,
			    &dRoundF_roundRhoA, &dRoundF_roundRhoB);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(vRhoA, dRoundF_roundRhoA, ED_CHECK(vRhoA));
  CPPUNIT_ASSERT_DOUBLES_EQUAL(vRhoB, dRoundF_roundRhoB, ED_CHECK(vRhoB));
}

void DfFunctional_VWNTest::testPointwise1_RKS(){
  // input
  const double dRhoA = 0.17E+01;
  //const double dRhoB = 0.17E+01;

  // expected value
  const double zk = -0.278978177367E+00;
  const double vRhoA = -0.907896301530E-01;
//   const double vRhoB = -0.907896301530E-01;

  // execute test
  DfFunctional_VWN f;
  
  double dFunctionalValue = f.getFunctional(dRhoA);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(zk, dFunctionalValue, ED_CHECK(zk));

  double dRoundF_roundRhoA;
  f.getDerivativeFunctional(dRhoA, &dRoundF_roundRhoA);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(vRhoA, dRoundF_roundRhoA, ED_CHECK(vRhoA));
//   CPPUNIT_ASSERT_DOUBLES_EQUAL(vRhoB, dRoundF_roundRhoB, ED_CHECK(vRhoB));
}

// test2
//  rhoa= 0.17E+01 rhob= 0.17E+01 sigmaaa= 0.17E+01 sigmaab= 0.17E+01 sigmabb= 0.17E+01
//  zk            = -0.278978177367E+00
//  vrhoa         = -0.907896301530E-01
//  vrhob         = -0.907896301530E-01
//  vsigmaaa      =  0.000000000000E+00
//  vsigmaab      =  0.000000000000E+00
//  vsigmabb      =  0.000000000000E+00

//  v2rhoa2       =  0.129443214985E-01
//  v2rhoab       = -0.182559901422E-01
//  v2rhob2       =  0.129443214985E-01

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
void DfFunctional_VWNTest::testPointwise2(){
  // input
  const double dRhoA = 0.17E+01;
  const double dRhoB = 0.17E+01;

  // expected value
  const double zk = -0.278978177367E+00;
  const double vRhoA = -0.907896301530E-01;
  const double vRhoB = -0.907896301530E-01;

  // execute test
  DfFunctional_VWN f;
  
  double dFunctionalValue = f.getFunctional(dRhoA, dRhoB);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(zk, dFunctionalValue, ED_CHECK(zk));

  double dRoundF_roundRhoA, dRoundF_roundRhoB;
  f.getDerivativeFunctional(dRhoA, dRhoB,
				&dRoundF_roundRhoA, &dRoundF_roundRhoB);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(vRhoA, dRoundF_roundRhoA, ED_CHECK(vRhoA));
  CPPUNIT_ASSERT_DOUBLES_EQUAL(vRhoB, dRoundF_roundRhoB, ED_CHECK(vRhoB));
}

void DfFunctional_VWNTest::testPointwise2_RKS(){
  // input
  const double dRhoA = 0.17E+01;
  //const double dRhoB = 0.17E+01;

  // expected value
  const double zk = -0.278978177367E+00;
  const double vRhoA = -0.907896301530E-01;
//   const double vRhoB = -0.907896301530E-01;

  // execute test
  DfFunctional_VWN f;
  
  double dFunctionalValue = f.getFunctional(dRhoA);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(zk, dFunctionalValue, ED_CHECK(zk));

  double dRoundF_roundRhoA;
  f.getDerivativeFunctional(dRhoA, &dRoundF_roundRhoA);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(vRhoA, dRoundF_roundRhoA, ED_CHECK(vRhoA));
//   CPPUNIT_ASSERT_DOUBLES_EQUAL(vRhoB, dRoundF_roundRhoB, ED_CHECK(vRhoB));
}

// test3
//  rhoa= 0.15E+01 rhob= 0.15E+01 sigmaaa= 0.36E+02 sigmaab= 0.36E+02 sigmabb= 0.36E+02
//  zk            = -0.242883397986E+00
//  vrhoa         = -0.896613951966E-01
//  vrhob         = -0.896613951966E-01
//  vsigmaaa      =  0.000000000000E+00
//  vsigmaab      =  0.000000000000E+00
//  vsigmabb      =  0.000000000000E+00

//  v2rhoa2       =  0.144649744464E-01
//  v2rhoab       = -0.204638346911E-01
//  v2rhob2       =  0.144649744464E-01

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
void DfFunctional_VWNTest::testPointwise3(){
  // input
  const double dRhoA = 0.15E+01;
  const double dRhoB = 0.15E+01;

  // expected value
  const double zk = -0.242883397986E+00;
  const double vRhoA = -0.896613951966E-01;
  const double vRhoB = -0.896613951966E-01;

  // execute test
  DfFunctional_VWN f;
  
  double dFunctionalValue = f.getFunctional(dRhoA, dRhoB);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(zk, dFunctionalValue, ED_CHECK(zk));

  double dRoundF_roundRhoA, dRoundF_roundRhoB;
  f.getDerivativeFunctional(dRhoA, dRhoB,
			    &dRoundF_roundRhoA, &dRoundF_roundRhoB);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(vRhoA, dRoundF_roundRhoA, ED_CHECK(vRhoA));
  CPPUNIT_ASSERT_DOUBLES_EQUAL(vRhoB, dRoundF_roundRhoB, ED_CHECK(vRhoB));
}

// test4
//  rhoa= 0.88E-01 rhob= 0.88E-01 sigmaaa= 0.87E-01 sigmaab= 0.87E-01 sigmabb= 0.87E-01
//  zk            = -0.101483720780E-01
//  vrhoa         = -0.653289336535E-01
//  vrhob         = -0.653289336535E-01
//  vsigmaaa      =  0.000000000000E+00
//  vsigmaab      =  0.000000000000E+00
//  vsigmabb      =  0.000000000000E+00

//  v2rhoa2       =  0.171179303519E+00
//  v2rhoab       = -0.263237442473E+00
//  v2rhob2       =  0.171179303519E+00

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
void DfFunctional_VWNTest::testPointwise4(){
  // input
  const double dRhoA = 0.88E-01;
  const double dRhoB = 0.88E-01;

  // expected value
  const double zk = -0.101483720780E-01;
  const double vRhoA = -0.653289336535E-01;
  const double vRhoB = -0.653289336535E-01;

  // execute test
  DfFunctional_VWN f;
  
  double dFunctionalValue = f.getFunctional(dRhoA, dRhoB);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(zk, dFunctionalValue, ED_CHECK(zk));

  double dRoundF_roundRhoA, dRoundF_roundRhoB;
  f.getDerivativeFunctional(dRhoA, dRhoB,
			    &dRoundF_roundRhoA, &dRoundF_roundRhoB);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(vRhoA, dRoundF_roundRhoA, ED_CHECK(vRhoA));
  CPPUNIT_ASSERT_DOUBLES_EQUAL(vRhoB, dRoundF_roundRhoB, ED_CHECK(vRhoB));
}

// test5
//  rhoa= 0.18E+04 rhob= 0.18E+04 sigmaaa= 0.55E+00 sigmaab= 0.55E+00 sigmabb= 0.55E+00
//  zk            = -0.532741477023E+03
//  vrhoa         = -0.157944671704E+00
//  vrhob         = -0.157944671704E+00
//  vsigmaaa      =  0.000000000000E+00
//  vsigmaab      =  0.000000000000E+00
//  vsigmabb      =  0.000000000000E+00

//  v2rhoa2       =  0.223669588268E-04
//  v2rhoab       = -0.279504750065E-04
//  v2rhob2       =  0.223669588268E-04

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
void DfFunctional_VWNTest::testPointwise5(){
  // input
  const double dRhoA = 0.18E+04;
  const double dRhoB = 0.18E+04;

  // expected value
  const double zk = -0.532741477023E+03;
  const double vRhoA = -0.157944671704E+00;
  const double vRhoB = -0.157944671704E+00;

  // execute test
  DfFunctional_VWN f;
  
  double dFunctionalValue = f.getFunctional(dRhoA, dRhoB);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(zk, dFunctionalValue, ED_CHECK(zk));

  double dRoundF_roundRhoA, dRoundF_roundRhoB;
  f.getDerivativeFunctional(dRhoA, dRhoB,
			    &dRoundF_roundRhoA, &dRoundF_roundRhoB);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(vRhoA, dRoundF_roundRhoA, ED_CHECK(vRhoA));
  CPPUNIT_ASSERT_DOUBLES_EQUAL(vRhoB, dRoundF_roundRhoB, ED_CHECK(vRhoB));
}

// test6
//  rhoa= 0.18E+04 rhob= 0.18E+04 sigmaaa= 0.86E+04 sigmaab= 0.86E+04 sigmabb= 0.86E+04
//  zk            = -0.532741477023E+03
//  vrhoa         = -0.157944671704E+00
//  vrhob         = -0.157944671704E+00
//  vsigmaaa      =  0.000000000000E+00
//  vsigmaab      =  0.000000000000E+00
//  vsigmabb      =  0.000000000000E+00

//  v2rhoa2       =  0.223669588268E-04
//  v2rhoab       = -0.279504750065E-04
//  v2rhob2       =  0.223669588268E-04

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
void DfFunctional_VWNTest::testPointwise6(){
  // input
  const double dRhoA = 0.18E+04;
  const double dRhoB = 0.18E+04;

  // expected value
  const double zk = -0.532741477023E+03;
  const double vRhoA = -0.157944671704E+00;
  const double vRhoB = -0.157944671704E+00;

  // execute test
  DfFunctional_VWN f;
  
  double dFunctionalValue = f.getFunctional(dRhoA, dRhoB);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(zk, dFunctionalValue, ED_CHECK(zk));

  double dRoundF_roundRhoA, dRoundF_roundRhoB;
  f.getDerivativeFunctional(dRhoA, dRhoB,
			    &dRoundF_roundRhoA, &dRoundF_roundRhoB);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(vRhoA, dRoundF_roundRhoA, ED_CHECK(vRhoA));
  CPPUNIT_ASSERT_DOUBLES_EQUAL(vRhoB, dRoundF_roundRhoB, ED_CHECK(vRhoB));
}

// test7
//  rhoa= 0.16E+04 rhob= 0.16E+04 sigmaaa= 0.37E+10 sigmaab= 0.37E+10 sigmabb= 0.37E+10
//  zk            = -0.469795648279E+03
//  vrhoa         = -0.156761418492E+00
//  vrhob         = -0.156761418492E+00
//  vsigmaaa      =  0.000000000000E+00
//  vsigmaab      =  0.000000000000E+00
//  vsigmabb      =  0.000000000000E+00

//  v2rhoa2       =  0.249624818738E-04
//  v2rhoab       = -0.312385564000E-04
//  v2rhob2       =  0.249624818738E-04

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
void DfFunctional_VWNTest::testPointwise7(){
  // input
  const double dRhoA = 0.16E+04;
  const double dRhoB = 0.16E+04;

  // expected value
  const double zk = -0.469795648279E+03;
  const double vRhoA = -0.156761418492E+00;
  const double vRhoB = -0.156761418492E+00;

  // execute test
  DfFunctional_VWN f;
  
  double dFunctionalValue = f.getFunctional(dRhoA, dRhoB);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(zk, dFunctionalValue, ED_CHECK(zk));

  double dRoundF_roundRhoA, dRoundF_roundRhoB;
  f.getDerivativeFunctional(dRhoA, dRhoB,
			    &dRoundF_roundRhoA, &dRoundF_roundRhoB);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(vRhoA, dRoundF_roundRhoA, ED_CHECK(vRhoA));
  CPPUNIT_ASSERT_DOUBLES_EQUAL(vRhoB, dRoundF_roundRhoB, ED_CHECK(vRhoB));
}

//test8
//  rhoa= 0.26E+00 rhob= 0.26E+00 sigmaaa= 0.28E+00 sigmaab= 0.28E+00 sigmabb= 0.28E+00
//  zk            = -0.344300981310E-01
//  vrhoa         = -0.743196778205E-01
//  vrhob         = -0.743196778205E-01
//  vsigmaaa      =  0.000000000000E+00
//  vsigmaab      =  0.000000000000E+00
//  vsigmabb      =  0.000000000000E+00

//  v2rhoa2       =  0.673517043262E-01
//  v2rhoab       = -0.999958453856E-01
//  v2rhob2       =  0.673517043262E-01

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
void DfFunctional_VWNTest::testPointwise8(){
  // input
  const double dRhoA = 0.26E+00;
  const double dRhoB = 0.26E+00;

  // expected value
  const double zk = -0.344300981310E-01;
  const double vRhoA = -0.743196778205E-01;
  const double vRhoB = -0.743196778205E-01;

  // execute test
  DfFunctional_VWN f;
  
  double dFunctionalValue = f.getFunctional(dRhoA, dRhoB);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(zk, dFunctionalValue, ED_CHECK(zk));

  double dRoundF_roundRhoA, dRoundF_roundRhoB;
  f.getDerivativeFunctional(dRhoA, dRhoB,
			    &dRoundF_roundRhoA, &dRoundF_roundRhoB);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(vRhoA, dRoundF_roundRhoA, ED_CHECK(vRhoA));
  CPPUNIT_ASSERT_DOUBLES_EQUAL(vRhoB, dRoundF_roundRhoB, ED_CHECK(vRhoB));
}

//test9
//  rhoa= 0.53E+05 rhob= 0.53E+05 sigmaaa= 0.96E+05 sigmaab= 0.96E+05 sigmabb= 0.96E+05
//  zk            = -0.193016440590E+05
//  vrhoa         = -0.192271893744E+00
//  vrhob         = -0.192271893744E+00
//  vsigmaaa      =  0.000000000000E+00
//  vsigmaab      =  0.000000000000E+00
//  vsigmabb      =  0.000000000000E+00

//  v2rhoa2       =  0.934975367062E-06
//  v2rhoab       = -0.112791439004E-05
//  v2rhob2       =  0.934975367062E-06

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
void DfFunctional_VWNTest::testPointwise9(){
  // input
  const double dRhoA = 0.53E+05;
  const double dRhoB = 0.53E+05;

  // expected value
  const double zk = -0.193016440590E+05;
  const double vRhoA = -0.192271893744E+00;
  const double vRhoB = -0.192271893744E+00;

  // execute test
  DfFunctional_VWN f;
  
  double dFunctionalValue = f.getFunctional(dRhoA, dRhoB);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(zk, dFunctionalValue, ED_CHECK(zk));

  double dRoundF_roundRhoA, dRoundF_roundRhoB;
  f.getDerivativeFunctional(dRhoA, dRhoB,
			    &dRoundF_roundRhoA, &dRoundF_roundRhoB);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(vRhoA, dRoundF_roundRhoA, ED_CHECK(vRhoA));
  CPPUNIT_ASSERT_DOUBLES_EQUAL(vRhoB, dRoundF_roundRhoB, ED_CHECK(vRhoB));
}

//test10
//  rhoa= 0.47E+05 rhob= 0.47E+05 sigmaaa= 0.29E+14 sigmaab= 0.29E+14 sigmabb= 0.29E+14
//  zk            = -0.170016041806E+05
//  vrhoa         = -0.191043581788E+00
//  vrhob         = -0.191043581788E+00
//  vsigmaaa      =  0.000000000000E+00
//  vsigmaab      =  0.000000000000E+00
//  vsigmabb      =  0.000000000000E+00

//  v2rhoa2       =  0.104726097228E-05
//  v2rhoab       = -0.126473945017E-05
//  v2rhob2       =  0.104726097228E-05

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
void DfFunctional_VWNTest::testPointwise10(){
  // input
  const double dRhoA = 0.47E+05;
  const double dRhoB = 0.47E+05;

  // expected value
  const double zk = -0.170016041806E+05;
  const double vRhoA = -0.191043581788E+00;
  const double vRhoB = -0.191043581788E+00;

  // execute test
  DfFunctional_VWN f;
  
  double dFunctionalValue = f.getFunctional(dRhoA, dRhoB);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(zk, dFunctionalValue, ED_CHECK(zk));

  double dRoundF_roundRhoA, dRoundF_roundRhoB;
  f.getDerivativeFunctional(dRhoA, dRhoB,
			    &dRoundF_roundRhoA, &dRoundF_roundRhoB);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(vRhoA, dRoundF_roundRhoA, ED_CHECK(vRhoA));
  CPPUNIT_ASSERT_DOUBLES_EQUAL(vRhoB, dRoundF_roundRhoB, ED_CHECK(vRhoB));
}

//test11
//  rhoa= 0.15E+00 rhob= 0.15E+00 sigmaaa= 0.16E+00 sigmaab= 0.16E+00 sigmabb= 0.16E+00
//  zk            = -0.185432270230E-01
//  vrhoa         = -0.697024933328E-01
//  vrhob         = -0.697024933328E-01
//  vsigmaaa      =  0.000000000000E+00
//  vsigmaab      =  0.000000000000E+00
//  vsigmabb      =  0.000000000000E+00

//  v2rhoa2       =  0.108356325631E+00
//  v2rhoab       = -0.163679368941E+00
//  v2rhob2       =  0.108356325631E+00

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
void DfFunctional_VWNTest::testPointwise11(){
  // input
  const double dRhoA = 0.15E+00;
  const double dRhoB = 0.15E+00;

  // expected value
  const double zk = -0.185432270230E-01;
  const double vRhoA = -0.697024933328E-01;
  const double vRhoB = -0.697024933328E-01;

  // execute test
  DfFunctional_VWN f;
  
  double dFunctionalValue = f.getFunctional(dRhoA, dRhoB);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(zk, dFunctionalValue, ED_CHECK(zk));

  double dRoundF_roundRhoA, dRoundF_roundRhoB;
  f.getDerivativeFunctional(dRhoA, dRhoB,
			    &dRoundF_roundRhoA, &dRoundF_roundRhoB);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(vRhoA, dRoundF_roundRhoA, ED_CHECK(vRhoA));
  CPPUNIT_ASSERT_DOUBLES_EQUAL(vRhoB, dRoundF_roundRhoB, ED_CHECK(vRhoB));
}

// test12
//  rhoa= 0.35E+01 rhob= 0.00E+00 sigmaaa= 0.46E-10 sigmaab= 0.00E+00 sigmabb= 0.00E+00
//  zk            = -0.149673920800E+00
//  vrhoa         = -0.471789714951E-01
//  vrhob         =  0.000000000000E+00
//  vsigmaaa      =  0.000000000000E+00
//  vsigmaab      =  0.000000000000E+00
//  vsigmabb      =  0.000000000000E+00

//  v2rhoa2       = -0.130468710291E-02
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
void DfFunctional_VWNTest::testPointwise12(){
  // input
  const double dRhoA = 0.35E+01;
  const double dRhoB = 0.00E+00;

  // expected value
  const double zk = -0.149673920800E+00;
  const double vRhoA = -0.471789714951E-01;
  const double vRhoB = 0.000000000000E+00;

  // execute test
  DfFunctional_VWN f;
  
  double dFunctionalValue = f.getFunctional(dRhoA, dRhoB);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(zk, dFunctionalValue, ED_CHECK(zk));

  double dRoundF_roundRhoA, dRoundF_roundRhoB;
  f.getDerivativeFunctional(dRhoA, dRhoB,
			    &dRoundF_roundRhoA, &dRoundF_roundRhoB);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(vRhoA, dRoundF_roundRhoA, ED_CHECK(vRhoA));
  CPPUNIT_ASSERT_DOUBLES_EQUAL(vRhoB, dRoundF_roundRhoB, ED_CHECK(vRhoB));
}

//test13
//  rhoa= 0.35E+01 rhob= 0.00E+00 sigmaaa= 0.34E+01 sigmaab= 0.00E+00 sigmabb= 0.00E+00
//  zk            = -0.149673920800E+00
//  vrhoa         = -0.471789714951E-01
//  vrhob         =  0.000000000000E+00
//  vsigmaaa      =  0.000000000000E+00
//  vsigmaab      =  0.000000000000E+00
//  vsigmabb      =  0.000000000000E+00

//  v2rhoa2       = -0.130468710291E-02
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
void DfFunctional_VWNTest::testPointwise13(){
  // input
  const double dRhoA = 0.35E+01;
  const double dRhoB = 0.00E+00;

  // expected value
  const double zk = -0.149673920800E+00;
  const double vRhoA = -0.471789714951E-01;
  const double vRhoB = 0.000000000000E+00;

  // execute test
  DfFunctional_VWN f;
  
  double dFunctionalValue = f.getFunctional(dRhoA, dRhoB);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(zk, dFunctionalValue, ED_CHECK(zk));

  double dRoundF_roundRhoA, dRoundF_roundRhoB;
  f.getDerivativeFunctional(dRhoA, dRhoB,
			    &dRoundF_roundRhoA, &dRoundF_roundRhoB);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(vRhoA, dRoundF_roundRhoA, ED_CHECK(vRhoA));
  CPPUNIT_ASSERT_DOUBLES_EQUAL(vRhoB, dRoundF_roundRhoB, ED_CHECK(vRhoB));
}

// test14
//  rhoa= 0.30E+01 rhob= 0.00E+00 sigmaaa= 0.20E+03 sigmaab= 0.00E+00 sigmabb= 0.00E+00
//  zk            = -0.126255646553E+00
//  vrhoa         = -0.464766057364E-01
//  vrhob         =  0.000000000000E+00
//  vsigmaaa      =  0.000000000000E+00
//  vsigmaab      =  0.000000000000E+00
//  vsigmabb      =  0.000000000000E+00

//  v2rhoa2       = -0.151540949651E-02
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
void DfFunctional_VWNTest::testPointwise14(){
  // input
  const double dRhoA = 0.30E+01;
  const double dRhoB = 0.00E+00;

  // expected value
  const double zk = -0.126255646553E+00;
  const double vRhoA = -0.464766057364E-01;
  const double vRhoB = 0.000000000000E+00;

  // execute test
  DfFunctional_VWN f;
  
  double dFunctionalValue = f.getFunctional(dRhoA, dRhoB);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(zk, dFunctionalValue, ED_CHECK(zk));

  double dRoundF_roundRhoA, dRoundF_roundRhoB;
  f.getDerivativeFunctional(dRhoA, dRhoB,
				&dRoundF_roundRhoA, &dRoundF_roundRhoB);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(vRhoA, dRoundF_roundRhoA, ED_CHECK(vRhoA));
  CPPUNIT_ASSERT_DOUBLES_EQUAL(vRhoB, dRoundF_roundRhoB, ED_CHECK(vRhoB));
}

//test15
//  rhoa= 0.58E-01 rhob= 0.00E+00 sigmaaa= 0.47E-01 sigmaab= 0.00E+00 sigmabb= 0.00E+00
//  zk            = -0.151940115037E-02
//  vrhoa         = -0.297993118612E-01
//  vrhob         =  0.000000000000E+00
//  vsigmaaa      =  0.000000000000E+00
//  vsigmaab      =  0.000000000000E+00
//  vsigmabb      =  0.000000000000E+00

//  v2rhoa2       = -0.663179036789E-01
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
void DfFunctional_VWNTest::testPointwise15(){
  // input
  const double dRhoA = 0.58E-01;
  const double dRhoB = 0.00E+00;

  // expected value
  const double zk = -0.151940115037E-02;
  const double vRhoA = -0.297993118612E-01;
  const double vRhoB = 0.000000000000E+00;

  // execute test
  DfFunctional_VWN f;
  
  double dFunctionalValue = f.getFunctional(dRhoA, dRhoB);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(zk, dFunctionalValue, ED_CHECK(zk));

  double dRoundF_roundRhoA, dRoundF_roundRhoB;
  f.getDerivativeFunctional(dRhoA, dRhoB,
			    &dRoundF_roundRhoA, &dRoundF_roundRhoB);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(vRhoA, dRoundF_roundRhoA, ED_CHECK(vRhoA));
  CPPUNIT_ASSERT_DOUBLES_EQUAL(vRhoB, dRoundF_roundRhoB, ED_CHECK(vRhoB));
}

// test16
//  rhoa= 0.82E+02 rhob= 0.81E+02 sigmaaa= 0.49E+07 sigmaab= 0.49E+07 sigmabb= 0.49E+07
//  zk            = -0.191815088538E+02
//  vrhoa         = -0.126816376070E+00
//  vrhob         = -0.127719732599E+00
//  vsigmaaa      =  0.000000000000E+00
//  vsigmaab      =  0.000000000000E+00
//  vsigmabb      =  0.000000000000E+00

//  v2rhoa2       =  0.386783435451E-03
//  v2rhoab       = -0.511451913589E-03
//  v2rhob2       =  0.397050991054E-03

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
void DfFunctional_VWNTest::testPointwise16(){
  // input
  const double dRhoA = 0.82E+02;
  const double dRhoB = 0.81E+02;

  // expected value
  const double zk = -0.191815088538E+02;
  const double vRhoA = -0.126816376070E+00;
  const double vRhoB = -0.127719732599E+00;

  // execute test
  DfFunctional_VWN f;
  
  double dFunctionalValue = f.getFunctional(dRhoA, dRhoB);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(zk, dFunctionalValue, ED_CHECK(zk));

  double dRoundF_roundRhoA, dRoundF_roundRhoB;
  f.getDerivativeFunctional(dRhoA, dRhoB,
			    &dRoundF_roundRhoA, &dRoundF_roundRhoB);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(vRhoA, dRoundF_roundRhoA, ED_CHECK(vRhoA));
  CPPUNIT_ASSERT_DOUBLES_EQUAL(vRhoB, dRoundF_roundRhoB, ED_CHECK(vRhoB));
}

// test17
//  rhoa= 0.39E+02 rhob= 0.38E+02 sigmaaa= 0.81E+06 sigmaab= 0.82E+06 sigmabb= 0.82E+06
//  zk            = -0.851077910672E+01
//  vrhoa         = -0.119099058995E+00
//  vrhob         = -0.120906044904E+00
//  vsigmaaa      =  0.000000000000E+00
//  vsigmaab      =  0.000000000000E+00
//  vsigmabb      =  0.000000000000E+00

//  v2rhoa2       =  0.756836181702E-03
//  v2rhoab       = -0.102861281830E-02
//  v2rhob2       =  0.800136175083E-03

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
void DfFunctional_VWNTest::testPointwise17(){
  // input
  const double dRhoA = 0.39E+02;
  const double dRhoB = 0.38E+02;

  // expected value
  const double zk = -0.851077910672E+01;
  const double vRhoA = -0.119099058995E+00;
  const double vRhoB = -0.120906044904E+00;

  // execute test
  DfFunctional_VWN f;
  
  double dFunctionalValue = f.getFunctional(dRhoA, dRhoB);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(zk, dFunctionalValue, ED_CHECK(zk));

  double dRoundF_roundRhoA, dRoundF_roundRhoB;
  f.getDerivativeFunctional(dRhoA, dRhoB,
			    &dRoundF_roundRhoA, &dRoundF_roundRhoB);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(vRhoA, dRoundF_roundRhoA, ED_CHECK(vRhoA));
  CPPUNIT_ASSERT_DOUBLES_EQUAL(vRhoB, dRoundF_roundRhoB, ED_CHECK(vRhoB));
}

// test 18
//  rhoa= 0.13E+00 rhob= 0.95E-01 sigmaaa= 0.15E+00 sigmaab= 0.18E+00 sigmabb= 0.22E+00
//  zk            = -0.132928938310E-01
//  vrhoa         = -0.615913447654E-01
//  vrhob         = -0.739125478228E-01
//  vsigmaaa      =  0.000000000000E+00
//  vsigmaab      =  0.000000000000E+00
//  vsigmabb      =  0.000000000000E+00

//  v2rhoa2       =  0.963811005993E-01
//  v2rhoab       = -0.210790984367E+00
//  v2rhob2       =  0.193655360428E+00

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
void DfFunctional_VWNTest::testPointwise18(){
  // input
  const double dRhoA = 0.13E+00;
  const double dRhoB = 0.95E-01;

  // expected value
  const double zk = -0.132928938310E-01;
  const double vRhoA = -0.615913447654E-01;
  const double vRhoB = -0.739125478228E-01;

  // execute test
  DfFunctional_VWN f;
  
  double dFunctionalValue = f.getFunctional(dRhoA, dRhoB);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(zk, dFunctionalValue, ED_CHECK(zk));

  double dRoundF_roundRhoA, dRoundF_roundRhoB;
  f.getDerivativeFunctional(dRhoA, dRhoB,
			    &dRoundF_roundRhoA, &dRoundF_roundRhoB);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(vRhoA, dRoundF_roundRhoA, ED_CHECK(vRhoA));
  CPPUNIT_ASSERT_DOUBLES_EQUAL(vRhoB, dRoundF_roundRhoB, ED_CHECK(vRhoB));
}

// test19
//  rhoa= 0.78E-01 rhob= 0.31E-01 sigmaaa= 0.41E-02 sigmaab= 0.38E-02 sigmabb= 0.36E-02
//  zk            = -0.551636081757E-02
//  vrhoa         = -0.482590171567E-01
//  vrhob         = -0.811494055217E-01
//  vsigmaaa      =  0.000000000000E+00
//  vsigmaab      =  0.000000000000E+00
//  vsigmabb      =  0.000000000000E+00

//  v2rhoa2       =  0.864399166736E-01
//  v2rhoab       = -0.417632287561E+00
//  v2rhob2       =  0.710228207528E+00

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
void DfFunctional_VWNTest::testPointwise19(){
  // input
  const double dRhoA = 0.78E-01;
  const double dRhoB = 0.31E-01;

  // expected value
  const double zk = -0.551636081757E-02;
  const double vRhoA = -0.482590171567E-01;
  const double vRhoB = -0.811494055217E-01;

  // execute test
  DfFunctional_VWN f;
  
  double dFunctionalValue = f.getFunctional(dRhoA, dRhoB);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(zk, dFunctionalValue, ED_CHECK(zk));

  double dRoundF_roundRhoA, dRoundF_roundRhoB;
  f.getDerivativeFunctional(dRhoA, dRhoB,
			    &dRoundF_roundRhoA, &dRoundF_roundRhoB);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(vRhoA, dRoundF_roundRhoA, ED_CHECK(vRhoA));
  CPPUNIT_ASSERT_DOUBLES_EQUAL(vRhoB, dRoundF_roundRhoB, ED_CHECK(vRhoB));
}

// test20
//  rhoa= 0.50E+02 rhob= 0.49E+02 sigmaaa= 0.11E+06 sigmaab= 0.11E+06 sigmabb= 0.11E+06
//  zk            = -0.111786110384E+02
//  vrhoa         = -0.121711427226E+00
//  vrhob         = -0.123144247973E+00
//  vsigmaaa      =  0.000000000000E+00
//  vsigmaab      =  0.000000000000E+00
//  vsigmabb      =  0.000000000000E+00

//  v2rhoa2       =  0.605386382334E-03
//  v2rhoab       = -0.814117737252E-03
//  v2rhob2       =  0.632128050135E-03

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
void DfFunctional_VWNTest::testPointwise20(){
  // input
  const double dRhoA = 0.50E+02;
  const double dRhoB = 0.49E+02;

  // expected value
  const double zk = -0.111786110384E+02;
  const double vRhoA = -0.121711427226E+00;
  const double vRhoB = -0.123144247973E+00;

  // execute test
  DfFunctional_VWN f;
  
  double dFunctionalValue = f.getFunctional(dRhoA, dRhoB);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(zk, dFunctionalValue, ED_CHECK(zk));

  double dRoundF_roundRhoA, dRoundF_roundRhoB;
  f.getDerivativeFunctional(dRhoA, dRhoB,
			    &dRoundF_roundRhoA, &dRoundF_roundRhoB);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(vRhoA, dRoundF_roundRhoA, ED_CHECK(vRhoA));
  CPPUNIT_ASSERT_DOUBLES_EQUAL(vRhoB, dRoundF_roundRhoB, ED_CHECK(vRhoB));
}

//  rhoa= 0.40E+02 rhob= 0.40E+02 sigmaaa= 0.99E+05 sigmaab= 0.98E+05 sigmabb= 0.98E+05
//  zk            = -0.887177858837E+01
//  vrhoa         = -0.120365709995E+00
//  vrhob         = -0.120365709995E+00
//  vsigmaaa      =  0.000000000000E+00
//  vsigmaab      =  0.000000000000E+00
//  vsigmabb      =  0.000000000000E+00

//  v2rhoa2       =  0.751587360250E-03
//  v2rhoab       = -0.992735036010E-03
//  v2rhob2       =  0.751587360250E-03

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
void DfFunctional_VWNTest::testPointwise21(){
  // input
  const double dRhoA = 0.40E+02;
  const double dRhoB = 0.40E+02;

  // expected value
  const double zk = -0.887177858837E+01;
  const double vRhoA = -0.120365709995E+00;
  const double vRhoB = -0.120365709995E+00;

  // execute test
  DfFunctional_VWN f;
  
  double dFunctionalValue = f.getFunctional(dRhoA, dRhoB);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(zk, dFunctionalValue, ED_CHECK(zk));

  double dRoundF_roundRhoA, dRoundF_roundRhoB;
  f.getDerivativeFunctional(dRhoA, dRhoB,
			    &dRoundF_roundRhoA, &dRoundF_roundRhoB);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(vRhoA, dRoundF_roundRhoA, ED_CHECK(vRhoA));
  CPPUNIT_ASSERT_DOUBLES_EQUAL(vRhoB, dRoundF_roundRhoB, ED_CHECK(vRhoB));
}

//  rhoa= 0.12E+00 rhob= 0.10E+00 sigmaaa= 0.12E+00 sigmaab= 0.13E+00 sigmabb= 0.14E+00
//  zk            = -0.130284832590E-01
//  vrhoa         = -0.637104594188E-01
//  vrhob         = -0.708674029459E-01
//  vsigmaaa      =  0.000000000000E+00
//  vsigmaab      =  0.000000000000E+00
//  vsigmabb      =  0.000000000000E+00

//  v2rhoa2       =  0.114871866539E+00
//  v2rhoab       = -0.215499991137E+00
//  v2rhob2       =  0.172160006882E+00

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
void DfFunctional_VWNTest::testPointwise22(){
  // input
  const double dRhoA = 0.12E+00;
  const double dRhoB = 0.10E+00;

  // expected value
  const double zk = -0.130284832590E-01;
  const double vRhoA = -0.637104594188E-01;
  const double vRhoB = -0.708674029459E-01;

  // execute test
  DfFunctional_VWN f;
  
  double dFunctionalValue = f.getFunctional(dRhoA, dRhoB);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(zk, dFunctionalValue, ED_CHECK(zk));

  double dRoundF_roundRhoA, dRoundF_roundRhoB;
  f.getDerivativeFunctional(dRhoA, dRhoB,
			    &dRoundF_roundRhoA, &dRoundF_roundRhoB);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(vRhoA, dRoundF_roundRhoA, ED_CHECK(vRhoA));
  CPPUNIT_ASSERT_DOUBLES_EQUAL(vRhoB, dRoundF_roundRhoB, ED_CHECK(vRhoB));
}

//  rhoa= 0.48E-01 rhob= 0.25E-01 sigmaaa= 0.46E-02 sigmaab= 0.44E-02 sigmabb= 0.41E-02
//  zk            = -0.360439923610E-02
//  vrhoa         = -0.488659901421E-01
//  vrhob         = -0.708954619404E-01
//  vsigmaaa      =  0.000000000000E+00
//  vsigmaab      =  0.000000000000E+00
//  vsigmabb      =  0.000000000000E+00

//  v2rhoa2       =  0.164372114971E+00
//  v2rhoab       = -0.574640089589E+00
//  v2rhob2       =  0.724506769540E+00

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
void DfFunctional_VWNTest::testPointwise23(){
  // input
  const double dRhoA = 0.48E-01;
  const double dRhoB = 0.25E-01;

  // expected value
  const double zk = -0.360439923610E-02;
  const double vRhoA = -0.488659901421E-01;
  const double vRhoB = -0.708954619404E-01;

  // execute test
  DfFunctional_VWN f;
  
  double dFunctionalValue = f.getFunctional(dRhoA, dRhoB);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(zk, dFunctionalValue, ED_CHECK(zk));

  double dRoundF_roundRhoA, dRoundF_roundRhoB;
  f.getDerivativeFunctional(dRhoA, dRhoB,
			    &dRoundF_roundRhoA, &dRoundF_roundRhoB);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(vRhoA, dRoundF_roundRhoA, ED_CHECK(vRhoA));
  CPPUNIT_ASSERT_DOUBLES_EQUAL(vRhoB, dRoundF_roundRhoB, ED_CHECK(vRhoB));
}


