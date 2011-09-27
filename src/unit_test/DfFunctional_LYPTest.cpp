#include <vector>
#include <string>
#include "digit.h"
#include "DfFunctional_LYPTest.h"

#define EFFCTIVE_DIGITS 12
#define ED_CHECK(x) (std::pow(10.0, -EFFCTIVE_DIGITS + digit(x)))
const double DfFunctional_LYPTest::EPS = 1.0E-10;

// CPPUNIT_ASSERT( condition );
// condition$B$,56(B(false,0)$B$G$"$C$?$H$-!"<:GT$7$^$9!#(B
//
// CPPUNIT_ASSERT_MESSAGE( message, condition );
// condition$B$,56$G$"$C$?$H$-!"<:GT$7$^$9!#$3$N$H$-(Bmessage$B$r=PNO$7$^$9!#(B
//
// CPPUNIT_FAIL( message );
// $BI,$:<:GT$7$^$9!#(Bmessage$B$r=PNO$7$^$9!#(B
//
// CPPUNIT_ASSERT_EQUAL( expected, actual );
// $BF@$i$l$?7k2L(Bactual$B$,4|BT$9$kCM(Bexpected$B$G$J$+$C$?$H$-!"$9$J$o$A(Bexpected != actual$B$N$H$-$K<:GT$7$^$9!#(B
//
// CPPUNIT_ASSERT_EQUAL_MESSAGE( message, expected, actual );
// $BF@$i$l$?7k2L(Bactual$B$,4|BT$9$kCM(Bexpected$B$G$J$+$C$?$H$-!"$9$J$o$A(Bexpected != actual$B$N$H$-$K<:GT$7$^$9!#$3$N$H$-(Bmessage$B$r=PNO$7$^$9!#(B
//
// CPPUNIT_ASSERT_DOUBLES_EQUAL( expected, actual, delta );
// $BF@$i$l$?7k2L(Bactual$B$H4|BT$9$kCM(Bexpected$B$H$N:9$,(Bdelta$B$h$jBg$-$$$H$-!"<:GT$7$^$9!#(B

// =====================================================================
// data used from: http://www.cse.scitech.ac.uk/ccg/dft/data_pt_x_lda.html
//

void DfFunctional_LYPTest::testConstructer(){
  DfFunctional_LYP a;
}

// test1
//  rhoa= 0.17E+01 rhob= 0.17E+01 sigmaaa= 0.81E-11 sigmaab= 0.81E-11 sigmabb= 0.81E-11
//  zk            = -0.179175399535E+00
//  vrhoa         = -0.567254370239E-01
//  vrhob         = -0.567254370239E-01
//  vsigmaaa      =  0.603063052247E-04
//  vsigmaab      =  0.562668577012E-04
//  vsigmabb      =  0.603063052247E-04

//  v2rhoa2       =  0.133393844920E-01
//  v2rhoab       = -0.152396936370E-01
//  v2rhob2       =  0.133393844920E-01

//  v2rhoasigmaaa = -0.139567528486E-03
//  v2rhoasigmaab = -0.291609267169E-04
//  v2rhoasigmabb =  0.811155700026E-04
//  v2rhobsigmaaa =  0.811155700026E-04
//  v2rhobsigmaab = -0.291609267169E-04
//  v2rhobsigmabb = -0.139567528486E-03

//  v2sigmaaa2    =  0.000000000000E+00
//  v2sigmaaaab   =  0.000000000000E+00
//  v2sigmaaabb   =  0.000000000000E+00
//  v2sigmaab2    =  0.000000000000E+00
//  v2sigmaabbb   =  0.000000000000E+00
//  v2sigmabb2    =  0.000000000000E+00
void DfFunctional_LYPTest::testPointwise1(){
  // input
  const double dRhoA = 0.17E+01;
  const double dRhoB = 0.17E+01;
  const double dGammaAA = 0.81E-11;
  const double dGammaAB = 0.81E-11;
  const double dGammaBB = 0.81E-11;

  // expected value
  const double zk = -0.179175399535E+00;
  const double vRhoA = -0.567254370239E-01;
  const double vRhoB = -0.567254370239E-01;
  const double vGammaAA = 0.603063052247E-04;
  const double vGammaAB = 0.562668577012E-04;
  const double vGammaBB = 0.603063052247E-04;

  // execute test
  DfFunctional_LYP f;
  
  double dFunctionalValue = f.getFunctional(dRhoA, dRhoB, dGammaAA, dGammaAB, dGammaBB);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(zk, dFunctionalValue, ED_CHECK(zk));

  double dRoundF_roundRhoA, dRoundF_roundRhoB;
  double dRoundF_roundGammaAA, dRoundF_roundGammaAB, dRoundF_roundGammaBB;
  f.getDerivativeFunctional(dRhoA, dRhoB, dGammaAA, dGammaAB, dGammaBB,
			    &dRoundF_roundRhoA, &dRoundF_roundRhoB,
			    &dRoundF_roundGammaAA, &dRoundF_roundGammaAB, &dRoundF_roundGammaBB);

  CPPUNIT_ASSERT_DOUBLES_EQUAL(vRhoA, dRoundF_roundRhoA, ED_CHECK(vRhoA));
  CPPUNIT_ASSERT_DOUBLES_EQUAL(vRhoB, dRoundF_roundRhoB, ED_CHECK(vRhoB));
  CPPUNIT_ASSERT_DOUBLES_EQUAL(vGammaAA, dRoundF_roundGammaAA, ED_CHECK(vGammaAA));
  CPPUNIT_ASSERT_DOUBLES_EQUAL(vGammaAB, dRoundF_roundGammaAB, ED_CHECK(vGammaAB));
  CPPUNIT_ASSERT_DOUBLES_EQUAL(vGammaBB, dRoundF_roundGammaBB, ED_CHECK(vGammaBB));
}

void DfFunctional_LYPTest::testPointwise1_RKS(){
  // input
  const double dRhoA = 0.17E+01;
  //const double dRhoB = 0.17E+01;
  const double dGammaAA = 0.81E-11;
  //const double dGammaAB = 0.81E-11;
  //const double dGammaBB = 0.81E-11;

  // expected value
  const double zk = -0.179175399535E+00;
  const double vRhoA = -0.567254370239E-01;
  //const double vRhoB = -0.567254370239E-01;
  const double vGammaAA = 0.603063052247E-04;
  const double vGammaAB = 0.562668577012E-04;
  //const double vGammaBB = 0.603063052247E-04;

  // execute test
  DfFunctional_LYP f;
  
  double dFunctionalValue = f.getFunctional(dRhoA, dGammaAA);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(zk, dFunctionalValue, ED_CHECK(zk));

  double dRoundF_roundRhoA;
  double dRoundF_roundGammaAA, dRoundF_roundGammaAB;
  f.getDerivativeFunctional(dRhoA, dGammaAA,
			    &dRoundF_roundRhoA, &dRoundF_roundGammaAA, &dRoundF_roundGammaAB);

  CPPUNIT_ASSERT_DOUBLES_EQUAL(vRhoA, dRoundF_roundRhoA, ED_CHECK(vRhoA));
  //CPPUNIT_ASSERT_DOUBLES_EQUAL(vRhoB, dRoundF_roundRhoB, ED_CHECK(vRhoB));
  CPPUNIT_ASSERT_DOUBLES_EQUAL(vGammaAA, dRoundF_roundGammaAA, ED_CHECK(vGammaAA));
  CPPUNIT_ASSERT_DOUBLES_EQUAL(vGammaAB, dRoundF_roundGammaAB, ED_CHECK(vGammaAB));
  //CPPUNIT_ASSERT_DOUBLES_EQUAL(vGammaBB, dRoundF_roundGammaBB, ED_CHECK(vGammaBB));
}

// test2
//  rhoa= 0.17E+01 rhob= 0.17E+01 sigmaaa= 0.17E+01 sigmaab= 0.17E+01 sigmabb= 0.17E+01
//  zk            = -0.178874704439E+00
//  vrhoa         = -0.568743789287E-01
//  vrhob         = -0.568743789287E-01
//  vsigmaaa      =  0.603063052247E-04
//  vsigmaab      =  0.562668577012E-04
//  vsigmabb      =  0.603063052247E-04

//  v2rhoa2       =  0.138019299001E-01
//  v2rhoab       = -0.154679607816E-01
//  v2rhob2       =  0.138019299001E-01

//  v2rhoasigmaaa = -0.139567528486E-03
//  v2rhoasigmaab = -0.291609267169E-04
//  v2rhoasigmabb =  0.811155700026E-04
//  v2rhobsigmaaa =  0.811155700026E-04
//  v2rhobsigmaab = -0.291609267169E-04
//  v2rhobsigmabb = -0.139567528486E-03

//  v2sigmaaa2    =  0.000000000000E+00
//  v2sigmaaaab   =  0.000000000000E+00
//  v2sigmaaabb   =  0.000000000000E+00
//  v2sigmaab2    =  0.000000000000E+00
//  v2sigmaabbb   =  0.000000000000E+00
//  v2sigmabb2    =  0.000000000000E+00
void DfFunctional_LYPTest::testPointwise2(){
  // input
  const double dRhoA = 0.17E+01;
  const double dRhoB = 0.17E+01;
  const double dGammaAA = 0.17E+01;
  const double dGammaAB = 0.17E+01;
  const double dGammaBB = 0.17E+01;

  // expected value
  const double zk = -0.178874704439E+00;
  const double vRhoA = -0.568743789287E-01;
  const double vRhoB = -0.568743789287E-01;
  const double vGammaAA = 0.603063052247E-04;
  const double vGammaAB = 0.562668577012E-04;
  const double vGammaBB = 0.603063052247E-04;

  // execute test
  DfFunctional_LYP f;
  
  double dFunctionalValue = f.getFunctional(dRhoA, dRhoB, dGammaAA, dGammaAB, dGammaBB);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(zk, dFunctionalValue, ED_CHECK(zk));

  double dRoundF_roundRhoA, dRoundF_roundRhoB;
  double dRoundF_roundGammaAA, dRoundF_roundGammaAB, dRoundF_roundGammaBB;
  f.getDerivativeFunctional(dRhoA, dRhoB, dGammaAA, dGammaAB, dGammaBB,
			    &dRoundF_roundRhoA, &dRoundF_roundRhoB,
			    &dRoundF_roundGammaAA, &dRoundF_roundGammaAB, &dRoundF_roundGammaBB);

  CPPUNIT_ASSERT_DOUBLES_EQUAL(vRhoA, dRoundF_roundRhoA, ED_CHECK(vRhoA));
  CPPUNIT_ASSERT_DOUBLES_EQUAL(vRhoB, dRoundF_roundRhoB, ED_CHECK(vRhoB));
  CPPUNIT_ASSERT_DOUBLES_EQUAL(vGammaAA, dRoundF_roundGammaAA, ED_CHECK(vGammaAA));
  CPPUNIT_ASSERT_DOUBLES_EQUAL(vGammaAB, dRoundF_roundGammaAB, ED_CHECK(vGammaAB));
  CPPUNIT_ASSERT_DOUBLES_EQUAL(vGammaBB, dRoundF_roundGammaBB, ED_CHECK(vGammaBB));
}

void DfFunctional_LYPTest::testPointwise2_RKS(){
  // input
  const double dRhoA = 0.17E+01;
  //const double dRhoB = 0.17E+01;
  const double dGammaAA = 0.17E+01;
  //const double dGammaAB = 0.17E+01;
  //const double dGammaBB = 0.17E+01;

  // expected value
  const double zk = -0.178874704439E+00;
  const double vRhoA = -0.568743789287E-01;
  //const double vRhoB = -0.568743789287E-01;
  const double vGammaAA = 0.603063052247E-04;
  const double vGammaAB = 0.562668577012E-04;
  //const double vGammaBB = 0.603063052247E-04;

  // execute test
  DfFunctional_LYP f;
  
  double dFunctionalValue = f.getFunctional(dRhoA, dGammaAA);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(zk, dFunctionalValue, ED_CHECK(zk));

  double dRoundF_roundRhoA;
  double dRoundF_roundGammaAA, dRoundF_roundGammaAB;
  f.getDerivativeFunctional(dRhoA, dGammaAA,
			    &dRoundF_roundRhoA, &dRoundF_roundGammaAA, &dRoundF_roundGammaAB);

  CPPUNIT_ASSERT_DOUBLES_EQUAL(vRhoA, dRoundF_roundRhoA, ED_CHECK(vRhoA));
  //CPPUNIT_ASSERT_DOUBLES_EQUAL(vRhoB, dRoundF_roundRhoB, ED_CHECK(vRhoB));
  CPPUNIT_ASSERT_DOUBLES_EQUAL(vGammaAA, dRoundF_roundGammaAA, ED_CHECK(vGammaAA));
  CPPUNIT_ASSERT_DOUBLES_EQUAL(vGammaAB, dRoundF_roundGammaAB, ED_CHECK(vGammaAB));
  //CPPUNIT_ASSERT_DOUBLES_EQUAL(vGammaBB, dRoundF_roundGammaBB, ED_CHECK(vGammaBB));
}

// test3
//  rhoa= 0.15E+01 rhob= 0.15E+01 sigmaaa= 0.36E+02 sigmaab= 0.36E+02 sigmabb= 0.36E+02
//  zk            = -0.148704448499E+00
//  vrhoa         = -0.607218516385E-01
//  vrhob         = -0.607218516385E-01
//  vsigmaaa      =  0.741080708899E-04
//  vsigmaab      =  0.701383136250E-04
//  vsigmabb      =  0.741080708899E-04

//  v2rhoa2       =  0.303235212020E-01
//  v2rhoab       = -0.246964871909E-01
//  v2rhob2       =  0.303235212020E-01

//  v2rhoasigmaaa = -0.192341863313E-03
//  v2rhoasigmaab = -0.411266458866E-04
//  v2rhoasigmabb =  0.111052541885E-03
//  v2rhobsigmaaa =  0.111052541885E-03
//  v2rhobsigmaab = -0.411266458866E-04
//  v2rhobsigmabb = -0.192341863313E-03

//  v2sigmaaa2    =  0.000000000000E+00
//  v2sigmaaaab   =  0.000000000000E+00
//  v2sigmaaabb   =  0.000000000000E+00
//  v2sigmaab2    =  0.000000000000E+00
//  v2sigmaabbb   =  0.000000000000E+00
//  v2sigmabb2    =  0.000000000000E+00
void DfFunctional_LYPTest::testPointwise3(){
  // input
  const double dRhoA = 0.15E+01;
  const double dRhoB = 0.15E+01;
  const double dGammaAA = 0.36E+02;
  const double dGammaAB = 0.36E+02;
  const double dGammaBB = 0.36E+02;

  // expected value
  const double zk = -0.148704448499E+00;
  const double vRhoA = -0.607218516385E-01;
  const double vRhoB = -0.607218516385E-01;
  const double vGammaAA = 0.741080708899E-04;
  const double vGammaAB = 0.701383136250E-04;
  const double vGammaBB = 0.741080708899E-04;

  // execute test
  DfFunctional_LYP f;
  
  double dFunctionalValue = f.getFunctional(dRhoA, dRhoB, dGammaAA, dGammaAB, dGammaBB);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(zk, dFunctionalValue, ED_CHECK(zk));

  double dRoundF_roundRhoA, dRoundF_roundRhoB;
  double dRoundF_roundGammaAA, dRoundF_roundGammaAB, dRoundF_roundGammaBB;
  f.getDerivativeFunctional(dRhoA, dRhoB, dGammaAA, dGammaAB, dGammaBB,
			    &dRoundF_roundRhoA, &dRoundF_roundRhoB,
			    &dRoundF_roundGammaAA, &dRoundF_roundGammaAB, &dRoundF_roundGammaBB);

  CPPUNIT_ASSERT_DOUBLES_EQUAL(vRhoA, dRoundF_roundRhoA, ED_CHECK(vRhoA));
  CPPUNIT_ASSERT_DOUBLES_EQUAL(vRhoB, dRoundF_roundRhoB, ED_CHECK(vRhoB));
  CPPUNIT_ASSERT_DOUBLES_EQUAL(vGammaAA, dRoundF_roundGammaAA, ED_CHECK(vGammaAA));
  CPPUNIT_ASSERT_DOUBLES_EQUAL(vGammaAB, dRoundF_roundGammaAB, ED_CHECK(vGammaAB));
  CPPUNIT_ASSERT_DOUBLES_EQUAL(vGammaBB, dRoundF_roundGammaBB, ED_CHECK(vGammaBB));
}

void DfFunctional_LYPTest::testPointwise3_RKS(){
  // input
  const double dRhoA = 0.15E+01;
  //const double dRhoB = 0.15E+01;
  const double dGammaAA = 0.36E+02;
  //const double dGammaAB = 0.36E+02;
  //const double dGammaBB = 0.36E+02;

  // expected value
  const double zk = -0.148704448499E+00;
  const double vRhoA = -0.607218516385E-01;
  //const double vRhoB = -0.607218516385E-01;
  const double vGammaAA = 0.741080708899E-04;
  const double vGammaAB = 0.701383136250E-04;
  //const double vGammaBB = 0.741080708899E-04;

  // execute test
  DfFunctional_LYP f;
  
  double dFunctionalValue = f.getFunctional(dRhoA, dGammaAA);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(zk, dFunctionalValue, ED_CHECK(zk));

  double dRoundF_roundRhoA;
  double dRoundF_roundGammaAA, dRoundF_roundGammaAB;
  f.getDerivativeFunctional(dRhoA, dGammaAA,
			    &dRoundF_roundRhoA, &dRoundF_roundGammaAA, &dRoundF_roundGammaAB);

  CPPUNIT_ASSERT_DOUBLES_EQUAL(vRhoA, dRoundF_roundRhoA, ED_CHECK(vRhoA));
  //CPPUNIT_ASSERT_DOUBLES_EQUAL(vRhoB, dRoundF_roundRhoB, ED_CHECK(vRhoB));
  CPPUNIT_ASSERT_DOUBLES_EQUAL(vGammaAA, dRoundF_roundGammaAA, ED_CHECK(vGammaAA));
  CPPUNIT_ASSERT_DOUBLES_EQUAL(vGammaAB, dRoundF_roundGammaAB, ED_CHECK(vGammaAB));
  //CPPUNIT_ASSERT_DOUBLES_EQUAL(vGammaBB, dRoundF_roundGammaBB, ED_CHECK(vGammaBB));
}

// test4
//  rhoa= 0.88E-01 rhob= 0.88E-01 sigmaaa= 0.87E-01 sigmaab= 0.87E-01 sigmabb= 0.87E-01
//  zk            = -0.465024100803E-02
//  vrhoa         = -0.610987646298E-01
//  vrhob         = -0.610987646298E-01
//  vsigmaaa      =  0.694062751783E-02
//  vsigmaab      =  0.876388739167E-02
//  vsigmabb      =  0.694062751783E-02

//  v2rhoa2       =  0.939780887292E+00
//  v2rhoab       = -0.498868805878E+00
//  v2rhob2       =  0.939780887292E+00

//  v2rhoasigmaaa = -0.228361811495E+00
//  v2rhoasigmaab = -0.807965355076E-01
//  v2rhoasigmabb =  0.107205512562E+00
//  v2rhobsigmaaa =  0.107205512562E+00
//  v2rhobsigmaab = -0.807965355076E-01
//  v2rhobsigmabb = -0.228361811495E+00

//  v2sigmaaa2    =  0.000000000000E+00
//  v2sigmaaaab   =  0.000000000000E+00
//  v2sigmaaabb   =  0.000000000000E+00
//  v2sigmaab2    =  0.000000000000E+00
//  v2sigmaabbb   =  0.000000000000E+00
//  v2sigmabb2    =  0.000000000000E+00
void DfFunctional_LYPTest::testPointwise4(){
  // input
  const double dRhoA = 0.88E-01;
  const double dRhoB = 0.88E-01;
  const double dGammaAA = 0.87E-01;
  const double dGammaAB = 0.87E-01;
  const double dGammaBB = 0.87E-01;

  // expected value
  const double zk = -0.465024100803E-02;
  const double vRhoA = -0.610987646298E-01;
  const double vRhoB = -0.610987646298E-01;
  const double vGammaAA = 0.694062751783E-02;
  const double vGammaAB = 0.876388739167E-02;
  const double vGammaBB = 0.694062751783E-02;

  // execute test
  DfFunctional_LYP f;
  
  double dFunctionalValue = f.getFunctional(dRhoA, dRhoB, dGammaAA, dGammaAB, dGammaBB);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(zk, dFunctionalValue, ED_CHECK(zk));

  double dRoundF_roundRhoA, dRoundF_roundRhoB;
  double dRoundF_roundGammaAA, dRoundF_roundGammaAB, dRoundF_roundGammaBB;
  f.getDerivativeFunctional(dRhoA, dRhoB, dGammaAA, dGammaAB, dGammaBB,
			    &dRoundF_roundRhoA, &dRoundF_roundRhoB,
			    &dRoundF_roundGammaAA, &dRoundF_roundGammaAB, &dRoundF_roundGammaBB);

  CPPUNIT_ASSERT_DOUBLES_EQUAL(vRhoA, dRoundF_roundRhoA, ED_CHECK(vRhoA));
  CPPUNIT_ASSERT_DOUBLES_EQUAL(vRhoB, dRoundF_roundRhoB, ED_CHECK(vRhoB));
  CPPUNIT_ASSERT_DOUBLES_EQUAL(vGammaAA, dRoundF_roundGammaAA, ED_CHECK(vGammaAA));
  CPPUNIT_ASSERT_DOUBLES_EQUAL(vGammaAB, dRoundF_roundGammaAB, ED_CHECK(vGammaAB));
  CPPUNIT_ASSERT_DOUBLES_EQUAL(vGammaBB, dRoundF_roundGammaBB, ED_CHECK(vGammaBB));
}

void DfFunctional_LYPTest::testPointwise4_RKS(){
  // input
  const double dRhoA = 0.88E-01;
  //const double dRhoB = 0.88E-01;
  const double dGammaAA = 0.87E-01;
  //const double dGammaAB = 0.87E-01;
  //const double dGammaBB = 0.87E-01;

  // expected value
  const double zk = -0.465024100803E-02;
  const double vRhoA = -0.610987646298E-01;
  //const double vRhoB = -0.610987646298E-01;
  const double vGammaAA = 0.694062751783E-02;
  const double vGammaAB = 0.876388739167E-02;
  //const double vGammaBB = 0.694062751783E-02;

  // execute test
  DfFunctional_LYP f;
  
  double dFunctionalValue = f.getFunctional(dRhoA, dGammaAA);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(zk, dFunctionalValue, ED_CHECK(zk));

  double dRoundF_roundRhoA;
  double dRoundF_roundGammaAA, dRoundF_roundGammaAB;
  f.getDerivativeFunctional(dRhoA, dGammaAA,
			    &dRoundF_roundRhoA, &dRoundF_roundGammaAA, &dRoundF_roundGammaAB);

  CPPUNIT_ASSERT_DOUBLES_EQUAL(vRhoA, dRoundF_roundRhoA, ED_CHECK(vRhoA));
  //CPPUNIT_ASSERT_DOUBLES_EQUAL(vRhoB, dRoundF_roundRhoB, ED_CHECK(vRhoB));
  CPPUNIT_ASSERT_DOUBLES_EQUAL(vGammaAA, dRoundF_roundGammaAA, ED_CHECK(vGammaAA));
  CPPUNIT_ASSERT_DOUBLES_EQUAL(vGammaAB, dRoundF_roundGammaAB, ED_CHECK(vGammaAB));
  //CPPUNIT_ASSERT_DOUBLES_EQUAL(vGammaBB, dRoundF_roundGammaBB, ED_CHECK(vGammaBB));
}

//test5
//  rhoa= 0.18E+04 rhob= 0.18E+04 sigmaaa= 0.55E+00 sigmaab= 0.55E+00 sigmabb= 0.55E+00
//  zk            = -0.237638430952E+03
//  vrhoa         = -0.665993271139E-01
//  vrhob         = -0.665993271139E-01
//  vsigmaaa      =  0.540555769454E-09
//  vsigmaab      =  0.260773591409E-09
//  vsigmabb      =  0.540555769454E-09

//  v2rhoa2       =  0.144305866044E-04
//  v2rhoab       = -0.146537621392E-04
//  v2rhob2       =  0.144305866044E-04

//  v2rhoasigmaaa = -0.158850198364E-11
//  v2rhoasigmaab = -0.124882367057E-12
//  v2rhoasigmabb =  0.108678004324E-11
//  v2rhobsigmaaa =  0.108678004324E-11
//  v2rhobsigmaab = -0.124882367057E-12
//  v2rhobsigmabb = -0.158850198364E-11

//  v2sigmaaa2    =  0.000000000000E+00
//  v2sigmaaaab   =  0.000000000000E+00
//  v2sigmaaabb   =  0.000000000000E+00
//  v2sigmaab2    =  0.000000000000E+00
//  v2sigmaabbb   =  0.000000000000E+00
//  v2sigmabb2    =  0.000000000000E+00
void DfFunctional_LYPTest::testPointwise5(){
  // input
  const double dRhoA = 0.18E+04;
  const double dRhoB = 0.18E+04;
  const double dGammaAA = 0.55E+00;
  const double dGammaAB = 0.55E+00;
  const double dGammaBB = 0.55E+00;

  // expected value
  const double zk = -0.237638430952E+03;
  const double vRhoA = -0.665993271139E-01;
  const double vRhoB = -0.665993271139E-01;
  const double vGammaAA = 0.540555769454E-09;
  const double vGammaAB = 0.260773591409E-09;
  const double vGammaBB = 0.540555769454E-09;

  // execute test
  DfFunctional_LYP f;
  
  double dFunctionalValue = f.getFunctional(dRhoA, dRhoB, dGammaAA, dGammaAB, dGammaBB);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(zk, dFunctionalValue, ED_CHECK(zk));

  double dRoundF_roundRhoA, dRoundF_roundRhoB;
  double dRoundF_roundGammaAA, dRoundF_roundGammaAB, dRoundF_roundGammaBB;
  f.getDerivativeFunctional(dRhoA, dRhoB, dGammaAA, dGammaAB, dGammaBB,
			    &dRoundF_roundRhoA, &dRoundF_roundRhoB,
			    &dRoundF_roundGammaAA, &dRoundF_roundGammaAB, &dRoundF_roundGammaBB);

  CPPUNIT_ASSERT_DOUBLES_EQUAL(vRhoA, dRoundF_roundRhoA, ED_CHECK(vRhoA));
  CPPUNIT_ASSERT_DOUBLES_EQUAL(vRhoB, dRoundF_roundRhoB, ED_CHECK(vRhoB));
  CPPUNIT_ASSERT_DOUBLES_EQUAL(vGammaAA, dRoundF_roundGammaAA, ED_CHECK(vGammaAA));
  CPPUNIT_ASSERT_DOUBLES_EQUAL(vGammaAB, dRoundF_roundGammaAB, ED_CHECK(vGammaAB));
  CPPUNIT_ASSERT_DOUBLES_EQUAL(vGammaBB, dRoundF_roundGammaBB, ED_CHECK(vGammaBB));
}


void DfFunctional_LYPTest::testPointwise5_RKS(){
  // input
  const double dRhoA = 0.18E+04;
  //const double dRhoB = 0.18E+04;
  const double dGammaAA = 0.55E+00;
  //const double dGammaAB = 0.55E+00;
  //const double dGammaBB = 0.55E+00;

  // expected value
  const double zk = -0.237638430952E+03;
  const double vRhoA = -0.665993271139E-01;
  //const double vRhoB = -0.665993271139E-01;
  const double vGammaAA = 0.540555769454E-09;
  const double vGammaAB = 0.260773591409E-09;
  //const double vGammaBB = 0.540555769454E-09;

  // execute test
  DfFunctional_LYP f;
  
  double dFunctionalValue = f.getFunctional(dRhoA, dGammaAA);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(zk, dFunctionalValue, ED_CHECK(zk));

  double dRoundF_roundRhoA;
  double dRoundF_roundGammaAA, dRoundF_roundGammaAB;
  f.getDerivativeFunctional(dRhoA, dGammaAA,
			    &dRoundF_roundRhoA, &dRoundF_roundGammaAA, &dRoundF_roundGammaAB);

  CPPUNIT_ASSERT_DOUBLES_EQUAL(vRhoA, dRoundF_roundRhoA, ED_CHECK(vRhoA));
  //CPPUNIT_ASSERT_DOUBLES_EQUAL(vRhoB, dRoundF_roundRhoB, ED_CHECK(vRhoB));
  CPPUNIT_ASSERT_DOUBLES_EQUAL(vGammaAA, dRoundF_roundGammaAA, ED_CHECK(vGammaAA));
  CPPUNIT_ASSERT_DOUBLES_EQUAL(vGammaAB, dRoundF_roundGammaAB, ED_CHECK(vGammaAB));
  //CPPUNIT_ASSERT_DOUBLES_EQUAL(vGammaBB, dRoundF_roundGammaBB, ED_CHECK(vGammaBB));
}

//test6
//  rhoa= 0.18E+04 rhob= 0.18E+04 sigmaaa= 0.86E+04 sigmaab= 0.86E+04 sigmabb= 0.86E+04
//  zk            = -0.237638419413E+03
//  vrhoa         = -0.665993325023E-01
//  vrhob         = -0.665993325023E-01
//  vsigmaaa      =  0.540555769454E-09
//  vsigmaab      =  0.260773591409E-09
//  vsigmabb      =  0.540555769454E-09

//  v2rhoa2       =  0.144306117017E-04
//  v2rhoab       = -0.146537792037E-04
//  v2rhob2       =  0.144306117017E-04

//  v2rhoasigmaaa = -0.158850198364E-11
//  v2rhoasigmaab = -0.124882367057E-12
//  v2rhoasigmabb =  0.108678004324E-11
//  v2rhobsigmaaa =  0.108678004324E-11
//  v2rhobsigmaab = -0.124882367057E-12
//  v2rhobsigmabb = -0.158850198364E-11

//  v2sigmaaa2    =  0.000000000000E+00
//  v2sigmaaaab   =  0.000000000000E+00
//  v2sigmaaabb   =  0.000000000000E+00
//  v2sigmaab2    =  0.000000000000E+00
//  v2sigmaabbb   =  0.000000000000E+00
//  v2sigmabb2    =  0.000000000000E+00

//test7
//  rhoa= 0.16E+04 rhob= 0.16E+04 sigmaaa= 0.37E+10 sigmaab= 0.37E+10 sigmabb= 0.37E+10
//  zk            = -0.204955819468E+03
//  vrhoa         = -0.697313972428E-01
//  vrhob         = -0.697313972428E-01
//  vsigmaaa      =  0.658117255270E-09
//  vsigmaab      =  0.319515846955E-09
//  vsigmabb      =  0.658117255270E-09

//  v2rhoa2       =  0.328352356048E-04
//  v2rhoab       = -0.277614126197E-04
//  v2rhob2       =  0.328352356048E-04

//  v2rhoasigmaaa = -0.217196289305E-11
//  v2rhoasigmaab = -0.172301824814E-12
//  v2rhoasigmabb =  0.148472344173E-11
//  v2rhobsigmaaa =  0.148472344173E-11
//  v2rhobsigmaab = -0.172301824814E-12
//  v2rhobsigmabb = -0.217196289305E-11

//  v2sigmaaa2    =  0.000000000000E+00
//  v2sigmaaaab   =  0.000000000000E+00
//  v2sigmaaabb   =  0.000000000000E+00
//  v2sigmaab2    =  0.000000000000E+00
//  v2sigmaabbb   =  0.000000000000E+00
//  v2sigmabb2    =  0.000000000000E+00

//test8
//  rhoa= 0.26E+00 rhob= 0.26E+00 sigmaaa= 0.28E+00 sigmaab= 0.28E+00 sigmabb= 0.28E+00
//  zk            = -0.216471443930E-01
//  vrhoa         = -0.526995644881E-01
//  vrhob         = -0.526995644881E-01
//  vsigmaaa      =  0.127229458316E-02
//  vsigmaab      =  0.145290124256E-02
//  vsigmabb      =  0.127229458316E-02

//  v2rhoa2       =  0.134806861660E+00
//  v2rhoab       = -0.118221268828E+00
//  v2rhob2       =  0.134806861660E+00

//  v2rhoasigmaaa = -0.160681338920E-01
//  v2rhoasigmaab = -0.472571671370E-02
//  v2rhoasigmabb =  0.827564425842E-02
//  v2rhobsigmaaa =  0.827564425842E-02
//  v2rhobsigmaab = -0.472571671370E-02
//  v2rhobsigmabb = -0.160681338920E-01

//  v2sigmaaa2    =  0.000000000000E+00
//  v2sigmaaaab   =  0.000000000000E+00
//  v2sigmaaabb   =  0.000000000000E+00
//  v2sigmaab2    =  0.000000000000E+00
//  v2sigmaabbb   =  0.000000000000E+00
//  v2sigmabb2    =  0.000000000000E+00

//test9
//  rhoa= 0.53E+05 rhob= 0.53E+05 sigmaaa= 0.96E+05 sigmaab= 0.96E+05 sigmabb= 0.96E+05
//  zk            = -0.712575807209E+04
//  vrhoa         = -0.674210020682E-01
//  vrhob         = -0.674210020682E-01
//  vsigmaaa      =  0.190814948140E-11
//  vsigmaab      =  0.816481365419E-12
//  vsigmabb      =  0.190814948140E-11

//  v2rhoa2       =  0.495469678644E-06
//  v2rhoab       = -0.497965254694E-06
//  v2rhob2       =  0.495469678644E-06

//  v2rhoasigmaaa = -0.196249114596E-15
//  v2rhoasigmaab = -0.130135394800E-16
//  v2rhoasigmabb =  0.136188177073E-15
//  v2rhobsigmaaa =  0.136188177073E-15
//  v2rhobsigmaab = -0.130135394800E-16
//  v2rhobsigmabb = -0.196249114596E-15

//  v2sigmaaa2    =  0.000000000000E+00
//  v2sigmaaaab   =  0.000000000000E+00
//  v2sigmaaabb   =  0.000000000000E+00
//  v2sigmaab2    =  0.000000000000E+00
//  v2sigmaabbb   =  0.000000000000E+00
//  v2sigmabb2    =  0.000000000000E+00

//test10
//  rhoa= 0.47E+05 rhob= 0.47E+05 sigmaaa= 0.29E+14 sigmaab= 0.29E+14 sigmabb= 0.29E+14
//  zk            = -0.615255831041E+04
//  vrhoa         = -0.703265205307E-01
//  vrhob         = -0.703265205307E-01
//  vsigmaaa      =  0.233162530668E-11
//  vsigmaab      =  0.100028081622E-11
//  vsigmabb      =  0.233162530668E-11

//  v2rhoa2       =  0.110943737062E-05
//  v2rhoab       = -0.946187279564E-06
//  v2rhob2       =  0.110943737062E-05

//  v2rhoasigmaaa = -0.270252437366E-15
//  v2rhoasigmaab = -0.179871585341E-16
//  v2rhoasigmabb =  0.187490309515E-15
//  v2rhobsigmaaa =  0.187490309515E-15
//  v2rhobsigmaab = -0.179871585341E-16
//  v2rhobsigmabb = -0.270252437366E-15

//  v2sigmaaa2    =  0.000000000000E+00
//  v2sigmaaaab   =  0.000000000000E+00
//  v2sigmaaabb   =  0.000000000000E+00
//  v2sigmaab2    =  0.000000000000E+00
//  v2sigmaabbb   =  0.000000000000E+00
//  v2sigmabb2    =  0.000000000000E+00

//test11
//  rhoa= 0.15E+00 rhob= 0.15E+00 sigmaaa= 0.16E+00 sigmaab= 0.16E+00 sigmabb= 0.16E+00
//  zk            = -0.106605459023E-01
//  vrhoa         = -0.547299644084E-01
//  vrhob         = -0.547299644084E-01
//  vsigmaaa      =  0.303386607339E-02
//  vsigmaab      =  0.365209975632E-02
//  vsigmabb      =  0.303386607339E-02

//  v2rhoa2       =  0.346499910453E+00
//  v2rhoab       = -0.239377906462E+00
//  v2rhob2       =  0.346499910453E+00

//  v2rhoasigmaaa = -0.624221561613E-01
//  v2rhoasigmaab = -0.201987433700E-01
//  v2rhoasigmabb =  0.307385926277E-01
//  v2rhobsigmaaa =  0.307385926277E-01
//  v2rhobsigmaab = -0.201987433700E-01
//  v2rhobsigmabb = -0.624221561613E-01

//  v2sigmaaa2    =  0.000000000000E+00
//  v2sigmaaaab   =  0.000000000000E+00
//  v2sigmaaabb   =  0.000000000000E+00
//  v2sigmaab2    =  0.000000000000E+00
//  v2sigmaabbb   =  0.000000000000E+00
//  v2sigmabb2    =  0.000000000000E+00

//test12
//  rhoa= 0.35E+01 rhob= 0.00E+00 sigmaaa= 0.46E-10 sigmaab= 0.00E+00 sigmabb= 0.00E+00
//  zk            =  0.000000000000E+00
//  vrhoa         =  0.000000000000E+00
//  vrhob         =  0.000000000000E+00
//  vsigmaaa      =  0.000000000000E+00
//  vsigmaab      =  0.000000000000E+00
//  vsigmabb      =  0.000000000000E+00

//  v2rhoa2       =  0.000000000000E+00
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
void DfFunctional_LYPTest::testPointwise12(){
  // input
  const double dRhoA = 0.35E+01;
  const double dRhoB = 0.00E+00;
  const double dGammaAA = 0.46E-10;
  const double dGammaAB = 0.00E+00;
  const double dGammaBB = 0.00E+00;

  // expected value
  const double zk = 0.000000000000E+00;
  const double vRhoA = 0.000000000000E+00;
  const double vRhoB = 0.000000000000E+00;
  const double vGammaAA = 0.000000000000E+00;
  const double vGammaAB = 0.000000000000E+00;
  const double vGammaBB = 0.000000000000E+00;

  // execute test
  DfFunctional_LYP f;
  
  double dFunctionalValue = f.getFunctional(dRhoA, dRhoB, dGammaAA, dGammaAB, dGammaBB);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(zk, dFunctionalValue, ED_CHECK(zk));

  double dRoundF_roundRhoA, dRoundF_roundRhoB;
  double dRoundF_roundGammaAA, dRoundF_roundGammaAB, dRoundF_roundGammaBB;
  f.getDerivativeFunctional(dRhoA, dRhoB, dGammaAA, dGammaAB, dGammaBB,
			    &dRoundF_roundRhoA, &dRoundF_roundRhoB,
			    &dRoundF_roundGammaAA, &dRoundF_roundGammaAB, &dRoundF_roundGammaBB);

  CPPUNIT_ASSERT_DOUBLES_EQUAL(vRhoA, dRoundF_roundRhoA, ED_CHECK(vRhoA));
  CPPUNIT_ASSERT_DOUBLES_EQUAL(vRhoB, dRoundF_roundRhoB, ED_CHECK(vRhoB));
  CPPUNIT_ASSERT_DOUBLES_EQUAL(vGammaAA, dRoundF_roundGammaAA, ED_CHECK(vGammaAA));
  CPPUNIT_ASSERT_DOUBLES_EQUAL(vGammaAB, dRoundF_roundGammaAB, ED_CHECK(vGammaAB));
  CPPUNIT_ASSERT_DOUBLES_EQUAL(vGammaBB, dRoundF_roundGammaBB, ED_CHECK(vGammaBB));
}

//test13
//  rhoa= 0.35E+01 rhob= 0.00E+00 sigmaaa= 0.34E+01 sigmaab= 0.00E+00 sigmabb= 0.00E+00
//  zk            =  0.000000000000E+00
//  vrhoa         =  0.000000000000E+00
//  vrhob         =  0.000000000000E+00
//  vsigmaaa      =  0.000000000000E+00
//  vsigmaab      =  0.000000000000E+00
//  vsigmabb      =  0.000000000000E+00

//  v2rhoa2       =  0.000000000000E+00
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
void DfFunctional_LYPTest::testPointwise13(){
  // input
  const double dRhoA = 0.35E+01;
  const double dRhoB = 0.00E+00;
  const double dGammaAA = 0.34E+01;
  const double dGammaAB = 0.00E+00;
  const double dGammaBB = 0.00E+00;

  // expected value
  const double zk = 0.000000000000E+00;
  const double vRhoA = 0.000000000000E+00;
  const double vRhoB = 0.000000000000E+00;
  const double vGammaAA = 0.000000000000E+00;
  const double vGammaAB = 0.000000000000E+00;
  const double vGammaBB = 0.000000000000E+00;

  // execute test
  DfFunctional_LYP f;
  
  double dFunctionalValue = f.getFunctional(dRhoA, dRhoB, dGammaAA, dGammaAB, dGammaBB);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(zk, dFunctionalValue, ED_CHECK(zk));

  double dRoundF_roundRhoA, dRoundF_roundRhoB;
  double dRoundF_roundGammaAA, dRoundF_roundGammaAB, dRoundF_roundGammaBB;
  f.getDerivativeFunctional(dRhoA, dRhoB, dGammaAA, dGammaAB, dGammaBB,
			    &dRoundF_roundRhoA, &dRoundF_roundRhoB,
			    &dRoundF_roundGammaAA, &dRoundF_roundGammaAB, &dRoundF_roundGammaBB);

  CPPUNIT_ASSERT_DOUBLES_EQUAL(vRhoA, dRoundF_roundRhoA, ED_CHECK(vRhoA));
  CPPUNIT_ASSERT_DOUBLES_EQUAL(vRhoB, dRoundF_roundRhoB, ED_CHECK(vRhoB));
  CPPUNIT_ASSERT_DOUBLES_EQUAL(vGammaAA, dRoundF_roundGammaAA, ED_CHECK(vGammaAA));
  CPPUNIT_ASSERT_DOUBLES_EQUAL(vGammaAB, dRoundF_roundGammaAB, ED_CHECK(vGammaAB));
  CPPUNIT_ASSERT_DOUBLES_EQUAL(vGammaBB, dRoundF_roundGammaBB, ED_CHECK(vGammaBB));
}

//test14
//  rhoa= 0.30E+01 rhob= 0.00E+00 sigmaaa= 0.20E+03 sigmaab= 0.00E+00 sigmabb= 0.00E+00
//  zk            =  0.000000000000E+00
//  vrhoa         =  0.000000000000E+00
//  vrhob         =  0.000000000000E+00
//  vsigmaaa      =  0.000000000000E+00
//  vsigmaab      =  0.000000000000E+00
//  vsigmabb      =  0.000000000000E+00

//  v2rhoa2       =  0.000000000000E+00
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
void DfFunctional_LYPTest::testPointwise14(){
  // input
  const double dRhoA = 0.30E+01;
  const double dRhoB = 0.00E+00;
  const double dGammaAA = 0.20E+03;
  const double dGammaAB = 0.00E+00;
  const double dGammaBB = 0.00E+00;

  // expected value
  const double zk = 0.000000000000E+00;
  const double vRhoA = 0.000000000000E+00;
  const double vRhoB = 0.000000000000E+00;
  const double vGammaAA = 0.000000000000E+00;
  const double vGammaAB = 0.000000000000E+00;
  const double vGammaBB = 0.000000000000E+00;

  // execute test
  DfFunctional_LYP f;
  
  double dFunctionalValue = f.getFunctional(dRhoA, dRhoB, dGammaAA, dGammaAB, dGammaBB);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(zk, dFunctionalValue, ED_CHECK(zk));

  double dRoundF_roundRhoA, dRoundF_roundRhoB;
  double dRoundF_roundGammaAA, dRoundF_roundGammaAB, dRoundF_roundGammaBB;
  f.getDerivativeFunctional(dRhoA, dRhoB, dGammaAA, dGammaAB, dGammaBB,
			    &dRoundF_roundRhoA, &dRoundF_roundRhoB,
			    &dRoundF_roundGammaAA, &dRoundF_roundGammaAB, &dRoundF_roundGammaBB);

  CPPUNIT_ASSERT_DOUBLES_EQUAL(vRhoA, dRoundF_roundRhoA, ED_CHECK(vRhoA));
  CPPUNIT_ASSERT_DOUBLES_EQUAL(vRhoB, dRoundF_roundRhoB, ED_CHECK(vRhoB));
  CPPUNIT_ASSERT_DOUBLES_EQUAL(vGammaAA, dRoundF_roundGammaAA, ED_CHECK(vGammaAA));
  CPPUNIT_ASSERT_DOUBLES_EQUAL(vGammaAB, dRoundF_roundGammaAB, ED_CHECK(vGammaAB));
  CPPUNIT_ASSERT_DOUBLES_EQUAL(vGammaBB, dRoundF_roundGammaBB, ED_CHECK(vGammaBB));
}

//  rhoa= 0.58E-01 rhob= 0.00E+00 sigmaaa= 0.47E-01 sigmaab= 0.00E+00 sigmabb= 0.00E+00
//  zk            =  0.000000000000E+00
//  vrhoa         =  0.000000000000E+00
//  vrhob         =  0.000000000000E+00
//  vsigmaaa      =  0.000000000000E+00
//  vsigmaab      =  0.000000000000E+00
//  vsigmabb      =  0.000000000000E+00

//  v2rhoa2       =  0.000000000000E+00
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
void DfFunctional_LYPTest::testPointwise15(){
  // input
  const double dRhoA = 0.58E-01;
  const double dRhoB = 0.00E+00;
  const double dGammaAA = 0.47E-01;
  const double dGammaAB = 0.00E+00;
  const double dGammaBB = 0.00E+00;

  // expected value
  const double zk = 0.000000000000E+00;
  const double vRhoA = 0.000000000000E+00;
  const double vRhoB = 0.000000000000E+00;
  const double vGammaAA = 0.000000000000E+00;
  const double vGammaAB = 0.000000000000E+00;
  const double vGammaBB = 0.000000000000E+00;

  // execute test
  DfFunctional_LYP f;
  
  double dFunctionalValue = f.getFunctional(dRhoA, dRhoB, dGammaAA, dGammaAB, dGammaBB);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(zk, dFunctionalValue, ED_CHECK(zk));

  double dRoundF_roundRhoA, dRoundF_roundRhoB;
  double dRoundF_roundGammaAA, dRoundF_roundGammaAB, dRoundF_roundGammaBB;
  f.getDerivativeFunctional(dRhoA, dRhoB, dGammaAA, dGammaAB, dGammaBB,
			    &dRoundF_roundRhoA, &dRoundF_roundRhoB,
			    &dRoundF_roundGammaAA, &dRoundF_roundGammaAB, &dRoundF_roundGammaBB);

  CPPUNIT_ASSERT_DOUBLES_EQUAL(vRhoA, dRoundF_roundRhoA, ED_CHECK(vRhoA));
  CPPUNIT_ASSERT_DOUBLES_EQUAL(vRhoB, dRoundF_roundRhoB, ED_CHECK(vRhoB));
  CPPUNIT_ASSERT_DOUBLES_EQUAL(vGammaAA, dRoundF_roundGammaAA, ED_CHECK(vGammaAA));
  CPPUNIT_ASSERT_DOUBLES_EQUAL(vGammaAB, dRoundF_roundGammaAB, ED_CHECK(vGammaAB));
  CPPUNIT_ASSERT_DOUBLES_EQUAL(vGammaBB, dRoundF_roundGammaBB, ED_CHECK(vGammaBB));
}

//  rhoa= 0.82E+02 rhob= 0.81E+02 sigmaaa= 0.49E+07 sigmaab= 0.49E+07 sigmabb= 0.49E+07
//  zk            = -0.903966286601E+01
//  vrhoa         = -0.759264713358E-01
//  vrhob         = -0.784370481382E-01
//  vsigmaaa      =  0.907992031805E-07
//  vsigmaab      =  0.581118313258E-07
//  vsigmabb      =  0.100370432217E-06

//  v2rhoa2       =  0.141362773226E-02
//  v2rhoab       = -0.105147646057E-02
//  v2rhob2       =  0.150485268125E-02

//  v2rhoasigmaaa = -0.565167905921E-08
//  v2rhoasigmaab = -0.515378386013E-09
//  v2rhoasigmabb =  0.376492504615E-08
//  v2rhobsigmaaa =  0.384361243685E-08
//  v2rhobsigmaab = -0.747058939790E-09
//  v2rhobsigmabb = -0.588190678666E-08

//  v2sigmaaa2    =  0.000000000000E+00
//  v2sigmaaaab   =  0.000000000000E+00
//  v2sigmaaabb   =  0.000000000000E+00
//  v2sigmaab2    =  0.000000000000E+00
//  v2sigmaabbb   =  0.000000000000E+00
//  v2sigmabb2    =  0.000000000000E+00
void DfFunctional_LYPTest::testPointwise16(){
  // input
  const double dRhoA = 0.82E+02;
  const double dRhoB = 0.81E+02;
  const double dGammaAA = 0.49E+07;
  const double dGammaAB = 0.49E+07;
  const double dGammaBB = 0.49E+07;

  // expected value
  const double zk = -0.903966286601E+01;
  const double vRhoA = -0.759264713358E-01;
  const double vRhoB = -0.784370481382E-01;
  const double vGammaAA = 0.907992031805E-07;
  const double vGammaAB = 0.581118313258E-07;
  const double vGammaBB = 0.100370432217E-06;

  // execute test
  DfFunctional_LYP f;
  
  double dFunctionalValue = f.getFunctional(dRhoA, dRhoB, dGammaAA, dGammaAB, dGammaBB);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(zk, dFunctionalValue, ED_CHECK(zk));

  double dRoundF_roundRhoA, dRoundF_roundRhoB;
  double dRoundF_roundGammaAA, dRoundF_roundGammaAB, dRoundF_roundGammaBB;
  f.getDerivativeFunctional(dRhoA, dRhoB, dGammaAA, dGammaAB, dGammaBB,
			    &dRoundF_roundRhoA, &dRoundF_roundRhoB,
			    &dRoundF_roundGammaAA, &dRoundF_roundGammaAB, &dRoundF_roundGammaBB);

  CPPUNIT_ASSERT_DOUBLES_EQUAL(vRhoA, dRoundF_roundRhoA, ED_CHECK(vRhoA));
  CPPUNIT_ASSERT_DOUBLES_EQUAL(vRhoB, dRoundF_roundRhoB, ED_CHECK(vRhoB));
  CPPUNIT_ASSERT_DOUBLES_EQUAL(vGammaAA, dRoundF_roundGammaAA, ED_CHECK(vGammaAA));
  CPPUNIT_ASSERT_DOUBLES_EQUAL(vGammaAB, dRoundF_roundGammaAB, ED_CHECK(vGammaAB));
  CPPUNIT_ASSERT_DOUBLES_EQUAL(vGammaBB, dRoundF_roundGammaBB, ED_CHECK(vGammaBB));
}

//  rhoa= 0.39E+02 rhob= 0.38E+02 sigmaaa= 0.81E+06 sigmaab= 0.82E+06 sigmabb= 0.82E+06
//  zk            = -0.402158795173E+01
//  vrhoa         = -0.762734644914E-01
//  vrhob         = -0.830226435821E-01
//  vsigmaaa      =  0.301052145436E-06
//  vsigmaab      =  0.220298633297E-06
//  vsigmabb      =  0.369624286402E-06

//  v2rhoa2       =  0.331769729999E-02
//  v2rhoab       = -0.248438749270E-02
//  v2rhob2       =  0.384280348438E-02

//  v2rhoasigmaaa = -0.398359773843E-07
//  v2rhoasigmaab = -0.335415277613E-08
//  v2rhoasigmabb =  0.263970784129E-07
//  v2rhobsigmaaa =  0.275886078235E-07
//  v2rhobsigmaab = -0.685474898360E-08
//  v2rhobsigmabb = -0.433118929134E-07

//  v2sigmaaa2    =  0.000000000000E+00
//  v2sigmaaaab   =  0.000000000000E+00
//  v2sigmaaabb   =  0.000000000000E+00
//  v2sigmaab2    =  0.000000000000E+00
//  v2sigmaabbb   =  0.000000000000E+00
//  v2sigmabb2    =  0.000000000000E+00
void DfFunctional_LYPTest::testPointwise17(){
  // input
  const double dRhoA = 0.39E+02;
  const double dRhoB = 0.38E+02;
  const double dGammaAA = 0.81E+06;
  const double dGammaAB = 0.82E+06;
  const double dGammaBB = 0.82E+06;

  // expected value
  const double zk = -0.402158795173E+01;
  const double vRhoA = -0.762734644914E-01;
  const double vRhoB = -0.830226435821E-01;
  const double vGammaAA = 0.301052145436E-06;
  const double vGammaAB = 0.220298633297E-06;
  const double vGammaBB = 0.369624286402E-06;

  // execute test
  DfFunctional_LYP f;
  
  double dFunctionalValue = f.getFunctional(dRhoA, dRhoB, dGammaAA, dGammaAB, dGammaBB);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(zk, dFunctionalValue, ED_CHECK(zk));

  double dRoundF_roundRhoA, dRoundF_roundRhoB;
  double dRoundF_roundGammaAA, dRoundF_roundGammaAB, dRoundF_roundGammaBB;
  f.getDerivativeFunctional(dRhoA, dRhoB, dGammaAA, dGammaAB, dGammaBB,
			    &dRoundF_roundRhoA, &dRoundF_roundRhoB,
			    &dRoundF_roundGammaAA, &dRoundF_roundGammaAB, &dRoundF_roundGammaBB);

  CPPUNIT_ASSERT_DOUBLES_EQUAL(vRhoA, dRoundF_roundRhoA, ED_CHECK(vRhoA));
  CPPUNIT_ASSERT_DOUBLES_EQUAL(vRhoB, dRoundF_roundRhoB, ED_CHECK(vRhoB));
  CPPUNIT_ASSERT_DOUBLES_EQUAL(vGammaAA, dRoundF_roundGammaAA, ED_CHECK(vGammaAA));
  CPPUNIT_ASSERT_DOUBLES_EQUAL(vGammaAB, dRoundF_roundGammaAB, ED_CHECK(vGammaAB));
  CPPUNIT_ASSERT_DOUBLES_EQUAL(vGammaBB, dRoundF_roundGammaBB, ED_CHECK(vGammaBB));
}

//  rhoa= 0.13E+00 rhob= 0.95E-01 sigmaaa= 0.15E+00 sigmaab= 0.18E+00 sigmabb= 0.22E+00
//  zk            = -0.535812609335E-02
//  vrhoa         = -0.428946532738E-01
//  vrhob         = -0.100396349615E+00
//  vsigmaaa      =  0.178546727213E-02
//  vsigmaab      =  0.678262565697E-02
//  vsigmabb      =  0.827650024238E-02

//  v2rhoa2       =  0.309477814856E+00
//  v2rhoab       = -0.329371433841E+00
//  v2rhob2       =  0.162055387526E+01

//  v2rhoasigmaaa = -0.781343180939E-01
//  v2rhoasigmaab = -0.432850540037E-02
//  v2rhoasigmabb =  0.363514849034E-01
//  v2rhobsigmaaa =  0.732832606483E-01
//  v2rhobsigmaab = -0.108581862510E+00
//  v2rhobsigmabb = -0.179156926286E+00

//  v2sigmaaa2    =  0.000000000000E+00
//  v2sigmaaaab   =  0.000000000000E+00
//  v2sigmaaabb   =  0.000000000000E+00
//  v2sigmaab2    =  0.000000000000E+00
//  v2sigmaabbb   =  0.000000000000E+00
//  v2sigmabb2    =  0.000000000000E+00
void DfFunctional_LYPTest::testPointwise18(){
  // input
  const double dRhoA = 0.13E+00;
  const double dRhoB = 0.95E-01;
  const double dGammaAA = 0.15E+00;
  const double dGammaAB = 0.18E+00;
  const double dGammaBB = 0.22E+00;

  // expected value
  const double zk = -0.535812609335E-02;
  const double vRhoA = -0.428946532738E-01;
  const double vRhoB = -0.100396349615E+00;
  const double vGammaAA = 0.178546727213E-02;
  const double vGammaAB = 0.678262565697E-02;
  const double vGammaBB = 0.827650024238E-02;

  // execute test
  DfFunctional_LYP f;
  
  double dFunctionalValue = f.getFunctional(dRhoA, dRhoB, dGammaAA, dGammaAB, dGammaBB);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(zk, dFunctionalValue, ED_CHECK(zk));

  double dRoundF_roundRhoA, dRoundF_roundRhoB;
  double dRoundF_roundGammaAA, dRoundF_roundGammaAB, dRoundF_roundGammaBB;
  f.getDerivativeFunctional(dRhoA, dRhoB, dGammaAA, dGammaAB, dGammaBB,
			    &dRoundF_roundRhoA, &dRoundF_roundRhoB,
			    &dRoundF_roundGammaAA, &dRoundF_roundGammaAB, &dRoundF_roundGammaBB);

  CPPUNIT_ASSERT_DOUBLES_EQUAL(vRhoA, dRoundF_roundRhoA, ED_CHECK(vRhoA));
  CPPUNIT_ASSERT_DOUBLES_EQUAL(vRhoB, dRoundF_roundRhoB, ED_CHECK(vRhoB));
  CPPUNIT_ASSERT_DOUBLES_EQUAL(vGammaAA, dRoundF_roundGammaAA, ED_CHECK(vGammaAA));
  CPPUNIT_ASSERT_DOUBLES_EQUAL(vGammaAB, dRoundF_roundGammaAB, ED_CHECK(vGammaAB));
  CPPUNIT_ASSERT_DOUBLES_EQUAL(vGammaBB, dRoundF_roundGammaBB, ED_CHECK(vGammaBB));
}

//  rhoa= 0.78E-01 rhob= 0.31E-01 sigmaaa= 0.41E-02 sigmaab= 0.38E-02 sigmabb= 0.36E-02
//  zk            = -0.303154148293E-02
//  vrhoa         = -0.244020499979E-01
//  vrhob         = -0.783890652209E-01
//  vsigmaaa      = -0.352365411226E-02
//  vsigmaab      =  0.373879285446E-01
//  vsigmabb      =  0.434388664083E-01

//  v2rhoa2       =  0.202330890903E+00
//  v2rhoab       = -0.536912334084E+00
//  v2rhob2       =  0.182166278214E+01

//  v2rhoasigmaaa = -0.992049673308E-01
//  v2rhoasigmaab = -0.485780496011E-01
//  v2rhoasigmabb = -0.158164241941E+00
//  v2rhobsigmaaa =  0.340893479361E+00
//  v2rhobsigmaab = -0.162192597781E+01
//  v2rhobsigmabb = -0.154746697775E+01

//  v2sigmaaa2    =  0.000000000000E+00
//  v2sigmaaaab   =  0.000000000000E+00
//  v2sigmaaabb   =  0.000000000000E+00
//  v2sigmaab2    =  0.000000000000E+00
//  v2sigmaabbb   =  0.000000000000E+00
//  v2sigmabb2    =  0.000000000000E+00
void DfFunctional_LYPTest::testPointwise19(){
  // input
  const double dRhoA = 0.78E-01;
  const double dRhoB = 0.31E-01;
  const double dGammaAA = 0.41E-02;
  const double dGammaAB = 0.38E-02;
  const double dGammaBB = 0.36E-02;

  // expected value
  const double zk = -0.303154148293E-02;
  const double vRhoA = -0.244020499979E-01;
  const double vRhoB = -0.783890652209E-01;
  const double vGammaAA = -0.352365411226E-02;
  const double vGammaAB =  0.373879285446E-01;
  const double vGammaBB =  0.434388664083E-01;

  // execute test
  DfFunctional_LYP f;
  
  double dFunctionalValue = f.getFunctional(dRhoA, dRhoB, dGammaAA, dGammaAB, dGammaBB);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(zk, dFunctionalValue, ED_CHECK(zk));

  double dRoundF_roundRhoA, dRoundF_roundRhoB;
  double dRoundF_roundGammaAA, dRoundF_roundGammaAB, dRoundF_roundGammaBB;
  f.getDerivativeFunctional(dRhoA, dRhoB, dGammaAA, dGammaAB, dGammaBB,
			    &dRoundF_roundRhoA, &dRoundF_roundRhoB,
			    &dRoundF_roundGammaAA, &dRoundF_roundGammaAB, &dRoundF_roundGammaBB);

  CPPUNIT_ASSERT_DOUBLES_EQUAL(vRhoA, dRoundF_roundRhoA, ED_CHECK(vRhoA));
  CPPUNIT_ASSERT_DOUBLES_EQUAL(vRhoB, dRoundF_roundRhoB, ED_CHECK(vRhoB));
  CPPUNIT_ASSERT_DOUBLES_EQUAL(vGammaAA, dRoundF_roundGammaAA, ED_CHECK(vGammaAA));
  CPPUNIT_ASSERT_DOUBLES_EQUAL(vGammaAB, dRoundF_roundGammaAB, ED_CHECK(vGammaAB));
  CPPUNIT_ASSERT_DOUBLES_EQUAL(vGammaBB, dRoundF_roundGammaBB, ED_CHECK(vGammaBB));
}

//  rhoa= 0.50E+02 rhob= 0.49E+02 sigmaaa= 0.11E+06 sigmaab= 0.11E+06 sigmabb= 0.11E+06
//  zk            = -0.608750998131E+01
//  vrhoa         = -0.643435848234E-01
//  vrhob         = -0.656377677365E-01
//  vsigmaaa      =  0.202448437183E-06
//  vsigmaab      =  0.140789565804E-06
//  vsigmabb      =  0.237924154032E-06

//  v2rhoa2       =  0.644965086508E-03
//  v2rhoab       = -0.630065421347E-03
//  v2rhob2       =  0.684494690431E-03

//  v2rhoasigmaaa = -0.207590458741E-07
//  v2rhoasigmaab = -0.182473154038E-08
//  v2rhoasigmabb =  0.137746587331E-07
//  v2rhobsigmaaa =  0.142543841876E-07
//  v2rhobsigmaab = -0.323521094196E-08
//  v2rhobsigmabb = -0.221599858116E-07

//  v2sigmaaa2    =  0.000000000000E+00
//  v2sigmaaaab   =  0.000000000000E+00
//  v2sigmaaabb   =  0.000000000000E+00
//  v2sigmaab2    =  0.000000000000E+00
//  v2sigmaabbb   =  0.000000000000E+00
//  v2sigmabb2    =  0.000000000000E+00
void DfFunctional_LYPTest::testPointwise20(){
  // input
  const double dRhoA = 0.50E+02;
  const double dRhoB = 0.49E+02;
  const double dGammaAA = 0.11E+06;
  const double dGammaAB = 0.11E+06;
  const double dGammaBB = 0.11E+06;

  // expected value
  const double zk = -0.608750998131E+01;
  const double vRhoA = -0.643435848234E-01;
  const double vRhoB = -0.656377677365E-01;
  const double vGammaAA = 0.202448437183E-06;
  const double vGammaAB = 0.140789565804E-06;
  const double vGammaBB = 0.237924154032E-06;

  // execute test
  DfFunctional_LYP f;
  
  double dFunctionalValue = f.getFunctional(dRhoA, dRhoB, dGammaAA, dGammaAB, dGammaBB);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(zk, dFunctionalValue, ED_CHECK(zk));

  double dRoundF_roundRhoA, dRoundF_roundRhoB;
  double dRoundF_roundGammaAA, dRoundF_roundGammaAB, dRoundF_roundGammaBB;
  f.getDerivativeFunctional(dRhoA, dRhoB, dGammaAA, dGammaAB, dGammaBB,
			    &dRoundF_roundRhoA, &dRoundF_roundRhoB,
			    &dRoundF_roundGammaAA, &dRoundF_roundGammaAB, &dRoundF_roundGammaBB);

  CPPUNIT_ASSERT_DOUBLES_EQUAL(vRhoA, dRoundF_roundRhoA, ED_CHECK(vRhoA));
  CPPUNIT_ASSERT_DOUBLES_EQUAL(vRhoB, dRoundF_roundRhoB, ED_CHECK(vRhoB));
  CPPUNIT_ASSERT_DOUBLES_EQUAL(vGammaAA, dRoundF_roundGammaAA, ED_CHECK(vGammaAA));
  CPPUNIT_ASSERT_DOUBLES_EQUAL(vGammaAB, dRoundF_roundGammaAB, ED_CHECK(vGammaAB));
  CPPUNIT_ASSERT_DOUBLES_EQUAL(vGammaBB, dRoundF_roundGammaBB, ED_CHECK(vGammaBB));
}

//  rhoa= 0.40E+02 rhob= 0.40E+02 sigmaaa= 0.99E+05 sigmaab= 0.98E+05 sigmabb= 0.98E+05
//  zk            = -0.485826625474E+01
//  vrhoa         = -0.653859880589E-01
//  vrhob         = -0.653239498026E-01
//  vsigmaaa      =  0.314305582522E-06
//  vsigmaab      =  0.205017935750E-06
//  vsigmabb      =  0.314305582522E-06

//  v2rhoa2       =  0.932011907377E-03
//  v2rhoab       = -0.845990806989E-03
//  v2rhob2       =  0.927946564209E-03

//  v2rhoasigmaaa = -0.375898792310E-07
//  v2rhoasigmaab = -0.455264241776E-08
//  v2rhoasigmabb =  0.244483769850E-07
//  v2rhobsigmaaa =  0.244483769850E-07
//  v2rhobsigmaab = -0.455264241776E-08
//  v2rhobsigmabb = -0.375898792310E-07

//  v2sigmaaa2    =  0.000000000000E+00
//  v2sigmaaaab   =  0.000000000000E+00
//  v2sigmaaabb   =  0.000000000000E+00
//  v2sigmaab2    =  0.000000000000E+00
//  v2sigmaabbb   =  0.000000000000E+00
//  v2sigmabb2    =  0.000000000000E+00
void DfFunctional_LYPTest::testPointwise21(){
  // input
  const double dRhoA = 0.40E+02;
  const double dRhoB = 0.40E+02;
  const double dGammaAA = 0.99E+05;
  const double dGammaAB = 0.98E+05;
  const double dGammaBB = 0.98E+05;

  // expected value
  const double zk = -0.485826625474E+01;
  const double vRhoA = -0.653859880589E-01;
  const double vRhoB = -0.653239498026E-01;
  const double vGammaAA = 0.314305582522E-06;
  const double vGammaAB = 0.205017935750E-06;
  const double vGammaBB = 0.314305582522E-06;

  // execute test
  DfFunctional_LYP f;
  
  double dFunctionalValue = f.getFunctional(dRhoA, dRhoB, dGammaAA, dGammaAB, dGammaBB);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(zk, dFunctionalValue, ED_CHECK(zk));

  double dRoundF_roundRhoA, dRoundF_roundRhoB;
  double dRoundF_roundGammaAA, dRoundF_roundGammaAB, dRoundF_roundGammaBB;
  f.getDerivativeFunctional(dRhoA, dRhoB, dGammaAA, dGammaAB, dGammaBB,
			    &dRoundF_roundRhoA, &dRoundF_roundRhoB,
			    &dRoundF_roundGammaAA, &dRoundF_roundGammaAB, &dRoundF_roundGammaBB);

  CPPUNIT_ASSERT_DOUBLES_EQUAL(vRhoA, dRoundF_roundRhoA, ED_CHECK(vRhoA));
  CPPUNIT_ASSERT_DOUBLES_EQUAL(vRhoB, dRoundF_roundRhoB, ED_CHECK(vRhoB));
  CPPUNIT_ASSERT_DOUBLES_EQUAL(vGammaAA, dRoundF_roundGammaAA, ED_CHECK(vGammaAA));
  CPPUNIT_ASSERT_DOUBLES_EQUAL(vGammaAB, dRoundF_roundGammaAB, ED_CHECK(vGammaAB));
  CPPUNIT_ASSERT_DOUBLES_EQUAL(vGammaBB, dRoundF_roundGammaBB, ED_CHECK(vGammaBB));
}

//  rhoa= 0.12E+00 rhob= 0.10E+00 sigmaaa= 0.12E+00 sigmaab= 0.13E+00 sigmabb= 0.14E+00
//  zk            = -0.634688010938E-02
//  vrhoa         = -0.486068357646E-01
//  vrhob         = -0.739949355207E-01
//  vsigmaaa      =  0.305337986175E-02
//  vsigmaab      =  0.641200968640E-02
//  vsigmabb      =  0.698202739086E-02

//  v2rhoa2       =  0.433585531118E+00
//  v2rhoab       = -0.365108886112E+00
//  v2rhob2       =  0.100038106732E+01

//  v2rhoasigmaaa = -0.102735565156E+00
//  v2rhoasigmaab = -0.181538750084E-01
//  v2rhoasigmabb =  0.501849138748E-01
//  v2rhobsigmaaa =  0.732361795124E-01
//  v2rhobsigmaab = -0.824126803195E-01
//  v2rhobsigmabb = -0.165270268130E+00

//  v2sigmaaa2    =  0.000000000000E+00
//  v2sigmaaaab   =  0.000000000000E+00
//  v2sigmaaabb   =  0.000000000000E+00
//  v2sigmaab2    =  0.000000000000E+00
//  v2sigmaabbb   =  0.000000000000E+00
//  v2sigmabb2    =  0.000000000000E+00
void DfFunctional_LYPTest::testPointwise22(){
  // input
  const double dRhoA = 0.12E+00;
  const double dRhoB = 0.10E+00;
  const double dGammaAA = 0.12E+00;
  const double dGammaAB = 0.13E+00;
  const double dGammaBB = 0.14E+00;

  // expected value
  const double zk = -0.634688010938E-02;
  const double vRhoA = -0.486068357646E-01;
  const double vRhoB = -0.739949355207E-01;
  const double vGammaAA = 0.305337986175E-02;
  const double vGammaAB = 0.641200968640E-02;
  const double vGammaBB = 0.698202739086E-02;

  // execute test
  DfFunctional_LYP f;
  
  double dFunctionalValue = f.getFunctional(dRhoA, dRhoB, dGammaAA, dGammaAB, dGammaBB);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(zk, dFunctionalValue, ED_CHECK(zk));

  double dRoundF_roundRhoA, dRoundF_roundRhoB;
  double dRoundF_roundGammaAA, dRoundF_roundGammaAB, dRoundF_roundGammaBB;
  f.getDerivativeFunctional(dRhoA, dRhoB, dGammaAA, dGammaAB, dGammaBB,
			    &dRoundF_roundRhoA, &dRoundF_roundRhoB,
			    &dRoundF_roundGammaAA, &dRoundF_roundGammaAB, &dRoundF_roundGammaBB);

  CPPUNIT_ASSERT_DOUBLES_EQUAL(vRhoA, dRoundF_roundRhoA, ED_CHECK(vRhoA));
  CPPUNIT_ASSERT_DOUBLES_EQUAL(vRhoB, dRoundF_roundRhoB, ED_CHECK(vRhoB));
  CPPUNIT_ASSERT_DOUBLES_EQUAL(vGammaAA, dRoundF_roundGammaAA, ED_CHECK(vGammaAA));
  CPPUNIT_ASSERT_DOUBLES_EQUAL(vGammaAB, dRoundF_roundGammaAB, ED_CHECK(vGammaAB));
  CPPUNIT_ASSERT_DOUBLES_EQUAL(vGammaBB, dRoundF_roundGammaBB, ED_CHECK(vGammaBB));
}

//  rhoa= 0.48E-01 rhob= 0.25E-01 sigmaaa= 0.46E-02 sigmaab= 0.44E-02 sigmabb= 0.41E-02
//  zk            = -0.172301075975E-02
//  vrhoa         = -0.294822900234E-01
//  vrhob         = -0.766946918598E-01
//  vsigmaaa      =  0.127140367254E-02
//  vsigmaab      =  0.519446045760E-01
//  vsigmabb      =  0.608118658277E-01

//  v2rhoa2       =  0.427000673980E+00
//  v2rhoab       = -0.671140471911E+00
//  v2rhob2       =  0.300158916489E+01

//  v2rhoasigmaaa = -0.631566175836E+00
//  v2rhoasigmaab = -0.540623931547E-01
//  v2rhoasigmabb = -0.252818541902E-02
//  v2rhobsigmaaa =  0.988123041938E+00
//  v2rhobsigmaab = -0.292607572225E+01
//  v2rhobsigmabb = -0.333300175870E+01

//  v2sigmaaa2    =  0.000000000000E+00
//  v2sigmaaaab   =  0.000000000000E+00
//  v2sigmaaabb   =  0.000000000000E+00
//  v2sigmaab2    =  0.000000000000E+00
//  v2sigmaabbb   =  0.000000000000E+00
//  v2sigmabb2    =  0.000000000000E+00
void DfFunctional_LYPTest::testPointwise23(){
  // input
  const double dRhoA = 0.48E-01;
  const double dRhoB = 0.25E-01;
  const double dGammaAA = 0.46E-02;
  const double dGammaAB = 0.44E-02;
  const double dGammaBB = 0.41E-02;

  // expected value
  const double zk = -0.172301075975E-02;
  const double vRhoA = -0.294822900234E-01;
  const double vRhoB = -0.766946918598E-01;
  const double vGammaAA = 0.127140367254E-02;
  const double vGammaAB = 0.519446045760E-01;
  const double vGammaBB = 0.608118658277E-01;

  // execute test
  DfFunctional_LYP f;
  
  double dFunctionalValue = f.getFunctional(dRhoA, dRhoB, dGammaAA, dGammaAB, dGammaBB);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(zk, dFunctionalValue, ED_CHECK(zk));

  double dRoundF_roundRhoA, dRoundF_roundRhoB;
  double dRoundF_roundGammaAA, dRoundF_roundGammaAB, dRoundF_roundGammaBB;
  f.getDerivativeFunctional(dRhoA, dRhoB, dGammaAA, dGammaAB, dGammaBB,
			    &dRoundF_roundRhoA, &dRoundF_roundRhoB,
			    &dRoundF_roundGammaAA, &dRoundF_roundGammaAB, &dRoundF_roundGammaBB);

  CPPUNIT_ASSERT_DOUBLES_EQUAL(vRhoA, dRoundF_roundRhoA, ED_CHECK(vRhoA));
  CPPUNIT_ASSERT_DOUBLES_EQUAL(vRhoB, dRoundF_roundRhoB, ED_CHECK(vRhoB));
  CPPUNIT_ASSERT_DOUBLES_EQUAL(vGammaAA, dRoundF_roundGammaAA, ED_CHECK(vGammaAA));
  CPPUNIT_ASSERT_DOUBLES_EQUAL(vGammaAB, dRoundF_roundGammaAB, ED_CHECK(vGammaAB));
  CPPUNIT_ASSERT_DOUBLES_EQUAL(vGammaBB, dRoundF_roundGammaBB, ED_CHECK(vGammaBB));
}
