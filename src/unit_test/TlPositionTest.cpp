#include <vector>
#include <string>
#include <limits>
#include "TlPositionTest.h"

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

const double TlPositionTest::threshold = std::numeric_limits<double>::epsilon();

void TlPositionTest::testConstructer(){
  TlPosition a;
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, a.x(), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, a.y(), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, a.z(), threshold);

  TlPosition b(1.0, 2.0, 3.0);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, b.x(), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(2.0, b.y(), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(3.0, b.z(), threshold);
}

void TlPositionTest::testCopyConstructer(){
  TlPosition a(1.0, 2.0, 3.0);
  TlPosition b(a);

  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, b.x(), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(2.0, b.y(), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(3.0, b.z(), threshold);
}

void TlPositionTest::testOperatorEqual(){
  TlPosition a(1.0, 2.0, 3.0);
  TlPosition b(2.0, 3.0, 4.0);

  b = a;
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, b.x(), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(2.0, b.y(), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(3.0, b.z(), threshold);
}

void TlPositionTest::testSquareDistanceFrom(){
  TlPosition a(1.0, 2.0, 3.0);
  TlPosition b(4.0, 5.0, 6.0);

  double squareDistance = a.squareDistanceFrom(b);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(155.0, squareDistance, threshold);
}


