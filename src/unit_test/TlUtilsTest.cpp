#include <vector>
#include <string>
#include "TlUtilsTest.h"

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

void TlUtilsTest::testPad(){
  std::string str1 = "line";
  TlUtils::pad(str1, 10, '=');

  CPPUNIT_ASSERT_EQUAL(std::string("line======"), str1);
}

void TlUtilsTest::testTrim(){
  std::string str1 = "aaaaaThis is a pen.";
  TlUtils::trim(str1, 'a');

  CPPUNIT_ASSERT_EQUAL(std::string("This is a pen."), str1);
}

void TlUtilsTest::testTrim_ws(){
  std::string str1 = "     This is a pen.";
  TlUtils::trim_ws(str1);

  CPPUNIT_ASSERT_EQUAL(std::string("This is a pen."), str1);
}



