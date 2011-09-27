#include <vector>
#include <string>
#include "TlParameterTest.h"

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

void TlParameterTest::testConstructer(){
  TlParameter a;
  a["group"]["keyword"] = "value";

  CPPUNIT_ASSERT_EQUAL(std::string("value"), a["group"]["keyword"]);
}

void TlParameterTest::testCopyConstructer(){
  TlParameter a;
  a["group"]["keyword"] = "value";

  TlParameter b = a;
  CPPUNIT_ASSERT_EQUAL(std::string("value"), b["group"]["keyword"]);
}

void TlParameterTest::testOperatorEqual(){
  TlParameter a;
  a["group"]["keyword"] = "hoge";

  TlParameter b;
  b["group"]["keyword"] = "foo";
  CPPUNIT_ASSERT_EQUAL(std::string("foo"), b["group"]["keyword"]);

  b = a;
  CPPUNIT_ASSERT_EQUAL(std::string("hoge"), b["group"]["keyword"]);
}
