#include <limits>
#include "Fl_UserinputXTest.h"
#include "TlParameter.h"

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

void Fl_UserinputXTest::testConstructer(){
  Fl_UserinputX flUserinput;
}

void Fl_UserinputXTest::testGetFlGlobalinputX(){
  Fl_UserinputX flUserinput("Fl_UserinputXTest.txt");
  flUserinput.load();
  TlParameter param = flUserinput.getParameter();
  
  CPPUNIT_ASSERT_EQUAL(std::string("create integral guess scf"), param["MAIN"]    ["step-control"]);
  CPPUNIT_ASSERT_EQUAL(std::string("1.0e-16"),                   param["INTEGRAL"]["cut-value"]);
}



