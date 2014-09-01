#include <limits>
#include "Fl_UserinputXTest.h"
#include "TlParameter.h"

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



