#include <vector>
#include <string>
#include "TlUtilsTest.h"

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



