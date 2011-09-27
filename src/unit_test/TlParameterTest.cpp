#include <vector>
#include <string>
#include "TlParameterTest.h"

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
