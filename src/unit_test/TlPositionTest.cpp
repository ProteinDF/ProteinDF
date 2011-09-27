#include <vector>
#include <string>
#include <limits>
#include "TlPositionTest.h"

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


