#include <vector>
#include <string>
#include "TlStringTokenizerTest.h"

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

void TlStringTokenizerTest::testConstructer(){
  TlStringTokenizer st("This is a test string.");
}

void TlStringTokenizerTest::testCountTokens(){
  TlStringTokenizer st("This is a test string.");

  int count = st.countTokens();

  CPPUNIT_ASSERT_EQUAL(5, count);
}

void TlStringTokenizerTest::testGetTokens(){
  TlStringTokenizer st("This is a test string.");
  
  std::vector<std::string> tokens;
  while (st.hasMoreTokens() == true){
    std::string s = st.nextToken();
    tokens.push_back(s);
  }

  CPPUNIT_ASSERT_EQUAL(size_t(5), tokens.size());
  CPPUNIT_ASSERT_EQUAL(std::string("This"),    tokens[0]);
  CPPUNIT_ASSERT_EQUAL(std::string("is"),      tokens[1]);
  CPPUNIT_ASSERT_EQUAL(std::string("a"),       tokens[2]);
  CPPUNIT_ASSERT_EQUAL(std::string("test"),    tokens[3]);
  CPPUNIT_ASSERT_EQUAL(std::string("string."), tokens[4]);
}


