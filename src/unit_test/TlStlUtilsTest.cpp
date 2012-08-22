#include <vector>
#include <string>
#include "TlStlUtilsTest.h"

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

void TlStlUtilsTest::testCache() {
    TlCache<int, std::string> cache(10);
    CPPUNIT_ASSERT_EQUAL(std::size_t(10), cache.getMaxItems());

    cache.setMaxItems(8);
    CPPUNIT_ASSERT_EQUAL(std::size_t(8), cache.getMaxItems());

    cache.set(1, "one");
    cache.set(2, "two");
    cache.set(3, "three");
    cache.set(4, "four");
    cache.set(5, "five");

    CPPUNIT_ASSERT_EQUAL(std::size_t(5), cache.getNumOfItems());
    CPPUNIT_ASSERT_EQUAL(std::string("three"), *(cache.get(3)));

    cache.set(6, "six");
    cache.set(7, "seven");
    cache.set(8, "eight");
    CPPUNIT_ASSERT_EQUAL(std::size_t(8), cache.getNumOfItems());

    cache.set(9, "nine");
    CPPUNIT_ASSERT_EQUAL(std::size_t(8), cache.getNumOfItems());

    cache.set(10, "ten");
    CPPUNIT_ASSERT_EQUAL(std::string("five"), *(cache.get(5)));

    // empty
    CPPUNIT_ASSERT_EQUAL(std::string(""), *(cache.get(1)));
}

