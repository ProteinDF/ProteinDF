#include <vector>
#include <string>
#include "TlSerializeDataTest.h"

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

void TlSerializeDataTest::testConstructor()
{
    TlSerializeData so("Hello World!");

    const std::string tmp = so.getStr();
    
    CPPUNIT_ASSERT_EQUAL(std::string("Hello World!"), tmp);
}

void TlSerializeDataTest::testCopyConstructor()
{
    TlSerializeData root;
    TlSerializeData chiyoda_line;
    chiyoda_line.pushBack("yoyogi-uehara");
    chiyoda_line.pushBack("yoyogi-kouen");
    chiyoda_line.pushBack("meiji-jinguumae");
    chiyoda_line.pushBack("omotesandou");
    root.add(TlSerializeData("chiyoda-line"), chiyoda_line);
    
    TlSerializeData root2 = root;
    
}


void TlSerializeDataTest::testMapAccess()
{
    TlSerializeData root;
    root["chiyoda-line"] = "unknown";
    root["marunouchi-line"] = "red";
    root["ginza-line"] = "orange";

    TlSerializeData chiyoda;
    chiyoda["yoyogi-uehara"] = "C1";
    chiyoda["yoyogi-kouen"] = "C2";
    chiyoda["meiji-jinguumae"] = "C3";
    root["chiyoda-line"] = chiyoda; // overwrite
    
    CPPUNIT_ASSERT_EQUAL(std::string("red"), root["marunouchi-line"].getStr());
    CPPUNIT_ASSERT_EQUAL(std::string("orange"), root["ginza-line"].getStr());
    CPPUNIT_ASSERT_EQUAL(std::string("C2"), root["chiyoda-line"]["yoyogi-kouen"].getStr());
}
