#include <vector>
#include <string>
#include "TlMsgPackTest.h"
#include "TlSerializeData.h"

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

void TlMsgPackTest::testLoad()
{
    TlMsgPack mpack;
    mpack.load("sample.mpac");

//     TlSerializeData so = mpack.getSerializeObject();
//     std::string str = so.str();
//     std::cout << str << std::endl;
    
    //CPPUNIT_ASSERT_EQUAL(std::string("value"), a["group"]["keyword"]);
}

void TlMsgPackTest::testDumpAndPack()
{
    TlMsgPack mpack;
    mpack.load("sample.mpac");

//     TlSerializeData so = mpack.getSerializeObject();
//     std::string str = so.str();
//     std::cout << str << std::endl;

    std::string packStr = mpack.dump();

    TlMsgPack mpack2;
    mpack2.pack(packStr);

    TlSerializeData so2 = mpack2.getSerializeData();
    std::string str2 = so2.str();
//     std::cout << str2 << std::endl;
    
    //CPPUNIT_ASSERT_EQUAL(std::string("value"), a["group"]["keyword"]);
}

