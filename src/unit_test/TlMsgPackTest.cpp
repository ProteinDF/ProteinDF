#include <string>
#include <vector>
#include "TlMsgPack.h"
#include "TlSerializeData.h"
#include "gtest/gtest.h"

TEST(TlMsgPack, load) {
    TlMsgPack mpack;
    mpack.load("sample.mpac");

    //     TlSerializeData so = mpack.getSerializeObject();
    //     std::string str = so.str();
    //     std::cout << str << std::endl;

    // EXPECT_EQ(std::string("value"), a["group"]["keyword"]);
}

TEST(TlMsgPack, dumpAndPack) {
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

    // EXPECT_EQ(std::string("value"), a["group"]["keyword"]);
}

TEST(TlMsgPack, emptyStringInMapAndArray) {
    TlSerializeData root;
    root["empty_key"] = std::string("");
    root["non_empty_key"] = std::string("hello");

    TlSerializeData arr;
    arr.pushBack(std::string(""));
    arr.pushBack(std::string("world"));
    root["empty_array"] = arr;

    TlMsgPack packer(root);
    const std::string packStr = packer.dump();

    TlMsgPack unpacker;
    unpacker.pack(packStr);
    TlSerializeData result = unpacker.getSerializeData();

    EXPECT_EQ(TlSerializeData::STRING, result["empty_key"].getType());
    EXPECT_EQ(std::string(""), result["empty_key"].getStr());

    EXPECT_EQ(TlSerializeData::STRING, result["non_empty_key"].getType());
    EXPECT_EQ(std::string("hello"), result["non_empty_key"].getStr());

    EXPECT_EQ(TlSerializeData::ARRAY, result["empty_array"].getType());
    EXPECT_EQ(TlSerializeData::STRING, result["empty_array"].getAt(0).getType());
    EXPECT_EQ(std::string(""), result["empty_array"].getAt(0).getStr());
    EXPECT_EQ(TlSerializeData::STRING, result["empty_array"].getAt(1).getType());
    EXPECT_EQ(std::string("world"), result["empty_array"].getAt(1).getStr());
}
