#include <string>
#include <vector>
#include "TlSerializeData.h"
#include "gtest/gtest.h"

TEST(TlSerializeData, constructor) {
    TlSerializeData so("Hello World!");

    const std::string tmp = so.getStr();

    EXPECT_EQ(std::string("Hello World!"), tmp);
}

TEST(TlSerializeData, copyConstructor) {
    TlSerializeData root;
    TlSerializeData chiyoda_line;
    chiyoda_line.pushBack("yoyogi-uehara");
    chiyoda_line.pushBack("yoyogi-kouen");
    chiyoda_line.pushBack("meiji-jinguumae");
    chiyoda_line.pushBack("omotesandou");
    root.add(TlSerializeData("chiyoda-line"), chiyoda_line);

    TlSerializeData root2 = root;
}

TEST(TlSerializeData, mapAccess) {
    TlSerializeData root;
    root["chiyoda-line"] = "unknown";
    root["marunouchi-line"] = "red";
    root["ginza-line"] = "orange";

    TlSerializeData chiyoda;
    chiyoda["yoyogi-uehara"] = "C1";
    chiyoda["yoyogi-kouen"] = "C2";
    chiyoda["meiji-jinguumae"] = "C3";
    root["chiyoda-line"] = chiyoda;  // overwrite

    EXPECT_EQ(std::string("red"), root["marunouchi-line"].getStr());
    EXPECT_EQ(std::string("orange"), root["ginza-line"].getStr());
    EXPECT_EQ(std::string("C2"), root["chiyoda-line"]["yoyogi-kouen"].getStr());
}

TEST(TlSerializeData, converVectorInt) {
    const int size = 10;
    std::vector<int> v(size);
    for (int i = 0; i < size; ++i) {
        v[i] = i;
    }

    TlSerializeData root = v;
    EXPECT_EQ(size, root.getSize());
    for (int i = 0; i < size; ++i) {
        EXPECT_EQ(i, root.getAt(i).getInt());
    }
}
