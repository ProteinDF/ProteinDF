#include <string>
#include <vector>

#include "TlUtils.h"
#include "gtest/gtest.h"

TEST(TlUtils, join) {
    std::vector<std::string> strings(4);
    strings[0] = "This";
    strings[1] = "is";
    strings[2] = "a";
    strings[3] = "pen";

    const std::string sentense = TlUtils::join(strings, " ");

    EXPECT_EQ(std::string("This is a pen"), sentense);
}

TEST(TlUtils, pad) {
    std::string str1 = "line";
    TlUtils::pad(str1, 10, '=');

    EXPECT_EQ(std::string("line======"), str1);
}

TEST(TlUtils, trim) {
    std::string str1 = "aaaaaThis is a pen.";
    TlUtils::trim(str1, 'a');

    EXPECT_EQ(std::string("This is a pen."), str1);
}

TEST(TlUtils, trim_ws) {
    std::string str1 = "     This is a pen.";
    TlUtils::trim_ws(str1);

    EXPECT_EQ(std::string("This is a pen."), str1);
}
