#include "TlFile.h"
#include "gtest/gtest.h"

TEST(TlFile, dir) {
    const std::string dir = TlFile::dir("/home/abc/def.ghi");
    EXPECT_EQ(std::string("/home/abc"), dir);
}

TEST(TlFile, filename) {
    const std::string filename = TlFile::filename("/home/abc/def.ghi");
    EXPECT_EQ(std::string("def.ghi"), filename);
}

TEST(TlFile, stem) {
    const std::string stem = TlFile::stem("/home/abc/def.ghi");
    EXPECT_EQ(std::string("def"), stem);
}

TEST(TlFile, extension) {
    const std::string ext = TlFile::extension("/home/abc/def.ghi");
    EXPECT_EQ(std::string(".ghi"), ext);
}

TEST(TlFile, isExistFile) {
    bool a = TlFile::isExistFile("/bin/bash");
    bool b = TlFile::isExistFile("unexisted_file");
    // bool c = TlFile::isExistFile("/tmp");

    EXPECT_TRUE(a);
    EXPECT_FALSE(b);
    // EXPECT_FALSE(c);
}

TEST(TlFile, getFileSize) {
    std::size_t length = TlFile::getFileSize("/bin/bash");
    std::cout << "length = " << length << " byte." << std::endl;
}
