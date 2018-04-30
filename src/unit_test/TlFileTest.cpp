#include "TlFile.h"
#include "gtest/gtest.h"

TEST(TlFIle, isExistFile) {
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
