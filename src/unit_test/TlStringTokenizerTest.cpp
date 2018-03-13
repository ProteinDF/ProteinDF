#include "TlStringTokenizer.h"
#include <string>
#include <vector>
#include "gtest/gtest.h"

TEST(TlStringTokenizer, countTokens) {
  TlStringTokenizer st("This is a test string.");

  int count = st.countTokens();

  EXPECT_EQ(5, count);
}

TEST(TlStringTokenizer, testGetTokens) {
  TlStringTokenizer st("This is a test string.");

  std::vector<std::string> tokens;
  while (st.hasMoreTokens() == true) {
    std::string s = st.nextToken();
    tokens.push_back(s);
  }

  EXPECT_EQ(size_t(5), tokens.size());
  EXPECT_EQ(std::string("This"), tokens[0]);
  EXPECT_EQ(std::string("is"), tokens[1]);
  EXPECT_EQ(std::string("a"), tokens[2]);
  EXPECT_EQ(std::string("test"), tokens[3]);
  EXPECT_EQ(std::string("string."), tokens[4]);
}
