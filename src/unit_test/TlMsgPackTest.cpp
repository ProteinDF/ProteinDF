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
