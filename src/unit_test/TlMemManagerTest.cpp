#include <vector>
#include <string>
#include "TlMemManagerTest.h"

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

void TlMemManagerTest::testAllocate() {
  const std::size_t size = 1024 * 1024 * 1024; // 1 GB
  const std::string filePath = "temp.map";
  TlMemManager::setParam(size, filePath);

  TlMemManager& rMem = TlMemManager::getInstance();

  const std::size_t p1_size = 100;
  char* p1 = rMem.allocate(p1_size);
  //rMem.debugOutFreeMemList();

  const std::size_t p2_size = 100;
  char* p2 = rMem.allocate(p2_size);
  //rMem.debugOutFreeMemList();
  
  CPPUNIT_ASSERT_EQUAL(100, int(p2 - p1));

  // 再確保
  rMem.deallocate(p2, p2_size);
  //rMem.debugOutFreeMemList();

  p2 = rMem.allocate(p2_size);
  //rMem.debugOutFreeMemList();

  CPPUNIT_ASSERT_EQUAL(100, int(p2 - p1));

  // 追加確保
  const std::size_t p3_size = 100;
  char* p3 = rMem.allocate(p3_size);
  //rMem.debugOutFreeMemList();

  CPPUNIT_ASSERT_EQUAL(100, int(p3 - p2));
  CPPUNIT_ASSERT_EQUAL(200, int(p3 - p1));

  // 中間解放 -> 空間内確保
  rMem.deallocate(p2, p2_size);
  //rMem.debugOutFreeMemList();

  const std::size_t p4_size = 50;
  char* p4 = rMem.allocate(p4_size);
  CPPUNIT_ASSERT_EQUAL(100, int(p4 - p1)); // p2の場所にできるはず
  //rMem.debugOutFreeMemList();

  // 空間より大きな領域確保
  const std::size_t p5_size = 100;
  char* p5 = rMem.allocate(p5_size);
  CPPUNIT_ASSERT_EQUAL(300, int(p5 - p1));
  //rMem.debugOutFreeMemList();
}


