#include <string>
#include <vector>
#include "TlMemManager.h"
#include "gtest/gtest.h"

TEST(TlMemManager, allocat) {
    const std::size_t size = 1024 * 1024 * 1024;  // 1 GB
    const std::string filePath = "/tmp/temp.map";
    TlMemManager::setParam(size, filePath);

    TlMemManager& rMem = TlMemManager::getInstance();

    const std::size_t p1_size = 100;
    char* p1 = rMem.allocate(p1_size);
    // rMem.debugOutFreeMemList();

    const std::size_t p2_size = 100;
    char* p2 = rMem.allocate(p2_size);
    // rMem.debugOutFreeMemList();

    ASSERT_EQ(100, int(p2 - p1));

    // 再確保
    rMem.deallocate(p2);
    // rMem.debugOutFreeMemList();

    p2 = rMem.allocate(p2_size);
    // rMem.debugOutFreeMemList();

    ASSERT_EQ(100, int(p2 - p1));

    // 追加確保
    const std::size_t p3_size = 100;
    char* p3 = rMem.allocate(p3_size);
    // rMem.debugOutFreeMemList();

    ASSERT_EQ(100, int(p3 - p2));
    ASSERT_EQ(200, int(p3 - p1));

    // 中間解放 -> 空間内確保
    rMem.deallocate(p2);
    // rMem.debugOutFreeMemList();

    const std::size_t p4_size = 50;
    char* p4 = rMem.allocate(p4_size);
    ASSERT_EQ(100, int(p4 - p1));  // p2の場所にできるはず
    // rMem.debugOutFreeMemList();

    // 空間より大きな領域確保
    const std::size_t p5_size = 100;
    char* p5 = rMem.allocate(p5_size);
    ASSERT_EQ(300, int(p5 - p1));
    // rMem.debugOutFreeMemList();
}
