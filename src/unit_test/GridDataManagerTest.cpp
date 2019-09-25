#include "GridDataManagerTest.h"
#include <limits>
#include "TlFile.h"

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
// 得られた結果actualが期待する値expectedでなかったとき、すなわちexpected !=
// actualのときに失敗します。
//
// CPPUNIT_ASSERT_EQUAL_MESSAGE( message, expected, actual );
// 得られた結果actualが期待する値expectedでなかったとき、すなわちexpected !=
// actualのときに失敗します。このときmessageを出力します。
//
// CPPUNIT_ASSERT_DOUBLES_EQUAL( expected, actual, delta );
// 得られた結果actualと期待する値expectedとの差がdeltaより大きいとき、失敗します。

void GridDataManagerTest::testConstructer() {
    if (TlFile::isExistFile("GridDataManager.dat") == true) {
        TlFile::remove("GridDataManager.dat");
    }
    GridDataManager gdm("GridDataManager.dat");
}

void GridDataManagerTest::testSetDensity() {
    GridDataManager gdm("GridDataManager.dat");

    std::vector<double> v1(10);
    for (int i = 0; i < 10; ++i) {
        v1[i] = double(i);
    }
    gdm.setData(1, GridDataManager::DENSITY, v1);

    std::vector<double> v2(10);
    for (int i = 0; i < 10; ++i) {
        v2[i] = double(i * 2);
    }
    gdm.setData(2, GridDataManager::DENSITY, v2);

    std::vector<double> v3(10);
    for (int i = 0; i < 10; ++i) {
        v3[i] = double(i * 3);
    }
    gdm.setData(3, GridDataManager::DENSITY, v3);
}

void GridDataManagerTest::testGetDensity() {
    GridDataManager gdm("GridDataManager.dat");

    std::vector<double> v1 = gdm.getData(1, GridDataManager::DENSITY);
    const std::size_t v1_size = 10;
    CPPUNIT_ASSERT_EQUAL(v1_size, v1.size());
    for (std::size_t i = 0; i < v1_size; ++i) {
        CPPUNIT_ASSERT_DOUBLES_EQUAL(v1[i], double(i), 1.0E-16);
    }

    std::vector<double> v2 = gdm.getData(2, GridDataManager::DENSITY);
    const std::size_t v2_size = 10;
    CPPUNIT_ASSERT_EQUAL(v2_size, v2.size());
    for (std::size_t i = 0; i < v2_size; ++i) {
        CPPUNIT_ASSERT_DOUBLES_EQUAL(v2[i], double(i * 2), 1.0E-16);
    }
}

void GridDataManagerTest::testUpdate() {
    // prepare
    {
        GridDataManager gdm("GridDataManager.dat");

        std::vector<double> v1(10);
        for (int i = 0; i < 10; ++i) {
            v1[i] = double(i);
        }
        gdm.setData(1, GridDataManager::DENSITY, v1);

        std::vector<double> v2(10);
        for (int i = 0; i < 10; ++i) {
            v2[i] = double(i * 2);
        }
        gdm.setData(2, GridDataManager::DENSITY, v2);

        std::vector<double> v3(10);
        for (int i = 0; i < 10; ++i) {
            v3[i] = double(i * 3);
        }
        gdm.setData(3, GridDataManager::DENSITY, v3);

        // update v2
        v2.resize(20, 0.0);
        for (int i = 0; i < 20; ++i) {
            v2[i] = double(i * 2);
        }
        gdm.setData(2, GridDataManager::DENSITY, v2);
    }

    // check
    {
        GridDataManager gdm("GridDataManager.dat");

        std::vector<double> v1 = gdm.getData(1, GridDataManager::DENSITY);
        const std::size_t v1_size = 10;
        CPPUNIT_ASSERT_EQUAL(v1_size, v1.size());
        for (std::size_t i = 0; i < v1_size; ++i) {
            CPPUNIT_ASSERT_DOUBLES_EQUAL(v1[i], double(i), 1.0E-16);
        }

        std::vector<double> v2 = gdm.getData(2, GridDataManager::DENSITY);
        const std::size_t v2_size = 20;
        CPPUNIT_ASSERT_EQUAL(v2_size, v2.size());
        for (std::size_t i = 0; i < v2_size; ++i) {
            CPPUNIT_ASSERT_DOUBLES_EQUAL(v2[i], double(i * 2), 1.0E-16);
        }

        std::vector<double> v3 = gdm.getData(3, GridDataManager::DENSITY);
        const std::size_t v3_size = 10;
        CPPUNIT_ASSERT_EQUAL(v3_size, v3.size());
        for (std::size_t i = 0; i < v3_size; ++i) {
            CPPUNIT_ASSERT_DOUBLES_EQUAL(v3[i], double(i * 3), 1.0E-16);
        }
    }
}
