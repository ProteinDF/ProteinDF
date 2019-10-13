#include "gtest/gtest.h"

#include "TlFile.h"
#include "TlMatrix.h"
#include "TlMmapMatrix.h"
#include "common.h"

static const std::string matPath = "temp.mmapmat.mat";

static void cleanup() {
    if (TlFile::isExist(matPath)) {
        TlFile::remove(matPath);
    }
}

TEST(TlMmapMatrix, construct) {
    const int row = 100;
    const int col = 200;
    {
        cleanup();
        TlMmapMatrix m(matPath, row, col);
    }

    {
        TlMatrix a;
        a.load(matPath);

        EXPECT_EQ(row, a.getNumOfRows());
        EXPECT_EQ(col, a.getNumOfCols());
    }
}

TEST(TlMmapMatrix, resize) {
    cleanup();
    const int row = 100;
    const int col = 200;
    TlMatrix ref = getTlMatrix(row, col);
    ref.save(matPath);
    TlMmapMatrix a(matPath);

    const int newRow = 200;
    const int newCol = 200;
    ref.resize(newRow, newCol);
    a.resize(newRow, newCol);

    EXPECT_EQ(newRow, a.getNumOfRows());
    EXPECT_EQ(newCol, a.getNumOfCols());
    for (int i = 0; i < newRow; ++i) {
        for (int j = 0; j < newCol; ++j) {
            EXPECT_DOUBLE_EQ(ref.get(i, j), a.get(i, j));
        }
    }
}

TEST(TlMmapMatrix, set_matrix) {
    cleanup();
    const int row = 100;
    const int col = 200;
    TlMatrix ref = getTlMatrix(row, col);

    {
        TlMmapMatrix m(matPath, ref.getNumOfRows(), ref.getNumOfCols());
        for (int r = 0; r < ref.getNumOfRows(); ++r) {
            for (int c = 0; c < ref.getNumOfCols(); ++c) {
                m.set(r, c, ref.get(r, c));
            }
        }
    }

    {
        TlMatrix a;
        a.load(matPath);

        EXPECT_EQ(ref.getNumOfRows(), a.getNumOfRows());
        EXPECT_EQ(ref.getNumOfCols(), a.getNumOfCols());
        for (int r = 0; r < ref.getNumOfRows(); ++r) {
            for (int c = 0; c < ref.getNumOfCols(); ++c) {
                EXPECT_DOUBLE_EQ(ref.get(r, c), a.get(r, c));
            }
        }
    }
}

TEST(TlMmapMatrix, get_matrix) {
    cleanup();
    const int row = 200;
    const int col = 100;
    TlMatrix ref = getTlMatrix(row, col);
    ref.save(matPath);

    {
        TlMmapMatrix m(matPath);
        EXPECT_EQ(ref.getNumOfRows(), m.getNumOfRows());
        EXPECT_EQ(ref.getNumOfCols(), m.getNumOfCols());
        for (int r = 0; r < ref.getNumOfRows(); ++r) {
            for (int c = 0; c < ref.getNumOfCols(); ++c) {
                EXPECT_DOUBLE_EQ(ref.get(r, c), m.get(r, c));
            }
        }
    }
}
