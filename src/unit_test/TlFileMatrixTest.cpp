#include <stdio.h>
#include <time.h>

#include "TlFile.h"
#include "TlFileMatrix.h"
#include "TlMatrix.h"
#include "gtest/gtest.h"

static const std::string mat_path = "temp.filemat.mat";

static void cleanup() {
  if (TlFile::isExistFile(mat_path)) {
    TlFile::remove(mat_path);
  }
}

TlMatrix getMatrix(int row, int col) {
  TlMatrix m(row, col);

  std::srand((unsigned int)time(NULL));
  for (int r = 0; r < row; ++r) {
    for (int c = 0; c < col; ++c) {
      m.set(r, c, double(double(std::rand()) / double(RAND_MAX)));
    }
  }

  return m;
}

TEST(TlFileMatrix, construct) {
  const int row = 100;
  const int col = 200;
  {
    cleanup();
    TlFileMatrix m(mat_path, row, col);

    EXPECT_EQ(row, m.getNumOfRows());
    EXPECT_EQ(col, m.getNumOfCols());
  }
  EXPECT_EQ(true, TlFile::isExistFile(mat_path));

  {
    TlMatrix a;
    a.load(mat_path);

    EXPECT_EQ(row, a.getNumOfRows());
    EXPECT_EQ(col, a.getNumOfCols());
  }
}

TEST(TlFileMatrix, construct_by_existed) {
  cleanup();

  TlMatrix m(100, 200);
  {
    m(3, 17) = 51.0;
    m(0, 28) = -1.0;
    m.save(mat_path);
  }

  {
    TlFileMatrix a(mat_path);
    EXPECT_EQ(100, a.getNumOfRows());
    EXPECT_EQ(200, a.getNumOfCols());
    for (int r = 0; r < 100; ++r) {
      for (int c = 0; c <= r; ++c) {
        EXPECT_DOUBLE_EQ(m(r, c), a.get(r, c));
      }
    }
  }
}

TEST(TlFileMatrix, get_from_saved_matrix) {
  cleanup();

  TlMatrix m(100, 200);
  {
    m(3, 17) = 51.0;
    m(0, 28) = -1.0;
    m.save(mat_path);
  }

  {
    TlFileMatrix fm(mat_path);
    EXPECT_DOUBLE_EQ(0.0, fm.get(0, 0));
    EXPECT_DOUBLE_EQ(0.0, fm.get(0, 99));
    EXPECT_DOUBLE_EQ(0.0, fm.get(99, 0));
    EXPECT_DOUBLE_EQ(0.0, fm.get(99, 99));
    EXPECT_DOUBLE_EQ(51.0, fm.get(3, 17));
    EXPECT_DOUBLE_EQ(-1.0, fm.get(0, 28));
  }
}

TEST(TlFileMatrix, set) {
  cleanup();

  {
    TlFileMatrix fm(mat_path, 100, 200);
    fm.set(3, 17, 51.0);
    fm.set(0, 28, -1.0);
  }

  {
    TlMatrix m;
    m.load(mat_path);

    EXPECT_EQ(100, m.getNumOfRows());
    EXPECT_EQ(200, m.getNumOfCols());
    EXPECT_DOUBLE_EQ(0.0, m(0, 0));
    EXPECT_DOUBLE_EQ(0.0, m(0, 99));
    EXPECT_DOUBLE_EQ(0.0, m(99, 0));
    EXPECT_DOUBLE_EQ(0.0, m(99, 99));
    EXPECT_DOUBLE_EQ(51.0, m(3, 17));
    EXPECT_DOUBLE_EQ(-1.0, m(0, 28));
  }
}

TEST(TlFileMatrix, add) {
  cleanup();

  // prepare
  {
    TlMatrix m(100, 200);
    m(3, 17) = 51.0;
    m(0, 28) = -1.0;
    m.save(mat_path);
  }

  {
    TlFileMatrix fsm(mat_path);
    fsm.add(3, 17, 12.3);
  }

  // check
  {
    TlMatrix m;
    m.load(mat_path);
    EXPECT_EQ(100, m.getNumOfRows());
    EXPECT_EQ(200, m.getNumOfCols());
    EXPECT_DOUBLE_EQ(0.0, m(0, 0));
    EXPECT_DOUBLE_EQ(0.0, m(0, 99));
    EXPECT_DOUBLE_EQ(0.0, m(99, 0));
    EXPECT_DOUBLE_EQ(0.0, m(99, 99));
    EXPECT_DOUBLE_EQ(63.3, m(3, 17));
    EXPECT_DOUBLE_EQ(-1.0, m(0, 28));
  }
}

TEST(TlFileMatrix, resize) {
  cleanup();

  // prepare
  TlFileMatrix m(mat_path, 100, 200);
  m.set(30, 40, 3.14);
  EXPECT_EQ(100, m.getNumOfRows());
  EXPECT_EQ(200, m.getNumOfCols());

  // resize
  m.resize(150, 300);
  EXPECT_EQ(150, m.getNumOfRows());
  EXPECT_EQ(300, m.getNumOfCols());
  EXPECT_DOUBLE_EQ(3.14, m.get(30, 40));
}

TEST(TlFileMatrix, getRowVector) {
  cleanup();

  // prepare
  const int row = 1000;
  const int col = 2000;
  TlMatrix m = getMatrix(row, col);
  m.save(mat_path);

  TlFileMatrix fm(mat_path);

  for (int i = 0; i < row; ++i) {
    const TlVector v1 = m.getRowVector(i);
    const TlVector v2 = fm.getRowVector(i);

    EXPECT_EQ(col, v1.getSize());
    EXPECT_EQ(col, v2.getSize());
    for (int j = 0; j < col; ++j) {
      EXPECT_DOUBLE_EQ(v1[j], v2[j]);
    }
  }
}

TEST(TlFileMatrix, getColVector) {
  cleanup();

  // prepare
  const int row = 1000;
  const int col = 2000;
  TlMatrix m = getMatrix(row, col);
  m.save(mat_path);

  TlFileMatrix fm(mat_path);

  for (int i = 0; i < col; ++i) {
    const TlVector v1 = m.getColVector(i);
    const TlVector v2 = fm.getColVector(i);

    EXPECT_EQ(row, v1.getSize());
    EXPECT_EQ(row, v2.getSize());
    for (int j = 0; j < row; ++j) {
      EXPECT_DOUBLE_EQ(v1[j], v2[j]);
    }
  }
}

TEST(TlFileMatrix, setRowVector) {
  cleanup();

  // prepare
  const int row = 1000;
  const int col = 2000;
  TlMatrix m = getMatrix(row, col);

  TlFileMatrix fm(mat_path, row, col);
  for (int i = 0; i < row; ++i) {
    const TlVector v1 = m.getRowVector(i);
    fm.setRowVector(i, v1);
  }

  // check
  EXPECT_EQ(row, fm.getNumOfRows());
  EXPECT_EQ(col, fm.getNumOfCols());
  for (int i = 0; i < row; ++i) {
    for (int j = 0; j < col; ++j) {
      EXPECT_DOUBLE_EQ(m.get(i, j), fm.get(i, j));
    }
  }
}

TEST(TlFileMatrix, setColVector) {
  cleanup();

  // prepare
  const int row = 1000;
  const int col = 2000;
  TlMatrix m = getMatrix(row, col);

  TlFileMatrix fm(mat_path, row, col);
  for (int i = 0; i < col; ++i) {
    const TlVector v1 = m.getColVector(i);
    fm.setColVector(i, v1);
  }

  // check
  EXPECT_EQ(row, fm.getNumOfRows());
  EXPECT_EQ(col, fm.getNumOfCols());
  for (int i = 0; i < row; ++i) {
    for (int j = 0; j < col; ++j) {
      EXPECT_DOUBLE_EQ(m.get(i, j), fm.get(i, j));
    }
  }
}
