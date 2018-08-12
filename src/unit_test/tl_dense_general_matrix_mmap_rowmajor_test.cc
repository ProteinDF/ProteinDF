#include "gtest/gtest.h"

#include "TlFile.h"
#include "tl_dense_general_matrix_lapack.h"
#include "tl_matrix_mmap_row_major.h"

static const std::string matPath = "temp.mmap.rsfd.mat";

static void cleanup() {
  if (TlFile::isExistFile(matPath)) {
    TlFile::remove(matPath);
  }
}

TEST(TlMatrixMmapRowMajor, construct) {
  const int row = 100;
  const int col = 200;
  {
    cleanup();
    TlMatrixMmapRowMajor m(matPath, row, col);
  }

  {
    TlDenseGeneralMatrix_Lapack a;
    a.load(matPath);

    EXPECT_EQ(row, a.getNumOfRows());
    EXPECT_EQ(col, a.getNumOfCols());
  }
}

// TEST(TlMatrixMmapRowMajor, resize) {
//   cleanup();
//   const int row = 100;
//   const int col = 200;
//   TlDenseGeneralMatrix_Lapack_rsfd ref = getTlMatrix(row, col);
//   ref.save(matPath);
//   TlMatrixMmapRowMajor a(matPath);
//
//   const int newRow = 200;
//   const int newCol = 200;
//   ref.resize(newRow, newCol);
//   a.resize(newRow, newCol);
//
//   EXPECT_EQ(newRow, a.getNumOfRows());
//   EXPECT_EQ(newCol, a.getNumOfCols());
//   for (int i = 0; i < newRow; ++i) {
//     for (int j = 0; j < newCol; ++j) {
//       EXPECT_DOUBLE_EQ(ref.get(i, j), a.get(i, j));
//     }
//   }
// }
//
// TEST(TlMatrixMmapRowMajor, set_matrix) {
//   cleanup();
//   const int row = 100;
//   const int col = 200;
//   TlDenseGeneralMatrix_Lapack ref = getTlMatrix(row, col);
//
//   {
//     TlMatrixMmapRowMajor m(matPath, ref.getNumOfRows(), ref.getNumOfCols());
//     for (int r = 0; r < ref.getNumOfRows(); ++r) {
//       for (int c = 0; c < ref.getNumOfCols(); ++c) {
//         m.set(r, c, ref.get(r, c));
//       }
//     }
//   }
//
//   {
//     TlDenseGeneralMatrix_Lapack a;
//     a.load(matPath);
//
//     EXPECT_EQ(ref.getNumOfRows(), a.getNumOfRows());
//     EXPECT_EQ(ref.getNumOfCols(), a.getNumOfCols());
//     for (int r = 0; r < ref.getNumOfRows(); ++r) {
//       for (int c = 0; c < ref.getNumOfCols(); ++c) {
//         EXPECT_DOUBLE_EQ(ref.get(r, c), a.get(r, c));
//       }
//     }
//   }
// }
//
// TEST(TlMatrixMmapRowMajor, get_matrix) {
//   cleanup();
//   const int row = 200;
//   const int col = 100;
//   TlDenseGeneralMatrix_Lapack ref = getTlMatrix(row, col);
//   ref.save(matPath);
//
//   {
//     TlMatrixMmapRowMajor m(matPath);
//     EXPECT_EQ(ref.getNumOfRows(), m.getNumOfRows());
//     EXPECT_EQ(ref.getNumOfCols(), m.getNumOfCols());
//     for (int r = 0; r < ref.getNumOfRows(); ++r) {
//       for (int c = 0; c < ref.getNumOfCols(); ++c) {
//         EXPECT_DOUBLE_EQ(ref.get(r, c), m.get(r, c));
//       }
//     }
//   }
// }
