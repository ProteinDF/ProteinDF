#include <cstdlib>

#include "TlCommunicate.h"
#include "TlDistributeMatrix.h"
#include "TlMatrix.h"
#include "config.h"
#include "gtest/gtest.h"

static const double EPS = 1.0E-10;
static const std::string mat_path = "temp.mat";
static const std::string distmat_path = "temp.dist.mat";
static const std::string distmat_h5 = "temp.distmat.h5";

static void getMatrix(const int row, const int col, TlMatrix* pRefMat,
                      TlDistributeMatrix* pDistMat) {
  TlCommunicate& rComm = TlCommunicate::getInstance();
  srand((unsigned int)time(NULL));

  pRefMat->resize(row, col);
  if (rComm.isMaster()) {
    for (int i = 0; i < row; ++i) {
      for (int j = 0; j < col; ++j) {
        if (i != j) {
          (*pRefMat)(i, j) = double(rand() / RAND_MAX);
        } else {
          (*pRefMat)(i, i) = double(i);
        }
      }
    }
  }
  rComm.broadcast(*pRefMat);

  pDistMat->resize(row, col);
  for (int i = 0; i < row; ++i) {
    for (int j = 0; j < col; ++j) {
      (*pDistMat)(i, j) = (*pRefMat)(i, j);
    }
  }
}

static void getSymmetricMatrix(const int dim, TlMatrix* pRefMat,
                               TlDistributeMatrix* pDistMat) {
  TlCommunicate& rComm = TlCommunicate::getInstance();
  srand((unsigned int)time(NULL));

  pRefMat->resize(dim, dim);
  if (rComm.isMaster()) {
    for (int i = 0; i < dim; ++i) {
      for (int j = 0; j < i; ++j) {
        const double v = 0.01 * double(i);
        (*pRefMat)(i, j) = v;
        (*pRefMat)(j, i) = v;
      }

      (*pRefMat)(i, i) = double(i);
    }
  }
  rComm.broadcast(*pRefMat);

  pDistMat->resize(dim, dim);
  for (int i = 0; i < dim; ++i) {
    for (int j = 0; j < i; ++j) {
      const double v = (*pRefMat)(i, j);
      (*pDistMat)(i, j) = v;
      (*pDistMat)(j, i) = v;
    }
    (*pDistMat)(i, i) = (*pRefMat)(i, i);
  }
}

TEST(TlDistributeMatrix, constructer) {
  std::size_t row = 1000;
  std::size_t col = 800;
  TlDistributeMatrix A(row, col);
  for (std::size_t i = 0; i < row; ++i) {
    for (std::size_t j = 0; j < col; ++j) {
      EXPECT_DOUBLE_EQ(0.0, A(i, j));
    }
  }

  A(10, 10) = 10.0;
  A(17, 4) = 56.0;
  A(21, 18) = -2.5;
  A(0, 10) = 12.0;

  EXPECT_EQ(row, A.getNumOfRows());
  EXPECT_EQ(col, A.getNumOfCols());
  EXPECT_DOUBLE_EQ(10.0, A(10, 10));
  EXPECT_DOUBLE_EQ(56.0, A(17, 4));
  EXPECT_DOUBLE_EQ(-2.5, A(21, 18));
  EXPECT_DOUBLE_EQ(12.0, A(0, 10));
}

TEST(TlDistributeMatrix, copy_constructer) {
  TlDistributeMatrix A(100, 80);
  A(10, 10) = 10.0;
  A(17, 4) = 56.0;
  A(21, 18) = -2.5;
  A(0, 10) = 12.0;

  const TlDistributeMatrix B = A;
  EXPECT_EQ(100, B.getNumOfRows());
  EXPECT_EQ(80, B.getNumOfCols());
  EXPECT_DOUBLE_EQ(10.0, B(10, 10));
  EXPECT_DOUBLE_EQ(56.0, B(17, 4));
  EXPECT_DOUBLE_EQ(-2.5, B(21, 18));
  EXPECT_DOUBLE_EQ(12.0, B(0, 10));
}

TEST(TlDistributeMatrix, resize) {
  TlDistributeMatrix A(100, 80);
  A(10, 10) = 10.0;
  A(17, 4) = 56.0;
  A(21, 18) = -2.5;
  A(0, 10) = 12.0;

  A.resize(1000, 500);

  EXPECT_EQ(1000, A.getNumOfRows());
  EXPECT_EQ(500, A.getNumOfCols());
  EXPECT_DOUBLE_EQ(10.0, A(10, 10));
  EXPECT_DOUBLE_EQ(56.0, A(17, 4));
  EXPECT_DOUBLE_EQ(-2.5, A(21, 18));
  EXPECT_DOUBLE_EQ(12.0, A(0, 10));
  EXPECT_DOUBLE_EQ(0.0, A(800, 450));
  EXPECT_DOUBLE_EQ(0.0, A(999, 499));
}

TEST(TlDistributeMatrix, trace) {
  TlDistributeMatrix A(100, 80);
  A(10, 10) = 10.0;
  A(75, 75) = 13.0;
  A(17, 4) = 56.0;
  A(21, 18) = -2.5;
  A(0, 10) = 12.0;

  const double trace = A.trace();
  EXPECT_EQ(100, A.getNumOfRows());
  EXPECT_EQ(80, A.getNumOfCols());
  EXPECT_DOUBLE_EQ(23.0, trace);
}

TEST(TlDistributeMatrix, getMaxAbsoluteElement) {
  TlDistributeMatrix A(100, 80);
  A(10, 10) = 10.0;
  A(75, 75) = 13.0;
  A(17, 4) = -56.0;
  A(21, 18) = -2.5;
  A(0, 10) = 12.0;

  TlMatrix::index_type row, col;
  const double v = A.getMaxAbsoluteElement(&row, &col);
  EXPECT_EQ(100, A.getNumOfRows());
  EXPECT_EQ(80, A.getNumOfCols());
  EXPECT_DOUBLE_EQ(56.0, v);
  EXPECT_DOUBLE_EQ(17, row);
  EXPECT_DOUBLE_EQ(4, col);
}

TEST(TlDistributeMatrix, getRMS) {
  TlDistributeMatrix A(100, 80);
  A(10, 10) = 10.0;
  A(75, 75) = 13.0;
  A(17, 4) = -56.0;
  A(21, 18) = -2.5;
  A(0, 10) = 12.0;

  double rms = 0.0;
  rms += 10.0 * 10.0;
  rms += 13.0 * 13.0;
  rms += -56.0 * -56.0;
  rms += -2.5 * -2.5;
  rms += 12.0 * 12.0;
  rms = std::sqrt(rms / (100.0 * 80.0));

  const double v = A.getRMS();
  EXPECT_EQ(100, A.getNumOfRows());
  EXPECT_EQ(80, A.getNumOfCols());
  EXPECT_DOUBLE_EQ(rms, v);
}

TEST(TlDistributeMatrix, operator_eq) {
  TlDistributeMatrix A(100, 120);
  A(10, 10) = 10.0;
  A(17, 4) = 56.0;
  A(21, 18) = -2.5;
  A(0, 10) = 12.0;

  TlDistributeMatrix B;
  B = A;

  EXPECT_EQ(100, B.getNumOfRows());
  EXPECT_EQ(120, B.getNumOfCols());
  EXPECT_DOUBLE_EQ(10.0, B(10, 10));
  EXPECT_DOUBLE_EQ(56.0, B(17, 4));
  EXPECT_DOUBLE_EQ(-2.5, B(21, 18));
  EXPECT_DOUBLE_EQ(12.0, B(0, 10));
  EXPECT_DOUBLE_EQ(0.0, B(80, 90));
}

TEST(TlDistributeMatrix, operator_iadd) {
  TlDistributeMatrix A(100, 120);
  A(10, 10) = 10.0;
  A(17, 4) = 56.0;
  A(21, 18) = -2.5;
  A(0, 10) = 12.0;

  TlDistributeMatrix B(100, 120);
  B(10, 10) = 5.0;
  B(0, 10) = -12.0;
  B(3, 17) = 51.0;

  A += B;

  EXPECT_EQ(100, A.getNumOfRows());
  EXPECT_EQ(120, A.getNumOfCols());
  EXPECT_DOUBLE_EQ(10.0 + 5.0, A(10, 10));
  EXPECT_DOUBLE_EQ(12.0 - 12.0, A(0, 10));
  EXPECT_DOUBLE_EQ(56.0, A(17, 4));
  EXPECT_DOUBLE_EQ(-2.5, A(21, 18));
  EXPECT_DOUBLE_EQ(51.0, A(3, 17));
}

TEST(TlDistributeMatrix, operator_mul_matrix) {
  TlCommunicate& rComm = TlCommunicate::getInstance();

  const int rowA = 100;
  const int colA = 120;
  const int& rowB = colA;
  const int colB = 80;
  TlMatrix refA, refB;
  TlDistributeMatrix A, B;

  getMatrix(rowA, colA, &refA, &A);
  getMatrix(rowB, colB, &refB, &B);

  EXPECT_EQ(rowA, A.getNumOfRows());
  EXPECT_EQ(colA, A.getNumOfCols());
  EXPECT_EQ(rowB, B.getNumOfRows());
  EXPECT_EQ(colB, B.getNumOfCols());

  TlMatrix refC = refA * refB;
  TlDistributeMatrix C = A * B;

  EXPECT_EQ(rowA, C.getNumOfRows());
  EXPECT_EQ(colB, C.getNumOfCols());
  for (int i = 0; i < rowA; ++i) {
    for (int j = 0; j < colB; ++j) {
      EXPECT_NEAR(refC(i, j), C(i, j), EPS);
    }
  }
}

TEST(TlDistributeMatrix, operator_mul_double) {
  TlCommunicate& rComm = TlCommunicate::getInstance();

  const int rowA = 100;
  const int colA = 120;
  TlMatrix refA;
  TlDistributeMatrix A;

  getMatrix(rowA, colA, &refA, &A);

  EXPECT_EQ(rowA, A.getNumOfRows());
  EXPECT_EQ(colA, A.getNumOfCols());

  TlMatrix refC = 3.0 * refA;
  TlDistributeMatrix C = 3.0 * A;

  EXPECT_EQ(rowA, C.getNumOfRows());
  EXPECT_EQ(colA, C.getNumOfCols());
  for (int i = 0; i < rowA; ++i) {
    for (int j = 0; j < colA; ++j) {
      EXPECT_NEAR(refC(i, j), C(i, j), EPS);
    }
  }
}

TEST(TlDistributeMatrix, operator_imul_double) {
  const int rowA = 100;
  const int colA = 120;
  const int& rowB = colA;
  const int colB = 80;
  TlMatrix refA;
  TlDistributeMatrix A;

  getMatrix(rowA, colA, &refA, &A);

  EXPECT_EQ(rowA, A.getNumOfRows());
  EXPECT_EQ(colA, A.getNumOfCols());

  refA *= -2.5;
  A *= -2.5;

  EXPECT_EQ(rowA, A.getNumOfRows());
  EXPECT_EQ(colA, A.getNumOfCols());
  for (int i = 0; i < rowA; ++i) {
    for (int j = 0; j < colA; ++j) {
      EXPECT_NEAR(refA(i, j), A(i, j), EPS);
    }
  }
}

TEST(TlDistributeMatrix, save) {
  const int rowA = 1000;
  const int colA = 1200;
  TlMatrix refA;
  TlDistributeMatrix A;
  getMatrix(rowA, colA, &refA, &A);

  A.save(distmat_path);

  TlMatrix B;
  B.load(distmat_path);

  EXPECT_EQ(rowA, B.getNumOfRows());
  EXPECT_EQ(colA, B.getNumOfCols());
  for (int i = 0; i < rowA; ++i) {
    for (int j = 0; j < colA; ++j) {
      EXPECT_DOUBLE_EQ(refA(i, j), B(i, j));
    }
  }
}

TEST(TlDistributeMatrix, load) {
  const int rowA = 1000;
  const int colA = 1200;
  TlMatrix refA;
  TlDistributeMatrix A;
  getMatrix(rowA, colA, &refA, &A);

  refA.save(mat_path);

  TlDistributeMatrix B;
  B.load(mat_path);

  EXPECT_EQ(rowA, B.getNumOfRows());
  EXPECT_EQ(colA, B.getNumOfCols());
  for (int i = 0; i < rowA; ++i) {
    for (int j = 0; j < colA; ++j) {
      EXPECT_DOUBLE_EQ(refA(i, j), B(i, j));
    }
  }
}

#ifdef HAVE_HDF5
TEST(TlDistributeMatrix, save_hdf5) {
  const int rowA = 100;
  const int colA = 120;
  TlMatrix refA;
  TlDistributeMatrix A;
  getMatrix(rowA, colA, &refA, &A);

  A.saveHdf5(distmat_h5, "save_A");

  TlMatrix B;
  B.loadHdf5(distmat_h5, "save_A");

  EXPECT_EQ(rowA, B.getNumOfRows());
  EXPECT_EQ(colA, B.getNumOfCols());
  for (int i = 0; i < rowA; ++i) {
    for (int j = 0; j < colA; ++j) {
      EXPECT_DOUBLE_EQ(refA(i, j), B(i, j));
    }
  }
}

TEST(TlDistributeMatrix, load_hdf5) {
  const int rowA = 100;
  const int colA = 120;
  TlMatrix refA;
  TlDistributeMatrix A;
  getMatrix(rowA, colA, &refA, &A);

  refA.saveHdf5(distmat_h5, "load_A");

  TlDistributeMatrix B;
  B.loadHdf5(distmat_h5, "load_A");

  EXPECT_EQ(rowA, B.getNumOfRows());
  EXPECT_EQ(colA, B.getNumOfCols());
  for (int i = 0; i < rowA; ++i) {
    for (int j = 0; j < colA; ++j) {
      EXPECT_DOUBLE_EQ(refA(i, j), B(i, j));
    }
  }
}

#endif  // HAVE_HDF5

TEST(TlDistributeMatrix, transpose) {
  const int rowA = 100;
  const int colA = 120;
  TlMatrix refA;
  TlDistributeMatrix A;

  getMatrix(rowA, colA, &refA, &A);

  refA.transpose();
  A.transpose();

  EXPECT_EQ(colA, A.getNumOfRows());
  EXPECT_EQ(rowA, A.getNumOfCols());
  for (int i = 0; i < colA; ++i) {
    for (int j = 0; j < rowA; ++j) {
      EXPECT_DOUBLE_EQ(refA(i, j), A(i, j));
    }
  }
}

TEST(TlDistributeMatrix, inverse) {
  const int dim = 100;
  TlMatrix refA;
  TlDistributeMatrix A;

  getSymmetricMatrix(dim, &refA, &A);

  TlMatrix refB = refA;
  TlDistributeMatrix B = A;

  refB.inverse();
  B.inverse();

  EXPECT_EQ(refB.getNumOfRows(), B.getNumOfRows());
  EXPECT_EQ(refB.getNumOfCols(), B.getNumOfCols());
  for (int i = 0; i < B.getNumOfRows(); ++i) {
    for (int j = 0; j < B.getNumOfCols(); ++j) {
      EXPECT_NEAR(refB(i, j), B(i, j), EPS);
    }
  }
}
