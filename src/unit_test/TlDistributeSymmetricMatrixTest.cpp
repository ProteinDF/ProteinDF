#include <cstdlib>

#include "TlCommunicate.h"
#include "TlDistributeSymmetricMatrix.h"
#include "TlSymmetricMatrix.h"
#include "config.h"
#include "gtest/gtest.h"

static const double EPS = 1.0E-5;
static const std::string mat_path = "temp.sym.mat";
static const std::string distmat_path = "temp.dist.sym.mat";
static const std::string distmat_h5 = "temp.dist.sym.h5";

static void getMatrix(const int dim, TlSymmetricMatrix* pRefMat,
                      TlDistributeSymmetricMatrix* pDistMat) {
  TlCommunicate& rComm = TlCommunicate::getInstance();
  srand((unsigned int)time(NULL));

  pRefMat->resize(dim);
  if (rComm.isMaster()) {
    int count = 0;
    for (int i = 0; i < dim; ++i) {
      for (int j = 0; j <= i; ++j) {
        // (*pRefMat)(i, j) = double(rand() / RAND_MAX);
        (*pRefMat)(i, j) = double(count);
        ++count;
      }
    }
  }
  rComm.broadcast(*pRefMat);

  pDistMat->resize(dim);
  for (int i = 0; i < dim; ++i) {
    for (int j = 0; j <= i; ++j) {
      (*pDistMat)(i, j) = (*pRefMat)(i, j);
    }
  }
}

TEST(TlDistributeSymmetricMatrix, constructer) {
  const std::size_t dim = 1000;
  TlDistributeSymmetricMatrix A(dim);
  for (std::size_t i = 0; i < dim; ++i) {
    for (std::size_t j = 0; j < dim; ++j) {
      EXPECT_DOUBLE_EQ(0.0, A(i, j));
    }
  }

  A(10, 10) = 10.0;
  A(17, 4) = 56.0;
  A(21, 18) = -2.5;
  A(0, 10) = 12.0;

  EXPECT_EQ(dim, A.getNumOfRows());
  EXPECT_EQ(dim, A.getNumOfCols());
  EXPECT_DOUBLE_EQ(10.0, A(10, 10));
  EXPECT_DOUBLE_EQ(56.0, A(17, 4));
  EXPECT_DOUBLE_EQ(56.0, A(4, 17));
  EXPECT_DOUBLE_EQ(-2.5, A(21, 18));
  EXPECT_DOUBLE_EQ(-2.5, A(18, 21));
  EXPECT_DOUBLE_EQ(12.0, A(0, 10));
  EXPECT_DOUBLE_EQ(12.0, A(10, 0));
}

TEST(TlDistributeSymmetricMatrix, copy_constructer) {
  TlDistributeSymmetricMatrix A(100);
  A(10, 10) = 10.0;
  A(17, 4) = 56.0;
  A(21, 18) = -2.5;
  A(0, 10) = 12.0;

  const TlDistributeMatrix B = A;
  EXPECT_EQ(100, B.getNumOfRows());
  EXPECT_EQ(100, B.getNumOfCols());
  EXPECT_DOUBLE_EQ(10.0, B(10, 10));
  EXPECT_DOUBLE_EQ(56.0, B(17, 4));
  EXPECT_DOUBLE_EQ(56.0, B(4, 17));
  EXPECT_DOUBLE_EQ(-2.5, B(21, 18));
  EXPECT_DOUBLE_EQ(-2.5, B(18, 21));
  EXPECT_DOUBLE_EQ(12.0, B(0, 10));
  EXPECT_DOUBLE_EQ(12.0, B(10, 0));
}

TEST(TlDistributeSymmetricMatrix, resize) {
  TlDistributeSymmetricMatrix A(100);
  A(10, 10) = 10.0;
  A(17, 4) = 56.0;
  A(21, 18) = -2.5;
  A(0, 10) = 12.0;
  A(72, 80) = 3.14;
  A(99, 9) = -1.23;
  A(3, 99) = 2.34;

  const int newDim = 1000;
  A.resize(newDim);

  EXPECT_EQ(newDim, A.getNumOfRows());
  EXPECT_EQ(newDim, A.getNumOfCols());
  EXPECT_DOUBLE_EQ(10.0, A(10, 10));
  EXPECT_DOUBLE_EQ(56.0, A(17, 4));
  EXPECT_DOUBLE_EQ(56.0, A(4, 17));
  EXPECT_DOUBLE_EQ(-2.5, A(21, 18));
  EXPECT_DOUBLE_EQ(-2.5, A(18, 21));
  EXPECT_DOUBLE_EQ(12.0, A(0, 10));
  EXPECT_DOUBLE_EQ(12.0, A(10, 0));
  EXPECT_DOUBLE_EQ(3.14, A(72, 80));
  EXPECT_DOUBLE_EQ(3.14, A(80, 72));
  EXPECT_DOUBLE_EQ(-1.23, A(99, 9));
  EXPECT_DOUBLE_EQ(-1.23, A(9, 99));
  EXPECT_DOUBLE_EQ(2.34, A(3, 99));
  EXPECT_DOUBLE_EQ(2.34, A(99, 3));
  EXPECT_DOUBLE_EQ(0.0, A(800, 450));
  EXPECT_DOUBLE_EQ(0.0, A(999, 999));
}

TEST(TlDistributeSymmetricMatrix, trace) {
  TlDistributeSymmetricMatrix A(100);
  A(10, 10) = 10.0;
  A(75, 75) = 13.0;
  A(17, 4) = 56.0;
  A(21, 18) = -2.5;
  A(0, 10) = 12.0;

  const double trace = A.trace();
  EXPECT_EQ(100, A.getNumOfRows());
  EXPECT_EQ(100, A.getNumOfCols());
  EXPECT_DOUBLE_EQ(23.0, trace);
}

TEST(TlDistributeSymmetricMatrix, getMaxAbsoluteElement) {
  TlDistributeSymmetricMatrix A(100);
  A(10, 10) = 10.0;
  A(75, 75) = 13.0;
  A(17, 4) = -56.0;
  A(21, 18) = -2.5;
  A(0, 10) = 12.0;

  TlMatrix::index_type row, col;
  const double v = A.getMaxAbsoluteElement(&row, &col);
  EXPECT_EQ(100, A.getNumOfRows());
  EXPECT_EQ(100, A.getNumOfCols());
  EXPECT_DOUBLE_EQ(56.0, v);
  EXPECT_DOUBLE_EQ(17, row);
  EXPECT_DOUBLE_EQ(4, col);
}

TEST(TlDistributeSymmetricMatrix, getRMS) {
  TlDistributeSymmetricMatrix A(100);
  A(10, 10) = 10.0;
  A(75, 75) = 13.0;
  A(17, 4) = -56.0;
  A(21, 18) = -2.5;
  A(0, 10) = 12.0;

  double rms = 0.0;
  rms += 10.0 * 10.0;    // (10, 10)
  rms += 13.0 * 13.0;    // (75, 75)
  rms += -56.0 * -56.0;  // (17,  4)
  rms += -56.0 * -56.0;  // ( 4, 17)
  rms += -2.5 * -2.5;    // (21, 18)
  rms += -2.5 * -2.5;    // (18, 21)
  rms += 12.0 * 12.0;    // ( 0, 10)
  rms += 12.0 * 12.0;    // (10,  0)
  rms = std::sqrt(rms / (100.0 * 100.0));

  const double v = A.getRMS();
  EXPECT_EQ(100, A.getNumOfRows());
  EXPECT_EQ(100, A.getNumOfCols());
  EXPECT_DOUBLE_EQ(rms, v);
}

TEST(TlDistributeSymmetricMatrix, operator_eq) {
  TlDistributeSymmetricMatrix A(100);
  A(10, 10) = 10.0;
  A(17, 4) = 56.0;
  A(21, 18) = -2.5;
  A(0, 10) = 12.0;

  TlDistributeSymmetricMatrix B;
  B = A;

  EXPECT_EQ(100, B.getNumOfRows());
  EXPECT_EQ(100, B.getNumOfCols());
  EXPECT_DOUBLE_EQ(10.0, B(10, 10));
  EXPECT_DOUBLE_EQ(56.0, B(17, 4));
  EXPECT_DOUBLE_EQ(56.0, B(4, 17));
  EXPECT_DOUBLE_EQ(-2.5, B(21, 18));
  EXPECT_DOUBLE_EQ(-2.5, B(18, 21));
  EXPECT_DOUBLE_EQ(12.0, B(0, 10));
  EXPECT_DOUBLE_EQ(12.0, B(10, 0));
  EXPECT_DOUBLE_EQ(0.0, B(80, 90));
}

TEST(TlDistributeSymmetricMatrix, operator_iadd) {
  TlDistributeSymmetricMatrix A(100);
  A(10, 10) = 10.0;
  A(17, 4) = 56.0;
  A(21, 18) = -2.5;
  A(0, 10) = 12.0;

  TlDistributeSymmetricMatrix B(100);
  B(10, 10) = 5.0;
  B(0, 10) = -12.0;
  B(3, 17) = 51.0;

  A += B;

  EXPECT_EQ(100, A.getNumOfRows());
  EXPECT_EQ(100, A.getNumOfCols());
  EXPECT_DOUBLE_EQ(10.0 + 5.0, A(10, 10));
  EXPECT_DOUBLE_EQ(12.0 - 12.0, A(0, 10));
  EXPECT_DOUBLE_EQ(12.0 - 12.0, A(10, 0));
  EXPECT_DOUBLE_EQ(56.0, A(17, 4));
  EXPECT_DOUBLE_EQ(56.0, A(4, 17));
  EXPECT_DOUBLE_EQ(-2.5, A(21, 18));
  EXPECT_DOUBLE_EQ(-2.5, A(18, 21));
  EXPECT_DOUBLE_EQ(51.0, A(3, 17));
  EXPECT_DOUBLE_EQ(51.0, A(17, 3));
}

TEST(TlDistributeSymmetricMatrix, operator_mul_matrix) {
  TlCommunicate& rComm = TlCommunicate::getInstance();

  const int dim = 1000;
  TlSymmetricMatrix refA, refB;
  TlDistributeSymmetricMatrix A, B;

  getMatrix(dim, &refA, &A);
  getMatrix(dim, &refB, &B);

  EXPECT_EQ(dim, A.getNumOfRows());
  EXPECT_EQ(dim, A.getNumOfCols());
  EXPECT_EQ(dim, B.getNumOfRows());
  EXPECT_EQ(dim, B.getNumOfCols());

  TlMatrix refC = refA * refB;
  TlDistributeMatrix C = A * B;

  EXPECT_EQ(dim, C.getNumOfRows());
  EXPECT_EQ(dim, C.getNumOfCols());
  for (int i = 0; i < dim; ++i) {
    for (int j = 0; j < dim; ++j) {
      EXPECT_NEAR(refC(i, j), C(i, j), EPS);
    }
  }
}

TEST(TlDistributeSymmetricMatrix, operator_mul_double) {
  TlCommunicate& rComm = TlCommunicate::getInstance();

  const int dim = 1000;
  TlSymmetricMatrix refA;
  TlDistributeSymmetricMatrix A;

  getMatrix(dim, &refA, &A);

  EXPECT_EQ(dim, A.getNumOfRows());
  EXPECT_EQ(dim, A.getNumOfCols());

  TlSymmetricMatrix refC = 3.0 * refA;
  TlDistributeSymmetricMatrix C = 3.0 * A;

  EXPECT_EQ(dim, C.getNumOfRows());
  EXPECT_EQ(dim, C.getNumOfCols());
  for (int i = 0; i < dim; ++i) {
    for (int j = 0; j < dim; ++j) {
      EXPECT_NEAR(refC(i, j), C(i, j), EPS);
    }
  }
}

TEST(TlDistributeSymmetricMatrix, operator_imul_double) {
  const int dim = 1000;
  TlSymmetricMatrix refA;
  TlDistributeSymmetricMatrix A;

  getMatrix(dim, &refA, &A);

  EXPECT_EQ(dim, A.getNumOfRows());
  EXPECT_EQ(dim, A.getNumOfCols());

  refA *= -2.5;
  A *= -2.5;

  EXPECT_EQ(dim, A.getNumOfRows());
  EXPECT_EQ(dim, A.getNumOfCols());
  for (int i = 0; i < dim; ++i) {
    for (int j = 0; j < dim; ++j) {
      EXPECT_NEAR(refA(i, j), A(i, j), EPS);
    }
  }
}

TEST(TlDistributeSymmetricMatrix, save) {
  const int dim = 1000;
  TlSymmetricMatrix refA;
  TlDistributeSymmetricMatrix A;

  getMatrix(dim, &refA, &A);

  A.save(distmat_path);

  TlSymmetricMatrix B;
  B.load(distmat_path);

  EXPECT_EQ(dim, B.getNumOfRows());
  EXPECT_EQ(dim, B.getNumOfCols());
  for (int i = 0; i < dim; ++i) {
    for (int j = 0; j < dim; ++j) {
      EXPECT_DOUBLE_EQ(refA(i, j), B(i, j));
    }
  }
}

TEST(TlDistributeSymmetricMatrix, load) {
  const int dim = 1000;
  TlSymmetricMatrix refA;
  TlDistributeSymmetricMatrix A;

  getMatrix(dim, &refA, &A);

  refA.save(mat_path);

  TlDistributeSymmetricMatrix B;
  B.load(mat_path);

  EXPECT_EQ(dim, B.getNumOfRows());
  EXPECT_EQ(dim, B.getNumOfCols());
  for (int i = 0; i < dim; ++i) {
    for (int j = 0; j < dim; ++j) {
      EXPECT_DOUBLE_EQ(refA(i, j), B(i, j));
    }
  }
}

#ifdef HAVE_HDF5
TEST(TlDistributeSymmetricMatrix, save_hdf5) {
  TlCommunicate& rComm = TlCommunicate::getInstance();

  const int dim = 1000;
  TlSymmetricMatrix refA;
  TlDistributeSymmetricMatrix A;
  getMatrix(dim, &refA, &A);

  A.saveHdf5(distmat_h5, "save_A");

  if (rComm.isMaster()) {
    TlSymmetricMatrix B;
    B.loadHdf5(distmat_h5, "save_A");
    // B.save("temp.B.mat");

    EXPECT_EQ(dim, B.getNumOfRows());
    EXPECT_EQ(dim, B.getNumOfCols());
    for (int i = 0; i < dim; ++i) {
      for (int j = 0; j < dim; ++j) {
        EXPECT_DOUBLE_EQ(refA(i, j), B(i, j));
      }
    }
  }
}

TEST(TlDistributeSymmetricMatrix, load_hdf5) {
  TlCommunicate& rComm = TlCommunicate::getInstance();

  const int dim = 1000;
  TlSymmetricMatrix refA;
  TlDistributeSymmetricMatrix A;
  getMatrix(dim, &refA, &A);

  if (rComm.isMaster()) {
    refA.saveHdf5(distmat_h5, "load_A");
  }

  TlDistributeSymmetricMatrix B;
  B.loadHdf5(distmat_h5, "load_A");

  EXPECT_EQ(dim, B.getNumOfRows());
  EXPECT_EQ(dim, B.getNumOfCols());
  for (int i = 0; i < dim; ++i) {
    for (int j = 0; j < dim; ++j) {
      EXPECT_DOUBLE_EQ(refA(i, j), B(i, j));
    }
  }
}

#endif  // HAVE_HDF5

// TEST(TlDistributeSymmetricMatrix, inverse)
// {
//     const int dim = 1000;
//     TlSymmetricMatrix refA;
//     TlDistributeSymmetricMatrix A;
//
//     getMatrix(dim, &refA, &A);
//
//     TlSymmetricMatrix refB = refA;
//     TlDistributeSymmetricMatrix B = A;
//     // EXPECT_EQ(refB.getNumOfRows(), B.getNumOfRows());
//     // EXPECT_EQ(refB.getNumOfCols(), B.getNumOfCols());
//     // for (int i = 0; i < B.getNumOfRows(); ++i) {
//     //     for (int j = 0; j < B.getNumOfCols(); ++j) {
//     //         EXPECT_NEAR(refB(i, j), B(i, j), EPS);
//     //     }
//     // }
//
//     refB.inverse();
//     B.inverse();
//
//     EXPECT_EQ(refB.getNumOfRows(), B.getNumOfRows());
//     EXPECT_EQ(refB.getNumOfCols(), B.getNumOfCols());
//     for (int i = 0; i < B.getNumOfRows(); ++i) {
//         for (int j = 0; j < B.getNumOfCols(); ++j) {
//             EXPECT_NEAR(refB(i, j), B(i, j), EPS);
//         }
//     }
// }
