#include <cstdlib>

#include "TlCommunicate.h"
#include "config.h"
#include "gtest/gtest.h"
#include "tl_dense_general_matrix_lapack.h"
#include "tl_dense_general_matrix_scalapack.h"
#include "tl_dense_symmetric_matrix_lapack.h"
#include "tl_dense_symmetric_matrix_scalapack.h"

static const double EPS = 1.0E-5;
static const std::string mat_save_path = "temp.sym.scalapack.save.mat";
static const std::string refmat_save_path = "temp.sym.scalapack.save.ref.mat";
static const std::string mat_load_path = "temp.sym.scalapack.load.mat";
static const std::string distmat_h5 = "temp.dist.sym.h5";

static void getMatrix(const int dim, TlDenseSymmetricMatrix_Lapack* pRefMat,
                      TlDenseSymmetricMatrix_Scalapack* pDistMat) {
  TlCommunicate& rComm = TlCommunicate::getInstance();
  srand((unsigned int)time(NULL));

  pRefMat->resize(dim);
  pDistMat->resize(dim);
  if (rComm.isMaster()) {
    int count = 0;

    const double coef = 1.0 / (dim * dim);
    for (int i = 0; i < dim; ++i) {
      for (int j = 0; j <= i; ++j) {
        // const double v = double(rand() / RAND_MAX);
        const double v = double(count) * coef;
        pRefMat->set(i, j, v);
        // pDistMat->set(i, j, v);
        ++count;
      }
    }
  }
  rComm.broadcast(pRefMat);

  for (int i = 0; i < dim; ++i) {
    for (int j = 0; j <= i; ++j) {
      pDistMat->set(i, j, pRefMat->get(i, j));
    }
  }
}

TEST(TlDenseSymmetricMatrix_Scalapack, constructor) {
  const TlMatrixObject::index_type dim = 1000;
  TlDenseSymmetricMatrix_Scalapack A(dim);
  for (std::size_t i = 0; i < dim; ++i) {
    for (std::size_t j = 0; j < dim; ++j) {
      EXPECT_DOUBLE_EQ(0.0, A.get(i, j));
    }
  }

  A.set(10, 10, 10.0);
  A.set(17, 4, 56.0);
  A.set(21, 18, -2.5);
  A.set(0, 10, 12.0);

  EXPECT_EQ(dim, A.getNumOfRows());
  EXPECT_EQ(dim, A.getNumOfCols());
  EXPECT_DOUBLE_EQ(10.0, A.get(10, 10));
  EXPECT_DOUBLE_EQ(56.0, A.get(17, 4));
  EXPECT_DOUBLE_EQ(56.0, A.get(4, 17));
  EXPECT_DOUBLE_EQ(-2.5, A.get(21, 18));
  EXPECT_DOUBLE_EQ(-2.5, A.get(18, 21));
  EXPECT_DOUBLE_EQ(12.0, A.get(0, 10));
  EXPECT_DOUBLE_EQ(12.0, A.get(10, 0));
}

TEST(TlDenseSymmetricMatrix_Scalapack, copy_constructor) {
  TlDenseSymmetricMatrix_Scalapack A(100);
  A.set(10, 10, 10.0);
  A.set(17, 4, 56.0);
  A.set(21, 18, -2.5);
  A.set(0, 10, 12.0);

  const TlDenseSymmetricMatrix_Scalapack B = A;
  EXPECT_EQ(100, B.getNumOfRows());
  EXPECT_EQ(100, B.getNumOfCols());
  EXPECT_DOUBLE_EQ(10.0, B.get(10, 10));
  EXPECT_DOUBLE_EQ(56.0, B.get(17, 4));
  EXPECT_DOUBLE_EQ(56.0, B.get(4, 17));
  EXPECT_DOUBLE_EQ(-2.5, B.get(21, 18));
  EXPECT_DOUBLE_EQ(-2.5, B.get(18, 21));
  EXPECT_DOUBLE_EQ(12.0, B.get(0, 10));
  EXPECT_DOUBLE_EQ(12.0, B.get(10, 0));
}

TEST(TlDenseSymmetricMatrix_Scalapack, resize) {
  TlDenseSymmetricMatrix_Scalapack A(100);
  A.set(10, 10, 10.0);
  A.set(17, 4, 56.0);
  A.set(21, 18, -2.5);
  A.set(0, 10, 12.0);
  A.set(72, 80, 3.14);
  A.set(99, 9, -1.23);
  A.set(3, 99, 2.34);

  const int newDim = 1000;
  A.resize(newDim);

  EXPECT_EQ(newDim, A.getNumOfRows());
  EXPECT_EQ(newDim, A.getNumOfCols());
  EXPECT_DOUBLE_EQ(10.0, A.get(10, 10));
  EXPECT_DOUBLE_EQ(56.0, A.get(17, 4));
  EXPECT_DOUBLE_EQ(56.0, A.get(4, 17));
  EXPECT_DOUBLE_EQ(-2.5, A.get(21, 18));
  EXPECT_DOUBLE_EQ(-2.5, A.get(18, 21));
  EXPECT_DOUBLE_EQ(12.0, A.get(0, 10));
  EXPECT_DOUBLE_EQ(12.0, A.get(10, 0));
  EXPECT_DOUBLE_EQ(3.14, A.get(72, 80));
  EXPECT_DOUBLE_EQ(3.14, A.get(80, 72));
  EXPECT_DOUBLE_EQ(-1.23, A.get(99, 9));
  EXPECT_DOUBLE_EQ(-1.23, A.get(9, 99));
  EXPECT_DOUBLE_EQ(2.34, A.get(3, 99));
  EXPECT_DOUBLE_EQ(2.34, A.get(99, 3));
  EXPECT_DOUBLE_EQ(0.0, A.get(800, 450));
  EXPECT_DOUBLE_EQ(0.0, A.get(999, 999));
}

TEST(TlDenseSymmetricMatrix_Scalapack, trace) {
  TlDenseSymmetricMatrix_Scalapack A(100);
  A.set(10, 10, 10.0);
  A.set(75, 75, 13.0);
  A.set(17, 4, 56.0);
  A.set(21, 18, -2.5);
  A.set(0, 10, 12.0);

  const double trace = A.trace();
  EXPECT_EQ(100, A.getNumOfRows());
  EXPECT_EQ(100, A.getNumOfCols());
  EXPECT_DOUBLE_EQ(23.0, trace);
}

TEST(TlDenseSymmetricMatrix_Scalapack, getMaxAbsoluteElement) {
  TlDenseSymmetricMatrix_Scalapack A(100);
  A.set(10, 10, 10.0);
  A.set(75, 75, 13.0);
  A.set(17, 4, -56.0);
  A.set(21, 18, -2.5);
  A.set(0, 10, 12.0);

  TlMatrixObject::index_type row, col;
  const double v = A.getMaxAbsoluteElement(&row, &col);
  EXPECT_EQ(100, A.getNumOfRows());
  EXPECT_EQ(100, A.getNumOfCols());
  EXPECT_DOUBLE_EQ(56.0, v);
  EXPECT_DOUBLE_EQ(17, row);
  EXPECT_DOUBLE_EQ(4, col);
}

TEST(TlDenseSymmetricMatrix_Scalapack, getRMS) {
  TlDenseSymmetricMatrix_Scalapack A(100);
  A.set(10, 10, 10.0);
  A.set(75, 75, 13.0);
  A.set(17, 4, -56.0);
  A.set(21, 18, -2.5);
  A.set(0, 10, 12.0);

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

TEST(TlDenseSymmetricMatrix_Scalapack, operator_eq) {
  TlDenseSymmetricMatrix_Scalapack A(100);
  A.set(10, 10, 10.0);
  A.set(17, 4, 56.0);
  A.set(21, 18, -2.5);
  A.set(0, 10, 12.0);

  TlDenseSymmetricMatrix_Scalapack B;
  B = A;

  EXPECT_EQ(100, B.getNumOfRows());
  EXPECT_EQ(100, B.getNumOfCols());
  EXPECT_DOUBLE_EQ(10.0, B.get(10, 10));
  EXPECT_DOUBLE_EQ(56.0, B.get(17, 4));
  EXPECT_DOUBLE_EQ(56.0, B.get(4, 17));
  EXPECT_DOUBLE_EQ(-2.5, B.get(21, 18));
  EXPECT_DOUBLE_EQ(-2.5, B.get(18, 21));
  EXPECT_DOUBLE_EQ(12.0, B.get(0, 10));
  EXPECT_DOUBLE_EQ(12.0, B.get(10, 0));
  EXPECT_DOUBLE_EQ(0.0, B.get(80, 90));
}

TEST(TlDenseSymmetricMatrix_Scalapack, operator_iadd) {
  TlDenseSymmetricMatrix_Scalapack A(100);
  A.set(10, 10, 10.0);
  A.set(17, 4, 56.0);
  A.set(21, 18, -2.5);
  A.set(0, 10, 12.0);

  TlDenseSymmetricMatrix_Scalapack B(100);
  B.set(10, 10, 5.0);
  B.set(0, 10, -12.0);
  B.set(3, 17, 51.0);

  A += B;

  EXPECT_EQ(100, A.getNumOfRows());
  EXPECT_EQ(100, A.getNumOfCols());
  EXPECT_DOUBLE_EQ(10.0 + 5.0, A.get(10, 10));
  EXPECT_DOUBLE_EQ(12.0 - 12.0, A.get(0, 10));
  EXPECT_DOUBLE_EQ(12.0 - 12.0, A.get(10, 0));
  EXPECT_DOUBLE_EQ(56.0, A.get(17, 4));
  EXPECT_DOUBLE_EQ(56.0, A.get(4, 17));
  EXPECT_DOUBLE_EQ(-2.5, A.get(21, 18));
  EXPECT_DOUBLE_EQ(-2.5, A.get(18, 21));
  EXPECT_DOUBLE_EQ(51.0, A.get(3, 17));
  EXPECT_DOUBLE_EQ(51.0, A.get(17, 3));
}

TEST(TlDenseSymmetricMatrix_Scalapack, operator_isub) {
  TlDenseSymmetricMatrix_Scalapack A(100);
  A.set(10, 10, 10.0);
  A.set(17, 4, 56.0);
  A.set(21, 18, -2.5);
  A.set(0, 10, 12.0);

  TlDenseSymmetricMatrix_Scalapack B(100);
  B.set(10, 10, 5.0);
  B.set(0, 10, -12.0);
  B.set(3, 17, 51.0);

  A -= B;

  EXPECT_EQ(100, A.getNumOfRows());
  EXPECT_EQ(100, A.getNumOfCols());
  EXPECT_DOUBLE_EQ(10.0 - (+5.0), A.get(10, 10));
  EXPECT_DOUBLE_EQ(12.0 - (-12.0), A.get(0, 10));
  EXPECT_DOUBLE_EQ(12.0 - (-12.0), A.get(10, 0));
  EXPECT_DOUBLE_EQ(56.0, A.get(17, 4));
  EXPECT_DOUBLE_EQ(56.0, A.get(4, 17));
  EXPECT_DOUBLE_EQ(-2.5, A.get(21, 18));
  EXPECT_DOUBLE_EQ(-2.5, A.get(18, 21));
  EXPECT_DOUBLE_EQ(-51.0, A.get(3, 17));
  EXPECT_DOUBLE_EQ(-51.0, A.get(17, 3));
}

TEST(TlDenseSymmetricMatrix_Scalapack, operator_mul_matrix) {
  const int dim = 1000;
  TlDenseSymmetricMatrix_Lapack refA, refB;
  TlDenseSymmetricMatrix_Scalapack A, B;

  getMatrix(dim, &refA, &A);
  getMatrix(dim, &refB, &B);

  EXPECT_EQ(dim, A.getNumOfRows());
  EXPECT_EQ(dim, A.getNumOfCols());
  for (int i = 0; i < dim; ++i) {
    for (int j = 0; j < dim; ++j) {
      EXPECT_NEAR(refA.get(i, j), A.get(i, j), EPS);
    }
  }

  EXPECT_EQ(dim, A.getNumOfRows());
  EXPECT_EQ(dim, A.getNumOfCols());
  EXPECT_EQ(dim, B.getNumOfRows());
  EXPECT_EQ(dim, B.getNumOfCols());

  TlDenseGeneralMatrix_Lapack refC = refA * refB;
  TlDenseGeneralMatrix_Scalapack C = A * B;

  EXPECT_EQ(dim, C.getNumOfRows());
  EXPECT_EQ(dim, C.getNumOfCols());
  for (int i = 0; i < dim; ++i) {
    for (int j = 0; j < dim; ++j) {
      EXPECT_NEAR(refC.get(i, j), C.get(i, j), EPS);
    }
  }
}

TEST(TlDenseSymmetricMatrix_Scalapack, operator_mul_double) {
  const int dim = 1000;
  TlDenseSymmetricMatrix_Lapack refA;
  TlDenseSymmetricMatrix_Scalapack A;

  getMatrix(dim, &refA, &A);

  EXPECT_EQ(dim, A.getNumOfRows());
  EXPECT_EQ(dim, A.getNumOfCols());

  TlDenseSymmetricMatrix_Lapack refC = 3.0 * refA;
  TlDenseSymmetricMatrix_Scalapack C = 3.0 * A;

  EXPECT_EQ(dim, C.getNumOfRows());
  EXPECT_EQ(dim, C.getNumOfCols());
  for (int i = 0; i < dim; ++i) {
    for (int j = 0; j < dim; ++j) {
      EXPECT_NEAR(refC.get(i, j), C.get(i, j), EPS);
    }
  }
}

TEST(TlDenseSymmetricMatrix_Scalapack, operator_imul_double) {
  const int dim = 1000;
  TlDenseSymmetricMatrix_Lapack refA;
  TlDenseSymmetricMatrix_Scalapack A;

  getMatrix(dim, &refA, &A);

  EXPECT_EQ(dim, A.getNumOfRows());
  EXPECT_EQ(dim, A.getNumOfCols());

  refA *= -2.5;
  A *= -2.5;

  EXPECT_EQ(dim, A.getNumOfRows());
  EXPECT_EQ(dim, A.getNumOfCols());
  for (int i = 0; i < dim; ++i) {
    for (int j = 0; j < dim; ++j) {
      EXPECT_NEAR(refA.get(i, j), A.get(i, j), EPS);
    }
  }
}

TEST(TlDenseSymmetricMatrix_Scalapack, save) {
  TlCommunicate& rComm = TlCommunicate::getInstance();

  const int dim = 1000;
  TlDenseSymmetricMatrix_Lapack refA;
  TlDenseSymmetricMatrix_Scalapack A;

  getMatrix(dim, &refA, &A);
  if (rComm.isMaster()) {
      refA.save(refmat_save_path);
  }
  A.save(mat_save_path);

  if (rComm.isMaster()) {
    TlDenseSymmetricMatrix_Lapack B;
    B.load(mat_save_path);

    EXPECT_EQ(dim, B.getNumOfRows());
    EXPECT_EQ(dim, B.getNumOfCols());
    for (int i = 0; i < dim; ++i) {
      for (int j = 0; j < dim; ++j) {
        EXPECT_DOUBLE_EQ(refA.get(i, j), B.get(i, j));
      }
    }
  }
}

TEST(TlDenseSymmetricMatrix_Scalapack, load) {
  TlCommunicate& rComm = TlCommunicate::getInstance();

  const int dim = 1000;
  TlDenseSymmetricMatrix_Lapack refA;
  TlDenseSymmetricMatrix_Scalapack A;

  getMatrix(dim, &refA, &A);
  if (rComm.isMaster()) {
    refA.save(mat_load_path);
  }
  rComm.barrier();

  TlDenseSymmetricMatrix_Scalapack B;
  B.load(mat_load_path);

  EXPECT_EQ(dim, B.getNumOfRows());
  EXPECT_EQ(dim, B.getNumOfCols());
  for (int i = 0; i < dim; ++i) {
    for (int j = 0; j < dim; ++j) {
      const double v = B.get(i, j);
      if (rComm.isMaster()) {
        EXPECT_DOUBLE_EQ(refA.get(i, j), v);
      }
    }
  }
}

// #ifdef HAVE_HDF5
// TEST(TlDenseSymmetricMatrix_Scalapack, save_hdf5) {
//   TlCommunicate& rComm = TlCommunicate::getInstance();
//
//   const int dim = 1000;
//   TlDenseSymmetricMatrix_Lapack refA;
//   TlDenseSymmetricMatrix_Scalapack A;
//   getMatrix(dim, &refA, &A);
//
//   A.saveHdf5(distmat_h5, "save_A");
//
//   if (rComm.isMaster()) {
//     TlDenseSymmetricMatrix_Lapack B;
//     B.loadHdf5(distmat_h5, "save_A");
//     // B.save("temp.B.mat");
//
//     EXPECT_EQ(dim, B.getNumOfRows());
//     EXPECT_EQ(dim, B.getNumOfCols());
//     for (int i = 0; i < dim; ++i) {
//       for (int j = 0; j < dim; ++j) {
//         EXPECT_DOUBLE_EQ(refA.get(i, j), B.get(i, j));
//       }
//     }
//   }
// }
//
// TEST(TlDenseSymmetricMatrix_Scalapack, load_hdf5) {
//   TlCommunicate& rComm = TlCommunicate::getInstance();
//
//   const int dim = 1000;
//   TlDenseSymmetricMatrix_Lapack refA;
//   TlDenseSymmetricMatrix_Scalapack A;
//   getMatrix(dim, &refA, &A);
//
//   if (rComm.isMaster()) {
//     refA.saveHdf5(distmat_h5, "load_A");
//   }
//
//   TlDenseSymmetricMatrix_Scalapack B;
//   B.loadHdf5(distmat_h5, "load_A");
//
//   EXPECT_EQ(dim, B.getNumOfRows());
//   EXPECT_EQ(dim, B.getNumOfCols());
//   for (int i = 0; i < dim; ++i) {
//     for (int j = 0; j < dim; ++j) {
//       EXPECT_DOUBLE_EQ(refA.get(i, j), B.get(i, j));
//     }
//   }
// }
//
// #endif  // HAVE_HDF5

// TEST(TlDenseSymmetricMatrix_Scalapack, inverse)
// {
//     const int dim = 1000;
//     TlDenseSymmetricMatrix_Lapack refA;
//     TlDenseSymmetricMatrix_Scalapack A;
//
//     getMatrix(dim, &refA, &A);
//
//     TlDenseSymmetricMatrix_Lapack refB = refA;
//     TlDenseSymmetricMatrix_Scalapack B = A;
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
