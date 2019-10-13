#include <cstdlib>

#include "TlCommunicate.h"
#include "config.h"
#include "gtest/gtest.h"
#include "tl_dense_general_matrix_blacs.h"
#include "tl_dense_general_matrix_blas_old.h"

static const double EPS = 1.0E-10;
static const std::string mat_path = "temp.mat";
static const std::string distmat_path = "temp.dist.mat";
static const std::string distmat_h5 = "temp.distmat.h5";

static void getMatrix(const int row, const int col,
                      TlDenseGeneralMatrix_BLAS_old* pRefMat,
                      TlDenseGeneralMatrix_blacs* pDistMat) {
    TlCommunicate& rComm = TlCommunicate::getInstance();
    srand((unsigned int)time(NULL));

    pRefMat->resize(row, col);
    if (rComm.isMaster()) {
        for (int i = 0; i < row; ++i) {
            for (int j = 0; j < col; ++j) {
                if (i != j) {
                    pRefMat->set(i, j, double(rand() / RAND_MAX));
                } else {
                    pRefMat->set(i, i, double(i));
                }
            }
        }
    }
    rComm.broadcast(*pRefMat);

    pDistMat->resize(row, col);
    for (int i = 0; i < row; ++i) {
        for (int j = 0; j < col; ++j) {
            pDistMat->set(i, j, pRefMat->get(i, j));
        }
    }
}

static void getSymmetricMatrix(const int dim,
                               TlDenseGeneralMatrix_BLAS_old* pRefMat,
                               TlDenseGeneralMatrix_blacs* pDistMat) {
    TlCommunicate& rComm = TlCommunicate::getInstance();
    srand((unsigned int)time(NULL));

    pRefMat->resize(dim, dim);
    if (rComm.isMaster()) {
        for (int i = 0; i < dim; ++i) {
            for (int j = 0; j < i; ++j) {
                const double v = 0.01 * double(i);
                pRefMat->set(i, j, v);
                pRefMat->set(j, i, v);
            }

            pRefMat->set(i, i, double(i));
        }
    }
    rComm.broadcast(*pRefMat);

    pDistMat->resize(dim, dim);
    for (int i = 0; i < dim; ++i) {
        for (int j = 0; j < i; ++j) {
            const double v = pRefMat->get(i, j);
            pDistMat->set(i, j, v);
            pDistMat->set(j, i, v);
        }
        pDistMat->set(i, i, pRefMat->get(i, i));
    }
}

TEST(TlDenseGeneralMatrix_blacs, constructer) {
    std::size_t row = 1000;
    std::size_t col = 800;
    TlDenseGeneralMatrix_blacs A(row, col);
    for (std::size_t i = 0; i < row; ++i) {
        for (std::size_t j = 0; j < col; ++j) {
            EXPECT_DOUBLE_EQ(0.0, A.get(i, j));
        }
    }

    A.set(10, 10, 10.0);
    A.set(17, 4, 56.0);
    A.set(21, 18, -2.5);
    A.set(0, 10, 12.0);

    EXPECT_EQ(row, A.getNumOfRows());
    EXPECT_EQ(col, A.getNumOfCols());
    EXPECT_DOUBLE_EQ(10.0, A.get(10, 10));
    EXPECT_DOUBLE_EQ(56.0, A.get(17, 4));
    EXPECT_DOUBLE_EQ(-2.5, A.get(21, 18));
    EXPECT_DOUBLE_EQ(12.0, A.get(0, 10));
}

TEST(TlDenseGeneralMatrix_blacs, copy_constructer) {
    TlDenseGeneralMatrix_blacs A(100, 80);
    A.set(10, 10, 10.0);
    A.set(17, 4, 56.0);
    A.set(21, 18, -2.5);
    A.set(0, 10, 12.0);

    const TlDenseGeneralMatrix_blacs B = A;
    EXPECT_EQ(100, B.getNumOfRows());
    EXPECT_EQ(80, B.getNumOfCols());
    EXPECT_DOUBLE_EQ(10.0, B.get(10, 10));
    EXPECT_DOUBLE_EQ(56.0, B.get(17, 4));
    EXPECT_DOUBLE_EQ(-2.5, B.get(21, 18));
    EXPECT_DOUBLE_EQ(12.0, B.get(0, 10));
}

TEST(TlDenseGeneralMatrix_blacs, resize) {
    TlDenseGeneralMatrix_blacs A(100, 80);
    A.set(10, 10, 10.0);
    A.set(17, 4, 56.0);
    A.set(21, 18, -2.5);
    A.set(0, 10, 12.0);

    A.resize(1000, 500);

    EXPECT_EQ(1000, A.getNumOfRows());
    EXPECT_EQ(500, A.getNumOfCols());
    EXPECT_DOUBLE_EQ(10.0, A.get(10, 10));
    EXPECT_DOUBLE_EQ(56.0, A.get(17, 4));
    EXPECT_DOUBLE_EQ(-2.5, A.get(21, 18));
    EXPECT_DOUBLE_EQ(12.0, A.get(0, 10));
    EXPECT_DOUBLE_EQ(0.0, A.get(800, 450));
    EXPECT_DOUBLE_EQ(0.0, A.get(999, 499));
}

TEST(TlDenseGeneralMatrix_blacs, trace) {
    TlDenseGeneralMatrix_blacs A(100, 80);
    A.set(10, 10, 10.0);
    A.set(75, 75, 13.0);
    A.set(17, 4, 56.0);
    A.set(21, 18, -2.5);
    A.set(0, 10, 12.0);

    const double trace = A.trace();
    EXPECT_EQ(100, A.getNumOfRows());
    EXPECT_EQ(80, A.getNumOfCols());
    EXPECT_DOUBLE_EQ(23.0, trace);
}

TEST(TlDenseGeneralMatrix_blacs, getMaxAbsoluteElement) {
    TlDenseGeneralMatrix_blacs A(100, 80);
    A.set(10, 10, 10.0);
    A.set(75, 75, 13.0);
    A.set(17, 4, -56.0);
    A.set(21, 18, -2.5);
    A.set(0, 10, 12.0);

    TlDenseGeneralMatrix_BLAS_old::index_type row, col;
    const double v = A.getMaxAbsoluteElement(&row, &col);
    EXPECT_EQ(100, A.getNumOfRows());
    EXPECT_EQ(80, A.getNumOfCols());
    EXPECT_DOUBLE_EQ(56.0, v);
    EXPECT_DOUBLE_EQ(17, row);
    EXPECT_DOUBLE_EQ(4, col);
}

TEST(TlDenseGeneralMatrix_blacs, getRMS) {
    TlDenseGeneralMatrix_blacs A(100, 80);
    A.set(10, 10, 10.0);
    A.set(75, 75, 13.0);
    A.set(17, 4, -56.0);
    A.set(21, 18, -2.5);
    A.set(0, 10, 12.0);

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

TEST(TlDenseGeneralMatrix_blacs, operator_eq) {
    TlDenseGeneralMatrix_blacs A(100, 120);
    A.set(10, 10, 10.0);
    A.set(17, 4, 56.0);
    A.set(21, 18, -2.5);
    A.set(0, 10, 12.0);

    TlDenseGeneralMatrix_blacs B;
    B = A;

    EXPECT_EQ(100, B.getNumOfRows());
    EXPECT_EQ(120, B.getNumOfCols());
    EXPECT_DOUBLE_EQ(10.0, B.get(10, 10));
    EXPECT_DOUBLE_EQ(56.0, B.get(17, 4));
    EXPECT_DOUBLE_EQ(-2.5, B.get(21, 18));
    EXPECT_DOUBLE_EQ(12.0, B.get(0, 10));
    EXPECT_DOUBLE_EQ(0.0, B.get(80, 90));
}

TEST(TlDenseGeneralMatrix_blacs, operator_iadd) {
    TlDenseGeneralMatrix_blacs A(100, 120);
    A.set(10, 10, 10.0);
    A.set(17, 4, 56.0);
    A.set(21, 18, -2.5);
    A.set(0, 10, 12.0);

    TlDenseGeneralMatrix_blacs B(100, 120);
    B.set(10, 10, 5.0);
    B.set(0, 10, -12.0);
    B.set(3, 17, 51.0);

    A += B;

    EXPECT_EQ(100, A.getNumOfRows());
    EXPECT_EQ(120, A.getNumOfCols());
    EXPECT_DOUBLE_EQ(10.0 + 5.0, A.get(10, 10));
    EXPECT_DOUBLE_EQ(12.0 - 12.0, A.get(0, 10));
    EXPECT_DOUBLE_EQ(56.0, A.get(17, 4));
    EXPECT_DOUBLE_EQ(-2.5, A.get(21, 18));
    EXPECT_DOUBLE_EQ(51.0, A.get(3, 17));
}

TEST(TlDenseGeneralMatrix_blacs, operator_mul_matrix) {
    const int rowA = 100;
    const int colA = 120;
    const int& rowB = colA;
    const int colB = 80;
    TlDenseGeneralMatrix_BLAS_old refA, refB;
    TlDenseGeneralMatrix_blacs A, B;

    getMatrix(rowA, colA, &refA, &A);
    getMatrix(rowB, colB, &refB, &B);

    EXPECT_EQ(rowA, A.getNumOfRows());
    EXPECT_EQ(colA, A.getNumOfCols());
    EXPECT_EQ(rowB, B.getNumOfRows());
    EXPECT_EQ(colB, B.getNumOfCols());

    TlDenseGeneralMatrix_BLAS_old refC = refA * refB;
    TlDenseGeneralMatrix_blacs C = A * B;

    EXPECT_EQ(rowA, C.getNumOfRows());
    EXPECT_EQ(colB, C.getNumOfCols());
    for (int i = 0; i < rowA; ++i) {
        for (int j = 0; j < colB; ++j) {
            EXPECT_NEAR(refC.get(i, j), C.get(i, j), EPS);
        }
    }
}

TEST(TlDenseGeneralMatrix_blacs, operator_mul_double) {
    const int rowA = 100;
    const int colA = 120;
    TlDenseGeneralMatrix_BLAS_old refA;
    TlDenseGeneralMatrix_blacs A;

    getMatrix(rowA, colA, &refA, &A);

    EXPECT_EQ(rowA, A.getNumOfRows());
    EXPECT_EQ(colA, A.getNumOfCols());

    TlDenseGeneralMatrix_BLAS_old refC = 3.0 * refA;
    TlDenseGeneralMatrix_blacs C = 3.0 * A;

    EXPECT_EQ(rowA, C.getNumOfRows());
    EXPECT_EQ(colA, C.getNumOfCols());
    for (int i = 0; i < rowA; ++i) {
        for (int j = 0; j < colA; ++j) {
            EXPECT_NEAR(refC.get(i, j), C.get(i, j), EPS);
        }
    }
}

TEST(TlDenseGeneralMatrix_blacs, operator_imul_double) {
    const int rowA = 100;
    const int colA = 120;
    // const int& rowB = colA;
    // const int colB = 80;
    TlDenseGeneralMatrix_BLAS_old refA;
    TlDenseGeneralMatrix_blacs A;

    getMatrix(rowA, colA, &refA, &A);

    EXPECT_EQ(rowA, A.getNumOfRows());
    EXPECT_EQ(colA, A.getNumOfCols());

    refA *= -2.5;
    A *= -2.5;

    EXPECT_EQ(rowA, A.getNumOfRows());
    EXPECT_EQ(colA, A.getNumOfCols());
    for (int i = 0; i < rowA; ++i) {
        for (int j = 0; j < colA; ++j) {
            EXPECT_NEAR(refA.get(i, j), A.get(i, j), EPS);
        }
    }
}

TEST(TlDenseGeneralMatrix_blacs, save) {
    const int rowA = 1000;
    const int colA = 1200;
    TlDenseGeneralMatrix_BLAS_old refA;
    TlDenseGeneralMatrix_blacs A;
    getMatrix(rowA, colA, &refA, &A);

    A.save(distmat_path);

    TlDenseGeneralMatrix_BLAS_old B;
    B.load(distmat_path);

    EXPECT_EQ(rowA, B.getNumOfRows());
    EXPECT_EQ(colA, B.getNumOfCols());
    for (int i = 0; i < rowA; ++i) {
        for (int j = 0; j < colA; ++j) {
            EXPECT_DOUBLE_EQ(refA.get(i, j), B.get(i, j));
        }
    }
}

TEST(TlDenseGeneralMatrix_blacs, load) {
    const int rowA = 1000;
    const int colA = 1200;
    TlDenseGeneralMatrix_BLAS_old refA;
    TlDenseGeneralMatrix_blacs A;
    getMatrix(rowA, colA, &refA, &A);

    refA.save(mat_path);

    TlDenseGeneralMatrix_blacs B;
    B.load(mat_path);

    EXPECT_EQ(rowA, B.getNumOfRows());
    EXPECT_EQ(colA, B.getNumOfCols());
    for (int i = 0; i < rowA; ++i) {
        for (int j = 0; j < colA; ++j) {
            EXPECT_DOUBLE_EQ(refA.get(i, j), B.get(i, j));
        }
    }
}

#ifdef HAVE_HDF5
TEST(TlDenseGeneralMatrix_blacs, save_hdf5) {
    const int rowA = 100;
    const int colA = 120;
    TlDenseGeneralMatrix_BLAS_old refA;
    TlDenseGeneralMatrix_blacs A;
    getMatrix(rowA, colA, &refA, &A);

    A.saveHdf5(distmat_h5, "save_A");

    TlDenseGeneralMatrix_BLAS_old B;
    B.loadHdf5(distmat_h5, "save_A");

    EXPECT_EQ(rowA, B.getNumOfRows());
    EXPECT_EQ(colA, B.getNumOfCols());
    for (int i = 0; i < rowA; ++i) {
        for (int j = 0; j < colA; ++j) {
            EXPECT_DOUBLE_EQ(refA.get(i, j), B.get(i, j));
        }
    }
}

TEST(TlDenseGeneralMatrix_blacs, load_hdf5) {
    const int rowA = 100;
    const int colA = 120;
    TlDenseGeneralMatrix_BLAS_old refA;
    TlDenseGeneralMatrix_blacs A;
    getMatrix(rowA, colA, &refA, &A);

    refA.saveHdf5(distmat_h5, "load_A");

    TlDenseGeneralMatrix_blacs B;
    B.loadHdf5(distmat_h5, "load_A");

    EXPECT_EQ(rowA, B.getNumOfRows());
    EXPECT_EQ(colA, B.getNumOfCols());
    for (int i = 0; i < rowA; ++i) {
        for (int j = 0; j < colA; ++j) {
            EXPECT_DOUBLE_EQ(refA.get(i, j), B.get(i, j));
        }
    }
}

#endif  // HAVE_HDF5

TEST(TlDenseGeneralMatrix_blacs, transpose) {
    const int rowA = 100;
    const int colA = 120;
    TlDenseGeneralMatrix_BLAS_old refA;
    TlDenseGeneralMatrix_blacs A;

    getMatrix(rowA, colA, &refA, &A);

    refA.transposeInPlace();
    A.transposeInPlace();

    EXPECT_EQ(colA, A.getNumOfRows());
    EXPECT_EQ(rowA, A.getNumOfCols());
    for (int i = 0; i < colA; ++i) {
        for (int j = 0; j < rowA; ++j) {
            EXPECT_DOUBLE_EQ(refA.get(i, j), A.get(i, j));
        }
    }
}

TEST(TlDenseGeneralMatrix_blacs, inverse) {
    const int dim = 100;
    TlDenseGeneralMatrix_BLAS_old refA;
    TlDenseGeneralMatrix_blacs A;

    getSymmetricMatrix(dim, &refA, &A);

    TlDenseGeneralMatrix_BLAS_old refB = refA;
    TlDenseGeneralMatrix_blacs B = A;

    refB.inverse();
    B.inverse();

    EXPECT_EQ(refB.getNumOfRows(), B.getNumOfRows());
    EXPECT_EQ(refB.getNumOfCols(), B.getNumOfCols());
    for (int i = 0; i < B.getNumOfRows(); ++i) {
        for (int j = 0; j < B.getNumOfCols(); ++j) {
            EXPECT_NEAR(refB.get(i, j), B.get(i, j), EPS);
        }
    }
}
