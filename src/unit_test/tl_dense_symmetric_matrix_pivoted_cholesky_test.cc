#include <cmath>
#include <vector>

#include "gtest/gtest.h"
#include "tl_dense_general_matrix_lapack.h"
#include "tl_dense_symmetric_matrix_lapack.h"

namespace {
    const double THRESHOLD = 1.0E-8;

    // A = B * B^T (B: N x rank), a symmetric positive semi-definite matrix of the given rank
    TlDenseSymmetricMatrix_Lapack makeLowRank(const int N, const int rank) {
        TlDenseGeneralMatrix_Lapack B(N, rank);
        for (int i = 0; i < N; ++i) {
            for (int k = 0; k < rank; ++k) {
                B.set(i, k, std::sin(1.0 + 0.7 * i * (k + 1)) + 0.3 * (k + 1) + (i == k ? 2.0 : 0.0));
            }
        }
        TlDenseGeneralMatrix_Lapack BBt = B * TlDenseGeneralMatrix_Lapack(B).transpose();
        TlDenseSymmetricMatrix_Lapack A(N);
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j <= i; ++j) {
                A.set(i, j, BBt.get(i, j));
            }
        }
        return A;
    }

    void checkDecomposition(const TlDenseSymmetricMatrix_Lapack& A, const int expectedRank) {
        const int N = A.getNumOfRows();
        TlDenseGeneralMatrix_Lapack L;
        A.pivotedCholeskyDecomposition(&L, THRESHOLD);

        ASSERT_EQ(N, L.getNumOfRows());
        ASSERT_EQ(expectedRank, L.getNumOfCols());
        for (int i = 0; i < N; ++i) {
            for (int k = 0; k < L.getNumOfCols(); ++k) {
                ASSERT_FALSE(std::isnan(L.get(i, k)));
            }
        }

        TlDenseGeneralMatrix_Lapack LLt = L * TlDenseGeneralMatrix_Lapack(L).transpose();
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < N; ++j) {
                EXPECT_NEAR(A.get(i, j), LLt.get(i, j), THRESHOLD);
            }
        }
    }

    // subtract a small amount from the diagonals so that the residual becomes negative once the rank is exhausted
    TlDenseSymmetricMatrix_Lapack makeNoisy(const int N, const int rank) {
        TlDenseSymmetricMatrix_Lapack A = makeLowRank(N, rank);
        for (int i = 0; i < N; ++i) {
            A.add(i, i, -1.0E-12);
        }
        return A;
    }
}  // namespace

TEST(TlDenseSymmetricMatrix_Lapack_PivotedCholesky, rank1) {
    checkDecomposition(makeLowRank(10, 1), 1);
}

TEST(TlDenseSymmetricMatrix_Lapack_PivotedCholesky, rank3) {
    checkDecomposition(makeLowRank(10, 3), 3);
}

TEST(TlDenseSymmetricMatrix_Lapack_PivotedCholesky, fullRank) {
    checkDecomposition(makeLowRank(10, 10), 10);
}

TEST(TlDenseSymmetricMatrix_Lapack_PivotedCholesky, rank1NegativeNoise) {
    checkDecomposition(makeNoisy(10, 1), 1);
}

TEST(TlDenseSymmetricMatrix_Lapack_PivotedCholesky, rank3NegativeNoise) {
    checkDecomposition(makeNoisy(10, 3), 3);
}
