#include "tl_sparse_matrix.h"
#include <limits>
#include "gtest/gtest.h"

static const double EPS = 1.0E-10;  // std::numeric_limits<double>::epsilon();

TEST(TlSparseMatrix, constructor) {
    TlSparseMatrix a(3, 3);

    EXPECT_NEAR(0.0, a.get(0, 0), EPS);
    EXPECT_NEAR(0.0, a.get(0, 1), EPS);
    EXPECT_NEAR(0.0, a.get(0, 2), EPS);
    EXPECT_NEAR(0.0, a.get(1, 0), EPS);
    EXPECT_NEAR(0.0, a.get(1, 1), EPS);
    EXPECT_NEAR(0.0, a.get(1, 2), EPS);
    EXPECT_NEAR(0.0, a.get(2, 0), EPS);
    EXPECT_NEAR(0.0, a.get(2, 1), EPS);
    EXPECT_NEAR(0.0, a.get(2, 2), EPS);
}

TEST(TlSparseMatrix, copyConstructor) {
    TlSparseMatrix b(5, 5);
    b.set(0, 0, 1.0);
    b.set(1, 0, 2.0);
    b.set(2, 0, 3.0);
    b.set(3, 3, 4.0);

    TlSparseMatrix a;
    a = b;

    EXPECT_NEAR(1.0, a.get(0, 0), EPS);
    EXPECT_NEAR(0.0, a.get(0, 1), EPS);
    EXPECT_NEAR(0.0, a.get(0, 2), EPS);
    EXPECT_NEAR(0.0, a.get(0, 3), EPS);
    EXPECT_NEAR(0.0, a.get(0, 4), EPS);

    EXPECT_NEAR(2.0, a.get(1, 0), EPS);
    EXPECT_NEAR(0.0, a.get(1, 1), EPS);
    EXPECT_NEAR(0.0, a.get(1, 2), EPS);
    EXPECT_NEAR(0.0, a.get(1, 3), EPS);
    EXPECT_NEAR(0.0, a.get(1, 4), EPS);

    EXPECT_NEAR(3.0, a.get(2, 0), EPS);
    EXPECT_NEAR(0.0, a.get(2, 1), EPS);
    EXPECT_NEAR(0.0, a.get(2, 2), EPS);
    EXPECT_NEAR(0.0, a.get(2, 3), EPS);
    EXPECT_NEAR(0.0, a.get(2, 4), EPS);

    EXPECT_NEAR(0.0, a.get(3, 0), EPS);
    EXPECT_NEAR(0.0, a.get(3, 1), EPS);
    EXPECT_NEAR(0.0, a.get(3, 2), EPS);
    EXPECT_NEAR(4.0, a.get(3, 3), EPS);
    EXPECT_NEAR(0.0, a.get(3, 4), EPS);

    EXPECT_NEAR(0.0, a.get(4, 0), EPS);
    EXPECT_NEAR(0.0, a.get(4, 1), EPS);
    EXPECT_NEAR(0.0, a.get(4, 2), EPS);
    EXPECT_NEAR(0.0, a.get(4, 3), EPS);
    EXPECT_NEAR(0.0, a.get(4, 4), EPS);
}

TEST(TlSparseMatrix, merge) {
    TlSparseMatrix a(5, 5);
    TlSparseMatrix b(5, 5);

    a.set(0, 0, 1.0);
    a.set(1, 0, 2.0);
    b.set(2, 0, 3.0);
    b.set(3, 3, 4.0);

    a.merge(b);

    EXPECT_NEAR(1.0, a.get(0, 0), EPS);
    EXPECT_NEAR(0.0, a.get(0, 1), EPS);
    EXPECT_NEAR(0.0, a.get(0, 2), EPS);
    EXPECT_NEAR(0.0, a.get(0, 3), EPS);
    EXPECT_NEAR(0.0, a.get(0, 4), EPS);

    EXPECT_NEAR(2.0, a.get(1, 0), EPS);
    EXPECT_NEAR(0.0, a.get(1, 1), EPS);
    EXPECT_NEAR(0.0, a.get(1, 2), EPS);
    EXPECT_NEAR(0.0, a.get(1, 3), EPS);
    EXPECT_NEAR(0.0, a.get(1, 4), EPS);

    EXPECT_NEAR(3.0, a.get(2, 0), EPS);
    EXPECT_NEAR(0.0, a.get(2, 1), EPS);
    EXPECT_NEAR(0.0, a.get(2, 2), EPS);
    EXPECT_NEAR(0.0, a.get(2, 3), EPS);
    EXPECT_NEAR(0.0, a.get(2, 4), EPS);

    EXPECT_NEAR(0.0, a.get(3, 0), EPS);
    EXPECT_NEAR(0.0, a.get(3, 1), EPS);
    EXPECT_NEAR(0.0, a.get(3, 2), EPS);
    EXPECT_NEAR(4.0, a.get(3, 3), EPS);
    EXPECT_NEAR(0.0, a.get(3, 4), EPS);

    EXPECT_NEAR(0.0, a.get(4, 0), EPS);
    EXPECT_NEAR(0.0, a.get(4, 1), EPS);
    EXPECT_NEAR(0.0, a.get(4, 2), EPS);
    EXPECT_NEAR(0.0, a.get(4, 3), EPS);
    EXPECT_NEAR(0.0, a.get(4, 4), EPS);
}
