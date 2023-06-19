#ifdef HAVE_CONFIG_H
#include "config.h"
#endif  // HAVE_CONFIG_H

#include "TlCommunicate.h"
#include "gtest/gtest.h"
#include "tl_dense_general_matrix_eigen.h"
#include "tl_dense_symmetric_matrix_eigen.h"
#include "tl_dense_vector_eigen.h"

// -----------------------------------------------------------------------------
// broadcast
// -----------------------------------------------------------------------------
TEST(TlCommunicate_eigen, broadcastGeneralMatrix) {
    TlCommunicate& rComm = TlCommunicate::getInstance();

    const int row = 20;
    const int col = 30;
    TlDenseGeneralMatrix_Eigen mat(row, col);

    if (rComm.isMaster()) {
        int index = 0;
        for (int r = 0; r < row; ++r) {
            for (int c = 0; c < col; ++c) {
                mat.set(r, c, index * 0.01);
                ++index;
            }
        }
    }

    rComm.broadcast(&mat);
    {
        EXPECT_EQ(row, mat.getNumOfRows());
        EXPECT_EQ(col, mat.getNumOfCols());

        int index = 0;
        for (int r = 0; r < row; ++r) {
            for (int c = 0; c < col; ++c) {
                EXPECT_DOUBLE_EQ(index * 0.01, mat.get(r, c));
                ++index;
            }
        }
    }
}

TEST(TlCommunicate_eigen, broadcastSymmetricMatrix) {
    TlCommunicate& rComm = TlCommunicate::getInstance();

    const int dim = 30;
    TlDenseSymmetricMatrix_Eigen mat(dim);

    if (rComm.isMaster()) {
        int index = 0;
        for (int r = 0; r < dim; ++r) {
            for (int c = 0; c <= r; ++c) {
                mat.set(r, c, index * 0.01);
                ++index;
            }
        }
    }

    rComm.broadcast(&mat);
    {
        EXPECT_EQ(dim, mat.getNumOfRows());
        EXPECT_EQ(dim, mat.getNumOfCols());

        int index = 0;
        for (int r = 0; r < dim; ++r) {
            for (int c = 0; c <= r; ++c) {
                EXPECT_DOUBLE_EQ(index * 0.01, mat.get(r, c));
                EXPECT_DOUBLE_EQ(index * 0.01, mat.get(c, r));
                ++index;
            }
        }
    }
}

TEST(TlCommunicate_eigen, broadcastDenseVector) {
    TlCommunicate& rComm = TlCommunicate::getInstance();

    const int dim = 30;
    TlDenseVector_Eigen vtr(dim);

    if (rComm.isMaster()) {
        int index = 0;
        for (int i = 0; i < dim; ++i) {
            vtr.set(i, index * 0.01);
            ++index;
        }
    }

    rComm.broadcast(&vtr);
    {
        EXPECT_EQ(dim, vtr.getSize());
        int index = 0;
        for (int i = 0; i < dim; ++i) {
            EXPECT_DOUBLE_EQ(index * 0.01, vtr.get(i));
            ++index;
        }
    }
}

// -----------------------------------------------------------------------------
// all_reduce(SUM)
// -----------------------------------------------------------------------------
TEST(TlCommunicate_eigen, allReduceSumGeneralMatrix) {
    TlCommunicate& rComm = TlCommunicate::getInstance();
    const int proc = rComm.getNumOfProc();

    const int row = 20;
    const int col = 30;
    TlDenseGeneralMatrix_Eigen mat(row, col);

    {
        int index = 0;
        for (int r = 0; r < row; ++r) {
            for (int c = 0; c < col; ++c) {
                mat.set(r, c, index * 0.01);
                ++index;
            }
        }
    }

    rComm.allReduce_SUM(&mat);
    {
        EXPECT_EQ(row, mat.getNumOfRows());
        EXPECT_EQ(col, mat.getNumOfCols());

        int index = 0;
        for (int r = 0; r < row; ++r) {
            for (int c = 0; c < col; ++c) {
                EXPECT_DOUBLE_EQ(index * 0.01 * proc, mat.get(r, c));
                ++index;
            }
        }
    }
}

TEST(TlCommunicate_eigen, allReduceSumSymmetricMatrix) {
    TlCommunicate& rComm = TlCommunicate::getInstance();
    const int proc = rComm.getNumOfProc();

    const int dim = 100;
    TlDenseSymmetricMatrix_Eigen mat(dim);

    {
        int index = 0;
        for (int r = 0; r < dim; ++r) {
            for (int c = 0; c <= r; ++c) {
                mat.set(r, c, index * 0.01);
                ++index;
            }
        }
    }

    rComm.allReduce_SUM(&mat);
    {
        EXPECT_EQ(dim, mat.getNumOfRows());
        EXPECT_EQ(dim, mat.getNumOfCols());

        int index = 0;
        for (int r = 0; r < dim; ++r) {
            for (int c = 0; c <= r; ++c) {
                EXPECT_DOUBLE_EQ(index * 0.01 * proc, mat.get(r, c));
                EXPECT_DOUBLE_EQ(index * 0.01 * proc, mat.get(c, r));
                ++index;
            }
        }
    }
}
