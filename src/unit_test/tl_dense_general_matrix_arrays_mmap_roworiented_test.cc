#include "tl_dense_general_matrix_arrays_mmap_roworiented.h"

#include <iostream>
#include <limits>

#include "TlFile.h"
#include "TlUtils.h"
#include "gtest/gtest.h"
#include "tl_dense_general_matrix_lapack.h"

static const double EPS = 1.0E-10;  // std::numeric_limits<double>::epsilon();

static const std::string matPath1 = "temp.array_mmap1.mat";
static const std::string matPath2 = "temp.array_mmap2.mat";
static const std::string csfdPath = "temp.csfd.mat";

static void cleanup(const std::string& baseFilePath, const int numOfSubunits = 0) {
    for (int i = 0; i < numOfSubunits; ++i) {
        std::string path = TlDenseMatrix_arrays_mmap_Object::getFileName(baseFilePath, i);
        if (TlFile::isExistFile(path)) {
            TlFile::remove(path);
        }
    }
}

TEST(TlDenseGeneralMatrix_arrays_mmap_RowOriented, constructor) {
    const int numOfRows = 1000;
    const int numOfCols = 2000;
    const int numOfSubunits = 1;
    cleanup(matPath1, numOfSubunits);

    TlDenseGeneralMatrix_arrays_mmap_RowOriented A(matPath1, numOfRows, numOfCols);

    EXPECT_EQ(numOfRows, A.getNumOfRows());
    EXPECT_EQ(numOfCols, A.getNumOfCols());
    EXPECT_NEAR(0.0, A.get(0, 0), EPS);
}

TEST(TlDenseGeneralMatrix_arrays_mmap_RowOriented, constructor_group) {
    const int numOfRows = 1000;
    const int numOfCols = 2000;
    const int numOfSubunits = 4;
    cleanup(matPath1, numOfSubunits);

    TlDenseGeneralMatrix_arrays_mmap_RowOriented A0(matPath1, numOfRows, numOfCols, numOfSubunits, 0);
    TlDenseGeneralMatrix_arrays_mmap_RowOriented A1(matPath1, numOfRows, numOfCols, numOfSubunits, 1);

    EXPECT_EQ(numOfRows, A0.getNumOfRows());
    EXPECT_EQ(numOfCols, A0.getNumOfCols());
    EXPECT_EQ(numOfSubunits, A0.getNumOfSubunits());
    EXPECT_EQ(0, A0.getSubunitID());

    EXPECT_EQ(numOfRows, A1.getNumOfRows());
    EXPECT_EQ(numOfCols, A1.getNumOfCols());
    EXPECT_EQ(numOfSubunits, A1.getNumOfSubunits());
    EXPECT_EQ(1, A1.getSubunitID());
}

TEST(TlDenseGeneralMatrix_arrays_mmap_RowOriented, contents) {
    const int numOfRows = 1000;
    const int numOfCols = 2000;
    const int numOfSubunits = 1;
    cleanup(matPath1, numOfSubunits);

    TlDenseGeneralMatrix_Lapack matA(numOfRows, numOfCols);
    TlDenseGeneralMatrix_arrays_mmap_RowOriented arrayMatA(matPath1, numOfRows, numOfCols);

    // setup
    int count = 0;
    for (int r = 0; r < numOfRows; ++r) {
        for (int c = 0; c < numOfCols; ++c) {
            double v = double(count);
            matA.set(r, c, v);
            arrayMatA.set(r, c, v);

            ++count;
        }
    }

    // test
    for (int r = 0; r < numOfRows; ++r) {
        for (int c = 0; c < numOfCols; ++c) {
            // std::cerr << TlUtils::format("(%d, %d) %f, %f", r, c, matA.get(r, c), arrayMatA.get(r, c)) << std::endl;
            EXPECT_NEAR(matA.get(r, c), arrayMatA.get(r, c), EPS);
        }
    }
}

TEST(TlDenseGeneralMatrix_arrays_mmap_RowOriented, contents_reserved) {
    const int numOfRows = 1000;
    const int numOfCols = 2000;
    const int numOfSubunits = 1;
    cleanup(matPath1, numOfSubunits);

    TlDenseGeneralMatrix_Lapack matA(numOfRows, numOfCols);
    TlDenseGeneralMatrix_arrays_mmap_RowOriented arrayMatA(matPath1, numOfRows, numOfCols, numOfSubunits, 0, 3000);

    // setup
    int count = 0;
    for (int r = 0; r < numOfRows; ++r) {
        for (int c = 0; c < numOfCols; ++c) {
            double v = double(count);
            matA.set(r, c, v);
            arrayMatA.set(r, c, v);

            ++count;
        }
    }

    // test
    for (int r = 0; r < numOfRows; ++r) {
        for (int c = 0; c < numOfCols; ++c) {
            // std::cerr << TlUtils::format("(%d, %d) %f, %f", r, c, matA.get(r, c), arrayMatA.get(r, c)) << std::endl;
            EXPECT_NEAR(matA.get(r, c), arrayMatA.get(r, c), EPS);
        }
    }
}

TEST(TlDenseGeneralMatrix_arrays_mmap_RowOriented, save_load) {
    const int numOfRows = 1000;
    const int numOfCols = 2000;
    const int numOfSubunits = 1;
    cleanup(matPath1, numOfSubunits);

    TlDenseGeneralMatrix_Lapack matA(numOfRows, numOfCols);

    // setup
    {
        TlDenseGeneralMatrix_arrays_mmap_RowOriented arrayMatA(matPath1, numOfRows, numOfCols);
        int count = 0;
        for (int r = 0; r < numOfRows; ++r) {
            for (int c = 0; c < numOfCols; ++c) {
                double v = double(count);
                matA.set(r, c, v);
                arrayMatA.set(r, c, v);

                ++count;
            }
        }
    }

    // load new matrix
    const std::string path = matPath1 + ".part0.mat";
    TlDenseGeneralMatrix_arrays_mmap_RowOriented arrayMatB(path);

    // test
    EXPECT_EQ(numOfRows, arrayMatB.getNumOfRows());
    EXPECT_EQ(numOfCols, arrayMatB.getNumOfCols());
    for (int r = 0; r < numOfRows; ++r) {
        for (int c = 0; c < numOfCols; ++c) {
            // std::cerr << TlUtils::format("(%d, %d) %f, %f", r, c, matA.get(r, c), arrayMatB.get(r, c)) << std::endl;
            EXPECT_NEAR(matA.get(r, c), arrayMatB.get(r, c), EPS);
        }
    }
}

TEST(TlDenseGeneralMatrix_arrays_mmap_RowOriented, resize_same) {
    // resize (same size)
    const int numOfRows = 1000;
    const int numOfCols = 2000;
    const int numOfSubunits = 1;
    cleanup(matPath1, numOfSubunits);

    TlDenseGeneralMatrix_arrays_mmap_RowOriented A(matPath1, numOfRows, numOfCols);
    EXPECT_EQ(numOfRows, A.getNumOfRows());
    EXPECT_EQ(numOfCols, A.getNumOfCols());

    const int newNumOfRows = 1000;
    const int newNumOfCols = 2000;

    A.resize(newNumOfRows, newNumOfCols);
    EXPECT_EQ(newNumOfRows, A.getNumOfRows());
    EXPECT_EQ(newNumOfCols, A.getNumOfCols());
}

TEST(TlDenseGeneralMatrix_arrays_mmap_RowOriented, resize_expand_row) {
    const int numOfRows = 1000;
    const int numOfCols = 2000;
    const int numOfSubunits = 1;
    cleanup(matPath1, numOfSubunits);

    TlDenseGeneralMatrix_arrays_mmap_RowOriented A(matPath1, numOfRows, numOfCols);
    EXPECT_EQ(numOfRows, A.getNumOfRows());
    EXPECT_EQ(numOfCols, A.getNumOfCols());

    const int newNumOfRows = 2000;
    const int newNumOfCols = 2000;

    A.resize(newNumOfRows, newNumOfCols);
    EXPECT_EQ(newNumOfRows, A.getNumOfRows());
    EXPECT_EQ(newNumOfCols, A.getNumOfCols());
}

TEST(TlDenseGeneralMatrix_arrays_mmap_RowOriented, resize_expand_col) {
    const int numOfRows = 1000;
    const int numOfCols = 2000;
    const int numOfSubunits = 1;
    cleanup(matPath1, numOfSubunits);

    TlDenseGeneralMatrix_arrays_mmap_RowOriented A(matPath1, numOfRows, numOfCols);
    EXPECT_EQ(numOfRows, A.getNumOfRows());
    EXPECT_EQ(numOfCols, A.getNumOfCols());

    const int newNumOfRows = 1000;
    const int newNumOfCols = 3000;

    A.resize(newNumOfRows, newNumOfCols);
    EXPECT_EQ(newNumOfRows, A.getNumOfRows());
    EXPECT_EQ(newNumOfCols, A.getNumOfCols());
}

TEST(TlDenseGeneralMatrix_arrays_mmap_RowOriented, resize_expand_row_col) {
    const int numOfRows = 1000;
    const int numOfCols = 2000;
    const int numOfSubunits = 1;
    cleanup(matPath1, numOfSubunits);

    TlDenseGeneralMatrix_arrays_mmap_RowOriented A(matPath1, numOfRows, numOfCols);
    EXPECT_EQ(numOfRows, A.getNumOfRows());
    EXPECT_EQ(numOfCols, A.getNumOfCols());

    const int newNumOfRows = 2000;
    const int newNumOfCols = 3000;

    A.resize(newNumOfRows, newNumOfCols);
    EXPECT_EQ(newNumOfRows, A.getNumOfRows());
    EXPECT_EQ(newNumOfCols, A.getNumOfCols());
}

TEST(TlDenseGeneralMatrix_arrays_mmap_RowOriented, resize_reduction_row) {
    const int numOfRows = 1000;
    const int numOfCols = 2000;
    const int numOfSubunits = 1;
    cleanup(matPath1, numOfSubunits);

    TlDenseGeneralMatrix_arrays_mmap_RowOriented A(matPath1, numOfRows, numOfCols);
    EXPECT_EQ(numOfRows, A.getNumOfRows());
    EXPECT_EQ(numOfCols, A.getNumOfCols());

    const int newNumOfRows = 500;
    const int newNumOfCols = 2000;

    A.resize(newNumOfRows, newNumOfCols);
    EXPECT_EQ(newNumOfRows, A.getNumOfRows());
    EXPECT_EQ(newNumOfCols, A.getNumOfCols());
}

TEST(TlDenseGeneralMatrix_arrays_mmap_RowOriented, resize_reduction_col) {
    const int numOfRows = 1000;
    const int numOfCols = 2000;
    const int numOfSubunits = 1;
    cleanup(matPath1, numOfSubunits);

    TlDenseGeneralMatrix_arrays_mmap_RowOriented A(matPath1, numOfRows, numOfCols);
    EXPECT_EQ(numOfRows, A.getNumOfRows());
    EXPECT_EQ(numOfCols, A.getNumOfCols());

    const int newNumOfRows = 1000;
    const int newNumOfCols = 1000;

    A.resize(newNumOfRows, newNumOfCols);
    EXPECT_EQ(newNumOfRows, A.getNumOfRows());
    EXPECT_EQ(newNumOfCols, A.getNumOfCols());
}

TEST(TlDenseGeneralMatrix_arrays_mmap_RowOriented, resize_reduction_row_col) {
    const int numOfRows = 1000;
    const int numOfCols = 2000;
    const int numOfSubunits = 1;
    cleanup(matPath1, numOfSubunits);

    TlDenseGeneralMatrix_arrays_mmap_RowOriented A(matPath1, numOfRows, numOfCols);
    EXPECT_EQ(numOfRows, A.getNumOfRows());
    EXPECT_EQ(numOfCols, A.getNumOfCols());

    const int newNumOfRows = 500;
    const int newNumOfCols = 1000;

    A.resize(newNumOfRows, newNumOfCols);
    EXPECT_EQ(newNumOfRows, A.getNumOfRows());
    EXPECT_EQ(newNumOfCols, A.getNumOfCols());
}

TEST(TlDenseGeneralMatrix_arrays_mmap_RowOriented, resize_same_group) {
    // resize (same size)
    const int numOfRows = 1000;
    const int numOfCols = 2000;
    const int numOfSubunits = 4;
    cleanup(matPath1, numOfSubunits);

    // construct
    std::vector<TlDenseGeneralMatrix_arrays_mmap_RowOriented*> pMatrices(numOfSubunits, NULL);
    for (int id = 0; id < numOfSubunits; ++id) {
        pMatrices[id] =
            new TlDenseGeneralMatrix_arrays_mmap_RowOriented(matPath1, numOfRows, numOfCols, numOfSubunits, id);
    }

    // check
    for (int id = 0; id < numOfSubunits; ++id) {
        EXPECT_EQ(numOfRows, pMatrices[id]->getNumOfRows());
        EXPECT_EQ(numOfCols, pMatrices[id]->getNumOfCols());
    }

    const int newNumOfRows = 1000;
    const int newNumOfCols = 2000;

    // resize
    for (int id = 0; id < numOfSubunits; ++id) {
        pMatrices[id]->resize(newNumOfRows, newNumOfCols);
    }

    // check
    for (int id = 0; id < numOfSubunits; ++id) {
        EXPECT_EQ(newNumOfRows, pMatrices[id]->getNumOfRows());
        EXPECT_EQ(newNumOfCols, pMatrices[id]->getNumOfCols());
    }
    for (int r = 0; r < newNumOfRows; ++r) {
        for (int c = 0; c < newNumOfCols; ++c) {
            const int id = pMatrices[0]->getSubunitID(r);
            EXPECT_NEAR(0.0, pMatrices[id]->get(r, c), EPS);
        }
    }
}

TEST(TlDenseGeneralMatrix_arrays_mmap_RowOriented, resize_expand_row_group) {
    const int numOfRows = 1000;
    const int numOfCols = 2000;
    const int numOfSubunits = 4;
    cleanup(matPath1, numOfSubunits);

    // construct
    std::vector<TlDenseGeneralMatrix_arrays_mmap_RowOriented*> pMatrices(numOfSubunits, NULL);
    for (int id = 0; id < numOfSubunits; ++id) {
        pMatrices[id] =
            new TlDenseGeneralMatrix_arrays_mmap_RowOriented(matPath1, numOfRows, numOfCols, numOfSubunits, id);
    }
    TlDenseGeneralMatrix_Lapack refMat(numOfRows, numOfCols);

    // check
    for (int id = 0; id < numOfSubunits; ++id) {
        EXPECT_EQ(numOfRows, pMatrices[id]->getNumOfRows());
        EXPECT_EQ(numOfCols, pMatrices[id]->getNumOfCols());
    }

    // prepare
    double val = 0.0;
    for (int i = 0; i < numOfRows; ++i) {
        for (int j = 0; j < numOfCols; ++j) {
            for (int id = 0; id < numOfSubunits; ++id) {
                pMatrices[id]->set(i, j, val);
            }
            refMat.set(i, j, val);
            val += 0.01;
        }
    }

    // resize
    const int newNumOfRows = 2000;
    const int newNumOfCols = 2000;
    for (int id = 0; id < numOfSubunits; ++id) {
        pMatrices[id]->resize(newNumOfRows, newNumOfCols);
    }
    refMat.resize(newNumOfRows, newNumOfCols);

    for (int id = 0; id < numOfSubunits; ++id) {
        EXPECT_EQ(newNumOfRows, pMatrices[id]->getNumOfRows());
        EXPECT_EQ(newNumOfCols, pMatrices[id]->getNumOfCols());
    }

    // check
    for (int i = 0; i < newNumOfRows; ++i) {
        for (int j = 0; j < newNumOfCols; ++j) {
            const double val = refMat.get(i, j);

            const int subunit = pMatrices[0]->getSubunitID(i);
            for (int id = 0; id < numOfSubunits; ++id) {
                const int subunit2 = pMatrices[id]->getSubunitID(i);
                EXPECT_EQ(subunit, subunit2);
            }
            const double val2 = pMatrices[subunit]->get(i, j);
            // std::cerr << TlUtils::format("(%5d, %5d): %f, %f", i, j, val, val2) << std::endl;
            EXPECT_DOUBLE_EQ(val, val2);
        }
    }

    // destroy
    for (int id = 0; id < numOfSubunits; ++id) {
        delete pMatrices[id];
        pMatrices[id] = NULL;
    }
}

TEST(TlDenseGeneralMatrix_arrays_mmap_RowOriented, resize_expand_col_group) {
    const int numOfRows = 1000;
    const int numOfCols = 2000;
    const int numOfSubunits = 4;
    cleanup(matPath1, numOfSubunits);

    // construct
    std::vector<TlDenseGeneralMatrix_arrays_mmap_RowOriented*> pMatrices(numOfSubunits, NULL);
    for (int id = 0; id < numOfSubunits; ++id) {
        pMatrices[id] =
            new TlDenseGeneralMatrix_arrays_mmap_RowOriented(matPath1, numOfRows, numOfCols, numOfSubunits, id);
    }
    TlDenseGeneralMatrix_Lapack refMat(numOfRows, numOfCols);

    // check
    for (int id = 0; id < numOfSubunits; ++id) {
        EXPECT_EQ(numOfRows, pMatrices[id]->getNumOfRows());
        EXPECT_EQ(numOfCols, pMatrices[id]->getNumOfCols());
    }

    // prepare
    double val = 0.0;
    for (int i = 0; i < numOfRows; ++i) {
        for (int j = 0; j < numOfCols; ++j) {
            for (int id = 0; id < numOfSubunits; ++id) {
                pMatrices[id]->set(i, j, val);
            }
            refMat.set(i, j, val);
            val += 0.01;
        }
    }

    // resize
    const int newNumOfRows = 1000;
    const int newNumOfCols = 3000;
    for (int id = 0; id < numOfSubunits; ++id) {
        pMatrices[id]->resize(newNumOfRows, newNumOfCols);
    }
    refMat.resize(newNumOfRows, newNumOfCols);

    for (int id = 0; id < numOfSubunits; ++id) {
        EXPECT_EQ(newNumOfRows, pMatrices[id]->getNumOfRows());
        EXPECT_EQ(newNumOfCols, pMatrices[id]->getNumOfCols());
    }

    // check
    for (int i = 0; i < newNumOfRows; ++i) {
        for (int j = 0; j < newNumOfCols; ++j) {
            const double val = refMat.get(i, j);

            const int subunit = pMatrices[0]->getSubunitID(i);
            for (int id = 0; id < numOfSubunits; ++id) {
                const int subunit2 = pMatrices[id]->getSubunitID(i);
                EXPECT_EQ(subunit, subunit2);
            }
            const double val2 = pMatrices[subunit]->get(i, j);
            // std::cerr << TlUtils::format("(%5d, %5d): %f, %f", i, j, val, val2) << std::endl;
            EXPECT_DOUBLE_EQ(val, val2);
        }
    }

    // destroy
    for (int id = 0; id < numOfSubunits; ++id) {
        delete pMatrices[id];
        pMatrices[id] = NULL;
    }
}

TEST(TlDenseGeneralMatrix_arrays_mmap_RowOriented, resize_expand_row_col_group) {
    // resize (same size)
    const int numOfRows = 1000;
    const int numOfCols = 2000;
    const int numOfSubunits = 4;
    cleanup(matPath1, numOfSubunits);

    // construct
    std::vector<TlDenseGeneralMatrix_arrays_mmap_RowOriented*> pMatrices(numOfSubunits, NULL);
    for (int id = 0; id < numOfSubunits; ++id) {
        pMatrices[id] =
            new TlDenseGeneralMatrix_arrays_mmap_RowOriented(matPath1, numOfRows, numOfCols, numOfSubunits, id);
    }

    // check
    for (int id = 0; id < numOfSubunits; ++id) {
        EXPECT_EQ(numOfRows, pMatrices[id]->getNumOfRows());
        EXPECT_EQ(numOfCols, pMatrices[id]->getNumOfCols());
    }

    const int newNumOfRows = 2000;
    const int newNumOfCols = 3000;

    // resize
    for (int id = 0; id < numOfSubunits; ++id) {
        pMatrices[id]->resize(newNumOfRows, newNumOfCols);
    }

    // check
    for (int id = 0; id < numOfSubunits; ++id) {
        EXPECT_EQ(newNumOfRows, pMatrices[id]->getNumOfRows());
        EXPECT_EQ(newNumOfCols, pMatrices[id]->getNumOfCols());
    }
    for (int r = 0; r < newNumOfRows; ++r) {
        for (int c = 0; c < newNumOfCols; ++c) {
            const int id = pMatrices[0]->getSubunitID(r);
            EXPECT_NEAR(0.0, pMatrices[id]->get(r, c), EPS);
        }
    }
}

TEST(TlDenseGeneralMatrix_arrays_mmap_RowOriented, resize_reduce_row_col_group) {
    const int numOfRows = 1000;
    const int numOfCols = 2000;
    const int numOfSubunits = 4;
    cleanup(matPath1, numOfSubunits);

    // construct
    std::vector<TlDenseGeneralMatrix_arrays_mmap_RowOriented*> pMatrices(numOfSubunits, NULL);
    for (int id = 0; id < numOfSubunits; ++id) {
        pMatrices[id] =
            new TlDenseGeneralMatrix_arrays_mmap_RowOriented(matPath1, numOfRows, numOfCols, numOfSubunits, id);
    }
    TlDenseGeneralMatrix_Lapack refMat(numOfRows, numOfCols);

    // check
    for (int id = 0; id < numOfSubunits; ++id) {
        EXPECT_EQ(numOfRows, pMatrices[id]->getNumOfRows());
        EXPECT_EQ(numOfCols, pMatrices[id]->getNumOfCols());
    }

    // prepare
    double val = 0.0;
    for (int i = 0; i < numOfRows; ++i) {
        for (int j = 0; j < numOfCols; ++j) {
            for (int id = 0; id < numOfSubunits; ++id) {
                pMatrices[id]->set(i, j, val);
            }
            refMat.set(i, j, val);
            val += 0.01;
        }
    }

    // resize
    const int newNumOfRows = 500;
    const int newNumOfCols = 1000;
    for (int id = 0; id < numOfSubunits; ++id) {
        pMatrices[id]->resize(newNumOfRows, newNumOfCols);
    }
    refMat.resize(newNumOfRows, newNumOfCols);

    for (int id = 0; id < numOfSubunits; ++id) {
        EXPECT_EQ(newNumOfRows, pMatrices[id]->getNumOfRows());
        EXPECT_EQ(newNumOfCols, pMatrices[id]->getNumOfCols());
    }

    // check
    for (int i = 0; i < newNumOfRows; ++i) {
        for (int j = 0; j < newNumOfCols; ++j) {
            const double val = refMat.get(i, j);

            const int subunit = pMatrices[0]->getSubunitID(i);
            for (int id = 0; id < numOfSubunits; ++id) {
                const int subunit2 = pMatrices[id]->getSubunitID(i);
                EXPECT_EQ(subunit, subunit2);
            }
            const double val2 = pMatrices[subunit]->get(i, j);
            // std::cerr << TlUtils::format("(%5d, %5d): %f, %f", i, j, val, val2) << std::endl;
            EXPECT_DOUBLE_EQ(val, val2);
        }
    }

    // destroy
    for (int id = 0; id < numOfSubunits; ++id) {
        delete pMatrices[id];
        pMatrices[id] = NULL;
    }
}

TEST(TlDenseGeneralMatrix_arrays_mmap_RowOriented, reserve) {
    const int numOfRows = 1000;
    const int numOfCols = 2000;
    const int numOfSubunits = 1;
    cleanup(matPath1, numOfSubunits);

    // construct
    std::vector<TlDenseGeneralMatrix_arrays_mmap_RowOriented*> pMatrices(numOfSubunits, NULL);
    for (int id = 0; id < numOfSubunits; ++id) {
        pMatrices[id] =
            new TlDenseGeneralMatrix_arrays_mmap_RowOriented(matPath1, numOfRows, numOfCols, numOfSubunits, id);
    }
    TlDenseGeneralMatrix_Lapack refMat(numOfRows, numOfCols);

    for (int id = 0; id < numOfSubunits; ++id) {
        EXPECT_EQ(numOfRows, pMatrices[id]->getNumOfRows());
        EXPECT_EQ(numOfCols, pMatrices[id]->getNumOfCols());
    }

    // prepare
    double val = 0.0;
    for (int i = 0; i < numOfRows; ++i) {
        for (int j = 0; j < numOfCols; ++j) {
            for (int id = 0; id < numOfSubunits; ++id) {
                pMatrices[id]->set(i, j, val);
            }
            refMat.set(i, j, val);
            val += 0.01;
        }
    }

    // reserve
    for (int id = 0; id < numOfSubunits; ++id) {
        pMatrices[id]->reserveColSize(numOfCols + 1000);
    }

    // check
    for (int i = 0; i < numOfRows; ++i) {
        for (int j = 0; j < numOfCols; ++j) {
            const double val = refMat.get(i, j);

            const int subunit = pMatrices[0]->getSubunitID(i);
            for (int id = 0; id < numOfSubunits; ++id) {
                const int subunit2 = pMatrices[id]->getSubunitID(i);
                EXPECT_EQ(subunit, subunit2);
            }
            const double val2 = pMatrices[subunit]->get(i, j);
            EXPECT_DOUBLE_EQ(val, val2);
            // std::cerr << TlUtils::format("M(%d, %d)[%d]", i, j, subunit);
        }
    }

    // destroy
    for (int id = 0; id < numOfSubunits; ++id) {
        delete pMatrices[id];
        pMatrices[id] = NULL;
    }
}

TEST(TlDenseGeneralMatrix_arrays_mmap_RowOriented, reserve_group) {
    const int numOfRows = 1000;
    const int numOfCols = 2000;
    const int numOfSubunits = 4;
    cleanup(matPath1, numOfSubunits);

    // construct
    std::vector<TlDenseGeneralMatrix_arrays_mmap_RowOriented*> pMatrices(numOfSubunits, NULL);
    for (int id = 0; id < numOfSubunits; ++id) {
        pMatrices[id] =
            new TlDenseGeneralMatrix_arrays_mmap_RowOriented(matPath1, numOfRows, numOfCols, numOfSubunits, id);
    }
    TlDenseGeneralMatrix_Lapack refMat(numOfRows, numOfCols);

    for (int id = 0; id < numOfSubunits; ++id) {
        EXPECT_EQ(numOfRows, pMatrices[id]->getNumOfRows());
        EXPECT_EQ(numOfCols, pMatrices[id]->getNumOfCols());
    }

    // prepare
    double val = 0.0;
    for (int i = 0; i < numOfRows; ++i) {
        for (int j = 0; j < numOfCols; ++j) {
            for (int id = 0; id < numOfSubunits; ++id) {
                pMatrices[id]->set(i, j, val);
            }
            refMat.set(i, j, val);
            val += 0.01;
        }
    }

    // reserve
    for (int id = 0; id < numOfSubunits; ++id) {
        pMatrices[id]->reserveColSize(numOfCols + 1000);
    }

    // check
    for (int i = 0; i < numOfRows; ++i) {
        for (int j = 0; j < numOfCols; ++j) {
            const double val = refMat.get(i, j);

            const int subunit = pMatrices[0]->getSubunitID(i);
            for (int id = 0; id < numOfSubunits; ++id) {
                const int subunit2 = pMatrices[id]->getSubunitID(i);
                EXPECT_EQ(subunit, subunit2);
            }
            const double val2 = pMatrices[subunit]->get(i, j);
            EXPECT_DOUBLE_EQ(val, val2);
            // std::cerr << TlUtils::format("M(%d, %d)[%d]", i, j, subunit);
        }
    }

    // destroy
    for (int id = 0; id < numOfSubunits; ++id) {
        delete pMatrices[id];
        pMatrices[id] = NULL;
    }
}

TEST(TlDenseGeneralMatrix_arrays_mmap_RowOriented, reserve_resize) {
    const int numOfRows = 50;
    const int numOfCols = 10;
    const int numOfSubunits = 1;
    cleanup(matPath1, numOfSubunits);

    // construct
    std::vector<TlDenseGeneralMatrix_arrays_mmap_RowOriented*> pMatrices(numOfSubunits, NULL);
    for (int id = 0; id < numOfSubunits; ++id) {
        pMatrices[id] =
            new TlDenseGeneralMatrix_arrays_mmap_RowOriented(matPath1, numOfRows, numOfCols, numOfSubunits, id);
    }
    TlDenseGeneralMatrix_Lapack refMat(numOfRows, numOfCols);

    for (int id = 0; id < numOfSubunits; ++id) {
        EXPECT_EQ(numOfRows, pMatrices[id]->getNumOfRows());
        EXPECT_EQ(numOfCols, pMatrices[id]->getNumOfCols());
    }

    // prepare
    double val = 0.0;
    for (int i = 0; i < numOfRows; ++i) {
        for (int j = 0; j < numOfCols; ++j) {
            for (int id = 0; id < numOfSubunits; ++id) {
                pMatrices[id]->set(i, j, val);
            }
            refMat.set(i, j, val);

            val += 0.01;
        }
    }

    // reserve
    for (int id = 0; id < numOfSubunits; ++id) {
        pMatrices[id]->reserveColSize(numOfCols + 10);
    }

    // check
    for (int i = 0; i < numOfRows; ++i) {
        for (int j = 0; j < numOfCols; ++j) {
            const double val = refMat.get(i, j);

            const int subunit = pMatrices[0]->getSubunitID(i);
            for (int id = 0; id < numOfSubunits; ++id) {
                const int subunit2 = pMatrices[id]->getSubunitID(i);
                EXPECT_EQ(subunit, subunit2);
            }
            const double val2 = pMatrices[subunit]->get(i, j);
            EXPECT_DOUBLE_EQ(val, val2);
            // std::cerr << TlUtils::format("M(%d, %d)[%d] %f, %f", i, j, subunit, val, val2) << std::endl;
        }
    }

    // resize
    const int newNumOfRows = numOfRows;
    const int newNumOfCols = numOfCols + 5;
    for (int id = 0; id < numOfSubunits; ++id) {
        pMatrices[id]->resize(newNumOfRows, newNumOfCols);
    }
    refMat.resize(newNumOfRows, newNumOfCols);

    // check
    for (int i = 0; i < newNumOfRows; ++i) {
        for (int j = 0; j < newNumOfCols; ++j) {
            const double val = refMat.get(i, j);

            const int subunit = pMatrices[0]->getSubunitID(i);
            for (int id = 0; id < numOfSubunits; ++id) {
                const int subunit2 = pMatrices[id]->getSubunitID(i);
                EXPECT_EQ(subunit, subunit2);
            }
            const double val2 = pMatrices[subunit]->get(i, j);
            EXPECT_DOUBLE_EQ(val, val2);
            // std::cerr << TlUtils::format("M(%d, %d)[%d] %f, %f", i, j, subunit, val, val2) << std::endl;
        }
    }

    // destroy
    for (int id = 0; id < numOfSubunits; ++id) {
        delete pMatrices[id];
        pMatrices[id] = NULL;
    }
}

TEST(TlDenseGeneralMatrix_arrays_mmap_RowOriented, reserve_resize_group) {
    const int numOfRows = 1000;
    const int numOfCols = 100;
    const int numOfSubunits = 4;
    cleanup(matPath1, numOfSubunits);

    // construct
    std::vector<TlDenseGeneralMatrix_arrays_mmap_RowOriented*> pMatrices(numOfSubunits, NULL);
    for (int id = 0; id < numOfSubunits; ++id) {
        pMatrices[id] =
            new TlDenseGeneralMatrix_arrays_mmap_RowOriented(matPath1, numOfRows, numOfCols, numOfSubunits, id);
    }
    TlDenseGeneralMatrix_Lapack refMat(numOfRows, numOfCols);

    for (int id = 0; id < numOfSubunits; ++id) {
        EXPECT_EQ(numOfRows, pMatrices[id]->getNumOfRows());
        EXPECT_EQ(numOfCols, pMatrices[id]->getNumOfCols());
    }

    // prepare
    double val = 0.0;
    for (int i = 0; i < numOfRows; ++i) {
        for (int j = 0; j < numOfCols; ++j) {
            for (int id = 0; id < numOfSubunits; ++id) {
                pMatrices[id]->set(i, j, val);
            }
            refMat.set(i, j, val);

            val += 0.01;
        }
    }

    // reserve
    for (int id = 0; id < numOfSubunits; ++id) {
        pMatrices[id]->reserveColSize(numOfCols + 500);
    }

    // check
    for (int i = 0; i < numOfRows; ++i) {
        for (int j = 0; j < numOfCols; ++j) {
            const double val = refMat.get(i, j);

            const int subunit = pMatrices[0]->getSubunitID(i);
            for (int id = 0; id < numOfSubunits; ++id) {
                const int subunit2 = pMatrices[id]->getSubunitID(i);
                EXPECT_EQ(subunit, subunit2);
            }
            const double val2 = pMatrices[subunit]->get(i, j);
            EXPECT_DOUBLE_EQ(val, val2);
            // std::cerr << TlUtils::format("M(%d, %d)[%d] %f, %f", i, j, subunit, val, val2) << std::endl;
        }
    }

    // resize
    const int newNumOfRows = numOfRows;
    const int newNumOfCols = numOfCols + 200;
    for (int id = 0; id < numOfSubunits; ++id) {
        pMatrices[id]->resize(newNumOfRows, newNumOfCols);
    }
    refMat.resize(newNumOfRows, newNumOfCols);

    // check
    for (int i = 0; i < newNumOfRows; ++i) {
        for (int j = 0; j < newNumOfCols; ++j) {
            const double val = refMat.get(i, j);

            const int subunit = pMatrices[0]->getSubunitID(i);
            for (int id = 0; id < numOfSubunits; ++id) {
                const int subunit2 = pMatrices[id]->getSubunitID(i);
                EXPECT_EQ(subunit, subunit2);
            }
            const double val2 = pMatrices[subunit]->get(i, j);
            EXPECT_DOUBLE_EQ(val, val2);
            // std::cerr << TlUtils::format("M(%d, %d)[%d] %f, %f", i, j, subunit, val, val2) << std::endl;
        }
    }

    // destroy
    for (int id = 0; id < numOfSubunits; ++id) {
        delete pMatrices[id];
        pMatrices[id] = NULL;
    }
}

// TEST(TlDenseGeneralMatrix_arrays_mmap_RowOriented, resize_group1) {
//     const int numOfRows = 1000;
//     const int numOfCols = 2000;
//     const int numOfSubunits = 1;
//     std::vector<TlDenseGeneralMatrix_arrays_mmap_RowOriented*> pMatrices(numOfSubunits, NULL);
//     TlDenseGeneralMatrix_Lapack refMat(numOfRows, numOfCols);

//     // construct
//     for (int id = 0; id < numOfSubunits; ++id) {
//         pMatrices[id] = new TlDenseGeneralMatrix_arrays_mmap_RowOriented(numOfRows, numOfCols, numOfSubunits, id);
//     }

//     for (int id = 0; id < numOfSubunits; ++id) {
//         EXPECT_EQ(numOfRows, pMatrices[id]->getNumOfRows());
//         EXPECT_EQ(numOfCols, pMatrices[id]->getNumOfCols());
//     }

//     // prepare
//     double val = 0.0;
//     for (int i = 0; i < numOfRows; ++i) {
//         for (int j = 0; j < numOfCols; ++j) {
//             for (int id = 0; id < numOfSubunits; ++id) {
//                 pMatrices[id]->set(i, j, val);
//             }
//             refMat.set(i, j, val);
//         }
//         val += 0.01;
//     }

//     // resize
//     const int newNumOfRows = 4000;
//     const int newNumOfCols = 3000;
//     for (int id = 0; id < numOfSubunits; ++id) {
//         pMatrices[id]->resize(newNumOfRows, newNumOfCols);
//     }
//     refMat.resize(newNumOfRows, newNumOfCols);

//     for (int id = 0; id < numOfSubunits; ++id) {
//         EXPECT_EQ(newNumOfRows, pMatrices[id]->getNumOfRows());
//         EXPECT_EQ(newNumOfCols, pMatrices[id]->getNumOfCols());
//     }

//     // check
//     for (int i = 0; i < newNumOfRows; ++i) {
//         for (int j = 0; j < newNumOfCols; ++j) {
//             const double val = refMat.get(i, j);

//             const int subunit = pMatrices[0]->getSubunitID(i);
//             for (int id = 0; id < numOfSubunits; ++id) {
//                 const int subunit2 = pMatrices[id]->getSubunitID(i);
//                 EXPECT_EQ(subunit, subunit2);
//             }
//             const double val2 = pMatrices[subunit]->get(i, j);
//             EXPECT_DOUBLE_EQ(val, val2) << TlUtils::format("M(%d, %d)[%d]", i, j, subunit);
//         }
//     }

//     // destroy
//     for (int id = 0; id < numOfSubunits; ++id) {
//         delete pMatrices[id];
//         pMatrices[id] = NULL;
//     }
// }

// TEST(TlDenseGeneralMatrix_arrays_mmap_RowOriented, resize_row_group) {
//     const int numOfRows = 1000;
//     const int numOfCols = 2000;
//     const int numOfSubunits = 4;
//     std::vector<TlDenseGeneralMatrix_arrays_mmap_RowOriented*> pMatrices(numOfSubunits, NULL);
//     TlDenseGeneralMatrix_Lapack refMat(numOfRows, numOfCols);

//     // construct
//     for (int id = 0; id < numOfSubunits; ++id) {
//         pMatrices[id] = new TlDenseGeneralMatrix_arrays_mmap_RowOriented(numOfRows, numOfCols, numOfSubunits, id);
//     }

//     for (int id = 0; id < numOfSubunits; ++id) {
//         EXPECT_EQ(numOfRows, pMatrices[id]->getNumOfRows());
//         EXPECT_EQ(numOfCols, pMatrices[id]->getNumOfCols());
//     }

//     // prepare
//     double val = 0.0;
//     for (int i = 0; i < numOfRows; ++i) {
//         for (int j = 0; j < numOfCols; ++j) {
//             for (int id = 0; id < numOfSubunits; ++id) {
//                 pMatrices[id]->set(i, j, val);
//             }
//             refMat.set(i, j, val);
//         }
//         val += 0.01;
//     }

//     // resize
//     const int newNumOfRows = 2000;
//     const int newNumOfCols = numOfCols;
//     for (int id = 0; id < numOfSubunits; ++id) {
//         pMatrices[id]->resize(newNumOfRows, newNumOfCols);
//     }
//     refMat.resize(newNumOfRows, newNumOfCols);

//     for (int id = 0; id < numOfSubunits; ++id) {
//         EXPECT_EQ(newNumOfRows, pMatrices[id]->getNumOfRows());
//         EXPECT_EQ(newNumOfCols, pMatrices[id]->getNumOfCols());
//     }

//     // check
//     for (int i = 0; i < newNumOfRows; ++i) {
//         for (int j = 0; j < newNumOfCols; ++j) {
//             const double val = refMat.get(i, j);

//             const int subunit = pMatrices[0]->getSubunitID(i);
//             for (int id = 0; id < numOfSubunits; ++id) {
//                 const int subunit2 = pMatrices[id]->getSubunitID(i);
//                 EXPECT_EQ(subunit, subunit2);
//             }
//             const double val2 = pMatrices[subunit]->get(i, j);
//             EXPECT_DOUBLE_EQ(val, val2);
//         }
//     }

//     // destroy
//     for (int id = 0; id < numOfSubunits; ++id) {
//         delete pMatrices[id];
//         pMatrices[id] = NULL;
//     }
// }

// TEST(TlDenseGeneralMatrix_arrays_mmap_RowOriented, resize_col_group) {
//     const int numOfRows = 1000;
//     const int numOfCols = 2000;
//     const int numOfSubunits = 4;
//     std::vector<TlDenseGeneralMatrix_arrays_mmap_RowOriented*> pMatrices(numOfSubunits, NULL);
//     TlDenseGeneralMatrix_Lapack refMat(numOfRows, numOfCols);

//     // construct
//     for (int id = 0; id < numOfSubunits; ++id) {
//         pMatrices[id] = new TlDenseGeneralMatrix_arrays_mmap_RowOriented(numOfRows, numOfCols, numOfSubunits, id);
//     }

//     for (int id = 0; id < numOfSubunits; ++id) {
//         EXPECT_EQ(numOfRows, pMatrices[id]->getNumOfRows());
//         EXPECT_EQ(numOfCols, pMatrices[id]->getNumOfCols());
//     }

//     // prepare
//     double val = 0.0;
//     for (int i = 0; i < numOfRows; ++i) {
//         for (int j = 0; j < numOfCols; ++j) {
//             for (int id = 0; id < numOfSubunits; ++id) {
//                 pMatrices[id]->set(i, j, val);
//             }
//             refMat.set(i, j, val);
//         }
//         val += 0.01;
//     }

//     // resize
//     const int newNumOfRows = numOfRows;
//     const int newNumOfCols = 3000;
//     for (int id = 0; id < numOfSubunits; ++id) {
//         pMatrices[id]->resize(newNumOfRows, newNumOfCols);
//     }
//     refMat.resize(newNumOfRows, newNumOfCols);

//     for (int id = 0; id < numOfSubunits; ++id) {
//         EXPECT_EQ(newNumOfRows, pMatrices[id]->getNumOfRows());
//         EXPECT_EQ(newNumOfCols, pMatrices[id]->getNumOfCols());
//     }

//     // check
//     for (int i = 0; i < newNumOfRows; ++i) {
//         for (int j = 0; j < newNumOfCols; ++j) {
//             const double val = refMat.get(i, j);

//             const int subunit = pMatrices[0]->getSubunitID(i);
//             for (int id = 0; id < numOfSubunits; ++id) {
//                 const int subunit2 = pMatrices[id]->getSubunitID(i);
//                 EXPECT_EQ(subunit, subunit2);
//             }
//             const double val2 = pMatrices[subunit]->get(i, j);
//             EXPECT_DOUBLE_EQ(val, val2);
//         }
//     }

//     // destroy
//     for (int id = 0; id < numOfSubunits; ++id) {
//         delete pMatrices[id];
//         pMatrices[id] = NULL;
//     }
// }

// TEST(TlDenseGeneralMatrix_arrays_mmap_RowOriented, resize_group) {
//     const int numOfRows = 1000;
//     const int numOfCols = 2000;
//     const int numOfSubunits = 4;
//     std::vector<TlDenseGeneralMatrix_arrays_mmap_RowOriented*> pMatrices(numOfSubunits, NULL);
//     TlDenseGeneralMatrix_Lapack refMat(numOfRows, numOfCols);

//     // construct
//     for (int id = 0; id < numOfSubunits; ++id) {
//         pMatrices[id] = new TlDenseGeneralMatrix_arrays_mmap_RowOriented(numOfRows, numOfCols, numOfSubunits, id);
//     }

//     for (int id = 0; id < numOfSubunits; ++id) {
//         EXPECT_EQ(numOfRows, pMatrices[id]->getNumOfRows());
//         EXPECT_EQ(numOfCols, pMatrices[id]->getNumOfCols());
//     }

//     // prepare
//     double val = 0.0;
//     for (int i = 0; i < numOfRows; ++i) {
//         for (int j = 0; j < numOfCols; ++j) {
//             for (int id = 0; id < numOfSubunits; ++id) {
//                 pMatrices[id]->set(i, j, val);
//             }
//             refMat.set(i, j, val);
//         }
//         val += 0.01;
//     }

//     // resize
//     const int newNumOfRows = 4000;
//     const int newNumOfCols = 3000;
//     for (int id = 0; id < numOfSubunits; ++id) {
//         pMatrices[id]->resize(newNumOfRows, newNumOfCols);
//     }
//     refMat.resize(newNumOfRows, newNumOfCols);

//     for (int id = 0; id < numOfSubunits; ++id) {
//         EXPECT_EQ(newNumOfRows, pMatrices[id]->getNumOfRows());
//         EXPECT_EQ(newNumOfCols, pMatrices[id]->getNumOfCols());
//     }

//     // check
//     for (int i = 0; i < newNumOfRows; ++i) {
//         for (int j = 0; j < newNumOfCols; ++j) {
//             const double val = refMat.get(i, j);

//             const int subunit = pMatrices[0]->getSubunitID(i);
//             for (int id = 0; id < numOfSubunits; ++id) {
//                 const int subunit2 = pMatrices[id]->getSubunitID(i);
//                 EXPECT_EQ(subunit, subunit2);
//             }
//             const double val2 = pMatrices[subunit]->get(i, j);
//             EXPECT_DOUBLE_EQ(val, val2) << TlUtils::format("M(%d, %d)[%d]", i, j, subunit);
//         }
//     }

//     // destroy
//     for (int id = 0; id < numOfSubunits; ++id) {
//         delete pMatrices[id];
//         pMatrices[id] = NULL;
//     }
// }

// TEST(TlDenseGeneralMatrix_arrays_mmap_RowOriented, resize_group_decrese) {
//     const int numOfRows = 1000;
//     const int numOfCols = 2000;
//     const int numOfSubunits = 4;
//     std::vector<TlDenseGeneralMatrix_arrays_mmap_RowOriented*> pMatrices(numOfSubunits, NULL);
//     TlDenseGeneralMatrix_Lapack refMat(numOfRows, numOfCols);

//     // construct
//     for (int id = 0; id < numOfSubunits; ++id) {
//         pMatrices[id] = new TlDenseGeneralMatrix_arrays_mmap_RowOriented(numOfRows, numOfCols, numOfSubunits, id);
//     }

//     for (int id = 0; id < numOfSubunits; ++id) {
//         EXPECT_EQ(numOfRows, pMatrices[id]->getNumOfRows());
//         EXPECT_EQ(numOfCols, pMatrices[id]->getNumOfCols());
//     }

//     // prepare
//     double val = 0.0;
//     for (int i = 0; i < numOfRows; ++i) {
//         for (int j = 0; j < numOfCols; ++j) {
//             for (int id = 0; id < numOfSubunits; ++id) {
//                 pMatrices[id]->set(i, j, val);
//             }
//             refMat.set(i, j, val);
//         }
//         val += 0.01;
//     }

//     // resize
//     const int newNumOfRows = 500;
//     const int newNumOfCols = 1000;
//     for (int id = 0; id < numOfSubunits; ++id) {
//         pMatrices[id]->resize(newNumOfRows, newNumOfCols);
//     }
//     refMat.resize(newNumOfRows, newNumOfCols);

//     for (int id = 0; id < numOfSubunits; ++id) {
//         EXPECT_EQ(newNumOfRows, pMatrices[id]->getNumOfRows());
//         EXPECT_EQ(newNumOfCols, pMatrices[id]->getNumOfCols());
//     }

//     // check
//     for (int i = 0; i < newNumOfRows; ++i) {
//         for (int j = 0; j < newNumOfCols; ++j) {
//             const double val = refMat.get(i, j);

//             const int subunit = pMatrices[0]->getSubunitID(i);
//             for (int id = 0; id < numOfSubunits; ++id) {
//                 const int subunit2 = pMatrices[id]->getSubunitID(i);
//                 EXPECT_EQ(subunit, subunit2);
//             }
//             const double val2 = pMatrices[subunit]->get(i, j);
//             EXPECT_DOUBLE_EQ(val, val2) << TlUtils::format("M(%d, %d)[%d]", i, j, subunit);
//         }
//     }

//     // destroy
//     for (int id = 0; id < numOfSubunits; ++id) {
//         delete pMatrices[id];
//         pMatrices[id] = NULL;
//     }
// }

// TEST(TlDenseGeneralMatrix_arrays_mmap_RowOriented, toTlMatrix) {
//     const int maxRow = 100;
//     const int maxCol = 80;
//     TlDenseGeneralMatrix_arrays_mmap_RowOriented vecA(maxRow, maxCol);
//     TlDenseGeneralMatrix_Lapack matA(maxRow, maxCol);

//     // setup
//     int count = 0;
//     for (int r = 0; r < maxRow; ++r) {
//         for (int c = 0; c < maxCol; ++c) {
//             double v = double(count);
//             matA.set(r, c, v);
//             vecA.set(r, c, v);

//             ++count;
//         }
//     }

//     TlDenseGeneralMatrix_Lapack matB = vecA.getTlMatrixObject();

//     // test
//     for (int r = 0; r < maxRow; ++r) {
//         for (int c = 0; c < maxCol; ++c) {
//             EXPECT_NEAR(matA.get(r, c), matB.get(r, c), EPS);
//         }
//     }
// }

TEST(TlDenseGeneralMatrix_arrays_mmap_RowOriented, RowVectorMatrix2CSFD) {
    const int numOfRows = 1000;
    const int numOfCols = 2000;
    const int numOfSubunits = 1;
    cleanup(matPath1, numOfSubunits);

    std::vector<TlDenseGeneralMatrix_arrays_mmap_RowOriented*> pMatrices(numOfSubunits, NULL);
    TlDenseGeneralMatrix_Lapack refMat(numOfRows, numOfCols);

    // construct
    for (int id = 0; id < numOfSubunits; ++id) {
        pMatrices[id] =
            new TlDenseGeneralMatrix_arrays_mmap_RowOriented(matPath1, numOfRows, numOfCols, numOfSubunits, id);
    }

    // prepare
    double val = 0.0;
    for (int i = 0; i < numOfRows; ++i) {
        for (int j = 0; j < numOfCols; ++j) {
            for (int id = 0; id < numOfSubunits; ++id) {
                pMatrices[id]->set(i, j, val);
            }
            refMat.set(i, j, val);
        }
        val += 0.01;
    }

    // transform
    std::cerr << "transform begin." << std::endl;
    transpose2CSFD(matPath1, csfdPath, true);
    std::cerr << "transform done." << std::endl;

    // load
    TlDenseGeneralMatrix_Lapack chkMat;
    chkMat.load(csfdPath);

    // check
    EXPECT_EQ(chkMat.getNumOfRows(), numOfRows);
    EXPECT_EQ(chkMat.getNumOfCols(), numOfCols);
    for (int i = 0; i < numOfRows; ++i) {
        for (int j = 0; j < numOfCols; ++j) {
            const double val = refMat.get(i, j);
            const double val2 = chkMat.get(i, j);
            EXPECT_DOUBLE_EQ(val, val2);
        }
    }

    // destroy
    for (int id = 0; id < numOfSubunits; ++id) {
        delete pMatrices[id];
        pMatrices[id] = NULL;
    }
}
