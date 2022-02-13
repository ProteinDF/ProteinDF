#include "tl_dense_general_matrix_arrays_roworiented.h"

#include <iostream>
#include <limits>

#include "TlUtils.h"
#include "gtest/gtest.h"
#include "tl_dense_general_matrix_arrays_coloriented.h"
#include "tl_dense_general_matrix_lapack.h"

static const double EPS = 1.0E-10;  // std::numeric_limits<double>::epsilon();

TEST(TlDenseGeneralMatrix_arrays_RowOriented, constructer) {
    TlDenseGeneralMatrix_arrays_RowOriented A(100, 600);

    EXPECT_EQ(100, A.getNumOfRows());
    EXPECT_EQ(600, A.getNumOfCols());
    EXPECT_NEAR(0.0, A.get(0, 0), EPS);
}

TEST(TlDenseGeneralMatrix_arrays_RowOriented, constructer2) {
    TlDenseGeneralMatrix_arrays_RowOriented A0(100, 600, 10, 0);
    TlDenseGeneralMatrix_arrays_RowOriented A1(100, 600, 10, 1);

    EXPECT_EQ(100, A0.getNumOfRows());
    EXPECT_EQ(600, A0.getNumOfCols());
    EXPECT_EQ(10, A0.getNumOfSubunits());
    EXPECT_EQ(0, A0.getSubunitID());

    EXPECT_EQ(100, A1.getNumOfRows());
    EXPECT_EQ(600, A1.getNumOfCols());
    EXPECT_EQ(10, A1.getNumOfSubunits());
    EXPECT_EQ(1, A1.getSubunitID());
}

TEST(TlDenseGeneralMatrix_arrays_RowOriented, resize) {
    TlDenseGeneralMatrix_arrays_RowOriented A(100, 100);
    A.resize(200, 100);
    EXPECT_EQ(200, A.getNumOfRows());
    EXPECT_EQ(100, A.getNumOfCols());

    TlDenseGeneralMatrix_arrays_RowOriented B(100, 100);
    B.resize(100, 200);
    EXPECT_EQ(100, B.getNumOfRows());
    EXPECT_EQ(200, B.getNumOfCols());

    TlDenseGeneralMatrix_arrays_RowOriented C(50, 100);
    C.resize(50, 100);
    EXPECT_EQ(50, C.getNumOfRows());
    EXPECT_EQ(100, C.getNumOfCols());

    TlDenseGeneralMatrix_arrays_RowOriented D(100, 50);
    D.resize(100, 50);
    EXPECT_EQ(100, D.getNumOfRows());
    EXPECT_EQ(50, D.getNumOfCols());
}

TEST(TlDenseGeneralMatrix_arrays_RowOriented, contents) {
    const int maxRow = 100;
    const int maxCol = 80;
    TlDenseGeneralMatrix_arrays_RowOriented vecA(maxRow, maxCol);
    TlDenseGeneralMatrix_Lapack matA(maxRow, maxCol);

    // setup
    int count = 0;
    for (int r = 0; r < maxRow; ++r) {
        for (int c = 0; c < maxCol; ++c) {
            double v = double(count);
            matA.set(r, c, v);
            vecA.set(r, c, v);

            ++count;
        }
    }

    // test
    for (int r = 0; r < maxRow; ++r) {
        for (int c = 0; c < maxCol; ++c) {
            EXPECT_NEAR(matA.get(r, c), vecA.get(r, c), EPS);
        }
    }
}

TEST(TlDenseGeneralMatrix_arrays_RowOriented, contents_group) {
    const int numOfRows = 1000;
    const int numOfCols = 2000;
    const int numOfSubunits = 4;
    std::vector<TlDenseGeneralMatrix_arrays_RowOriented*> pMatrices(
        numOfSubunits, NULL);
    TlDenseGeneralMatrix_Lapack refMat(numOfRows, numOfCols);

    // construct
    for (int id = 0; id < numOfSubunits; ++id) {
        pMatrices[id] = new TlDenseGeneralMatrix_arrays_RowOriented(
            numOfRows, numOfCols, numOfSubunits, id);
    }

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
        }
        val += 0.01;
    }

    // resize
    const int newNumOfRows = numOfRows;
    const int newNumOfCols = numOfCols;
    // for (int id = 0; id < numOfSubunits; ++id) {
    //     pMatrices[id]->resize(newNumOfRows, newNumOfCols);
    // }
    // refMat.resize(newNumOfRows, newNumOfCols);

    // for (int id = 0; id < numOfSubunits; ++id) {
    //     EXPECT_EQ(newNumOfRows, pMatrices[id]->getNumOfRows());
    //     EXPECT_EQ(newNumOfCols, pMatrices[id]->getNumOfCols());
    // }

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
        }
    }

    // destroy
    for (int id = 0; id < numOfSubunits; ++id) {
        delete pMatrices[id];
        pMatrices[id] = NULL;
    }
}

TEST(TlDenseGeneralMatrix_arrays_RowOriented, reserve1) {
    const int numOfRows = 1000;
    const int numOfCols = 2000;
    const int numOfSubunits = 1;
    std::vector<TlDenseGeneralMatrix_arrays_RowOriented*> pMatrices(
        numOfSubunits, NULL);
    TlDenseGeneralMatrix_Lapack refMat(numOfRows, numOfCols);

    // construct
    for (int id = 0; id < numOfSubunits; ++id) {
        pMatrices[id] = new TlDenseGeneralMatrix_arrays_RowOriented(
            numOfRows, numOfCols, numOfSubunits, id);
    }

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
        }
        val += 0.01;
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
            EXPECT_DOUBLE_EQ(val, val2)
                << TlUtils::format("M(%d, %d)[%d]", i, j, subunit);
        }
    }

    // destroy
    for (int id = 0; id < numOfSubunits; ++id) {
        delete pMatrices[id];
        pMatrices[id] = NULL;
    }
}

TEST(TlDenseGeneralMatrix_arrays_RowOriented, resize_group1) {
    const int numOfRows = 1000;
    const int numOfCols = 2000;
    const int numOfSubunits = 1;
    std::vector<TlDenseGeneralMatrix_arrays_RowOriented*> pMatrices(
        numOfSubunits, NULL);
    TlDenseGeneralMatrix_Lapack refMat(numOfRows, numOfCols);

    // construct
    for (int id = 0; id < numOfSubunits; ++id) {
        pMatrices[id] = new TlDenseGeneralMatrix_arrays_RowOriented(
            numOfRows, numOfCols, numOfSubunits, id);
    }

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
        }
        val += 0.01;
    }

    // resize
    const int newNumOfRows = 4000;
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
            EXPECT_DOUBLE_EQ(val, val2)
                << TlUtils::format("M(%d, %d)[%d]", i, j, subunit);
        }
    }

    // destroy
    for (int id = 0; id < numOfSubunits; ++id) {
        delete pMatrices[id];
        pMatrices[id] = NULL;
    }
}

TEST(TlDenseGeneralMatrix_arrays_RowOriented, resize_row_group) {
    const int numOfRows = 1000;
    const int numOfCols = 2000;
    const int numOfSubunits = 4;
    std::vector<TlDenseGeneralMatrix_arrays_RowOriented*> pMatrices(
        numOfSubunits, NULL);
    TlDenseGeneralMatrix_Lapack refMat(numOfRows, numOfCols);

    // construct
    for (int id = 0; id < numOfSubunits; ++id) {
        pMatrices[id] = new TlDenseGeneralMatrix_arrays_RowOriented(
            numOfRows, numOfCols, numOfSubunits, id);
    }

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
        }
        val += 0.01;
    }

    // resize
    const int newNumOfRows = 2000;
    const int newNumOfCols = numOfCols;
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
            EXPECT_DOUBLE_EQ(val, val2);
        }
    }

    // destroy
    for (int id = 0; id < numOfSubunits; ++id) {
        delete pMatrices[id];
        pMatrices[id] = NULL;
    }
}

TEST(TlDenseGeneralMatrix_arrays_RowOriented, resize_col_group) {
    const int numOfRows = 1000;
    const int numOfCols = 2000;
    const int numOfSubunits = 4;
    std::vector<TlDenseGeneralMatrix_arrays_RowOriented*> pMatrices(
        numOfSubunits, NULL);
    TlDenseGeneralMatrix_Lapack refMat(numOfRows, numOfCols);

    // construct
    for (int id = 0; id < numOfSubunits; ++id) {
        pMatrices[id] = new TlDenseGeneralMatrix_arrays_RowOriented(
            numOfRows, numOfCols, numOfSubunits, id);
    }

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
        }
        val += 0.01;
    }

    // resize
    const int newNumOfRows = numOfRows;
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
            EXPECT_DOUBLE_EQ(val, val2);
        }
    }

    // destroy
    for (int id = 0; id < numOfSubunits; ++id) {
        delete pMatrices[id];
        pMatrices[id] = NULL;
    }
}

TEST(TlDenseGeneralMatrix_arrays_RowOriented, resize_group) {
    const int numOfRows = 1000;
    const int numOfCols = 2000;
    const int numOfSubunits = 4;
    std::vector<TlDenseGeneralMatrix_arrays_RowOriented*> pMatrices(
        numOfSubunits, NULL);
    TlDenseGeneralMatrix_Lapack refMat(numOfRows, numOfCols);

    // construct
    for (int id = 0; id < numOfSubunits; ++id) {
        pMatrices[id] = new TlDenseGeneralMatrix_arrays_RowOriented(
            numOfRows, numOfCols, numOfSubunits, id);
    }

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
        }
        val += 0.01;
    }

    // resize
    const int newNumOfRows = 4000;
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
            EXPECT_DOUBLE_EQ(val, val2)
                << TlUtils::format("M(%d, %d)[%d]", i, j, subunit);
        }
    }

    // destroy
    for (int id = 0; id < numOfSubunits; ++id) {
        delete pMatrices[id];
        pMatrices[id] = NULL;
    }
}

TEST(TlDenseGeneralMatrix_arrays_RowOriented, resize_group_decrese) {
    const int numOfRows = 1000;
    const int numOfCols = 2000;
    const int numOfSubunits = 4;
    std::vector<TlDenseGeneralMatrix_arrays_RowOriented*> pMatrices(
        numOfSubunits, NULL);
    TlDenseGeneralMatrix_Lapack refMat(numOfRows, numOfCols);

    // construct
    for (int id = 0; id < numOfSubunits; ++id) {
        pMatrices[id] = new TlDenseGeneralMatrix_arrays_RowOriented(
            numOfRows, numOfCols, numOfSubunits, id);
    }

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
        }
        val += 0.01;
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
            EXPECT_DOUBLE_EQ(val, val2)
                << TlUtils::format("M(%d, %d)[%d]", i, j, subunit);
        }
    }

    // destroy
    for (int id = 0; id < numOfSubunits; ++id) {
        delete pMatrices[id];
        pMatrices[id] = NULL;
    }
}

TEST(TlDenseGeneralMatrix_arrays_RowOriented, save_load) {
    const int maxRow = 100;
    const int maxCol = 80;
    TlDenseGeneralMatrix_arrays_RowOriented vecA(maxRow, maxCol);
    TlDenseGeneralMatrix_Lapack matA(maxRow, maxCol);

    // setup
    int count = 0;
    for (int r = 0; r < maxRow; ++r) {
        for (int c = 0; c < maxCol; ++c) {
            double v = double(count);
            matA.set(r, c, v);
            vecA.set(r, c, v);

            ++count;
        }
    }

    // save
    vecA.save("/tmp/vecA.mat");

    // load
    TlDenseGeneralMatrix_arrays_RowOriented vecB;
    vecB.load("/tmp/vecA.mat", 0);

    // test
    for (int r = 0; r < maxRow; ++r) {
        for (int c = 0; c < maxCol; ++c) {
            EXPECT_NEAR(matA.get(r, c), vecB.get(r, c), EPS);
        }
    }
}

TEST(TlDenseGeneralMatrix_arrays_RowOriented, toTlMatrix) {
    const int maxRow = 100;
    const int maxCol = 80;
    TlDenseGeneralMatrix_arrays_RowOriented vecA(maxRow, maxCol);
    TlDenseGeneralMatrix_Lapack matA(maxRow, maxCol);

    // setup
    int count = 0;
    for (int r = 0; r < maxRow; ++r) {
        for (int c = 0; c < maxCol; ++c) {
            double v = double(count);
            matA.set(r, c, v);
            vecA.set(r, c, v);

            ++count;
        }
    }

    TlDenseGeneralMatrix_Lapack matB = vecA.getTlMatrixObject();

    // test
    for (int r = 0; r < maxRow; ++r) {
        for (int c = 0; c < maxCol; ++c) {
            EXPECT_NEAR(matA.get(r, c), matB.get(r, c), EPS);
        }
    }
}

// TEST(TlDenseGeneralMatrix_arrays_RowOriented, saveAsColOriented1) {
//     const int numOfRows = 1000;
//     const int numOfCols = 2000;
//     const int numOfSubunits = 1;
//     std::vector<TlDenseGeneralMatrix_arrays_RowOriented*> pMatrices_row(numOfSubunits, NULL);
//     TlDenseGeneralMatrix_Lapack refMat(numOfRows, numOfCols);

//     // construct
//     for (int id = 0; id < numOfSubunits; ++id) {
//         pMatrices_row[id] = new TlDenseGeneralMatrix_arrays_RowOriented(numOfRows, numOfCols, numOfSubunits, id);
//     }

//     for (int id = 0; id < numOfSubunits; ++id) {
//         EXPECT_EQ(numOfRows, pMatrices_row[id]->getNumOfRows());
//         EXPECT_EQ(numOfCols, pMatrices_row[id]->getNumOfCols());
//     }

//     // prepare
//     double val = 0.0;
//     for (int i = 0; i < numOfRows; ++i) {
//         for (int j = 0; j < numOfCols; ++j) {
//             for (int id = 0; id < numOfSubunits; ++id) {
//                 pMatrices_row[id]->set(i, j, val);
//             }
//             refMat.set(i, j, val);
//         }
//         val += 0.01;
//     }

//     // transform
//     std::cerr << "saveAsColOriented1 transform" << std::endl;
//     const std::string basename = "/tmp/test_array_col";
//     for (int id = 0; id < numOfSubunits; ++id) {
//         pMatrices_row[id]->saveByTlDenseGeneralMatrix_arrays_ColOriented(basename);
//     }

//     // load
//     std::cerr << "saveAsColOriented1 load" << std::endl;
//     std::vector<TlDenseGeneralMatrix_arrays_ColOriented*> pMatrices_col(numOfSubunits, NULL);
//     for (int id = 0; id < numOfSubunits; ++id) {
//         pMatrices_col[id] = new TlDenseGeneralMatrix_arrays_ColOriented(1, 1, numOfSubunits, id);
//         pMatrices_col[id]->load(basename);
//     }

//     // check
//     std::cerr << "saveAsColOriented1 check" << std::endl;
//     for (int id = 0; id < numOfSubunits; ++id) {
//         EXPECT_EQ(pMatrices_col[id]->getNumOfRows(), numOfRows);
//         EXPECT_EQ(pMatrices_col[id]->getNumOfCols(), numOfCols);
//     }
//     for (int i = 0; i < numOfRows; ++i) {
//         for (int j = 0; j < numOfCols; ++j) {
//             const double val = refMat.get(i, j);

//             const int subunit = pMatrices_col[0]->getSubunitID(i);
//             for (int id = 0; id < numOfSubunits; ++id) {
//                 const int subunit2 = pMatrices_col[id]->getSubunitID(i);
//                 EXPECT_EQ(subunit, subunit2);
//             }
//             const double val2 = pMatrices_col[subunit]->get(i, j);
//             EXPECT_DOUBLE_EQ(val, val2);
//         }
//     }

//     // destroy
//     for (int id = 0; id < numOfSubunits; ++id) {
//         delete pMatrices_row[id];
//         pMatrices_row[id] = NULL;

//         delete pMatrices_col[id];
//         pMatrices_col[id] = NULL;
//     }
// }
