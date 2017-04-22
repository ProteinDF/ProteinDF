#include <cstdlib>
#include <iostream>
#include <string>
#include <stdlib.h>
#include <time.h>

#include "TlMatrix.h"
#include "TlSymmetricMatrix.h"
#include "TlDistributeSymmetricMatrix.h"
#include "TlCommunicate.h"

static const int MATRIX_SIZE = 100;

void showResultMessage(const std::string& sFunction, bool bIsPassed);
void showResultMessageAll(const std::string& sFunction, bool bIsPassed);

void makeMatrixA(TlDistributeSymmetricMatrix* pDistA,
                 TlSymmetricMatrix* pA);
void makeMatrixB(TlDistributeSymmetricMatrix* pDistB,
                 TlSymmetricMatrix* pB);

void testConstructer();
void testCopyConstructer();
void testCopyConstructer2();
void testSet();
void testOperatorPlusEqual();
void testOperatorMultiEqual();
void testMulti1();
void testMulti2();
void testMulti3();
void testSave();
void testLoad();
void testInverse();
void testInverse2();
void testDiagonal();
void testDot();
void testSum();
void testMergePartialMatrix();
// void testGetPartialMatrix();
void testGetSparseMatrix2();
void testGetPartialMatrix2();

int main(int argc, char *argv[])
{
    const int systemBlockSize = 4;

    // initialize
    TlCommunicate& rComm = TlCommunicate::getInstance(argc, argv);
    TlDistributeMatrix::setSystemBlockSize(systemBlockSize);

    // ===================================================================
    testConstructer();
    testSet();
    testCopyConstructer();
    testCopyConstructer2();
    testOperatorPlusEqual();
    testOperatorMultiEqual();
    testMulti1();
    testMulti2();
    testMulti3();
    testSave();
    testLoad();
    testDot();
    testSum();
    testInverse();
    testInverse2();
    testDiagonal();
    testMergePartialMatrix();
    // testGetPartialMatrix();
    testGetSparseMatrix2();
    testGetPartialMatrix2();
    // ===================================================================

    // finalize
    rComm.finalize();
    return EXIT_SUCCESS;
}

void showResultMessage(const std::string& sFunction, bool bIsPassed)
{
    TlCommunicate& rComm = TlCommunicate::getInstance();
    if (rComm.isMaster() == true) {
        std::cout << "TEST: " << sFunction << "() ";
        std::cout << ((bIsPassed == true) ? "." : "F") << std::endl;
    }
    rComm.barrier();
}

void showResultMessageAll(const std::string& sFunction, bool bIsPassed)
{
    TlCommunicate& rComm = TlCommunicate::getInstance();
    int nProc = rComm.getNumOfProc();
    int nRank = rComm.getRank();

    rComm.barrier();
    if (rComm.isMaster() == true) {
        std::cout << "TEST: " << sFunction << "() ";
    }
    rComm.barrier();
    for (int i = 0; i < nProc; ++i) {
        if (i == nRank) {
            std::cout << ((bIsPassed == true) ? "." : "F");
            std::cout.flush();
        }
        rComm.barrier();
    }
    rComm.barrier();
    if (rComm.isMaster() == true) {
        std::cout << std::endl;
    }

    rComm.barrier();
}


void makeMatrixA(TlDistributeSymmetricMatrix* pDistA,
                 TlSymmetricMatrix* pA)
{
    const int size = MATRIX_SIZE;
    pDistA->resize(size);
    pA->resize(size);

    double count = 1.0;
    for (int i = 0; i < size; ++i) {
        for (int j = 0; j <= i; ++j) {
            pDistA->set(i, j, count);
            pA->set(i, j, count);
            //count += 1.0;
        }
    }
}


void makeMatrixB(TlDistributeSymmetricMatrix* pDistB,
                 TlSymmetricMatrix* pB)
{
    srand((unsigned int)time(NULL));
    
    const int size = MATRIX_SIZE;
    pDistB->resize(size);
    pB->resize(size);

    for (int i = 0; i < size; ++i) {
        for (int j = 0; j <= i; ++j) {
            double v = (rand() / ((double)RAND_MAX+1.0)) * 100.0;
            pDistB->set(i, j, v);
            pB->set(i, j, v);
        }
    }
}


void testConstructer()
{
    TlDistributeSymmetricMatrix matrix(100);

    int nRow = matrix.getNumOfRows();
    int nCol = matrix.getNumOfCols();

    bool bIsPassed = ((nRow == 100) && (nCol == 100)) ? true : false;

    showResultMessageAll("testConstructer", bIsPassed);
}


void testSet()
{
    TlDistributeSymmetricMatrix A(25);
    A(10, 10) = 10.0;
    A(17,  4) = 56.0;
    A(21, 18) = -2.5;
    A(0, 10) = 12.0;

    // A.print(std::cout);

    bool bIsPassed = true;
    if ((A.getNumOfRows() != 25) || (A.getNumOfCols() != 25)) {
        bIsPassed = false;
    }
    if (std::fabs(A(10, 10) -  10.0) > 1.0E-16) {
        bIsPassed = false;
    }
    if (std::fabs(A(17,  4) -  56.0) > 1.0E-16) {
        bIsPassed = false;
    }
    if (std::fabs(A(4,  17) -  56.0) > 1.0E-16) {
        bIsPassed = false;
    }
    if (std::fabs(A(21, 18) - (-2.5)) > 1.0E-16) {
        bIsPassed = false;
    }
    if (std::fabs(A(18, 21) - (-2.5)) > 1.0E-16) {
        bIsPassed = false;
    }
    if (std::fabs(A(0, 10) -  12.0) > 1.0E-16) {
        bIsPassed = false;
    }
    if (std::fabs(A(10,  0) -  12.0) > 1.0E-16) {
        bIsPassed = false;
    }

    showResultMessageAll("testSet", bIsPassed);
}


void testCopyConstructer()
{
    TlDistributeSymmetricMatrix A(100);
    A(10, 10) = 10.0;
    A(17,  4) = 56.0;
    A(21, 18) = -2.5;
    A(0, 10) = 12.0;

    const TlDistributeSymmetricMatrix B(A);

    //A.print(std::cout);
    //B.print(std::cout);

    bool bIsPassed = true;
    if ((B.getNumOfRows() != 100) || (B.getNumOfCols() != 100)) {
        bIsPassed = false;
    }
    if (std::fabs(B(10, 10) -  10.0) > 1.0E-16) {
        bIsPassed = false;
        std::cout << "B(10, 10) = " << B(10, 10) << ", exact = 10.0" << std::endl;
    }
    if (std::fabs(B(17,  4) -  56.0) > 1.0E-16) {
        bIsPassed = false;
        std::cout << "B(17,  4) = " << B(17,  4) << ", exact = 56.0" << std::endl;
    }
    if (std::fabs(B(4,  17) -  56.0) > 1.0E-16) {
        bIsPassed = false;
        std::cout << "B( 4, 17) = " << B(4, 17) << ", exact = 56.0" << std::endl;
    }
    if (std::fabs(B(21, 18) - (-2.5)) > 1.0E-16) {
        bIsPassed = false;
        std::cout << "B(21, 18) = " << B(21, 18) << ", exact = -2.5" << std::endl;
    }
    if (std::fabs(B(18, 21) - (-2.5)) > 1.0E-16) {
        bIsPassed = false;
        std::cout << "B(18, 21) = " << B(18, 21) << ", exact = -2.5" << std::endl;
    }
    if (std::fabs(B(0, 10) -  12.0) > 1.0E-16) {
        bIsPassed = false;
        std::cout << "B(0, 10) = " << B(0, 10) << ", exact = 12.0" << std::endl;
    }
    if (std::fabs(B(10,  0) -  12.0) > 1.0E-16) {
        bIsPassed = false;
        std::cout << "B(10, 0) = " << B(10, 0) << ", exact = 12.0" << std::endl;
    }

    showResultMessageAll("testCopyConstructer", bIsPassed);
}


void testCopyConstructer2()
{
    TlDistributeSymmetricMatrix A(100);
    A(10, 10) = 10.0;
    A(17,  4) = 56.0;
    A(21, 18) = -2.5;
    A(0, 10) = 12.0;
    //A.print(std::cout);

    const TlDistributeMatrix B(A);
    //B.print(std::cout);

    bool bIsPassed = true;
    if ((B.getNumOfRows() != 100) || (B.getNumOfCols() != 100)) {
        bIsPassed = false;
    }
    if (std::fabs(B(10, 10) -  10.0) > 1.0E-16) {
        bIsPassed = false;
    }
    if (std::fabs(B(17,  4) -  56.0) > 1.0E-16) {
        bIsPassed = false;
    }
    if (std::fabs(B(4,  17) -  56.0) > 1.0E-16) {
        bIsPassed = false;
    }
    if (std::fabs(B(21, 18) - (-2.5)) > 1.0E-16) {
        bIsPassed = false;
    }
    if (std::fabs(B(18, 21) - (-2.5)) > 1.0E-16) {
        bIsPassed = false;
    }
    if (std::fabs(B(0, 10) -  12.0) > 1.0E-16) {
        bIsPassed = false;
    }
    if (std::fabs(B(10,  0) -  12.0) > 1.0E-16) {
        bIsPassed = false;
    }

    showResultMessageAll("testCopyConstructer2", bIsPassed);
}


void testOperatorPlusEqual()
{
    TlDistributeSymmetricMatrix A(25);
    A(10, 10) = 10.0;
    A(17,  4) = 56.0;
    A(21, 18) = -2.5;
    A(0, 10) = 12.0;

    TlDistributeSymmetricMatrix B(25);
    B(10, 10) =  5.0;
    B(0, 10) = -12.0;
    B(3, 17) =  51.0;

    A += B;

    bool bIsPassed = true;
    if ((A.getNumOfRows() != 25) || (A.getNumOfCols() != 25)) {
        bIsPassed = false;
    }
    if (std::fabs(A(10, 10) -  15.0) > 1.0E-16) {
        bIsPassed = false;
    }
    if (std::fabs(A(17,  4) -  56.0) > 1.0E-16) {
        bIsPassed = false;
    }
    if (std::fabs(A(4,  17) -  56.0) > 1.0E-16) {
        bIsPassed = false;
    }
    if (std::fabs(A(21, 18) - (-2.5)) > 1.0E-16) {
        bIsPassed = false;
    }
    if (std::fabs(A(18, 21) - (-2.5)) > 1.0E-16) {
        bIsPassed = false;
    }
    if (std::fabs(A(0, 10) -  0.0) > 1.0E-16) {
        bIsPassed = false;
    }
    if (std::fabs(A(10,  0) -  0.0) > 1.0E-16) {
        bIsPassed = false;
    }
    if (std::fabs(A(3, 17) -  51.0) > 1.0E-16) {
        bIsPassed = false;
    }
    if (std::fabs(A(17, 3) -  51.0) > 1.0E-16) {
        bIsPassed = false;
    }

    showResultMessageAll("testOperatorPlusEqual", bIsPassed);
}

void testOperatorMultiEqual()
{
    TlDistributeSymmetricMatrix A(25);
    A(10, 10) = 10.0;
    A(17,  4) = 56.0;
    A(21, 18) = -2.5;
    A(0, 10) = 12.0;

    A *= 3.0;

    bool bIsPassed = true;
    if ((A.getNumOfRows() != 25) || (A.getNumOfCols() != 25)) {
        bIsPassed = false;
    }
    if (std::fabs(A(10, 10) -  10.0*3.0) > 1.0E-16) {
        bIsPassed = false;
    }
    if (std::fabs(A(17,  4) -  56.0*3.0) > 1.0E-16) {
        bIsPassed = false;
    }
    if (std::fabs(A(4, 17) -  56.0*3.0) > 1.0E-16) {
        bIsPassed = false;
    }
    if (std::fabs(A(21, 18) - (-2.5)*3.0) > 1.0E-16) {
        bIsPassed = false;
    }
    if (std::fabs(A(18, 21) - (-2.5)*3.0) > 1.0E-16) {
        bIsPassed = false;
    }
    if (std::fabs(A(0, 10) -  12.0*3.0) > 1.0E-16) {
        bIsPassed = false;
    }
    if (std::fabs(A(10,  0) -  12.0*3.0) > 1.0E-16) {
        bIsPassed = false;
    }

    showResultMessageAll("testOperatorMultiEqual", bIsPassed);
}

void testMulti1()
{
    bool bIsPassed = true;
    const int dim = 100;
    
    TlDistributeSymmetricMatrix A(dim);
    TlSymmetricMatrix exactA(dim);
    {
        double count = 1.0;
        for (int i = 0; i < dim; ++i) {
            for (int j = 0; j <= i; ++j) {
                A(i, j) = count;
                exactA(i, j) = count;
                count += 1.0 / 100;
            }
        }
    }

    TlDistributeSymmetricMatrix B(dim);
    TlSymmetricMatrix exactB(dim);
    {
        double count = 1.82;
        for (int i = 0; i < dim; ++i) {
            for (int j = 0; j <= i; ++j) {
                B(i, j) = count;
                exactB(i, j) = count;
                count -= 1.0 / 100;
            }
        }
    }

    TlDistributeMatrix C = A * B;
    TlMatrix exactC = exactA * exactB;

    // judge
    for (int r = 0; r < dim; ++r) {
        for (int s = 0; s < dim; ++s) {
            double c = C(r, s);
            double exact_c = exactC(r, s);

            if (fabs(c - exact_c) > 1.0E-8) {
                bIsPassed = false;
                std::cout << TlUtils::format("(%2d, %2d): actual= %f, expect= %f, err= %e", r, s, c, exact_c, fabs(c - exact_c)) << std::endl;
            }
        }
    }

    showResultMessageAll("testMulti1", bIsPassed);
}


void testMulti2()
{
    bool bIsPassed = true;
    const int dim = 100;
    
    TlDistributeSymmetricMatrix A(dim);
    TlSymmetricMatrix exactA(dim);
    {
        double count = 1.0;
        for (int i = 0; i < dim; ++i) {
            for (int j = 0; j <= i; ++j) {
                A(i, j) = count;
                exactA(i, j) = count;
                count += 1.0 / 100;
            }
        }
    }

    TlDistributeMatrix B(dim, dim);
    TlMatrix exactB(dim, dim);
    {
        double count = 1.82;
        for (int i = 0; i < dim; ++i) {
            for (int j = 0; j < dim; ++j) {
                B(i, j) = count;
                exactB(i, j) = count;
                count -= 1.0 / 100;
            }
        }
    }

    TlDistributeMatrix C = A * B;
    TlMatrix exactC = exactA * exactB;

    // judge
    for (int r = 0; r < dim; ++r) {
        for (int s = 0; s < dim; ++s) {
            double c = C(r, s);
            double exact_c = exactC(r, s);

            if (fabs(c - exact_c) > 1.0E-8) {
                bIsPassed = false;
                std::cout << TlUtils::format("(%2d, %2d): actual= %f, expect= %f, err= %e", r, s, c, exact_c, fabs(c - exact_c)) << std::endl;
            }
        }
    }

    showResultMessageAll("testMulti2", bIsPassed);
}

void testMulti3()
{
    bool bIsPassed = true;
    const int dim = 100;

    TlDistributeSymmetricMatrix A(dim);
    TlSymmetricMatrix exactA(dim);
    {
        double count = 1.0;
        for (int i = 0; i < dim; ++i) {
            for (int j = 0; j < i; ++j) {
                A(i, j) = count;
                exactA(i, j) = count;
                count += 1.0 / 100;
            }
        }
    }

    TlDistributeMatrix B(dim, dim);
    TlMatrix exactB(dim, dim);
    {
        double count = 1.82;
        for (int i = 0; i < dim; ++i) {
            for (int j = 0; j < dim; ++j) {
                B(i, j) = count;
                exactB(i, j) = count;
                count -= 1.0 / 100;
            }
        }
    }

    TlDistributeMatrix C = B * A;
    TlMatrix exactC = exactB * exactA;

    // judge
    for (int r = 0; r < dim; ++r) {
        for (int s = 0; s < dim; ++s) {
            double c = C(r, s);
            double exact_c = exactC(r, s);

            if (fabs(c - exact_c) > 1.0E-8) {
                bIsPassed = false;
                std::cout << TlUtils::format("(%2d, %2d): actual= %f, expect= %f, err= %e", r, s, c, exact_c, fabs(c - exact_c)) << std::endl;
            }
        }
    }

    showResultMessageAll("testMulti3", bIsPassed);
}

void testSave()
{
    bool bIsPassed = true;

    const int nSize = 100;
    TlDistributeSymmetricMatrix A(nSize);
    TlSymmetricMatrix exactA(nSize);

    for (int r = 0; r < nSize; ++r) {
        for (int c = 0; c <= r; ++c) {
            const double v = double(rand()) / double(rand());
            A(r, c) = v;
            exactA(r, c) = v;
        }
    }

    A.save("test.symmetricmatrix_palallel");
    exactA.save("test.symmetricmatrix_exact");

    showResultMessageAll("testSave", bIsPassed);
}

void testLoad()
{
    bool bIsPassed = true;

    TlDistributeSymmetricMatrix A;
    A.load("Sab.mtx");
    TlSymmetricMatrix exactA;
    exactA.load("Sab.mtx");

    if ((A.getNumOfRows() != exactA.getNumOfRows()) || (A.getNumOfCols() != exactA.getNumOfCols())) {
        std::cout << "row = " << A.getNumOfRows() << std::endl;
        std::cout << "col = " << A.getNumOfRows() << std::endl;
        bIsPassed = false;
    }

    //A.print(std::cout);
    //exactA.print(std::cout);

    const int nSize = A.getNumOfRows();
    for (int r = 0; r < nSize; ++r) {
        for (int c = 0; c <= r; ++c) {
            if (fabs(A(r, c) - exactA(r, c)) > 1.0E-5) {
                std::cout << TlUtils::format("(%2d, %2d): actual= %f, expect= %f, err= %e",
                                             r, c, A(r,c), exactA(r,c), fabs(A(r,c) - exactA(r,c))) << std::endl;
                bIsPassed = false;
                break;
            }
        }
    }

    showResultMessageAll("testLoad", bIsPassed);
}

void testInverse()
{
    TlCommunicate& rComm = TlCommunicate::getInstance();
    bool bIsPassed = true;

    TlDistributeSymmetricMatrix A(9);
    //TlSymmetricMatrix exactA(9);
    {
        double count = 1.0;
        for (int i = 0; i < 9; ++i) {
            for (int j = 0; j <= i; ++j) {
                A(i, j) = count;
                //exactA(i, j) = count;
                count += 1.0 / 100;
            }
        }
    }

    TlDistributeSymmetricMatrix B = A;
    //TlSymmetricMatrix exactB = exactA;

    bool check = B.inverse();
    if (check == false) {
        std::cout << "inverse function returns false!" << std::endl;
    }
    //exactB.inverse();

    TlDistributeMatrix C = A * B;
//   C.print(std::cout);
//   {
//     TlDistributeMatrix X = A;
//     A.print(std::cout);
//     X.print(std::cout);
//   }

    // judge
    for (int r = 0; r < 9; ++r) {
        for (int s = 0; s < 9; ++s) {
            double c = C(r, s);
            double exact_c = (r == s) ? 1.0 : 0.0;
            //double c = B(r, s);
            //double exact_c = exactB(r, s);

            if (fabs(c - exact_c) > 1.0E-5) {
                bIsPassed = false;
                if (rComm.isMaster() == true) {
                    std::cout << TlUtils::format("(%2d, %2d): actual= %+e, expect= %+e, err= %+e", r, s, c, exact_c, fabs(c - exact_c)) << std::endl;
                }
                rComm.barrier();
            }
        }
    }

    showResultMessageAll("testInverse", bIsPassed);
}

void testInverse2()
{
    TlCommunicate& rComm = TlCommunicate::getInstance();
    bool bIsPassed = true;
    const int nSize = 100;

    TlDistributeSymmetricMatrix A(nSize);
    TlSymmetricMatrix exactA(nSize);
    {
        //double count = 1.0;
        for (int i = 0; i < nSize; ++i) {
            for (int j = 0; j <= i; ++j) {
                //const double value = double(rand())/double(RAND_MAX) -0.5;
                const double value = 1.0 + 0.1*i;
                A(i, j) = value;
                exactA(i, j) = value;
            }

            // 対角要素はノンゼロ
            //A(i, i) = 1.0;
            //exactA(i, i) = 1.0;
        }
    }

    //A.print(std::cout);
    //exactA.print(std::cout);

    TlDistributeSymmetricMatrix Ainv = A;
    bool check = Ainv.inverse();
    if (check == false) {
        std::cout << "inverse function returns false!" << std::endl;
    }

    TlSymmetricMatrix exactAinv = exactA;
    exactA.inverse();

    const TlDistributeMatrix C = A * Ainv;
    const TlMatrix exactC = exactA * exactAinv;

    //C.print(std::cout);

    // judge
    for (int r = 0; r < nSize; ++r) {
        for (int s = 0; s < nSize; ++s) {
            double c = C(r, s);
            double exact_c = exactC(r, s);

            if (fabs(c - exact_c) > 1.0E-9) {
                bIsPassed = false;
                if (rComm.isMaster() == true) {
                    std::cout
                        << TlUtils::format("(%2d, %2d): actual= %f, expect= %f, err= %e", r, s, c, exact_c, fabs(c - exact_c))
                        << std::endl;
                }
                rComm.barrier();
            }
        }
    }

    showResultMessageAll("testInverse2", bIsPassed);
}


void testDiagonal()
{
    TlCommunicate& rComm = TlCommunicate::getInstance();
    bool bIsPassed = true;

    const int size = 100;
    TlDistributeSymmetricMatrix A(size);
    TlSymmetricMatrix exactA(size);
    {
        double count = 1.0;
        for (int i = 0; i < size; ++i) {
            for (int j = 0; j <= i; ++j) {
                A(i, j) = count;
                exactA(i, j) = count;
                count += 1.0 / 100;
            }
        }
    }

    TlVector eigVal;
    TlDistributeMatrix eigVec;
    A.diagonal(&eigVal, &eigVec);

    TlVector exactEigVal;
    TlMatrix exactEigVec;
    exactA.diagonal(&exactEigVal, &exactEigVec);

    // judge
    int nSize = exactEigVal.getSize();
    for (int i = 0; i < nSize; ++i) {
        if (fabs(eigVal[i] - exactEigVal[i]) > 1.0E-9) {
            bIsPassed = false;
            if (rComm.isMaster() == true) {
                std::cout << TlUtils::format("eigval[%2d]: actual= %f, expect= %f, err= %e",
                                             i, eigVal[i], exactEigVal[i], fabs(eigVal[i] - exactEigVal[i])) << std::endl;
            }
            rComm.barrier();
        }
    }

    for (int r = 0; r < size; ++r) {
        TlVector eigvec(size);
        TlVector exact_eigvec(size);
        for (int s = 0; s < size; ++s) {
            eigvec[s] = eigVec(r, s);
            exact_eigvec[s] = exactEigVec(r, s);
        }
        TlVector eigvec2 = -1.0 * eigvec;

        for (int s = 0; s < size; ++s) {
            if ((fabs(eigvec[s] - exact_eigvec[s]) > 1.0E-9) &&
                    (fabs(eigvec2[s] - exact_eigvec[s]) > 1.0E-9)) {
                bIsPassed = false;
                if (rComm.isMaster() == true) {
                    std::cout
                        << TlUtils::format("(%2d, %2d): actual= %f, expect= %f, err= %e",
                                           r, s, eigvec[s], exact_eigvec[s], fabs(eigvec[s] - exact_eigvec[s]))
                        << std::endl;
                }
                rComm.barrier();
            }
        }
    }

    showResultMessageAll("testDiagonal", bIsPassed);
}


void testDot()
{
    // TlCommunicate& rComm = TlCommunicate::getInstance();
    bool bIsPassed = true;

    TlDistributeSymmetricMatrix A;
    TlSymmetricMatrix exactA;
    makeMatrixA(&A, &exactA);

    TlDistributeSymmetricMatrix B;
    TlSymmetricMatrix exactB;
    makeMatrixB(&B, &exactB);
   
    A.dot(B);
    exactA.dot(exactB);
    
    // judge
    const int dim = A.getNumOfRows();
    for (int r = 0; r < dim; ++r) {
        for (int s = 0; s <= r; ++s) {
            double a = A.get(r, s);
            double exact_a = exactA.get(r, s);

            if (std::fabs(a - exact_a) > 1.0E-12) {
                bIsPassed = false;
                std::cout << TlUtils::format("(%2d, %2d): actual= %f, expect= %f, err= %e", r, s,
                                             a, exact_a, fabs(a - exact_a))
                          << std::endl;
            }
        }
    }

    showResultMessageAll("testDot", bIsPassed);
}


void testSum()
{
    // TlCommunicate& rComm = TlCommunicate::getInstance();
    bool bIsPassed = true;

    TlDistributeSymmetricMatrix A;
    TlSymmetricMatrix exactA;
    makeMatrixA(&A, &exactA);
    
    double sum = A.sum();
    double exactSum = exactA.sum();
    
    // judge
    if (std::fabs(sum - exactSum) > 1.0E-12) {
        bIsPassed = false;
        std::cout << TlUtils::format("actual= %f, expect= %f, err= %e",
                                     sum, exactSum,
                                     fabs(sum - exactSum)) << std::endl;
    }
    
    showResultMessageAll("testSum", bIsPassed);
}


void testMergePartialMatrix()
{
    TlCommunicate& rComm = TlCommunicate::getInstance();
    bool bIsPassed = true;

    int dim = 10;
    
    TlSparseSymmetricMatrix sm(dim);
    TlSymmetricMatrix exactA(dim);
    {
        double count = 1.0;
        for (int i = 0; i < dim; ++i) {
            for (int j = 0; j <= i; ++j) {
                if ((j % rComm.getNumOfProc()) == rComm.getRank()) {
                    sm(i, j) = count;
                }
                exactA(i, j) = count;
                count += 1.0 / 100;
            }
        }
    }
//     if (rComm.getRank() == 0) {
//         sm.set(0, 0, 1.0);
//     }
//     exactA.add(0, 0, 1.0);

//     if (rComm.getRank() == 1) {
//         sm.add(0, 0, 2.0);
//     }
//     exactA.set(0, 0, 2.0);

//     if (rComm.getRank() == 2) {
//         sm.set(0, 0, 4.0);
//     }
//     exactA.add(0, 0, 4.0);
    
    TlDistributeSymmetricMatrix A(dim);
    A.mergeSparseMatrix(sm);

//     rComm.barrier();
    //exactA.print(std::cout);
    //rComm.barrier();
//     A.print(std::cout);
//     rComm.barrier();

    for (int r = 0; r < dim; ++r) {
        for (int s = 0; s < dim; ++s) {
            double c = A(r, s);
            double exact_c = exactA(r, s);
            if (fabs(c - exact_c) > 1.0E-9) {
                bIsPassed = false;
                if (rComm.isMaster() == true) {
                    std::cout << TlUtils::format("(%2d, %2d): actual= %f, expect= %f, err= %e",
                                                 r, s, c, exact_c, fabs(c - exact_c))
                              << std::endl;
                }
                rComm.barrier();
            }
        }
    }

    showResultMessageAll("testMergePartialMatrix", bIsPassed);
}

// void testGetPartialMatrix()
// {
//     bool bIsPassed = true;

//     TlDistributeSymmetricMatrix A(100);
//     A(3, 17) = -51.0;
//     A(3,  9) = 33.0;

//     TlSparseSymmetricMatrix B(100);
//     B(1,  0) = 1.0;
//     B(3, 17) = 2.0;
//     B(3,  9) = 3.0;
//     B(51, 17) = 4.0;

//     A.getPartialMatrix(B);

//     if (fabs(B(1,  0) - 0.0) > 1.0E-10) {
//         bIsPassed = false;
//     }
//     if (fabs(B(3, 17) - (-51.0)) > 1.0E-10) {
//         bIsPassed = false;
//     }
//     if (fabs(B(3,  9) - 33.0) > 1.0E-10) {
//         bIsPassed = false;
//     }
//     if (fabs(B(51, 17) - 0.0) > 1.0E-10) {
//         bIsPassed = false;
//     }

//     showResultMessageAll("testGetPartialMatrix", bIsPassed);
// }


void testGetSparseMatrix2()
{
    TlCommunicate& rComm = TlCommunicate::getInstance();
    bool bIsPassed = true;

    TlDistributeSymmetricMatrix A(1000);
    A(3, 17) = -51.0;
    A(3,  9) = 33.0;

    TlSparseSymmetricMatrix B(1000);
    B(1,  0) = 1.0;
    B(3, 17) = 2.0;
    B(3,  9) = 3.0;
    B(51, 17) = 4.0;

    int fin = 0;
    const int target = 1;
    const int END_LOOP = 99;
    if (rComm.getRank() != target) {
        rComm.iReceiveData(fin, target, END_LOOP);
    }
    //std::cerr << "start loop." << std::endl;
    while (true) {
        if (rComm.getRank() == target) {
            bool isComplete = A.getSparseMatrixX(&B);
            if (isComplete == true) {
                for (int i = 0; i < rComm.getNumOfProc(); ++i) {
                    if (i != target) {
                        rComm.sendData(fin, i, END_LOOP);
                    }
                }
                break;
            }
        } else {
             A.getSparseMatrixX(NULL);
            if (rComm.test(fin) == true) {
                rComm.wait(fin);
                break;
            }
        }
    }
    //std::cerr << "end loop." << std::endl;
    A.getSparseMatrixX(NULL, true);

    if (rComm.getRank() == target) {
        if (fabs(B(1,  0) - 0.0) > 1.0E-10) {
            bIsPassed = false;
        }
        if (fabs(B(3, 17) - (-51.0)) > 1.0E-10) {
            bIsPassed = false;
        }
        if (fabs(B(3,  9) - 33.0) > 1.0E-10) {
            bIsPassed = false;
        }
        if (fabs(B(51, 17) - 0.0) > 1.0E-10) {
            bIsPassed = false;
        }
        // symmetric
        if (fabs(B(0,  1) - 0.0) > 1.0E-10) {
            bIsPassed = false;
        }
        if (fabs(B(17, 3) - (-51.0)) > 1.0E-10) {
            bIsPassed = false;
        }
        if (fabs(B(9, 3) - 33.0) > 1.0E-10) {
            bIsPassed = false;
        }
        if (fabs(B(17, 51) - 0.0) > 1.0E-10) {
            bIsPassed = false;
        }
        
        if (target != 0) {
            rComm.sendData(bIsPassed, 0);
        }
    } else if (rComm.isMaster() == true) {
        if (target != 0) {
            rComm.receiveData(bIsPassed, target);
        }
    }

    showResultMessageAll("testSparseMatrix2", bIsPassed);
}


void testGetPartialMatrix2()
{
    TlCommunicate& rComm = TlCommunicate::getInstance();
    bool bIsPassed = true;

    TlDistributeSymmetricMatrix A(1000);
    A(915, 915) =  30.0;
    A(977, 977) =   7.7;
    A(903, 917) = -51.0;
    A(903, 909) =  33.0;

    TlPartialSymmetricMatrix B(1000, 900, 900, 200);
    B.set(901, 900, 1.0);
    B.set(903, 917, 2.0);
    B.set(903, 909, 3.0);
    B.set(951, 917, 4.0);

    int fin = 0;
    const int target = 1;
    const int END_LOOP = 99;
    if (rComm.getRank() != target) {
        rComm.iReceiveData(fin, target, END_LOOP);
    }
    //std::cerr << "start loop." << std::endl;
    while (true) {
        if (rComm.getRank() == target) {
            bool isComplete = A.getPartialMatrixX(&B);
            if (isComplete == true) {
                for (int i = 0; i < rComm.getNumOfProc(); ++i) {
                    if (i != target) {
                        rComm.sendData(fin, i, END_LOOP);
                    }
                }
                break;
            }
        } else {
             A.getPartialMatrixX(NULL);
            if (rComm.test(fin) == true) {
                rComm.wait(fin);
                break;
            }
        }
    }
    //std::cerr << "end loop." << std::endl;
    A.getPartialMatrixX(NULL, true);

    if (rComm.getRank() == target) {
        if (fabs(B.get(915, 915) - 30.0) > 1.0E-10) {
            bIsPassed = false;
        }
        if (fabs(B.get(977, 977) -  7.7) > 1.0E-10) {
            bIsPassed = false;
        }
        if (fabs(B.get(901, 900) - 0.0) > 1.0E-10) {
            bIsPassed = false;
        }
        if (fabs(B.get(903, 917) - (-51.0)) > 1.0E-10) {
            bIsPassed = false;
        }
        if (fabs(B.get(903, 909) - 33.0) > 1.0E-10) {
            bIsPassed = false;
        }
        if (fabs(B.get(951, 917) - 0.0) > 1.0E-10) {
            bIsPassed = false;
        }
        // symmetric
        if (fabs(B.get(900, 901) - 0.0) > 1.0E-10) {
            bIsPassed = false;
        }
        if (fabs(B.get(917, 903) - (-51.0)) > 1.0E-10) {
            bIsPassed = false;
        }
        if (fabs(B.get(909, 903) - 33.0) > 1.0E-10) {
            bIsPassed = false;
        }
        if (fabs(B.get(917, 951) - 0.0) > 1.0E-10) {
            bIsPassed = false;
        }
        
        if (target != 0) {
            rComm.sendData(bIsPassed, 0);
        }
    } else if (rComm.isMaster() == true) {
        if (target != 0) {
            rComm.receiveData(bIsPassed, target);
        }
    }

    showResultMessageAll("testPartialMatrix2", bIsPassed);
}


