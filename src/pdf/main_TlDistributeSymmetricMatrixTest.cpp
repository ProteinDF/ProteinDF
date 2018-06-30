// Copyright (C) 2002-2014 The ProteinDF project
// see also AUTHORS and README.
//
// This file is part of ProteinDF.
//
// ProteinDF is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// ProteinDF is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with ProteinDF.  If not, see <http://www.gnu.org/licenses/>.

#include <cstdlib>
#include <iostream>
#include <string>

#include "TlCommunicate.h"
#include "TlDenseSymmetricMatrix_BLAS_Old.h"
#include "TlDistributeSymmetricMatrix.h"
#include "TlMatrix.h"

void showResultMessage(const std::string& sFunction, bool bIsPassed);
void showResultMessageAll(const std::string& sFunction, bool bIsPassed);

void testConstructer();
void testCopyConstructer();
void testCopyConstructer2();
void testSet();
void testOperatorPlusEqual();
void testOperatorMultiEqual();
void testMulti1();
void testMulti2();
void testSave();
void testLoad();
void testInverse();
void testInverse2();
void testDiagonal();

int main(int argc, char* argv[]) {
  // initialize
  TlCommunicate& rComm = TlCommunicate::getInstance(argc, argv);

  // ===================================================================
  //   testConstructer();
  //   testCopyConstructer();
  //   testCopyConstructer2();
  //   testSet();
  //   testOperatorPlusEqual();
  //   testOperatorMultiEqual();
  //   testMulti1();
  //   testMulti2();
  testSave();
  testLoad();
  //   testInverse();
  //   testInverse2();
  //   testDiagonal();
  // ===================================================================

  // finalize
  rComm.finalize();
  return EXIT_SUCCESS;
}

void showResultMessage(const std::string& sFunction, bool bIsPassed) {
  TlCommunicate& rComm = TlCommunicate::getInstance();
  if (rComm.isMaster() == true) {
    std::cout << "TEST: " << sFunction << "() ";
    std::cout << ((bIsPassed == true) ? "." : "F") << std::endl;
  }
  rComm.barrier();
}

void showResultMessageAll(const std::string& sFunction, bool bIsPassed) {
  TlCommunicate& rComm = TlCommunicate::getInstance();
  int nProc = rComm.getNumOfProc();
  int nRank = rComm.getRank();

  rComm.barrier();
  if (rComm.isMaster() == true) {
    std::cout << "TEST: " << sFunction << "() ";
  }
  for (int i = 0; i < nProc; ++i) {
    if (i == nRank) {
      std::cout << ((bIsPassed == true) ? "." : "F");
      std::cout.flush();
    }
    rComm.barrier();
  }
  if (rComm.isMaster() == true) {
    std::cout << std::endl;
  }

  rComm.barrier();
}

void testConstructer() {
  TlDistributeSymmetricMatrix matrix(100);

  int nRow = matrix.getNumOfRows();
  int nCol = matrix.getNumOfCols();

  bool bIsPassed = ((nRow == 100) && (nCol == 100)) ? true : false;

  showResultMessageAll("testConstructer", bIsPassed);
}

void testCopyConstructer() {
  TlDistributeSymmetricMatrix A(100);
  A(10, 10) = 10.0;
  A(17, 4) = 56.0;
  A(21, 18) = -2.5;
  A(0, 10) = 12.0;

  const TlDistributeSymmetricMatrix B(A);

  // A.print(std::cout);
  // B.print(std::cout);

  bool bIsPassed = true;
  if ((B.getNumOfRows() != 100) || (B.getNumOfCols() != 100)) {
    bIsPassed = false;
  }
  if (std::fabs(B(10, 10) - 10.0) > 1.0E-16) {
    bIsPassed = false;
    std::cout << "B(10, 10) = " << B(10, 10) << ", exact = 10.0" << std::endl;
  }
  if (std::fabs(B(17, 4) - 56.0) > 1.0E-16) {
    bIsPassed = false;
    std::cout << "B(17,  4) = " << B(17, 4) << ", exact = 56.0" << std::endl;
  }
  if (std::fabs(B(4, 17) - 56.0) > 1.0E-16) {
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
  if (std::fabs(B(0, 10) - 12.0) > 1.0E-16) {
    bIsPassed = false;
    std::cout << "B(0, 10) = " << B(0, 10) << ", exact = 12.0" << std::endl;
  }
  if (std::fabs(B(10, 0) - 12.0) > 1.0E-16) {
    bIsPassed = false;
    std::cout << "B(10, 0) = " << B(10, 0) << ", exact = 12.0" << std::endl;
  }

  showResultMessageAll("testCopyConstructer", bIsPassed);
}

void testCopyConstructer2() {
  TlDistributeSymmetricMatrix A(100);
  A(10, 10) = 10.0;
  A(17, 4) = 56.0;
  A(21, 18) = -2.5;
  A(0, 10) = 12.0;

  // A.print(std::cout);

  const TlDistributeMatrix B(A);

  bool bIsPassed = true;
  if ((B.getNumOfRows() != 100) || (B.getNumOfCols() != 100)) {
    bIsPassed = false;
  }
  if (std::fabs(B(10, 10) - 10.0) > 1.0E-16) {
    bIsPassed = false;
  }
  if (std::fabs(B(17, 4) - 56.0) > 1.0E-16) {
    bIsPassed = false;
  }
  if (std::fabs(B(4, 17) - 56.0) > 1.0E-16) {
    bIsPassed = false;
  }
  if (std::fabs(B(21, 18) - (-2.5)) > 1.0E-16) {
    bIsPassed = false;
  }
  if (std::fabs(B(18, 21) - (-2.5)) > 1.0E-16) {
    bIsPassed = false;
  }
  if (std::fabs(B(0, 10) - 12.0) > 1.0E-16) {
    bIsPassed = false;
  }
  if (std::fabs(B(10, 0) - 12.0) > 1.0E-16) {
    bIsPassed = false;
  }

  showResultMessageAll("testCopyConstructer2", bIsPassed);
}

void testSet() {
  TlDistributeSymmetricMatrix A(25, 4);
  A(10, 10) = 10.0;
  A(17, 4) = 56.0;
  A(21, 18) = -2.5;
  A(0, 10) = 12.0;

  // A.print(std::cout);

  bool bIsPassed = true;
  if ((A.getNumOfRows() != 25) || (A.getNumOfCols() != 25)) {
    bIsPassed = false;
  }
  if (std::fabs(A(10, 10) - 10.0) > 1.0E-16) {
    bIsPassed = false;
  }
  if (std::fabs(A(17, 4) - 56.0) > 1.0E-16) {
    bIsPassed = false;
  }
  if (std::fabs(A(4, 17) - 56.0) > 1.0E-16) {
    bIsPassed = false;
  }
  if (std::fabs(A(21, 18) - (-2.5)) > 1.0E-16) {
    bIsPassed = false;
  }
  if (std::fabs(A(18, 21) - (-2.5)) > 1.0E-16) {
    bIsPassed = false;
  }
  if (std::fabs(A(0, 10) - 12.0) > 1.0E-16) {
    bIsPassed = false;
  }
  if (std::fabs(A(10, 0) - 12.0) > 1.0E-16) {
    bIsPassed = false;
  }

  showResultMessageAll("testSet", bIsPassed);
}

void testOperatorPlusEqual() {
  TlDistributeSymmetricMatrix A(25, 4);
  A(10, 10) = 10.0;
  A(17, 4) = 56.0;
  A(21, 18) = -2.5;
  A(0, 10) = 12.0;

  TlDistributeSymmetricMatrix B(25, 4);
  B(10, 10) = 5.0;
  B(0, 10) = -12.0;
  B(3, 17) = 51.0;

  A += B;

  bool bIsPassed = true;
  if ((A.getNumOfRows() != 25) || (A.getNumOfCols() != 25)) {
    bIsPassed = false;
  }
  if (std::fabs(A(10, 10) - 15.0) > 1.0E-16) {
    bIsPassed = false;
  }
  if (std::fabs(A(17, 4) - 56.0) > 1.0E-16) {
    bIsPassed = false;
  }
  if (std::fabs(A(4, 17) - 56.0) > 1.0E-16) {
    bIsPassed = false;
  }
  if (std::fabs(A(21, 18) - (-2.5)) > 1.0E-16) {
    bIsPassed = false;
  }
  if (std::fabs(A(18, 21) - (-2.5)) > 1.0E-16) {
    bIsPassed = false;
  }
  if (std::fabs(A(0, 10) - 0.0) > 1.0E-16) {
    bIsPassed = false;
  }
  if (std::fabs(A(10, 0) - 0.0) > 1.0E-16) {
    bIsPassed = false;
  }
  if (std::fabs(A(3, 17) - 51.0) > 1.0E-16) {
    bIsPassed = false;
  }
  if (std::fabs(A(17, 3) - 51.0) > 1.0E-16) {
    bIsPassed = false;
  }

  showResultMessageAll("testOperatorPlusEqual", bIsPassed);
}

void testOperatorMultiEqual() {
  TlDistributeSymmetricMatrix A(25, 4);
  A(10, 10) = 10.0;
  A(17, 4) = 56.0;
  A(21, 18) = -2.5;
  A(0, 10) = 12.0;

  A *= 3.0;

  bool bIsPassed = true;
  if ((A.getNumOfRows() != 25) || (A.getNumOfCols() != 25)) {
    bIsPassed = false;
  }
  if (std::fabs(A(10, 10) - 10.0 * 3.0) > 1.0E-16) {
    bIsPassed = false;
  }
  if (std::fabs(A(17, 4) - 56.0 * 3.0) > 1.0E-16) {
    bIsPassed = false;
  }
  if (std::fabs(A(4, 17) - 56.0 * 3.0) > 1.0E-16) {
    bIsPassed = false;
  }
  if (std::fabs(A(21, 18) - (-2.5) * 3.0) > 1.0E-16) {
    bIsPassed = false;
  }
  if (std::fabs(A(18, 21) - (-2.5) * 3.0) > 1.0E-16) {
    bIsPassed = false;
  }
  if (std::fabs(A(0, 10) - 12.0 * 3.0) > 1.0E-16) {
    bIsPassed = false;
  }
  if (std::fabs(A(10, 0) - 12.0 * 3.0) > 1.0E-16) {
    bIsPassed = false;
  }

  showResultMessageAll("testOperatorMultiEqual", bIsPassed);
}

void testMulti1() {
  bool bIsPassed = true;

  TlDistributeSymmetricMatrix A(9, 4);
  TlDenseSymmetricMatrix_BLAS_Old exactA(9);
  {
    double count = 1.0;
    for (int i = 0; i < 9; ++i) {
      for (int j = 0; j < i; ++j) {
        A(i, j) = count;
        exactA(i, j) = count;
        count += 1.0 / 100;
      }
    }
  }

  TlDistributeMatrix B(9, 9, 4);
  TlMatrix exactB(9, 9);
  {
    double count = 1.82;
    for (int i = 0; i < 9; ++i) {
      for (int j = 0; j < 9; ++j) {
        B(i, j) = count;
        exactB(i, j) = count;
        count -= 1.0 / 100;
      }
    }
  }

  TlDistributeMatrix C = A * B;
  TlMatrix exactC = exactA * exactB;

  // judge
  for (int r = 0; r < 9; ++r) {
    for (int s = 0; s < 9; ++s) {
      double c = C(r, s);
      double exact_c = exactC(r, s);

      if (fabs(c - exact_c) > 1.0E-12) {
        bIsPassed = false;
        std::cout << TlUtils::format(
                         "(%2d, %2d): actual= %f, expect= %f, err= %e", r, s, c,
                         exact_c, fabs(c - exact_c))
                  << std::endl;
      }
    }
  }

  showResultMessageAll("testMulti1", bIsPassed);
}

void testMulti2() {
  bool bIsPassed = true;

  TlDistributeSymmetricMatrix A(9, 4);
  TlDenseSymmetricMatrix_BLAS_Old exactA(9);
  {
    double count = 1.0;
    for (int i = 0; i < 9; ++i) {
      for (int j = 0; j < i; ++j) {
        A(i, j) = count;
        exactA(i, j) = count;
        count += 1.0 / 100;
      }
    }
  }

  TlDistributeMatrix B(9, 9, 4);
  TlMatrix exactB(9, 9);
  {
    double count = 1.82;
    for (int i = 0; i < 9; ++i) {
      for (int j = 0; j < 9; ++j) {
        B(i, j) = count;
        exactB(i, j) = count;
        count -= 1.0 / 100;
      }
    }
  }

  TlDistributeMatrix C = B * A;
  TlMatrix exactC = exactB * exactA;

  // judge
  for (int r = 0; r < 9; ++r) {
    for (int s = 0; s < 9; ++s) {
      double c = C(r, s);
      double exact_c = exactC(r, s);

      if (fabs(c - exact_c) > 1.0E-12) {
        bIsPassed = false;
        std::cout << TlUtils::format(
                         "(%2d, %2d): actual= %f, expect= %f, err= %e", r, s, c,
                         exact_c, fabs(c - exact_c))
                  << std::endl;
      }
    }
  }

  showResultMessageAll("testMulti2", bIsPassed);
}

void testSave() {
  bool bIsPassed = true;

  const int nSize = 100;
  TlDistributeSymmetricMatrix A(nSize);
  TlDenseSymmetricMatrix_BLAS_Old exactA(nSize);

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

void testLoad() {
  bool bIsPassed = true;

  TlDistributeSymmetricMatrix A;
  A.load("test.symmetricmatrix_exact");
  TlDenseSymmetricMatrix_BLAS_Old exactA;
  exactA.load("test.symmetricmatrix_exact");

  if ((A.getNumOfRows() != exactA.getNumOfRows()) ||
      (A.getNumOfCols() != exactA.getNumOfCols())) {
    std::cout << "row = " << A.getNumOfRows() << std::endl;
    std::cout << "col = " << A.getNumOfRows() << std::endl;
    bIsPassed = false;
  }

  const int nSize = A.getNumOfRows();
  for (int r = 0; r < nSize; ++r) {
    for (int c = 0; c <= r; ++c) {
      if (fabs(A(r, c) - exactA(r, c)) < 1.0E-5) {
        bIsPassed = false;
        break;
      }
    }
  }

  showResultMessageAll("testLoad", bIsPassed);
}

void testInverse() {
  TlCommunicate& rComm = TlCommunicate::getInstance();
  bool bIsPassed = true;

  TlDistributeSymmetricMatrix A(9, 4);
  {
    double count = 1.0;
    for (int i = 0; i < 9; ++i) {
      for (int j = 0; j <= i; ++j) {
        A(i, j) = count;
        count += 1.0 / 100;
      }
    }
  }

  TlDistributeSymmetricMatrix B = A;

  bool check = B.inverse();
  if (check == false) {
    std::cout << "inverse function returns false!" << std::endl;
  }

  TlDistributeMatrix C = A * B;

  // judge
  for (int r = 0; r < 9; ++r) {
    for (int s = 0; s < 9; ++s) {
      double c = C(r, s);
      double exact_c = (r == s) ? 1.0 : 0.0;

      if (fabs(c - exact_c) > 1.0E-9) {
        bIsPassed = false;
        if (rComm.isMaster() == true) {
          std::cout << TlUtils::format(
                           "(%2d, %2d): actual= %f, expect= %f, err= %e", r, s,
                           c, exact_c, fabs(c - exact_c))
                    << std::endl;
        }
        rComm.barrier();
      }
    }
  }

  showResultMessageAll("testInverse", bIsPassed);
}

void testInverse2() {
  TlCommunicate& rComm = TlCommunicate::getInstance();
  bool bIsPassed = true;
  const int nSize = 100;

  TlDistributeSymmetricMatrix A(nSize);
  TlDenseSymmetricMatrix_BLAS_Old exactA(nSize);
  {
    // double count = 1.0;
    for (int i = 0; i < nSize; ++i) {
      for (int j = 0; j < i; ++j) {
        const double value = double(rand()) / double(RAND_MAX) - 0.5;
        A(i, j) = value;
        exactA(i, j) = value;
      }

      // 対角要素はノンゼロ
      A(i, i) = 1.0;
      exactA(i, i) = 1.0;
    }
  }

  // A.print(std::cout);
  // exactA.print(std::cout);

  bool check = A.inverse();
  if (check == false) {
    std::cout << "inverse function returns false!" << std::endl;
  }

  exactA.inverse();

  // judge
  for (int r = 0; r < 9; ++r) {
    for (int s = 0; s < 9; ++s) {
      double c = A(r, s);
      double exact_c = exactA(r, s);

      if (fabs(c - exact_c) > 1.0E-9) {
        bIsPassed = false;
        if (rComm.isMaster() == true) {
          std::cout << TlUtils::format(
                           "(%2d, %2d): actual= %f, expect= %f, err= %e", r, s,
                           c, exact_c, fabs(c - exact_c))
                    << std::endl;
        }
        rComm.barrier();
      }
    }
  }

  showResultMessageAll("testInverse2", bIsPassed);
}

void testDiagonal() {
  TlCommunicate& rComm = TlCommunicate::getInstance();
  bool bIsPassed = true;

  TlDistributeSymmetricMatrix A(9, 4);
  TlDenseSymmetricMatrix_BLAS_Old exactA(9);
  {
    double count = 1.0;
    for (int i = 0; i < 9; ++i) {
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
        std::cout << TlUtils::format(
                         "eigval[%2d]: actual= %f, expect= %f, err= %e", i,
                         eigVal[i], exactEigVal[i],
                         fabs(eigVal[i] - exactEigVal[i]))
                  << std::endl;
      }
      rComm.barrier();
    }
  }

  for (int r = 0; r < 9; ++r) {
    for (int s = 0; s < 9; ++s) {
      double c = eigVec(r, s);
      double exact_c = exactEigVec(r, s);

      if (fabs(c - exact_c) > 1.0E-9) {
        bIsPassed = false;
        if (rComm.isMaster() == true) {
          std::cout << TlUtils::format(
                           "(%2d, %2d): actual= %f, expect= %f, err= %e", r, s,
                           c, exact_c, fabs(c - exact_c))
                    << std::endl;
        }
        rComm.barrier();
      }
    }
  }

  showResultMessageAll("testDiagonal", bIsPassed);
}
