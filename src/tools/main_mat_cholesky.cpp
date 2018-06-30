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

#include "TlGetopt.h"
#include "TlMatrix.h"
#include "tl_dense_symmetric_matrix_blas_old.h"

void showHelp(const std::string& progname) {
  std::cout << "cholesky [options] input_file_path" << std::endl;
  std::cout << " OPTIONS:" << std::endl;
  // std::cout << "  -l FILE: save vector for eigen values" << std::endl;
  // std::cout << "  -x FILE: save matrix for eigen vector" << std::endl;
  std::cout << "  -h:      show help" << std::endl;
  std::cout << "  -v:      verbose" << std::endl;
}

int main(int argc, char* argv[]) {
  TlGetopt opt(argc, argv, "hv");

  if (opt["h"] == "defined") {
    showHelp(opt[0]);
    return EXIT_SUCCESS;
  }

  const bool bVerbose = (opt["v"] == "defined");

  if (opt.getCount() <= 1) {
    showHelp(opt[0]);
    return EXIT_FAILURE;
  }
  std::string inputMatrixPath = opt[1];

  if (bVerbose == true) {
    std::cerr << "load matrix: " << inputMatrixPath << std::endl;
  }
  if (TlDenseSymmetricMatrix_BLAS_Old::isLoadable(inputMatrixPath) != true) {
    std::cerr << "can not open file: " << inputMatrixPath << std::endl;
    return EXIT_FAILURE;
  }

  TlDenseSymmetricMatrix_BLAS_Old A;
  A.load(inputMatrixPath);
  const int numOfDims = A.getNumOfRows();

  if (bVerbose == true) {
    std::cerr << "running..." << inputMatrixPath << std::endl;
  }
  TlMatrix L = A.choleskyFactorization2(1.0E-16);

  TlMatrix Lt = L;
  Lt.transpose();

  TlMatrix LL = L * Lt;

  std::cout << ">>>> L" << std::endl;
  L.print(std::cout);
  std::cout << ">>>> A" << std::endl;
  A.print(std::cout);
  std::cout << ">>>> LL" << std::endl;
  LL.print(std::cout);

  return EXIT_SUCCESS;
}
