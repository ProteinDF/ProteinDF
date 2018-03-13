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
#include "tl_dense_general_matrix_blas_old.h"

void showHelp() {
  std::cout << "pdf-mat-resize [options] input_path output_path" << std::endl;
  std::cout << " OPTIONS:" << std::endl;
  std::cout << "  -r rows: new number of rows" << std::endl;
  std::cout << "  -c cols: new number of cols" << std::endl;
  std::cout << "  -h:      show help" << std::endl;
  std::cout << "  -v:      verbose" << std::endl;
}

int main(int argc, char* argv[]) {
  TlGetopt opt(argc, argv, "r:c:hv");

  if (opt["h"] == "defined") {
    showHelp();
    return EXIT_SUCCESS;
  }

  const bool bVerbose = (opt["v"] == "defined");
  int newNumOfRows = 0;
  if (!opt["r"].empty()) {
    newNumOfRows = std::atoi(opt["r"].c_str());
  }
  int newNumOfCols = 0;
  if (!opt["c"].empty()) {
    newNumOfCols = std::atoi(opt["c"].c_str());
  }

  if (opt.getCount() <= 2) {
    showHelp();
    return EXIT_FAILURE;
  }
  std::string inputMatrixPath = opt[1];
  std::string outputMatrixPath = opt[2];

  if (bVerbose == true) {
    std::cerr << "load matrix: " << inputMatrixPath << std::endl;
  }

  TlDenseGeneralMatrix_BLAS_old A;
  A.load(inputMatrixPath);

  TlDenseGeneralMatrix_BLAS_old::index_type numOfRows =
      (newNumOfRows != 0) ? newNumOfRows : A.getNumOfRows();
  TlDenseGeneralMatrix_BLAS_old::index_type numOfCols =
      (newNumOfCols != 0) ? newNumOfCols : A.getNumOfCols();
  A.resize(numOfRows, numOfCols);

  if (bVerbose == true) {
    std::cerr << "save matrix: " << outputMatrixPath << std::endl;
  }
  A.save(outputMatrixPath);

  return EXIT_SUCCESS;
}
