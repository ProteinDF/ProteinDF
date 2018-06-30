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
#include "tl_dense_symmetric_matrix_blas_old.h"
#include "tl_matrix_utils.h"

void showHelp(const std::string& name) {
  std::cout << TlUtils::format("%s [options] INPUT OUTPUT", name.c_str())
            << std::endl;
  std::cout << "convert normal matrix to symmetric matrix" << std::endl;
  std::cout << " OPTIONS:" << std::endl;
  std::cout << "  -h:      show help" << std::endl;
}

int main(int argc, char* argv[]) {
  TlGetopt opt(argc, argv, "h");

  if ((opt["h"] == "defined") || (opt.getCount() < 3)) {
    showHelp(opt[0]);
    return EXIT_FAILURE;
  }

  std::string input_path = opt[1];
  std::string output_path = opt[2];

  if (TlMatrixUtils::isLoadable(input_path, TlMatrixObject::CSFD) != true) {
    std::cerr << TlUtils::format("type mismatch: %s", input_path.c_str())
              << std::endl;
    return EXIT_FAILURE;
  }
  TlDenseGeneralMatrix_BLAS_old in;
  in.load(input_path);

  if (in.getNumOfRows() != in.getNumOfCols()) {
    std::cerr << TlUtils::format("wrong matrix size: %d x %d",
                                 in.getNumOfRows(), in.getNumOfCols())
              << std::endl;
    return EXIT_FAILURE;
  }

  TlDenseSymmetricMatrix_BLAS_Old out(in);
  out.save(output_path);

  return EXIT_SUCCESS;
}
