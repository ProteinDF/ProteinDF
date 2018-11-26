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
#include "TlUtils.h"
#include "tl_dense_general_matrix_arrays_roworiented.h"

typedef TlMatrixObject::index_type index_type;

void showHelp() {
  std::cout << "pdf-mat-format [options] IN_MATRIX_FILE OUT_MATRIX_FILE"
            << std::endl;
  std::cout << " OPTIONS:" << std::endl;
  std::cout << "  -h:      show help" << std::endl;
}

int main(int argc, char* argv[]) {
  TlGetopt opt(argc, argv, "h");

  std::string in_path = opt[1];
  std::string out_path = opt[2];

  TlMatrix out;
  if (TlMatrix::isLoadable(out_path)) {
    std::cerr << "loading output matrix: " << out_path << std::endl;
    out.load(out_path);
  }

  index_type numOfVectors = 0;
  index_type sizeOfVector = 0;
  int numOfSubunits = 0;
  int subunitID = 0;
  if (TlDenseGeneralMatrix_arrays_RowOriented::isLoadable(
          in_path, &numOfVectors, &sizeOfVector, &numOfSubunits, &subunitID)) {
    TlDenseGeneralMatrix_arrays_RowOriented in;
    in.load(in_path);

    out.resize(numOfVectors, sizeOfVector);
    for (index_type i = 0; i < numOfVectors; ++i) {
      if (in.getSubunitID(i) == subunitID) {
        TlVector_BLAS v = in.getVector(i);
        assert(v.getSize() == sizeOfVector);

        for (index_type j = 0; j < sizeOfVector; ++j) {
          out.set(i, j, v[j]);
        }
      }
    }

    std::cerr << "save: " << out_path << std::endl;
    out.save(out_path);

  } else {
    std::cerr << "cannot load: " << in_path << std::endl;
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
