#include <cstdlib>
#include <iostream>

#include "TlGetopt.h"
#include "tl_matrix_utils.h"
#include "tl_dense_general_matrix_lapack.h"
#include "tl_dense_symmetric_matrix_lapack.h"

int main(int argc, char* argv[]) {
  TlGetopt opt(argc, argv, "dhv");

  const bool bVerbose = (opt["v"] == "defined");
  const bool bDiagonalMode = (opt["d"] == "defined");

  std::string sPath = opt[1];
  if (bVerbose) {
    std::cerr << "loading... " << sPath << std::endl;
  }

  double sum = 0.0;
  if (TlMatrixUtils::isLoadable(sPath, TlMatrixObject::RLHD) == true) {
    TlDenseSymmetricMatrix_Lapack M;
    M.load(sPath);

    TlMatrixObject::index_type dim = M.getNumOfRows();

    if (bDiagonalMode) {
      for (TlMatrixObject::index_type r = 0; r < dim; ++r) {
        sum += M.get(r, r);
      }
    } else {
      for (TlMatrixObject::index_type r = 0; r < dim; ++r) {
        for (TlMatrixObject::index_type c = 0; c < r; ++c) {
          sum += 2.0 * M.get(r, c);
        }

        sum += M.get(r, r);
      }
    }
  } else if (TlMatrixUtils::isLoadable(sPath, TlMatrixObject::CSFD) == true) {
    TlDenseGeneralMatrix_Lapack M;
    M.load(sPath);

    TlMatrixObject::index_type numOfRows = M.getNumOfRows();
    TlMatrixObject::index_type numOfCols = M.getNumOfCols();
    for (TlMatrixObject::index_type r = 0; r < numOfRows; ++r) {
      for (TlMatrixObject::index_type c = 0; c < numOfCols; ++c) {
        sum += M.get(r, c);
      }
    }

  } else {
    std::cerr << "unknown file type: " << sPath << std::endl;
    return EXIT_FAILURE;
  }

  std::cout << "sum = " << sum << std::endl;

  return EXIT_SUCCESS;
}
