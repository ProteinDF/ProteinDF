#include <cstdlib>
#include <iostream>

#include "TlGetopt.h"
#include "TlMatrix.h"
#include "TlSymmetricMatrix.h"

int main(int argc, char* argv[]) {
  TlGetopt opt(argc, argv, "dhv");

  const bool bVerbose = (opt["v"] == "defined");
  const bool bDiagonalMode = (opt["d"] == "defined");

  std::string sPath = opt[1];
  if (bVerbose) {
    std::cerr << "loading... " << sPath << std::endl;
  }

  double sum = 0.0;
  if (TlSymmetricMatrix::isLoadable(sPath) == true) {
    TlSymmetricMatrix M;
    M.load(sPath);

    TlMatrix::index_type dim = M.getNumOfRows();

    if (bDiagonalMode) {
      for (TlMatrix::index_type r = 0; r < dim; ++r) {
        sum += M.get(r, r);
      }
    } else {
      for (TlMatrix::index_type r = 0; r < dim; ++r) {
        for (TlMatrix::index_type c = 0; c < r; ++c) {
          sum += 2.0 * M.get(r, c);
        }

        sum += M.get(r, r);
      }
    }
  } else if (TlMatrix::isLoadable(sPath) == true) {
    TlMatrix M;
    M.load(sPath);

    TlMatrix::index_type numOfRows = M.getNumOfRows();
    TlMatrix::index_type numOfCols = M.getNumOfCols();
    for (TlMatrix::index_type r = 0; r < numOfRows; ++r) {
      for (TlMatrix::index_type c = 0; c < numOfCols; ++c) {
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
