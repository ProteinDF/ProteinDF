#include <cstdlib>
#include <iostream>
#include <string>

#include "TlGetopt.h"
#include "TlHdf5Utils.h"
#include "TlMatrix.h"
#include "TlSymmetricMatrix.h"

void showHelp(const std::string& progname) {
  std::cout << TlUtils::format("%s [options] MATRIX_FILE HDF5_FILE",
                               progname.c_str())
            << std::endl;
  std::cout << " OPTIONS:" << std::endl;
  std::cout << "  -h:      show help" << std::endl;
}

int main(int argc, char* argv[]) {
  TlGetopt opt(argc, argv, "h");

  if (opt["h"] == "defined") {
    showHelp(opt[0]);
    return EXIT_SUCCESS;
  }

  std::string mat_path = opt[1];
  std::string hdf5_path = opt[2];

  TlMatrix mat;
  mat.load(mat_path);

  TlHdf5Utils h5(hdf5_path);

  return EXIT_SUCCESS;
}
