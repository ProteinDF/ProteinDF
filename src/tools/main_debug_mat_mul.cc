#include <iostream>

#include "TlGetopt.h"
#include "TlTime.h"
#include "TlUtils.h"
#include "tool_common.h"

// -----------------------------------------------------------------------------
// prepare
// -----------------------------------------------------------------------------
void showHelp() {
  std::cerr << "a matrix(multiplication) debug of linear algebra packages"
            << std::endl;
  std::cerr << "OPTIONS:" << std::endl;
  std::cerr << "-d <device id>: switch device id" << std::endl;
  std::cerr << "-s <size>: matrix size" << std::endl;
}

// -----------------------------------------------------------------------------
// MAIN
// -----------------------------------------------------------------------------
int main(int argc, char* argv[]) {
  TlGetopt opt(argc, argv, "hd:s:");

  if (opt["h"] == "defined") {
    showHelp();
    return EXIT_SUCCESS;
  }

// check NDEBUG
#ifndef NDEBUG
  std::cout
      << "NDEBUG is not defined in this build. This may cause bad performance."
      << std::endl;
#endif  // NDEBUG

// check openmp
#ifdef _OPENMP
  std::cout << "# of Procs: " << omp_get_num_procs() << std::endl;
  std::cout << "# of OpenMP Max threads: " << omp_get_max_threads()
            << std::endl;
#endif  // _OPENMP

  int deviceId = 0;
  if (!opt["d"].empty()) {
    deviceId = std::atoi(opt["d"].c_str());
  }
  TlMatrixObject::index_type dim = 10;
  if (!opt["s"].empty()) {
    dim = std::atoi(opt["s"].c_str());
  }

  std::cout << "dim: " << dim << std::endl;

  initRand();

  // --------------------------------------------------------------------------
  // OpenCL device
  // --------------------------------------------------------------------------
#ifdef HAVE_VIENNACL
  {
    TlViennaCL vcl;
    vcl.setupAllAvailableDevices();
    std::cout << vcl.listDevices() << std::endl;

    vcl.switchDevice(deviceId);
    std::cout << vcl.listCurrentDevice() << std::endl;
  }
#endif  // HAVE_VIENNACL

  const TlMatrixObject::index_type row = dim;
  const TlMatrixObject::index_type col = dim;

  // --------------------------------------------------------------------------
  // general matrix
  {
    TlDenseGeneralMatrix_Eigen eigen_A =
        getGeneralMatrixA<TlDenseGeneralMatrix_Eigen>(dim, dim);
    TlDenseGeneralMatrix_Eigen eigen_B =
        getGeneralMatrixB<TlDenseGeneralMatrix_Eigen>(dim, dim);
    TlDenseGeneralMatrix_Eigen eigen_C = eigen_A * eigen_B;
    std::cout << eigen_C << std::endl;

    TlDenseGeneralMatrix_ViennaCL vcl_A =
        getGeneralMatrixA<TlDenseGeneralMatrix_ViennaCL>(dim, dim);
    TlDenseGeneralMatrix_ViennaCL vcl_B =
        getGeneralMatrixB<TlDenseGeneralMatrix_ViennaCL>(dim, dim);
    TlDenseGeneralMatrix_ViennaCL vcl_C = eigen_A * eigen_B;
    std::cout << vcl_C << std::endl;

    std::cout << ">>>> diff infomation >>>>" << std::endl;
    for (TlMatrixObject::index_type r = 0; r < row; ++r) {
      for (TlMatrixObject::index_type c = 0; c < col; ++c) {
        const double eigen_v = eigen_C.get(r, c);
        const double vcl_v = vcl_C.get(r, c);
        const double delta = eigen_v - vcl_v;
        if (std::fabs(delta) > 1.0E-5) {
          std::cout << TlUtils::format("! (%5d, %5d): % 8.3e != %8.3e (% 8.3e)",
                                       r, c, eigen_v, vcl_v, delta)
                    << std::endl;
        }
      }
    }
  }

  std::cout << "done." << std::endl;
  return EXIT_SUCCESS;
}
