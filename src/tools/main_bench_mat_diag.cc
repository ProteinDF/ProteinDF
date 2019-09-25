#include "main_bench_mat_common.h"

// -----------------------------------------------------------------------------
// Symmetric Matrix
// -----------------------------------------------------------------------------
template <typename SymmetricMatrixType, typename GeneralMatrixType,
          typename VectorType>
void benchDiagonal(TlMatrixObject::index_type dim) {
    TlTime time;
    time.start();
    SymmetricMatrixType A(dim);
    const double createTime = time.getElapseTime();
    time.reset();
    time.start();

    setupSymmetricMatrixForDiagonal(dim, &A);
    const double setupTime = time.getElapseTime();
    time.reset();
    time.start();

    GeneralMatrixType X;
    VectorType v;
    A.eig(&v, &X);
    const double calcTime = time.getElapseTime();

    std::cout << TlUtils::format("create:   %8.3e sec", createTime)
              << std::endl;
    std::cout << TlUtils::format("setup:    %8.3e sec", setupTime) << std::endl;
    std::cout << TlUtils::format("diagonal: %8.3e sec", calcTime) << std::endl;
}

#ifdef HAVE_EIGEN
template <typename SymmetricMatrixType, typename GeneralMatrixType,
          typename VectorType>
void benchDiagonalWithEigen(TlMatrixObject::index_type dim) {
    TlTime time;
    time.start();
    SymmetricMatrixType A(dim);
    const double createTime = time.getElapseTime();
    std::cout << TlUtils::format("create:   %8.3e sec", createTime)
              << std::endl;
    time.reset();
    time.start();

    TlDenseSymmetricMatrix_Eigen tmpA(dim);
    setupSymmetricMatrixForDiagonal(dim, &tmpA);
    A = tmpA;
    const double setupTime = time.getElapseTime();
    std::cout << TlUtils::format("setup:    %8.3e sec", setupTime) << std::endl;
    time.reset();
    time.start();

    GeneralMatrixType X;
    VectorType v;
    A.eig(&v, &X);
    const double calcTime = time.getElapseTime();

    std::cout << TlUtils::format("diagonal: %8.3e sec", calcTime) << std::endl;
}
#endif  // HAVE_EIGEN

// -----------------------------------------------------------------------------
// prepare
// -----------------------------------------------------------------------------
void showHelp() {
    std::cerr << "a benchmark (matrix diagonal) of linear algebra packages"
              << std::endl;
    std::cerr << "OPTIONS:" << std::endl;
    std::cerr << "-d <device id>: switch device id" << std::endl;
    std::cerr << "-s <size>: matrix size" << std::endl;
    std::cerr << "-c: test CPU code" << std::endl;
    std::cerr << "-g: test GPU code" << std::endl;
}

// -----------------------------------------------------------------------------
// MAIN
// -----------------------------------------------------------------------------
int main(int argc, char* argv[]) {
    TlGetopt opt(argc, argv, "hcd:s:g");

    if (opt["h"] == "defined") {
        showHelp();
        return EXIT_SUCCESS;
    }

// check NDEBUG
#ifndef NDEBUG
    std::cout << "NDEBUG is not defined in this build. This may cause bad "
                 "performance."
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
    TlMatrixObject::index_type dim = 10240;
    if (!opt["s"].empty()) {
        dim = std::atoi(opt["s"].c_str());
    }

    bool onCPU = (opt["c"] == "defined");
    bool onGPU = (opt["g"] == "defined");
    if (!(onCPU ^ onGPU)) {
        onCPU = true;
        onGPU = true;
    }

    initRand();

// ViennaCL
#ifdef HAVE_VIENNACL
    {
        TlViennaCL vcl;
        vcl.setupAllAvailableDevices();
        std::cout << vcl.listDevices() << std::endl;

        vcl.switchDevice(deviceId);
        std::cout << vcl.listCurrentDevice() << std::endl;
    }
#endif  // HAVE_VIENNACL

    // diagonal
    if (onCPU) {
        std::cout << ">>>> BLAS@SymmetricMatrix (diagonal)" << std::endl;
        benchDiagonal<TlDenseSymmetricMatrix_Lapack,
                      TlDenseGeneralMatrix_Lapack, TlDenseVector_Lapack>(dim);

#ifdef HAVE_EIGEN
        std::cout << ">>>> Eigen@SymmetricMatrix (diagonal)" << std::endl;
        benchDiagonal<TlDenseSymmetricMatrix_Eigen, TlDenseGeneralMatrix_Eigen,
                      TlDenseVector_Eigen>(dim);
#endif  // HAVE_EIGEN
    }

#ifdef HAVE_VIENNACL
    if (onGPU) {
        std::cout << ">>>> ViennaCL@SymmetricMatrix (diagonal)" << std::endl;
        benchDiagonalWithEigen<TlDenseSymmetricMatrix_ViennaCL,
                               TlDenseGeneralMatrix_ViennaCL,
                               TlDenseVector_ViennaCL>(dim);
    }
#endif  // HAVE_VIENNACL

    std::cout << "done." << std::endl;
    return EXIT_SUCCESS;
}
