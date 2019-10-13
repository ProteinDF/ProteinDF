#include "main_bench_mat_common.h"

// -----------------------------------------------------------------------------
// General Matrix
// -----------------------------------------------------------------------------
template <typename DenseMatrixType>
void benchGeneralMatrix(TlMatrixObject::index_type dim) {
    TlTime time;
    time.start();
    DenseMatrixType A(dim, dim);
    DenseMatrixType B(dim, dim);
    const double createTime = time.getElapseTime();
    std::cout << TlUtils::format("create: %8.3e sec", createTime) << std::endl;
    time.reset();
    time.start();

    setupGeneralMatrix(dim, &A);
    setupGeneralMatrix(dim, &B);
    const double setupTime = time.getElapseTime();
    std::cout << TlUtils::format("setup : %8.3e sec", setupTime) << std::endl;
    time.reset();
    time.start();

    DenseMatrixType C = A * B;
    const double mulTime = time.getElapseTime();

    // FLOP
    double FLOP = 0.0;
    {
        const double row1 = A.getNumOfRows();
        const double col1 = A.getNumOfCols();
        // const double row2 = B.getNumOfRows(); // col1 == row2
        const double col2 = B.getNumOfCols();
        FLOP = (col1 + (col1 - 1)) * row1 * col2;
        // std::cout << "FLOP: " << FLOP << std::endl;
    }

    std::cout << TlUtils::format("mul   : %8.3e sec (%e GFLOPS)", mulTime,
                                 FLOP / (mulTime * 1024.0 * 1024.0 * 1024.0))
              << std::endl;
}

#if defined(HAVE_EIGEN)
template <typename DenseMatrixType>
void benchGeneralMatrixWithEigen(TlMatrixObject::index_type dim) {
    TlTime time;
    time.start();
    DenseMatrixType A(dim, dim);
    DenseMatrixType B(dim, dim);
    const double createTime = time.getElapseTime();
    std::cout << TlUtils::format("create: %8.3e sec", createTime) << std::endl;
    time.reset();
    time.start();

    setupGeneralMatrixWithEigen(dim, &A);
    setupGeneralMatrixWithEigen(dim, &B);
    const double setupTime = time.getElapseTime();
    std::cout << TlUtils::format("setup : %8.3e sec", setupTime) << std::endl;
    time.reset();
    time.start();

    DenseMatrixType C = A * B;
    const double mulTime = time.getElapseTime();

    // FLOP
    double FLOP = 0.0;
    {
        const double row1 = A.getNumOfRows();
        const double col1 = A.getNumOfCols();
        // const double row2 = B.getNumOfRows(); // col1 == row2
        const double col2 = B.getNumOfCols();
        FLOP = (col1 + (col1 - 1)) * row1 * col2;
        // std::cout << "FLOP: " << FLOP << std::endl;
    }

    std::cout << TlUtils::format("mul   : %8.3e sec (%e GFLOPS)", mulTime,
                                 FLOP / (mulTime * 1024.0 * 1024.0 * 1024.0))
              << std::endl;
}
#endif  // defined(HAVE_EIGEN)

// -----------------------------------------------------------------------------
// Symmetric Matrix
// -----------------------------------------------------------------------------
template <typename SymmetricMatrixType, typename GeneralMatrixType>
void benchSymmetricMatrix(const TlMatrixObject::index_type dim) {
    TlTime time;
    time.start();
    SymmetricMatrixType A(dim);
    SymmetricMatrixType B(dim);
    const double createTime = time.getElapseTime();
    time.reset();
    time.start();

    setupSymmetricMatrix(dim, &A);
    setupSymmetricMatrix(dim, &B);
    const double setupTime = time.getElapseTime();
    time.reset();
    time.start();

    GeneralMatrixType C = A * B;
    const double mulTime = time.getElapseTime();

    // FLOP
    double FLOP = 0.0;
    {
        const double row1 = A.getNumOfRows();
        const double col1 = A.getNumOfCols();
        // const double row2 = B.getNumOfRows(); // col1 == row2
        const double col2 = B.getNumOfCols();
        FLOP = (col1 + (col1 - 1)) * row1 * col2;
        // std::cout << "FLOP: " << FLOP << std::endl;
    }

    std::cout << TlUtils::format("create: %8.3e sec", createTime) << std::endl;
    std::cout << TlUtils::format("setup : %8.3e sec", setupTime) << std::endl;
    std::cout << TlUtils::format("mul   : %8.3e sec (%e GFLOPS)", mulTime,
                                 FLOP / (mulTime * 1024.0 * 1024.0 * 1024.0))
              << std::endl;
}

#if defined(HAVE_EIGEN)
template <typename SymmetricMatrixType, typename GeneralMatrixType>
void benchSymmetricMatrixWithEigen(const TlMatrixObject::index_type dim) {
    TlTime time;
    time.start();
    SymmetricMatrixType A(dim);
    SymmetricMatrixType B(dim);
    const double createTime = time.getElapseTime();
    time.reset();
    time.start();

    setupSymmetricMatrixWithEigen(dim, &A);
    setupSymmetricMatrixWithEigen(dim, &B);
    const double setupTime = time.getElapseTime();
    time.reset();
    time.start();

    GeneralMatrixType C = A * B;
    const double mulTime = time.getElapseTime();

    // FLOP
    double FLOP = 0.0;
    {
        const double row1 = A.getNumOfRows();
        const double col1 = A.getNumOfCols();
        // const double row2 = B.getNumOfRows(); // col1 == row2
        const double col2 = B.getNumOfCols();
        FLOP = (col1 + (col1 - 1)) * row1 * col2;
        // std::cout << "FLOP: " << FLOP << std::endl;
    }

    std::cout << TlUtils::format("create: %8.3e sec", createTime) << std::endl;
    std::cout << TlUtils::format("setup : %8.3e sec", setupTime) << std::endl;
    std::cout << TlUtils::format("mul   : %8.3e sec (%e GFLOPS)", mulTime,
                                 FLOP / (mulTime * 1024.0 * 1024.0 * 1024.0))
              << std::endl;
}
#endif  // defined(HAVE_EIGEN)

// -----------------------------------------------------------------------------
// prepare
// -----------------------------------------------------------------------------
void showHelp() {
    std::cerr << "a matrix(multiplication) benchmark of linear algebra packages"
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

    std::cout << "dim: " << dim << std::endl;

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

    // general matrix
    if (onCPU) {
        std::cout << ">>>> BLAS@GeneralMatrix" << std::endl;
        benchGeneralMatrix<TlDenseGeneralMatrix_Lapack>(dim);

#ifdef HAVE_EIGEN
        std::cout << ">>>> Eigen@GeneralMatrix" << std::endl;
        benchGeneralMatrix<TlDenseGeneralMatrix_Eigen>(dim);
#endif  // HAVE_EIGEN
    }

#ifdef HAVE_VIENNACL
    if (onGPU) {
        // std::cout << ">>>> ViennaCL@GeneralMatrix" << std::endl;
        // benchGeneralMatrix<TlDenseGeneralMatrix_ViennaCL>(dim);

        std::cout << ">>>> ViennaCL(Eigen)@GeneralMatrix" << std::endl;
        benchGeneralMatrixWithEigen<TlDenseGeneralMatrix_ViennaCL>(dim);
    }
#endif  // HAVE_VIENNACL

    // symmetric matrix
    if (onCPU) {
        std::cout << ">>>> BLAS@SymmetricMatrix" << std::endl;
        benchSymmetricMatrix<TlDenseSymmetricMatrix_Lapack,
                             TlDenseGeneralMatrix_Lapack>(dim);

#ifdef HAVE_EIGEN
        std::cout << ">>>> Eigen@SymmetriclMatrix" << std::endl;
        benchSymmetricMatrix<TlDenseSymmetricMatrix_Eigen,
                             TlDenseGeneralMatrix_Eigen>(dim);
#endif  // HAVE_EIGEN
    }

#ifdef HAVE_VIENNACL
    if (onGPU) {
        // std::cout << ">>>> ViennaCL@SymmetricMatrix" << std::endl;
        // benchSymmetricMatrix<TlDenseSymmetricMatrix_ViennaCL,
        //                      TlDenseGeneralMatrix_ViennaCL>(dim);

        std::cout << ">>>> ViennaCL(Eigen)@SymmetricMatrix" << std::endl;
        benchSymmetricMatrixWithEigen<TlDenseSymmetricMatrix_ViennaCL,
                                      TlDenseGeneralMatrix_ViennaCL>(dim);
    }
#endif  // HAVE_VIENNACL

    std::cout << "done." << std::endl;
    return EXIT_SUCCESS;
}
