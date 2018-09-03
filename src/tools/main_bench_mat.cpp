#include <stdio.h>
#include <time.h>
#include <cmath>
#include <cstdlib>
#include <iostream>

#include "TlGetopt.h"
#include "TlTime.h"
#include "TlUtils.h"
#include "tl_dense_general_matrix_lapack.h"
#include "tl_dense_symmetric_matrix_lapack.h"
#include "tl_dense_vector_lapack.h"

#ifdef _OPENMP
#include <omp.h>
#endif  // _OPENMP

#ifdef HAVE_EIGEN
#include "tl_dense_general_matrix_eigen.h"
#include "tl_dense_symmetric_matrix_eigen.h"
#include "tl_dense_vector_eigen.h"
#endif // HAVE_EIGEN

#ifdef HAVE_VIENNACL
#include "tl_dense_general_matrix_viennacl.h"
#include "tl_dense_symmetric_matrix_viennacl.h"
#include "tl_dense_vector_viennacl.h"

#include "tl_viennacl.h"
#include "viennacl/ocl/backend.hpp"
#include "viennacl/ocl/context.hpp"
#include "viennacl/ocl/enqueue.hpp"
#include "viennacl/ocl/platform.hpp"
#endif  // HAVE_VIENNACL

// -----------------------------------------------------------------------------
// General Matrix
// -----------------------------------------------------------------------------
template <typename DenseMatrixType>
void setupGeneralMatrix(const TlMatrixObject::index_type dim,
                        DenseMatrixType* pMat) {
  for (TlMatrixObject::index_type r = 0; r < dim; ++r) {
    for (TlMatrixObject::index_type c = 0; c < dim; ++c) {
      const double v = double(rand()) / double((RAND_MAX));
      pMat->set(r, c, v);
    }
  }
}

#ifdef HAVE_EIGEN
template <typename DenseMatrixType>
void setupGeneralMatrixWithEigen(const TlMatrixObject::index_type dim,
                                 DenseMatrixType* pMat) {
  TlDenseGeneralMatrix_Eigen tmp(dim, dim);
  for (TlMatrixObject::index_type r = 0; r < dim; ++r) {
    for (TlMatrixObject::index_type c = 0; c < dim; ++c) {
      const double v = double(rand()) / double((RAND_MAX));
      tmp.set(r, c, v);
    }
  }

  *pMat = tmp;
}
#endif // HAVE_EIGEN

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
#endif // defined(HAVE_EIGEN)


// -----------------------------------------------------------------------------
// Symmetric Matrix
// -----------------------------------------------------------------------------
template <typename SymmetricMatrixType>
void setupSymmetricMatrix(const TlMatrixObject::index_type dim,
                          SymmetricMatrixType* pMat) {
  for (TlMatrixObject::index_type r = 0; r < dim; ++r) {
    for (TlMatrixObject::index_type c = 0; c < r; ++c) {
      const double v = double(rand()) / double((RAND_MAX));
      pMat->set(r, c, v);
    }
    pMat->set(r, r, double(r));
  }
}

#ifdef HAVE_EIGEN
template <typename SymmetricMatrixType>
void setupSymmetricMatrixWithEigen(const TlMatrixObject::index_type dim,
                                   SymmetricMatrixType* pMat) {
  TlDenseSymmetricMatrix_Eigen tmp(dim);
  for (TlMatrixObject::index_type r = 0; r < dim; ++r) {
    for (TlMatrixObject::index_type c = 0; c < r; ++c) {
      const double v = double(rand()) / double((RAND_MAX));
      tmp.set(r, c, v);
    }
    tmp.set(r, r, double(r));
  }

  *pMat = tmp;
}
#endif // HAVE_EIGEN

template <typename SymmetricMatrixType>
void setupSymmetricMatrixForDiagonal(const TlMatrixObject::index_type dim,
                                     SymmetricMatrixType* pMat) {
  pMat->resize(dim);
  for (TlMatrixObject::index_type r = 0; r < dim; ++r) {
    for (TlMatrixObject::index_type c = 0; c < r; ++c) {
      const double v = std::sqrt(0.3 * double(r) * double(c));
      pMat->set(r, c, v);
    }
    pMat->set(r, r, double(r));
  }
}

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
#endif // defined(HAVE_EIGEN)


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

  std::cout << TlUtils::format("create:   %8.3e sec", createTime) << std::endl;
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
  std::cout << TlUtils::format("create:   %8.3e sec", createTime) << std::endl;
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
#endif // HAVE_EIGEN

// -----------------------------------------------------------------------------
// prepare
// -----------------------------------------------------------------------------
void initRand() { srand((unsigned int)time(NULL)); }

void showHelp() {
  std::cerr << "a benchmark of linear algebra packages" << std::endl;
  std::cerr << "OPTIONS:" << std::endl;
  std::cerr << "-d <device id>: switch device id" << std::endl;
  std::cerr << "-s <size>: matrix size" << std::endl;
  std::cerr << "-g: test GPU code only" << std::endl;
}

// -----------------------------------------------------------------------------
// MAIN
// -----------------------------------------------------------------------------
int main(int argc, char* argv[]) {
  TlGetopt opt(argc, argv, "ehd:s:g");

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
  TlMatrixObject::index_type dim = 10240;
  if (!opt["s"].empty()) {
    dim = std::atoi(opt["s"].c_str());
  }

  bool gpuOnly = false;
  if (opt["g"] == "defined") {
    gpuOnly = true;
  }

  bool eigOnly = false;
  if (opt["e"] == "defined") {
    eigOnly = true;
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
  if (eigOnly != true) {
    if (gpuOnly != true) {
      std::cout << ">>>> BLAS@GeneralMatrix" << std::endl;
      benchGeneralMatrix<TlDenseGeneralMatrix_Lapack>(dim);

#ifdef HAVE_EIGEN
      std::cout << ">>>> Eigen@GeneralMatrix" << std::endl;
      benchGeneralMatrix<TlDenseGeneralMatrix_Eigen>(dim);
#endif // HAVE_EIGEN
    }
#ifdef HAVE_VIENNACL
    // std::cout << ">>>> ViennaCL@GeneralMatrix" << std::endl;
    // benchGeneralMatrix<TlDenseGeneralMatrix_ViennaCL>(dim);

    std::cout << ">>>> ViennaCL(Eigen)@GeneralMatrix" << std::endl;
    benchGeneralMatrixWithEigen<TlDenseGeneralMatrix_ViennaCL>(dim);
#endif  // HAVE_VIENNACL

    // symmetric matrix
    if (gpuOnly != true) {
      std::cout << ">>>> BLAS@SymmetricMatrix" << std::endl;
      benchSymmetricMatrix<TlDenseSymmetricMatrix_Lapack,
                           TlDenseGeneralMatrix_Lapack>(dim);

#ifdef HAVE_EIGEN                      
      std::cout << ">>>> Eigen@SymmetriclMatrix" << std::endl;
      benchSymmetricMatrix<TlDenseSymmetricMatrix_Eigen,
                           TlDenseGeneralMatrix_Eigen>(dim);
#endif // HAVE_EIGEN
    }
#ifdef HAVE_VIENNACL
    // std::cout << ">>>> ViennaCL@SymmetricMatrix" << std::endl;
    // benchSymmetricMatrix<TlDenseSymmetricMatrix_ViennaCL,
    //                      TlDenseGeneralMatrix_ViennaCL>(dim);

    std::cout << ">>>> ViennaCL(Eigen)@SymmetricMatrix" << std::endl;
    benchSymmetricMatrixWithEigen<TlDenseSymmetricMatrix_ViennaCL,
                                  TlDenseGeneralMatrix_ViennaCL>(dim);
#endif  // HAVE_VIENNACL
  }

  // diagonal
  if (gpuOnly != true) {
    std::cout << ">>>> BLAS@SymmetricMatrix (diagonal)" << std::endl;
    benchDiagonal<TlDenseSymmetricMatrix_Lapack, TlDenseGeneralMatrix_Lapack,
                  TlDenseVector_Lapack>(dim);

#ifdef HAVE_EIGEN                  
    std::cout << ">>>> Eigen@SymmetricMatrix (diagonal)" << std::endl;
    benchDiagonal<TlDenseSymmetricMatrix_Eigen, TlDenseGeneralMatrix_Eigen,
                  TlDenseVector_Eigen>(dim);
#endif // HAVE_EIGEN
  }

#ifdef HAVE_VIENNACL
  std::cout << ">>>> ViennaCL@SymmetricMatrix (diagonal)" << std::endl;
  benchDiagonalWithEigen<TlDenseSymmetricMatrix_ViennaCL,
                         TlDenseGeneralMatrix_ViennaCL, TlDenseVector_ViennaCL>(
      dim);
#endif  // HAVE_VIENNACL

  std::cout << "done." << std::endl;
  return EXIT_SUCCESS;
}
