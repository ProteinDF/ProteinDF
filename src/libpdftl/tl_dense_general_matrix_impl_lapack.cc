#include "tl_dense_general_matrix_impl_lapack.h"
#include <iostream>
#include "lapack.h"
#include "tl_dense_general_matrix_lapack.h"
#include "tl_dense_symmetric_matrix_impl_lapack.h"
#include "tl_dense_vector_impl_lapack.h"

TlDenseGeneralMatrix_ImplLapack::TlDenseGeneralMatrix_ImplLapack(
    const TlMatrixObject::index_type row, const TlMatrixObject::index_type col)
    : row_(row), col_(col), matrix_(NULL) {
  this->initialize();
};

TlDenseGeneralMatrix_ImplLapack::TlDenseGeneralMatrix_ImplLapack(
    const TlDenseGeneralMatrix_ImplLapack& rhs)
    : row_(rhs.getNumOfRows()), col_(rhs.getNumOfCols()), matrix_(NULL) {
  this->initialize(false);
  std::copy(rhs.matrix_, rhs.matrix_ + this->getNumOfElements(), this->matrix_);
}

TlDenseGeneralMatrix_ImplLapack::TlDenseGeneralMatrix_ImplLapack(
    const TlDenseSymmetricMatrix_ImplLapack& rhs)
    : row_(rhs.getNumOfRows()), col_(rhs.getNumOfCols()), matrix_(NULL) {
  this->initialize();

  const TlMatrixObject::index_type dim = rhs.getNumOfRows();

  const char UPLO = 'U';
  const int N = dim;
  const int LDA = std::max(1, N);
  int INFO;
  dtpttr_(&UPLO, &N, rhs.matrix_, this->matrix_, &LDA, &INFO);
  if (INFO != 0) {
    this->log_.critical(
        TlUtils::format("program error: %d@%s", __LINE__, __FILE__));
  }

  for (TlMatrixObject::index_type r = 0; r < dim; r++) {
    for (TlMatrixObject::index_type c = 0; c < r; c++) {
      this->set(r, c, this->get(c, r));
    }
  }
}

TlDenseGeneralMatrix_ImplLapack::~TlDenseGeneralMatrix_ImplLapack() {
  delete[] this->matrix_;
  this->matrix_ = NULL;
}

// ---------------------------------------------------------------------------
// static
// ---------------------------------------------------------------------------
TlDenseGeneralMatrix_ImplLapack TlDenseGeneralMatrix_ImplLapack::E(
    const TlMatrixObject::index_type dim) {
  assert(dim > 0);
  TlDenseGeneralMatrix_ImplLapack E(dim);
  for (TlMatrixObject::index_type i = 0; i < dim; ++i) {
    E.set(i, i, 1.0);
  }

  return E;
}

// ---------------------------------------------------------------------------
// properties
// ---------------------------------------------------------------------------
TlMatrixObject::index_type TlDenseGeneralMatrix_ImplLapack::getNumOfRows() const {
  return this->row_;
}

TlMatrixObject::index_type TlDenseGeneralMatrix_ImplLapack::getNumOfCols() const {
  return this->col_;
}

void TlDenseGeneralMatrix_ImplLapack::resize(
    const TlMatrixObject::index_type newRow,
    const TlMatrixObject::index_type newCol) {
  assert((newRow > 0) && (newCol > 0));
  const TlDenseGeneralMatrix_ImplLapack oldMatrix(*this);

  // destroy object
  delete[] this->matrix_;
  this->matrix_ = NULL;

  // initialize
  this->row_ = newRow;
  this->col_ = newCol;
  this->initialize(true);

  // copy old data
  const TlMatrixObject::index_type maxRowForCopy =
      std::min(oldMatrix.getNumOfRows(), newRow);
  const TlMatrixObject::index_type maxColForCopy =
      std::min(oldMatrix.getNumOfCols(), newCol);
#pragma omp parallel for
  for (TlMatrixObject::index_type c = 0; c < maxColForCopy; ++c) {
    for (TlMatrixObject::index_type r = 0; r < maxRowForCopy; ++r) {
      this->set(r, c, oldMatrix.get(r, c));
    }
  }
}

double TlDenseGeneralMatrix_ImplLapack::get(
    const TlMatrixObject::index_type row,
    const TlMatrixObject::index_type col) const {
  return this->matrix_[this->index(row, col)];
}

void TlDenseGeneralMatrix_ImplLapack::set(const TlMatrixObject::index_type row,
                                        const TlMatrixObject::index_type col,
                                        const double value) {
  const TlMatrixObject::size_type index = this->index(row, col);

#pragma omp critical(TlDenseGeneralMatrix_ImplLapack__set)
  { this->matrix_[index] = value; }
}

void TlDenseGeneralMatrix_ImplLapack::add(const TlMatrixObject::index_type row,
                                        const TlMatrixObject::index_type col,
                                        const double value) {
  const TlMatrixObject::size_type index = this->index(row, col);

#pragma omp atomic
  this->matrix_[index] += value;
}

// ---------------------------------------------------------------------------
// operators
// ---------------------------------------------------------------------------
TlDenseGeneralMatrix_ImplLapack& TlDenseGeneralMatrix_ImplLapack::operator=(
    const TlDenseGeneralMatrix_ImplLapack& rhs) {
  if (this != &rhs) {
    delete[] this->matrix_;
    this->matrix_ = NULL;

    this->row_ = rhs.getNumOfRows();
    this->col_ = rhs.getNumOfCols();
    this->initialize(false);
    const std::size_t size = this->getNumOfElements();
    std::copy(rhs.matrix_, rhs.matrix_ + size, this->matrix_);
  }

  return (*this);
}

TlDenseGeneralMatrix_ImplLapack TlDenseGeneralMatrix_ImplLapack::operator*(
    const TlDenseGeneralMatrix_ImplLapack& rhs) const {
  const TlMatrixObject::index_type row1 = this->getNumOfRows();
  const TlMatrixObject::index_type col1 = this->getNumOfCols();
  const TlMatrixObject::index_type row2 = rhs.getNumOfRows();
  const TlMatrixObject::index_type col2 = rhs.getNumOfCols();
  assert(col1 == row2);
  TlDenseGeneralMatrix_ImplLapack answer(row1, col2);

  const char TRANSA = 'N';
  const char TRANSB = 'N';
  const double alpha = 1.0;
  const double beta = 0.0;
  const int M = row1;
  const int N = col2;
  const int K = col1;
  const int LDA = M;
  const int LDB = K;
  const int LDC = M;

  // DGEMM  performs one of the matrix-matrix operations
  // C := alpha*op( A )*op( B ) + beta*C,
  // where  op( X ) is one of op( X ) = X   or   op( X ) = X',
  // alpha and beta are scalars, and A, B and C are matrices, with op( A )
  // an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix
  dgemm_(&TRANSA, &TRANSB, &M, &N, &K, &alpha,
         const_cast<TlDenseGeneralMatrix_ImplLapack*>(this)->matrix_, &LDA,
         const_cast<TlDenseGeneralMatrix_ImplLapack&>(rhs).matrix_, &LDB, &beta,
         answer.matrix_, &LDC);

  return answer;
}

TlDenseGeneralMatrix_ImplLapack& TlDenseGeneralMatrix_ImplLapack::operator+=(
    const TlDenseGeneralMatrix_ImplLapack& rhs) {
  const TlMatrixObject::index_type row1 = this->getNumOfRows();
  const TlMatrixObject::index_type col1 = this->getNumOfCols();
  const TlMatrixObject::index_type row2 = rhs.getNumOfRows();
  const TlMatrixObject::index_type col2 = rhs.getNumOfCols();
  assert(row1 == row2);
  assert(col1 == col2);

  const int n = this->getNumOfElements();
  const double alpha = 1.0;
  const int incx = 1;
  const int incy = 1;
  daxpy_(&n, &alpha, rhs.matrix_, &incx, this->matrix_, &incy);

  return *this;
}

TlDenseGeneralMatrix_ImplLapack& TlDenseGeneralMatrix_ImplLapack::operator-=(
    const TlDenseGeneralMatrix_ImplLapack& rhs) {
  const TlMatrixObject::index_type row1 = this->getNumOfRows();
  const TlMatrixObject::index_type col1 = this->getNumOfCols();
  const TlMatrixObject::index_type row2 = rhs.getNumOfRows();
  const TlMatrixObject::index_type col2 = rhs.getNumOfCols();
  assert(row1 == row2);
  assert(col1 == col2);

  const int n = this->getNumOfElements();
  const double alpha = -1.0;
  const int incx = 1;
  const int incy = 1;
  daxpy_(&n, &alpha, rhs.matrix_, &incx, this->matrix_, &incy);

  return *this;
}

TlDenseGeneralMatrix_ImplLapack& TlDenseGeneralMatrix_ImplLapack::operator*=(
    const double coef) {
  const int n = this->getNumOfElements();
  const double a = -1.0;
  const int incx = 1;
  dscal_(&n, &coef, this->matrix_, &incx);

  return *this;
}

TlDenseGeneralMatrix_ImplLapack& TlDenseGeneralMatrix_ImplLapack::operator/=(
    const double coef) {
  return this->operator*=(1.0 / coef);
}

TlDenseGeneralMatrix_ImplLapack& TlDenseGeneralMatrix_ImplLapack::operator*=(
    const TlDenseGeneralMatrix_ImplLapack& rhs) {
  const TlDenseGeneralMatrix_ImplLapack answer = *this * rhs;
  *this = answer;

  return *this;
}

// ---------------------------------------------------------------------------
// operations
// ---------------------------------------------------------------------------
const TlDenseGeneralMatrix_ImplLapack& TlDenseGeneralMatrix_ImplLapack::dotInPlace(
    const TlDenseGeneralMatrix_ImplLapack& rhs) {
  assert(this->getNumOfRows() == rhs.getNumOfRows());
  assert(this->getNumOfCols() == rhs.getNumOfCols());

  const TlMatrixObject::size_type length = this->getNumOfElements();
  for (TlMatrixObject::size_type i = 0; i < length; ++i) {
    this->matrix_[i] *= rhs.matrix_[i];
  }

  return *this;
}

TlDenseGeneralMatrix_ImplLapack TlDenseGeneralMatrix_ImplLapack::transpose() const {
  const TlMatrixObject::index_type row = this->getNumOfRows();
  const TlMatrixObject::index_type col = this->getNumOfCols();

  TlDenseGeneralMatrix_ImplLapack E = TlDenseGeneralMatrix_ImplLapack::E(row);

  TlDenseGeneralMatrix_ImplLapack answer(col, row);
  const char transa = 'N';
  const char transb = 'T';
  const int m = row;
  const int n = col;
  const int k = row;
  const double alpha = 1.0;
  const int ldA = m;
  const int ldB = k;
  const double beta = 0.0;
  const int ldC = col;

  dgemm_(&transa, &transb, &m, &n, &k, &alpha, E.matrix_, &ldA, this->matrix_,
         &ldB, &beta, answer.matrix_, &ldC);

  return answer;
}

TlDenseGeneralMatrix_ImplLapack TlDenseGeneralMatrix_ImplLapack::inverse() const {
  TlDenseGeneralMatrix_ImplLapack answer = *this;
  const int M = answer.getNumOfRows();
  const int N = answer.getNumOfCols();

  const int LDA = std::max(1, M);
  int* IPIV = new int[std::min(M, N)];
  int INFO = 0;

  double* A = answer.matrix_;
  int LWORK = std::max(1, N);
  double* WORK = new double[LWORK];

  dgetrf_(&M, &N, A, &LDA, IPIV, &INFO);
  if (INFO == 0) {
    dgetri_(&N, A, &LDA, IPIV, WORK, &LWORK, &INFO);
    if (INFO == 0) {
      ;
    } else {
      this->log_.critical(
          TlUtils::format("dgetri() failed. return code = %d (%d@%s)", INFO,
                          __LINE__, __FILE__));
    }
  } else {
    this->log_.critical(TlUtils::format(
        "dgetrf() failed. return code = %d (%d@%s)", INFO, __LINE__, __FILE__));
  }

  std::swap(answer.row_, answer.col_);

  delete[] WORK;
  WORK = NULL;

  delete[] IPIV;
  IPIV = NULL;

  return answer;
}

// ---------------------------------------------------------------------------
// protected
// ---------------------------------------------------------------------------
void TlDenseGeneralMatrix_ImplLapack::initialize(bool clearIfNeeded) {
  const TlMatrixObject::index_type row = this->getNumOfRows();
  const TlMatrixObject::index_type col = this->getNumOfCols();

  const std::size_t size = this->getNumOfElements();
  if (size > 0) {
    try {
      this->matrix_ = new double[size];
    } catch (std::bad_alloc& ba) {
      this->log_.critical(
          TlUtils::format("bad_alloc caught: %s: row=%d, col=%d, size=%ld",
                          ba.what(), row, col, size));
      throw;
    } catch (...) {
      this->log_.critical("unknown error.");
      throw;
    }
    assert(this->matrix_ != NULL);

    if (clearIfNeeded) {
      std::fill(this->matrix_, this->matrix_ + size, 0.0);
    }
  }
}

TlMatrixObject::size_type TlDenseGeneralMatrix_ImplLapack::getNumOfElements()
    const {
  const TlMatrixObject::size_type row = this->getNumOfRows();
  const TlMatrixObject::size_type col = this->getNumOfCols();
  const TlMatrixObject::size_type elements = row * col;

  return elements;
}

TlMatrixObject::size_type TlDenseGeneralMatrix_ImplLapack::index(
    TlMatrixObject::index_type row, TlMatrixObject::index_type col) const {
  // Column-major order
  // see. https://en.wikipedia.org/wiki/Row-_and_column-major_order
  assert((0 <= row) && (row < this->getNumOfRows()));
  assert((0 <= col) && (col < this->getNumOfCols()));

  const TlMatrixObject::size_type index = row + this->getNumOfRows() * col;
  return index;
}

// ---------------------------------------------------------------------------
// others
// ---------------------------------------------------------------------------
TlDenseVector_ImplLapack operator*(const TlDenseGeneralMatrix_ImplLapack& mat,
                                   const TlDenseVector_ImplLapack& vec) {
  TlLogging& logger = TlLogging::getInstance();
  if (mat.getNumOfCols() != vec.getSize()) {
    logger.critical(TlUtils::format("size mismatch: %d != %d (%d@%s)",
                                    mat.getNumOfCols(), vec.getSize(), __LINE__,
                                    __FILE__));
  }

  TlDenseVector_ImplLapack answer(vec.getSize());
  const char TRANS = 'N';
  const int M = mat.getNumOfRows();
  const int N = mat.getNumOfCols();
  const double alpha = 1.0;
  const int LDA = std::max(1, M);
  const int INCX = 1;
  const double beta = 1.0;
  const int INCY = 1;

  // *  DGEMV  performs one of the matrix-vector operations
  // *
  // *     y := alpha*A*x + beta*y,   or   y := alpha*A'*x + beta*y,
  // *
  // *  where alpha and beta are scalars, x and y are vectors and A is an
  // *  m by n matrix.
  dgemv_(&TRANS, &M, &N, &alpha, mat.matrix_, &LDA, vec.vector_, &INCX, &beta,
         answer.vector_, &INCY);

  return answer;
}

TlDenseVector_ImplLapack operator*(const TlDenseVector_ImplLapack& vec,
                                   const TlDenseGeneralMatrix_ImplLapack& mat) {
  TlLogging& logger = TlLogging::getInstance();
  if (mat.getNumOfRows() != vec.getSize()) {
    logger.critical(TlUtils::format("size mismatch: %d != %d (%d@%s)",
                                    mat.getNumOfCols(), vec.getSize(), __LINE__,
                                    __FILE__));
  }

  TlDenseVector_ImplLapack answer(vec.getSize());
  const char TRANS = 'T';
  const int M = mat.getNumOfRows();
  const int N = mat.getNumOfCols();
  const double alpha = 1.0;
  const int LDA = std::max(1, M);
  const int INCX = 1;
  const double beta = 1.0;
  const int INCY = 1;

  // *  DGEMV  performs one of the matrix-vector operations
  // *
  // *     y := alpha*A*x + beta*y,   or   y := alpha*A'*x + beta*y,
  // *
  // *  where alpha and beta are scalars, x and y are vectors and A is an
  // *  m by n matrix.
  dgemv_(&TRANS, &M, &N, &alpha, mat.matrix_, &LDA, vec.vector_, &INCX, &beta,
         answer.vector_, &INCY);

  return answer;
}
