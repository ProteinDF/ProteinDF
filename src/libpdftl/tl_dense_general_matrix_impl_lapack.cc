#include "tl_dense_general_matrix_impl_lapack.h"

#include <cassert>
#include <iostream>

#include "TlUtils.h"
#include "lapack.h"
#include "tl_dense_general_matrix_lapack.h"
#include "tl_dense_symmetric_matrix_impl_lapack.h"
#include "tl_dense_vector_impl_lapack.h"
#include "tl_dense_vector_impl_scalapack.h"

// ----------------------------------------------------------------------------
// constructor & destructor
// ----------------------------------------------------------------------------
TlDenseGeneralMatrix_ImplLapack::TlDenseGeneralMatrix_ImplLapack(const TlMatrixObject::index_type row,
                                                                 const TlMatrixObject::index_type col,
                                                                 double const* const pBuf)
    : row_(row), col_(col), matrix_(NULL) {
    this->initialize();

    if (pBuf != NULL) {
        this->vtr2mat(pBuf);
    }
};

TlDenseGeneralMatrix_ImplLapack::TlDenseGeneralMatrix_ImplLapack(const TlDenseGeneralMatrix_ImplLapack& rhs)
    : row_(rhs.getNumOfRows()), col_(rhs.getNumOfCols()), matrix_(NULL) {
    this->initialize(false);
    std::copy(rhs.matrix_, rhs.matrix_ + this->getNumOfElements(), this->matrix_);
}

TlDenseGeneralMatrix_ImplLapack::TlDenseGeneralMatrix_ImplLapack(const TlDenseSymmetricMatrix_ImplLapack& rhs)
    : row_(rhs.getNumOfRows()), col_(rhs.getNumOfCols()), matrix_(NULL) {
    this->initialize();

    const TlMatrixObject::index_type dim = rhs.getNumOfRows();

    const char UPLO = 'U';
    const int N = dim;
    const int LDA = std::max(1, N);
    int INFO;
    dtpttr_(&UPLO, &N, rhs.matrix_, this->matrix_, &LDA, &INFO);
    if (INFO != 0) {
        this->log_.critical(TlUtils::format("program error: %d@%s", __LINE__, __FILE__));
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

TlDenseGeneralMatrix_ImplLapack::operator std::vector<double>() const {
    std::vector<double> answer(this->matrix_, this->matrix_ + this->getNumOfElements());
    return answer;
}

// ---------------------------------------------------------------------------
// static
// ---------------------------------------------------------------------------
TlDenseGeneralMatrix_ImplLapack TlDenseGeneralMatrix_ImplLapack::E(const TlMatrixObject::index_type dim) {
    assert(dim > 0);
    TlDenseGeneralMatrix_ImplLapack E(dim, dim);
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

void TlDenseGeneralMatrix_ImplLapack::resize(const TlMatrixObject::index_type newRow,
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
    const TlMatrixObject::index_type maxRowForCopy = std::min(oldMatrix.getNumOfRows(), newRow);
    const TlMatrixObject::index_type maxColForCopy = std::min(oldMatrix.getNumOfCols(), newCol);
#pragma omp parallel for
    for (TlMatrixObject::index_type c = 0; c < maxColForCopy; ++c) {
        for (TlMatrixObject::index_type r = 0; r < maxRowForCopy; ++r) {
            this->set(r, c, oldMatrix.get(r, c));
        }
    }
}

double TlDenseGeneralMatrix_ImplLapack::get(const TlMatrixObject::index_type row,
                                            const TlMatrixObject::index_type col) const {
    return this->matrix_[this->index(row, col)];
}

void TlDenseGeneralMatrix_ImplLapack::set(const TlMatrixObject::index_type row, const TlMatrixObject::index_type col,
                                          const double value) {
    const TlMatrixObject::size_type index = this->index(row, col);

#pragma omp critical(TlDenseGeneralMatrix_ImplLapack__set)
    { this->matrix_[index] = value; }
}

void TlDenseGeneralMatrix_ImplLapack::add(const TlMatrixObject::index_type row, const TlMatrixObject::index_type col,
                                          const double value) {
    const TlMatrixObject::size_type index = this->index(row, col);

#pragma omp atomic
    this->matrix_[index] += value;
}

TlMatrixObject::index_type TlDenseGeneralMatrix_ImplLapack::getRowVector(const TlMatrixObject::index_type row,
                                                                         const TlMatrixObject::index_type length,
                                                                         double* pBuf) const {
    TlMatrixObject::index_type copiedLength = std::min(length, this->getNumOfCols());
    const TlMatrixObject::size_type base = this->index(row, 0);
    const TlMatrixObject::size_type rows = this->getNumOfRows();

#pragma omp parallel
    for (TlMatrixObject::index_type c = 0; c < copiedLength; ++c) {
        pBuf[c] = this->matrix_[base + rows * c];
    }

    return copiedLength;
}

TlMatrixObject::index_type TlDenseGeneralMatrix_ImplLapack::getColVector(const TlMatrixObject::index_type col,
                                                                         const TlMatrixObject::index_type length,
                                                                         double* pBuf) const {
    TlMatrixObject::index_type copiedLength = std::min(length, this->getNumOfRows());
    const TlMatrixObject::size_type base = this->index(0, col);

    std::copy(this->matrix_ + base, this->matrix_ + base + copiedLength, pBuf);

    return copiedLength;
}

TlMatrixObject::index_type TlDenseGeneralMatrix_ImplLapack::setRowVector(const TlMatrixObject::index_type row,
                                                                         const TlMatrixObject::index_type length,
                                                                         const double* pBuf) {
    TlMatrixObject::index_type copiedLength = std::min(length, this->getNumOfCols());
    const TlMatrixObject::size_type base = this->index(row, 0);
    const TlMatrixObject::size_type rows = this->getNumOfRows();

#pragma omp parallel
    for (TlMatrixObject::index_type c = 0; c < copiedLength; ++c) {
        this->matrix_[base + rows * c] = pBuf[c];
    }

    return copiedLength;
}

TlMatrixObject::index_type TlDenseGeneralMatrix_ImplLapack::setColVector(const TlMatrixObject::index_type col,
                                                                         const TlMatrixObject::index_type length,
                                                                         const double* pBuf) {
    TlMatrixObject::index_type copiedLength = std::min(length, this->getNumOfRows());
    const TlMatrixObject::size_type base = this->index(0, col);

    std::copy(pBuf, pBuf + copiedLength, this->matrix_ + base);

    return copiedLength;
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
    dgemm_(&TRANSA, &TRANSB, &M, &N, &K, &alpha, const_cast<TlDenseGeneralMatrix_ImplLapack*>(this)->matrix_, &LDA,
           const_cast<TlDenseGeneralMatrix_ImplLapack&>(rhs).matrix_, &LDB, &beta, answer.matrix_, &LDC);

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

TlDenseGeneralMatrix_ImplLapack& TlDenseGeneralMatrix_ImplLapack::operator*=(const double coef) {
    const int n = this->getNumOfElements();
    const int incx = 1;
    dscal_(&n, &coef, this->matrix_, &incx);

    return *this;
}

TlDenseGeneralMatrix_ImplLapack& TlDenseGeneralMatrix_ImplLapack::operator/=(const double coef) {
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
double* TlDenseGeneralMatrix_ImplLapack::data() {
    return this->matrix_;
}

const double* TlDenseGeneralMatrix_ImplLapack::data() const {
    return this->matrix_;
}

TlDenseGeneralMatrix_ImplLapack TlDenseGeneralMatrix_ImplLapack::transpose() const {
    const TlMatrixObject::index_type row = this->getNumOfRows();
    const TlMatrixObject::index_type col = this->getNumOfCols();

    TlDenseGeneralMatrix_ImplLapack E = TlDenseGeneralMatrix_ImplLapack::E(row);

    TlDenseGeneralMatrix_ImplLapack answer(col, row);
    const char transa = 'T';
    const char transb = 'N';
    const int m = col;
    const int n = row;
    const int k = row;
    const double alpha = 1.0;
    const int ldA = k;
    const int ldB = k;
    const double beta = 0.0;
    const int ldC = m;

    dgemm_(&transa, &transb, &m, &n, &k, &alpha, this->matrix_, &ldA, E.matrix_, &ldB, &beta, answer.matrix_, &ldC);

    return answer;
}

void TlDenseGeneralMatrix_ImplLapack::transposeInPlace() {
    TlDenseGeneralMatrix_ImplLapack tmp = this->transpose();
    std::swap(*this, tmp);
}

TlDenseGeneralMatrix_ImplLapack TlDenseGeneralMatrix_ImplLapack::dot(const TlDenseGeneralMatrix_ImplLapack& rhs) const {
    TlDenseGeneralMatrix_ImplLapack answer = *this;
    answer.dotInPlace(rhs);

    return answer;
}

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
            this->log_.critical(TlUtils::format("dgetri() failed. return code = %d (%d@%s)", INFO, __LINE__, __FILE__));
        }
    } else {
        this->log_.critical(TlUtils::format("dgetrf() failed. return code = %d (%d@%s)", INFO, __LINE__, __FILE__));
    }

    std::swap(answer.row_, answer.col_);

    delete[] WORK;
    WORK = NULL;

    delete[] IPIV;
    IPIV = NULL;

    return answer;
}

// Ax = B
TlDenseGeneralMatrix_ImplLapack TlDenseGeneralMatrix_ImplLapack::getLeastSquaresSolution(
    const TlDenseGeneralMatrix_ImplLapack& B) const {
    const int M = this->getNumOfRows();
    // if (M != B.getNumOfRows()) {
    //   this->log_.critical(
    //       TlUtils::format("the numbers are not consistent: %d != %d @%s.%d",
    //       M,
    //                       B.getNumOfRows(), __FILE__, __LINE__));
    // }
    const int N = this->getNumOfCols();
    const int NRHS = B.getNumOfCols();

    double* A = new double[M * N];
    std::copy(this->matrix_, this->matrix_ + (M * N), A);
    const int LDA = std::max(1, M);

    TlDenseGeneralMatrix_ImplLapack X = B;
    const int LDB = std::max(1, std::max(M, N));
    X.resize(LDB, NRHS);  // size extended.

    double* S = new double[std::min(M, N)];

    // If RCOND < 0, machine precision is used instead.
    const double RCOND = -1.0;

    int RANK = 0;
    const int LWORK = 3 * std::min(M, N) + std::max(std::max(2 * std::min(M, N), std::max(M, N)), NRHS);
    double* WORK = new double[std::max(1, LWORK)];
    int INFO = 0;

    dgelss_(&M, &N, &NRHS, A, &LDA, X.matrix_, &LDB, S, &RCOND, &RANK, WORK, &LWORK, &INFO);
    if (INFO != 0) {
        if (INFO > 0) {
            this->log_.critical(
                TlUtils::format("%d-th argument had an illegal value. @%s.%d", -INFO, __FILE__, __LINE__));
        } else {
            this->log_.critical(
                TlUtils::format("%d-th off-diagonal elements of an intermediate "
                                "bidiagonal form did not converge to zero. @%s.%d",
                                INFO, __FILE__, __LINE__));
        }
    }

    delete[] S;
    S = NULL;
    delete[] A;
    A = NULL;

    X.resize(N, NRHS);
    return X;
}

// ---------------------------------------------------------------------------
// I/O
// ---------------------------------------------------------------------------
void TlDenseGeneralMatrix_ImplLapack::dump(TlDenseVector_ImplLapack* v) const {
    v->resize(this->getNumOfElements());
    std::copy(this->matrix_, this->matrix_ + this->getNumOfElements(), v->vector_);
}

void TlDenseGeneralMatrix_ImplLapack::restore(const TlDenseVector_ImplLapack& v) {
    const std::size_t copySize = std::min<std::size_t>(this->getNumOfElements(), v.getSize());
    std::copy(v.vector_, v.vector_ + copySize, this->matrix_);
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
                TlUtils::format("bad_alloc caught: %s: row=%d, col=%d, size=%ld", ba.what(), row, col, size));
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

TlMatrixObject::size_type TlDenseGeneralMatrix_ImplLapack::getNumOfElements() const {
    const std::size_t row = this->getNumOfRows();
    const std::size_t col = this->getNumOfCols();
    const std::size_t elements = row * col;

    return elements;
}

TlMatrixObject::size_type TlDenseGeneralMatrix_ImplLapack::index(TlMatrixObject::index_type row,
                                                                 TlMatrixObject::index_type col) const {
    // Column-major order
    // see. https://en.wikipedia.org/wiki/Row-_and_column-major_order
    assert((0 <= row) && (row < this->getNumOfRows()));
    assert((0 <= col) && (col < this->getNumOfCols()));

    const std::size_t r = row;
    const std::size_t c = col;
    const std::size_t index = r + this->getNumOfRows() * c;

    return index;
}

void TlDenseGeneralMatrix_ImplLapack::vtr2mat(double const* const pBuf) {
    std::copy(pBuf, pBuf + this->getNumOfElements(), this->matrix_);
}

// ---------------------------------------------------------------------------
// others
// ---------------------------------------------------------------------------
// DV = DM(G) * DV
TlDenseVector_ImplLapack operator*(const TlDenseGeneralMatrix_ImplLapack& mat, const TlDenseVector_ImplLapack& vec) {
    TlLogging& logger = TlLogging::getInstance();
    if (mat.getNumOfCols() != vec.getSize()) {
        logger.critical(
            TlUtils::format("size mismatch: %d != %d (%d@%s)", mat.getNumOfCols(), vec.getSize(), __LINE__, __FILE__));
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
    dgemv_(&TRANS, &M, &N, &alpha, mat.matrix_, &LDA, vec.vector_, &INCX, &beta, answer.vector_, &INCY);

    return answer;
}

// DV = DV * DM(G)
TlDenseVector_ImplLapack operator*(const TlDenseVector_ImplLapack& vec, const TlDenseGeneralMatrix_ImplLapack& mat) {
    TlLogging& logger = TlLogging::getInstance();
    if (mat.getNumOfRows() != vec.getSize()) {
        logger.critical(
            TlUtils::format("size mismatch: %d != %d (%d@%s)", mat.getNumOfCols(), vec.getSize(), __LINE__, __FILE__));
    }

    TlDenseVector_ImplLapack answer(mat.getNumOfCols());
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
    dgemv_(&TRANS, &M, &N, &alpha, mat.matrix_, &LDA, vec.vector_, &INCX, &beta, answer.vector_, &INCY);

    return answer;
}
