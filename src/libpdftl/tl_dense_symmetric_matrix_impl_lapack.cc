#include <iostream>
#include <cassert>

#include "tl_dense_symmetric_matrix_impl_lapack.h"
#include "lapack.h"
#include "lapack.h"
#include "tl_dense_general_matrix_lapack.h"
#include "tl_dense_vector_impl_lapack.h"
#include "TlUtils.h"

// ---------------------------------------------------------------------------
// constructor & destructor
// ---------------------------------------------------------------------------
TlDenseSymmetricMatrix_ImplLapack::TlDenseSymmetricMatrix_ImplLapack(
    const TlMatrixObject::index_type dim)
    : TlDenseGeneralMatrix_ImplLapack(dim, dim) {}

TlDenseSymmetricMatrix_ImplLapack::TlDenseSymmetricMatrix_ImplLapack(
    const TlDenseSymmetricMatrix_ImplLapack& rhs)
    : TlDenseGeneralMatrix_ImplLapack(rhs.getNumOfRows(), rhs.getNumOfCols()) {
  std::copy(rhs.matrix_, rhs.matrix_ + rhs.getNumOfElements(), this->matrix_);
}

// dim = rows
TlDenseSymmetricMatrix_ImplLapack::TlDenseSymmetricMatrix_ImplLapack(
    const TlDenseGeneralMatrix_ImplLapack& rhs)
    : TlDenseGeneralMatrix_ImplLapack(rhs.getNumOfRows(), rhs.getNumOfRows()) {
  const TlMatrixObject::index_type dim = rhs.getNumOfRows();

  const char UPLO = 'U';
  const int N = dim;
  const int LDA = std::max(1, N);
  int INFO;
  dtrttp_(&UPLO, &N, rhs.matrix_, &LDA, this->matrix_, &INFO);
  if (INFO != 0) {
    this->log_.critical(
        TlUtils::format("program error: %d@%s", __LINE__, __FILE__));
  }
}

TlDenseSymmetricMatrix_ImplLapack::~TlDenseSymmetricMatrix_ImplLapack() {}

// ---------------------------------------------------------------------------
// properties
// ---------------------------------------------------------------------------
void TlDenseSymmetricMatrix_ImplLapack::resize(TlMatrixObject::index_type row,
                                               TlMatrixObject::index_type col) {
  const TlMatrixObject::index_type& dim = row;
  assert(row == col);
  assert(dim > 0);

  TlDenseSymmetricMatrix_ImplLapack oldMatrix(*this);

  this->row_ = row;
  this->col_ = col;
  this->initialize(true);

  const TlMatrixObject::index_type dimForCopy =
      std::min<TlMatrixObject::index_type>(oldMatrix.getNumOfRows(), dim);

#pragma omp parallel for
  for (TlMatrixObject::index_type i = 0; i < dimForCopy; ++i) {
    for (TlMatrixObject::index_type j = 0; j <= i; ++j) {
      this->set(i, j, oldMatrix.get(i, j));
    }
  }
}

// ---------------------------------------------------------------------------
// operators
// ---------------------------------------------------------------------------
TlDenseSymmetricMatrix_ImplLapack& TlDenseSymmetricMatrix_ImplLapack::
operator*=(const double coef) {
  const int n = this->getNumOfElements();
  const int incx = 1;
  dscal_(&n, &coef, this->matrix_, &incx);

  return *this;
}

// const TlDenseGeneralMatrix_ImplLapack
// TlDenseSymmetricMatrix_ImplLapack::operator*=(
//     const TlDenseSymmetricMatrix_ImplLapack& rhs) {
//   TlDenseGeneralMatrix_ImplLapack answer = *this;
//   TlDenseGeneralMatrix_ImplLapack tmp = rhs;
//   answer *= tmp;
//
//   return answer;
// }
// TlDenseSymmetricMatrix_ImplLapack
// TlDenseSymmetricMatrix_ImplLapack::operator*(
//     const TlDenseSymmetricMatrix_ImplLapack& rhs) const {
//   const TlMatrixObject::index_type row1 = this->getNumOfRows();
//   const TlMatrixObject::index_type col1 = this->getNumOfCols();
//   const TlMatrixObject::index_type row2 = rhs.getNumOfRows();
//   const TlMatrixObject::index_type col2 = rhs.getNumOfCols();
//   TlDenseGeneralMatrix_ImplLapack answer(row1, col2);
//
//   char TRANSA = 'N';
//   char TRANSB = 'N';
//   double alpha = 1.0;
//   double beta = 0.0;
//   int M = row1;
//   int N = col2;
//   int K = col1;
//   int LDA = M;
//   int LDB = K;
//   int LDC = M;
//
//   // DGEMM  performs one of the matrix-matrix operations
//   // C := alpha*op( A )*op( B ) + beta*C,
//   // where  op( X ) is one of op( X ) = X   or   op( X ) = X',
//   // alpha and beta are scalars, and A, B and C are matrices, with op( A )
//   // an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix
//   dgemm_(&TRANSA, &TRANSB, &M, &N, &K, &alpha,
//          const_cast<TlDenseGeneralMatrix_ImplLapack*>(this)->matrix_, &LDA,
//          const_cast<TlDenseGeneralMatrix_ImplLapack&>(rhs).matrix_, &LDB,
//          &beta,
//          answer.matrix_, &LDC);
//
//   return answer;
// }

// TlDenseSymmetricMatrix_ImplLapack operator*(
//     const TlDenseSymmetricMatrix_ImplLapack& A,
//     const TlDenseSymmetricMatrix_ImplLapack& B) {
//   TlDenseGeneralMatrix_ImplLapack ans = A;
//   ans *= TlDenseGeneralMatrix_ImplLapack(B);
//
//   return TlDenseSymmetricMatrix_ImplLapack(ans);
// }

// ---------------------------------------------------------------------------
// operations
// ---------------------------------------------------------------------------
TlDenseSymmetricMatrix_ImplLapack TlDenseSymmetricMatrix_ImplLapack::transpose()
    const {
  // do nothing
  return *this;
}

TlDenseSymmetricMatrix_ImplLapack TlDenseSymmetricMatrix_ImplLapack::inverse()
    const {
  TlDenseSymmetricMatrix_ImplLapack answer = *this;

  // (input)
  // 'U':  Upper triangle of A is stored
  // 'L':  Lower triangle of A is stored.
  char UPLO = 'U';

  //(input) The order of the matrix A.  N >= 0.
  const int N = this->getNumOfRows();

  // (input/output) The order of the matrix A.  N >= 0.
  // if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;
  // if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = A(i,j) for j<=i<=n.
  // On exit, if INFO = 0, the triangular factor U or L from the
  // Cholesky factorization A = U**T*U or A = L*L**T, in the same
  // storage format as A.
  double* AP = answer.matrix_;

  // (output) for dsptrf_
  // INTEGER array, dimension (N)
  // Details of the interchanges and the block structure of D.
  // If IPIV(k) > 0, then rows and columns k and IPIV(k) were
  // interchanged and D(k,k) is a 1-by-1 diagonal block.
  // If UPLO = 'U' and IPIV(k) = IPIV(k-1) < 0, then rows and
  // columns k-1 and -IPIV(k) were interchanged and D(k-1:k,k-1:k)
  // is a 2-by-2 diagonal block.  If UPLO = 'L' and IPIV(k) =
  // IPIV(k+1) < 0, then rows and columns k+1 and -IPIV(k) were
  // interchanged and D(k:k+1,k:k+1) is a 2-by-2 diagonal block.
  //
  // (input)  for dsptri_
  // INTEGER array, dimension (N)
  // Details of the interchanges and the block structure of D
  // as determined by DSPTRF.
  int* IPIV = new int[N];

  // (workspace) DOUBLE PRECISION array, dimension (N)
  double* WORK = new double[N];

  // (output)
  // = 0:  successful exit
  // < 0:  if INFO = -i, the i-th argument had an illegal value
  // > 0:  if INFO = i, the leading minor of order i is not
  // positive definite, and the factorization could not be completed.
  int INFO = 0;

  // execute LAPACK
  dsptrf_(&UPLO, &N, AP, IPIV, &INFO);
  if (INFO == 0) {
    dsptri_(&UPLO, &N, AP, IPIV, WORK, &INFO);
    if (INFO != 0) {
      this->log_.critical(TlUtils::format("dsptri() return code = %d. (%d@%s)",
                                          INFO, __LINE__, __FILE__));
    }
  } else {
    this->log_.critical(TlUtils::format("dsptrf() return code = %d. (%d@%s)",
                                        INFO, __LINE__, __FILE__));
  }

  // finalize
  delete[] IPIV;
  IPIV = NULL;
  delete[] WORK;
  WORK = NULL;

  return answer;
}

bool TlDenseSymmetricMatrix_ImplLapack::eig(
    TlDenseVector_ImplLapack* pEigVal,
    TlDenseGeneralMatrix_ImplLapack* pEigVec) const {
  bool answer = false;
  const char JOBZ = 'V';  // 固有値と固有ベクトルを計算する。
  const char UPLO = 'U';
  const int N = this->getNumOfRows();  // 行列Aの次数(N>=0)

  assert(this->getNumOfElements() == N * (N + 1) / 2);
  double* AP = new double[this->getNumOfElements()];
  std::copy(this->matrix_, this->matrix_ + this->getNumOfElements(), AP);

  pEigVal->resize(N);
  // (出力用) INFO=0のとき, Wに固有値が昇順で入る。大きさN
  double* W = pEigVal->vector_;

  const int LDZ = N;
  pEigVec->resize(LDZ, N);
  double* Z = pEigVec->matrix_;

  double* WORK = new double[3 * N];  // (作業/出力用)
  int INFO = 0;  // (出力用) =0: 正常終了, >0: 収束しなかった

  dspev_(&JOBZ, &UPLO, &N, AP, W, Z, &LDZ, WORK, &INFO);

  if (INFO == 0) {
    answer = true;
  } else {
    this->log_.critical(TlUtils::format(
        "dspev calculation faild: INFO=%d (%d@%s)", INFO, __LINE__, __FILE__));
    answer = false;
  }

  delete[] WORK;
  WORK = NULL;
  delete[] AP;
  AP = NULL;

  return answer;
}

// ---------------------------------------------------------------------------
// I/O
// ---------------------------------------------------------------------------
// void TlDenseSymmetricMatrix_ImplLapack::dump(double* buf, const std::size_t
// size) const {
//     const std::size_t copySize = std::min<std::size_t>(size,
//     this->getNumOfElements());
//
//     std::copy(this->matrix_, this->matrix_ + copySize, buf);
// }
//
// void TlDenseSymmetricMatrix_ImplLapack::restore(const double* buf, const
// std::size_t size) {
//     const std::size_t copySize = std::min<std::size_t>(size,
//     this->getNumOfElements());
//
//     std::copy(buf, buf + copySize, this->matrix_);
// }

// ---------------------------------------------------------------------------
// protected
// ---------------------------------------------------------------------------
TlMatrixObject::size_type TlDenseSymmetricMatrix_ImplLapack::getNumOfElements()
    const {
  const TlMatrixObject::size_type dim = this->getNumOfRows();
  assert(dim == this->getNumOfCols());

  const TlMatrixObject::size_type elements = dim * (dim + 1) / 2;

  return elements;
}

TlMatrixObject::size_type TlDenseSymmetricMatrix_ImplLapack::index(
    TlMatrixObject::index_type row, TlMatrixObject::index_type col) const {
  // This class treats 'U' type matrix.
  // if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j; in Fortran
  // index = row + col * (col +1) /2; // in C/C++ (in row >= col)
  assert((0 <= row) && (row < this->getNumOfRows()));
  assert((0 <= col) && (col < this->getNumOfCols()));

  if (row < col) {
    std::swap(row, col);
  }
  const TlMatrixObject::size_type index = row * (row + 1) / 2 + col;
  assert(index < this->getNumOfElements());

  return index;
}

void TlDenseSymmetricMatrix_ImplLapack::vtr2mat(const std::vector<double>& vtr) {
  assert(vtr.size() == this->getNumOfElements());
  std::copy(vtr.begin(), vtr.end(), this->matrix_);
}


// ---------------------------------------------------------------------------
// private
// ---------------------------------------------------------------------------

// ---------------------------------------------------------------------------
// others
// ---------------------------------------------------------------------------
// DM(G) = DM(S) * DM(G)
TlDenseGeneralMatrix_ImplLapack operator*(
    const TlDenseSymmetricMatrix_ImplLapack& rhs1,
    const TlDenseGeneralMatrix_ImplLapack& rhs2) {
  TlDenseGeneralMatrix_ImplLapack gen_rhs1 = rhs1;

  return gen_rhs1 * rhs2;
}

// DM(G) = DM(G) * DM(S)
TlDenseGeneralMatrix_ImplLapack operator*(
    const TlDenseGeneralMatrix_ImplLapack& rhs1,
    const TlDenseSymmetricMatrix_ImplLapack& rhs2) {
  TlDenseGeneralMatrix_ImplLapack gen_rhs2 = rhs2;

  return rhs1 * gen_rhs2;
}

// DV = DM(S) * DV
TlDenseVector_ImplLapack operator*(const TlDenseSymmetricMatrix_ImplLapack& mat,
                                   const TlDenseVector_ImplLapack& vec) {
  TlLogging& logger = TlLogging::getInstance();
  if (mat.getNumOfCols() != vec.getSize()) {
    logger.critical(TlUtils::format("size mismatch: %d != %d (%d@%s)",
                                    mat.getNumOfCols(), vec.getSize(), __LINE__,
                                    __FILE__));
  }

  TlDenseVector_ImplLapack answer(vec.getSize());

  const char UPLO =
      'U';  // L means the lower triangular part of the symmetric matrix
            // U means the upper triangular part of the symmetric matrix
  const int N = mat.getNumOfRows();
  const double ALPHA = 1.0;  // ALPHA specifies the scalar alpha
  const double* AP = mat.matrix_;
  const double* pX = vec.vector_;
  const int INCX = 1;
  double BETA = 1.0;
  double* Y = answer.vector_;
  const int INCY = 1;

  dspmv_(&UPLO, &N, &ALPHA, AP, pX, &INCX, &BETA, Y, &INCY);

  return answer;
}
