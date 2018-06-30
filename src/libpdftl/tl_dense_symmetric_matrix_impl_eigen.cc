#include "tl_dense_symmetric_matrix_impl_eigen.h"
#include <Eigen/Dense>
#include "tl_dense_vector_impl_eigen.h"

TlDenseSymmetricMatrix_ImplEigen::TlDenseSymmetricMatrix_ImplEigen(
    const TlMatrixObject::index_type dim)
    : TlDenseGeneralMatrix_ImplEigen(dim, dim) {}

TlDenseSymmetricMatrix_ImplEigen::TlDenseSymmetricMatrix_ImplEigen(
    const TlDenseSymmetricMatrix_ImplEigen& rhs)
    : TlDenseGeneralMatrix_ImplEigen(rhs) {}

TlDenseSymmetricMatrix_ImplEigen::TlDenseSymmetricMatrix_ImplEigen(
    const TlDenseGeneralMatrix_ImplEigen& rhs) {
  this->matrix_ = rhs.matrix_;
  this->resize(rhs.getNumOfRows(), rhs.getNumOfRows());
  this->matrix_ = this->matrix_.selfadjointView<Eigen::Upper>();
}

TlDenseSymmetricMatrix_ImplEigen::~TlDenseSymmetricMatrix_ImplEigen() {}

void TlDenseSymmetricMatrix_ImplEigen::set(const TlMatrixObject::index_type row,
                                           const TlMatrixObject::index_type col,
                                           const double value) {
  this->matrix_(row, col) = value;
  if (row != col) {
    this->matrix_(col, row) = value;
  }
}

void TlDenseSymmetricMatrix_ImplEigen::add(const TlMatrixObject::index_type row,
                                           const TlMatrixObject::index_type col,
                                           const double value) {
  this->matrix_(row, col) += value;
  if (row != col) {
    this->matrix_(col, row) += value;
  }
}

// ---------------------------------------------------------------------------
// operators
// ---------------------------------------------------------------------------

// ---------------------------------------------------------------------------
// operations
// ---------------------------------------------------------------------------
TlDenseSymmetricMatrix_ImplEigen TlDenseSymmetricMatrix_ImplEigen::transpose()
    const {
  // do nothing
  return *this;
}

TlDenseSymmetricMatrix_ImplEigen TlDenseSymmetricMatrix_ImplEigen::inverse()
    const {
#if __cplusplus >= 201103L
  std::lock_guard<std::mutex> lock(this->matrix_mutex_);
#endif

  this->matrix_ = this->matrix_.selfadjointView<Eigen::Lower>();

  TlDenseSymmetricMatrix_ImplEigen answer;
  answer.matrix_ = this->matrix_.inverse();

  return answer;
}

bool TlDenseSymmetricMatrix_ImplEigen::eig(
    TlDenseVector_ImplEigen* pEigVal,
    TlDenseGeneralMatrix_ImplEigen* pEigVec) const {
#if __cplusplus >= 201103L
  std::lock_guard<std::mutex> lock(this->matrix_mutex_);
#endif
  bool answer = false;
  this->matrix_ = this->matrix_.selfadjointView<Eigen::Lower>();

  Eigen::SelfAdjointEigenSolver<MatrixDataType> es(this->matrix_);
  if (es.info() == Eigen::Success) {
    if (pEigVal != NULL) {
      *pEigVal = TlDenseVector_ImplEigen(es.eigenvalues());
    }
    if (pEigVec != NULL) {
      *pEigVec = TlDenseGeneralMatrix_ImplEigen(es.eigenvectors());
    }

    answer = true;
  }

  return answer;
}

// ---------------------------------------------------------------------------
// private
// ---------------------------------------------------------------------------
