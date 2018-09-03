#include <Eigen/Sparse>

#include "tl_dense_general_matrix_impl_eigen.h"
#include "tl_dense_symmetric_matrix_impl_eigen.h"
#include "tl_dense_vector_impl_eigen.h"
#include "tl_sparse_symmetric_matrix_impl_eigen.h"

#ifdef HAVE_VIENNACL
#include "tl_sparse_symmetric_matrix_impl_viennacl.h"
#endif  // HAVE_VIENNACL

TlSparseSymmetricMatrix_ImplEigen::TlSparseSymmetricMatrix_ImplEigen(
    const TlMatrixObject::index_type dim)
    : TlSparseGeneralMatrix_ImplEigen(dim, dim) {}

TlSparseSymmetricMatrix_ImplEigen::TlSparseSymmetricMatrix_ImplEigen(
    const TlSparseSymmetricMatrix_ImplEigen& rhs)
    : TlSparseGeneralMatrix_ImplEigen(rhs) {}

TlSparseSymmetricMatrix_ImplEigen::TlSparseSymmetricMatrix_ImplEigen(
    const TlSparseGeneralMatrix_ImplEigen& rhs) {
  assert(rhs.getNumOfRows() == rhs.getNumOfCols());
  this->resize(rhs.getNumOfRows(), rhs.getNumOfRows());
  this->matrix_ = rhs.matrix_.selfadjointView<Eigen::Upper>();
  // std::cout << this->matrix_ << std::endl;
}

TlSparseSymmetricMatrix_ImplEigen::TlSparseSymmetricMatrix_ImplEigen(
    const TlDenseSymmetricMatrix_ImplEigen& rhs) {
  this->matrix_ =
      rhs.matrix_.sparseView(TlSparseGeneralMatrix_ImplEigen::reference_,
                             TlSparseGeneralMatrix_ImplEigen::epsilon_);
}

#ifdef HAVE_VIENNACL
TlSparseSymmetricMatrix_ImplEigen::TlSparseSymmetricMatrix_ImplEigen(
    const TlSparseSymmetricMatrix_ImplViennaCL& rhs) : TlSparseGeneralMatrix_ImplEigen(rhs.getNumOfRows(), rhs.getNumOfCols()) {
  viennacl::copy(rhs.matrix_, this->matrix_);
}
#endif  // HAVE_VIENNACL

TlSparseSymmetricMatrix_ImplEigen::~TlSparseSymmetricMatrix_ImplEigen() {}

double TlSparseSymmetricMatrix_ImplEigen::get(
    const TlMatrixObject::index_type row,
    const TlMatrixObject::index_type col) const {
  double answer;
  // access elements in top-right angular matrix
  if (row >= col) {
    answer = this->matrix_.coeffRef(row, col);
  } else {
    answer = this->matrix_.coeffRef(col, row);
  }

  return answer;
}

void TlSparseSymmetricMatrix_ImplEigen::set(
    const TlMatrixObject::index_type row, const TlMatrixObject::index_type col,
    const double value) {
  // set elements in top-right angular matrix
  if (row >= col) {
    this->matrix_.coeffRef(row, col) = value;
  } else {
    this->matrix_.coeffRef(col, row) = value;
  }
}

void TlSparseSymmetricMatrix_ImplEigen::add(
    const TlMatrixObject::index_type row, const TlMatrixObject::index_type col,
    const double value) {
  // set elements in top-right angular matrix
  if (row >= col) {
    this->matrix_.coeffRef(row, col) += value;
  } else {
    this->matrix_.coeffRef(col, row) += value;
  }
}

void TlSparseSymmetricMatrix_ImplEigen::mul(
    const TlMatrixObject::index_type row, const TlMatrixObject::index_type col,
    const double value) {
  // set elements in top-right angular matrix
  if (row >= col) {
    this->matrix_.coeffRef(row, col) *= value;
  } else {
    this->matrix_.coeffRef(col, row) *= value;
  }
}

// ----------------------------------------------------------------------------
// operators
// ----------------------------------------------------------------------------
TlSparseSymmetricMatrix_ImplEigen& TlSparseSymmetricMatrix_ImplEigen::
operator+=(const TlSparseSymmetricMatrix_ImplEigen& sm) {
  this->matrix_ += sm.matrix_;
  return *this;
}

TlSparseSymmetricMatrix_ImplEigen& TlSparseSymmetricMatrix_ImplEigen::
operator-=(const TlSparseSymmetricMatrix_ImplEigen& sm) {
  this->matrix_ -= sm.matrix_;
  return *this;
}

// ----------------------------------------------------------------------------
// others
// ----------------------------------------------------------------------------
// DM(G) = DM(G) * SM(S)
TlDenseGeneralMatrix_ImplEigen operator*(
    const TlDenseGeneralMatrix_ImplEigen& mat1,
    const TlSparseSymmetricMatrix_ImplEigen& mat2) {
  assert(mat1.getNumOfCols() == mat2.getNumOfRows());
  TlDenseGeneralMatrix_ImplEigen answer;
  answer.matrix_ = mat1.matrix_ * mat2.matrix_.selfadjointView<Eigen::Lower>();

  return answer;
}

// DM(G) = SM(S) * DM(G)
TlDenseGeneralMatrix_ImplEigen operator*(
    const TlSparseSymmetricMatrix_ImplEigen& mat1,
    const TlDenseGeneralMatrix_ImplEigen& mat2) {
  assert(mat1.getNumOfCols() == mat2.getNumOfRows());
  TlDenseGeneralMatrix_ImplEigen answer;
  answer.matrix_ = mat1.matrix_.selfadjointView<Eigen::Lower>() * mat2.matrix_;

  return answer;
}

// SM(G) = SM(G) * SM(S)
TlSparseGeneralMatrix_ImplEigen operator*(
    const TlSparseGeneralMatrix_ImplEigen& sm1,
    const TlSparseSymmetricMatrix_ImplEigen& sm2) {
  assert(sm1.getNumOfCols() == sm2.getNumOfRows());
  TlSparseGeneralMatrix_ImplEigen answer;
  answer.matrix_ =
      sm1.matrix_ * TlSparseGeneralMatrix_ImplEigen::MatrixDataType(
                        sm2.matrix_.selfadjointView<Eigen::Lower>());

  return answer;
}

// SM(G) = SM(S) * SM(G)
TlSparseGeneralMatrix_ImplEigen operator*(
    const TlSparseSymmetricMatrix_ImplEigen& sm1,
    const TlSparseGeneralMatrix_ImplEigen& sm2) {
  assert(sm1.getNumOfCols() == sm2.getNumOfRows());
  TlSparseGeneralMatrix_ImplEigen answer;
  answer.matrix_ = TlSparseGeneralMatrix_ImplEigen::MatrixDataType(
                       sm1.matrix_.selfadjointView<Eigen::Lower>()) *
                   sm2.matrix_;

  return answer;
}

// SM(G) = SM(S) * SM(S)
TlSparseGeneralMatrix_ImplEigen operator*(
    const TlSparseSymmetricMatrix_ImplEigen& sm1,
    const TlSparseSymmetricMatrix_ImplEigen& sm2) {
  assert(sm1.getNumOfCols() == sm2.getNumOfRows());
  TlSparseGeneralMatrix_ImplEigen answer;
  answer.matrix_ = TlSparseGeneralMatrix_ImplEigen::MatrixDataType(
                       sm1.matrix_.selfadjointView<Eigen::Lower>()) *
                   TlSparseGeneralMatrix_ImplEigen::MatrixDataType(
                       sm2.matrix_.selfadjointView<Eigen::Lower>());

  return answer;
}

TlDenseVector_ImplEigen operator*(const TlSparseSymmetricMatrix_ImplEigen& mat,
                                  const TlDenseVector_ImplEigen& vtr) {
  assert(mat.getNumOfCols() == vtr.getSize());
  TlDenseVector_ImplEigen answer;
  answer.vector_ = mat.matrix_.selfadjointView<Eigen::Lower>() * vtr.vector_;

  return answer;
}
