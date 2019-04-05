#ifdef HAVE_CONFIG_H
#include "config.h"
#endif  // HAVE_CONFIG_H

#include <Eigen/Core>
#include <Eigen/LU>
#include <iostream>

#include "tl_dense_general_matrix_impl_eigen.h"
#include "tl_dense_symmetric_matrix_impl_eigen.h"
#include "tl_dense_vector_impl_eigen.h"
#include "tl_sparse_general_matrix_impl_eigen.h"

#ifdef HAVE_VIENNACL
#define VIENNACL_HAVE_EIGEN 1
#include <viennacl/matrix.hpp>
#include <viennacl/matrix_proxy.hpp>
#include "tl_dense_general_matrix_impl_viennacl.h"
#endif  //

TlDenseGeneralMatrix_ImplEigen::TlDenseGeneralMatrix_ImplEigen(
    const TlMatrixObject::index_type row, const TlMatrixObject::index_type col)
    : matrix_(MatrixDataType::Zero(row, col)){};

TlDenseGeneralMatrix_ImplEigen::TlDenseGeneralMatrix_ImplEigen(
    const TlDenseGeneralMatrix_ImplEigen& rhs) {
  this->matrix_ = rhs.matrix_;
}

TlDenseGeneralMatrix_ImplEigen::TlDenseGeneralMatrix_ImplEigen(
    const TlDenseSymmetricMatrix_ImplEigen& rhs) {
  this->matrix_ = rhs.matrix_;
}

TlDenseGeneralMatrix_ImplEigen::TlDenseGeneralMatrix_ImplEigen(
    const MatrixDataType& rhs) {
  this->matrix_ = rhs;
}

TlDenseGeneralMatrix_ImplEigen::TlDenseGeneralMatrix_ImplEigen(
    const TlSparseGeneralMatrix_ImplEigen& sm)
    : matrix_(sm.matrix_) {}

#ifdef HAVE_VIENNACL
TlDenseGeneralMatrix_ImplEigen::TlDenseGeneralMatrix_ImplEigen(
    const TlDenseGeneralMatrix_ImplViennaCL& rhs)
    : matrix_(rhs.getNumOfRows(), rhs.getNumOfCols()) {
  viennacl::copy(rhs.matrix_, this->matrix_);
}
#endif  // HAVE_VIENNACL

TlDenseGeneralMatrix_ImplEigen::~TlDenseGeneralMatrix_ImplEigen() {}

void TlDenseGeneralMatrix_ImplEigen::vtr2mat(const std::vector<double>& vtr) {
  assert(vtr.size() == this->getNumOfRows() * this->getNumOfCols());
  this->matrix_ =
      MapTypeConst(&(vtr[0]), this->getNumOfRows(), this->getNumOfCols());
}

// ----------------------
TlMatrixObject::index_type TlDenseGeneralMatrix_ImplEigen::getNumOfRows()
    const {
  return this->matrix_.rows();
}

TlMatrixObject::index_type TlDenseGeneralMatrix_ImplEigen::getNumOfCols()
    const {
  return this->matrix_.cols();
}

void TlDenseGeneralMatrix_ImplEigen::resize(
    const TlMatrixObject::index_type newRow,
    const TlMatrixObject::index_type newCol) {
  this->matrix_.conservativeResizeLike(MatrixDataType::Zero(newRow, newCol));
}

double TlDenseGeneralMatrix_ImplEigen::get(
    const TlMatrixObject::index_type row,
    const TlMatrixObject::index_type col) const {
  return this->matrix_(row, col);
}

void TlDenseGeneralMatrix_ImplEigen::set(const TlMatrixObject::index_type row,
                                         const TlMatrixObject::index_type col,
                                         const double value) {
  this->matrix_(row, col) = value;
}

void TlDenseGeneralMatrix_ImplEigen::add(const TlMatrixObject::index_type row,
                                         const TlMatrixObject::index_type col,
                                         const double value) {
  this->matrix_(row, col) += value;
}

// ---------------------------------------------------------------------------
// operators
// ---------------------------------------------------------------------------
TlDenseGeneralMatrix_ImplEigen& TlDenseGeneralMatrix_ImplEigen::operator=(
    const TlDenseGeneralMatrix_ImplEigen& rhs) {
  this->matrix_ = rhs.matrix_;

  return *this;
}

TlDenseGeneralMatrix_ImplEigen& TlDenseGeneralMatrix_ImplEigen::operator+=(
    const TlDenseGeneralMatrix_ImplEigen& rhs) {
  const TlMatrixObject::index_type row1 = this->getNumOfRows();
  const TlMatrixObject::index_type col1 = this->getNumOfCols();
  const TlMatrixObject::index_type row2 = rhs.getNumOfRows();
  const TlMatrixObject::index_type col2 = rhs.getNumOfCols();
  assert(row1 == row2);
  assert(col1 == col2);

  this->matrix_ += rhs.matrix_;

  return *this;
}

TlDenseGeneralMatrix_ImplEigen& TlDenseGeneralMatrix_ImplEigen::operator-=(
    const TlDenseGeneralMatrix_ImplEigen& rhs) {
  const TlMatrixObject::index_type row1 = this->getNumOfRows();
  const TlMatrixObject::index_type col1 = this->getNumOfCols();
  const TlMatrixObject::index_type row2 = rhs.getNumOfRows();
  const TlMatrixObject::index_type col2 = rhs.getNumOfCols();
  assert(row1 == row2);
  assert(col1 == col2);

  this->matrix_ -= rhs.matrix_;

  return *this;
}

TlDenseGeneralMatrix_ImplEigen& TlDenseGeneralMatrix_ImplEigen::operator*=(
    const double coef) {
  this->matrix_ *= coef;

  return *this;
}

TlDenseGeneralMatrix_ImplEigen& TlDenseGeneralMatrix_ImplEigen::operator/=(
    const double coef) {
  this->matrix_ *= (1.0 / coef);

  return *this;
}

TlDenseGeneralMatrix_ImplEigen& TlDenseGeneralMatrix_ImplEigen::operator*=(
    const TlDenseGeneralMatrix_ImplEigen& rhs) {
  const TlMatrixObject::index_type row1 = this->getNumOfRows();
  const TlMatrixObject::index_type col1 = this->getNumOfCols();
  const TlMatrixObject::index_type row2 = rhs.getNumOfRows();
  const TlMatrixObject::index_type col2 = rhs.getNumOfCols();
  assert(col1 == row2);
  const MatrixDataType tmp = this->matrix_;

  this->resize(row1, col2);
  this->matrix_ = tmp * rhs.matrix_;

  return *this;
}

// -----------------------------------------------------------------------------
// operations
// -----------------------------------------------------------------------------
double TlDenseGeneralMatrix_ImplEigen::sum() const {
  return this->matrix_.array().sum();
}

double TlDenseGeneralMatrix_ImplEigen::getRMS() const {
  const double elements = this->getNumOfRows() * this->getNumOfCols();
  const double sum2 = (this->matrix_.array() * this->matrix_.array()).sum();
  const double rms = std::sqrt(sum2 / elements);

  return rms;
}

double TlDenseGeneralMatrix_ImplEigen::getMaxAbsoluteElement(
    TlMatrixObject::index_type* outRow,
    TlMatrixObject::index_type* outCol) const {
  TlMatrixObject::index_type max_row, max_col;
  TlMatrixObject::index_type min_row, min_col;
  const double max_v = this->matrix_.maxCoeff(&max_row, &max_col);
  const double min_v = this->matrix_.minCoeff(&min_row, &min_col);

  double answer;
  if (std::fabs(min_v) < max_v) {
    answer = max_v;
    if (outRow != NULL) {
      *outRow = max_row;
    }
    if (outCol != NULL) {
      *outCol = max_col;
    }
  } else {
    answer = min_v;
    if (outRow != NULL) {
      *outRow = min_row;
    }
    if (outCol != NULL) {
      *outCol = min_col;
    }
  }

  return answer;
}

void TlDenseGeneralMatrix_ImplEigen::transposeInPlace() {
  TlDenseGeneralMatrix_ImplEigen tmp = this->transpose();
  std::swap(*this, tmp);
}

const TlDenseGeneralMatrix_ImplEigen&
TlDenseGeneralMatrix_ImplEigen::dotInPlace(
    const TlDenseGeneralMatrix_ImplEigen& rhs) {
  this->matrix_.array() *= rhs.matrix_.array();
  return *this;
}

TlDenseGeneralMatrix_ImplEigen TlDenseGeneralMatrix_ImplEigen::transpose()
    const {
  TlDenseGeneralMatrix_ImplEigen answer;
  answer.matrix_ = this->matrix_.transpose();
  return answer;
}

TlDenseGeneralMatrix_ImplEigen TlDenseGeneralMatrix_ImplEigen::inverse() const {
  TlDenseGeneralMatrix_ImplEigen answer;
  answer.matrix_ = this->matrix_.inverse();
  return answer;
}

// ---------------------------------------------------------------------------
// protected
// ---------------------------------------------------------------------------

// ---------------------------------------------------------------------------
// others
// ---------------------------------------------------------------------------
TlDenseVector_ImplEigen operator*(const TlDenseGeneralMatrix_ImplEigen& mat,
                                  const TlDenseVector_ImplEigen& vec) {
  TlDenseVector_ImplEigen answer;
  answer.vector_ = mat.matrix_ * vec.vector_;

  return answer;
}

TlDenseVector_ImplEigen operator*(const TlDenseVector_ImplEigen& vec,
                                  const TlDenseGeneralMatrix_ImplEigen& mat) {
  TlDenseVector_ImplEigen answer;
  answer.vector_ = vec.vector_.transpose() * mat.matrix_;

  return answer;
}

TlDenseGeneralMatrix_ImplEigen operator*(const double coef, const TlDenseGeneralMatrix_ImplEigen& DM) {
  TlDenseGeneralMatrix_ImplEigen answer = DM;
  answer *= coef;
  return answer;
}
