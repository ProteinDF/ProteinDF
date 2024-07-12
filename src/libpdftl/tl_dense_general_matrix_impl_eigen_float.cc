#ifdef HAVE_CONFIG_H
#include "config.h"
#endif  // HAVE_CONFIG_H

#include <Eigen/Core>
#include <Eigen/LU>
#include <iostream>

#include "tl_dense_general_matrix_impl_eigen_float.h"
#include "tl_dense_symmetric_matrix_impl_eigen_float.h"
#include "tl_dense_vector_impl_eigen_float.h"
#include "tl_sparse_general_matrix_impl_eigen_float.h"

#ifdef HAVE_VIENNACL
#define VIENNACL_HAVE_EIGEN 1
#include <viennacl/matrix.hpp>
#include <viennacl/matrix_proxy.hpp>

#include "tl_dense_general_matrix_impl_viennacl.h"
#endif  //

// ----------------------------------------------------------------------------
// constructor & destructor
// ----------------------------------------------------------------------------
TlDenseGeneralMatrix_ImplEigenFloat::TlDenseGeneralMatrix_ImplEigenFloat(
    const TlMatrixObject::index_type row, const TlMatrixObject::index_type col,
    double const* const pBuf)
    : matrix_(MatrixDataType::Zero(row, col)) {
    if (pBuf != NULL) {
        this->vtr2mat(pBuf);
    }
};

TlDenseGeneralMatrix_ImplEigenFloat::TlDenseGeneralMatrix_ImplEigenFloat(
    const TlDenseGeneralMatrix_ImplEigenFloat& rhs) {
    this->matrix_ = rhs.matrix_;
}

TlDenseGeneralMatrix_ImplEigenFloat::TlDenseGeneralMatrix_ImplEigenFloat(
    const TlDenseSymmetricMatrix_ImplEigenFloat& rhs) {
    this->matrix_ = rhs.matrix_;
}

TlDenseGeneralMatrix_ImplEigenFloat::TlDenseGeneralMatrix_ImplEigenFloat(
    const MatrixDataType& rhs) {
    this->matrix_ = rhs;
}

TlDenseGeneralMatrix_ImplEigenFloat::TlDenseGeneralMatrix_ImplEigenFloat(
    const TlSparseGeneralMatrix_ImplEigenFloat& sm)
    : matrix_(sm.matrix_) {}

#ifdef HAVE_VIENNACL
TlDenseGeneralMatrix_ImplEigenFloat::TlDenseGeneralMatrix_ImplEigenFloat(
    const TlDenseGeneralMatrix_ImplViennaCL& rhs)
    : matrix_(rhs.getNumOfRows(), rhs.getNumOfCols()) {
    viennacl::copy(rhs.matrix_, this->matrix_);
}
#endif  // HAVE_VIENNACL

TlDenseGeneralMatrix_ImplEigenFloat::~TlDenseGeneralMatrix_ImplEigenFloat() {}

TlDenseGeneralMatrix_ImplEigenFloat::operator std::vector<double>() const {
    const std::size_t row = this->getNumOfRows();
    const std::size_t col = this->getNumOfCols();
    std::vector<float> tmp(row * col);
    MapType(&(tmp[0]), row, col) = this->matrix_;

    return std::vector<double>(tmp.begin(), tmp.end());
}

// ---------------------------------------------------------------------------
// properties
// ---------------------------------------------------------------------------
TlMatrixObject::index_type TlDenseGeneralMatrix_ImplEigenFloat::getNumOfRows()
    const {
    return this->matrix_.rows();
}

TlMatrixObject::index_type TlDenseGeneralMatrix_ImplEigenFloat::getNumOfCols()
    const {
    return this->matrix_.cols();
}

void TlDenseGeneralMatrix_ImplEigenFloat::resize(
    const TlMatrixObject::index_type newRow,
    const TlMatrixObject::index_type newCol) {
    this->matrix_.conservativeResizeLike(MatrixDataType::Zero(newRow, newCol));
}

double TlDenseGeneralMatrix_ImplEigenFloat::get(
    const TlMatrixObject::index_type row,
    const TlMatrixObject::index_type col) const {
    return this->matrix_(row, col);
}

void TlDenseGeneralMatrix_ImplEigenFloat::set(const TlMatrixObject::index_type row,
                                              const TlMatrixObject::index_type col,
                                              const double value) {
    this->matrix_(row, col) = value;
}

void TlDenseGeneralMatrix_ImplEigenFloat::add(const TlMatrixObject::index_type row,
                                              const TlMatrixObject::index_type col,
                                              const double value) {
    this->matrix_(row, col) += value;
}

// ---------------------------------------------------------------------------
// operators
// ---------------------------------------------------------------------------
TlDenseGeneralMatrix_ImplEigenFloat& TlDenseGeneralMatrix_ImplEigenFloat::operator=(
    const TlDenseGeneralMatrix_ImplEigenFloat& rhs) {
    this->matrix_ = rhs.matrix_;

    return *this;
}

TlDenseGeneralMatrix_ImplEigenFloat& TlDenseGeneralMatrix_ImplEigenFloat::operator+=(
    const TlDenseGeneralMatrix_ImplEigenFloat& rhs) {
    assert(this->getNumOfRows() == rhs.getNumOfRows());
    assert(this->getNumOfCols() == rhs.getNumOfCols());

    this->matrix_ += rhs.matrix_;

    return *this;
}

TlDenseGeneralMatrix_ImplEigenFloat& TlDenseGeneralMatrix_ImplEigenFloat::operator-=(
    const TlDenseGeneralMatrix_ImplEigenFloat& rhs) {
    assert(this->getNumOfRows() == rhs.getNumOfRows());
    assert(this->getNumOfCols() == rhs.getNumOfCols());

    this->matrix_ -= rhs.matrix_;

    return *this;
}

TlDenseGeneralMatrix_ImplEigenFloat& TlDenseGeneralMatrix_ImplEigenFloat::operator*=(
    const double coef) {
    this->matrix_ *= coef;

    return *this;
}

TlDenseGeneralMatrix_ImplEigenFloat& TlDenseGeneralMatrix_ImplEigenFloat::operator/=(
    const double coef) {
    this->matrix_ *= (1.0 / coef);

    return *this;
}

TlDenseGeneralMatrix_ImplEigenFloat& TlDenseGeneralMatrix_ImplEigenFloat::operator*=(
    const TlDenseGeneralMatrix_ImplEigenFloat& rhs) {
    assert(this->getNumOfCols() == rhs.getNumOfRows());
    const MatrixDataType tmp = this->matrix_;

    this->resize(this->getNumOfRows(), rhs.getNumOfCols());
    this->matrix_ = tmp * rhs.matrix_;

    return *this;
}

// -----------------------------------------------------------------------------
// operations
// -----------------------------------------------------------------------------
double TlDenseGeneralMatrix_ImplEigenFloat::sum() const {
    return this->matrix_.array().sum();
}

double TlDenseGeneralMatrix_ImplEigenFloat::getRMS() const {
    const double elements = this->getNumOfRows() * this->getNumOfCols();
    const double sum2 = (this->matrix_.array() * this->matrix_.array()).sum();
    const double rms = std::sqrt(sum2 / elements);

    return rms;
}

double TlDenseGeneralMatrix_ImplEigenFloat::getMaxAbsoluteElement(
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

void TlDenseGeneralMatrix_ImplEigenFloat::transposeInPlace() {
    TlDenseGeneralMatrix_ImplEigenFloat tmp = this->transpose();
    std::swap(*this, tmp);
}

TlDenseGeneralMatrix_ImplEigenFloat TlDenseGeneralMatrix_ImplEigenFloat::dot(
    const TlDenseGeneralMatrix_ImplEigenFloat& rhs) const {
    TlDenseGeneralMatrix_ImplEigenFloat answer = *this;
    answer.dotInPlace(rhs);

    return answer;
}

const TlDenseGeneralMatrix_ImplEigenFloat&
TlDenseGeneralMatrix_ImplEigenFloat::dotInPlace(
    const TlDenseGeneralMatrix_ImplEigenFloat& rhs) {
    this->matrix_.array() *= rhs.matrix_.array();
    return *this;
}

TlDenseGeneralMatrix_ImplEigenFloat TlDenseGeneralMatrix_ImplEigenFloat::transpose()
    const {
    TlDenseGeneralMatrix_ImplEigenFloat answer;
    answer.matrix_ = this->matrix_.transpose();
    return answer;
}

TlDenseGeneralMatrix_ImplEigenFloat TlDenseGeneralMatrix_ImplEigenFloat::inverse() const {
    TlDenseGeneralMatrix_ImplEigenFloat answer;
    answer.matrix_ = this->matrix_.inverse();
    return answer;
}

// ---------------------------------------------------------------------------
// protected
// ---------------------------------------------------------------------------
void TlDenseGeneralMatrix_ImplEigenFloat::vtr2mat(double const* const pBuf) {
    const std::size_t row = this->getNumOfRows();
    const std::size_t col = this->getNumOfCols();

    std::vector<float> tmp(pBuf, pBuf + (row * col));
    this->matrix_ = MapTypeConst(tmp.data(), row, col);
}

TlMatrixObject::size_type TlDenseGeneralMatrix_ImplEigenFloat::getNumOfElements() const {
    return this->matrix_.size();
}

float* TlDenseGeneralMatrix_ImplEigenFloat::data() {
    return this->matrix_.data();
}

const float* TlDenseGeneralMatrix_ImplEigenFloat::data() const {
    return this->matrix_.data();
}

// ---------------------------------------------------------------------------
// others
// ---------------------------------------------------------------------------
TlDenseVector_ImplEigenFloat operator*(const TlDenseGeneralMatrix_ImplEigenFloat& mat,
                                       const TlDenseVector_ImplEigenFloat& vec) {
    TlDenseVector_ImplEigenFloat answer;
    answer.vector_ = mat.matrix_ * vec.vector_;

    return answer;
}

TlDenseVector_ImplEigenFloat operator*(const TlDenseVector_ImplEigenFloat& vec,
                                       const TlDenseGeneralMatrix_ImplEigenFloat& mat) {
    TlDenseVector_ImplEigenFloat answer;
    answer.vector_ = vec.vector_.transpose() * mat.matrix_;

    return answer;
}

TlDenseGeneralMatrix_ImplEigenFloat operator*(
    const double coef, const TlDenseGeneralMatrix_ImplEigenFloat& DM) {
    TlDenseGeneralMatrix_ImplEigenFloat answer = DM;
    answer *= coef;
    return answer;
}

TlDenseGeneralMatrix_ImplEigenFloat operator*(
    const TlDenseGeneralMatrix_ImplEigenFloat& DM, const double coef) {
    TlDenseGeneralMatrix_ImplEigenFloat answer = DM;
    answer *= coef;
    return answer;
}
