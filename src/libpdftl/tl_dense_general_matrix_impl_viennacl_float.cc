#ifdef HAVE_CONFIG_H
#include "config.h"
#endif  // HAVE_CONFIG_H

#ifdef HAVE_EIGEN
#define VIENNACL_HAVE_EIGEN
#include <Eigen/Core>
#include <Eigen/LU>
#endif  // HAVE_EIGEN

#include <viennacl/linalg/cg.hpp>
#include <viennacl/linalg/direct_solve.hpp>
#include <viennacl/linalg/fft_operations.hpp>
#include <viennacl/linalg/lu.hpp>
#include <viennacl/linalg/sum.hpp>
#include <viennacl/matrix.hpp>
#include <viennacl/matrix_proxy.hpp>

#include "TlUtils.h"
#include "tl_dense_general_matrix_impl_viennacl.h"
#include "tl_dense_general_matrix_impl_viennacl_float.h"
#include "tl_dense_symmetric_matrix_impl_viennacl_float.h"
#include "tl_dense_vector_impl_viennacl_float.h"
#include "tl_sparse_general_matrix_impl_viennacl_float.h"

#ifdef HAVE_EIGEN
#include "tl_dense_general_matrix_impl_eigen_float.h"
#include "tl_sparse_general_matrix_impl_eigen_float.h"
#endif  // HAVE_EIGEN

// ---------------------------------------------------------------------------
// constructor & destructor
// ---------------------------------------------------------------------------
TlDenseGeneralMatrix_ImplViennaCLFloat::TlDenseGeneralMatrix_ImplViennaCLFloat(const TlMatrixObject::index_type row, const TlMatrixObject::index_type col,
                                                                               double const* const pBuf)
    : matrix_(row, col) {
    if (pBuf != NULL) {
        this->vtr2mat(pBuf);
    }
}

TlDenseGeneralMatrix_ImplViennaCLFloat::TlDenseGeneralMatrix_ImplViennaCLFloat(const TlDenseGeneralMatrix_ImplViennaCLFloat& rhs) {
    this->matrix_ = rhs.matrix_;
}

TlDenseGeneralMatrix_ImplViennaCLFloat::TlDenseGeneralMatrix_ImplViennaCLFloat(const TlDenseGeneralMatrix_ImplViennaCL& rhs)
    : matrix_(rhs.getNumOfRows(), rhs.getNumOfCols()) {
#ifdef HAVE_EIGEN
    {
        Eigen::MatrixXd eigenMd(rhs.getNumOfRows(), rhs.getNumOfCols());
        viennacl::copy(rhs.matrix_, eigenMd);
        Eigen::MatrixXf eigenMf = eigenMd.cast<float>();
        viennacl::copy(eigenMf, this->matrix_);
    }
#else
    {
        const TlMatrixObject::index_type row = rhs.getNumOfRows();
        const TlMatrixObject::index_type col = rhs.getNumOfCols();
        for (TlMatrixObject::index_type r = 0; r < row; ++r) {
            for (TlMatrixObject::index_type c = 0; c < col; ++c) {
                this->matrix_(r, c) = static_cast<float>(rhs.matrix_(r, c));
            }
        }
    }
#endif  // HAVE_EIGEN
}

TlDenseGeneralMatrix_ImplViennaCLFloat::TlDenseGeneralMatrix_ImplViennaCLFloat(const TlDenseSymmetricMatrix_ImplViennaCLFloat& rhs) {
    this->matrix_ = rhs.matrix_;
}

TlDenseGeneralMatrix_ImplViennaCLFloat::TlDenseGeneralMatrix_ImplViennaCLFloat(const TlSparseGeneralMatrix_ImplViennaCLFloat& rhs)
    : matrix_(rhs.getNumOfRows(), rhs.getNumOfCols()) {
    TlDenseGeneralMatrix_ImplEigenFloat DM = TlSparseGeneralMatrix_ImplEigenFloat(rhs);
    viennacl::copy(DM.matrix_, this->matrix_);
}

#ifdef HAVE_EIGEN
TlDenseGeneralMatrix_ImplViennaCLFloat::TlDenseGeneralMatrix_ImplViennaCLFloat(const TlDenseGeneralMatrix_ImplEigenFloat& rhs)
    : matrix_(rhs.getNumOfRows(), rhs.getNumOfCols()) {
    viennacl::copy(rhs.matrix_, this->matrix_);
}
#endif  // HAVE_EIGEN

TlDenseGeneralMatrix_ImplViennaCLFloat::~TlDenseGeneralMatrix_ImplViennaCLFloat() {
}

TlDenseGeneralMatrix_ImplViennaCLFloat::operator std::vector<double>() const {
    const std::size_t row = this->getNumOfRows();
    const std::size_t col = this->getNumOfCols();
    std::vector<double> v(row * col);

#ifdef HAVE_EIGEN
    {
        std::vector<float> v_tmp(row * col);
        {
            EigenMatrixDataType tmp(row, col);
            viennacl::copy(this->matrix_, tmp);
            Eigen::Map<EigenMatrixDataType>(&(v_tmp[0]), row, col) = tmp;
        }
        v = std::vector<double>(v_tmp.begin(), v_tmp.end());
    }
#else
    {
        std::size_t i = 0;
        for (std::size_t c = 0; c < col; ++c) {
            for (std::size_t r = 0; r < row; ++r) {
                v[i] = this->get(r, c);
                ++i;
            }
        }
    }
#endif  // HAVE_EIGEN

    return v;
}

// ---------------------------------------------------------------------------
// properties
// ---------------------------------------------------------------------------
TlMatrixObject::index_type TlDenseGeneralMatrix_ImplViennaCLFloat::getNumOfRows() const {
    return this->matrix_.size1();
}

TlMatrixObject::index_type TlDenseGeneralMatrix_ImplViennaCLFloat::getNumOfCols() const {
    return this->matrix_.size2();
}

void TlDenseGeneralMatrix_ImplViennaCLFloat::resize(
    const TlMatrixObject::index_type newRow,
    const TlMatrixObject::index_type newCol) {
    this->matrix_.resize(newRow, newCol, true);
}

double TlDenseGeneralMatrix_ImplViennaCLFloat::get(const TlMatrixObject::index_type row, const TlMatrixObject::index_type col) const {
    return static_cast<double>(this->matrix_(row, col));
}

void TlDenseGeneralMatrix_ImplViennaCLFloat::set(const TlMatrixObject::index_type row, const TlMatrixObject::index_type col,
                                                 const double value) {
    this->matrix_(row, col) = static_cast<float>(value);
}

void TlDenseGeneralMatrix_ImplViennaCLFloat::add(const TlMatrixObject::index_type row, const TlMatrixObject::index_type col,
                                                 const double value) {
    this->matrix_(row, col) += static_cast<float>(value);
}

// ---------------------------------------------------------------------------
// operators
// ---------------------------------------------------------------------------
TlDenseGeneralMatrix_ImplViennaCLFloat& TlDenseGeneralMatrix_ImplViennaCLFloat::operator=(const TlDenseGeneralMatrix_ImplViennaCLFloat& rhs) {
    if (this != &rhs) {
        this->matrix_ = rhs.matrix_;
    }

    return *this;
}

TlDenseGeneralMatrix_ImplViennaCLFloat& TlDenseGeneralMatrix_ImplViennaCLFloat::operator+=(const TlDenseGeneralMatrix_ImplViennaCLFloat& rhs) {
    const TlMatrixObject::index_type row1 = this->getNumOfRows();
    const TlMatrixObject::index_type col1 = this->getNumOfCols();
    const TlMatrixObject::index_type row2 = rhs.getNumOfRows();
    const TlMatrixObject::index_type col2 = rhs.getNumOfCols();
    assert(row1 == row2);
    assert(col1 == col2);

    this->matrix_ += rhs.matrix_;

    return *this;
}

TlDenseGeneralMatrix_ImplViennaCLFloat& TlDenseGeneralMatrix_ImplViennaCLFloat::operator-=(const TlDenseGeneralMatrix_ImplViennaCLFloat& rhs) {
    const TlMatrixObject::index_type row1 = this->getNumOfRows();
    const TlMatrixObject::index_type col1 = this->getNumOfCols();
    const TlMatrixObject::index_type row2 = rhs.getNumOfRows();
    const TlMatrixObject::index_type col2 = rhs.getNumOfCols();
    assert(row1 == row2);
    assert(col1 == col2);

    this->matrix_ -= rhs.matrix_;

    return *this;
}

TlDenseGeneralMatrix_ImplViennaCLFloat& TlDenseGeneralMatrix_ImplViennaCLFloat::operator*=(const float coef) {
    this->matrix_ *= coef;

    return *this;
}

TlDenseGeneralMatrix_ImplViennaCLFloat& TlDenseGeneralMatrix_ImplViennaCLFloat::operator/=(const float coef) {
    this->matrix_ *= (1.0 / coef);

    return *this;
}

TlDenseGeneralMatrix_ImplViennaCLFloat& TlDenseGeneralMatrix_ImplViennaCLFloat::operator*=(const TlDenseGeneralMatrix_ImplViennaCLFloat& rhs) {
    const TlMatrixObject::index_type row1 = this->getNumOfRows();
    const TlMatrixObject::index_type col1 = this->getNumOfCols();
    const TlMatrixObject::index_type row2 = rhs.getNumOfRows();
    const TlMatrixObject::index_type col2 = rhs.getNumOfCols();
    assert(col1 == row2);
    const MatrixDataType tmp = this->matrix_;

    this->resize(row1, col2);
    this->matrix_ = viennacl::linalg::prod(tmp, rhs.matrix_);

    return *this;
}

// -----------------------------------------------------------------------------
// operations
// -----------------------------------------------------------------------------
double TlDenseGeneralMatrix_ImplViennaCLFloat::sum() const {
    const VectorDataType v = viennacl::linalg::row_sum(this->matrix_);
    const float sum = viennacl::linalg::sum(v);

    return static_cast<double>(sum);
}

double TlDenseGeneralMatrix_ImplViennaCLFloat::getRMS() const {
    const float elements = this->getNumOfRows() * this->getNumOfCols();

    const MatrixDataType mat2 =
        viennacl::linalg::element_prod(this->matrix_, this->matrix_);
    const VectorDataType rows = viennacl::linalg::row_sum(mat2);
    const float sum2 = viennacl::linalg::sum(rows);
    const float rms = std::sqrt(sum2 / elements);

    return static_cast<double>(rms);
}

double TlDenseGeneralMatrix_ImplViennaCLFloat::getMaxAbsoluteElement(TlMatrixObject::index_type* outRow,
                                                                     TlMatrixObject::index_type* outCol) const {
    TlMatrixObject::index_type max_row = 0, max_col = 0;
    float answer = 0.0;
    const unsigned int numOfRows = this->getNumOfRows();
    const unsigned int numOfCols = this->getNumOfCols();
    for (unsigned int r = 0; r < numOfRows; ++r) {
        VectorDataType vec_col(numOfCols);
        viennacl::linalg::matrix_row(this->matrix_, r, vec_col);
        const unsigned int col = viennacl::linalg::index_norm_inf(vec_col);
        float value = vec_col[col];
        if (std::fabs(answer) < std::fabs(value)) {
            max_row = r;
            max_col = col;
            answer = value;
        }
    }

    if (outRow != NULL) {
        *outRow = max_row;
    }
    if (outCol != NULL) {
        *outCol = max_col;
    }

    return static_cast<double>(answer);
}

void TlDenseGeneralMatrix_ImplViennaCLFloat::transposeInPlace() {
    this->matrix_ = viennacl::trans(this->matrix_);
}

TlDenseGeneralMatrix_ImplViennaCLFloat& TlDenseGeneralMatrix_ImplViennaCLFloat::dotInPlace(const TlDenseGeneralMatrix_ImplViennaCLFloat& rhs) {
    const MatrixDataType tmp =
        viennacl::linalg::element_prod(this->matrix_, rhs.matrix_);
    this->matrix_ = tmp;
    return *this;
}

TlDenseGeneralMatrix_ImplViennaCLFloat TlDenseGeneralMatrix_ImplViennaCLFloat::transpose() const {
    TlDenseGeneralMatrix_ImplViennaCLFloat answer(this->getNumOfCols(),
                                                  this->getNumOfRows());
    answer.matrix_ = viennacl::trans(this->matrix_);

    return answer;
}

TlDenseGeneralMatrix_ImplViennaCLFloat TlDenseGeneralMatrix_ImplViennaCLFloat::inverse() const {
    // const TlMatrixObject::index_type dim = this->getNumOfRows();
    // const VectorDataType v = viennacl::scalar_vector<float>(dim, 1.0);
    // MatrixDataType E = viennacl::diag(v);

    TlDenseGeneralMatrix_ImplViennaCLFloat answer(this->getNumOfCols(),
                                                  this->getNumOfRows());
    // answer.matrix_ = viennacl::linalg::solve(this->matrix_, E,
    // viennacl::linalg::cg_tag());

    // LU factorization
    // MatrixDataType tmp = this->matrix_;
    // viennacl::linalg::lu_factorize(tmp);
    // viennacl::linalg::lu_substitute(tmp, E);

#ifdef HAVE_EIGEN
    {
        EigenMatrixDataType eigenMat(this->getNumOfRows(), this->getNumOfCols());
        copy(this->matrix_, eigenMat);
        const EigenMatrixDataType eigenInvMat = eigenMat.inverse();
        answer.resize(eigenInvMat.rows(), eigenInvMat.cols());
        copy(eigenInvMat, answer.matrix_);
    }
#endif  // HAVE_EIGEN

    return answer;
}

TlDenseGeneralMatrix_ImplViennaCLFloat& TlDenseGeneralMatrix_ImplViennaCLFloat::reverseColumns() {
    viennacl::slice sr(0, 1, this->getNumOfRows());
    viennacl::slice sc(this->getNumOfCols() - 1, -1, this->getNumOfCols());

    viennacl::matrix_slice<MatrixDataType> s(this->matrix_, sr, sc);
    const MatrixDataType tmp = s;
    this->matrix_ = tmp;

    return *this;
}

// ---------------------------------------------------------------------------
// protected
// ---------------------------------------------------------------------------
void TlDenseGeneralMatrix_ImplViennaCLFloat::vtr2mat(const double* pBuf) {
    const std::size_t row = this->getNumOfRows();
    const std::size_t col = this->getNumOfCols();

#ifdef HAVE_EIGEN
    {
        std::vector<float> buf_fp32(pBuf, pBuf + (row * col));
        const Eigen::MatrixXf tmp = Eigen::Map<const Eigen::MatrixXf>(&(buf_fp32[0]), row, col);
        viennacl::copy(tmp, this->matrix_);
    }
#else
    {
        std::size_t i = 0;
        for (std::size_t c = 0; c < col; ++c) {
            for (std::size_t r = 0; r < row; ++r) {
                this->set(r, c, pBuf[i]);
                ++i;
            }
        }
    }
#endif  // HAVE_EIGEN
}

// ---------------------------------------------------------------------------
// others
// ---------------------------------------------------------------------------
// DV = DM(G) * DV
TlDenseVector_ImplViennaCLFloat operator*(const TlDenseGeneralMatrix_ImplViennaCLFloat& mat,
                                          const TlDenseVector_ImplViennaCLFloat& vec) {
    assert(mat.getNumOfCols() == vec.getSize());
    TlDenseVector_ImplViennaCLFloat answer(mat.getNumOfRows());
    answer.vector_ = viennacl::linalg::prod(mat.matrix_, vec.vector_);

    return answer;
}

// DV = DV * DM(G)
TlDenseVector_ImplViennaCLFloat operator*(const TlDenseVector_ImplViennaCLFloat& vec,
                                          const TlDenseGeneralMatrix_ImplViennaCLFloat& mat) {
    assert(mat.getNumOfRows() == vec.getSize());
    TlDenseVector_ImplViennaCLFloat answer(mat.getNumOfCols());
    answer.vector_ =
        viennacl::linalg::prod(viennacl::trans(mat.matrix_), vec.vector_);

    return answer;
}

// DM(G) = float * DM(G)
TlDenseGeneralMatrix_ImplViennaCLFloat operator*(const float coef, const TlDenseGeneralMatrix_ImplViennaCLFloat& DM) {
    TlDenseGeneralMatrix_ImplViennaCLFloat answer = DM;
    answer *= coef;
    return answer;
}
