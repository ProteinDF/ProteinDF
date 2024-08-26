#include <Eigen/Sparse>
#include <iostream>
#include <limits>

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif  // HAVE_CONFIG_H

#include "tl_dense_general_matrix_impl_eigen_float.h"
#include "tl_dense_vector_impl_eigen_float.h"
#include "tl_matrix_utils.h"
#include "tl_sparse_general_matrix_impl_eigen_float.h"
#include "tl_sparse_symmetric_matrix_impl_eigen_float.h"

#ifdef HAVE_VIENNACL
#include "tl_sparse_general_matrix_impl_viennacl_float.h"
#endif  // HAVE_VIENNACL

const double TlSparseGeneralMatrix_ImplEigenFloat::reference_ = 1.0;
const double TlSparseGeneralMatrix_ImplEigenFloat::epsilon_ = std::numeric_limits<double>::epsilon();

TlSparseGeneralMatrix_ImplEigenFloat::TlSparseGeneralMatrix_ImplEigenFloat(const TlMatrixObject::index_type row,
                                                                           const TlMatrixObject::index_type col)
    : matrix_(row, col) {};

TlSparseGeneralMatrix_ImplEigenFloat::TlSparseGeneralMatrix_ImplEigenFloat(const TlSparseGeneralMatrix_ImplEigenFloat& rhs) {
    this->matrix_ = rhs.matrix_;
}

TlSparseGeneralMatrix_ImplEigenFloat::TlSparseGeneralMatrix_ImplEigenFloat(const TlSparseSymmetricMatrix_ImplEigenFloat& rhs) {
    this->matrix_ = rhs.matrix_.selfadjointView<Eigen::Lower>();
}

TlSparseGeneralMatrix_ImplEigenFloat::TlSparseGeneralMatrix_ImplEigenFloat(const MatrixDataType& rhs) {
    this->matrix_ = rhs;
}

TlSparseGeneralMatrix_ImplEigenFloat::TlSparseGeneralMatrix_ImplEigenFloat(const TlDenseGeneralMatrix_ImplEigenFloat& rhs)
    : matrix_(rhs.matrix_.sparseView(TlSparseGeneralMatrix_ImplEigenFloat::reference_,
                                     TlSparseGeneralMatrix_ImplEigenFloat::epsilon_)) {
}

#ifdef HAVE_VIENNACL
TlSparseGeneralMatrix_ImplEigenFloat::TlSparseGeneralMatrix_ImplEigenFloat(const TlSparseGeneralMatrix_ImplViennaCLFloat& rhs)
    : matrix_(rhs.getNumOfRows(), rhs.getNumOfCols()) {
    // std::cout << "TlSparseGeneralMatrix_ImplEigenFloat::TlSparseGeneralMatrix_"
    //              "ImplEigen(const TlSparseGeneralMatrix_ImplViennaCLFloat& rhs)"
    //           << std::endl;
    viennacl::copy(rhs.matrix_, this->matrix_);
}
#endif  // HAVE_VIENNACL

TlSparseGeneralMatrix_ImplEigenFloat::~TlSparseGeneralMatrix_ImplEigenFloat() {
}

// ----------------------------------------------------------------------------
// properties
// ----------------------------------------------------------------------------
TlMatrixObject::index_type TlSparseGeneralMatrix_ImplEigenFloat::getNumOfRows() const {
    return this->matrix_.rows();
}

TlMatrixObject::index_type TlSparseGeneralMatrix_ImplEigenFloat::getNumOfCols() const {
    return this->matrix_.cols();
}

void TlSparseGeneralMatrix_ImplEigenFloat::resize(const TlMatrixObject::index_type newRow,
                                                  const TlMatrixObject::index_type newCol) {
    this->matrix_.resize(newRow, newCol);
}

double TlSparseGeneralMatrix_ImplEigenFloat::get(const TlMatrixObject::index_type row,
                                                 const TlMatrixObject::index_type col) const {
    return this->matrix_.coeff(row, col);
}

void TlSparseGeneralMatrix_ImplEigenFloat::set(const TlMatrixObject::index_type row, const TlMatrixObject::index_type col,
                                               const double value) {
    this->matrix_.coeffRef(row, col) = value;
}

void TlSparseGeneralMatrix_ImplEigenFloat::add(const TlMatrixObject::index_type row, const TlMatrixObject::index_type col,
                                               const double value) {
    this->matrix_.coeffRef(row, col) += value;
}

void TlSparseGeneralMatrix_ImplEigenFloat::mul(const TlMatrixObject::index_type row, const TlMatrixObject::index_type col,
                                               const double value) {
    this->matrix_.coeffRef(row, col) *= value;
}

// ----------------------------------------------------------------------------
// operators
// ----------------------------------------------------------------------------
TlSparseGeneralMatrix_ImplEigenFloat& TlSparseGeneralMatrix_ImplEigenFloat::operator=(
    const TlSparseGeneralMatrix_ImplEigenFloat& rhs) {
    if (this != &rhs) {
        this->matrix_ = rhs.matrix_;
    }

    return *this;
}

TlSparseGeneralMatrix_ImplEigenFloat& TlSparseGeneralMatrix_ImplEigenFloat::operator+=(
    const TlSparseGeneralMatrix_ImplEigenFloat& rhs) {
    this->matrix_ += rhs.matrix_;
    return *this;
}

TlSparseGeneralMatrix_ImplEigenFloat& TlSparseGeneralMatrix_ImplEigenFloat::operator-=(
    const TlSparseGeneralMatrix_ImplEigenFloat& rhs) {
    this->matrix_ -= rhs.matrix_;
    return *this;
}

TlSparseGeneralMatrix_ImplEigenFloat& TlSparseGeneralMatrix_ImplEigenFloat::operator*=(const double coef) {
    this->matrix_ *= coef;
    return *this;
}

// ----------------------------------------------------------------------------
// I/O
// ----------------------------------------------------------------------------
bool TlSparseGeneralMatrix_ImplEigenFloat::load(const std::string& path) {
    bool answer = false;
    TlMatrixObject::HeaderInfo headerInfo;

    const bool isLoadable = TlMatrixUtils::getHeaderInfo(path, &headerInfo);
    if (isLoadable == true) {
        this->resize(headerInfo.numOfRows, headerInfo.numOfCols);
        this->matrix_.setZero();

        std::fstream fs;
        fs.open(path.c_str(), std::ios::in | std::ios::binary);
        if (!fs.fail()) {
            fs.seekg(headerInfo.headerSize);
            const std::size_t numOfItems = headerInfo.numOfItems;
            for (std::size_t i = 0; i < numOfItems; ++i) {
                TlMatrixObject::index_type row;
                TlMatrixObject::index_type col;
                double value;
                fs.read(reinterpret_cast<char*>(&row), sizeof(TlMatrixObject::index_type));
                fs.read(reinterpret_cast<char*>(&col), sizeof(TlMatrixObject::index_type));
                fs.read(reinterpret_cast<char*>(&value), sizeof(double));
                this->set(row, col, value);
            }
        }
        fs.close();
        answer = true;
    }

    return answer;
}

bool TlSparseGeneralMatrix_ImplEigenFloat::save(const std::string& path) const {
    bool answer = false;

    std::fstream fs;
    fs.open(path.c_str(), std::ios::out | std::ios::binary);
    if (!fs.fail()) {
        // header
        char matrixType = static_cast<char>(TlMatrixObject::COOF);
        const TlMatrixObject::index_type numOfRows = this->getNumOfRows();
        const TlMatrixObject::index_type numOfCols = this->getNumOfCols();
        const std::size_t numOfItems = this->matrix_.nonZeros();
        fs.write(&matrixType, sizeof(char));
        fs.write(reinterpret_cast<const char*>(&numOfRows), sizeof(TlMatrixObject::index_type));
        fs.write(reinterpret_cast<const char*>(&numOfCols), sizeof(TlMatrixObject::index_type));
        fs.write(reinterpret_cast<const char*>(numOfItems), sizeof(std::size_t));

        // contents
        for (int k = 0; k < this->matrix_.outerSize(); ++k) {
            for (MatrixDataType::InnerIterator it(this->matrix_, k); it; ++it) {
                const TlMatrixObject::index_type row = it.row();  // row index
                const TlMatrixObject::index_type col = it.col();  // col index (here it is equal to k)
                const double value = it.value();
                fs.write(reinterpret_cast<const char*>(&row), sizeof(TlMatrixObject::index_type));
                fs.write(reinterpret_cast<const char*>(&col), sizeof(TlMatrixObject::index_type));
                fs.write(reinterpret_cast<const char*>(&value), sizeof(double));
            }
        }
        answer = true;
    }
    fs.close();

    return answer;
}

// ----------------------------------------------------------------------------
// others
// ----------------------------------------------------------------------------
// DM(G) = DM(G) * SM(G)
TlDenseGeneralMatrix_ImplEigenFloat operator*(const TlDenseGeneralMatrix_ImplEigenFloat& mat1,
                                         const TlSparseGeneralMatrix_ImplEigenFloat& mat2) {
    assert(mat1.getNumOfCols() == mat2.getNumOfRows());
    TlDenseGeneralMatrix_ImplEigenFloat answer;
    answer.matrix_ = mat1.matrix_ * mat2.matrix_;

    return answer;
}

// DM(G) = SM(G) * DM(G)
TlDenseGeneralMatrix_ImplEigenFloat operator*(const TlSparseGeneralMatrix_ImplEigenFloat& mat1,
                                         const TlDenseGeneralMatrix_ImplEigenFloat& mat2) {
    assert(mat1.getNumOfCols() == mat2.getNumOfRows());
    TlDenseGeneralMatrix_ImplEigenFloat answer;
    answer.matrix_ = mat1.matrix_ * mat2.matrix_;

    return answer;
}

// SM(G) = SM(G) * SM(G)
TlSparseGeneralMatrix_ImplEigenFloat operator*(const TlSparseGeneralMatrix_ImplEigenFloat& mat1,
                                               const TlSparseGeneralMatrix_ImplEigenFloat& mat2) {
    TlSparseGeneralMatrix_ImplEigenFloat answer;
    answer.matrix_ = mat1.matrix_ * mat2.matrix_;

    return answer;
}

// DV = SM(G) * DV
TlDenseVector_ImplEigenFloat operator*(const TlSparseGeneralMatrix_ImplEigenFloat& mat, const TlDenseVector_ImplEigenFloat& vtr) {
    // std::cout << "TlDenseVector_ImplEigenFloat operator*(const
    // TlSparseGeneralMatrix_ImplEigenFloat& mat, const TlDenseVector_ImplEigenFloat& vtr)"
    // << std::endl;
    assert(mat.getNumOfCols() == vtr.getSize());
    TlDenseVector_ImplEigenFloat answer;
    answer.vector_ = mat.matrix_ * vtr.vector_;

    return answer;
}

// DV = DV * SM(G)
TlDenseVector_ImplEigenFloat operator*(const TlDenseVector_ImplEigenFloat& vtr, const TlSparseGeneralMatrix_ImplEigenFloat& mat) {
    // std::cout << "TlDenseVector_ImplEigenFloat operator*(const
    // TlDenseVector_ImplEigenFloat& vtr, const TlSparseGeneralMatrix_ImplEigenFloat& mat)"
    // << std::endl;
    assert(mat.getNumOfRows() == vtr.getSize());
    TlDenseVector_ImplEigenFloat answer;
    answer.vector_ = mat.matrix_.transpose() * vtr.vector_;

    return answer;
}
