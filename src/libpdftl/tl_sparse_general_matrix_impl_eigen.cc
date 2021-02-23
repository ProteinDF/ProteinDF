#include <Eigen/Sparse>
#include <iostream>
#include <limits>

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif  // HAVE_CONFIG_H

#include "tl_dense_general_matrix_impl_eigen.h"
#include "tl_dense_vector_impl_eigen.h"
#include "tl_matrix_utils.h"
#include "tl_sparse_general_matrix_impl_eigen.h"
#include "tl_sparse_symmetric_matrix_impl_eigen.h"

#ifdef HAVE_VIENNACL
#include "tl_sparse_general_matrix_impl_viennacl.h"
#endif  // HAVE_VIENNACL

const double TlSparseGeneralMatrix_ImplEigen::reference_ = 1.0;
const double TlSparseGeneralMatrix_ImplEigen::epsilon_ = std::numeric_limits<double>::epsilon();

TlSparseGeneralMatrix_ImplEigen::TlSparseGeneralMatrix_ImplEigen(const TlMatrixObject::index_type row,
                                                                 const TlMatrixObject::index_type col)
    : matrix_(row, col){};

TlSparseGeneralMatrix_ImplEigen::TlSparseGeneralMatrix_ImplEigen(const TlSparseGeneralMatrix_ImplEigen& rhs) {
    this->matrix_ = rhs.matrix_;
}

TlSparseGeneralMatrix_ImplEigen::TlSparseGeneralMatrix_ImplEigen(const TlSparseSymmetricMatrix_ImplEigen& rhs) {
    this->matrix_ = rhs.matrix_.selfadjointView<Eigen::Lower>();
}

TlSparseGeneralMatrix_ImplEigen::TlSparseGeneralMatrix_ImplEigen(const MatrixDataType& rhs) {
    this->matrix_ = rhs;
}

TlSparseGeneralMatrix_ImplEigen::TlSparseGeneralMatrix_ImplEigen(const TlDenseGeneralMatrix_ImplEigen& rhs)
    : matrix_(rhs.matrix_.sparseView(TlSparseGeneralMatrix_ImplEigen::reference_,
                                     TlSparseGeneralMatrix_ImplEigen::epsilon_)) {
}

#ifdef HAVE_VIENNACL
TlSparseGeneralMatrix_ImplEigen::TlSparseGeneralMatrix_ImplEigen(const TlSparseGeneralMatrix_ImplViennaCL& rhs)
    : matrix_(rhs.getNumOfRows(), rhs.getNumOfCols()) {
    // std::cout << "TlSparseGeneralMatrix_ImplEigen::TlSparseGeneralMatrix_"
    //              "ImplEigen(const TlSparseGeneralMatrix_ImplViennaCL& rhs)"
    //           << std::endl;
    viennacl::copy(rhs.matrix_, this->matrix_);
}
#endif  // HAVE_VIENNACL

TlSparseGeneralMatrix_ImplEigen::~TlSparseGeneralMatrix_ImplEigen() {
}

// ----------------------------------------------------------------------------
// properties
// ----------------------------------------------------------------------------
TlMatrixObject::index_type TlSparseGeneralMatrix_ImplEigen::getNumOfRows() const {
    return this->matrix_.rows();
}

TlMatrixObject::index_type TlSparseGeneralMatrix_ImplEigen::getNumOfCols() const {
    return this->matrix_.cols();
}

void TlSparseGeneralMatrix_ImplEigen::resize(const TlMatrixObject::index_type newRow,
                                             const TlMatrixObject::index_type newCol) {
    this->matrix_.resize(newRow, newCol);
}

double TlSparseGeneralMatrix_ImplEigen::get(const TlMatrixObject::index_type row,
                                            const TlMatrixObject::index_type col) const {
    return this->matrix_.coeff(row, col);
}

void TlSparseGeneralMatrix_ImplEigen::set(const TlMatrixObject::index_type row, const TlMatrixObject::index_type col,
                                          const double value) {
    this->matrix_.coeffRef(row, col) = value;
}

void TlSparseGeneralMatrix_ImplEigen::add(const TlMatrixObject::index_type row, const TlMatrixObject::index_type col,
                                          const double value) {
    this->matrix_.coeffRef(row, col) += value;
}

void TlSparseGeneralMatrix_ImplEigen::mul(const TlMatrixObject::index_type row, const TlMatrixObject::index_type col,
                                          const double value) {
    this->matrix_.coeffRef(row, col) *= value;
}

// ----------------------------------------------------------------------------
// operators
// ----------------------------------------------------------------------------
TlSparseGeneralMatrix_ImplEigen& TlSparseGeneralMatrix_ImplEigen::operator=(
    const TlSparseGeneralMatrix_ImplEigen& rhs) {
    if (this != &rhs) {
        this->matrix_ = rhs.matrix_;
    }

    return *this;
}

TlSparseGeneralMatrix_ImplEigen& TlSparseGeneralMatrix_ImplEigen::operator+=(
    const TlSparseGeneralMatrix_ImplEigen& rhs) {
    this->matrix_ += rhs.matrix_;
    return *this;
}

TlSparseGeneralMatrix_ImplEigen& TlSparseGeneralMatrix_ImplEigen::operator-=(
    const TlSparseGeneralMatrix_ImplEigen& rhs) {
    this->matrix_ -= rhs.matrix_;
    return *this;
}

TlSparseGeneralMatrix_ImplEigen& TlSparseGeneralMatrix_ImplEigen::operator*=(const double coef) {
    this->matrix_ *= coef;
    return *this;
}

// ----------------------------------------------------------------------------
// I/O
// ----------------------------------------------------------------------------
bool TlSparseGeneralMatrix_ImplEigen::load(const std::string& path) {
    bool answer = false;
    TlMatrixObject::HeaderInfo headerInfo;
    // TlMatrixObject::MatrixType matrixType;

    const TlMatrixUtils::FileSize headerSize = TlMatrixUtils::getHeaderInfo(path, &headerInfo);
    if (headerSize > 0) {
        this->resize(headerInfo.numOfRows, headerInfo.numOfCols);
        this->matrix_.setZero();

        std::fstream fs;
        fs.open(path.c_str(), std::ios::in | std::ios::binary);
        if (!fs.fail()) {
            fs.seekg(headerSize);
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

bool TlSparseGeneralMatrix_ImplEigen::save(const std::string& path) const {
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
TlDenseGeneralMatrix_ImplEigen operator*(const TlDenseGeneralMatrix_ImplEigen& mat1,
                                         const TlSparseGeneralMatrix_ImplEigen& mat2) {
    assert(mat1.getNumOfCols() == mat2.getNumOfRows());
    TlDenseGeneralMatrix_ImplEigen answer;
    answer.matrix_ = mat1.matrix_ * mat2.matrix_;

    return answer;
}

// DM(G) = SM(G) * DM(G)
TlDenseGeneralMatrix_ImplEigen operator*(const TlSparseGeneralMatrix_ImplEigen& mat1,
                                         const TlDenseGeneralMatrix_ImplEigen& mat2) {
    assert(mat1.getNumOfCols() == mat2.getNumOfRows());
    TlDenseGeneralMatrix_ImplEigen answer;
    answer.matrix_ = mat1.matrix_ * mat2.matrix_;

    return answer;
}

// SM(G) = SM(G) * SM(G)
TlSparseGeneralMatrix_ImplEigen operator*(const TlSparseGeneralMatrix_ImplEigen& mat1,
                                          const TlSparseGeneralMatrix_ImplEigen& mat2) {
    TlSparseGeneralMatrix_ImplEigen answer;
    answer.matrix_ = mat1.matrix_ * mat2.matrix_;

    return answer;
}

// DV = SM(G) * DV
TlDenseVector_ImplEigen operator*(const TlSparseGeneralMatrix_ImplEigen& mat, const TlDenseVector_ImplEigen& vtr) {
    // std::cout << "TlDenseVector_ImplEigen operator*(const
    // TlSparseGeneralMatrix_ImplEigen& mat, const TlDenseVector_ImplEigen& vtr)"
    // << std::endl;
    assert(mat.getNumOfCols() == vtr.getSize());
    TlDenseVector_ImplEigen answer;
    answer.vector_ = mat.matrix_ * vtr.vector_;

    return answer;
}

// DV = DV * SM(G)
TlDenseVector_ImplEigen operator*(const TlDenseVector_ImplEigen& vtr, const TlSparseGeneralMatrix_ImplEigen& mat) {
    // std::cout << "TlDenseVector_ImplEigen operator*(const
    // TlDenseVector_ImplEigen& vtr, const TlSparseGeneralMatrix_ImplEigen& mat)"
    // << std::endl;
    assert(mat.getNumOfRows() == vtr.getSize());
    TlDenseVector_ImplEigen answer;
    answer.vector_ = mat.matrix_.transpose() * vtr.vector_;

    return answer;
}
