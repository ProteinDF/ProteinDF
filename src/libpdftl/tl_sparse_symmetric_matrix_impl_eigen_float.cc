#include "tl_sparse_symmetric_matrix_impl_eigen_float.h"

#include <Eigen/Sparse>

#include "tl_dense_general_matrix_impl_eigen_float.h"
#include "tl_dense_symmetric_matrix_impl_eigen_float.h"
#include "tl_dense_vector_impl_eigen_float.h"

#ifdef HAVE_VIENNACL
#include "tl_sparse_symmetric_matrix_impl_viennacl.h"
#endif  // HAVE_VIENNACL

TlSparseSymmetricMatrix_ImplEigenFloat::TlSparseSymmetricMatrix_ImplEigenFloat(const TlMatrixObject::index_type dim)
    : TlSparseGeneralMatrix_ImplEigenFloat(dim, dim) {}

TlSparseSymmetricMatrix_ImplEigenFloat::TlSparseSymmetricMatrix_ImplEigenFloat(const TlSparseSymmetricMatrix_ImplEigenFloat& rhs)
    : TlSparseGeneralMatrix_ImplEigenFloat(rhs) {}

TlSparseSymmetricMatrix_ImplEigenFloat::TlSparseSymmetricMatrix_ImplEigenFloat(const TlSparseGeneralMatrix_ImplEigenFloat& rhs) {
    assert(rhs.getNumOfRows() == rhs.getNumOfCols());
    this->resize(rhs.getNumOfRows(), rhs.getNumOfRows());
    this->matrix_ = rhs.matrix_.selfadjointView<Eigen::Upper>();
    // std::cout << this->matrix_ << std::endl;
}

TlSparseSymmetricMatrix_ImplEigenFloat::TlSparseSymmetricMatrix_ImplEigenFloat(const TlDenseSymmetricMatrix_ImplEigenFloat& rhs) {
    this->matrix_ =
        rhs.matrix_.sparseView(TlSparseGeneralMatrix_ImplEigenFloat::reference_,
                               TlSparseGeneralMatrix_ImplEigenFloat::epsilon_);
}

#ifdef HAVE_VIENNACL
TlSparseSymmetricMatrix_ImplEigenFloat::TlSparseSymmetricMatrix_ImplEigenFloat(
    const TlSparseSymmetricMatrix_ImplViennaCL& rhs)
    : TlSparseGeneralMatrix_ImplEigenFloat(rhs.getNumOfRows(), rhs.getNumOfCols()) {
    viennacl::copy(rhs.matrix_, this->matrix_);
}
#endif  // HAVE_VIENNACL

TlSparseSymmetricMatrix_ImplEigenFloat::~TlSparseSymmetricMatrix_ImplEigenFloat() {}

double TlSparseSymmetricMatrix_ImplEigenFloat::get(
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

void TlSparseSymmetricMatrix_ImplEigenFloat::set(
    const TlMatrixObject::index_type row, const TlMatrixObject::index_type col,
    const double value) {
    // set elements in top-right angular matrix
    if (row >= col) {
        this->matrix_.coeffRef(row, col) = value;
    } else {
        this->matrix_.coeffRef(col, row) = value;
    }
}

void TlSparseSymmetricMatrix_ImplEigenFloat::add(
    const TlMatrixObject::index_type row, const TlMatrixObject::index_type col,
    const double value) {
    // set elements in top-right angular matrix
    if (row >= col) {
        this->matrix_.coeffRef(row, col) += value;
    } else {
        this->matrix_.coeffRef(col, row) += value;
    }
}

void TlSparseSymmetricMatrix_ImplEigenFloat::mul(
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
TlSparseSymmetricMatrix_ImplEigenFloat& TlSparseSymmetricMatrix_ImplEigenFloat::
operator+=(const TlSparseSymmetricMatrix_ImplEigenFloat& sm) {
    this->matrix_ += sm.matrix_;
    return *this;
}

TlSparseSymmetricMatrix_ImplEigenFloat& TlSparseSymmetricMatrix_ImplEigenFloat::
operator-=(const TlSparseSymmetricMatrix_ImplEigenFloat& sm) {
    this->matrix_ -= sm.matrix_;
    return *this;
}

// ----------------------------------------------------------------------------
// others
// ----------------------------------------------------------------------------
// DM(G) = DM(G) * SM(S)
TlDenseGeneralMatrix_ImplEigenFloat operator*(const TlDenseGeneralMatrix_ImplEigenFloat& mat1,
                                              const TlSparseSymmetricMatrix_ImplEigenFloat& mat2) {
    assert(mat1.getNumOfCols() == mat2.getNumOfRows());
    TlDenseGeneralMatrix_ImplEigenFloat answer;

    // answer.matrix_ = mat1.matrix_ *
    // (mat2.matrix_.selfadjointView<Eigen::Lower>());
    const TlSparseGeneralMatrix_ImplEigenFloat::MatrixDataType tmp = mat2.matrix_.selfadjointView<Eigen::Lower>();
    answer.matrix_ = mat1.matrix_ * tmp;

    return answer;
}

// DM(G) = SM(S) * DM(G)
TlDenseGeneralMatrix_ImplEigenFloat operator*(const TlSparseSymmetricMatrix_ImplEigenFloat& mat1,
                                              const TlDenseGeneralMatrix_ImplEigenFloat& mat2) {
    assert(mat1.getNumOfCols() == mat2.getNumOfRows());
    TlDenseGeneralMatrix_ImplEigenFloat answer;
    answer.matrix_ = mat1.matrix_.selfadjointView<Eigen::Lower>() * mat2.matrix_;

    return answer;
}

// SM(G) = SM(G) * SM(S)
TlSparseGeneralMatrix_ImplEigenFloat operator*(const TlSparseGeneralMatrix_ImplEigenFloat& sm1,
                                               const TlSparseSymmetricMatrix_ImplEigenFloat& sm2) {
    assert(sm1.getNumOfCols() == sm2.getNumOfRows());
    TlSparseGeneralMatrix_ImplEigenFloat answer;
    answer.matrix_ = sm1.matrix_ * TlSparseGeneralMatrix_ImplEigenFloat::MatrixDataType(sm2.matrix_.selfadjointView<Eigen::Lower>());

    return answer;
}

// SM(G) = SM(S) * SM(G)
TlSparseGeneralMatrix_ImplEigenFloat operator*(const TlSparseSymmetricMatrix_ImplEigenFloat& sm1,
                                               const TlSparseGeneralMatrix_ImplEigenFloat& sm2) {
    assert(sm1.getNumOfCols() == sm2.getNumOfRows());
    TlSparseGeneralMatrix_ImplEigenFloat answer;
    answer.matrix_ = TlSparseGeneralMatrix_ImplEigenFloat::MatrixDataType(sm1.matrix_.selfadjointView<Eigen::Lower>()) * sm2.matrix_;

    return answer;
}

// SM(G) = SM(S) * SM(S)
TlSparseGeneralMatrix_ImplEigenFloat operator*(const TlSparseSymmetricMatrix_ImplEigenFloat& sm1,
                                               const TlSparseSymmetricMatrix_ImplEigenFloat& sm2) {
    assert(sm1.getNumOfCols() == sm2.getNumOfRows());
    TlSparseGeneralMatrix_ImplEigenFloat answer;
    answer.matrix_ = TlSparseGeneralMatrix_ImplEigenFloat::MatrixDataType(sm1.matrix_.selfadjointView<Eigen::Lower>()) *
                     TlSparseGeneralMatrix_ImplEigenFloat::MatrixDataType(sm2.matrix_.selfadjointView<Eigen::Lower>());

    return answer;
}

TlDenseVector_ImplEigenFloat operator*(const TlSparseSymmetricMatrix_ImplEigenFloat& mat,
                                  const TlDenseVector_ImplEigenFloat& vtr) {
    assert(mat.getNumOfCols() == vtr.getSize());
    TlDenseVector_ImplEigenFloat answer;
    answer.vector_ = mat.matrix_.selfadjointView<Eigen::Lower>() * vtr.vector_;

    return answer;
}
