#include "tl_dense_symmetric_matrix_impl_eigen_float.h"

#include <Eigen/Dense>
#include <iostream>

#include "tl_dense_vector_impl_eigen_float.h"
#include "tl_sparse_symmetric_matrix_impl_eigen_float.h"

#ifdef HAVE_VIENNACL
#include "tl_dense_symmetric_matrix_impl_viennacl_float.h"
#endif  // HAVE_VIENNACL

TlDenseSymmetricMatrix_ImplEigenFloat::TlDenseSymmetricMatrix_ImplEigenFloat(
    const TlMatrixObject::index_type dim, double const* const pBuf)
    : TlDenseGeneralMatrix_ImplEigenFloat(dim, dim) {
    if (pBuf != NULL) {
        this->vtr2mat(pBuf);
    }
}

TlDenseSymmetricMatrix_ImplEigenFloat::TlDenseSymmetricMatrix_ImplEigenFloat(
    const TlDenseSymmetricMatrix_ImplEigenFloat& rhs)
    : TlDenseGeneralMatrix_ImplEigenFloat(rhs) {}

TlDenseSymmetricMatrix_ImplEigenFloat::TlDenseSymmetricMatrix_ImplEigenFloat(
    const TlDenseGeneralMatrix_ImplEigenFloat& rhs) {
    this->matrix_ = rhs.matrix_;
    this->resize(rhs.getNumOfRows(), rhs.getNumOfRows());
    this->matrix_ = this->matrix_.selfadjointView<Eigen::Upper>();
}

TlDenseSymmetricMatrix_ImplEigenFloat::TlDenseSymmetricMatrix_ImplEigenFloat(
    const TlSparseSymmetricMatrix_ImplEigenFloat& sm)
    : TlDenseGeneralMatrix_ImplEigenFloat(sm.getNumOfRows(), sm.getNumOfCols()) {
    // this->matrix_ = sm.matrix_.selfadjointView<Eigen::Upper>();
    this->matrix_ = sm.matrix_;
}

#ifdef HAVE_VIENNACL
TlDenseSymmetricMatrix_ImplEigenFloat::TlDenseSymmetricMatrix_ImplEigenFloat(const TlDenseSymmetricMatrix_ImplViennaCLFloat& rhs)
    : TlDenseGeneralMatrix_ImplEigenFloat(rhs.getNumOfRows(), rhs.getNumOfCols()) {
    viennacl::copy(rhs.matrix_, this->matrix_);
}
#endif  // HAVE_VIENNACL

TlDenseSymmetricMatrix_ImplEigenFloat::~TlDenseSymmetricMatrix_ImplEigenFloat() {}

TlDenseSymmetricMatrix_ImplEigenFloat::operator std::vector<double>() const {
    const std::size_t dim = this->getNumOfRows();
    const std::size_t size = dim * (dim + 1) / 2;
    std::cout << "TlDenseSymmetricMatrix_ImplEigenFloat::operator std::vector<double>() const: size=" << size << std::endl;
    std::vector<double> buf(size);

    // column-major
    std::size_t i = 0;
    for (std::size_t c = 0; c < dim; ++c) {
        for (std::size_t r = 0; r <= c; ++r) {
            const double v = this->get(r, c);
            assert(i < size);
            buf[i] = v;
            ++i;
        }
    }

    std::cout << "TlDenseSymmetricMatrix_ImplEigenFloat::operator std::vector<double>() return" << std::endl;
    return buf;
}

// ---------------------------------------------------------------------------
// properties
// ---------------------------------------------------------------------------
void TlDenseSymmetricMatrix_ImplEigenFloat::set(const TlMatrixObject::index_type row,
                                                const TlMatrixObject::index_type col,
                                                const double value) {
    this->matrix_(row, col) = value;
    if (row != col) {
        this->matrix_(col, row) = value;
    }
}

void TlDenseSymmetricMatrix_ImplEigenFloat::add(const TlMatrixObject::index_type row,
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
TlDenseSymmetricMatrix_ImplEigenFloat TlDenseSymmetricMatrix_ImplEigenFloat::transpose()
    const {
    // do nothing
    return *this;
}

TlDenseSymmetricMatrix_ImplEigenFloat TlDenseSymmetricMatrix_ImplEigenFloat::inverse()
    const {
#if __cplusplus >= 201103L
    std::lock_guard<std::mutex> lock(this->matrix_mutex_);
#endif

    this->matrix_ = this->matrix_.selfadjointView<Eigen::Lower>();

    TlDenseSymmetricMatrix_ImplEigenFloat answer;
    answer.matrix_ = this->matrix_.inverse();

    return answer;
}

bool TlDenseSymmetricMatrix_ImplEigenFloat::eig(TlDenseVector_ImplEigenFloat* pEigVal, TlDenseGeneralMatrix_ImplEigenFloat* pEigVec) const {
#if __cplusplus >= 201103L
    std::lock_guard<std::mutex> lock(this->matrix_mutex_);
#endif
    bool answer = false;
    this->matrix_ = this->matrix_.selfadjointView<Eigen::Lower>();

    Eigen::SelfAdjointEigenSolver<MatrixDataType> es(this->matrix_);
    if (es.info() == Eigen::Success) {
        if (pEigVal != NULL) {
            *pEigVal = TlDenseVector_ImplEigenFloat(es.eigenvalues());
        }
        if (pEigVec != NULL) {
            *pEigVec = TlDenseGeneralMatrix_ImplEigenFloat(es.eigenvectors());
        }

        answer = true;
    }

    return answer;
}

// ---------------------------------------------------------------------------
// I/O
// ---------------------------------------------------------------------------
void TlDenseSymmetricMatrix_ImplEigenFloat::vtr2mat(double const* const pBuf) {
    const TlMatrixObject::index_type dim = this->getNumOfRows();

    std::size_t i = 0;
    // column-major
    for (TlMatrixObject::index_type c = 0; c < dim; ++c) {
        for (TlMatrixObject::index_type r = 0; r <= c; ++r) {
            double v = pBuf[i];
            this->set(r, c, v);
            ++i;
        }
    }
}

TlMatrixObject::size_type TlDenseSymmetricMatrix_ImplEigenFloat::getNumOfElements() const {
    return this->matrix_.size();
}

float* TlDenseSymmetricMatrix_ImplEigenFloat::data() {
    return this->matrix_.data();
}

const float* TlDenseSymmetricMatrix_ImplEigenFloat::data() const {
    return this->matrix_.data();
}

// ---------------------------------------------------------------------------
// others
// ---------------------------------------------------------------------------
// DM(G) * DM(S)
TlDenseGeneralMatrix_ImplEigenFloat operator*(
    const TlDenseGeneralMatrix_ImplEigenFloat& mat1,
    const TlDenseSymmetricMatrix_ImplEigenFloat& mat2) {
    TlDenseGeneralMatrix_ImplEigenFloat answer;
    answer.matrix_ = mat1.matrix_ * mat2.matrix_;
    return answer;
}

// DM(S) * DM(G)
TlDenseGeneralMatrix_ImplEigenFloat operator*(
    const TlDenseSymmetricMatrix_ImplEigenFloat& mat1,
    const TlDenseGeneralMatrix_ImplEigenFloat& mat2) {
    TlDenseGeneralMatrix_ImplEigenFloat answer;
    answer.matrix_ = mat1.matrix_ * mat2.matrix_;
    return answer;
}

// DM(S) * DV
TlDenseVector_ImplEigenFloat operator*(const TlDenseSymmetricMatrix_ImplEigenFloat& dms,
                                       const TlDenseVector_ImplEigenFloat& dv) {
    TlDenseVector_ImplEigenFloat answer;
    answer.vector_ = dms.matrix_ * dv.vector_;
    return answer;
}

// DV * DM(S)
TlDenseVector_ImplEigenFloat operator*(const TlDenseVector_ImplEigenFloat& dv,
                                       const TlDenseSymmetricMatrix_ImplEigenFloat& dms) {
    TlDenseVector_ImplEigenFloat answer;
    answer.vector_ = dv.vector_ * dms.matrix_;
    return answer;
}

TlDenseSymmetricMatrix_ImplEigenFloat operator*(
    const double coef, const TlDenseSymmetricMatrix_ImplEigenFloat& DM) {
    TlDenseSymmetricMatrix_ImplEigenFloat answer = DM;
    answer *= coef;
    return answer;
}

TlDenseSymmetricMatrix_ImplEigenFloat operator*(
    const TlDenseSymmetricMatrix_ImplEigenFloat& DM, const double coef) {
    TlDenseSymmetricMatrix_ImplEigenFloat answer = DM;
    answer *= coef;
    return answer;
}
