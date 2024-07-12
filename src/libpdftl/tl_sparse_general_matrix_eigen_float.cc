#include "tl_sparse_general_matrix_eigen.h"
#include "tl_dense_general_matrix_eigen.h"
#include "tl_dense_general_matrix_impl_eigen.h"
#include "tl_dense_vector_eigen.h"
#include "tl_dense_vector_impl_eigen.h"
#include "tl_sparse_general_matrix_impl_eigen.h"
#include "tl_sparse_symmetric_matrix_eigen.h"
#include "tl_sparse_symmetric_matrix_impl_eigen.h"

#ifdef HAVE_VIENNACL
#include "tl_sparse_general_matrix_impl_viennacl.h"
#include "tl_sparse_general_matrix_viennacl.h"
#endif  // HAVE_VIENNACL

// ---------------------------------------------------------------------------
// constructor & destructor
// ---------------------------------------------------------------------------
TlSparseGeneralMatrix_Eigen::TlSparseGeneralMatrix_Eigen(
    const TlMatrixObject::index_type row,
    const TlMatrixObject::index_type col) {
    this->pImpl_ = new TlSparseGeneralMatrix_ImplEigen(row, col);
}

TlSparseGeneralMatrix_Eigen::TlSparseGeneralMatrix_Eigen(
    const TlSparseGeneralMatrix_Eigen& rhs) {
    this->pImpl_ = new TlSparseGeneralMatrix_ImplEigen(
        *(dynamic_cast<const TlSparseGeneralMatrix_ImplEigen*>(rhs.pImpl_)));
}

TlSparseGeneralMatrix_Eigen::TlSparseGeneralMatrix_Eigen(
    const TlSparseSymmetricMatrix_Eigen& rhs) {
    this->pImpl_ = new TlSparseGeneralMatrix_ImplEigen(
        *(dynamic_cast<const TlSparseSymmetricMatrix_ImplEigen*>(rhs.pImpl_)));
}

TlSparseGeneralMatrix_Eigen::TlSparseGeneralMatrix_Eigen(
    const TlSparseGeneralMatrix_ImplEigen& rhs) {
    this->pImpl_ = new TlSparseGeneralMatrix_ImplEigen(rhs);
}

TlSparseGeneralMatrix_Eigen::TlSparseGeneralMatrix_Eigen(
    const TlDenseGeneralMatrix_Eigen& rhs) {
    this->pImpl_ = new TlSparseGeneralMatrix_ImplEigen(
        *(dynamic_cast<TlDenseGeneralMatrix_ImplEigen*>(rhs.pImpl_)));
}

#ifdef HAVE_VIENNACL
TlSparseGeneralMatrix_Eigen::TlSparseGeneralMatrix_Eigen(
    const TlSparseGeneralMatrix_ViennaCL& rhs) {
    std::cout << "TlSparseGeneralMatrix_Eigen::TlSparseGeneralMatrix_Eigen("
                 "const TlSparseGeneralMatrix_ViennaCL& rhs)"
              << std::endl;
    this->pImpl_ = new TlSparseGeneralMatrix_ImplEigen(
        *(dynamic_cast<TlSparseGeneralMatrix_ImplViennaCL*>(rhs.pImpl_)));
}
#endif  // HAVE_VIENNACL

TlSparseGeneralMatrix_Eigen::~TlSparseGeneralMatrix_Eigen() {
    delete this->pImpl_;
    this->pImpl_ = NULL;
}

// ---------------------------------------------------------------------------
// operators
// ---------------------------------------------------------------------------
TlSparseGeneralMatrix_Eigen& TlSparseGeneralMatrix_Eigen::operator=(
    const TlSparseGeneralMatrix_Eigen& rhs) {
    if (this != &rhs) {
        delete this->pImpl_;
        this->pImpl_ = new TlSparseGeneralMatrix_ImplEigen(
            *(dynamic_cast<TlSparseGeneralMatrix_ImplEigen*>(rhs.pImpl_)));
    }

    return *this;
}

TlSparseGeneralMatrix_Eigen& TlSparseGeneralMatrix_Eigen::operator+=(
    const TlSparseGeneralMatrix_Eigen& rhs) {
    *(dynamic_cast<TlSparseGeneralMatrix_ImplEigen*>(this->pImpl_)) +=
        *(dynamic_cast<TlSparseGeneralMatrix_ImplEigen*>(rhs.pImpl_));
    return *this;
}

TlSparseGeneralMatrix_Eigen& TlSparseGeneralMatrix_Eigen::operator-=(
    const TlSparseGeneralMatrix_Eigen& rhs) {
    *(dynamic_cast<TlSparseGeneralMatrix_ImplEigen*>(this->pImpl_)) -=
        *(dynamic_cast<TlSparseGeneralMatrix_ImplEigen*>(rhs.pImpl_));
    return *this;
}

TlSparseGeneralMatrix_Eigen& TlSparseGeneralMatrix_Eigen::operator*=(
    const double coef) {
    *(dynamic_cast<TlSparseGeneralMatrix_ImplEigen*>(this->pImpl_)) *= coef;
    return *this;
}

// ----------------------------------------------------------------------------
// others
// ----------------------------------------------------------------------------
// DM(G) = DM(G) * SM(G)
TlDenseGeneralMatrix_Eigen operator*(const TlDenseGeneralMatrix_Eigen& mat1,
                                     const TlSparseGeneralMatrix_Eigen& mat2) {
    return TlDenseGeneralMatrix_Eigen(
        *(dynamic_cast<TlDenseGeneralMatrix_ImplEigen*>(mat1.pImpl_)) *
        *(dynamic_cast<TlSparseGeneralMatrix_ImplEigen*>(mat2.pImpl_)));
}

// DM(G) = SM(G) * DM(G)
TlDenseGeneralMatrix_Eigen operator*(const TlSparseGeneralMatrix_Eigen& mat1,
                                     const TlDenseGeneralMatrix_Eigen& mat2) {
    return TlDenseGeneralMatrix_Eigen(
        *(dynamic_cast<TlSparseGeneralMatrix_ImplEigen*>(mat1.pImpl_)) *
        *(dynamic_cast<TlDenseGeneralMatrix_ImplEigen*>(mat2.pImpl_)));
}

// SM(G) = SM(S) * SM(G)
TlSparseGeneralMatrix_Eigen operator*(const TlSparseGeneralMatrix_Eigen& sm1,
                                      const TlSparseGeneralMatrix_Eigen& sm2) {
    return TlSparseGeneralMatrix_Eigen(
        *(dynamic_cast<TlSparseGeneralMatrix_ImplEigen*>(sm1.pImpl_)) *
        *(dynamic_cast<TlSparseGeneralMatrix_ImplEigen*>(sm2.pImpl_)));
}

// DV = SM(G) * DV
TlDenseVector_Eigen operator*(const TlSparseGeneralMatrix_Eigen& mat,
                              const TlDenseVector_Eigen& vtr) {
    return TlDenseVector_Eigen(
        *(dynamic_cast<TlSparseGeneralMatrix_ImplEigen*>(mat.pImpl_)) *
        *(dynamic_cast<TlDenseVector_ImplEigen*>(vtr.pImpl_)));
}
