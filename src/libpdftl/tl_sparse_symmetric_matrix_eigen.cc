#include "tl_sparse_symmetric_matrix_eigen.h"
#include "tl_dense_general_matrix_eigen.h"
#include "tl_dense_general_matrix_impl_eigen.h"
#include "tl_dense_symmetric_matrix_eigen.h"
#include "tl_dense_symmetric_matrix_impl_eigen.h"
#include "tl_dense_vector_eigen.h"
#include "tl_dense_vector_impl_eigen.h"
#include "tl_sparse_general_matrix_eigen.h"
#include "tl_sparse_general_matrix_impl_eigen.h"
#include "tl_sparse_symmetric_matrix_impl_eigen.h"

#ifdef HAVE_VIENNACL
#include "tl_sparse_symmetric_matrix_impl_viennacl.h"
#include "tl_sparse_symmetric_matrix_viennacl.h"
#endif  // HAVE_VIENNACL

// ---------------------------------------------------------------------------
// constructor & destructor
// ---------------------------------------------------------------------------
TlSparseSymmetricMatrix_Eigen::TlSparseSymmetricMatrix_Eigen(
    const TlMatrixObject::index_type dim) {
    this->pImpl_ = new TlSparseSymmetricMatrix_ImplEigen(dim);
}

TlSparseSymmetricMatrix_Eigen::TlSparseSymmetricMatrix_Eigen(
    const TlSparseSymmetricMatrix_Eigen& rhs) {
    this->pImpl_ = new TlSparseSymmetricMatrix_ImplEigen(
        *(dynamic_cast<const TlSparseSymmetricMatrix_ImplEigen*>(rhs.pImpl_)));
}

TlSparseSymmetricMatrix_Eigen::TlSparseSymmetricMatrix_Eigen(
    const TlSparseGeneralMatrix_Eigen& rhs) {
    this->pImpl_ = new TlSparseSymmetricMatrix_ImplEigen(
        *(dynamic_cast<const TlSparseGeneralMatrix_ImplEigen*>(rhs.pImpl_)));
}

TlSparseSymmetricMatrix_Eigen::TlSparseSymmetricMatrix_Eigen(
    const TlDenseSymmetricMatrix_Eigen& rhs) {
    this->pImpl_ = new TlSparseSymmetricMatrix_ImplEigen(
        *(dynamic_cast<TlDenseSymmetricMatrix_ImplEigen*>(rhs.pImpl_)));
}

#ifdef HAVE_VIENNACL
TlSparseSymmetricMatrix_Eigen::TlSparseSymmetricMatrix_Eigen(
    const TlSparseSymmetricMatrix_ViennaCL& rhs) {
    this->pImpl_ = new TlSparseSymmetricMatrix_ImplEigen(
        *(dynamic_cast<TlSparseSymmetricMatrix_ImplViennaCL*>(rhs.pImpl_)));
}
#endif  // HAVE_VIENNACL

TlSparseSymmetricMatrix_Eigen::~TlSparseSymmetricMatrix_Eigen() {
    delete this->pImpl_;
    this->pImpl_ = NULL;
}

// ---------------------------------------------------------------------------
// operators
// ---------------------------------------------------------------------------
TlSparseSymmetricMatrix_Eigen& TlSparseSymmetricMatrix_Eigen::operator=(
    const TlSparseSymmetricMatrix_Eigen& rhs) {
    if (this != &rhs) {
        delete this->pImpl_;
        this->pImpl_ = new TlSparseSymmetricMatrix_ImplEigen(
            *(dynamic_cast<TlSparseSymmetricMatrix_ImplEigen*>(rhs.pImpl_)));
    }

    return *this;
}

TlSparseSymmetricMatrix_Eigen& TlSparseSymmetricMatrix_Eigen::operator+=(
    const TlSparseSymmetricMatrix_Eigen& sm) {
    *(dynamic_cast<TlSparseSymmetricMatrix_ImplEigen*>(this->pImpl_)) +=
        *(dynamic_cast<TlSparseSymmetricMatrix_ImplEigen*>(sm.pImpl_));
    return *this;
}

TlSparseSymmetricMatrix_Eigen& TlSparseSymmetricMatrix_Eigen::operator-=(
    const TlSparseSymmetricMatrix_Eigen& sm) {
    *(dynamic_cast<TlSparseSymmetricMatrix_ImplEigen*>(this->pImpl_)) -=
        *(dynamic_cast<TlSparseSymmetricMatrix_ImplEigen*>(sm.pImpl_));
    return *this;
}

// ---------------------------------------------------------------------------
// others
// ---------------------------------------------------------------------------
// DM(G) = DM(G) * SM(S)
TlDenseGeneralMatrix_Eigen operator*(const TlDenseGeneralMatrix_Eigen& dm1,
                                     const TlSparseSymmetricMatrix_Eigen& sm2) {
    return TlDenseGeneralMatrix_Eigen(
        *(dynamic_cast<TlDenseGeneralMatrix_ImplEigen*>(dm1.pImpl_)) *
        *(dynamic_cast<TlSparseSymmetricMatrix_ImplEigen*>(sm2.pImpl_)));
}

// DM(G) = SM(S) * DM(G)
TlDenseGeneralMatrix_Eigen operator*(const TlSparseSymmetricMatrix_Eigen& sm1,
                                     const TlDenseGeneralMatrix_Eigen& dm2) {
    return TlDenseGeneralMatrix_Eigen(
        *(dynamic_cast<TlSparseSymmetricMatrix_ImplEigen*>(sm1.pImpl_)) *
        *(dynamic_cast<TlDenseGeneralMatrix_ImplEigen*>(dm2.pImpl_)));
}

// SM(G) = SM(G) * SM(S)
TlSparseGeneralMatrix_Eigen operator*(
    const TlSparseGeneralMatrix_Eigen& sm1,
    const TlSparseSymmetricMatrix_Eigen& sm2) {
    return TlSparseGeneralMatrix_Eigen(
        *(dynamic_cast<TlSparseGeneralMatrix_ImplEigen*>(sm1.pImpl_)) *
        *(dynamic_cast<TlSparseSymmetricMatrix_ImplEigen*>(sm2.pImpl_)));
}

// SM(G) = SM(S) * SM(G)
TlSparseGeneralMatrix_Eigen operator*(const TlSparseSymmetricMatrix_Eigen& sm1,
                                      const TlSparseGeneralMatrix_Eigen& sm2) {
    return TlSparseGeneralMatrix_Eigen(
        *(dynamic_cast<TlSparseSymmetricMatrix_ImplEigen*>(sm1.pImpl_)) *
        *(dynamic_cast<TlSparseGeneralMatrix_ImplEigen*>(sm2.pImpl_)));
}

// SM(G) = SM(S) * SM(S)
TlSparseGeneralMatrix_Eigen operator*(
    const TlSparseSymmetricMatrix_Eigen& sm1,
    const TlSparseSymmetricMatrix_Eigen& sm2) {
    return TlSparseGeneralMatrix_Eigen(
        *(dynamic_cast<TlSparseSymmetricMatrix_ImplEigen*>(sm1.pImpl_)) *
        *(dynamic_cast<TlSparseSymmetricMatrix_ImplEigen*>(sm2.pImpl_)));
}

// SM(S) * DV
TlDenseVector_Eigen operator*(const TlSparseSymmetricMatrix_Eigen& mat,
                              const TlDenseVector_Eigen& vtr) {
    return TlDenseVector_Eigen(
        *(dynamic_cast<TlSparseSymmetricMatrix_ImplEigen*>(mat.pImpl_)) *
        *(dynamic_cast<TlDenseVector_ImplEigen*>(vtr.pImpl_)));
}
