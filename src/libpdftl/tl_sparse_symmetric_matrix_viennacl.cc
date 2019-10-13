#ifdef HAVE_CONFIG_H
#include "config.h"
#endif  // HAVE_CONFIG_H

#include "tl_dense_symmetric_matrix_impl_viennacl.h"
#include "tl_dense_symmetric_matrix_viennacl.h"
#include "tl_sparse_general_matrix_impl_viennacl.h"
#include "tl_sparse_general_matrix_viennacl.h"
#include "tl_sparse_symmetric_matrix_impl_viennacl.h"
#include "tl_sparse_symmetric_matrix_viennacl.h"

#ifdef HAVE_EIGEN
#include "tl_sparse_symmetric_matrix_eigen.h"
#include "tl_sparse_symmetric_matrix_impl_eigen.h"
#endif  // HAVE_EIGEN

// ---------------------------------------------------------------------------
// constructor & destructor
// ---------------------------------------------------------------------------
TlSparseSymmetricMatrix_ViennaCL::TlSparseSymmetricMatrix_ViennaCL(
    const TlMatrixObject::index_type dim) {
    this->pImpl_ = new TlSparseSymmetricMatrix_ImplViennaCL(dim);
}

TlSparseSymmetricMatrix_ViennaCL::TlSparseSymmetricMatrix_ViennaCL(
    const TlSparseSymmetricMatrix_ViennaCL& rhs) {
    this->pImpl_ = new TlSparseSymmetricMatrix_ImplViennaCL(*(
        dynamic_cast<const TlSparseSymmetricMatrix_ImplViennaCL*>(rhs.pImpl_)));
}

TlSparseSymmetricMatrix_ViennaCL::TlSparseSymmetricMatrix_ViennaCL(
    const TlSparseGeneralMatrix_ViennaCL& rhs) {
    this->pImpl_ = new TlSparseSymmetricMatrix_ImplViennaCL(
        *(dynamic_cast<const TlSparseGeneralMatrix_ImplViennaCL*>(rhs.pImpl_)));
}

TlSparseSymmetricMatrix_ViennaCL::TlSparseSymmetricMatrix_ViennaCL(
    const TlDenseSymmetricMatrix_ViennaCL& rhs) {
    this->pImpl_ = new TlSparseSymmetricMatrix_ImplViennaCL(
        *(dynamic_cast<TlDenseSymmetricMatrix_ImplViennaCL*>(rhs.pImpl_)));
}

#ifdef HAVE_EIGEN
TlSparseSymmetricMatrix_ViennaCL::TlSparseSymmetricMatrix_ViennaCL(
    const TlSparseSymmetricMatrix_Eigen& rhs) {
    this->pImpl_ = new TlSparseSymmetricMatrix_ImplViennaCL(
        *(dynamic_cast<TlSparseSymmetricMatrix_ImplEigen*>(rhs.pImpl_)));
}
#endif  // HAVE_EIGEN

TlSparseSymmetricMatrix_ViennaCL::~TlSparseSymmetricMatrix_ViennaCL() {
    delete this->pImpl_;
    this->pImpl_ = NULL;
}

// ---------------------------------------------------------------------------
// operators
// ---------------------------------------------------------------------------
TlSparseSymmetricMatrix_ViennaCL& TlSparseSymmetricMatrix_ViennaCL::operator=(
    const TlSparseSymmetricMatrix_ViennaCL& rhs) {
    if (this != &rhs) {
        delete this->pImpl_;
        this->pImpl_ = new TlSparseSymmetricMatrix_ImplViennaCL(
            *(dynamic_cast<TlSparseSymmetricMatrix_ImplViennaCL*>(rhs.pImpl_)));
    }

    return *this;
}

TlSparseSymmetricMatrix_ViennaCL& TlSparseSymmetricMatrix_ViennaCL::operator+=(
    const TlSparseSymmetricMatrix_ViennaCL& rhs) {
    *(dynamic_cast<TlSparseSymmetricMatrix_ImplViennaCL*>(this->pImpl_)) +=
        *(dynamic_cast<TlSparseSymmetricMatrix_ImplViennaCL*>(rhs.pImpl_));
    return *this;
}

TlSparseSymmetricMatrix_ViennaCL& TlSparseSymmetricMatrix_ViennaCL::operator-=(
    const TlSparseSymmetricMatrix_ViennaCL& rhs) {
    *(dynamic_cast<TlSparseSymmetricMatrix_ImplViennaCL*>(this->pImpl_)) -=
        *(dynamic_cast<TlSparseSymmetricMatrix_ImplViennaCL*>(rhs.pImpl_));
    return *this;
}

TlSparseSymmetricMatrix_ViennaCL& TlSparseSymmetricMatrix_ViennaCL::operator*=(
    const double coef) {
    *(dynamic_cast<TlSparseSymmetricMatrix_ImplViennaCL*>(this->pImpl_)) *=
        coef;
    return *this;
}

// ----------------------------------------------------------------------------
// others
// ----------------------------------------------------------------------------
// SM(G) = SM(G) * SM(S)
TlSparseGeneralMatrix_ViennaCL operator*(
    const TlSparseGeneralMatrix_ViennaCL& sm1,
    const TlSparseSymmetricMatrix_ViennaCL& sm2) {
    return TlSparseGeneralMatrix_ViennaCL(
        *(dynamic_cast<TlSparseGeneralMatrix_ImplViennaCL*>(sm1.pImpl_)) *
        *(dynamic_cast<TlSparseSymmetricMatrix_ImplViennaCL*>(sm2.pImpl_)));
}

// SM(G) = SM(S) * SM(G)
TlSparseGeneralMatrix_ViennaCL operator*(
    const TlSparseSymmetricMatrix_ViennaCL& sm1,
    const TlSparseGeneralMatrix_ViennaCL& sm2) {
    return TlSparseGeneralMatrix_ViennaCL(
        *(dynamic_cast<TlSparseSymmetricMatrix_ImplViennaCL*>(sm1.pImpl_)) *
        *(dynamic_cast<TlSparseGeneralMatrix_ImplViennaCL*>(sm2.pImpl_)));
}

// SM(G) = SM(S) * SM(S)
TlSparseGeneralMatrix_ViennaCL operator*(
    const TlSparseSymmetricMatrix_ViennaCL& sm1,
    const TlSparseSymmetricMatrix_ViennaCL& sm2) {
    return TlSparseGeneralMatrix_ViennaCL(
        *(dynamic_cast<TlSparseSymmetricMatrix_ImplViennaCL*>(sm1.pImpl_)) *
        *(dynamic_cast<TlSparseSymmetricMatrix_ImplViennaCL*>(sm2.pImpl_)));
}
