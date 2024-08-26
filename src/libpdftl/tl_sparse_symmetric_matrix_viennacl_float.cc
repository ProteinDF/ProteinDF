#ifdef HAVE_CONFIG_H
#include "config.h"
#endif  // HAVE_CONFIG_H

#include "tl_dense_symmetric_matrix_impl_viennacl_float.h"
#include "tl_dense_symmetric_matrix_viennacl_float.h"
#include "tl_sparse_general_matrix_impl_viennacl_float.h"
#include "tl_sparse_general_matrix_viennacl_float.h"
#include "tl_sparse_symmetric_matrix_impl_viennacl_float.h"
#include "tl_sparse_symmetric_matrix_viennacl_float.h"

#ifdef HAVE_EIGEN
#include "tl_sparse_symmetric_matrix_eigen_float.h"
#include "tl_sparse_symmetric_matrix_impl_eigen_float.h"
#endif  // HAVE_EIGEN

// ---------------------------------------------------------------------------
// constructor & destructor
// ---------------------------------------------------------------------------
TlSparseSymmetricMatrix_ViennaCLFloat::TlSparseSymmetricMatrix_ViennaCLFloat(
    const TlMatrixObject::index_type dim) {
    this->pImpl_ = new TlSparseSymmetricMatrix_ImplViennaCLFloat(dim);
}

TlSparseSymmetricMatrix_ViennaCLFloat::TlSparseSymmetricMatrix_ViennaCLFloat(
    const TlSparseSymmetricMatrix_ViennaCLFloat& rhs) {
    this->pImpl_ = new TlSparseSymmetricMatrix_ImplViennaCLFloat(*(
        dynamic_cast<const TlSparseSymmetricMatrix_ImplViennaCLFloat*>(rhs.pImpl_)));
}

TlSparseSymmetricMatrix_ViennaCLFloat::TlSparseSymmetricMatrix_ViennaCLFloat(
    const TlSparseGeneralMatrix_ViennaCLFloat& rhs) {
    this->pImpl_ = new TlSparseSymmetricMatrix_ImplViennaCLFloat(
        *(dynamic_cast<const TlSparseGeneralMatrix_ImplViennaCLFloat*>(rhs.pImpl_)));
}

TlSparseSymmetricMatrix_ViennaCLFloat::TlSparseSymmetricMatrix_ViennaCLFloat(
    const TlDenseSymmetricMatrix_ViennaCLFloat& rhs) {
    this->pImpl_ = new TlSparseSymmetricMatrix_ImplViennaCLFloat(
        *(dynamic_cast<TlDenseSymmetricMatrix_ImplViennaCLFloat*>(rhs.pImpl_)));
}

#ifdef HAVE_EIGEN
TlSparseSymmetricMatrix_ViennaCLFloat::TlSparseSymmetricMatrix_ViennaCLFloat(
    const TlSparseSymmetricMatrix_EigenFloat& rhs) {
    this->pImpl_ = new TlSparseSymmetricMatrix_ImplViennaCLFloat(
        *(dynamic_cast<TlSparseSymmetricMatrix_ImplEigenFloat*>(rhs.pImpl_)));
}
#endif  // HAVE_EIGEN

TlSparseSymmetricMatrix_ViennaCLFloat::~TlSparseSymmetricMatrix_ViennaCLFloat() {
    delete this->pImpl_;
    this->pImpl_ = NULL;
}

// ---------------------------------------------------------------------------
// operators
// ---------------------------------------------------------------------------
TlSparseSymmetricMatrix_ViennaCLFloat& TlSparseSymmetricMatrix_ViennaCLFloat::operator=(
    const TlSparseSymmetricMatrix_ViennaCLFloat& rhs) {
    if (this != &rhs) {
        delete this->pImpl_;
        this->pImpl_ = new TlSparseSymmetricMatrix_ImplViennaCLFloat(
            *(dynamic_cast<TlSparseSymmetricMatrix_ImplViennaCLFloat*>(rhs.pImpl_)));
    }

    return *this;
}

TlSparseSymmetricMatrix_ViennaCLFloat& TlSparseSymmetricMatrix_ViennaCLFloat::operator+=(
    const TlSparseSymmetricMatrix_ViennaCLFloat& rhs) {
    *(dynamic_cast<TlSparseSymmetricMatrix_ImplViennaCLFloat*>(this->pImpl_)) +=
        *(dynamic_cast<TlSparseSymmetricMatrix_ImplViennaCLFloat*>(rhs.pImpl_));
    return *this;
}

TlSparseSymmetricMatrix_ViennaCLFloat& TlSparseSymmetricMatrix_ViennaCLFloat::operator-=(
    const TlSparseSymmetricMatrix_ViennaCLFloat& rhs) {
    *(dynamic_cast<TlSparseSymmetricMatrix_ImplViennaCLFloat*>(this->pImpl_)) -=
        *(dynamic_cast<TlSparseSymmetricMatrix_ImplViennaCLFloat*>(rhs.pImpl_));
    return *this;
}

TlSparseSymmetricMatrix_ViennaCLFloat& TlSparseSymmetricMatrix_ViennaCLFloat::operator*=(
    const double coef) {
    *(dynamic_cast<TlSparseSymmetricMatrix_ImplViennaCLFloat*>(this->pImpl_)) *=
        coef;
    return *this;
}

// ----------------------------------------------------------------------------
// others
// ----------------------------------------------------------------------------
// SM(G) = SM(G) * SM(S)
TlSparseGeneralMatrix_ViennaCLFloat operator*(
    const TlSparseGeneralMatrix_ViennaCLFloat& sm1,
    const TlSparseSymmetricMatrix_ViennaCLFloat& sm2) {
    return TlSparseGeneralMatrix_ViennaCLFloat(
        *(dynamic_cast<TlSparseGeneralMatrix_ImplViennaCLFloat*>(sm1.pImpl_)) *
        *(dynamic_cast<TlSparseSymmetricMatrix_ImplViennaCLFloat*>(sm2.pImpl_)));
}

// SM(G) = SM(S) * SM(G)
TlSparseGeneralMatrix_ViennaCLFloat operator*(
    const TlSparseSymmetricMatrix_ViennaCLFloat& sm1,
    const TlSparseGeneralMatrix_ViennaCLFloat& sm2) {
    return TlSparseGeneralMatrix_ViennaCLFloat(
        *(dynamic_cast<TlSparseSymmetricMatrix_ImplViennaCLFloat*>(sm1.pImpl_)) *
        *(dynamic_cast<TlSparseGeneralMatrix_ImplViennaCLFloat*>(sm2.pImpl_)));
}

// SM(G) = SM(S) * SM(S)
TlSparseGeneralMatrix_ViennaCLFloat operator*(
    const TlSparseSymmetricMatrix_ViennaCLFloat& sm1,
    const TlSparseSymmetricMatrix_ViennaCLFloat& sm2) {
    return TlSparseGeneralMatrix_ViennaCLFloat(
        *(dynamic_cast<TlSparseSymmetricMatrix_ImplViennaCLFloat*>(sm1.pImpl_)) *
        *(dynamic_cast<TlSparseSymmetricMatrix_ImplViennaCLFloat*>(sm2.pImpl_)));
}
