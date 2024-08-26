#ifdef HAVE_CONFIG_H
#include "config.h"
#endif  // HAVE_CONFIG_H

#include "tl_dense_general_matrix_impl_viennacl_float.h"
#include "tl_dense_general_matrix_viennacl_float.h"
#include "tl_dense_vector_impl_viennacl_float.h"
#include "tl_dense_vector_viennacl_float.h"
#include "tl_sparse_general_matrix_impl_viennacl_float.h"
#include "tl_sparse_general_matrix_viennacl_float.h"
#include "tl_sparse_symmetric_matrix_impl_viennacl_float.h"
#include "tl_sparse_symmetric_matrix_viennacl_float.h"

#ifdef HAVE_EIGEN
#include "tl_sparse_general_matrix_eigen_float.h"
#include "tl_sparse_general_matrix_impl_eigen_float.h"
#endif  // HAVE_EIGEN

// ---------------------------------------------------------------------------
// constructor & destructor
// ---------------------------------------------------------------------------A
TlSparseGeneralMatrix_ViennaCLFloat::TlSparseGeneralMatrix_ViennaCLFloat(
    const TlMatrixObject::index_type row,
    const TlMatrixObject::index_type col) {
    this->pImpl_ = new TlSparseGeneralMatrix_ImplViennaCLFloat(row, col);
}

TlSparseGeneralMatrix_ViennaCLFloat::TlSparseGeneralMatrix_ViennaCLFloat(
    const TlSparseGeneralMatrix_ViennaCLFloat& rhs) {
    this->pImpl_ = new TlSparseGeneralMatrix_ImplViennaCLFloat(
        *(dynamic_cast<const TlSparseGeneralMatrix_ImplViennaCLFloat*>(rhs.pImpl_)));
}

TlSparseGeneralMatrix_ViennaCLFloat::TlSparseGeneralMatrix_ViennaCLFloat(
    const TlSparseSymmetricMatrix_ViennaCLFloat& rhs) {
    this->pImpl_ = new TlSparseGeneralMatrix_ImplViennaCLFloat(*(
        dynamic_cast<const TlSparseSymmetricMatrix_ImplViennaCLFloat*>(rhs.pImpl_)));
}

TlSparseGeneralMatrix_ViennaCLFloat::TlSparseGeneralMatrix_ViennaCLFloat(
    const TlSparseGeneralMatrix_ImplViennaCLFloat& rhs) {
    this->pImpl_ = new TlSparseGeneralMatrix_ImplViennaCLFloat(rhs);
}

TlSparseGeneralMatrix_ViennaCLFloat::TlSparseGeneralMatrix_ViennaCLFloat(const TlDenseGeneralMatrix_ViennaCLFloat& rhs) {
    this->pImpl_ = new TlSparseGeneralMatrix_ImplViennaCLFloat(*(dynamic_cast<TlDenseGeneralMatrix_ImplViennaCLFloat*>(rhs.pImpl_)));
}

#ifdef HAVE_VIENNACL
TlSparseGeneralMatrix_ViennaCLFloat::TlSparseGeneralMatrix_ViennaCLFloat(const TlSparseGeneralMatrix_EigenFloat& rhs) {
    this->pImpl_ = new TlSparseGeneralMatrix_ImplViennaCLFloat(*(dynamic_cast<TlSparseGeneralMatrix_ImplEigenFloat*>(rhs.pImpl_)));
}
#endif  // HAVE_VIENNACL

TlSparseGeneralMatrix_ViennaCLFloat::~TlSparseGeneralMatrix_ViennaCLFloat() {
    delete this->pImpl_;
    this->pImpl_ = NULL;
}

// ---------------------------------------------------------------------------
// operators
// ---------------------------------------------------------------------------
TlSparseGeneralMatrix_ViennaCLFloat& TlSparseGeneralMatrix_ViennaCLFloat::operator=(
    const TlSparseGeneralMatrix_ViennaCLFloat& rhs) {
    if (this != &rhs) {
        delete this->pImpl_;
        this->pImpl_ = new TlSparseGeneralMatrix_ImplViennaCLFloat(
            *(dynamic_cast<TlSparseGeneralMatrix_ImplViennaCLFloat*>(rhs.pImpl_)));
    }

    return *this;
}

TlSparseGeneralMatrix_ViennaCLFloat& TlSparseGeneralMatrix_ViennaCLFloat::operator+=(
    const TlSparseGeneralMatrix_ViennaCLFloat& rhs) {
    *(dynamic_cast<TlSparseGeneralMatrix_ImplViennaCLFloat*>(this->pImpl_)) +=
        *(dynamic_cast<TlSparseGeneralMatrix_ImplViennaCLFloat*>(rhs.pImpl_));
    return *this;
}

TlSparseGeneralMatrix_ViennaCLFloat& TlSparseGeneralMatrix_ViennaCLFloat::operator-=(
    const TlSparseGeneralMatrix_ViennaCLFloat& rhs) {
    *(dynamic_cast<TlSparseGeneralMatrix_ImplViennaCLFloat*>(this->pImpl_)) -=
        *(dynamic_cast<TlSparseGeneralMatrix_ImplViennaCLFloat*>(rhs.pImpl_));
    return *this;
}

TlSparseGeneralMatrix_ViennaCLFloat& TlSparseGeneralMatrix_ViennaCLFloat::operator*=(
    const double coef) {
    *(dynamic_cast<TlSparseGeneralMatrix_ImplViennaCLFloat*>(this->pImpl_)) *= coef;
    return *this;
}

// ----------------------------------------------------------------------------
// others
// ----------------------------------------------------------------------------
TlSparseGeneralMatrix_ViennaCLFloat operator*(
    const TlSparseGeneralMatrix_ViennaCLFloat& sm1,
    const TlSparseGeneralMatrix_ViennaCLFloat& sm2) {
    return TlSparseGeneralMatrix_ViennaCLFloat(
        *(dynamic_cast<TlSparseGeneralMatrix_ImplViennaCLFloat*>(sm1.pImpl_)) *
        *(dynamic_cast<TlSparseGeneralMatrix_ImplViennaCLFloat*>(sm2.pImpl_)));
}

TlDenseVector_ViennaCLFloat operator*(const TlSparseGeneralMatrix_ViennaCLFloat& mat,
                                      const TlDenseVector_ViennaCLFloat& vtr) {
    // std::cout << "operator*(VCL)" << std::endl;
    return TlDenseVector_ViennaCLFloat(
        *(dynamic_cast<TlSparseGeneralMatrix_ImplViennaCLFloat*>(mat.pImpl_)) *
        *(dynamic_cast<TlDenseVector_ImplViennaCLFloat*>(vtr.pImpl_)));
}
