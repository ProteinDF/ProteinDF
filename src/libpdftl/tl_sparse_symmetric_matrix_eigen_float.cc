#include "tl_sparse_symmetric_matrix_eigen_float.h"

#include "tl_dense_general_matrix_eigen_float.h"
#include "tl_dense_general_matrix_impl_eigen_float.h"
#include "tl_dense_symmetric_matrix_eigen_float.h"
#include "tl_dense_symmetric_matrix_impl_eigen_float.h"
#include "tl_dense_vector_eigen_float.h"
#include "tl_dense_vector_impl_eigen_float.h"
#include "tl_sparse_general_matrix_eigen_float.h"
#include "tl_sparse_general_matrix_impl_eigen_float.h"
#include "tl_sparse_symmetric_matrix_impl_eigen_float.h"

#ifdef HAVE_VIENNACL
#include "tl_sparse_symmetric_matrix_impl_viennacl_float.h"
#include "tl_sparse_symmetric_matrix_viennacl_float.h"
#endif  // HAVE_VIENNACL

// ---------------------------------------------------------------------------
// constructor & destructor
// ---------------------------------------------------------------------------
TlSparseSymmetricMatrix_EigenFloat::TlSparseSymmetricMatrix_EigenFloat(const TlMatrixObject::index_type dim) {
    this->pImpl_ = new TlSparseSymmetricMatrix_ImplEigenFloat(dim);
}

TlSparseSymmetricMatrix_EigenFloat::TlSparseSymmetricMatrix_EigenFloat(const TlSparseSymmetricMatrix_EigenFloat& rhs) {
    this->pImpl_ = new TlSparseSymmetricMatrix_ImplEigenFloat(
        *(dynamic_cast<const TlSparseSymmetricMatrix_ImplEigenFloat*>(rhs.pImpl_)));
}

TlSparseSymmetricMatrix_EigenFloat::TlSparseSymmetricMatrix_EigenFloat(const TlSparseGeneralMatrix_EigenFloat& rhs) {
    this->pImpl_ = new TlSparseSymmetricMatrix_ImplEigenFloat(
        *(dynamic_cast<const TlSparseGeneralMatrix_ImplEigenFloat*>(rhs.pImpl_)));
}

TlSparseSymmetricMatrix_EigenFloat::TlSparseSymmetricMatrix_EigenFloat(const TlDenseSymmetricMatrix_EigenFloat& rhs) {
    this->pImpl_ = new TlSparseSymmetricMatrix_ImplEigenFloat(
        *(dynamic_cast<TlDenseSymmetricMatrix_ImplEigenFloat*>(rhs.pImpl_)));
}

#ifdef HAVE_VIENNACL
TlSparseSymmetricMatrix_EigenFloat::TlSparseSymmetricMatrix_EigenFloat(const TlSparseSymmetricMatrix_ViennaCLFloat& rhs) {
    this->pImpl_ = new TlSparseSymmetricMatrix_ImplEigenFloat(*(dynamic_cast<TlSparseSymmetricMatrix_ImplViennaCLFloat*>(rhs.pImpl_)));
}
#endif  // HAVE_VIENNACL

TlSparseSymmetricMatrix_EigenFloat::~TlSparseSymmetricMatrix_EigenFloat() {
    delete this->pImpl_;
    this->pImpl_ = NULL;
}

// ---------------------------------------------------------------------------
// operators
// ---------------------------------------------------------------------------
TlSparseSymmetricMatrix_EigenFloat& TlSparseSymmetricMatrix_EigenFloat::operator=(
    const TlSparseSymmetricMatrix_EigenFloat& rhs) {
    if (this != &rhs) {
        delete this->pImpl_;
        this->pImpl_ = new TlSparseSymmetricMatrix_ImplEigenFloat(
            *(dynamic_cast<TlSparseSymmetricMatrix_ImplEigenFloat*>(rhs.pImpl_)));
    }

    return *this;
}

TlSparseSymmetricMatrix_EigenFloat& TlSparseSymmetricMatrix_EigenFloat::operator+=(
    const TlSparseSymmetricMatrix_EigenFloat& sm) {
    *(dynamic_cast<TlSparseSymmetricMatrix_ImplEigenFloat*>(this->pImpl_)) +=
        *(dynamic_cast<TlSparseSymmetricMatrix_ImplEigenFloat*>(sm.pImpl_));
    return *this;
}

TlSparseSymmetricMatrix_EigenFloat& TlSparseSymmetricMatrix_EigenFloat::operator-=(
    const TlSparseSymmetricMatrix_EigenFloat& sm) {
    *(dynamic_cast<TlSparseSymmetricMatrix_ImplEigenFloat*>(this->pImpl_)) -=
        *(dynamic_cast<TlSparseSymmetricMatrix_ImplEigenFloat*>(sm.pImpl_));
    return *this;
}

// ---------------------------------------------------------------------------
// others
// ---------------------------------------------------------------------------
// DM(G) = DM(G) * SM(S)
TlDenseGeneralMatrix_EigenFloat operator*(const TlDenseGeneralMatrix_EigenFloat& dm1,
                                          const TlSparseSymmetricMatrix_EigenFloat& sm2) {
    return TlDenseGeneralMatrix_EigenFloat(*(dynamic_cast<TlDenseGeneralMatrix_ImplEigenFloat*>(dm1.pImpl_)) *
                                           *(dynamic_cast<TlSparseSymmetricMatrix_ImplEigenFloat*>(sm2.pImpl_)));
}

// DM(G) = SM(S) * DM(G)
TlDenseGeneralMatrix_EigenFloat operator*(const TlSparseSymmetricMatrix_EigenFloat& sm1,
                                          const TlDenseGeneralMatrix_EigenFloat& dm2) {
    return TlDenseGeneralMatrix_EigenFloat(*(dynamic_cast<TlSparseSymmetricMatrix_ImplEigenFloat*>(sm1.pImpl_)) *
                                           *(dynamic_cast<TlDenseGeneralMatrix_ImplEigenFloat*>(dm2.pImpl_)));
}

// SM(G) = SM(G) * SM(S)
TlSparseGeneralMatrix_EigenFloat operator*(
    const TlSparseGeneralMatrix_EigenFloat& sm1,
    const TlSparseSymmetricMatrix_EigenFloat& sm2) {
    return TlSparseGeneralMatrix_EigenFloat(*(dynamic_cast<TlSparseGeneralMatrix_ImplEigenFloat*>(sm1.pImpl_)) *
                                            *(dynamic_cast<TlSparseSymmetricMatrix_ImplEigenFloat*>(sm2.pImpl_)));
}

// SM(G) = SM(S) * SM(G)
TlSparseGeneralMatrix_EigenFloat operator*(const TlSparseSymmetricMatrix_EigenFloat& sm1,
                                           const TlSparseGeneralMatrix_EigenFloat& sm2) {
    return TlSparseGeneralMatrix_EigenFloat(*(dynamic_cast<TlSparseSymmetricMatrix_ImplEigenFloat*>(sm1.pImpl_)) *
                                            *(dynamic_cast<TlSparseGeneralMatrix_ImplEigenFloat*>(sm2.pImpl_)));
}

// SM(G) = SM(S) * SM(S)
TlSparseGeneralMatrix_EigenFloat operator*(
    const TlSparseSymmetricMatrix_EigenFloat& sm1,
    const TlSparseSymmetricMatrix_EigenFloat& sm2) {
    return TlSparseGeneralMatrix_EigenFloat(*(dynamic_cast<TlSparseSymmetricMatrix_ImplEigenFloat*>(sm1.pImpl_)) *
                                            *(dynamic_cast<TlSparseSymmetricMatrix_ImplEigenFloat*>(sm2.pImpl_)));
}

// SM(S) * DV
TlDenseVector_EigenFloat operator*(const TlSparseSymmetricMatrix_EigenFloat& mat,
                                   const TlDenseVector_EigenFloat& vtr) {
    return TlDenseVector_EigenFloat(*(dynamic_cast<TlSparseSymmetricMatrix_ImplEigenFloat*>(mat.pImpl_)) *
                                    *(dynamic_cast<TlDenseVector_ImplEigenFloat*>(vtr.pImpl_)));
}
