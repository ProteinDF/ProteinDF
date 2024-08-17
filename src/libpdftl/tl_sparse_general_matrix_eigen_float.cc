#include "tl_sparse_general_matrix_eigen_float.h"

#include "tl_dense_general_matrix_eigen_float.h"
#include "tl_dense_general_matrix_impl_eigen_float.h"
#include "tl_dense_vector_eigen_float.h"
#include "tl_dense_vector_impl_eigen_float.h"
#include "tl_sparse_general_matrix_impl_eigen_float.h"
#include "tl_sparse_symmetric_matrix_eigen_float.h"
#include "tl_sparse_symmetric_matrix_impl_eigen_float.h"

#ifdef HAVE_VIENNACL
#include "tl_sparse_general_matrix_impl_viennacl_float.h"
#include "tl_sparse_general_matrix_viennacl_float.h"
#endif  // HAVE_VIENNACL

// ---------------------------------------------------------------------------
// constructor & destructor
// ---------------------------------------------------------------------------
TlSparseGeneralMatrix_EigenFloat::TlSparseGeneralMatrix_EigenFloat(const TlMatrixObject::index_type row,
                                                                   const TlMatrixObject::index_type col) {
    this->pImpl_ = new TlSparseGeneralMatrix_ImplEigenFloat(row, col);
}

TlSparseGeneralMatrix_EigenFloat::TlSparseGeneralMatrix_EigenFloat(const TlSparseGeneralMatrix_EigenFloat& rhs) {
    this->pImpl_ = new TlSparseGeneralMatrix_ImplEigenFloat(*(dynamic_cast<const TlSparseGeneralMatrix_ImplEigenFloat*>(rhs.pImpl_)));
}

TlSparseGeneralMatrix_EigenFloat::TlSparseGeneralMatrix_EigenFloat(const TlSparseSymmetricMatrix_EigenFloat& rhs) {
    this->pImpl_ = new TlSparseGeneralMatrix_ImplEigenFloat(*(dynamic_cast<const TlSparseSymmetricMatrix_ImplEigenFloat*>(rhs.pImpl_)));
}

TlSparseGeneralMatrix_EigenFloat::TlSparseGeneralMatrix_EigenFloat(const TlSparseGeneralMatrix_ImplEigenFloat& rhs) {
    this->pImpl_ = new TlSparseGeneralMatrix_ImplEigenFloat(rhs);
}

TlSparseGeneralMatrix_EigenFloat::TlSparseGeneralMatrix_EigenFloat(const TlDenseGeneralMatrix_EigenFloat& rhs) {
    this->pImpl_ = new TlSparseGeneralMatrix_ImplEigenFloat(*(dynamic_cast<TlDenseGeneralMatrix_ImplEigenFloat*>(rhs.pImpl_)));
}

#ifdef HAVE_VIENNACL
TlSparseGeneralMatrix_EigenFloat::TlSparseGeneralMatrix_EigenFloat(const TlSparseGeneralMatrix_ViennaCLFloat& rhs) {
    this->pImpl_ = new TlSparseGeneralMatrix_ImplEigenFloat(*(dynamic_cast<TlSparseGeneralMatrix_ImplViennaCLFloat*>(rhs.pImpl_)));
}
#endif  // HAVE_VIENNACL

TlSparseGeneralMatrix_EigenFloat::~TlSparseGeneralMatrix_EigenFloat() {
    delete this->pImpl_;
    this->pImpl_ = NULL;
}

// ---------------------------------------------------------------------------
// operators
// ---------------------------------------------------------------------------
TlSparseGeneralMatrix_EigenFloat& TlSparseGeneralMatrix_EigenFloat::operator=(const TlSparseGeneralMatrix_EigenFloat& rhs) {
    if (this != &rhs) {
        delete this->pImpl_;
        this->pImpl_ = new TlSparseGeneralMatrix_ImplEigenFloat(*(dynamic_cast<TlSparseGeneralMatrix_ImplEigenFloat*>(rhs.pImpl_)));
    }

    return *this;
}

TlSparseGeneralMatrix_EigenFloat& TlSparseGeneralMatrix_EigenFloat::operator+=(const TlSparseGeneralMatrix_EigenFloat& rhs) {
    *(dynamic_cast<TlSparseGeneralMatrix_ImplEigenFloat*>(this->pImpl_)) += *(dynamic_cast<TlSparseGeneralMatrix_ImplEigenFloat*>(rhs.pImpl_));
    return *this;
}

TlSparseGeneralMatrix_EigenFloat& TlSparseGeneralMatrix_EigenFloat::operator-=(const TlSparseGeneralMatrix_EigenFloat& rhs) {
    *(dynamic_cast<TlSparseGeneralMatrix_ImplEigenFloat*>(this->pImpl_)) -= *(dynamic_cast<TlSparseGeneralMatrix_ImplEigenFloat*>(rhs.pImpl_));
    return *this;
}

TlSparseGeneralMatrix_EigenFloat& TlSparseGeneralMatrix_EigenFloat::operator*=(const double coef) {
    *(dynamic_cast<TlSparseGeneralMatrix_ImplEigenFloat*>(this->pImpl_)) *= coef;
    return *this;
}

// ----------------------------------------------------------------------------
// others
// ----------------------------------------------------------------------------
// DM(G) = DM(G) * SM(G)
TlDenseGeneralMatrix_EigenFloat operator*(const TlDenseGeneralMatrix_EigenFloat& mat1,
                                          const TlSparseGeneralMatrix_EigenFloat& mat2) {
    return TlDenseGeneralMatrix_EigenFloat(*(dynamic_cast<TlDenseGeneralMatrix_ImplEigenFloat*>(mat1.pImpl_)) * *(dynamic_cast<TlSparseGeneralMatrix_ImplEigenFloat*>(mat2.pImpl_)));
}

// DM(G) = SM(G) * DM(G)
TlDenseGeneralMatrix_EigenFloat operator*(const TlSparseGeneralMatrix_EigenFloat& mat1,
                                          const TlDenseGeneralMatrix_EigenFloat& mat2) {
    return TlDenseGeneralMatrix_EigenFloat(
        *(dynamic_cast<TlSparseGeneralMatrix_ImplEigenFloat*>(mat1.pImpl_)) *
        *(dynamic_cast<TlDenseGeneralMatrix_ImplEigenFloat*>(mat2.pImpl_)));
}

// SM(G) = SM(S) * SM(G)
TlSparseGeneralMatrix_EigenFloat operator*(const TlSparseGeneralMatrix_EigenFloat& sm1,
                                           const TlSparseGeneralMatrix_EigenFloat& sm2) {
    return TlSparseGeneralMatrix_EigenFloat(
        *(dynamic_cast<TlSparseGeneralMatrix_ImplEigenFloat*>(sm1.pImpl_)) *
        *(dynamic_cast<TlSparseGeneralMatrix_ImplEigenFloat*>(sm2.pImpl_)));
}

// DV = SM(G) * DV
TlDenseVector_EigenFloat operator*(const TlSparseGeneralMatrix_EigenFloat& mat,
                                   const TlDenseVector_EigenFloat& vtr) {
    return TlDenseVector_EigenFloat(
        *(dynamic_cast<TlSparseGeneralMatrix_ImplEigenFloat*>(mat.pImpl_)) *
        *(dynamic_cast<TlDenseVector_ImplEigenFloat*>(vtr.pImpl_)));
}
