#ifdef HAVE_CONFIG_H
#include "config.h"  // this file created by autotools
#endif               // HAVE_CONFIG_H

#include "tl_dense_general_matrix_eigen_float.h"
#include "tl_dense_symmetric_matrix_eigen.h"
#include "tl_dense_symmetric_matrix_impl_eigen.h"
#include "tl_dense_symmetric_matrix_eigen_float.h"
#include "tl_dense_symmetric_matrix_impl_eigen_float.h"
#include "tl_dense_vector_eigen_float.h"
#include "tl_dense_vector_impl_eigen_float.h"
#include "tl_sparse_symmetric_matrix_eigen_float.h"
#include "tl_sparse_symmetric_matrix_impl_eigen_float.h"

#ifdef HAVE_VIENNACL
#include "tl_dense_symmetric_matrix_impl_viennacl_float.h"
#include "tl_dense_symmetric_matrix_viennacl_float.h"
#endif  // HAVE_VIENNACL

TlDenseSymmetricMatrix_EigenFloat::TlDenseSymmetricMatrix_EigenFloat(const TlMatrixObject::index_type dim, double const* const pBuf) {
    this->pImpl_ = new TlDenseSymmetricMatrix_ImplEigenFloat(dim, pBuf);
}

TlDenseSymmetricMatrix_EigenFloat::TlDenseSymmetricMatrix_EigenFloat(const TlDenseSymmetricMatrix_EigenFloat& rhs) {
    this->pImpl_ = new TlDenseSymmetricMatrix_ImplEigenFloat(*(dynamic_cast<const TlDenseSymmetricMatrix_ImplEigenFloat*>(rhs.pImpl_)));
}

TlDenseSymmetricMatrix_EigenFloat::TlDenseSymmetricMatrix_EigenFloat(const TlDenseSymmetricMatrix_Eigen& rhs) {
    this->pImpl_ = new TlDenseSymmetricMatrix_ImplEigenFloat(*(dynamic_cast<const TlDenseSymmetricMatrix_ImplEigen*>(rhs.pImpl_)));
}

TlDenseSymmetricMatrix_EigenFloat::TlDenseSymmetricMatrix_EigenFloat(const TlDenseGeneralMatrix_EigenFloat& rhs) {
    this->pImpl_ = new TlDenseSymmetricMatrix_ImplEigenFloat(*(dynamic_cast<const TlDenseGeneralMatrix_ImplEigenFloat*>(rhs.pImpl_)));
}

TlDenseSymmetricMatrix_EigenFloat::TlDenseSymmetricMatrix_EigenFloat(const TlSparseSymmetricMatrix_EigenFloat& sm) {
    this->pImpl_ = new TlDenseSymmetricMatrix_ImplEigenFloat(*(dynamic_cast<TlSparseSymmetricMatrix_ImplEigenFloat*>(sm.pImpl_)));
}

#ifdef HAVE_VIENNACL
TlDenseSymmetricMatrix_EigenFloat::TlDenseSymmetricMatrix_EigenFloat(const TlDenseSymmetricMatrix_ViennaCLFloat& rhs) {
    this->pImpl_ = new TlDenseSymmetricMatrix_ImplEigenFloat(*(dynamic_cast<TlDenseSymmetricMatrix_ImplViennaCLFloat*>(rhs.pImpl_)));
}
#endif  // HAVE_VIENNACL

TlDenseSymmetricMatrix_EigenFloat::~TlDenseSymmetricMatrix_EigenFloat() {
    delete this->pImpl_;
    this->pImpl_ = NULL;
}

TlDenseSymmetricMatrix_EigenFloat::operator std::vector<double>() const {
    return *(dynamic_cast<TlDenseSymmetricMatrix_ImplEigenFloat*>(this->pImpl_));
}
// ---------------------------------------------------------------------------
// operators
// ---------------------------------------------------------------------------
TlDenseSymmetricMatrix_EigenFloat& TlDenseSymmetricMatrix_EigenFloat::operator=(
    const TlDenseSymmetricMatrix_EigenFloat& rhs) {
    if (this != &rhs) {
        delete this->pImpl_;
        this->pImpl_ = new TlDenseSymmetricMatrix_ImplEigenFloat(
            *(dynamic_cast<TlDenseSymmetricMatrix_ImplEigenFloat*>(rhs.pImpl_)));
    }

    return *this;
}

const TlDenseSymmetricMatrix_EigenFloat TlDenseSymmetricMatrix_EigenFloat::operator+(
    const TlDenseSymmetricMatrix_EigenFloat& rhs) const {
    TlDenseSymmetricMatrix_EigenFloat answer = *this;
    answer += rhs;
    return answer;
}

const TlDenseSymmetricMatrix_EigenFloat TlDenseSymmetricMatrix_EigenFloat::operator-(
    const TlDenseSymmetricMatrix_EigenFloat& rhs) const {
    TlDenseSymmetricMatrix_EigenFloat answer = *this;
    answer -= rhs;
    return answer;
}

const TlDenseSymmetricMatrix_EigenFloat TlDenseSymmetricMatrix_EigenFloat::operator*(
    const TlDenseSymmetricMatrix_EigenFloat& rhs) const {
    TlDenseSymmetricMatrix_EigenFloat answer = *this;
    answer *= rhs;
    return answer;
}

TlDenseSymmetricMatrix_EigenFloat& TlDenseSymmetricMatrix_EigenFloat::operator+=(
    const TlDenseSymmetricMatrix_EigenFloat& rhs) {
    *(dynamic_cast<TlDenseSymmetricMatrix_ImplEigenFloat*>(this->pImpl_)) +=
        *(dynamic_cast<TlDenseSymmetricMatrix_ImplEigenFloat*>(rhs.pImpl_));

    return *this;
}

TlDenseSymmetricMatrix_EigenFloat& TlDenseSymmetricMatrix_EigenFloat::operator-=(
    const TlDenseSymmetricMatrix_EigenFloat& rhs) {
    *(dynamic_cast<TlDenseSymmetricMatrix_ImplEigenFloat*>(this->pImpl_)) -=
        *(dynamic_cast<TlDenseSymmetricMatrix_ImplEigenFloat*>(rhs.pImpl_));

    return *this;
}

TlDenseSymmetricMatrix_EigenFloat& TlDenseSymmetricMatrix_EigenFloat::operator*=(
    const double coef) {
    *(dynamic_cast<TlDenseSymmetricMatrix_ImplEigenFloat*>(this->pImpl_)) *= coef;

    return *this;
}

TlDenseSymmetricMatrix_EigenFloat& TlDenseSymmetricMatrix_EigenFloat::operator/=(
    const double coef) {
    *(dynamic_cast<TlDenseSymmetricMatrix_ImplEigenFloat*>(this->pImpl_)) /= coef;

    return *this;
}

TlDenseSymmetricMatrix_EigenFloat& TlDenseSymmetricMatrix_EigenFloat::operator*=(
    const TlDenseSymmetricMatrix_EigenFloat& rhs) {
    *(dynamic_cast<TlDenseSymmetricMatrix_ImplEigenFloat*>(this->pImpl_)) *=
        *(dynamic_cast<TlDenseSymmetricMatrix_ImplEigenFloat*>(rhs.pImpl_));

    return *this;
}

// ---------------------------------------------------------------------------
// operations
// ---------------------------------------------------------------------------
const TlDenseSymmetricMatrix_EigenFloat& TlDenseSymmetricMatrix_EigenFloat::dotInPlace(
    const TlDenseSymmetricMatrix_EigenFloat& rhs) {
    dynamic_cast<TlDenseSymmetricMatrix_ImplEigenFloat*>(this->pImpl_)
        ->dotInPlace(
            *dynamic_cast<TlDenseSymmetricMatrix_ImplEigenFloat*>(rhs.pImpl_));

    return *this;
}

bool TlDenseSymmetricMatrix_EigenFloat::eig(TlDenseVector_EigenFloat* pEigVal, TlDenseGeneralMatrix_EigenFloat* pEigVec) const {
    TlDenseVector_ImplEigenFloat* pImpl_eigval = dynamic_cast<TlDenseVector_ImplEigenFloat*>(pEigVal->pImpl_);
    TlDenseGeneralMatrix_ImplEigenFloat* pImpl_eigvec = dynamic_cast<TlDenseGeneralMatrix_ImplEigenFloat*>(pEigVec->pImpl_);

    const bool answer = dynamic_cast<TlDenseSymmetricMatrix_ImplEigenFloat*>(this->pImpl_)->eig(pImpl_eigval, pImpl_eigvec);
    return answer;
}

TlDenseSymmetricMatrix_EigenFloat TlDenseSymmetricMatrix_EigenFloat::inverse() const {
    TlDenseSymmetricMatrix_ImplEigenFloat tmp = dynamic_cast<const TlDenseSymmetricMatrix_ImplEigenFloat*>(this->pImpl_)->inverse();
    TlDenseSymmetricMatrix_EigenFloat answer(tmp);
    return answer;
}

// ---------------------------------------------------------------------------
// I/O
// ---------------------------------------------------------------------------
TlMatrixObject::size_type TlDenseSymmetricMatrix_EigenFloat::getNumOfElements() const {
    return dynamic_cast<TlDenseSymmetricMatrix_ImplEigenFloat*>(this->pImpl_)->getNumOfElements();
}

float* TlDenseSymmetricMatrix_EigenFloat::data() {
    return dynamic_cast<TlDenseSymmetricMatrix_ImplEigenFloat*>(this->pImpl_)->data();
}

const float* TlDenseSymmetricMatrix_EigenFloat::data() const {
    return dynamic_cast<TlDenseSymmetricMatrix_ImplEigenFloat*>(this->pImpl_)->data();
}

// ---------------------------------------------------------------------------
// friend functions
// ---------------------------------------------------------------------------
TlDenseGeneralMatrix_EigenFloat operator*(const TlDenseGeneralMatrix_EigenFloat& mat1,
                                          const TlDenseSymmetricMatrix_EigenFloat& mat2) {
    const TlDenseGeneralMatrix_ImplEigenFloat tmp = *(dynamic_cast<TlDenseGeneralMatrix_ImplEigenFloat*>(mat1.pImpl_)) *
                                                    *(dynamic_cast<TlDenseSymmetricMatrix_ImplEigenFloat*>(mat2.pImpl_));
    TlDenseGeneralMatrix_EigenFloat answer(tmp);
    return answer;
}

TlDenseGeneralMatrix_EigenFloat operator*(const TlDenseSymmetricMatrix_EigenFloat& mat1,
                                          const TlDenseGeneralMatrix_EigenFloat& mat2) {
    const TlDenseGeneralMatrix_ImplEigenFloat tmp = *(dynamic_cast<TlDenseSymmetricMatrix_ImplEigenFloat*>(mat1.pImpl_)) *
                                                    *(dynamic_cast<TlDenseGeneralMatrix_ImplEigenFloat*>(mat2.pImpl_));
    TlDenseGeneralMatrix_EigenFloat answer(tmp);
    return answer;
}

TlDenseVector_EigenFloat operator*(const TlDenseSymmetricMatrix_EigenFloat& dms1,
                                   const TlDenseVector_EigenFloat& dv) {
    return TlDenseVector_EigenFloat(*(dynamic_cast<TlDenseSymmetricMatrix_ImplEigenFloat*>(dms1.pImpl_)) *
                                    *(dynamic_cast<TlDenseVector_ImplEigenFloat*>(dv.pImpl_)));
}

TlDenseVector_EigenFloat operator*(const TlDenseVector_EigenFloat& dv,
                                   const TlDenseSymmetricMatrix_EigenFloat& dms1) {
    return TlDenseVector_EigenFloat(*(dynamic_cast<TlDenseVector_ImplEigenFloat*>(dv.pImpl_)) *
                                    *(dynamic_cast<TlDenseSymmetricMatrix_ImplEigenFloat*>(dms1.pImpl_)));
}

TlDenseSymmetricMatrix_EigenFloat operator*(const double coef,
                                            const TlDenseSymmetricMatrix_EigenFloat& DM) {
    TlDenseSymmetricMatrix_EigenFloat answer = DM;
    answer *= coef;
    return answer;
}

TlDenseSymmetricMatrix_EigenFloat operator*(const TlDenseSymmetricMatrix_EigenFloat& DM,
                                            const double coef) {
    TlDenseSymmetricMatrix_EigenFloat answer = DM;
    answer *= coef;
    return answer;
}
