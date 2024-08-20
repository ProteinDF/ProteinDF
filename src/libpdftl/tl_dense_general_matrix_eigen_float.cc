#ifdef HAVE_CONFIG_H
#include "config.h"  // this file created by autotools
#endif               // HAVE_CONFIG_H

#include "tl_dense_general_matrix_eigen_float.h"

#include "tl_dense_general_matrix_eigen.h"
#include "tl_dense_general_matrix_impl_eigen.h"
#include "tl_dense_general_matrix_impl_eigen_float.h"
#include "tl_dense_symmetric_matrix_eigen_float.h"
#include "tl_dense_symmetric_matrix_impl_eigen_float.h"
#include "tl_dense_vector_eigen_float.h"
#include "tl_dense_vector_impl_eigen_float.h"
#include "tl_sparse_general_matrix_eigen_float.h"

#ifdef HAVE_VIENNACL
#include "tl_dense_general_matrix_impl_viennacl_float.h"
#include "tl_dense_general_matrix_viennacl_float.h"
#endif  // HAVE_VIENNACL

TlDenseGeneralMatrix_EigenFloat::TlDenseGeneralMatrix_EigenFloat(const TlMatrixObject::index_type row, const TlMatrixObject::index_type col, double const* const pBuf) {
    this->pImpl_ = new TlDenseGeneralMatrix_ImplEigenFloat(row, col, pBuf);
}

TlDenseGeneralMatrix_EigenFloat::TlDenseGeneralMatrix_EigenFloat(const TlDenseGeneralMatrix_EigenFloat& rhs) {
    this->pImpl_ = new TlDenseGeneralMatrix_ImplEigenFloat(*(dynamic_cast<const TlDenseGeneralMatrix_ImplEigenFloat*>(rhs.pImpl_)));
}

TlDenseGeneralMatrix_EigenFloat::TlDenseGeneralMatrix_EigenFloat(const TlDenseGeneralMatrix_Eigen& rhs) {
    this->pImpl_ = new TlDenseGeneralMatrix_ImplEigenFloat(*(dynamic_cast<const TlDenseGeneralMatrix_ImplEigen*>(rhs.pImpl_)));
}

TlDenseGeneralMatrix_EigenFloat::TlDenseGeneralMatrix_EigenFloat(const TlDenseSymmetricMatrix_EigenFloat& rhs) {
    this->pImpl_ = new TlDenseGeneralMatrix_ImplEigenFloat(*(dynamic_cast<const TlDenseSymmetricMatrix_ImplEigenFloat*>(rhs.pImpl_)));
}

TlDenseGeneralMatrix_EigenFloat::TlDenseGeneralMatrix_EigenFloat(const TlDenseGeneralMatrix_ImplEigenFloat& rhs) {
    this->pImpl_ = new TlDenseGeneralMatrix_ImplEigenFloat(rhs);
}

TlDenseGeneralMatrix_EigenFloat::TlDenseGeneralMatrix_EigenFloat(const TlSparseGeneralMatrix_EigenFloat& sm) {
    this->pImpl_ = new TlDenseGeneralMatrix_ImplEigenFloat(*(dynamic_cast<TlSparseGeneralMatrix_ImplEigenFloat*>(sm.pImpl_)));
}

#ifdef HAVE_VIENNACL
TlDenseGeneralMatrix_EigenFloat::TlDenseGeneralMatrix_EigenFloat(const TlDenseGeneralMatrix_ViennaCLFloat& rhs) {
    this->pImpl_ = new TlDenseGeneralMatrix_ImplEigenFloat(*(dynamic_cast<TlDenseGeneralMatrix_ImplViennaCLFloat*>(rhs.pImpl_)));
}
#endif  // HAVE_VIENNACL

TlDenseGeneralMatrix_EigenFloat::~TlDenseGeneralMatrix_EigenFloat() {
    delete this->pImpl_;
    this->pImpl_ = NULL;
}

TlDenseGeneralMatrix_EigenFloat::operator std::vector<double>() const {
    return *(dynamic_cast<TlDenseGeneralMatrix_ImplEigenFloat*>(this->pImpl_));
}

// ---------------------------------------------------------------------------
// operators
// ---------------------------------------------------------------------------
TlDenseGeneralMatrix_EigenFloat& TlDenseGeneralMatrix_EigenFloat::operator=(
    const TlDenseGeneralMatrix_EigenFloat& rhs) {
    delete this->pImpl_;
    this->pImpl_ = new TlDenseGeneralMatrix_ImplEigenFloat(
        *(dynamic_cast<TlDenseGeneralMatrix_ImplEigenFloat*>(rhs.pImpl_)));

    return *this;
}

const TlDenseGeneralMatrix_EigenFloat TlDenseGeneralMatrix_EigenFloat::operator+(
    const TlDenseGeneralMatrix_EigenFloat& rhs) const {
    TlDenseGeneralMatrix_EigenFloat answer = *this;
    answer += rhs;
    return answer;
}

const TlDenseGeneralMatrix_EigenFloat TlDenseGeneralMatrix_EigenFloat::operator-(
    const TlDenseGeneralMatrix_EigenFloat& rhs) const {
    TlDenseGeneralMatrix_EigenFloat answer = *this;
    answer -= rhs;
    return answer;
}

const TlDenseGeneralMatrix_EigenFloat TlDenseGeneralMatrix_EigenFloat::operator*(
    const TlDenseGeneralMatrix_EigenFloat& rhs) const {
    TlDenseGeneralMatrix_EigenFloat answer = *this;
    answer *= rhs;
    return answer;
}

TlDenseGeneralMatrix_EigenFloat& TlDenseGeneralMatrix_EigenFloat::operator+=(
    const TlDenseGeneralMatrix_EigenFloat& rhs) {
    *(dynamic_cast<TlDenseGeneralMatrix_ImplEigenFloat*>(this->pImpl_)) +=
        *(dynamic_cast<TlDenseGeneralMatrix_ImplEigenFloat*>(rhs.pImpl_));

    return *this;
}

TlDenseGeneralMatrix_EigenFloat& TlDenseGeneralMatrix_EigenFloat::operator-=(
    const TlDenseGeneralMatrix_EigenFloat& rhs) {
    *(dynamic_cast<TlDenseGeneralMatrix_ImplEigenFloat*>(this->pImpl_)) -=
        *(dynamic_cast<TlDenseGeneralMatrix_ImplEigenFloat*>(rhs.pImpl_));

    return *this;
}

TlDenseGeneralMatrix_EigenFloat& TlDenseGeneralMatrix_EigenFloat::operator*=(
    const double coef) {
    *(dynamic_cast<TlDenseGeneralMatrix_ImplEigenFloat*>(this->pImpl_)) *= coef;

    return *this;
}

TlDenseGeneralMatrix_EigenFloat& TlDenseGeneralMatrix_EigenFloat::operator/=(
    const double coef) {
    *(dynamic_cast<TlDenseGeneralMatrix_ImplEigenFloat*>(this->pImpl_)) /= coef;

    return *this;
}

TlDenseGeneralMatrix_EigenFloat& TlDenseGeneralMatrix_EigenFloat::operator*=(
    const TlDenseGeneralMatrix_EigenFloat& rhs) {
    *(dynamic_cast<TlDenseGeneralMatrix_ImplEigenFloat*>(this->pImpl_)) *=
        *(dynamic_cast<TlDenseGeneralMatrix_ImplEigenFloat*>(rhs.pImpl_));

    return *this;
}

// ---------------------------------------------------------------------------
// operations
// ---------------------------------------------------------------------------
double TlDenseGeneralMatrix_EigenFloat::sum() const {
    return this->pImpl_->sum();
}

double TlDenseGeneralMatrix_EigenFloat::getRMS() const {
    return this->pImpl_->getRMS();
}

TlDenseGeneralMatrix_EigenFloat TlDenseGeneralMatrix_EigenFloat::dot(
    const TlDenseGeneralMatrix_EigenFloat& rhs) const {
    TlDenseGeneralMatrix_EigenFloat answer(
        dynamic_cast<TlDenseGeneralMatrix_ImplEigenFloat*>(this->pImpl_)
            ->dot(*(dynamic_cast<const TlDenseGeneralMatrix_ImplEigenFloat*>(
                rhs.pImpl_))));

    return answer;
}

const TlDenseGeneralMatrix_EigenFloat& TlDenseGeneralMatrix_EigenFloat::dotInPlace(
    const TlDenseGeneralMatrix_EigenFloat& rhs) {
    dynamic_cast<TlDenseGeneralMatrix_ImplEigenFloat*>(this->pImpl_)
        ->dotInPlace(
            *(dynamic_cast<const TlDenseGeneralMatrix_ImplEigenFloat*>(rhs.pImpl_)));
    return *this;
}

TlDenseGeneralMatrix_EigenFloat TlDenseGeneralMatrix_EigenFloat::transpose() const {
    return TlDenseGeneralMatrix_EigenFloat(
        dynamic_cast<TlDenseGeneralMatrix_ImplEigenFloat*>(this->pImpl_)
            ->transpose());
}

TlDenseGeneralMatrix_EigenFloat TlDenseGeneralMatrix_EigenFloat::inverse() const {
    return TlDenseGeneralMatrix_EigenFloat(
        dynamic_cast<TlDenseGeneralMatrix_ImplEigenFloat*>(this->pImpl_)->inverse());
}

// ---------------------------------------------------------------------------
// I/O
// ---------------------------------------------------------------------------
TlMatrixObject::size_type TlDenseGeneralMatrix_EigenFloat::getNumOfElements() const {
    return dynamic_cast<TlDenseGeneralMatrix_ImplEigenFloat*>(this->pImpl_)->getNumOfElements();
}

float* TlDenseGeneralMatrix_EigenFloat::data() {
    return dynamic_cast<TlDenseGeneralMatrix_ImplEigenFloat*>(this->pImpl_)->data();
}

const float* TlDenseGeneralMatrix_EigenFloat::data() const {
    return dynamic_cast<TlDenseGeneralMatrix_ImplEigenFloat*>(this->pImpl_)->data();
}

// ---------------------------------------------------------------------------
// friend functions
// ---------------------------------------------------------------------------
TlDenseVector_EigenFloat operator*(const TlDenseGeneralMatrix_EigenFloat& rhs1,
                              const TlDenseVector_EigenFloat& rhs2) {
    return TlDenseVector_EigenFloat(
        *(dynamic_cast<TlDenseGeneralMatrix_ImplEigenFloat*>(rhs1.pImpl_)) *
        *(dynamic_cast<TlDenseVector_ImplEigenFloat*>(rhs2.pImpl_)));
}

TlDenseVector_EigenFloat operator*(const TlDenseVector_EigenFloat& rhs1,
                              const TlDenseGeneralMatrix_EigenFloat& rhs2) {
    return TlDenseVector_EigenFloat(
        *(dynamic_cast<TlDenseVector_ImplEigenFloat*>(rhs1.pImpl_)) *
        *(dynamic_cast<TlDenseGeneralMatrix_ImplEigenFloat*>(rhs2.pImpl_)));
}

TlDenseGeneralMatrix_EigenFloat operator*(const double coef,
                                          const TlDenseGeneralMatrix_EigenFloat& DM) {
    TlDenseGeneralMatrix_EigenFloat answer = DM;
    answer *= coef;
    return answer;
}

TlDenseGeneralMatrix_EigenFloat operator*(const TlDenseGeneralMatrix_EigenFloat& DM,
                                          const double coef) {
    TlDenseGeneralMatrix_EigenFloat answer = DM;
    answer *= coef;
    return answer;
}
