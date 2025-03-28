#ifdef HAVE_CONFIG_H
#include "config.h"  // this file created by autotools
#endif               // HAVE_CONFIG_H

#include "tl_dense_general_matrix_eigen.h"
#include "tl_dense_general_matrix_eigen_float.h"
#include "tl_dense_general_matrix_impl_eigen.h"
#include "tl_dense_general_matrix_impl_eigen_float.h"
#include "tl_dense_symmetric_matrix_eigen.h"
#include "tl_dense_symmetric_matrix_impl_eigen.h"
#include "tl_dense_vector_eigen.h"
#include "tl_dense_vector_impl_eigen.h"
#include "tl_sparse_general_matrix_eigen.h"

#ifdef HAVE_VIENNACL
#include "tl_dense_general_matrix_impl_viennacl.h"
#include "tl_dense_general_matrix_viennacl.h"
#endif  // HAVE_VIENNACL

TlDenseGeneralMatrix_Eigen::TlDenseGeneralMatrix_Eigen(const TlMatrixObject::index_type row, const TlMatrixObject::index_type col,
                                                       double const* const pBuf) {
    this->pImpl_ = new TlDenseGeneralMatrix_ImplEigen(row, col, pBuf);
}

TlDenseGeneralMatrix_Eigen::TlDenseGeneralMatrix_Eigen(const TlDenseGeneralMatrix_Eigen& rhs) {
    this->pImpl_ = new TlDenseGeneralMatrix_ImplEigen(*(dynamic_cast<const TlDenseGeneralMatrix_ImplEigen*>(rhs.pImpl_)));
}

TlDenseGeneralMatrix_Eigen::TlDenseGeneralMatrix_Eigen(const TlDenseGeneralMatrix_EigenFloat& rhs) {
    this->pImpl_ = new TlDenseGeneralMatrix_ImplEigen(*(dynamic_cast<const TlDenseGeneralMatrix_ImplEigenFloat*>(rhs.pImpl_)));
}

TlDenseGeneralMatrix_Eigen::TlDenseGeneralMatrix_Eigen(const TlDenseSymmetricMatrix_Eigen& rhs) {
    this->pImpl_ = new TlDenseGeneralMatrix_ImplEigen(*(dynamic_cast<const TlDenseSymmetricMatrix_ImplEigen*>(rhs.pImpl_)));
}

TlDenseGeneralMatrix_Eigen::TlDenseGeneralMatrix_Eigen(const TlDenseGeneralMatrix_ImplEigen& rhs) {
    this->pImpl_ = new TlDenseGeneralMatrix_ImplEigen(rhs);
}

TlDenseGeneralMatrix_Eigen::TlDenseGeneralMatrix_Eigen(
    const TlSparseGeneralMatrix_Eigen& sm) {
    this->pImpl_ = new TlDenseGeneralMatrix_ImplEigen(
        *(dynamic_cast<TlSparseGeneralMatrix_ImplEigen*>(sm.pImpl_)));
}

#ifdef HAVE_VIENNACL
TlDenseGeneralMatrix_Eigen::TlDenseGeneralMatrix_Eigen(
    const TlDenseGeneralMatrix_ViennaCL& rhs) {
    this->pImpl_ = new TlDenseGeneralMatrix_ImplEigen(
        *(dynamic_cast<TlDenseGeneralMatrix_ImplViennaCL*>(rhs.pImpl_)));
}
#endif  // HAVE_VIENNACL

TlDenseGeneralMatrix_Eigen::~TlDenseGeneralMatrix_Eigen() {
    delete this->pImpl_;
    this->pImpl_ = NULL;
}

TlDenseGeneralMatrix_Eigen::operator std::vector<double>() const {
    return *(dynamic_cast<TlDenseGeneralMatrix_ImplEigen*>(this->pImpl_));
}

// ---------------------------------------------------------------------------
// operators
// ---------------------------------------------------------------------------
TlDenseGeneralMatrix_Eigen& TlDenseGeneralMatrix_Eigen::operator=(
    const TlDenseGeneralMatrix_Eigen& rhs) {
    delete this->pImpl_;
    this->pImpl_ = new TlDenseGeneralMatrix_ImplEigen(
        *(dynamic_cast<TlDenseGeneralMatrix_ImplEigen*>(rhs.pImpl_)));

    return *this;
}

const TlDenseGeneralMatrix_Eigen TlDenseGeneralMatrix_Eigen::operator+(
    const TlDenseGeneralMatrix_Eigen& rhs) const {
    TlDenseGeneralMatrix_Eigen answer = *this;
    answer += rhs;
    return answer;
}

const TlDenseGeneralMatrix_Eigen TlDenseGeneralMatrix_Eigen::operator-(
    const TlDenseGeneralMatrix_Eigen& rhs) const {
    TlDenseGeneralMatrix_Eigen answer = *this;
    answer -= rhs;
    return answer;
}

const TlDenseGeneralMatrix_Eigen TlDenseGeneralMatrix_Eigen::operator*(
    const TlDenseGeneralMatrix_Eigen& rhs) const {
    TlDenseGeneralMatrix_Eigen answer = *this;
    answer *= rhs;
    return answer;
}

TlDenseGeneralMatrix_Eigen& TlDenseGeneralMatrix_Eigen::operator+=(
    const TlDenseGeneralMatrix_Eigen& rhs) {
    *(dynamic_cast<TlDenseGeneralMatrix_ImplEigen*>(this->pImpl_)) +=
        *(dynamic_cast<TlDenseGeneralMatrix_ImplEigen*>(rhs.pImpl_));

    return *this;
}

TlDenseGeneralMatrix_Eigen& TlDenseGeneralMatrix_Eigen::operator-=(
    const TlDenseGeneralMatrix_Eigen& rhs) {
    *(dynamic_cast<TlDenseGeneralMatrix_ImplEigen*>(this->pImpl_)) -=
        *(dynamic_cast<TlDenseGeneralMatrix_ImplEigen*>(rhs.pImpl_));

    return *this;
}

TlDenseGeneralMatrix_Eigen& TlDenseGeneralMatrix_Eigen::operator*=(
    const double coef) {
    *(dynamic_cast<TlDenseGeneralMatrix_ImplEigen*>(this->pImpl_)) *= coef;

    return *this;
}

TlDenseGeneralMatrix_Eigen& TlDenseGeneralMatrix_Eigen::operator/=(
    const double coef) {
    *(dynamic_cast<TlDenseGeneralMatrix_ImplEigen*>(this->pImpl_)) /= coef;

    return *this;
}

TlDenseGeneralMatrix_Eigen& TlDenseGeneralMatrix_Eigen::operator*=(
    const TlDenseGeneralMatrix_Eigen& rhs) {
    *(dynamic_cast<TlDenseGeneralMatrix_ImplEigen*>(this->pImpl_)) *=
        *(dynamic_cast<TlDenseGeneralMatrix_ImplEigen*>(rhs.pImpl_));

    return *this;
}

// ---------------------------------------------------------------------------
// operations
// ---------------------------------------------------------------------------
double TlDenseGeneralMatrix_Eigen::sum() const {
    return this->pImpl_->sum();
}

double TlDenseGeneralMatrix_Eigen::getRMS() const {
    return this->pImpl_->getRMS();
}

TlDenseGeneralMatrix_Eigen TlDenseGeneralMatrix_Eigen::dot(
    const TlDenseGeneralMatrix_Eigen& rhs) const {
    TlDenseGeneralMatrix_Eigen answer(
        dynamic_cast<TlDenseGeneralMatrix_ImplEigen*>(this->pImpl_)
            ->dot(*(dynamic_cast<const TlDenseGeneralMatrix_ImplEigen*>(
                rhs.pImpl_))));

    return answer;
}

const TlDenseGeneralMatrix_Eigen& TlDenseGeneralMatrix_Eigen::dotInPlace(
    const TlDenseGeneralMatrix_Eigen& rhs) {
    dynamic_cast<TlDenseGeneralMatrix_ImplEigen*>(this->pImpl_)
        ->dotInPlace(
            *(dynamic_cast<const TlDenseGeneralMatrix_ImplEigen*>(rhs.pImpl_)));
    return *this;
}

TlDenseGeneralMatrix_Eigen TlDenseGeneralMatrix_Eigen::transpose() const {
    return TlDenseGeneralMatrix_Eigen(
        dynamic_cast<TlDenseGeneralMatrix_ImplEigen*>(this->pImpl_)
            ->transpose());
}

TlDenseGeneralMatrix_Eigen TlDenseGeneralMatrix_Eigen::inverse() const {
    return TlDenseGeneralMatrix_Eigen(
        dynamic_cast<TlDenseGeneralMatrix_ImplEigen*>(this->pImpl_)->inverse());
}

// ---------------------------------------------------------------------------
// I/O
// ---------------------------------------------------------------------------
TlMatrixObject::size_type TlDenseGeneralMatrix_Eigen::getNumOfElements() const {
    return dynamic_cast<TlDenseGeneralMatrix_ImplEigen*>(this->pImpl_)->getNumOfElements();
}

double* TlDenseGeneralMatrix_Eigen::data() {
    return dynamic_cast<TlDenseGeneralMatrix_ImplEigen*>(this->pImpl_)->data();
}

const double* TlDenseGeneralMatrix_Eigen::data() const {
    return dynamic_cast<TlDenseGeneralMatrix_ImplEigen*>(this->pImpl_)->data();
}

// ---------------------------------------------------------------------------
// friend functions
// ---------------------------------------------------------------------------
TlDenseVector_Eigen operator*(const TlDenseGeneralMatrix_Eigen& rhs1,
                              const TlDenseVector_Eigen& rhs2) {
    return TlDenseVector_Eigen(
        *(dynamic_cast<TlDenseGeneralMatrix_ImplEigen*>(rhs1.pImpl_)) *
        *(dynamic_cast<TlDenseVector_ImplEigen*>(rhs2.pImpl_)));
}

TlDenseVector_Eigen operator*(const TlDenseVector_Eigen& rhs1,
                              const TlDenseGeneralMatrix_Eigen& rhs2) {
    return TlDenseVector_Eigen(
        *(dynamic_cast<TlDenseVector_ImplEigen*>(rhs1.pImpl_)) *
        *(dynamic_cast<TlDenseGeneralMatrix_ImplEigen*>(rhs2.pImpl_)));
}

TlDenseGeneralMatrix_Eigen operator*(const double coef,
                                     const TlDenseGeneralMatrix_Eigen& DM) {
    TlDenseGeneralMatrix_Eigen answer = DM;
    answer *= coef;
    return answer;
}

TlDenseGeneralMatrix_Eigen operator*(const TlDenseGeneralMatrix_Eigen& DM,
                                     const double coef) {
    TlDenseGeneralMatrix_Eigen answer = DM;
    answer *= coef;
    return answer;
}
