#ifdef HAVE_CONFIG_H
#include "config.h"
#endif  // HAVE_CONFIG_H

#include "tl_dense_general_matrix_impl_viennacl.h"
#include "tl_dense_general_matrix_impl_viennacl_float.h"
#include "tl_dense_general_matrix_viennacl.h"
#include "tl_dense_general_matrix_viennacl_float.h"
#include "tl_dense_symmetric_matrix_impl_viennacl_float.h"
#include "tl_dense_symmetric_matrix_viennacl_float.h"
#include "tl_dense_vector_impl_viennacl_float.h"
#include "tl_dense_vector_viennacl_float.h"
#include "tl_sparse_general_matrix_impl_viennacl_float.h"
#include "tl_sparse_general_matrix_viennacl_float.h"

#ifdef HAVE_EIGEN
#include "tl_dense_general_matrix_eigen_float.h"
#include "tl_dense_general_matrix_impl_eigen_float.h"
#endif  // HAVE_EIGEN

TlDenseGeneralMatrix_ViennaCLFloat::TlDenseGeneralMatrix_ViennaCLFloat(const TlMatrixObject::index_type row, const TlMatrixObject::index_type col, double const* const pBuf) {
    this->pImpl_ = new TlDenseGeneralMatrix_ImplViennaCLFloat(row, col, pBuf);
}

TlDenseGeneralMatrix_ViennaCLFloat::TlDenseGeneralMatrix_ViennaCLFloat(const TlDenseGeneralMatrix_ViennaCLFloat& rhs) {
    this->pImpl_ = new TlDenseGeneralMatrix_ImplViennaCLFloat(*(dynamic_cast<const TlDenseGeneralMatrix_ImplViennaCLFloat*>(rhs.pImpl_)));
}

TlDenseGeneralMatrix_ViennaCLFloat::TlDenseGeneralMatrix_ViennaCLFloat(const TlDenseGeneralMatrix_ViennaCL& rhs) {
    this->pImpl_ = new TlDenseGeneralMatrix_ImplViennaCLFloat(*(dynamic_cast<const TlDenseGeneralMatrix_ImplViennaCL*>(rhs.pImpl_)));
}

TlDenseGeneralMatrix_ViennaCLFloat::TlDenseGeneralMatrix_ViennaCLFloat(const TlDenseSymmetricMatrix_ViennaCLFloat& rhs) {
    this->pImpl_ = new TlDenseGeneralMatrix_ImplViennaCLFloat(*(dynamic_cast<const TlDenseSymmetricMatrix_ImplViennaCLFloat*>(rhs.pImpl_)));
}

TlDenseGeneralMatrix_ViennaCLFloat::TlDenseGeneralMatrix_ViennaCLFloat(const TlDenseGeneralMatrix_ImplViennaCLFloat& rhs) {
    this->pImpl_ = new TlDenseGeneralMatrix_ImplViennaCLFloat(rhs);
}

TlDenseGeneralMatrix_ViennaCLFloat::TlDenseGeneralMatrix_ViennaCLFloat(const TlSparseGeneralMatrix_ViennaCLFloat& rhs) {
    this->pImpl_ = new TlDenseGeneralMatrix_ImplViennaCLFloat(*(dynamic_cast<TlSparseGeneralMatrix_ImplViennaCLFloat*>(rhs.pImpl_)));
}

#ifdef HAVE_EIGEN
TlDenseGeneralMatrix_ViennaCLFloat::TlDenseGeneralMatrix_ViennaCLFloat(const TlDenseGeneralMatrix_EigenFloat& rhs) {
    this->pImpl_ = new TlDenseGeneralMatrix_ImplViennaCLFloat(*dynamic_cast<const TlDenseGeneralMatrix_ImplEigenFloat*>(rhs.pImpl_));
}
#endif  // HAVE_EIGEN

TlDenseGeneralMatrix_ViennaCLFloat::~TlDenseGeneralMatrix_ViennaCLFloat() {
    delete this->pImpl_;
    this->pImpl_ = NULL;
}

TlDenseGeneralMatrix_ViennaCLFloat::operator std::vector<double>() const {
    return *(dynamic_cast<TlDenseGeneralMatrix_ImplViennaCLFloat*>(this->pImpl_));
}

// ---------------------------------------------------------------------------
// operators
// ---------------------------------------------------------------------------
TlDenseGeneralMatrix_ViennaCLFloat& TlDenseGeneralMatrix_ViennaCLFloat::operator=(const TlDenseGeneralMatrix_ViennaCLFloat& rhs) {
    delete this->pImpl_;
    this->pImpl_ = new TlDenseGeneralMatrix_ImplViennaCLFloat(*(dynamic_cast<TlDenseGeneralMatrix_ImplViennaCLFloat*>(rhs.pImpl_)));

    return *this;
}

// TlDenseGeneralMatrix_ViennaCLFloat& TlDenseGeneralMatrix_ViennaCLFloat::operator=(const
// TlDenseGeneralMatrix_EigenFloat_Old& rhs) {
//     *(this->pImpl_) = *(dynamic_cast<TlDenselMatrix_ImplEigen>(rhs.pImpl_));
// }

const TlDenseGeneralMatrix_ViennaCLFloat TlDenseGeneralMatrix_ViennaCLFloat::operator+(const TlDenseGeneralMatrix_ViennaCLFloat& rhs) const {
    TlDenseGeneralMatrix_ViennaCLFloat answer = *this;
    answer += rhs;
    return answer;
}

const TlDenseGeneralMatrix_ViennaCLFloat TlDenseGeneralMatrix_ViennaCLFloat::operator-(const TlDenseGeneralMatrix_ViennaCLFloat& rhs) const {
    TlDenseGeneralMatrix_ViennaCLFloat answer = *this;
    answer -= rhs;
    return answer;
}

const TlDenseGeneralMatrix_ViennaCLFloat TlDenseGeneralMatrix_ViennaCLFloat::operator*(const TlDenseGeneralMatrix_ViennaCLFloat& rhs) const {
    TlDenseGeneralMatrix_ViennaCLFloat answer = *this;
    answer *= rhs;
    return answer;
}

TlDenseGeneralMatrix_ViennaCLFloat& TlDenseGeneralMatrix_ViennaCLFloat::operator+=(const TlDenseGeneralMatrix_ViennaCLFloat& rhs) {
    *(dynamic_cast<TlDenseGeneralMatrix_ImplViennaCLFloat*>(this->pImpl_)) +=
        *(dynamic_cast<TlDenseGeneralMatrix_ImplViennaCLFloat*>(rhs.pImpl_));

    return *this;
}

TlDenseGeneralMatrix_ViennaCLFloat& TlDenseGeneralMatrix_ViennaCLFloat::operator-=(const TlDenseGeneralMatrix_ViennaCLFloat& rhs) {
    *(dynamic_cast<TlDenseGeneralMatrix_ImplViennaCLFloat*>(this->pImpl_)) -=
        *(dynamic_cast<TlDenseGeneralMatrix_ImplViennaCLFloat*>(rhs.pImpl_));

    return *this;
}

TlDenseGeneralMatrix_ViennaCLFloat& TlDenseGeneralMatrix_ViennaCLFloat::operator*=(const double coef) {
    *(dynamic_cast<TlDenseGeneralMatrix_ImplViennaCLFloat*>(this->pImpl_)) *= static_cast<float>(coef);

    return *this;
}

TlDenseGeneralMatrix_ViennaCLFloat& TlDenseGeneralMatrix_ViennaCLFloat::operator/=(const double coef) {
    *(dynamic_cast<TlDenseGeneralMatrix_ImplViennaCLFloat*>(this->pImpl_)) /= static_cast<float>(coef);

    return *this;
}

TlDenseGeneralMatrix_ViennaCLFloat& TlDenseGeneralMatrix_ViennaCLFloat::operator*=(const TlDenseGeneralMatrix_ViennaCLFloat& rhs) {
    *(dynamic_cast<TlDenseGeneralMatrix_ImplViennaCLFloat*>(this->pImpl_)) *=
        *(dynamic_cast<TlDenseGeneralMatrix_ImplViennaCLFloat*>(rhs.pImpl_));

    return *this;
}

// ---------------------------------------------------------------------------
// operations
// ---------------------------------------------------------------------------
double TlDenseGeneralMatrix_ViennaCLFloat::sum() const {
    return static_cast<double>(this->pImpl_->sum());
}

double TlDenseGeneralMatrix_ViennaCLFloat::getRMS() const {
    return static_cast<double>(this->pImpl_->getRMS());
}

const TlDenseGeneralMatrix_ViennaCLFloat& TlDenseGeneralMatrix_ViennaCLFloat::dotInPlace(
    const TlDenseGeneralMatrix_ViennaCLFloat& rhs) {
    dynamic_cast<TlDenseGeneralMatrix_ImplViennaCLFloat*>(this->pImpl_)->dotInPlace(*(dynamic_cast<TlDenseGeneralMatrix_ImplViennaCLFloat*>(rhs.pImpl_)));
    return *this;
}

TlDenseGeneralMatrix_ViennaCLFloat TlDenseGeneralMatrix_ViennaCLFloat::transpose() const {
    return TlDenseGeneralMatrix_ViennaCLFloat(dynamic_cast<const TlDenseGeneralMatrix_ImplViennaCLFloat*>(this->pImpl_)->transpose());
}

TlDenseGeneralMatrix_ViennaCLFloat TlDenseGeneralMatrix_ViennaCLFloat::inverse() const {
    return TlDenseGeneralMatrix_ViennaCLFloat(dynamic_cast<const TlDenseGeneralMatrix_ImplViennaCLFloat*>(this->pImpl_)->inverse());
}

TlDenseGeneralMatrix_ViennaCLFloat& TlDenseGeneralMatrix_ViennaCLFloat::reverseColumns() {
    dynamic_cast<TlDenseGeneralMatrix_ImplViennaCLFloat*>(this->pImpl_)->reverseColumns();

    return *this;
}

// ---------------------------------------------------------------------------
// friend functions
// ---------------------------------------------------------------------------
TlDenseVector_ViennaCLFloat operator*(const TlDenseGeneralMatrix_ViennaCLFloat& rhs1,
                                      const TlDenseVector_ViennaCLFloat& rhs2) {
    TlDenseVector_ViennaCLFloat answer;
    *(dynamic_cast<TlDenseVector_ImplViennaCLFloat*>(answer.pImpl_)) =
        *(dynamic_cast<TlDenseGeneralMatrix_ImplViennaCLFloat*>(rhs1.pImpl_)) * *(dynamic_cast<TlDenseVector_ImplViennaCLFloat*>(rhs2.pImpl_));

    return answer;
}

TlDenseVector_ViennaCLFloat operator*(const TlDenseVector_ViennaCLFloat& rhs1,
                                      const TlDenseGeneralMatrix_ViennaCLFloat& rhs2) {
    TlDenseVector_ViennaCLFloat answer;
    *(dynamic_cast<TlDenseVector_ImplViennaCLFloat*>(answer.pImpl_)) =
        *(dynamic_cast<TlDenseVector_ImplViennaCLFloat*>(rhs1.pImpl_)) * *(dynamic_cast<TlDenseGeneralMatrix_ImplViennaCLFloat*>(rhs2.pImpl_));

    return answer;
}

TlDenseGeneralMatrix_ViennaCLFloat operator*(const double coef, const TlDenseGeneralMatrix_ViennaCLFloat& DM) {
    TlDenseGeneralMatrix_ViennaCLFloat answer = DM;
    answer *= static_cast<float>(coef);
    return answer;
}
