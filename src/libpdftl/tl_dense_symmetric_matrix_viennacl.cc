#ifdef HAVE_CONFIG_H
#include "config.h"
#endif  // HAVE_CONFIG_H

#include "tl_dense_general_matrix_impl_viennacl.h"
#include "tl_dense_general_matrix_viennacl.h"
#include "tl_dense_symmetric_matrix_eigen.h"
#include "tl_dense_symmetric_matrix_impl_eigen.h"
#include "tl_dense_symmetric_matrix_impl_viennacl.h"
#include "tl_dense_symmetric_matrix_viennacl.h"
#include "tl_dense_vector_impl_viennacl.h"
#include "tl_dense_vector_viennacl.h"
#include "tl_sparse_symmetric_matrix_impl_viennacl.h"
#include "tl_sparse_symmetric_matrix_viennacl.h"

TlDenseSymmetricMatrix_ViennaCL::TlDenseSymmetricMatrix_ViennaCL(
    const TlMatrixObject::index_type dim) {
    this->pImpl_ = new TlDenseSymmetricMatrix_ImplViennaCL(dim);
}

TlDenseSymmetricMatrix_ViennaCL::TlDenseSymmetricMatrix_ViennaCL(
    const TlDenseSymmetricMatrix_ViennaCL& rhs) {
    this->pImpl_ = new TlDenseSymmetricMatrix_ImplViennaCL(*(
        dynamic_cast<const TlDenseSymmetricMatrix_ImplViennaCL*>(rhs.pImpl_)));
}

TlDenseSymmetricMatrix_ViennaCL::TlDenseSymmetricMatrix_ViennaCL(
    const TlDenseGeneralMatrix_ViennaCL& rhs) {
    this->pImpl_ = new TlDenseSymmetricMatrix_ImplViennaCL(
        *(dynamic_cast<const TlDenseGeneralMatrix_ImplViennaCL*>(rhs.pImpl_)));
}

TlDenseSymmetricMatrix_ViennaCL::TlDenseSymmetricMatrix_ViennaCL(
    const TlSparseSymmetricMatrix_ViennaCL& rhs) {
    this->pImpl_ = new TlDenseSymmetricMatrix_ImplViennaCL(
        *(dynamic_cast<TlSparseSymmetricMatrix_ImplViennaCL*>(rhs.pImpl_)));
}

#ifdef HAVE_EIGEN
TlDenseSymmetricMatrix_ViennaCL::TlDenseSymmetricMatrix_ViennaCL(
    const TlDenseSymmetricMatrix_Eigen& rhs) {
    this->pImpl_ = new TlDenseSymmetricMatrix_ImplViennaCL(
        *dynamic_cast<const TlDenseSymmetricMatrix_ImplEigen*>(rhs.pImpl_));
}
#endif  // HAVE_EIGEN

TlDenseSymmetricMatrix_ViennaCL::~TlDenseSymmetricMatrix_ViennaCL() {
    delete this->pImpl_;
    this->pImpl_ = NULL;
}

void TlDenseSymmetricMatrix_ViennaCL::vtr2mat(const std::vector<double>& vtr) {
    dynamic_cast<TlDenseSymmetricMatrix_ImplViennaCL*>(this->pImpl_)
        ->vtr2mat(vtr);
}

// ---------------------------------------------------------------------------
// operators
// ---------------------------------------------------------------------------
TlDenseSymmetricMatrix_ViennaCL& TlDenseSymmetricMatrix_ViennaCL::operator=(
    const TlDenseSymmetricMatrix_ViennaCL& rhs) {
    if (this != &rhs) {
        delete this->pImpl_;
        this->pImpl_ = new TlDenseSymmetricMatrix_ImplViennaCL(
            *(dynamic_cast<TlDenseSymmetricMatrix_ImplViennaCL*>(rhs.pImpl_)));
    }

    return *this;
}

TlDenseSymmetricMatrix_ViennaCL& TlDenseSymmetricMatrix_ViennaCL::operator=(
    const TlDenseSymmetricMatrix_Eigen& rhs) {
    delete this->pImpl_;
    this->pImpl_ = new TlDenseSymmetricMatrix_ImplViennaCL(
        *(dynamic_cast<TlDenseSymmetricMatrix_ImplEigen*>(rhs.pImpl_)));

    return *this;
}

const TlDenseSymmetricMatrix_ViennaCL TlDenseSymmetricMatrix_ViennaCL::
operator+(const TlDenseSymmetricMatrix_ViennaCL& rhs) const {
    TlDenseSymmetricMatrix_ViennaCL answer = *this;
    answer += rhs;
    return answer;
}

const TlDenseSymmetricMatrix_ViennaCL TlDenseSymmetricMatrix_ViennaCL::
operator-(const TlDenseSymmetricMatrix_ViennaCL& rhs) const {
    TlDenseSymmetricMatrix_ViennaCL answer = *this;
    answer -= rhs;
    return answer;
}

const TlDenseSymmetricMatrix_ViennaCL TlDenseSymmetricMatrix_ViennaCL::
operator*(const TlDenseSymmetricMatrix_ViennaCL& rhs) const {
    TlDenseSymmetricMatrix_ViennaCL answer = *this;
    answer *= rhs;
    return answer;
}

TlDenseSymmetricMatrix_ViennaCL& TlDenseSymmetricMatrix_ViennaCL::operator+=(
    const TlDenseSymmetricMatrix_ViennaCL& rhs) {
    *(dynamic_cast<TlDenseSymmetricMatrix_ImplViennaCL*>(this->pImpl_)) +=
        *(dynamic_cast<TlDenseSymmetricMatrix_ImplViennaCL*>(rhs.pImpl_));

    return *this;
}

TlDenseSymmetricMatrix_ViennaCL& TlDenseSymmetricMatrix_ViennaCL::operator-=(
    const TlDenseSymmetricMatrix_ViennaCL& rhs) {
    *(dynamic_cast<TlDenseSymmetricMatrix_ImplViennaCL*>(this->pImpl_)) -=
        *(dynamic_cast<TlDenseSymmetricMatrix_ImplViennaCL*>(rhs.pImpl_));

    return *this;
}

TlDenseSymmetricMatrix_ViennaCL& TlDenseSymmetricMatrix_ViennaCL::operator*=(
    const double coef) {
    *(dynamic_cast<TlDenseSymmetricMatrix_ImplViennaCL*>(this->pImpl_)) *= coef;

    return *this;
}

TlDenseSymmetricMatrix_ViennaCL& TlDenseSymmetricMatrix_ViennaCL::operator/=(
    const double coef) {
    *(dynamic_cast<TlDenseSymmetricMatrix_ImplViennaCL*>(this->pImpl_)) /= coef;

    return *this;
}

TlDenseSymmetricMatrix_ViennaCL& TlDenseSymmetricMatrix_ViennaCL::operator*=(
    const TlDenseSymmetricMatrix_ViennaCL& rhs) {
    *(dynamic_cast<TlDenseSymmetricMatrix_ImplViennaCL*>(this->pImpl_)) *=
        *(dynamic_cast<TlDenseSymmetricMatrix_ImplViennaCL*>(rhs.pImpl_));

    return *this;
}

// ---------------------------------------------------------------------------
// operations
// ---------------------------------------------------------------------------
const TlDenseSymmetricMatrix_ViennaCL&
TlDenseSymmetricMatrix_ViennaCL::dotInPlace(
    const TlDenseSymmetricMatrix_ViennaCL& rhs) {
    dynamic_cast<TlDenseSymmetricMatrix_ImplViennaCL*>(this->pImpl_)
        ->dotInPlace(
            *dynamic_cast<TlDenseSymmetricMatrix_ImplViennaCL*>(rhs.pImpl_));

    return *this;
}

bool TlDenseSymmetricMatrix_ViennaCL::eig(
    TlDenseVector_ViennaCL* pEigVal, TlDenseGeneralMatrix_ViennaCL* pEigVec,
    EIG_METHOD eigMethod) const {
    bool answer = false;
    switch (eigMethod) {
        case EIG_QR:
            answer = this->eig_QR(pEigVal, pEigVec);
            break;

        case EIG_POWERITERATION:
            answer = this->eig_powerIteration(pEigVal, pEigVec);
            break;

        default:
            answer = this->eig_QR(pEigVal, pEigVec);
            break;
    }

    return answer;
}

TlDenseSymmetricMatrix_ViennaCL TlDenseSymmetricMatrix_ViennaCL::inverse()
    const {
    TlDenseSymmetricMatrix_ViennaCL answer;
    answer.pImpl_ = new TlDenseSymmetricMatrix_ImplViennaCL(
        dynamic_cast<const TlDenseSymmetricMatrix_ImplViennaCL*>(this->pImpl_)
            ->inverse());
    return answer;
}

// ---------------------------------------------------------------------------
// I/O
// ---------------------------------------------------------------------------
bool TlDenseSymmetricMatrix_ViennaCL::load(const std::string& filePath) {
    TlDenseSymmetricMatrix_Eigen eigenmat;
    const bool answer = eigenmat.load(filePath);

    if (answer) {
        delete this->pImpl_;
        this->pImpl_ = NULL;

        this->pImpl_ = new TlDenseSymmetricMatrix_ImplViennaCL(
            *dynamic_cast<const TlDenseSymmetricMatrix_ImplEigen*>(
                eigenmat.pImpl_));
    }

    return answer;
}

bool TlDenseSymmetricMatrix_ViennaCL::save(const std::string& filePath) const {
    TlDenseSymmetricMatrix_Eigen eigenmat(*this);
    return eigenmat.save(filePath);
}

// ---------------------------------------------------------------------------
// protected
// ---------------------------------------------------------------------------
bool TlDenseSymmetricMatrix_ViennaCL::eig_powerIteration(
    TlDenseVector_ViennaCL* pEigVal,
    TlDenseGeneralMatrix_ViennaCL* pEigVec) const {
    TlDenseVector_ImplViennaCL* pImpl_eigval =
        dynamic_cast<TlDenseVector_ImplViennaCL*>(pEigVal->pImpl_);
    TlDenseGeneralMatrix_ImplViennaCL* pImpl_eigvec =
        dynamic_cast<TlDenseGeneralMatrix_ImplViennaCL*>(pEigVec->pImpl_);

    const bool answer =
        dynamic_cast<TlDenseSymmetricMatrix_ImplViennaCL*>(this->pImpl_)
            ->eig(pImpl_eigval, pImpl_eigvec);
    return answer;
}

bool TlDenseSymmetricMatrix_ViennaCL::eig_QR(
    TlDenseVector_ViennaCL* pEigVal,
    TlDenseGeneralMatrix_ViennaCL* pEigVec) const {
    TlDenseVector_ImplViennaCL* pImpl_eigval =
        dynamic_cast<TlDenseVector_ImplViennaCL*>(pEigVal->pImpl_);
    TlDenseGeneralMatrix_ImplViennaCL* pImpl_eigvec =
        dynamic_cast<TlDenseGeneralMatrix_ImplViennaCL*>(pEigVec->pImpl_);

    const bool answer =
        dynamic_cast<TlDenseSymmetricMatrix_ImplViennaCL*>(this->pImpl_)
            ->eig_QR(pImpl_eigval, pImpl_eigvec);
    return answer;
}

// ---------------------------------------------------------------------------
// friend functions
// ---------------------------------------------------------------------------
TlDenseGeneralMatrix_ViennaCL operator*(
    const TlDenseGeneralMatrix_ViennaCL& mat1,
    const TlDenseSymmetricMatrix_ViennaCL& mat2) {
    return TlDenseGeneralMatrix_ViennaCL(
        *dynamic_cast<TlDenseGeneralMatrix_ImplViennaCL*>(mat1.pImpl_) *
        *dynamic_cast<TlDenseSymmetricMatrix_ImplViennaCL*>(mat2.pImpl_));
}

TlDenseGeneralMatrix_ViennaCL operator*(
    const TlDenseSymmetricMatrix_ViennaCL& mat1,
    const TlDenseGeneralMatrix_ViennaCL& mat2) {
    return TlDenseGeneralMatrix_ViennaCL(
        *dynamic_cast<TlDenseSymmetricMatrix_ImplViennaCL*>(mat1.pImpl_) *
        *dynamic_cast<TlDenseGeneralMatrix_ImplViennaCL*>(mat2.pImpl_));
}

TlDenseSymmetricMatrix_ViennaCL operator*(
    const double coef, const TlDenseSymmetricMatrix_ViennaCL& DM) {
    TlDenseSymmetricMatrix_ViennaCL answer = DM;
    answer *= coef;
    return answer;
}
