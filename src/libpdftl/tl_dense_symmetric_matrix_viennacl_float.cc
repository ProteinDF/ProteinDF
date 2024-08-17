#ifdef HAVE_CONFIG_H
#include "config.h"
#endif  // HAVE_CONFIG_H

#include "tl_dense_general_matrix_impl_viennacl_float.h"
#include "tl_dense_general_matrix_viennacl_float.h"
#include "tl_dense_symmetric_matrix_eigen_float.h"
#include "tl_dense_symmetric_matrix_impl_eigen_float.h"
#include "tl_dense_symmetric_matrix_impl_viennacl_float.h"
#include "tl_dense_symmetric_matrix_viennacl_float.h"
#include "tl_dense_vector_impl_viennacl_float.h"
#include "tl_dense_vector_viennacl_float.h"
#include "tl_sparse_symmetric_matrix_impl_viennacl_float.h"
#include "tl_sparse_symmetric_matrix_viennacl_float.h"

TlDenseSymmetricMatrix_ViennaCLFloat::TlDenseSymmetricMatrix_ViennaCLFloat(const TlMatrixObject::index_type dim) {
    this->pImpl_ = new TlDenseSymmetricMatrix_ImplViennaCLFloat(dim);
}

TlDenseSymmetricMatrix_ViennaCLFloat::TlDenseSymmetricMatrix_ViennaCLFloat(const TlDenseSymmetricMatrix_ViennaCLFloat& rhs) {
    this->pImpl_ = new TlDenseSymmetricMatrix_ImplViennaCLFloat(*(
        dynamic_cast<const TlDenseSymmetricMatrix_ImplViennaCLFloat*>(rhs.pImpl_)));
}

TlDenseSymmetricMatrix_ViennaCLFloat::TlDenseSymmetricMatrix_ViennaCLFloat(const TlDenseGeneralMatrix_ViennaCLFloat& rhs) {
    this->pImpl_ = new TlDenseSymmetricMatrix_ImplViennaCLFloat(
        *(dynamic_cast<const TlDenseGeneralMatrix_ImplViennaCLFloat*>(rhs.pImpl_)));
}

TlDenseSymmetricMatrix_ViennaCLFloat::TlDenseSymmetricMatrix_ViennaCLFloat(const TlSparseSymmetricMatrix_ViennaCLFloat& rhs) {
    this->pImpl_ = new TlDenseSymmetricMatrix_ImplViennaCLFloat(
        *(dynamic_cast<TlSparseSymmetricMatrix_ImplViennaCLFloat*>(rhs.pImpl_)));
}

#ifdef HAVE_EIGEN
TlDenseSymmetricMatrix_ViennaCLFloat::TlDenseSymmetricMatrix_ViennaCLFloat(const TlDenseSymmetricMatrix_EigenFloat& rhs) {
    this->pImpl_ = new TlDenseSymmetricMatrix_ImplViennaCLFloat(*dynamic_cast<const TlDenseSymmetricMatrix_ImplEigenFloat*>(rhs.pImpl_));
}
#endif  // HAVE_EIGEN

TlDenseSymmetricMatrix_ViennaCLFloat::~TlDenseSymmetricMatrix_ViennaCLFloat() {
    delete this->pImpl_;
    this->pImpl_ = NULL;
}

void TlDenseSymmetricMatrix_ViennaCLFloat::vtr2mat(const std::vector<double>& vtr) {
    dynamic_cast<TlDenseSymmetricMatrix_ImplViennaCLFloat*>(this->pImpl_)->vtr2mat(vtr);
}

// ---------------------------------------------------------------------------
// operators
// ---------------------------------------------------------------------------
TlDenseSymmetricMatrix_ViennaCLFloat& TlDenseSymmetricMatrix_ViennaCLFloat::operator=(
    const TlDenseSymmetricMatrix_ViennaCLFloat& rhs) {
    if (this != &rhs) {
        delete this->pImpl_;
        this->pImpl_ = new TlDenseSymmetricMatrix_ImplViennaCLFloat(
            *(dynamic_cast<TlDenseSymmetricMatrix_ImplViennaCLFloat*>(rhs.pImpl_)));
    }

    return *this;
}

TlDenseSymmetricMatrix_ViennaCLFloat& TlDenseSymmetricMatrix_ViennaCLFloat::operator=(
    const TlDenseSymmetricMatrix_EigenFloat& rhs) {
    delete this->pImpl_;
    this->pImpl_ = new TlDenseSymmetricMatrix_ImplViennaCLFloat(
        *(dynamic_cast<TlDenseSymmetricMatrix_ImplEigenFloat*>(rhs.pImpl_)));

    return *this;
}

const TlDenseSymmetricMatrix_ViennaCLFloat TlDenseSymmetricMatrix_ViennaCLFloat::
operator+(const TlDenseSymmetricMatrix_ViennaCLFloat& rhs) const {
    TlDenseSymmetricMatrix_ViennaCLFloat answer = *this;
    answer += rhs;
    return answer;
}

const TlDenseSymmetricMatrix_ViennaCLFloat TlDenseSymmetricMatrix_ViennaCLFloat::
operator-(const TlDenseSymmetricMatrix_ViennaCLFloat& rhs) const {
    TlDenseSymmetricMatrix_ViennaCLFloat answer = *this;
    answer -= rhs;
    return answer;
}

const TlDenseSymmetricMatrix_ViennaCLFloat TlDenseSymmetricMatrix_ViennaCLFloat::
operator*(const TlDenseSymmetricMatrix_ViennaCLFloat& rhs) const {
    TlDenseSymmetricMatrix_ViennaCLFloat answer = *this;
    answer *= rhs;
    return answer;
}

TlDenseSymmetricMatrix_ViennaCLFloat& TlDenseSymmetricMatrix_ViennaCLFloat::operator+=(
    const TlDenseSymmetricMatrix_ViennaCLFloat& rhs) {
    *(dynamic_cast<TlDenseSymmetricMatrix_ImplViennaCLFloat*>(this->pImpl_)) +=
        *(dynamic_cast<TlDenseSymmetricMatrix_ImplViennaCLFloat*>(rhs.pImpl_));

    return *this;
}

TlDenseSymmetricMatrix_ViennaCLFloat& TlDenseSymmetricMatrix_ViennaCLFloat::operator-=(
    const TlDenseSymmetricMatrix_ViennaCLFloat& rhs) {
    *(dynamic_cast<TlDenseSymmetricMatrix_ImplViennaCLFloat*>(this->pImpl_)) -=
        *(dynamic_cast<TlDenseSymmetricMatrix_ImplViennaCLFloat*>(rhs.pImpl_));

    return *this;
}

TlDenseSymmetricMatrix_ViennaCLFloat& TlDenseSymmetricMatrix_ViennaCLFloat::operator*=(
    const float coef) {
    *(dynamic_cast<TlDenseSymmetricMatrix_ImplViennaCLFloat*>(this->pImpl_)) *= coef;

    return *this;
}

TlDenseSymmetricMatrix_ViennaCLFloat& TlDenseSymmetricMatrix_ViennaCLFloat::operator/=(
    const float coef) {
    *(dynamic_cast<TlDenseSymmetricMatrix_ImplViennaCLFloat*>(this->pImpl_)) /= coef;

    return *this;
}

TlDenseSymmetricMatrix_ViennaCLFloat& TlDenseSymmetricMatrix_ViennaCLFloat::operator*=(
    const TlDenseSymmetricMatrix_ViennaCLFloat& rhs) {
    *(dynamic_cast<TlDenseSymmetricMatrix_ImplViennaCLFloat*>(this->pImpl_)) *=
        *(dynamic_cast<TlDenseSymmetricMatrix_ImplViennaCLFloat*>(rhs.pImpl_));

    return *this;
}

// ---------------------------------------------------------------------------
// operations
// ---------------------------------------------------------------------------
const TlDenseSymmetricMatrix_ViennaCLFloat&
TlDenseSymmetricMatrix_ViennaCLFloat::dotInPlace(const TlDenseSymmetricMatrix_ViennaCLFloat& rhs) {
    dynamic_cast<TlDenseSymmetricMatrix_ImplViennaCLFloat*>(this->pImpl_)->dotInPlace(*dynamic_cast<TlDenseSymmetricMatrix_ImplViennaCLFloat*>(rhs.pImpl_));

    return *this;
}

bool TlDenseSymmetricMatrix_ViennaCLFloat::eig(TlDenseVector_ViennaCLFloat* pEigVal, TlDenseGeneralMatrix_ViennaCLFloat* pEigVec,
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

TlDenseSymmetricMatrix_ViennaCLFloat TlDenseSymmetricMatrix_ViennaCLFloat::inverse()
    const {
    TlDenseSymmetricMatrix_ViennaCLFloat answer;
    answer.pImpl_ = new TlDenseSymmetricMatrix_ImplViennaCLFloat(
        dynamic_cast<const TlDenseSymmetricMatrix_ImplViennaCLFloat*>(this->pImpl_)
            ->inverse());
    return answer;
}

// ---------------------------------------------------------------------------
// I/O
// ---------------------------------------------------------------------------
bool TlDenseSymmetricMatrix_ViennaCLFloat::load(const std::string& filePath) {
    TlDenseSymmetricMatrix_EigenFloat eigenmat;
    const bool answer = eigenmat.load(filePath);

    if (answer) {
        delete this->pImpl_;
        this->pImpl_ = NULL;

        this->pImpl_ = new TlDenseSymmetricMatrix_ImplViennaCLFloat(
            *dynamic_cast<const TlDenseSymmetricMatrix_ImplEigenFloat*>(
                eigenmat.pImpl_));
    }

    return answer;
}

bool TlDenseSymmetricMatrix_ViennaCLFloat::save(const std::string& filePath) const {
    TlDenseSymmetricMatrix_EigenFloat eigenmat(*this);
    return eigenmat.save(filePath);
}

// ---------------------------------------------------------------------------
// protected
// ---------------------------------------------------------------------------
bool TlDenseSymmetricMatrix_ViennaCLFloat::eig_powerIteration(TlDenseVector_ViennaCLFloat* pEigVal,
                                                              TlDenseGeneralMatrix_ViennaCLFloat* pEigVec) const {
    TlDenseVector_ImplViennaCLFloat* pImpl_eigval = dynamic_cast<TlDenseVector_ImplViennaCLFloat*>(pEigVal->pImpl_);
    TlDenseGeneralMatrix_ImplViennaCLFloat* pImpl_eigvec = dynamic_cast<TlDenseGeneralMatrix_ImplViennaCLFloat*>(pEigVec->pImpl_);

    const bool answer = dynamic_cast<TlDenseSymmetricMatrix_ImplViennaCLFloat*>(this->pImpl_)->eig(pImpl_eigval, pImpl_eigvec);
    return answer;
}

bool TlDenseSymmetricMatrix_ViennaCLFloat::eig_QR(TlDenseVector_ViennaCLFloat* pEigVal,
                                                  TlDenseGeneralMatrix_ViennaCLFloat* pEigVec) const {
    TlDenseVector_ImplViennaCLFloat* pImpl_eigval = dynamic_cast<TlDenseVector_ImplViennaCLFloat*>(pEigVal->pImpl_);
    TlDenseGeneralMatrix_ImplViennaCLFloat* pImpl_eigvec = dynamic_cast<TlDenseGeneralMatrix_ImplViennaCLFloat*>(pEigVec->pImpl_);

    const bool answer = dynamic_cast<TlDenseSymmetricMatrix_ImplViennaCLFloat*>(this->pImpl_)->eig_QR(pImpl_eigval, pImpl_eigvec);
    return answer;
}

// ---------------------------------------------------------------------------
// friend functions
// ---------------------------------------------------------------------------
TlDenseGeneralMatrix_ViennaCLFloat operator*(const TlDenseGeneralMatrix_ViennaCLFloat& mat1,
                                             const TlDenseSymmetricMatrix_ViennaCLFloat& mat2) {
    return TlDenseGeneralMatrix_ViennaCLFloat(*dynamic_cast<TlDenseGeneralMatrix_ImplViennaCLFloat*>(mat1.pImpl_) *
                                              *dynamic_cast<TlDenseSymmetricMatrix_ImplViennaCLFloat*>(mat2.pImpl_));
}

TlDenseGeneralMatrix_ViennaCLFloat operator*(const TlDenseSymmetricMatrix_ViennaCLFloat& mat1,
                                             const TlDenseGeneralMatrix_ViennaCLFloat& mat2) {
    return TlDenseGeneralMatrix_ViennaCLFloat(*dynamic_cast<TlDenseSymmetricMatrix_ImplViennaCLFloat*>(mat1.pImpl_) *
                                              *dynamic_cast<TlDenseGeneralMatrix_ImplViennaCLFloat*>(mat2.pImpl_));
}

// ---------------------------------------------------------------------------
// arithmetic
// ---------------------------------------------------------------------------
TlDenseSymmetricMatrix_ViennaCLFloat operator*(const float coef, const TlDenseSymmetricMatrix_ViennaCLFloat& DM) {
    TlDenseSymmetricMatrix_ViennaCLFloat answer = DM;
    answer *= coef;
    return answer;
}

TlDenseSymmetricMatrix_ViennaCLFloat operator*(const TlDenseSymmetricMatrix_ViennaCLFloat& DM, const float coef) {
    return coef * DM;
}
