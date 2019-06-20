#include <cassert>
#include <vector>

#include "tl_dense_general_matrix_lapack.h"
#include "tl_dense_symmetric_matrix_impl_lapack.h"
#include "tl_dense_symmetric_matrix_lapack.h"
#include "tl_dense_vector_impl_lapack.h"
#include "tl_dense_vector_lapack.h"

// ---------------------------------------------------------------------------
// constructor & destructor
// ---------------------------------------------------------------------------
TlDenseGeneralMatrix_Lapack::TlDenseGeneralMatrix_Lapack(
    const TlMatrixObject::index_type row,
    const TlMatrixObject::index_type col) {
    this->pImpl_ = new TlDenseGeneralMatrix_ImplLapack(row, col);
}

TlDenseGeneralMatrix_Lapack::TlDenseGeneralMatrix_Lapack(
    const TlDenseGeneralMatrix_Lapack& rhs) {
    this->pImpl_ = new TlDenseGeneralMatrix_ImplLapack(
        *(dynamic_cast<TlDenseGeneralMatrix_ImplLapack*>(rhs.pImpl_)));
}

TlDenseGeneralMatrix_Lapack::TlDenseGeneralMatrix_Lapack(
    const TlDenseSymmetricMatrix_Lapack& rhs) {
    this->pImpl_ = new TlDenseGeneralMatrix_ImplLapack(
        *(dynamic_cast<const TlDenseSymmetricMatrix_ImplLapack*>(rhs.pImpl_)));
}

TlDenseGeneralMatrix_Lapack::TlDenseGeneralMatrix_Lapack(
    const TlDenseGeneralMatrix_ImplLapack& rhs) {
    this->pImpl_ = new TlDenseGeneralMatrix_ImplLapack(rhs);
}

void TlDenseGeneralMatrix_Lapack::vtr2mat(const std::vector<double>& vtr) {
    assert(vtr.size() == this->getNumOfRows() * this->getNumOfCols());
    dynamic_cast<TlDenseGeneralMatrix_ImplLapack*>(this->pImpl_)->vtr2mat(vtr);
}

TlDenseGeneralMatrix_Lapack::~TlDenseGeneralMatrix_Lapack() {
    delete this->pImpl_;
    this->pImpl_ = NULL;
}

// ---------------------------------------------------------------------------
// properties
// ---------------------------------------------------------------------------
TlMatrixObject::size_type TlDenseGeneralMatrix_Lapack::getNumOfElements()
    const {
    return dynamic_cast<TlDenseGeneralMatrix_ImplLapack*>(this->pImpl_)
        ->getNumOfElements();
}

// ---------------------------------------------------------------------------
// operators
// ---------------------------------------------------------------------------
TlDenseGeneralMatrix_Lapack& TlDenseGeneralMatrix_Lapack::operator=(
    const TlDenseGeneralMatrix_Lapack& rhs) {
    delete this->pImpl_;
    this->pImpl_ = new TlDenseGeneralMatrix_ImplLapack(
        *(dynamic_cast<TlDenseGeneralMatrix_ImplLapack*>(rhs.pImpl_)));

    return *this;
}

const TlDenseGeneralMatrix_Lapack TlDenseGeneralMatrix_Lapack::operator+(
    const TlDenseGeneralMatrix_Lapack& rhs) const {
    TlDenseGeneralMatrix_Lapack answer = *this;
    answer += rhs;
    return answer;
}

const TlDenseGeneralMatrix_Lapack TlDenseGeneralMatrix_Lapack::operator-(
    const TlDenseGeneralMatrix_Lapack& rhs) const {
    TlDenseGeneralMatrix_Lapack answer = *this;
    answer -= rhs;
    return answer;
}

const TlDenseGeneralMatrix_Lapack TlDenseGeneralMatrix_Lapack::operator*(
    const double coef) const {
    TlDenseGeneralMatrix_Lapack answer = *this;
    answer *= coef;
    return answer;
}

const TlDenseGeneralMatrix_Lapack TlDenseGeneralMatrix_Lapack::operator*(
    const TlDenseGeneralMatrix_Lapack& rhs) const {
    TlDenseGeneralMatrix_Lapack answer = *this;
    answer *= rhs;
    return answer;
}

TlDenseGeneralMatrix_Lapack& TlDenseGeneralMatrix_Lapack::operator+=(
    const TlDenseGeneralMatrix_Lapack& rhs) {
    *(dynamic_cast<TlDenseGeneralMatrix_ImplLapack*>(this->pImpl_)) +=
        *(dynamic_cast<TlDenseGeneralMatrix_ImplLapack*>(rhs.pImpl_));

    return *this;
}

TlDenseGeneralMatrix_Lapack& TlDenseGeneralMatrix_Lapack::operator-=(
    const TlDenseGeneralMatrix_Lapack& rhs) {
    *(dynamic_cast<TlDenseGeneralMatrix_ImplLapack*>(this->pImpl_)) -=
        *(dynamic_cast<TlDenseGeneralMatrix_ImplLapack*>(rhs.pImpl_));

    return *this;
}

TlDenseGeneralMatrix_Lapack& TlDenseGeneralMatrix_Lapack::operator*=(
    const double coef) {
    *(dynamic_cast<TlDenseGeneralMatrix_ImplLapack*>(this->pImpl_)) *= coef;

    return *this;
}

TlDenseGeneralMatrix_Lapack& TlDenseGeneralMatrix_Lapack::operator/=(
    const double coef) {
    *(dynamic_cast<TlDenseGeneralMatrix_ImplLapack*>(this->pImpl_)) /= coef;

    return *this;
}

TlDenseGeneralMatrix_Lapack& TlDenseGeneralMatrix_Lapack::operator*=(
    const TlDenseGeneralMatrix_Lapack& rhs) {
    *(dynamic_cast<TlDenseGeneralMatrix_ImplLapack*>(this->pImpl_)) *=
        *(dynamic_cast<TlDenseGeneralMatrix_ImplLapack*>(rhs.pImpl_));

    return *this;
}

// ---------------------------------------------------------------------------
// operations
// ---------------------------------------------------------------------------
double TlDenseGeneralMatrix_Lapack::sum() const { return this->pImpl_->sum(); }

double TlDenseGeneralMatrix_Lapack::getRMS() const {
    return this->pImpl_->getRMS();
}

TlDenseGeneralMatrix_Lapack TlDenseGeneralMatrix_Lapack::dot(
    const TlDenseGeneralMatrix_Lapack& rhs) const {
    TlDenseGeneralMatrix_Lapack answer(
        dynamic_cast<TlDenseGeneralMatrix_ImplLapack*>(this->pImpl_)
            ->dot(
                *(dynamic_cast<TlDenseGeneralMatrix_ImplLapack*>(rhs.pImpl_))));

    return answer;
}

const TlDenseGeneralMatrix_Lapack& TlDenseGeneralMatrix_Lapack::dotInPlace(
    const TlDenseGeneralMatrix_Lapack& rhs) {
    dynamic_cast<TlDenseGeneralMatrix_ImplLapack*>(this->pImpl_)
        ->dotInPlace(
            *(dynamic_cast<TlDenseGeneralMatrix_ImplLapack*>(rhs.pImpl_)));

    return *this;
}

TlDenseGeneralMatrix_Lapack TlDenseGeneralMatrix_Lapack::transpose() const {
    return TlDenseGeneralMatrix_Lapack(
        dynamic_cast<const TlDenseGeneralMatrix_ImplLapack*>(this->pImpl_)
            ->transpose());
}

TlDenseGeneralMatrix_Lapack TlDenseGeneralMatrix_Lapack::inverse() const {
    return TlDenseGeneralMatrix_Lapack(
        dynamic_cast<const TlDenseGeneralMatrix_ImplLapack*>(this->pImpl_)
            ->inverse());
}

TlDenseGeneralMatrix_Lapack
TlDenseGeneralMatrix_Lapack::getLeastSquaresSolution(
    const TlDenseGeneralMatrix_Lapack& B) const {
    return TlDenseGeneralMatrix_Lapack(
        dynamic_cast<const TlDenseGeneralMatrix_ImplLapack*>(this->pImpl_)
            ->getLeastSquaresSolution(
                *(dynamic_cast<const TlDenseGeneralMatrix_ImplLapack*>(
                    B.pImpl_))));
}

double* TlDenseGeneralMatrix_Lapack::data() {
    return dynamic_cast<TlDenseGeneralMatrix_ImplLapack*>(this->pImpl_)->data();
}

const double* TlDenseGeneralMatrix_Lapack::data() const {
    return dynamic_cast<TlDenseGeneralMatrix_ImplLapack*>(this->pImpl_)->data();
}

// ---------------------------------------------------------------------------
// I/O
// ---------------------------------------------------------------------------
void TlDenseGeneralMatrix_Lapack::dump(TlDenseVector_Lapack* v) const {
    dynamic_cast<TlDenseGeneralMatrix_ImplLapack*>(this->pImpl_)
        ->dump(dynamic_cast<TlDenseVector_ImplLapack*>(v->pImpl_));
}

void TlDenseGeneralMatrix_Lapack::restore(const TlDenseVector_Lapack& v) {
    dynamic_cast<TlDenseGeneralMatrix_ImplLapack*>(this->pImpl_)
        ->restore(*(dynamic_cast<TlDenseVector_ImplLapack*>(v.pImpl_)));
}

// ---------------------------------------------------------------------------
// friend functions
// ---------------------------------------------------------------------------
TlDenseVector_Lapack operator*(const TlDenseGeneralMatrix_Lapack& rhs1,
                               const TlDenseVector_Lapack& rhs2) {
    TlDenseVector_Lapack answer;
    *(dynamic_cast<TlDenseVector_ImplLapack*>(answer.pImpl_)) =
        *(dynamic_cast<TlDenseGeneralMatrix_ImplLapack*>(rhs1.pImpl_)) *
        *(dynamic_cast<TlDenseVector_ImplLapack*>(rhs2.pImpl_));

    return answer;
}

TlDenseVector_Lapack operator*(const TlDenseVector_Lapack& rhs1,
                               const TlDenseGeneralMatrix_Lapack& rhs2) {
    TlDenseVector_Lapack answer;
    *(dynamic_cast<TlDenseVector_ImplLapack*>(answer.pImpl_)) =
        *(dynamic_cast<TlDenseVector_ImplLapack*>(rhs1.pImpl_)) *
        *(dynamic_cast<TlDenseGeneralMatrix_ImplLapack*>(rhs2.pImpl_));

    return answer;
}
