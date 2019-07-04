#include "tl_dense_general_matrix_scalapack.h"
#include "tl_dense_general_matrix_impl_scalapack.h"
#include "tl_dense_symmetric_matrix_impl_scalapack.h"
#include "tl_dense_symmetric_matrix_scalapack.h"
#include "tl_dense_vector_impl_scalapack.h"
#include "tl_dense_vector_scalapack.h"

// ---------------------------------------------------------------------------
// constructor & destructor
// ---------------------------------------------------------------------------
TlDenseGeneralMatrix_Scalapack::TlDenseGeneralMatrix_Scalapack(
    const TlMatrixObject::index_type row,
    const TlMatrixObject::index_type col) {
    this->pImpl_ = new TlDenseGeneralMatrix_ImplScalapack(row, col);
}

TlDenseGeneralMatrix_Scalapack::TlDenseGeneralMatrix_Scalapack(
    const TlDenseGeneralMatrix_Scalapack& rhs) {
    this->pImpl_ = new TlDenseGeneralMatrix_ImplScalapack(
        *(dynamic_cast<TlDenseGeneralMatrix_ImplScalapack*>(rhs.pImpl_)));
}

TlDenseGeneralMatrix_Scalapack::TlDenseGeneralMatrix_Scalapack(
    const TlDenseSymmetricMatrix_Scalapack& rhs) {
    this->pImpl_ = new TlDenseGeneralMatrix_ImplScalapack(
        *(dynamic_cast<TlDenseSymmetricMatrix_ImplScalapack*>(rhs.pImpl_)));
}

TlDenseGeneralMatrix_Scalapack::TlDenseGeneralMatrix_Scalapack(
    const TlDenseGeneralMatrix_ImplScalapack& rhs) {
    this->pImpl_ = new TlDenseGeneralMatrix_ImplScalapack(rhs);
}

TlDenseGeneralMatrix_Scalapack::~TlDenseGeneralMatrix_Scalapack() {
    delete this->pImpl_;
    this->pImpl_ = NULL;
}

// ---------------------------------------------------------------------------
// properties
// ---------------------------------------------------------------------------
TlMatrixObject::size_type TlDenseGeneralMatrix_Scalapack::getNumOfElements()
    const {
    return dynamic_cast<TlDenseGeneralMatrix_ImplScalapack*>(this->pImpl_)
        ->getNumOfElements();
}

double TlDenseGeneralMatrix_Scalapack::getLocal(
    TlMatrixObject::index_type row, TlMatrixObject::index_type col) const {
    return dynamic_cast<TlDenseGeneralMatrix_ImplScalapack*>(this->pImpl_)
        ->getLocal(row, col);
}

// ---------------------------------------------------------------------------
// operators
// ---------------------------------------------------------------------------
TlDenseGeneralMatrix_Scalapack& TlDenseGeneralMatrix_Scalapack::operator=(
    const TlDenseGeneralMatrix_Scalapack& rhs) {
    delete this->pImpl_;
    this->pImpl_ = new TlDenseGeneralMatrix_ImplScalapack(
        *(dynamic_cast<TlDenseGeneralMatrix_ImplScalapack*>(rhs.pImpl_)));

    return *this;
}

const TlDenseGeneralMatrix_Scalapack TlDenseGeneralMatrix_Scalapack::operator+(
    const TlDenseGeneralMatrix_Scalapack& rhs) const {
    TlDenseGeneralMatrix_Scalapack answer = *this;
    answer += rhs;
    return answer;
}

const TlDenseGeneralMatrix_Scalapack TlDenseGeneralMatrix_Scalapack::operator-(
    const TlDenseGeneralMatrix_Scalapack& rhs) const {
    TlDenseGeneralMatrix_Scalapack answer = *this;
    answer -= rhs;
    return answer;
}

const TlDenseGeneralMatrix_Scalapack TlDenseGeneralMatrix_Scalapack::operator*(
    const double coef) const {
    TlDenseGeneralMatrix_Scalapack answer = *this;
    answer *= coef;
    return answer;
}

const TlDenseGeneralMatrix_Scalapack TlDenseGeneralMatrix_Scalapack::operator*(
    const TlDenseGeneralMatrix_Scalapack& rhs) const {
    TlDenseGeneralMatrix_Scalapack answer = *this;
    answer *= rhs;
    return answer;
}

TlDenseGeneralMatrix_Scalapack& TlDenseGeneralMatrix_Scalapack::operator+=(
    const TlDenseGeneralMatrix_Scalapack& rhs) {
    *(dynamic_cast<TlDenseGeneralMatrix_ImplScalapack*>(this->pImpl_)) +=
        *(dynamic_cast<TlDenseGeneralMatrix_ImplScalapack*>(rhs.pImpl_));

    return *this;
}

TlDenseGeneralMatrix_Scalapack& TlDenseGeneralMatrix_Scalapack::operator-=(
    const TlDenseGeneralMatrix_Scalapack& rhs) {
    *(dynamic_cast<TlDenseGeneralMatrix_ImplScalapack*>(this->pImpl_)) -=
        *(dynamic_cast<TlDenseGeneralMatrix_ImplScalapack*>(rhs.pImpl_));

    return *this;
}

TlDenseGeneralMatrix_Scalapack& TlDenseGeneralMatrix_Scalapack::operator*=(
    const double coef) {
    *(dynamic_cast<TlDenseGeneralMatrix_ImplScalapack*>(this->pImpl_)) *= coef;

    return *this;
}

TlDenseGeneralMatrix_Scalapack& TlDenseGeneralMatrix_Scalapack::operator/=(
    const double coef) {
    *(dynamic_cast<TlDenseGeneralMatrix_ImplScalapack*>(this->pImpl_)) /= coef;

    return *this;
}

TlDenseGeneralMatrix_Scalapack& TlDenseGeneralMatrix_Scalapack::operator*=(
    const TlDenseGeneralMatrix_Scalapack& rhs) {
    *(dynamic_cast<TlDenseGeneralMatrix_ImplScalapack*>(this->pImpl_)) *=
        *(dynamic_cast<TlDenseGeneralMatrix_ImplScalapack*>(rhs.pImpl_));

    return *this;
}

// ---------------------------------------------------------------------------
// operations
// ---------------------------------------------------------------------------
double TlDenseGeneralMatrix_Scalapack::sum() const {
    return this->pImpl_->sum();
}

double TlDenseGeneralMatrix_Scalapack::getRMS() const {
    return this->pImpl_->getRMS();
}

TlDenseGeneralMatrix_Scalapack TlDenseGeneralMatrix_Scalapack::dot(
    const TlDenseGeneralMatrix_Scalapack& rhs) const {
    TlDenseGeneralMatrix_Scalapack answer(
        dynamic_cast<TlDenseGeneralMatrix_ImplScalapack*>(this->pImpl_)
            ->dot(*(dynamic_cast<TlDenseGeneralMatrix_ImplScalapack*>(
                rhs.pImpl_))));

    return answer;
}

const TlDenseGeneralMatrix_Scalapack&
TlDenseGeneralMatrix_Scalapack::dotInPlace(
    const TlDenseGeneralMatrix_Scalapack& rhs) {
    dynamic_cast<TlDenseGeneralMatrix_ImplScalapack*>(this->pImpl_)
        ->dotInPlace(
            *(dynamic_cast<TlDenseGeneralMatrix_ImplScalapack*>(rhs.pImpl_)));

    return *this;
}

TlDenseGeneralMatrix_Scalapack TlDenseGeneralMatrix_Scalapack::transpose()
    const {
    return TlDenseGeneralMatrix_Scalapack(
        dynamic_cast<const TlDenseGeneralMatrix_ImplScalapack*>(this->pImpl_)
            ->transpose());
}

TlDenseGeneralMatrix_Scalapack TlDenseGeneralMatrix_Scalapack::inverse() const {
    return TlDenseGeneralMatrix_Scalapack(
        dynamic_cast<const TlDenseGeneralMatrix_ImplScalapack*>(this->pImpl_)
            ->inverse());
}

bool TlDenseGeneralMatrix_Scalapack::getSparseMatrix(TlSparseMatrix* pMatrix,
                                                     bool isFinalize) const {
    return dynamic_cast<TlDenseGeneralMatrix_ImplScalapack*>(this->pImpl_)
        ->getSparseMatrix(pMatrix, isFinalize);
}

void TlDenseGeneralMatrix_Scalapack::mergeSparseMatrix(
    const TlSparseMatrix& M) {
    dynamic_cast<TlDenseGeneralMatrix_ImplScalapack*>(this->pImpl_)
        ->mergeSparseMatrix(M);
}

std::vector<TlMatrixObject::index_type>
TlDenseGeneralMatrix_Scalapack::getRowIndexTable() const {
    return dynamic_cast<TlDenseGeneralMatrix_ImplScalapack*>(this->pImpl_)
        ->getRowIndexTable();
}

std::vector<TlMatrixObject::index_type>
TlDenseGeneralMatrix_Scalapack::getColIndexTable() const {
    return dynamic_cast<TlDenseGeneralMatrix_ImplScalapack*>(this->pImpl_)
        ->getColIndexTable();
}

void TlDenseGeneralMatrix_Scalapack::getLocalMatrix(
    TlDenseGeneralMatrixObject* pOutputMatrix) const {
    dynamic_cast<TlDenseGeneralMatrix_ImplScalapack*>(this->pImpl_)
        ->getLocalMatrix(pOutputMatrix);
}

// ---------------------------------------------------------------------------
// I/O
// ---------------------------------------------------------------------------
void TlDenseGeneralMatrix_Scalapack::setUsingPartialIO(const bool isUsePIO) {
    TlDenseGeneralMatrix_ImplScalapack::setUsingPartialIO(isUsePIO);
}

bool TlDenseGeneralMatrix_Scalapack::load(const std::string& filePath) {
    return dynamic_cast<TlDenseGeneralMatrix_ImplScalapack*>(this->pImpl_)
        ->load(filePath);
}

bool TlDenseGeneralMatrix_Scalapack::save(const std::string& filePath) const {
    return dynamic_cast<TlDenseGeneralMatrix_ImplScalapack*>(this->pImpl_)
        ->save(filePath);
}

// ---------------------------------------------------------------------------
// friend functions
// ---------------------------------------------------------------------------
TlDenseVector_Scalapack operator*(const TlDenseGeneralMatrix_Scalapack& rhs1,
                                  const TlDenseVector_Scalapack& rhs2) {
    TlDenseVector_Scalapack answer;
    *(dynamic_cast<TlDenseVector_ImplScalapack*>(answer.pImpl_)) =
        *(dynamic_cast<TlDenseGeneralMatrix_ImplScalapack*>(rhs1.pImpl_)) *
        *(dynamic_cast<TlDenseVector_ImplScalapack*>(rhs2.pImpl_));

    return answer;
}

TlDenseVector_Scalapack operator*(const TlDenseVector_Scalapack& rhs1,
                                  const TlDenseGeneralMatrix_Scalapack& rhs2) {
    TlDenseVector_Scalapack answer;
    *(dynamic_cast<TlDenseVector_ImplScalapack*>(answer.pImpl_)) =
        *(dynamic_cast<TlDenseVector_ImplScalapack*>(rhs1.pImpl_)) *
        *(dynamic_cast<TlDenseGeneralMatrix_ImplScalapack*>(rhs2.pImpl_));

    return answer;
}
