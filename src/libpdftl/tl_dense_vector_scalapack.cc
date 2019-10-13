#include "tl_dense_vector_scalapack.h"
#include "tl_dense_vector_impl_lapack.h"
#include "tl_dense_vector_impl_scalapack.h"
#include "tl_dense_vector_lapack.h"

// ---------------------------------------------------------------------------
// constructor & destructor
// ---------------------------------------------------------------------------
TlDenseVector_Scalapack::TlDenseVector_Scalapack(
    TlDenseVectorObject::index_type size) {
    this->pImpl_ = new TlDenseVector_ImplScalapack(size);
}

TlDenseVector_Scalapack::TlDenseVector_Scalapack(
    const TlDenseVector_Scalapack& rhs) {
    this->pImpl_ = new TlDenseVector_ImplScalapack(
        *dynamic_cast<const TlDenseVector_ImplScalapack*>(rhs.pImpl_));
}

TlDenseVector_Scalapack::TlDenseVector_Scalapack(
    const TlDenseVector_Lapack& rhs) {
    this->pImpl_ = new TlDenseVector_ImplScalapack(
        *dynamic_cast<const TlDenseVector_ImplLapack*>(rhs.pImpl_));
}

TlDenseVector_Scalapack::~TlDenseVector_Scalapack() {
    delete this->pImpl_;
    this->pImpl_ = NULL;
}

// ---------------------------------------------------------------------------
// properties
// ---------------------------------------------------------------------------
std::vector<double> TlDenseVector_Scalapack::getVector() const {
    return dynamic_cast<const TlDenseVector_ImplScalapack*>(this->pImpl_)
        ->getVector();
}

// ---------------------------------------------------------------------------
// operators
// ---------------------------------------------------------------------------
TlDenseVector_Scalapack& TlDenseVector_Scalapack::operator=(
    const TlDenseVector_Scalapack& rhs) {
    if (this != &rhs) {
        delete this->pImpl_;
        this->pImpl_ = NULL;

        this->pImpl_ = new TlDenseVector_ImplScalapack(
            *dynamic_cast<const TlDenseVector_ImplScalapack*>(rhs.pImpl_));
    }

    return *this;
}

TlDenseVector_Scalapack& TlDenseVector_Scalapack::operator+=(
    const TlDenseVector_Scalapack& rhs) {
    dynamic_cast<TlDenseVector_ImplScalapack*>(this->pImpl_)
        ->
        operator+=(*dynamic_cast<TlDenseVector_ImplScalapack*>(rhs.pImpl_));

    return *this;
}

TlDenseVector_Scalapack& TlDenseVector_Scalapack::operator-=(
    const TlDenseVector_Scalapack& rhs) {
    dynamic_cast<TlDenseVector_ImplScalapack*>(this->pImpl_)
        ->
        operator-=(*dynamic_cast<TlDenseVector_ImplScalapack*>(rhs.pImpl_));

    return *this;
}

TlDenseVector_Scalapack& TlDenseVector_Scalapack::operator*=(const double rhs) {
    dynamic_cast<TlDenseVector_ImplScalapack*>(this->pImpl_)->operator*=(rhs);

    return *this;
}

TlDenseVector_Scalapack& TlDenseVector_Scalapack::operator/=(const double rhs) {
    dynamic_cast<TlDenseVector_ImplScalapack*>(this->pImpl_)->operator/=(rhs);

    return *this;
}

double TlDenseVector_Scalapack::operator*(
    const TlDenseVector_Scalapack& rhs) const {
    return dynamic_cast<const TlDenseVector_ImplScalapack*>(this->pImpl_)
        ->
        operator*(
            *dynamic_cast<const TlDenseVector_ImplScalapack*>(rhs.pImpl_));
}

// ---------------------------------------------------------------------------
// operations
// ---------------------------------------------------------------------------
double TlDenseVector_Scalapack::dot(const TlDenseVector_Scalapack& rhs) const {
    double answer =
        dynamic_cast<TlDenseVector_ImplScalapack*>(this->pImpl_)
            ->dot(*dynamic_cast<TlDenseVector_ImplScalapack*>(rhs.pImpl_));

    return answer;
}

TlDenseVector_Scalapack& TlDenseVector_Scalapack::dotInPlace(
    const TlDenseVector_Scalapack& rhs) {
    dynamic_cast<TlDenseVector_ImplScalapack*>(this->pImpl_)
        ->dotInPlace(*dynamic_cast<TlDenseVector_ImplScalapack*>(rhs.pImpl_));

    return *this;
}

// double* TlDenseVector_Scalapack::data() {
//   return dynamic_cast<TlDenseVector_ImplScalapack*>(this->pImpl_)->data();
// }
//
// const double* TlDenseVector_Scalapack::data() const {
//   return dynamic_cast<TlDenseVector_ImplScalapack*>(this->pImpl_)->data();
// }

// ---------------------------------------------------------------------------
// I/O
// ---------------------------------------------------------------------------
bool TlDenseVector_Scalapack::load(const std::string& filePath) {
    bool answer = dynamic_cast<TlDenseVector_ImplScalapack*>(this->pImpl_)
                      ->load(filePath);
    return answer;
}

bool TlDenseVector_Scalapack::save(const std::string& filePath) const {
    bool answer = dynamic_cast<TlDenseVector_ImplScalapack*>(this->pImpl_)
                      ->save(filePath);
    return answer;
}

// ---------------------------------------------------------------------------
// others
// ---------------------------------------------------------------------------
TlDenseVector_Scalapack operator+(const TlDenseVector_Scalapack& rhs1,
                                  const TlDenseVector_Scalapack& rhs2) {
    TlDenseVector_Scalapack answer = rhs1;
    answer += rhs2;

    return answer;
}

TlDenseVector_Scalapack operator-(const TlDenseVector_Scalapack& rhs1,
                                  const TlDenseVector_Scalapack& rhs2) {
    TlDenseVector_Scalapack answer = rhs1;
    answer -= rhs2;

    return answer;
}

TlDenseVector_Scalapack operator*(const TlDenseVector_Scalapack& rhs1,
                                  const double rhs2) {
    TlDenseVector_Scalapack answer = rhs1;
    answer *= rhs2;

    return answer;
}

TlDenseVector_Scalapack operator*(const double rhs1,
                                  const TlDenseVector_Scalapack& rhs2) {
    return rhs2 * rhs1;
}
