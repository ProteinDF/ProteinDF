#include "tl_dense_vector_viennacl_float.h"

#include "tl_dense_vector_impl_viennacl_float.h"

#ifdef HAVE_EIGEN
#include "tl_dense_vector_eigen_float.h"
#include "tl_dense_vector_impl_eigen_float.h"
#endif  // HAVE_EIGEN

TlDenseVector_ViennaCLFloat::TlDenseVector_ViennaCLFloat(
    const TlDenseVectorObject::index_type size) {
    this->pImpl_ = new TlDenseVector_ImplViennaCLFloat(size);
}

TlDenseVector_ViennaCLFloat::TlDenseVector_ViennaCLFloat(
    const TlDenseVector_ViennaCLFloat& rhs) {
    this->pImpl_ = new TlDenseVector_ImplViennaCLFloat(
        *dynamic_cast<const TlDenseVector_ImplViennaCLFloat*>(rhs.pImpl_));
}

TlDenseVector_ViennaCLFloat::TlDenseVector_ViennaCLFloat(const std::vector<double>& rhs) {
    this->pImpl_ = new TlDenseVector_ImplViennaCLFloat(rhs);
}

#ifdef HAVE_EIGEN
TlDenseVector_ViennaCLFloat::TlDenseVector_ViennaCLFloat(const TlDenseVector_EigenFloat& rhs) {
    this->pImpl_ = new TlDenseVector_ImplViennaCLFloat(*(dynamic_cast<TlDenseVector_ImplEigenFloat*>(rhs.pImpl_)));
}
#endif  // HAVE_EIGEN

TlDenseVector_ViennaCLFloat::operator std::vector<double>() const {
    return std::vector<double>(*(dynamic_cast<TlDenseVector_ImplViennaCLFloat*>(this->pImpl_)));
}

TlDenseVector_ViennaCLFloat::~TlDenseVector_ViennaCLFloat() {
    delete this->pImpl_;
    this->pImpl_ = NULL;
}

// ---------------------------------------------------------------------------
// operators
// ---------------------------------------------------------------------------
TlDenseVector_ViennaCLFloat& TlDenseVector_ViennaCLFloat::operator=(
    const TlDenseVector_ViennaCLFloat& rhs) {
    if (this != &rhs) {
        delete this->pImpl_;
        this->pImpl_ = NULL;

        this->pImpl_ = new TlDenseVector_ImplViennaCLFloat(
            *dynamic_cast<const TlDenseVector_ImplViennaCLFloat*>(rhs.pImpl_));
    }

    return *this;
}

TlDenseVector_ViennaCLFloat& TlDenseVector_ViennaCLFloat::operator+=(
    const TlDenseVector_ViennaCLFloat& rhs) {
    dynamic_cast<TlDenseVector_ImplViennaCLFloat*>(this->pImpl_)
        ->
        operator+=(*dynamic_cast<TlDenseVector_ImplViennaCLFloat*>(rhs.pImpl_));

    return *this;
}

TlDenseVector_ViennaCLFloat& TlDenseVector_ViennaCLFloat::operator-=(
    const TlDenseVector_ViennaCLFloat& rhs) {
    dynamic_cast<TlDenseVector_ImplViennaCLFloat*>(this->pImpl_)
        ->
        operator-=(*dynamic_cast<TlDenseVector_ImplViennaCLFloat*>(rhs.pImpl_));

    return *this;
}

TlDenseVector_ViennaCLFloat& TlDenseVector_ViennaCLFloat::operator*=(const float rhs) {
    dynamic_cast<TlDenseVector_ImplViennaCLFloat*>(this->pImpl_)->operator*=(rhs);

    return *this;
}

TlDenseVector_ViennaCLFloat& TlDenseVector_ViennaCLFloat::operator/=(const float rhs) {
    dynamic_cast<TlDenseVector_ImplViennaCLFloat*>(this->pImpl_)->operator/=(rhs);

    return *this;
}

// ---------------------------------------------------------------------------
// operations
// ---------------------------------------------------------------------------
TlDenseVector_ViennaCLFloat& TlDenseVector_ViennaCLFloat::dotInPlace(
    const TlDenseVector_ViennaCLFloat& rhs) {
    dynamic_cast<TlDenseVector_ImplViennaCLFloat*>(this->pImpl_)
        ->dotInPlace(*dynamic_cast<TlDenseVector_ImplViennaCLFloat*>(rhs.pImpl_));

    return *this;
}

TlDenseVector_ViennaCLFloat& TlDenseVector_ViennaCLFloat::reverse() {
    dynamic_cast<TlDenseVector_ImplViennaCLFloat*>(this->pImpl_)->reverse();

    return *this;
}

// ---------------------------------------------------------------------------
// others
// ---------------------------------------------------------------------------
TlDenseVector_ViennaCLFloat operator+(const TlDenseVector_ViennaCLFloat& rhs1,
                                      const TlDenseVector_ViennaCLFloat& rhs2) {
    TlDenseVector_ViennaCLFloat answer = rhs1;
    answer += rhs2;

    return answer;
}

TlDenseVector_ViennaCLFloat operator-(const TlDenseVector_ViennaCLFloat& rhs1,
                                      const TlDenseVector_ViennaCLFloat& rhs2) {
    TlDenseVector_ViennaCLFloat answer = rhs1;
    answer -= rhs2;

    return answer;
}

TlDenseVector_ViennaCLFloat operator*(const TlDenseVector_ViennaCLFloat& rhs1,
                                      const float rhs2) {
    TlDenseVector_ViennaCLFloat answer = rhs1;
    answer *= rhs2;

    return answer;
}

TlDenseVector_ViennaCLFloat operator*(const float rhs1,
                                      const TlDenseVector_ViennaCLFloat& rhs2) {
    return rhs2 * rhs1;
}
